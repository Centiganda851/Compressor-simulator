"""Compressor simulation library converted from MATLAB model.

This module ports the provided MATLAB two-port reciprocating compressor model
into pure Python using SciPy's stiff ODE solver. The primary entry point is
`simulate(config)` which returns crank-angle-resolved cylinder pressure,
temperature, reed lift, and flow history.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict
import numpy as np
from scipy.integrate import solve_ivp


# Unit conversion constants.
PSI_TO_PA = 6894.75729  # 1 psi in pascals
PA_TO_PSI = 1 / PSI_TO_PA  # 1 pascal in psi


@dataclass
class Gas:
    """Properties of the working gas used inside the compressor."""

    gamma: float = 1.4  # heat-capacity ratio (Cp/Cv)
    M: float = 0.02896  # molar mass of air in kg/mol

    @property
    def Rspec(self) -> float:
        """Specific gas constant R = universal gas constant / molar mass."""
        return 8.314462618 / self.M

    @property
    def cv(self) -> float:
        """Specific heat capacity at constant volume."""
        return self.Rspec / (self.gamma - 1)

    @property
    def cp(self) -> float:
        """Specific heat capacity at constant pressure."""
        return self.cv + self.Rspec


@dataclass
class Boundaries:
    """Fixed boundary conditions used by the compressor model."""

    Psuc: float = 101325.0  # suction pressure in pascals
    Tsuc: float = 300.0  # suction temperature in kelvin
    Pdis: float = 175 * 6894.75729  # discharge pressure in pascals


@dataclass
class Kinematics:
    """Piston and crank geometry for the cylinder volume calculation."""

    RPM: float = 1902.0  # crankshaft speed in revolutions per minute
    a: float = 0.0379984  # crank radius in meters
    l: float = 0.040005  # connecting rod length in meters
    bore: float = 0.08001  # cylinder bore diameter in meters
    Vc: float = 0.0000100888  # clearance volume in cubic meters

    @property
    def omega(self) -> float:
        """Angular velocity in radians per second."""
        return 2 * np.pi * self.RPM / 60

    @property
    def Area(self) -> float:
        """Cylinder cross-sectional area in square meters."""
        return 0.25 * np.pi * self.bore**2


@dataclass
class Reed:
    """Dynamic and flow parameters for a reed valve."""

    M: float  # effective moving mass of the reed in kilograms
    k: float  # spring stiffness in newtons per meter
    c: float  # damping coefficient in newton-seconds per meter
    Ap: float  # maximum physical port area in square meters
    d: float  # effective opening diameter parameter in meters
    CD: float  # discharge coefficient for the orifice flow
    ymax: float  # maximum reed lift in meters
    tlock: float  # time before the reed is allowed to move after start
    y0: float = 1e-7  # initial reed lift position in meters
    v0: float = 0.0  # initial reed velocity in meters per second
    Pcrack: float = 0.0  # minimum pressure differential to open the reed


@dataclass
class Numerics:
    """Numerical tolerances and solver settings for the ODE integration."""

    epsA: float = 1e-10  # minimum effective area to consider a valve open
    epsm: float = 1e-9  # minimum cylinder mass to avoid division by zero
    nCycles: int = 6  # number of crank cycles to simulate
    rel_tol: float = 1e-7  # relative tolerance for the ODE solver
    abs_tol: float = 1e-9  # absolute tolerance for the ODE solver


def default_reeds() -> tuple[Reed, Reed]:
    """Create default suction and discharge reed valve parameters.

    Returns:
        tuple[Reed, Reed]: suction and discharge reed valve objects.
    """
    vM = 0.00049398  # base moving mass constant for reed valves

    s_holefactor = 6 * 2  # factor scaling the suction port geometry
    suction = Reed(
        M=vM * s_holefactor,
        k=9534,
        c=1.2725,
        Ap=0.000127235 * s_holefactor,
        d=0.009 * s_holefactor,
        CD=0.930,
        ymax=0.002,
        tlock=0.003,
    )

    d_holefactor = 3 * 2  # factor scaling the discharge port geometry
    discharge = Reed(
        M=vM * d_holefactor,
        k=10012,
        c=1.128,
        Ap=0.000127235 * d_holefactor,
        d=0.009 * d_holefactor,
        CD=0.953,
        ymax=0.002,
        tlock=suction.tlock,
        Pcrack=0.0,
    )
    return suction, discharge


def geom(theta: float, kin: Kinematics) -> tuple[float, float]:
    """Compute instantaneous cylinder volume and volume derivative.

    Args:
        theta: crank angle in radians.
        kin: kinematic parameters describing crank and connecting rod.

    Returns:
        tuple[float, float]: current cylinder volume V and its time derivative dV/dt.
    """
    a, L, A = kin.a, kin.l, kin.Area
    root = max(L**2 - (a * np.sin(theta)) ** 2, 0.0)
    s = a * np.cos(theta) + np.sqrt(root)
    x = (a + L) - s
    V = kin.Vc + A * x

    denom = max(np.sqrt(root), 1e-12)
    dx_dtheta = a * np.sin(theta) + (a**2 * np.sin(theta) * np.cos(theta)) / denom
    dVdt = A * dx_dtheta * kin.omega
    return V, dVdt


def orifice_mdot(pu: float, Tu: float, pd: float, Aeff: float, CD: float, gas: Gas) -> float:
    """Calculate mass flow through an orifice or reed valve.

    This function supports both subsonic and choked flow regimes.
    """
    if Aeff <= 0 or pu <= pd:
        return 0.0
    g = gas.gamma
    R = gas.Rspec
    pi_crit = (2 / (g + 1)) ** (g / (g - 1))
    Pi = pd / pu

    if Pi > pi_crit:
        rho_u = pu / (R * Tu)
        return CD * Aeff * np.sqrt(2 * rho_u * (pu - pd))

    return CD * Aeff * pu * np.sqrt(g / (R * Tu)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))


def reed_accel(p_up: float, p_dn: float, Aeff: float, reed: Reed, y: float, v: float) -> float:
    """Compute reed valve acceleration from pressure, spring, and damping forces."""
    F_p = max(p_up - p_dn, 0) * Aeff
    F_spr = -reed.k * y
    F_dmp = -reed.c * v

    F_seat = 0.0
    if y < 0:
        F_seat = -50 * reed.k * y - 10 * reed.c * v

    F_lim = 0.0
    if y > reed.ymax:
        F_lim = -50 * reed.k * (y - reed.ymax) - 10 * reed.c * v

    return (F_p + F_spr + F_dmp + F_seat + F_lim) / reed.M


def simulate(config: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """Run the compressor simulation for the requested configuration.

    Args:
        config: optional parameter dictionary. Supported keys include:
            RPM, Psuc, Tsuc, Pdis, kinL, kinA, kinBore, kinVc,
            suctionHolefactor, suctionAp, suctionD, suctionCD,
            dischargeHolefactor, dischargeAp, dischargeD, dischargeCD,
            nCycles.

    Returns:
        A dictionary containing time series results and final averaged metrics.
    """
    gas = Gas()
    bc = Boundaries()
    kin = Kinematics()
    s, d = default_reeds()
    num = Numerics()

    if config:
        kin.RPM = float(config.get("RPM", kin.RPM))
        kin.a = float(config.get("kinA", kin.a))
        kin.l = float(config.get("kinL", kin.l))
        kin.bore = float(config.get("kinBore", kin.bore))
        kin.Vc = float(config.get("kinVc", kin.Vc))

        # Convert boundary condition pressures from psi to pascals if needed.
        bc.Pdis = float(config.get("Pdis", bc.Pdis * PA_TO_PSI)) * PSI_TO_PA
        bc.Psuc = float(config.get("Psuc", bc.Psuc * PA_TO_PSI)) * PSI_TO_PA
        bc.Tsuc = float(config.get("Tsuc", bc.Tsuc))

        s_holefactor = float(config.get("suctionHolefactor", 12))
        d_holefactor = float(config.get("dischargeHolefactor", 6))
        s.M = 0.00049398 * s_holefactor
        d.M = 0.00049398 * d_holefactor
        s.Ap = float(config.get("suctionAp", 0.000127235)) * s_holefactor
        d.Ap = float(config.get("dischargeAp", 0.000127235)) * d_holefactor
        s.d = float(config.get("suctionD", 0.009)) * s_holefactor
        d.d = float(config.get("dischargeD", 0.009)) * d_holefactor
        s.CD = float(config.get("suctionCD", s.CD))
        d.CD = float(config.get("dischargeCD", d.CD))
        num.nCycles = int(config.get("nCycles", num.nCycles))

    period = 2 * np.pi / kin.omega
    t_end = num.nCycles * period
    t_eval = np.linspace(0, t_end, num.nCycles * 1200)

    V0, _ = geom(0.0, kin)
    m0 = max(bc.Psuc * V0 / (gas.Rspec * bc.Tsuc), 1e-6)
    x0 = np.array([m0, bc.Tsuc, s.y0, s.v0, d.y0, d.v0], dtype=float)

    def rhs(t: float, x: np.ndarray) -> np.ndarray:
        """Evaluate the ODE right-hand side for the compressor state vector."""
        m = max(x[0], num.epsm)
        T = x[1]
        ys, vs, yd, vd = x[2], x[3], x[4], x[5]

        theta = np.mod(kin.omega * t, 2 * np.pi)
        V, dVdt = geom(theta, kin)
        p_cyl = m * gas.Rspec * T / V

        Aeff_s = min(s.Ap, np.pi * s.d * max(ys, 0))
        open_s = (bc.Psuc > p_cyl) and (t > s.tlock) and (Aeff_s > num.epsA)
        mdot_s = orifice_mdot(bc.Psuc, bc.Tsuc, p_cyl, Aeff_s, s.CD, gas) if open_s else 0.0
        as_ = reed_accel(bc.Psuc, p_cyl, Aeff_s, s, ys, vs)

        Aeff_d = min(d.Ap, np.pi * d.d * max(yd, 0))
        open_d = (p_cyl > (bc.Pdis + d.Pcrack)) and (t > d.tlock) and (Aeff_d > num.epsA)
        mdot_d = orifice_mdot(p_cyl, T, bc.Pdis, Aeff_d, d.CD, gas) if open_d else 0.0
        ad_ = reed_accel(p_cyl, bc.Pdis, Aeff_d, d, yd, vd)

        dUdt = mdot_s * gas.cp * bc.Tsuc - mdot_d * gas.cp * T - p_cyl * dVdt
        mdot_cyl = mdot_s - mdot_d
        Tdot = (dUdt - gas.cv * T * mdot_cyl) / (m * gas.cv)

        return np.array([mdot_cyl, Tdot, vs, as_, vd, ad_], dtype=float)

    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        x0,
        method="BDF",
        t_eval=t_eval,
        rtol=num.rel_tol,
        atol=num.abs_tol,
        max_step=period / 2000,
    )

    t = sol.t
    m = sol.y[0]
    T = sol.y[1]
    ys = sol.y[2]
    yd = sol.y[4]

    theta = kin.omega * t
    V = np.array([geom(th, kin)[0] for th in theta])
    p_cyl = m * gas.Rspec * T / V

    # Compute the discharge volumetric flow rate from the discharge mass flow history.
    vdot_d = np.zeros_like(t)
    for i in range(len(t)):
        Aeff_d = min(d.Ap, np.pi * d.d * max(yd[i], 0))
        open_d = p_cyl[i] > (bc.Pdis + d.Pcrack)
        if open_d:
            mdot_d = orifice_mdot(p_cyl[i], T[i], bc.Pdis, Aeff_d, d.CD, gas)
            vdot_d[i] = mdot_d * gas.Rspec * T[i] / p_cyl[i]

    total_flow = np.trapezoid(vdot_d, t)
    average_flow_rate = float(total_flow / t[-1])
    time_estimate = float(0.22712461 / average_flow_rate) if average_flow_rate > 0 else float("inf")

    return {
        "time": t.tolist(),
        "pressure": (p_cyl * PA_TO_PSI).tolist(),
        "temperature": T.tolist(),
        "suction_lift": ys.tolist(),
        "discharge_lift": yd.tolist(),
        "average_flow_rate_m3s": average_flow_rate,
        "fill_time_seconds": time_estimate,
        "config": {
            "gas": asdict(gas),
            "boundaries": asdict(bc),
            "kinematics": asdict(kin),
            "suction": asdict(s),
            "discharge": asdict(d),
            "numerics": asdict(num),
        },
    }
