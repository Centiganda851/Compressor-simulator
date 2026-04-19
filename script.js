/**
 * Compressor simulation frontend logic.
 *
 * This script collects input values from the page, runs a compressor model
 * either in-browser or via the `/api/simulate` endpoint, and renders charts.
 */

// Unit conversion constants.
const PSI_TO_PA = 6894.75729; // pascals per psi
const PA_TO_PSI = 1 / PSI_TO_PA; // psi per pascal
const M3S_TO_FT3MIN = 2118.880003; // cubic meters per second to cubic feet per minute

// List of supported configuration keys from the UI.
const fields = [
  "RPM", "Psuc", "Tsuc", "Pdis", "kinL", "kinA", "kinBore", "kinVc",
  "suctionHolefactor", "suctionAp", "suctionD", "suctionCD",
  "dischargeHolefactor", "dischargeAp", "dischargeD", "dischargeCD", "nCycles",
];

// Active Chart.js instances keyed by series name.
const charts = {};

/**
 * Read current form input values and convert them to numbers.
 * @returns {Object<string, number>}
 */
function getConfig() {
  const cfg = {};
  fields.forEach((key) => {
    cfg[key] = Number(document.getElementById(key).value);
  });
  return cfg;
}

/**
 * Return the list of graph series the user wants displayed.
 * @returns {string[]}
 */
function getSelectedGraphs() {
  return [...document.querySelectorAll('.checks input[type="checkbox"]')]
    .filter((el) => el.checked)
    .map((el) => el.value);
}

/**
 * Convert seconds into a formatted minute/second string.
 * @param {number} sec
 * @returns {string}
 */
function formatDuration(sec) {
  if (!Number.isFinite(sec)) return "∞";
  const m = Math.floor(sec / 60);
  const s = sec % 60;
  return `${m} min ${s.toFixed(2)} s`;
}

/**
 * Build the default compressor model state from configuration values.
 * @param {Object} config
 * @returns {Object}
 */
function defaultModel(config) {
  const gas = {
    gamma: 1.4,
    M: 0.02896,
    get Rspec() { return 8.314462618 / this.M; },
    get cv() { return this.Rspec / (this.gamma - 1); },
    get cp() { return this.cv + this.Rspec; },
  };

  const bc = {
    Psuc: (config.Psuc ?? 14.6959) * PSI_TO_PA,
    Tsuc: config.Tsuc ?? 300,
    Pdis: (config.Pdis ?? 175) * PSI_TO_PA,
  };

  const kin = {
    RPM: config.RPM ?? 1902,
    a: config.kinA ?? 0.0379984,
    l: config.kinL ?? 0.040005,
    bore: config.kinBore ?? 0.08001,
    Vc: config.kinVc ?? 0.0000100888,
    get omega() { return 2 * Math.PI * this.RPM / 60; },
    get Area() { return 0.25 * Math.PI * this.bore ** 2; },
  };

  const valve_Mass = 0.00049398;
  const sHolefactor = config.suctionHolefactor ?? 12;
  const suction = {
    M: valve_Mass * sHolefactor,
    k: 9534,
    c: 1.2725,
    Ap: (config.suctionAp ?? 0.000127235) * sHolefactor,
    d: (config.suctionD ?? 0.009) * sHolefactor,
    CD: config.suctionCD ?? 0.93,
    ymax: 0.002,
    tlock: 0.003,
    y0: 1e-7,
    v0: 0,
  };

  const dHolefactor = config.dischargeHolefactor ?? 6;
  const discharge = {
    M: valve_Mass * dHolefactor,
    k: 10012,
    c: 1.128,
    Ap: (config.dischargeAp ?? 0.000127235) * dHolefactor,
    d: (config.dischargeD ?? 0.009) * dHolefactor,
    CD: config.dischargeCD ?? 0.953,
    ymax: 0.002,
    tlock: suction.tlock,
    y0: 1e-7,
    v0: 0,
    Pcrack: 0,
  };

  const numerics = {
    epsA: 1e-10,
    epsm: 1e-9,
    nCycles: Math.max(1, Math.floor(config.nCycles ?? 6)),
    samplesPerCycle: 1200,
  };

  return { gas, bc, kin, suction, discharge, numerics };
}

/**
 * Compute cylinder volume and volumetric rate from the crank angle.
 * @param {number} theta - crank angle in radians
 * @param {Object} kin - kinematic geometry parameters
 * @returns {{V: number, dVdt: number}}
 */
function geom(theta, kin) {
  const a = kin.a;
  const L = kin.l;
  const A = kin.Area;
  const root = Math.max(L ** 2 - (a * Math.sin(theta)) ** 2, 0);
  const s = a * Math.cos(theta) + Math.sqrt(root);
  const x = (a + L) - s;
  const V = kin.Vc + A * x;

  const denom = Math.max(Math.sqrt(root), 1e-12);
  const dx_dtheta = a * Math.sin(theta) + ((a ** 2 * Math.sin(theta) * Math.cos(theta)) / denom);
  const dVdt = A * dx_dtheta * kin.omega;
  return { V, dVdt };
}

/**
 * Estimate the mass flow through an orifice for the current pressure conditions.
 * @param {number} pu - upstream pressure
 * @param {number} Tu - upstream temperature
 * @param {number} pd - downstream pressure
 * @param {number} Aeff - effective orifice area
 * @param {number} CD - discharge coefficient
 * @param {Object} gas - gas thermodynamic properties
 * @returns {number}
 */
function orificeMdot(pu, Tu, pd, Aeff, CD, gas) {
  if (Aeff <= 0 || pu <= pd) return 0;
  const g = gas.gamma;
  const R = gas.Rspec;
  const piCrit = (2 / (g + 1)) ** (g / (g - 1));
  const Pi = pd / pu;

  if (Pi > piCrit) {
    const rhoU = pu / (R * Tu);
    return CD * Aeff * Math.sqrt(2 * rhoU * (pu - pd));
  }

  return CD * Aeff * pu * Math.sqrt(g / (R * Tu)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)));
}

/**
 * Compute reed acceleration from net force contributions.
 * @param {number} pUp - upstream pressure
 * @param {number} pDn - downstream pressure
 * @param {number} Aeff - effective valve area
 * @param {Object} reed - reed valve parameters
 * @param {number} y - current reed deflection
 * @param {number} v - current reed velocity
 * @returns {number}
 */
function reedAccel(pUp, pDn, Aeff, reed, y, v) {
  const Fp = Math.max(pUp - pDn, 0) * Aeff;
  const Fspr = -reed.k * y;
  const Fdmp = -reed.c * v;

  let Fseat = 0;
  if (y < 0) Fseat = -50 * reed.k * y - 10 * reed.c * v;

  let Flim = 0;
  if (y > reed.ymax) Flim = -50 * reed.k * (y - reed.ymax) - 10 * reed.c * v;

  return (Fp + Fspr + Fdmp + Fseat + Flim) / reed.M;
}

/**
 * Right-hand side of the ODE system describing the cylinder state.
 * @param {number} t - time in seconds
 * @param {number[]} x - state vector [mass, temperature, suction lift, suction velocity, discharge lift, discharge velocity]
 * @param {Object} model - current compressor model parameters
 * @returns {number[]}
 */
function rhs(t, x, model) {
  const { gas, bc, kin, suction: s, discharge: d, numerics: n } = model;

  const m = Math.max(x[0], n.epsm);
  const T = x[1];
  const ys = x[2];
  const vs = x[3];
  const yd = x[4];
  const vd = x[5];

  const theta = ((kin.omega * t) % (2 * Math.PI));
  const { V, dVdt } = geom(theta, kin);
  const pCyl = m * gas.Rspec * T / V;

  const AeffS = Math.min(s.Ap, Math.PI * s.d * Math.max(ys, 0));
  const openS = (bc.Psuc > pCyl) && (t > s.tlock) && (AeffS > n.epsA);
  const mdotS = openS ? orificeMdot(bc.Psuc, bc.Tsuc, pCyl, AeffS, s.CD, gas) : 0;
  const as = reedAccel(bc.Psuc, pCyl, AeffS, s, ys, vs);

  const AeffD = Math.min(d.Ap, Math.PI * d.d * Math.max(yd, 0));
  const openD = (pCyl > (bc.Pdis + d.Pcrack)) && (t > d.tlock) && (AeffD > n.epsA);
  const mdotD = openD ? orificeMdot(pCyl, T, bc.Pdis, AeffD, d.CD, gas) : 0;
  const ad = reedAccel(pCyl, bc.Pdis, AeffD, d, yd, vd);

  const dUdt = mdotS * gas.cp * bc.Tsuc - mdotD * gas.cp * T - pCyl * dVdt;
  const mdotCyl = mdotS - mdotD;
  const Tdot = (dUdt - gas.cv * T * mdotCyl) / (m * gas.cv);

  return [mdotCyl, Tdot, vs, as, vd, ad];
}

/**
 * Single fourth-order Runge-Kutta integration step.
 * @param {number} t - current time
 * @param {number[]} x - current state vector
 * @param {number} dt - timestep
 * @param {Object} model - compressor model parameters
 * @returns {number[]}
 */
function rk4Step(t, x, dt, model) {
  const k1 = rhs(t, x, model);
  const x2 = x.map((xi, i) => xi + 0.5 * dt * k1[i]);
  const k2 = rhs(t + 0.5 * dt, x2, model);
  const x3 = x.map((xi, i) => xi + 0.5 * dt * k2[i]);
  const k3 = rhs(t + 0.5 * dt, x3, model);
  const x4 = x.map((xi, i) => xi + dt * k3[i]);
  const k4 = rhs(t + dt, x4, model);
  return x.map((xi, i) => xi + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

/**
 * Run the compressor simulation inside the browser using explicit RK4 integration.
 * @param {Object} config
 * @returns {Object}
 */
function simulateInBrowser(config) {
  const model = defaultModel(config);
  const { gas, bc, kin, suction: s, discharge: d, numerics: n } = model;

  const period = 2 * Math.PI / kin.omega;
  const tEnd = n.nCycles * period;
  const count = n.nCycles * n.samplesPerCycle;
  const dt = tEnd / (count - 1);
  const rk4Substeps = 8;
  const dtInternal = dt / rk4Substeps;

  const V0 = geom(0, kin).V;
  const m0 = Math.max((bc.Psuc * V0) / (gas.Rspec * bc.Tsuc), 1e-6);

  let x = [m0, bc.Tsuc, s.y0, s.v0, d.y0, d.v0];
  const time = new Array(count);
  const pressure = new Array(count);
  const temperature = new Array(count);
  const suctionLift = new Array(count);
  const dischargeLift = new Array(count);
  const vdotD = new Array(count);

  let t = 0;
  for (let i = 0; i < count; i += 1) {
    const theta = kin.omega * t;
    const { V } = geom(theta, kin);
    const pCyl = x[0] * gas.Rspec * x[1] / V;

    time[i] = t;
    pressure[i] = pCyl * PA_TO_PSI;
    temperature[i] = x[1];
    suctionLift[i] = x[2];
    dischargeLift[i] = x[4];

    const AeffD = Math.min(d.Ap, Math.PI * d.d * Math.max(x[4], 0));
    if (pCyl > (bc.Pdis + d.Pcrack)) {
      const mdotD = orificeMdot(pCyl, x[1], bc.Pdis, AeffD, d.CD, gas);
      vdotD[i] = mdotD * gas.Rspec * x[1] / pCyl;
    } else {
      vdotD[i] = 0;
    }

    if (i < count - 1) {
      for (let j = 0; j < rk4Substeps; j += 1) {
        x = rk4Step(t, x, dtInternal, model);
        t += dtInternal;
      }
    }
  }

  let totalFlow = 0;
  for (let i = 1; i < count; i += 1) {
    totalFlow += 0.5 * (vdotD[i - 1] + vdotD[i]) * (time[i] - time[i - 1]);
  }

  const averageFlow = totalFlow / time[count - 1];
  const fillTime = averageFlow > 0 ? 0.22712461 / averageFlow : Infinity;

  return {
    time,
    pressure,
    temperature,
    suction_lift: suctionLift,
    discharge_lift: dischargeLift,
    average_flow_rate_m3s: averageFlow,
    fill_time_seconds: fillTime,
  };
}

/**
 * Call the backend API to perform the simulation.
 * @param {Object} config
 * @returns {Promise<Object>}
 */
async function simulateWithApi(config) {
  const response = await fetch('/api/simulate', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(config),
  });

  if (!response.ok) {
    throw new Error(`API failed with status ${response.status}`);
  }
  return response.json();
}

/**
 * Main simulation runner triggered by the UI.
 * Chooses API mode when available and falls back to browser mode.
 */
async function runSimulation() {
  const config = getConfig();
  let result;

  if (window.location.protocol === 'file:') {
    result = simulateInBrowser(config);
  } else {
    try {
      result = await simulateWithApi(config);
    } catch (_err) {
      result = simulateInBrowser(config);
    }
  }

  const avgFt3Min = result.average_flow_rate_m3s * M3S_TO_FT3MIN;
  document.getElementById('avgFlow').textContent = avgFt3Min.toFixed(4);
  document.getElementById('fillTime').textContent = formatDuration(result.fill_time_seconds);
  renderCharts(result);
}

/**
 * Render the selected charts using Chart.js.
 * @param {Object} result
 */
function renderCharts(result) {
  if (typeof Chart === 'undefined') {
    alert('Chart.js is unavailable. Results are still computed and shown as text.');
    return;
  }

  const selected = getSelectedGraphs();
  const container = document.getElementById('charts');

  Object.values(charts).forEach((instance) => instance.destroy());
  Object.keys(charts).forEach((key) => delete charts[key]);
  container.innerHTML = '';

  const colors = {
    pressure: '#ef4444',
    temperature: '#3b82f6',
    suction_lift: '#10b981',
    discharge_lift: '#a855f7',
  };

  const displayName = {
    pressure: 'Pressure (psi)',
    temperature: 'Temperature (K)',
    suction_lift: 'Suction Lift (m)',
    discharge_lift: 'Discharge Lift (m)',
  };

  selected.forEach((key) => {
    const figure = document.createElement('figure');
    figure.className = 'chart-figure';

    const title = document.createElement('figcaption');
    title.textContent = displayName[key];

    const canvas = document.createElement('canvas');
    figure.appendChild(title);
    figure.appendChild(canvas);
    container.appendChild(figure);

    charts[key] = new Chart(canvas.getContext('2d'), {
      type: 'line',
      data: {
        datasets: [{
          label: displayName[key],
          data: result.time.map((t, i) => ({ x: t, y: result[key][i] })),
          borderColor: colors[key],
          fill: false,
          pointRadius: 0,
          tension: 0.1,
        }],
      },
      options: {
        responsive: true,
        interaction: { mode: 'index', intersect: false },
        scales: {
          x: {
            type: 'linear',
            title: { display: true, text: 'Time (s)' },
            ticks: {
              maxTicksLimit: 8,
              callback: (value) => Number(value).toFixed(4),
            },
          },
          y: {
            title: { display: true, text: displayName[key] },
          },
        },
      },
    });
  });
}

/**
 * Load saved preset names from localStorage and populate the preset selector.
 */
function loadPresets() {
  const presets = JSON.parse(localStorage.getItem('compressor_presets') || '{}');
  const select = document.getElementById('presetSelect');
  select.innerHTML = '<option value="">Select preset</option>';
  Object.keys(presets).forEach((name) => {
    const opt = document.createElement('option');
    opt.value = name;
    opt.textContent = name;
    select.appendChild(opt);
  });
}

/**
 * Save the current configuration and selected graphs as a named preset.
 */
function savePreset() {
  const name = document.getElementById('presetName').value.trim();
  if (!name) return alert('Provide a preset name.');

  const presets = JSON.parse(localStorage.getItem('compressor_presets') || '{}');
  presets[name] = {
    config: getConfig(),
    graphs: getSelectedGraphs(),
  };
  localStorage.setItem('compressor_presets', JSON.stringify(presets));
  loadPresets();
}

/**
 * Load the selected preset values into the form and graph selection.
 */
function loadPreset() {
  const selected = document.getElementById('presetSelect').value;
  if (!selected) return;

  const presets = JSON.parse(localStorage.getItem('compressor_presets') || '{}');
  const preset = presets[selected];
  if (!preset) return;

  fields.forEach((key) => {
    if (preset.config[key] !== undefined) {
      document.getElementById(key).value = preset.config[key];
    }
  });

  document.querySelectorAll('.checks input[type="checkbox"]').forEach((cb) => {
    cb.checked = preset.graphs.includes(cb.value);
  });
}

/**
 * Delete the currently selected preset from localStorage.
 */
function deletePreset() {
  const selected = document.getElementById('presetSelect').value;
  if (!selected) return;
  const presets = JSON.parse(localStorage.getItem('compressor_presets') || '{}');
  delete presets[selected];
  localStorage.setItem('compressor_presets', JSON.stringify(presets));
  loadPresets();
}

document.getElementById('runBtn').addEventListener('click', runSimulation);
document.getElementById('savePresetBtn').addEventListener('click', savePreset);
document.getElementById('loadPresetBtn').addEventListener('click', loadPreset);
document.getElementById('deletePresetBtn').addEventListener('click', deletePreset);

document.addEventListener('DOMContentLoaded', loadPresets);
