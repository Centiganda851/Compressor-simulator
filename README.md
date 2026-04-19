# Reciprocating Compressor Simulator (MATLAB-to-Python)

This project converts the provided MATLAB compressor model into an open-source Python library and a browser UI.

## What this includes

- **Python simulation engine** (`compressor_sim.py`) that ports the MATLAB equations.
- **Web backend** (`app.py`) with one API endpoint: `POST /api/simulate`.
- **Frontend UI** (`index.html`, `script.js`, `style.css`) to:
  - Edit simulation parameters.
  - Choose graphs with checkboxes.
  - Save, load, and delete named presets in browser local storage.
  - Show average flow rate and estimated fill time.
  - Run in two modes:
    - **Download-and-run mode**: open `index.html` directly and run simulation in-browser (no Python/MATLAB required).
    - **Server mode**: run Flask backend and use SciPy-based simulation endpoint.

## Requirements

### For download-and-run mode

- Any modern browser.
- No local Python installation required.

### For server mode

- Python 3.10+
- Packages in `requirements.txt`

## Setup (server mode)


```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run (server mode)

```bash
python app.py
```

Open: `http://localhost:8000`

## Download-and-run mode

Just share these files together in one folder and open `index.html`:

- `index.html`
- `script.js`
- `style.css`

Then click **Run Simulation**.

## API example (server mode)

```bash
curl -X POST http://localhost:8000/api/simulate \
  -H "Content-Type: application/json" \
  -d '{"RPM":1902,"nCycles":6}'
```

## Notes for non-MATLAB users

- No MATLAB runtime is required.
- Server mode uses SciPy's `solve_ivp(..., method="BDF")`, suitable for stiff systems similar to MATLAB `ode15s`.
- Download-and-run mode uses an in-browser numerical solver for portability.
- Units are SI (Pa, K, m, kg, s).
