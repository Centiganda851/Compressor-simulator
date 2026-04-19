"""Web API entrypoint for the Shefa compressor simulation.

This Flask application serves the static frontend and exposes a single
simulation endpoint at `/api/simulate`.
"""

from flask import Flask, jsonify, request, send_from_directory
from compressor_sim import simulate

# Flask application configured to serve static files from the repository root.
app = Flask(__name__, static_folder=".", static_url_path="")


@app.get("/")
def root():
    """Serve the main HTML page for the compressor UI."""
    return send_from_directory(".", "index.html")


@app.post("/api/simulate")
def run_simulation():
    """Run the compressor simulation using JSON payload parameters.

    The request payload is parsed as JSON, passed through to the compressor
    model, and the returned result is serialized back to JSON.
    """
    payload = request.get_json(silent=True) or {}
    result = simulate(payload)
    return jsonify(result)


if __name__ == "__main__":
    # Start the Flask development server on port 8000.
    app.run(host="0.0.0.0", port=8000, debug=True)
