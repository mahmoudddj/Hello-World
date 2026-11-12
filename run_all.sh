#!/usr/bin/env bash
# Rerun everything from a clean environment.
#   ./run_all.sh              tests + synthetic validation only (no patient data needed)
#   GENIE_BPC_DATA_DIR=/secure/path ./run_all.sh     also runs G1 feasibility and, if the gate passes, the analysis
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv
.venv/bin/pip -q install -r requirements.lock
export PYTHONPATH=src
.venv/bin/python -m pytest -q tests
.venv/bin/python -m geniebpc.run --synthetic --out outputs/synthetic_validation --boot 200
if [ -n "${GENIE_BPC_DATA_DIR:-}" ]; then
  .venv/bin/python -m geniebpc.run --data-dir "$GENIE_BPC_DATA_DIR" --out outputs --boot 1000 || echo "Gate failed or analysis stopped: see outputs/feasibility_auto.md"
else
  echo "GENIE_BPC_DATA_DIR not set: real-data feasibility/analysis skipped."
fi
