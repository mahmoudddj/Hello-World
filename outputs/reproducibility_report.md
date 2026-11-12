# Reproducibility report (2026-11-12)

| Item | Status |
|---|---|
| Clean-environment rerun | `./run_all.sh` creates `.venv` from `requirements.lock`, runs tests and the synthetic validation |
| Tests | 29 passed (`pytest -q tests`) |
| Synthetic validation | `outputs/synthetic_validation/` — labelled synthetic in every file header/status; planted interaction recovered within bootstrap CIs |
| Real-file-layout smoke test | genieBPC simulated tables exported to the release file names: schema check passed, positivity gate failed as expected (6 eligible simulated patients); outputs not kept |
| Design precision | `outputs/diagnostics/design_precision.csv` (synthetic planning) |
| Real-data analysis | **Not run — no data access** |
| Seeds | Bootstrap 20261112; synthetic seeds fixed in code |
| Patient-level output | None written by design; only aggregate tables and figures |
