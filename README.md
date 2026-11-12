# STK11/KEAP1 and first-line immunotherapy strategy in advanced NSCLC (GENIE BPC) — independent study

**Status (2026-11-12): blocked at gate G0 — no GENIE BPC patient data are available to this workspace.
No real-data result exists. Everything numerical in this repository comes from synthetic validation data.**

Question: do pretreatment STK11 and/or KEAP1 alterations modify the comparative effectiveness of first-line
anti-PD-(L)1 monotherapy versus anti-PD-(L)1 plus platinum chemotherapy on overall survival?

| Where | What |
|---|---|
| `FINAL_REPORT.md` | Current verdict, what is done, what is blocked, rerun commands (interim, not frozen) |
| `docs/gate_G0.md` … `gate_G3.md` | Gate decisions with evidence |
| `docs/feasibility_report.md` | Variable-by-variable feasibility, design-stage precision, alternative questions |
| `docs/literature_matrix.csv`, `docs/literature_notes.md` | Literature (verification level stated per row) |
| `protocol/` | Target trial, statistical analysis plan, DAG (v0.1, pre-data) |
| `src/geniebpc/`, `tests/` | Pipeline and 29 automated checks on synthetic fixtures |
| `outputs/synthetic_validation/` | Pipeline output on synthetic data (validation only) |

Run: `./run_all.sh` (tests + synthetic validation); with data: `GENIE_BPC_DATA_DIR=/secure/path ./run_all.sh`.

Patient-level data, credentials and tokens are never committed (`.gitignore`); the pipeline writes aggregate outputs only.
