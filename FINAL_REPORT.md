# Report — INTERIM, NOT FROZEN (2026-11-12)

## Feasibility verdict
**Undetermined on real data; blocked at G0 (no data access).** Design-stage evidence suggests the primary question is unlikely to be answerable with useful precision in the public release (v2.0-public), and may only be answerable in the consortium release (v3.1).

## Gates
| Gate | Verdict |
|---|---|
| G0 data rights/files/dictionary | FAIL — no patient data reachable; dictionary only structurally verified |
| G1 variable feasibility | NOT RUN (generator ready) |
| G2 novelty/positivity/protocol | PROVISIONAL — protocol v0.1 written; novelty narrow (Elkrief 2023) |
| G3 SAP/fixtures/pipeline | PASS on synthetic data (29 tests) |
| G4 analyses | NOT STARTED |
| G5 audit/manuscript | NOT STARTED (introduction and methods drafts only) |

## Final estimand (pre-specified)
Between-stratum difference (STK11 and/or KEAP1 altered vs wild-type) in the 12-month RMST difference in overall survival between first-line anti-PD-(L)1 + platinum chemotherapy and anti-PD-(L)1 monotherapy.

## Sample size and numerical findings
None. No patient data were analysed.

## Main risks identified before data
1. Few first-line chemo-immunotherapy patients expected (cohort sequenced 2014–2017).
2. No performance status in the public release.
3. Overlapping prior work from the two largest contributing centres with a null genomic-interaction finding.

## Rerun
```
./run_all.sh                                        # tests + synthetic validation
GENIE_BPC_DATA_DIR=/secure/path ./run_all.sh        # adds G1 feasibility and, if the gate passes, the analysis
```
This report will be frozen with a commit hash only after G4/G5 on real data.
