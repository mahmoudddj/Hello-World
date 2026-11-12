# Decision log
| # | Date | Decision | Reason |
|---|---|---|---|
| 1 | 2026-11-12 | Python stack (pandas, lifelines, scikit-learn), single package `geniebpc` | R not installed; one language keeps the pipeline testable end to end |
| 2 | 2026-11-12 | Schema taken from genieBPC simulated data, flagged unverified | Only first-hand source of column names reachable; real data guide not readable |
| 3 | 2026-11-12 | Primary biomarker = STK11 and/or KEAP1 altered (pooled); each gene alone secondary | Frequent co-occurrence and shared biology; separate genes halve already small cells |
| 4 | 2026-11-12 | Primary cohort restricted to NGS reported on/before time zero | Biomarker must be a baseline characteristic; avoids immortal time. Any-time NGS with delayed entry is a sensitivity analysis |
| 5 | 2026-11-12 | Estimand on RMST difference scale at 12 and 24 months; hazard-ratio interaction secondary | Brief; additive scale is the clinically relevant one for treatment choice and does not rely on proportional hazards |
| 6 | 2026-11-12 | Hotspot panels (UHN-48/50, VICC-01-SOLIDTUMOR) treated as non-evaluable | Wild-type cannot be asserted without gene coverage; to be confirmed against panel gene lists |
| 7 | 2026-11-12 | Trial / investigational regimens excluded | Drug identity masked in the release |
| 8 | 2026-11-12 | Stop rule: ≥20 patients and ≥10 deaths in every treatment × biomarker cell | Below this, weights and bootstrap are unstable and the interaction is uninformative |
| 9 | 2026-11-12 | Missing PD-L1 handled by an explicit "unknown" category in the primary analysis | Expected high, era-dependent missingness; MI planned as sensitivity |
| 10 | 2026-11-12 | No pivot to alternative questions without researcher approval | Required by the brief |
| 11 | 2026-11-12 | Second agent's handoff archive opened; it contains only the shared brief | Supplied by the researcher; no cross-agent work product seen |
| 12 | 2026-11-12 | Isolation ended: researcher supplied the Codex R pipeline and asked for a review | Explicit researcher authorisation; findings in `codex_pipeline_review.md` |
| 13 | 2026-11-12 | Patient-level data stay on the researcher's PC; analyses delivered as R pipelines producing aggregate output | Researcher's decision on data handling |
| 14 | 2026-11-12 | Alternatives A, B, C implemented in R (`r_pipelines/`) ahead of formal approval of a pivot; feasibility counts decide which proceed | Researcher asked for pipelines for all three |
| 15 | 2026-11-12 | For A–C the primary NGS handling is delayed entry at report date; pretreatment-only is sensitivity | Requiring pretreatment NGS would discard most first-line chemotherapy patients in a 2014–2017 sequencing cohort |
