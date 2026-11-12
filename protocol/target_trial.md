# Target trial protocol v0.1 (written before data access, 2026-11-12)

| Component | Target trial | Emulation in GENIE BPC NSCLC |
|---|---|---|
| Eligibility | Adults with stage IV or recurrent metastatic NSCLC, no prior systemic therapy for advanced disease, no sensitising EGFR mutation or ALK/ROS1 fusion, tumour NGS result covering STK11 and KEAP1 available before treatment | Index NSCLC; stage IV at diagnosis or later distant metastasis; first regimen on/after that date; covering-panel NGS report dated ≤ regimen start; driver exclusions from mutation/fusion files |
| Strategies | (0) anti-PD-(L)1 monotherapy; (1) anti-PD-(L)1 + platinum-based chemotherapy (± bevacizumab) | First-line `regimen_drugs` classified by `classify_regimen`; CTLA-4 combinations, non-platinum combinations, trial and investigational regimens excluded |
| Assignment | Randomised, stratified by biomarker | Assumed conditionally exchangeable given measured covariates within biomarker stratum (not credible without performance status; see DAG) |
| Time zero | Randomisation = treatment start | First-line regimen start (`dx_reg_start_int`); eligibility, biomarker and covariates all assessed at or before it |
| Grace period | None | Agents grouped into one regimen by the release's curation rule; no additional grace period, so no immortal time is introduced by us |
| Outcome | Overall survival | `tt_os_g_days`, `os_g_status` |
| Follow-up | Until death, loss to follow-up or data cut-off; horizons 12 and 24 months | Same; 24-month horizon reported only if ≥10 patients per cell remain at risk beyond 24 months |
| Causal contrast | Intention-to-treat effect of starting strategy 1 vs 0 | Observational analogue (treatment initiated); later switching/cross-over is part of the strategy |
| Estimands | Δ(m,τ) = RMSTτ(strategy 1) − RMSTτ(strategy 0) within biomarker stratum m; interaction I(τ) = Δ(1,τ) − Δ(0,τ) | IPTW-standardised within stratum; see SAP |
| Effect modifier | STK11 and/or KEAP1 altered vs both wild-type (pretreatment) | As defined in `derive.py` |

Effect modification here means the treatment effect differs across biomarker strata; it is not a claim that intervening on the biomarker would change the effect. Survival differences by biomarker within one arm are prognostic and are reported separately, never as evidence of prediction.
