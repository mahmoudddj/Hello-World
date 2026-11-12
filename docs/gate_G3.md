# Gate G3 — analysis plan, fixtures, validated pipeline
**Verdict: PASS for the software, on synthetic data only.** (Does not imply the study is feasible.)

- SAP v0.1 fixed before data access (`protocol/SAP.md`).
- 29 automated checks pass (`pytest -q tests`): regimen classification (12 cases), KM/RMST agreement with lifelines (unweighted, weighted, delayed entry), chronology/leakage (post-baseline NGS and PD-L1 never used in the primary cohort), exclusions, planted-truth recovery for IPTW and AIPW, no manufactured interaction under the null, prognostic-only biomarker gives product-term HR ≈ 1, delayed entry removes immortal-time bias, gate stop rule, schema failure reporting, bootstrap determinism.
- Large-sample check (n=60,000 synthetic): IPTW interaction 1.39 / 3.28 months vs truth 1.34 / 3.10 at 12 / 24 months.
- Known weakness: with a main-effects propensity model, balance inside the smaller (altered) stratum can stay above |SMD| 0.1; a stratum-specific propensity model is the pre-specified sensitivity analysis.
- Not yet implemented: multiple imputation for PD-L1, inverse-probability-of-censoring weights, quantitative bias analysis for unmeasured performance status (listed in SAP as planned).
G4 and G5 not started (no data).
