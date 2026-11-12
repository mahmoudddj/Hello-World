# Statistical analysis plan v0.1 (fixed before data access, 2026-11-12)

1. **Population**: as in `target_trial.md`. Patient flow reported at every step.
2. **Primary estimand**: I(12) = Δ(1,12) − Δ(0,12), difference between biomarker strata in the 12-month RMST difference (chemo-IO minus monotherapy) for overall survival. Key secondary: I(24), stratum-specific Δ(m,τ), survival-probability differences at 12/24 months.
3. **Primary biomarker contrast**: STK11 and/or KEAP1 altered vs both wild-type. Secondary: STK11 alone, KEAP1 alone (each vs wild-type for that gene), loss-of-function-only definition. Exploratory: within KRAS-mutant tumours.
4. **Primary estimator**: stabilised inverse-probability-of-treatment weights from a logistic propensity model (covariates in `derive.design_matrix` plus biomarker), stabilised within biomarker stratum, truncated at the 1st/99th percentiles; weighted Kaplan–Meier → RMST. 95% CIs by patient-level bootstrap (1,000 replicates, nuisance models refitted). The interaction is tested with a Wald test using the bootstrap SE, two-sided α = 0.05, on I(12) only.
5. **Doubly robust analysis** (run if every cell has ≥40 patients): AIPW on RMST pseudo-observations computed within arm × biomarker.
6. **Secondary multiplicative scale**: weighted Cox model with treatment × biomarker product term, robust variance; reported only with a proportional-hazards check, never as the headline.
7. **Diagnostics**: standardised mean differences overall and within each stratum (target <0.1), propensity overlap, effective sample size, maximum weight.
8. **Missing data**: PD-L1, smoking and race use an explicit unknown category (primary). Planned sensitivity: multiple imputation by chained equations if any covariate is >10% missing (not yet implemented).
9. **Stop rules**: (a) any required table/column missing; (b) any treatment × biomarker cell with <20 patients or <10 deaths → no causal analysis, descriptive report only; (c) 24-month results suppressed if <10 at risk beyond 24 months in any cell.
10. **Sensitivity analyses**: any-time NGS with delayed entry; stratum-specific propensity models; loss-of-function-only biomarker; unadjusted; restriction to non-squamous; restriction to treatment start in the chemo-IO era; PD-L1 ≥50% and <50% separately (descriptive if cells are small); E-value for the stratum-specific effects; planned quantitative bias analysis for unmeasured performance status.
11. **Multiplicity**: one confirmatory test (I(12)). Everything else is estimation with 95% CIs, labelled secondary or exploratory.
12. **Interpretation rules**: no subgroup-significance comparisons; prognostic associations reported separately from the interaction; an inconclusive interaction (CI including both no interaction and clinically relevant values, pre-set at ±1.5 months) is reported as inconclusive.
13. **Reproducibility**: `run_all.sh`, `requirements.lock`, fixed seeds, aggregate outputs only.
Deviations after data access will be dated in `docs/decision_log.md`.
