# Review of the Codex R pipeline (GENIE_BPC_R_Pipeline v1.0.0)

Reviewed 2026-11-12 at the researcher's request (this ends the isolation between the two agents; logged in `decision_log.md`).
Everything below was checked by running the code in R 4.3.3 / survival 3.5.8 on artificial data. No patient data were used.

## What was checked and how
| Check | Result |
|---|---|
| Archive integrity (`SHA256SUMS.txt`) | All files match |
| Its own test suite (`tests/run_tests.R`) | Passes: 32 unit checks + end-to-end synthetic run, as claimed |
| Estimator correctness, independent cross-check | On one synthetic cohort (n ≈ 9,400) its IPTW RMST interaction is 1.116 / 2.784 months at 12 / 24 months; my separately written Python estimator gives 1.115 / 2.779; planted truth 1.29 / 3.07. The two implementations agree |
| Cohort logic, independent cross-check | Fed through the bridge below, its `make_cohort` reproduces exactly the four treatment × biomarker cells my own derivation gives on the same synthetic release-layout files (1285 / 1002 / 390 / 303) |
| RMST integration, weighted KM, bootstrap with refitting, Wald test | Correct |
| Export controls (allow-list, hash check, suppression) | Work as described; more thorough than mine |

## Problems found
| # | Severity | Finding | Fix |
|---|---|---|---|
| 1 | **Bug** | `secondary_cox()` uses `strata(biomarker)` without attaching `survival`. With survival 3.5.8 it fails with "could not find function strata"; the error is caught, so every run silently reports "Secondary Cox: NOT_ESTIMABLE" — including its own passing test run. Not tested on newer survival versions | Add `library(survival)` at the top of `R/common.R` (confirmed: the Cox table and PH test then appear) |
| 2 | **Blocks real data** | The release adapter is empty (`source_spec_template.R`); nothing runs on real files until ~50 columns are mapped by hand | Use `r_pipelines/codex_bridge/build_codex_canonical.R`, which builds its five input tables from the v2.0-public files |
| 3 | **Blocks real data** | `ecog` is a required complete-case covariate, but the public release has no performance-status variable → 100% complete-case loss → `EXCESS_COMPLETE_CASE_LOSS` | Remove `"ecog"` from `covariates` in `config/local.R` and state the limitation |
| 4 | **Likely blocks real data** | PD-L1 is a required complete-case covariate with a 180-day look-back and a 20% maximum loss. PD-L1 testing was uncommon for much of 2014–2017, so the gate will probably stop the run | Decide before looking at outcomes: either drop PD-L1 from the primary model (sensitivity: PD-L1-known subset) or use an "unknown" category; record the decision in the protocol |
| 5 | **Likely blocks real data** | Any patient with two pretreatment samples that disagree on STK11/KEAP1 aborts the whole run (`DISCORDANT_PRETREATMENT_ASSAYS`) | Pre-specify a rule (e.g. altered if any sample is altered) or exclude those patients and count them |
| 6 | Design | Strict same-day start of all components, no bevacizumab, NGS report strictly before treatment start, sample within 365 days: each rule is defensible but together they shrink an already small cohort | Keep, but report the flow; the bridge treats bevacizumab as "other systemic" |
| 7 | Design | No delayed-entry analysis: patients sequenced after treatment start are simply dropped | Fine for bias, costly for sample size; my pipeline offers the delayed-entry version as sensitivity |
| 8 | Minor | PD-L1 enters the propensity model as a linear percentage; `year` is only approximately recoverable from the release | Note in limitations |

## Overall
The statistics are sound and conservative; the package's main weakness is that it cannot be run on the actual release as delivered. With fixes 1–3 and the bridge it runs. Whether the original question is answerable still depends on the real cell counts, which `r_pipelines` step `feasibility` reports first.
