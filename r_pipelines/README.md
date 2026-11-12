# R pipelines for local execution (papers A, B, C + feasibility for the original question)

Runs entirely on your PC. Reads the GENIE BPC NSCLC v2.0-public files, writes **aggregate tables and figures only** (counts below 10 are blanked). Needs R ≥ 4.1 and the `survival` package (part of standard R installations). No other packages, no internet, no GPU.

## Steps
1. Put the release files in one folder (any sub-folders are fine): `patient_level_dataset.csv`, `cancer_level_dataset_index.csv`, `regimen_cancer_level_dataset.csv`, `pathology_report_level_dataset.csv`, `cancer_panel_test_level_dataset.csv`, `data_mutations_extended.txt`, and if available `data_CNA.txt`, `data_fusions.txt`.
2. Open a terminal in this folder and set the paths (Windows PowerShell: `$env:GENIE_BPC_DATA_DIR="D:\secure\genie"`; Linux/macOS: `export GENIE_BPC_DATA_DIR=/secure/genie`). Optionally set `GENIE_BPC_OUT_DIR` (default `results`).
3. `Rscript run_all.R test` — software checks on artificial data (about 1–2 minutes; should end with "46 checks passed").
4. `Rscript run_all.R feasibility` — counts only, no outcome modelling. **Send me these `00_*.csv` files first**; they decide which papers are viable.
5. `Rscript run_all.R A` (then `B`, `C`, or `all`). With 1,000 bootstrap replicates expect roughly 10–40 minutes per paper.
6. Look through the `results` folder before sharing it. It should contain only CSV/PNG/TXT files with aggregate numbers.

## What each step produces
| Step | Question | Main files |
|---|---|---|
| feasibility | Are there enough patients/deaths per cell for each question, when was NGS reported relative to treatment, which panels cover the genes | `00_cells_by_question.csv`, `00_regimen_classes_by_line.csv`, `00_ngs_timing_and_missingness.csv`, `00_panels_and_mutation_counts.csv` |
| A | Is the STK11/KEAP1 survival deficit larger on anti-PD-(L)1 monotherapy than on first-line platinum doublet chemotherapy? (OS and PFS; difference of RMST differences; Cox ratio of hazard ratios) | `A_status.csv`, `A_cells.csv`, `A_estimates.csv`, `A_cox.csv`, `A_balance.csv`, `A_curve_grid.csv`, `figures/` |
| B | How much do NGS-timing choices (ignore / pretreatment only / delayed entry) change survival and biomarker estimates? | `B_overall_survival_by_approach.csv`, `B_estimates.csv`, `B_quasi_independence.csv`, … |
| C | KRAS and STK11/KEAP1 co-alteration classes on first-line platinum doublet, non-squamous (prognostic only) | `C_class_counts.csv`, `C_estimates.csv`, `C_cox.csv`, … |

Reading `*_estimates.csv`: `rmst12|<stratum>|g1` is the adjusted 12-month restricted mean survival for group 1 (altered, or the named class) in that stratum; `diff12|<stratum>` is group 1 minus group 0; `did12` is the difference between strata of those differences (the interaction). CIs are patient-level bootstrap percentiles.

## Built-in rules
- An analysis stops (`STOPPED_SPARSE_CELLS` in `*_status.csv`) if any comparison cell has fewer than 20 rows or 10 events. Do not lower this to get a result.
- Primary analyses use every evaluable NGS report with delayed entry at the report date; "pretreatment only" is a sensitivity analysis. The "naive" approach appears only in paper B, as the thing being criticised.
- Biomarker = protein-altering STK11 or KEAP1 mutation or deep deletion on a covering panel; EGFR-sensitising / ALK / ROS1 excluded; trial and investigational regimens excluded.
- Performance status is not in the public release and is not adjusted for.

## Not yet verified on real data
Column names were checked against the genieBPC package's simulated NSCLC tables, not the real release. If step 4 reports missing columns in `00_schema_check.csv`, send me that file. Panel coverage of STK11/KEAP1 is assumed from panel names; `00_panels_and_mutation_counts.csv` lets us check it.
