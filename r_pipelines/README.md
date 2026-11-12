# R pipelines for local execution (papers A, B, C + feasibility for the original question)

Runs entirely on your PC. Reads the GENIE BPC NSCLC v2.0-public files, writes **aggregate tables and figures only** (counts below 10 are blanked). Needs R ≥ 4.1 and the `survival` package (part of standard R installations). No other packages, no internet, no GPU.

## Files to download (GENIE BPC NSCLC v2.0-public, Synapse project syn27056172)
| File | Synapse ID | Needed for |
|---|---|---|
| `patient_level_dataset.csv` | syn30358096 | required |
| `cancer_level_dataset_index.csv` | syn30358090 | required |
| `regimen_cancer_level_dataset.csv` | syn30358097 | required |
| `pathology_report_level_dataset.csv` | syn30358095 | required (PD-L1) |
| `cancer_panel_test_level_dataset.csv` | syn30358092 | required (NGS samples, panels, report timing) |
| `data_mutations_extended.txt` | syn30358120 | required |
| `data_CNA.txt` | syn30358099 | recommended (STK11/KEAP1 deletions) |
| `data_fusions.txt` | syn30358105 | recommended (ALK/ROS1/RET) |
| `cancer_level_dataset_non_index.csv` | syn30358091 | optional, not used yet |
| `imaging_level_dataset.csv` | syn30358093 | optional, for later imaging-based questions |
| `med_onc_note_level_dataset.csv` | syn30358094 | optional, for later questions |
| variable synopsis / Analytic Data Guide | syn30557312 / syn30557304 | documentation, no patient data; useful to share with me |

IDs are taken from the public `genieBPC` R package. In R you can fetch everything at once after `genieBPC::set_synapse_credentials()`:
`genieBPC::pull_data_synapse("NSCLC", "v2.0-public", download_location = "D:/secure/genie")`.

## Steps
1. Put the downloaded files in one folder (sub-folders are fine).
2. Open a terminal in this folder and set the paths (Windows PowerShell: `$env:GENIE_BPC_DATA_DIR="D:\secure\genie"`; Linux/macOS: `export GENIE_BPC_DATA_DIR=/secure/genie`). Optionally set `GENIE_BPC_OUT_DIR` (default `results`).
3. `Rscript run_all.R test` — software checks on artificial data (1–2 minutes; should end with "49 checks passed").
4. `Rscript run_all.R feasibility` — counts only, no outcome modelling, a minute or two. It writes the `00_*.csv` files plus the 30-question scan: **`F_verdicts.csv`** (one line per candidate paper: promising / marginal / insufficient) and `F_cells.csv` (group sizes and event counts, anything below 10 blanked). **Send me these files first.**
5. `Rscript run_all.R A` (then `B`, `C`, or `all`) for the three questions that already have full pipelines. With 1,000 bootstrap replicates expect roughly 10–40 minutes per paper.
6. Look through the `results` folder before sharing it. It should contain only CSV/PNG/TXT files with aggregate numbers.

### How to read `F_verdicts.csv`
"promising" = every key comparison group has at least 50 patients and 30 events; "marginal" = at least 20 and 10; otherwise "insufficient". This is a sample-size screen only. `novelty_prior_unverified` is my rough expectation before a literature check, not a finding. A "promising" row still needs that check and a protocol before any outcome analysis.

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
