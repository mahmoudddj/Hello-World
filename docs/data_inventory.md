# Data inventory

## Files physically available on 2026-11-12
| Item | Content | Patient data? |
|---|---|---|
| `GENIE_BPC_Claude_Code_Independent_Handoff.zip` | README, START_PROMPT, CLAUDE.md, SCIENTIFIC_BRIEF | No |
| `GENIE_BPC_Codex_Independent_Handoff.zip` | Same specification addressed to the other agent (identical brief; no code, decisions or results) | No |
| `genieBPC` R package (public GitHub, fetched for schema only) | Simulated `nsclc_test_data` (200 fake patients), Synapse ID table, panel list | No (simulated) |

Isolation note: the Codex archive was supplied by the researcher in the same message; it contains only the shared specification. No work product of the other agent has been seen.

## Expected source tables (GENIE BPC NSCLC v2.0-public) and how this study uses them
Column names verified against the genieBPC simulated data (structure only); to be re-verified on the real release by `python -m geniebpc.run` (writes `outputs/diagnostics/schema_check.csv`).

| Logical table | Release file | Key columns used | Used for |
|---|---|---|---|
| pt_char | patient_level_dataset.csv | record_id, institution, birth_year, naaccr_sex_code, naaccr_race_code_primary | demographics, site, approximate calendar year |
| ca_dx_index | cancer_level_dataset_index.csv | ca_seq, age_dx, dob_ca_dx_days, stage_dx_iv, dmets_stage_i_iii, dx_to_dmets_days, ca_hist_adeno_squamous, ca_lung_cigarette, dmets_brain/liver, dx_to_dmets_brain/liver_days | advanced-disease date, histology, smoking, baseline metastatic sites |
| ca_drugs | regimen_cancer_level_dataset.csv | regimen_number, regimen_drugs, drugs_ct_yn, dx_reg_start_int, os_g_status, tt_os_g_days | treatment strategy, time zero, overall survival |
| prissmm_pathology | pathology_report_level_dataset.csv | dx_path_proc_days, pdl1_perc | baseline PD-L1 |
| cpt | cancer_panel_test_level_dataset.csv | cpt_genie_sample_id, cpt_seq_assay_id, dx_cpt_rep_days | NGS panel and report timing |
| mutations_extended | data_mutations_extended.txt | Hugo_Symbol, Variant_Classification, HGVSp_Short, Tumor_Sample_Barcode | STK11, KEAP1, KRAS, TP53, EGFR |
| cna | data_CNA.txt | gene × sample matrix | STK11/KEAP1 deep deletion |
| fusions | data_fusions.txt | Hugo_Symbol, Tumor_Sample_Barcode | ALK/ROS1 exclusion |

## Variables requested by the brief that appear NOT to exist in the public release
- **ECOG / Karnofsky performance status**: no `md_ecog`/`md_karnof` column in the med-onc table of the 2.0-public structure. This is the single most important unmeasured confounder for monotherapy vs chemo-immunotherapy.
- **Exact calendar dates**: only intervals from birth/diagnosis plus `birth_year`; calendar year of treatment start can be reconstructed to about ±1 year.
- Structural variants file (`data_sv`) and radiation data start in later consortium releases.
Exact source-to-derived mapping: `src/geniebpc/derive.py` (single place where raw columns become analysis variables).
