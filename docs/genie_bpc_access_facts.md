# AACR Project GENIE BPC NSCLC: access and content facts

Compiled 2026-11-12. Evidence tiers: **[D]** read directly this session (GitHub source files); **[S]** search-engine snippet only (aacr.org, aacrjournals.org, PMC, synapse.org, cbioportal.org were not fetchable from this environment, so wording/numbers were not checked against the page itself).

## 1. Releases

| Release | Date | Access | Source |
|---|---|---|---|
| v1.1-consortium | 2020-10 | consortium only | [D] https://github.com/GENIE-BPC/genieBPC/blob/main/data-raw/synapse_tables.R |
| v2.1-consortium | 2021-08 | withdrawn (ages >89 included); replaced | [D] same file; https://github.com/GENIE-BPC/genieBPC/blob/main/NEWS.md |
| **v2.0-public** | 2022-05 (cBioPortal news: May 19, 2022) | **public** | [D] same file; https://github.com/cBioPortal/cbioportal/blob/master/docs/News-Genie.md |
| v2.2-consortium | 2024-02 | consortium only | [D] synapse_tables.R |
| v3.1-consortium | 2024-04 | consortium only | [D] synapse_tables.R |

- v2.0-public is the only public NSCLC release listed in genieBPC as of its last commit (2026-07-02, v2.2.0); later public additions were BrCa, PANC, Prostate only. [D] https://github.com/GENIE-BPC/genieBPC/blob/main/NEWS.md
- Consortium releases are restricted to GENIE consortium members and BPC pharma partners (concept sheet via genieinfo@aacr.org). [D] https://github.com/GENIE-BPC/genieBPC/blob/main/README.md
- v2.0-public: 1,846 patients, 4 institutions (MSKCC, DFCI, VICC, UHN); Synapse project syn27056172 (wiki 616601). [D] News-Genie.md (above). 2,004 tumor samples. [S] https://aacrjournals.org/clincancerres/article/29/17/3418/728542/The-GENIE-BPC-NSCLC-Cohort-A-Real-World-Repository
- Synapse folder IDs: consortium releases syn21241322, public releases syn27056700; consortium project syn21226493. [D] synapse_tables.R; README.md
- Analytic data guides: v2.0-public syn30557304; v1.1 syn23002641; v2.2 syn53463493; v3.1 syn58597690. [D] README.md
- Eligibility: sequencing report 2014-01-01 to 2017-12-31 at DFCI/MSK/VICC (UHN to 2015-12-31 per data guide; the cohort paper reportedly gives 2017 for UHN), age >=18 at sequencing, >=2 years possible follow-up, random selection; 4,982 of 6,152 screened were eligible. [S] https://www.aacr.org/wp-content/uploads/2022/05/GENIE-BPC-NSCLC-v2.0-public-Analytic-Data-Guide-1.pdf ; [S] https://aacrjournals.org/clincancerres/article-pdf/doi/10.1158/1078-0432.CCR-23-0580/3334888/ccr-23-0580.pdf

## 2. Access route

- Steps: Synapse account (not Google sign-in) -> enable 2FA (required since July 2025 for PAT use) -> accept Synapse terms -> open syn27056172 -> Data Releases -> Request Access, accept data-use terms. Synapse permission and per-release permission are separate. [D] README.md
- Programmatic access is by Personal Access Token only (username/password removed in genieBPC 2.1.0); env var `SYNAPSE_PAT`. [D] NEWS.md; https://github.com/GENIE-BPC/genieBPC/blob/main/R/pull_data_synapse.R
- Hosts hard-coded in genieBPC: `repo-prod.prod.sagebase.org` (entity bundle), `file-prod.prod.sagebase.org` (fileHandle/batch, returns a pre-signed URL), `auth-prod.prod.sagebase.org` (userProfile check). [D] pull_data_synapse.R; https://github.com/GENIE-BPC/genieBPC/blob/main/R/synapse_authentication.R
- Same three endpoints plus `www.synapse.org` are the Python client defaults. [D] https://github.com/Sage-Bionetworks/synapsePythonClient/blob/develop/synapseclient/client.py
- cBioPortal: a subset (genomics) is on the GENIE portal, study id `nsclc_public_genie_bpc`; clinical data only via Synapse. [S] https://rdrr.io/github/MSKCC-Epi-Bio/gnomeR/f/vignettes/genie-bpc-vignette.Rmd
- Terms (main GENIE data guides, not BPC-specific): no attempt to identify/contact participants; no redistribution without written permission of the Coordinating Center; cite the AACR Project GENIE Consortium, the 2017 Cancer Discovery paper and the dataset version. [S] https://www.aacr.org/wp-content/uploads/2024/07/data_guide.pdf

## 3. Files (v2.0-public), genieBPC name -> Synapse file

patient_level_dataset (pt_char, syn30358096); cancer_level_dataset_index (ca_dx_index, syn30358090); cancer_level_dataset_non_index (syn30358091); regimen_cancer_level_dataset (ca_drugs, syn30358097); imaging_level_dataset (syn30358093); pathology_report_level_dataset (syn30358095); med_onc_note_level_dataset (prissmm_md, syn30358094); cancer_panel_test_level_dataset (cpt, syn30358092); data_mutations_extended (syn30358120); data_fusions (syn30358105); data_CNA (syn30358099); variable synopsis (syn30557312). CSV for clinical, tab-delimited txt for genomic. [D] synapse_tables.R + `data/synapse_tables.rda` in https://github.com/GENIE-BPC/genieBPC
- No `data_sv` and no `ca_radtx` / `tm_level_dataset` for NSCLC v2.0-public (sv starts at v2.2; radiation at v3.1). [D] same

## 4. Variables

Confirmed as column names in genieBPC's simulated `nsclc_test_data` (structure mirrors NSCLC release; not the data guide itself) [D] https://github.com/GENIE-BPC/genieBPC/blob/main/data/nsclc_test_data.rda :
- PD-L1 (pathology dataset): pdl1_testing, pdl1_positive_any, pdl1_perc, pdl1_tclrange/tcurange, pdl1_icperc, pdl1_lcpsrange/ucpsrange, pdl1_sum (each with _2/_3).
- Smoking: ca_lung_cigarette. Stage: stage_dx, stage_dx_iv, best_ajcc_stage_cd, dmets_stage_i_iii, dx_to_dmets_days.
- OS: os_dx_status, tt_os_dx_days, os_adv_status, tt_os_adv_days (index dx); os_g_status, tt_os_g_days, os_d_status, tt_os_d1..5_days (regimen). PFS-I/PFS-M variables also present.
- Regimen: regimen_drugs, regimen_number, regimen_number_within_cancer, drugs_ct_yn, dx_reg_start_int, dx_reg_end_all_int, drugs_drug_1..5.
- NGS timing: dx_cpt_rep_days, dx_path_proc_cpt_days, path_proc_cpt_rep_days, cpt_seq_assay_id, cpt_genie_sample_id, sample_type.
- Investigational trial drugs: drug names set to "Investigational Drug", end interval = start interval. [D] https://github.com/GENIE-BPC/genieBPC/blob/main/vignettes/clinical_data_structure_vignette.Rmd
- Panels: MSK-IMPACT341/410/468, DFCI-ONCOPANEL-1/2/3/3.1, UHN-48-V1/50-V2, VICC-01-SOLIDTUMOR/T5A/T7 (UHN and VICC-SOLIDTUMOR panels are 31-50 genes; STK11/KEAP1 coverage per panel not checked). [D] `data/genie_panels.rda`

## 5. Feasibility and left truncation

- Cohort paper (Choudhury et al., Clin Cancer Res 2023;29:3418): first-line platinum without IO: KRAS alone n=113, KRAS+STK11 n=31, KRAS+KEAP1 n=4; STK11/KEAP1 worse outcomes on immunotherapy in KRAS-mutant but not KRAS-wild-type LUAD. [S] https://doi.org/10.1158/1078-0432.ccr-23-0580
- Brown et al., JAMA Oncol 2022: entry at sequencing causes left truncation; unadjusted median OS from diagnosis in stage IV NSCLC (N=727) overestimated by >1 year; trial-like subgroup 25.8 months adjusted vs 40.9 unadjusted; most common first-line regimen subgroup only N=137. [S] https://pmc.ncbi.nlm.nih.gov/articles/PMC9190030 ; code repo [D] https://github.com/slb2240/delayed_entry_clin_genom_studies
- Further methods papers listed by genieBPC: Kehl 2023 (informative cohort entry, CEBP), Kehl 2024 (NGS timing, JCO PO). [D] README.md

## Unverified / could not confirm

- Counts of first-line ICI monotherapy vs chemo-immunotherapy in BPC NSCLC: **not found**. A snippet giving 36.9% chemo / 30.5% IO mono / 27.6% chemo-IO could not be attributed to this cohort; do not use.
- Overall STK11 / KEAP1 prevalence (n, %) in the cohort: not found.
- md_ecog / md_karnof: absent from the simulated prissmm_md table (12 columns); presence in v2.0-public unverified.
- pdl1_* etc. in the actual v2.0-public data guide (only confirmed via simulated test data).
- Data cutoff / last follow-up date; diagnosis-year window; per-institution counts.
- Pre-signed download hostname (S3/CloudFront, e.g. proddata.sagebase.org.s3.amazonaws.com, data.prod.sagebase.org): not documented in any source seen; determined at runtime from the pre-signed URL.
- BPC-specific terms-of-use text and required acknowledgment wording (Synapse page not readable).
- "NSCLC 2.3-consortium" is named in NEWS.md but the package table lists v2.2-consortium; discrepancy unresolved.
- Whether any newer public NSCLC release was posted after July 2026.
