# Running the Codex pipeline on the real release

1. In the Codex folder, add one line at the top of `R/common.R`:  `library(survival)`  (fixes the secondary Cox model).
2. Copy `config/config.R` to `config/local.R`; set `release`, `secure_root`; remove `"ecog"` from `covariates`; decide on `"pdl1"` (see review, item 4).
3. From the `r_pipelines` folder: `Rscript codex_bridge/build_codex_canonical.R <secure_root>` (uses `GENIE_BPC_DATA_DIR`).
4. Open `<secure_root>/source_mapping.csv`, check the derivations in the bridge script, then set `release` and `verified = YES` yourself, and create `dictionary_reference.txt` as its instructions ask. These are your attestations; I have not made them for you.
5. In the Codex folder: `Rscript run_pipeline.R feasibility config/local.R`, then `analyze`.

Bridge choices to be aware of: one index NSCLC per patient (first `ca_seq`); survival from `tt_os_dx_days`/`os_dx_status`; metastases and age evaluated at first-line start; bevacizumab and trial/investigational drugs count as "other systemic"; panels treated as covering STK11/KEAP1 are MSK-IMPACT, DFCI-ONCOPANEL, VICC-01-T5A/T7.
