# Gate G0 — data rights, files and dictionary

**Verdict: FAIL (blocked on access). Date: 2026-11-12.**

| Check | Result | Evidence |
|---|---|---|
| Patient-level GENIE BPC files in the workspace | None | Handoff archives contain 4 text files each (specification only); repository contained only a README |
| Researcher's connected Drive | No GENIE/BPC files | Metadata search for GENIE/BPC/NSCLC/dataset file names returned one unrelated document |
| Public mirror (cBioPortal datahub on GitHub) | Not present | `public/nsclc_public_genie_bpc` → 404; directory listing has no BPC NSCLC study |
| Synapse / cBioPortal / AACR hosts | Unreachable | Session network policy refuses `repo-prod.prod.sagebase.org`, `genie.cbioportal.org`, `www.cbioportal.org` |
| Synapse credentials | None configured | No token, no `.synapseConfig`, no data-path variable |
| Data dictionary | Partly verified | Table and column names confirmed from the public `genieBPC` R package (simulated `nsclc_test_data`, release label 2.0-public); the official Analytic Data Guide (syn30557304) could not be read |
| Licence / terms of use | Not verified first-hand | Terms known only from search snippets (no re-identification, no redistribution, citation of AACR Project GENIE); researcher must have accepted them on Synapse |

## What exists (verified from genieBPC source, see `docs/genie_bpc_access_facts.md`)
- NSCLC **v2.0-public** (2022-05; 1,846 patients; MSK, DFCI, VICC, UHN; Synapse project syn27056172) is the only public NSCLC release.
- v2.2-consortium (2024-02) and v3.1-consortium (2024-04) exist but are restricted to GENIE consortium members / BPC partners.

## Precisely missing inputs
1. The eight v2.0-public files: `patient_level_dataset.csv`, `cancer_level_dataset_index.csv`, `regimen_cancer_level_dataset.csv`,
   `pathology_report_level_dataset.csv`, `cancer_panel_test_level_dataset.csv`, `data_mutations_extended.txt`, `data_CNA.txt`, `data_fusions.txt`
   (plus, ideally, the Analytic Data Guide PDF and `genomic_information.txt`/panel gene lists).
2. Confirmation that the researcher holds a Synapse account that has accepted the GENIE BPC NSCLC terms, and authorises use in this cloud session.
3. Either (a) the files uploaded to the session, or (b) a Synapse personal access token entered through the secure credential form plus network access to the Synapse hosts.

## Decision
G1 cannot be executed on real data. Work that does not need patient data (literature, protocol v0.1, pipeline, synthetic validation) was done and is labelled as such. No real-data claim is made anywhere.
