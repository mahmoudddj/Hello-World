# Feasibility report (pre-data; G1 not executable)

**No real GENIE BPC counts exist in this report.** Section 1 is structural (what the release layout can support), section 2 is design-stage planning on synthetic data, section 3 is the resulting risk assessment and alternatives.

## 1. Variable-by-variable feasibility matrix
"Structure verified" = column exists in the genieBPC simulated 2.0-public data. Coverage/missingness = unknown until real data are read; the pipeline fills these in automatically (`outputs/feasibility_auto.md`).

| Need | Source table.column | Mapping (src/geniebpc/derive.py) | Structure verified | Real coverage | Concern |
|---|---|---|---|---|---|
| Advanced disease date | ca_dx_index.stage_dx_iv, dmets_stage_i_iii, dx_to_dmets_days | day 0 if stage IV at diagnosis, else day of first distant metastasis | Yes | unknown | locally advanced unresectable stage III not captured as "advanced" |
| First-line regimen, time zero | ca_drugs.regimen_drugs, dx_reg_start_int, regimen_number | earliest regimen starting on/after advanced date | Yes | unknown | drugs grouped into regimens by curation rules; grace period inherited from curation |
| Treatment strategy | ca_drugs.regimen_drugs, drugs_ct_yn | `classify_regimen` | Yes | unknown | trial drugs masked as "Investigational Drug" → excluded |
| Overall survival from time zero | ca_drugs.tt_os_g_days, os_g_status | months = days/30.4375 | Yes | unknown | death ascertainment source varies by site |
| STK11 / KEAP1 / KRAS / TP53 | mutations_extended, cna | protein-altering variants or deep deletion; loss-of-function-only sensitivity | Yes | unknown | no OncoKB annotation; missense VUS included in primary definition |
| Panel covers both genes | cpt.cpt_seq_assay_id | panel-name prefix list (config) | names only | unknown | gene content per panel NOT yet verified; UHN-48/50 and VICC-01-SOLIDTUMOR assumed non-evaluable |
| NGS timing vs time zero | cpt.dx_cpt_rep_days − ca_drugs.dx_reg_start_int | report on/before t0 (primary); delayed entry otherwise | Yes | unknown | cohort was selected on sequencing in 2014–2017 → left truncation |
| PD-L1 | prissmm_pathology.pdl1_perc, dx_path_proc_days | last result on/before t0 → <1 / 1–49 / ≥50 / unknown | Yes | unknown, expected low before 2017 | key confounder and effect modifier; high missingness likely |
| EGFR/ALK/ROS1 exclusion | mutations_extended, fusions | sensitising EGFR, ALK/ROS1 fusion on any sample | Yes | unknown | — |
| Age, sex, race, site, histology, smoking | pt_char, ca_dx_index | direct | Yes | unknown | — |
| Brain / liver metastases at baseline | ca_dx_index.dmets_*, dx_to_dmets_*_days | present and dated ≤ t0 | Yes | unknown | — |
| Calendar year of t0 | birth_year + intervals | approximate (±1 y) | Yes | — | cannot enforce exact approval dates |
| **ECOG performance status** | — | — | **Absent** | — | **major unmeasured confounder** |
| Comorbidity, steroids, tumour burden | — | — | Absent | — | unmeasured confounders |

## 2. Design-stage precision (synthetic planning, not a result)
`python -m geniebpc.precision` → `outputs/diagnostics/design_precision.csv`. Assumes ~22–25% STK11/KEAP1-altered and median OS of roughly 15 (wild-type) and 9 (altered) months.

| Hypothetical first-line cohort | SD of 12-month RMST interaction | 95% CI half-width | Interaction detectable with 80% power |
|---|---|---|---|
| 150 mono + 40 chemo-IO | 2.0 months | ±3.9 | 5.5 months |
| 250 mono + 100 chemo-IO | 1.3 | ±2.5 | 3.6 |
| 400 + 400 | 0.7 | ±1.4 | 2.0 |
| 1000 + 1000 | 0.5 | ±0.9 | 1.3 |

A 12-month RMST can differ by at most 12 months; realistic interactions are about 1–2 months. Only the last two rows can inform the question.

## 3. Assessment
1. **Positivity is doubtful in the public release.** The cohort required an NGS report in 2014–2017 (search-snippet evidence, see `genie_bpc_access_facts.md`), whereas first-line chemo-immunotherapy entered practice from mid-2017. Patients with *pretreatment* NGS who then started first-line chemo-IO are expected to be few. The pre-specified gate (≥20 patients and ≥10 deaths in each of the four treatment × biomarker cells) will probably fail. This is an expectation, not an observed count.
2. **Novelty is limited.** Elkrief 2023 (matrix L25; MSK/DFCI, n=866, propensity-weighted mono vs chemo-IO) reports no genomic factor conferring differential benefit; MSK and DFCI are also the largest BPC contributors. A small BPC analysis would be a weaker partial replication.
3. **Confounding by performance status cannot be addressed** in the public release.

**Provisional verdict:** the primary question is unlikely to be answerable in v2.0-public with publishable precision. It may be answerable in v3.1-consortium (larger, later treatment era) — access requires consortium membership or an AACR-approved concept.

## 4. Alternative questions (need researcher approval before any pivot)
| # | Question | Type | Why feasible | Main risk |
|---|---|---|---|---|
| A | Among patients given anti-PD-(L)1 monotherapy in any line vs platinum-doublet chemotherapy first line, do STK11/KEAP1 alterations modify the ICI-vs-chemotherapy difference in OS and curated progression-free survival, with explicit left-truncation handling? | Predictive (interaction), with prognostic estimates reported separately | 2014–2017 era is dominated by platinum doublets and later-line nivolumab/pembrolizumab; PRISSMM PFS is unique to BPC | Different lines of therapy → needs sequential-trial design; prior null in Flatiron (L03) |
| B | How much do NGS-timing choices (ignore timing / restrict to pretreatment / delayed-entry risk sets) change STK11/KEAP1 prognostic and interaction estimates in a clinico-genomic cohort? | Methodological, prognostic | Needs no rare treatment cell; pipeline already implements all three | Overlap with Brown 2022 (L41) and Kehl; must add the biomarker-interaction angle to be novel |
| C | STK11/KEAP1/KRAS co-alteration classes and OS/PFS on first-line platinum doublet chemotherapy | Prognostic only | Largest treatment group in this era | Cohort paper already reports part of this (L38); small co-mutant cells (n=31 KRAS+STK11) |

Suggested pairing if data access is obtained: paper 1 = primary question on v3.1-consortium if accessible, otherwise alternative A on the public release; paper 2 = alternative B (methods), which reuses the same derived cohort.
