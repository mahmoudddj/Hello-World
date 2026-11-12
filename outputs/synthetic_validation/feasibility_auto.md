# G1 feasibility report (SYNTHETIC DATA - pipeline validation only, not study results)

## Schema check

Required tables/columns present: **True**

| table   | file   | column   | present   | optional   |
|:--------|:-------|:---------|:----------|:-----------|
| cna     |        | <table>  | False     | True       |
| fusions |        | <table>  | False     | True       |

## First-line regimen classes (all advanced patients)

| first_line_class                 |    n |
|:---------------------------------|-----:|
| ICI_mono                         | 1166 |
| ICI_chemo                        | 1092 |
| other:no_pd1_pdl1                |  873 |
| other:trial_or_investigational   |  305 |
| other:ici_plus_nonprotocol_agent |  290 |

## Cohort flow (primary: NGS report on/before time zero)

| step                                                                                                                    |    n |
|:------------------------------------------------------------------------------------------------------------------------|-----:|
| Patients with an index NSCLC diagnosis                                                                                  | 3726 |
| Advanced disease (stage IV at diagnosis or later distant metastasis)                                                    | 3726 |
| Received >=1 systemic regimen for advanced disease                                                                      | 3726 |
| First-line regimen is anti-PD-(L)1 monotherapy or anti-PD-(L)1 + platinum chemotherapy (non-trial, non-investigational) | 2258 |
| STK11 and KEAP1 evaluable on a covering NGS panel reported on/before treatment start                                    | 1645 |
| No sensitising EGFR mutation or ALK/ROS1 fusion detected                                                                | 1564 |
| Valid overall-survival follow-up after time zero                                                                        | 1563 |

## Cohort flow (sensitivity: any-time NGS with delayed entry)

| step                                                                                                                    |    n |
|:------------------------------------------------------------------------------------------------------------------------|-----:|
| Patients with an index NSCLC diagnosis                                                                                  | 3726 |
| Advanced disease (stage IV at diagnosis or later distant metastasis)                                                    | 3726 |
| Received >=1 systemic regimen for advanced disease                                                                      | 3726 |
| First-line regimen is anti-PD-(L)1 monotherapy or anti-PD-(L)1 + platinum chemotherapy (non-trial, non-investigational) | 2258 |
| STK11 and KEAP1 evaluable on a covering NGS panel (any time)                                                            | 2063 |
| No sensitising EGFR mutation or ALK/ROS1 fusion detected                                                                | 1965 |
| Valid overall-survival follow-up after time zero                                                                        | 1964 |
| Under observation after NGS report (delayed entry)                                                                      | 1964 |

## Treatment x biomarker cells (primary cohort)

| status                | arm       |   n |   deaths |   person_months |   median_followup_m |   deaths_by_12m |   at_risk_beyond_12m |   deaths_by_24m |   at_risk_beyond_24m |
|:----------------------|:----------|----:|---------:|----------------:|--------------------:|----------------:|---------------------:|----------------:|---------------------:|
| STK11/KEAP1 altered   | ICI_chemo | 190 |      124 |        2.71e+03 |               11    |              86 |                   88 |             112 |                   34 |
| STK11/KEAP1 altered   | ICI_mono  | 158 |      128 |        1.57e+03 |                7.69 |             102 |                   48 |             123 |                   14 |
| STK11/KEAP1 wild-type | ICI_chemo | 578 |      365 |        8.73e+03 |               12.5  |             220 |                  300 |             326 |                  122 |
| STK11/KEAP1 wild-type | ICI_mono  | 637 |      412 |        9.43e+03 |               11.4  |             269 |                  307 |             372 |                  136 |

## STK11 x KEAP1

| STK11     |   KEAP1 alt |   KEAP1 wt |
|:----------|------------:|-----------:|
| STK11 alt |         101 |        106 |
| STK11 wt  |         141 |       1215 |

## Covariate missingness (primary cohort)

| variable       |   missing_or_unknown_pct |
|:---------------|-------------------------:|
| age            |                      0   |
| year_t0        |                      0   |
| female         |                      0   |
| brain_mets     |                      0   |
| liver_mets     |                      0   |
| de_novo        |                      0   |
| prior_systemic |                      0   |
| KRAS           |                      0   |
| TP53           |                      0   |
| institution    |                      0   |
| race           |                      0   |
| histology      |                      0   |
| smoking        |                      0   |
| pdl1_cat       |                     20.7 |

## Chronology

| quantity                                                 |      value |
|:---------------------------------------------------------|-----------:|
| NGS report on/before t0 (any-time cohort), %             |  79.6      |
| median days from t0 to first evaluable NGS report        | -47        |
| PD-L1 result dated after t0 used as baseline (must be 0) |   0        |
| year of t0: min                                          |   2.02e+03 |
| year of t0: max                                          |   2.02e+03 |

## Positivity / precision gate

**PASS** - all four treatment x biomarker cells meet the pre-specified minimum.
