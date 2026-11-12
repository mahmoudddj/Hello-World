# Literature notes: STK11/KEAP1 as modifiers of ICI monotherapy vs chemo-IO (first-line advanced NSCLC)

Companion to `literature_matrix.csv` (53 rows; ids L01-L53). Searches run 2026-11-12. Verification was by web search results only: PubMed, PMC, Europe PMC, Crossref and publisher pages could not be fetched from this environment, so no full text was read.

## (a) Has the exact interaction been formally tested?

Partly, once, with a null result, and without verifiable gene-level detail.

- **Elkrief 2023 (L25)**: MSK/DFCI, 866 PD-L1 >=1% lung adenocarcinomas, IPTW Cox models, mono vs chemo-IO. The abstract states that no genomic factors conferred differential survival benefit to chemo-IO vs IO; only never-smokers with PD-L1 >=50% benefited. Which genes were tested, and the STK11/KEAP1 interaction estimates, could not be confirmed.
- **Papillon-Cavanagh 2020 (L03)**: Flatiron/FMI, formal interaction test, null (STK11 OS interaction HR 1.13, KEAP1 0.98), but the contrast was anti-PD-(L)1 vs other therapy, not mono vs chemo-IO.
- **Sun 2021 (L27)** tested the same contrast by KRAS only; **Ebot 2022 (L33)** contrasted CDKN2A effects across mono-IO and chemo-IO cohorts.
- Trial evidence is indirect only (L09, L10, L24).
- FDA pooled analyses (L28, L29) and the 2026 JAMA Oncology meta-analysis (L30) could not assess STK11/KEAP1. PERSEE (L31) is now a randomised head-to-head at PD-L1 >=50%; no genotype subgroups found. INSIGNA (L32) has no results.
- No GENIE BPC analysis of this question was found (L38, L42-L44).

## (b) Honest novelty gap

The gap is real but narrower than "never studied". Not found: a pre-specified, target-trial-emulated estimate of the STK11/KEAP1-by-regimen interaction on OS, on both additive (RMST) and multiplicative scales, with left-truncation handling, in a multi-institution cohort, including PD-L1 <50% and KRAS-stratified analyses. Caveats: Elkrief already reports a null in overlapping institutions (MSK and DFCI contribute to BPC, so this is not independent replication); PERSEE changes equipoise for the main effect; and BPC is small (1,846 total; KRAS+STK11 n=31 on first-line platinum, L38), so the study will probably yield wide intervals. Frame it as estimation, not hypothesis confirmation.

## (c) Effect sizes for precision planning

- **Prevalence**: STK11 20.6%, KEAP1 19.2%, KRAS 42.5% in ICI-treated LUAD (L02); deletions add ~14% each (L15); 44% had >=1 STK11/KEAP1/KRAS mutation in Flatiron (L17).
- **OS on first-line PD-(L)1-based therapy by genotype (regimens pooled)**: STK11 mutant 8.9 vs 17.6 months; KEAP1 9.9 vs 16.3; co-mutated 8.3 vs 18.2 with none (L17). STK11-mutant on single PD-(L)1 blockade in trials: 7.8 months (L22).
- **Mono vs chemo-IO main effect (PD-L1 >=50%)**: 20.9 vs 25.0 months, HR 0.82 (L28); 19.8 vs 29.2, HR 0.74 (L30); rwPFS aHR 1.04 (L26); all-PD-L1 real-world 15.7 vs 12.3 (L33).
- **By KRAS (mutant vs wild-type)**: mono 21.1 vs 13.6; chemo-IO 20.0 vs 19.3 (L27).
- **Trial HRs in mutants**: pembrolizumab vs chemo, STK11 0.37, KEAP1 0.75 (L09); chemo-IO vs chemo, STK11 0.80, KEAP1 0.76 (L23).
- **Median OS by genotype within each regimen (mono vs chemo-IO): not found in any verified source.** This is the key missing planning input.

## (d) Not verified

- "Akinboro, J Clin Oncol 2024": only the ASCO 2022 abstract was found (L28).
- A "Cordeiro de Lima" paper: not found; the pan-cancer STK11 paper is Krishnamurthy et al. (L07).
- The intended Aggarwal 2020 PD-L1/STK11 paper: L12 may be a different paper.
- PERSEE HR/medians (secondary coverage only) and DOI; JAMA Oncology 2026 DOI (L30); POSEIDON subgroup HRs (L21); Ricciuti 2022 HRs (L02).
- First authors of L15, L23, L24, L44; year/DOI of L34; database identity in L33; whether L42 used BPC or the main GENIE registry; L43 results (sources conflicting).
- Many PMIDs are marked "unverified" (e.g. L05, L06, L11, L26, L27, L48-L51, L53).
- VanderWeele 2009 (interaction vs effect modification) was not confirmed; L51 and L52 are substitutes.
