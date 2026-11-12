# Limitations (current)
1. **No patient data analysed.** All numbers are synthetic validation or planning output.
2. Performance status, comorbidity, steroid use and tumour burden are not in the public release; residual confounding of monotherapy vs chemo-immunotherapy is expected to be strong and directional (frailer patients receive monotherapy).
3. Cohort selection on sequencing in 2014–2017 limits first-line chemo-immunotherapy numbers and creates left truncation.
4. PD-L1 missingness is expected to be high and related to calendar time.
5. Calendar time is only approximately recoverable.
6. No functional annotation of variants; missense variants of unknown significance count as altered in the primary definition.
7. Panel coverage of STK11/KEAP1 is assumed from panel names until checked.
8. MSK and DFCI overlap with previously published cohorts addressing the same contrast.
9. Literature was verified through search results only; several PMIDs/DOIs are flagged unverified and no full text was read.
10. The synthetic generator is simpler than reality (proportional hazards, independent censoring), so passing tests shows the code is correct, not that assumptions hold.
