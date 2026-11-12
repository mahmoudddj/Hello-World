# Gate G1 — variable feasibility, strata, events, missingness, chronology
**Verdict: NOT RUN (blocked by G0).** No real cohort size, cross-tabulation, event count or missingness figure exists.

Done without patient data: structural feasibility matrix (`feasibility_report.md` §1); automated G1 report generator (`geniebpc.run.feasibility`) producing schema check, regimen classes, both cohort flows, treatment × biomarker cells with deaths and follow-up, STK11 × KEAP1 table, covariate missingness and chronology checks; tested on synthetic data and on the genieBPC simulated tables in the real file layout (schema check passed; gate correctly failed on 6 eligible simulated patients).

Pre-specified stop rule: any of the four treatment × biomarker cells with <20 patients or <10 deaths → primary causal analysis is not run (enforced in code; exit status 2).
Essential variable already known to be absent: ECOG performance status.
