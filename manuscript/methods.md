# Methods (draft of the pre-specified design; to be revised to match what is actually executed)

**Data source.** AACR Project GENIE Biopharma Collaborative NSCLC cohort [L38]; release and access date to be inserted once data are obtained. Curated with the PRISSMM framework [L39].

**Design.** Target trial emulation; protocol in `protocol/target_trial.md`. Eligible patients had stage IV or recurrent metastatic NSCLC, started first-line anti-PD-(L)1 monotherapy or anti-PD-(L)1 plus platinum chemotherapy outside a trial, had an NGS report from a panel covering STK11 and KEAP1 dated on or before treatment start, and had no sensitising EGFR mutation or ALK/ROS1 fusion. Time zero was the start of first-line therapy.

**Biomarker.** STK11 and/or KEAP1 protein-altering mutation or deep deletion versus both wild-type.

**Outcome and estimands.** Overall survival; difference in 12- and 24-month restricted mean survival time between strategies within each biomarker stratum, and the between-stratum difference of those differences (additive interaction) [L46, L52].

**Statistical analysis.** Stabilised inverse-probability-of-treatment weighting with weighted Kaplan–Meier estimation [L49], bootstrap confidence intervals, doubly robust pseudo-observation estimation [L48], and a secondary weighted Cox model with a product term; details, stop rules and sensitivity analyses in `protocol/SAP.md`. Analyses used Python (versions in `requirements.lock`); code and synthetic test fixtures are public.

**Ethics / data use.** To be completed from the release's terms of use.
