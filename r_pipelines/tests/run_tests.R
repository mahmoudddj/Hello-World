# Software checks on synthetic data only. Run: Rscript run_all.R test
checks <- 0L
expect <- function(x, what) { if (!isTRUE(x)) stop("TEST FAILED: ", what, call. = FALSE); checks <<- checks + 1L }
near <- function(a, b, tol) isTRUE(abs(a - b) <= tol)

# --- regimen classification
cases <- c("Pembrolizumab" = "ICI_mono", "Nivolumab(BMS936558,MDX1106,NIVO,ONO4538,Opdivo)" = "ICI_mono",
  "Carboplatin, Pembrolizumab, Pemetrexed Disodium" = "ICI_chemo", "Atezolizumab, Bevacizumab, Carboplatin, Paclitaxel" = "ICI_chemo",
  "Carboplatin, Pemetrexed Disodium" = "platinum_doublet", "Bevacizumab, Carboplatin, Paclitaxel" = "platinum_doublet",
  "Cisplatin, Etoposide" = "platinum_doublet", "Carboplatin" = "other:no_pd1_pdl1", "Docetaxel" = "other:no_pd1_pdl1",
  "Ipilimumab, Nivolumab" = "other:ici_plus_nonprotocol_agent", "Erlotinib Hydrochloride" = "other:no_pd1_pdl1",
  "Pembrolizumab, Pemetrexed Disodium" = "other:ici_nonplatinum_combo", "Investigational Drug, Pembrolizumab" = "other:trial_or_investigational")
for (k in names(cases)) expect(classify_regimen(k) == cases[[k]], paste("classify", k))
expect(classify_regimen("Pembrolizumab", "Yes") == "other:trial_or_investigational", "trial flag")
expect(classify_regimen(NA) == "other:missing", "missing regimen")

# --- Kaplan-Meier / RMST against survival::survfit (plain, weighted, delayed entry)
set.seed(1); n <- 400; tt <- round(rexp(n, .1), 1) + .1; ev <- rbinom(n, 1, .7); w <- runif(n, .5, 3); en <- pmin(runif(n, 0, 3), tt - .05); en[1:150] <- 0
chk <- function(k, f, what) { s <- summary(f, times = k$time)$surv; expect(max(abs(s - k$surv)) < 1e-8, what) }
chk(km(tt, ev), survfit(Surv(tt, ev) ~ 1), "KM plain")
chk(km(tt, ev, w), survfit(Surv(tt, ev) ~ 1, weights = w), "KM weighted")
chk(km(tt, ev, NULL, en), survfit(Surv(en, tt, ev) ~ 1), "KM delayed entry")
chk(km(tt, ev, w, en), survfit(Surv(en, tt, ev) ~ 1, weights = w), "KM weighted + delayed entry")
expect(near(rmst(c(1, 2, 3), c(1, 1, 0), 3), 2, 1e-12), "RMST analytic")
expect(near(rmst(c(1, 2, 3), c(1, 1, 0), 3, w = c(2, 1, 1)), 1.75, 1e-12), "RMST weighted analytic")
expect(near(rmst(c(5, 9), c(0, 0), 4), 4, 1e-12), "RMST no events")
expect(near(rmst(tt, ev, 12), summary(survfit(Surv(tt, ev) ~ 1), rmean = 12)$table[["rmean"]], 1e-6), "RMST vs survfit")

# --- privacy
x <- suppress(data.frame(n = c(100, 4), events = c(50, 2)), c("n", "events"), 10)
expect(is.na(x$n[2]) && is.na(x$events[2]) && x$n[1] == 100 && x$suppressed[2] == "<10", "small-cell suppression")
x <- suppress(data.frame(n = c(100, 40), events = c(50, 3)), c("n", "events"), 10); expect(is.na(x$n[2]), "whole row blanked")

# --- derivations and chronology on synthetic raw files (through the file reader)
s <- synth_tables(40000, seed = 5); dir <- file.path(tempdir(), paste0("syn_", Sys.getpid())); write_synth(s$tables, dir)
L <- load_tables(dir); expect(schema_ok(L$schema), "schema ok"); t <- L$tables
bad <- L$tables; bad$ca_drugs$tt_os_g_days <- NULL; write_synth(bad, file.path(dir, "bad"))
expect(!schema_ok(load_tables(file.path(dir, "bad"))$schema), "missing column detected")
r <- regimens(t, cfg)
expect(all(r$line[r$regimen_number == 2] == 2) && all(r$class[r$regimen_number == 2] == "ICI_mono"), "line numbering")
dp <- add_covariates(add_biomarker(r, t, cfg, "pretreatment"), t, cfg)
expect(all(dp$ngs_rel_m <= 0) && all(dp$entry_m == 0), "pretreatment mode uses no post-baseline NGS")
expect(all(is.na(dp$pdl1_day) | dp$pdl1_day <= dp$t0) && !any(dp$pdl1_cat == ">=50%" & is.na(dp$pdl1_day)), "no post-baseline PD-L1")
da <- add_biomarker(r, t, cfg, "anytime")
expect(any(da$entry_m > 0) && all(da$entry_m >= 0), "delayed entry present in anytime mode")
expect(all(add_biomarker(r, t, cfg, "naive")$entry_m == 0), "naive mode has no entry adjustment")
expect(!any(da$site == "UHN" | grepl("UHN", da$record_id)) || !any(t$cpt$cpt_seq_assay_id[match(da$record_id, t$cpt$record_id)] == "UHN-48-V1"), "hotspot panel not evaluable")
egfr <- sub("-S1$", "", t$mutations_extended$Tumor_Sample_Barcode[t$mutations_extended$Hugo_Symbol == "EGFR"])
expect(!any(da$record_id %in% egfr), "EGFR-mutant excluded")
fi <- first_ici(da); expect(!anyDuplicated(fi$record_id), "one ICI row per patient")
lof <- add_biomarker(r, t, cfg, "anytime", lof_only = TRUE); expect(sum(lof$KEAP1) == 0 && sum(lof$STK11) == sum(da$STK11), "loss-of-function definition")

# --- planted truth: paper A contrast, and the effect of ignoring NGS timing
cov <- c(COVARS, "KRAS", "line_gt1"); est <- list(); lvl <- list()
for (mode in c("naive", "pretreatment", "anytime")) { d <- set_outcome(paper_A_data(t, cfg, mode), "os")
  e <- cell_rmst(d, "biomarker", "tx", c(12, 24), ipw(d, "biomarker", cov, "tx")); est[[mode]] <- e; lvl[[mode]] <- e[["rmst12|1_platinum_doublet|g0"]] }
tr <- s$truth
cat(sprintf("  truth: chemo12 %.2f ici12 %.2f did12 %.2f did24 %.2f | delayed entry: %.2f %.2f %.2f %.2f\n", tr$chemo12, tr$ici12, tr$did12, tr$did24,
            est$anytime[["diff12|1_platinum_doublet"]], est$anytime[["diff12|2_ICI_mono"]], est$anytime[["did12"]], est$anytime[["did24"]]))
expect(near(est$anytime[["diff12|1_platinum_doublet"]], tr$chemo12, .35), "chemo contrast recovered (delayed entry)")
expect(near(est$anytime[["diff12|2_ICI_mono"]], tr$ici12, .35), "ICI contrast recovered (delayed entry)")
expect(near(est$anytime[["did12"]], tr$did12, .4) && near(est$anytime[["did24"]], tr$did24, .8), "difference-in-differences recovered")
expect(near(est$pretreatment[["did12"]], tr$did12, .5), "pretreatment-only also recovers it")
expect(lvl$naive > lvl$anytime + .4 && near(lvl$anytime, lvl$pretreatment, .25), "ignoring NGS timing inflates survival; delayed entry removes it")

# --- stop rule and end-to-end run with export checks
small <- synth_tables(300, seed = 2)$tables
x <- contrast_analysis(paper_A_data(small, cfg, "pretreatment"), "biomarker", "tx", "os", 12, cov, cfg, "small")
expect(x$status == "STOPPED_SPARSE_CELLS" && is.null(x$estimates), "sparse cells stop the analysis")
c2 <- cfg; c2$boot <- 25L; c2$out_dir <- file.path(dir, "out"); t3 <- synth_tables(5000, seed = 9)$tables
feasibility(t3, L$schema, c2); A <- paper_A(t3, c2); B <- paper_B(t3, c2); C <- paper_C(t3, c2)
expect(A$os_primary_delayed_entry$status == "COMPLETED", "paper A end-to-end")
expect(length(B) == 12 && length(C) == 8, "papers B and C end-to-end")
files <- list.files(c2$out_dir, recursive = TRUE, full.names = TRUE); csv <- files[grepl("csv$", files)]
expect(length(csv) >= 20 && !any(vapply(csv, function(f) any(grepl("GENIE-", readLines(f), fixed = TRUE)), logical(1))), "no patient identifiers in any exported table")
e1 <- boot_cluster(set_outcome(paper_A_data(t3, c2, "anytime"), "os"), function(z) cell_rmst(z, "biomarker", "tx", 12, rep(1, nrow(z))), 10, 1)
e2 <- boot_cluster(set_outcome(paper_A_data(t3, c2, "anytime"), "os"), function(z) cell_rmst(z, "biomarker", "tx", 12, rep(1, nrow(z))), 10, 1)
expect(identical(e1, e2), "bootstrap is reproducible")
unlink(dir, recursive = TRUE)
cat(checks, "checks passed (synthetic data only).\n")
