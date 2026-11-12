# Study-specific analyses. Each writes aggregate CSV/PNG files to cfg$out_dir and returns the result objects.

# ---- Step 0: feasibility for all four questions (counts only, no outcome modelling) ----
feasibility <- function(t, schema, cfg) {
  write_out(schema, "00_schema_check", cfg, character(0))
  r <- regimens(t, cfg)
  cl <- as.data.frame(table(line = pmin(r$line, 3), class = r$class), responseName = "n", stringsAsFactors = FALSE)
  cl$line <- ifelse(cl$line == "3", "3+", cl$line); write_out(cl, "00_regimen_classes_by_line", cfg)
  s <- sample_alterations(t)
  pan <- stats::aggregate(cbind(n = 1, STK11_alt = s$STK11_alt, KEAP1_alt = s$KEAP1_alt, KRAS_alt = s$KRAS_alt) ~ cpt_seq_assay_id + covers_both, data = s, FUN = sum)
  write_out(pan, "00_panels_and_mutation_counts", cfg, c("n", "STK11_alt", "KEAP1_alt", "KRAS_alt"))   # a panel with many samples and 0 KEAP1 is probably not covering it
  out <- list()
  for (mode in c("pretreatment", "anytime")) {
    d <- add_biomarker(r, t, cfg, mode)
    q <- list(P0_ICI_mono_vs_ICI_chemo_first_line = d[d$line == 1 & d$class %in% c("ICI_mono", "ICI_chemo"), ],
              A_chemo_first_line = d[d$line == 1 & d$class == "platinum_doublet", ],
              A_ICI_mono_first_use = first_ici(d),
              C_chemo_first_line_nonsquamous = { z <- add_covariates(d[d$line == 1 & d$class == "platinum_doublet", ], t, cfg); z[z$histology != "Squamous", ] })
    for (nm in names(q)) { z <- q[[nm]]; if (!nrow(z)) next
      a <- stats::aggregate(cbind(n = 1, deaths = z$death, kras = z$KRAS) ~ class + biomarker, data = z, FUN = sum)
      out[[length(out) + 1]] <- cbind(question = nm, ngs = mode, a) }
  }
  write_out(do.call(rbind, out), "00_cells_by_question", cfg, c("n", "deaths", "kras"))
  d <- add_covariates(add_biomarker(r[r$line == 1, ], t, cfg, "anytime"), t, cfg)
  qs <- stats::quantile(d$ngs_rel_m, c(.1, .25, .5, .75, .9), na.rm = TRUE)
  tim <- data.frame(quantity = c("first-line starts with evaluable NGS", "NGS report on/before first-line start (%)", paste0("months from first-line start to NGS report, p", c(10, 25, 50, 75, 90)),
                                 "calendar year of first-line start: p10", "p50", "p90", "PD-L1 unknown at first-line start (%)"),
                    value = c(nrow(d), 100 * mean(d$ngs_rel_m <= 0), qs, stats::quantile(d$year_t0, c(.1, .5, .9), na.rm = TRUE), 100 * mean(d$pdl1_cat == "Unknown")))
  write_out(tim, "00_ngs_timing_and_missingness", cfg, character(0))
  invisible(TRUE)
}
first_ici <- function(d) { z <- d[d$class == "ICI_mono", ]; z <- z[order(z$record_id, z$t0), ]; z[!duplicated(z$record_id), ] }

# ---- Paper A: is the STK11/KEAP1 survival deficit larger on anti-PD-(L)1 monotherapy than on platinum doublet chemotherapy? ----
paper_A_data <- function(t, cfg, ngs = "anytime", lof_only = FALSE) {
  d <- add_biomarker(regimens(t, cfg), t, cfg, ngs, lof_only)
  z <- rbind(d[d$line == 1 & d$class == "platinum_doublet", ], first_ici(d))
  z <- add_covariates(z, t, cfg); z$tx <- ifelse(z$class == "ICI_mono", "2_ICI_mono", "1_platinum_doublet"); z$line_gt1 <- as.integer(z$line > 1); z
}
paper_A <- function(t, cfg) {
  cov <- c(COVARS, "KRAS", "line_gt1"); res <- list()
  run <- function(tag, d, outcome, taus, adjust = TRUE) { res[[tag]] <<- contrast_analysis(d, "biomarker", "tx", outcome, taus, cov, cfg, tag, adjust) }
  d <- paper_A_data(t, cfg, "anytime")
  run("os_primary_delayed_entry", d, "os", cfg$taus_os)
  run("pfs_delayed_entry", d, "pfs", cfg$taus_pfs)
  run("os_unadjusted", d, "os", cfg$taus_os, adjust = FALSE)
  run("os_pretreatment_ngs_only", paper_A_data(t, cfg, "pretreatment"), "os", cfg$taus_os)
  run("os_loss_of_function_only", paper_A_data(t, cfg, "anytime", lof_only = TRUE), "os", cfg$taus_os)
  run("os_nonsquamous", d[d$histology != "Squamous", ], "os", cfg$taus_os)
  run("os_KRAS_mutant", d[d$KRAS == 1, ], "os", cfg$taus_os)
  for (g in c("STK11", "KEAP1")) { dg <- d; dg$biomarker <- dg[[g]]; run(paste0("os_", g, "_alone"), dg, "os", cfg$taus_os) }
  write_paper(res, "A", cfg); invisible(res)
}

# ---- Paper B: how much do NGS-timing choices change biomarker-survival estimates? ----
paper_B <- function(t, cfg) {
  cov <- c(COVARS, "KRAS"); res <- list(); r <- regimens(t, cfg)
  cohorts <- list(first_line_any_regimen = function(d) d[d$line == 1 & !grepl("trial", d$class), ],
                  first_line_platinum_doublet = function(d) d[d$line == 1 & d$class == "platinum_doublet", ],
                  ICI_mono_first_use = first_ici)
  desc <- list()
  for (cn in names(cohorts)) for (mode in c("naive", "pretreatment", "anytime")) {
    d <- add_covariates(cohorts[[cn]](add_biomarker(r, t, cfg, mode)), t, cfg); d$all <- cn
    tag <- paste(cn, c(naive = "1_naive_ignore_timing", pretreatment = "2_pretreatment_only", anytime = "3_delayed_entry")[mode], sep = "__")
    x <- contrast_analysis(d, "biomarker", "all", "os", cfg$taus_os, cov, cfg, tag); res[[tag]] <- x
    z <- set_outcome(d, "os"); k <- km(z$time, z$event, NULL, if (any(z$entry > 0)) z$entry else NULL)
    med <- if (any(k$surv <= .5)) k$time[which(k$surv <= .5)[1]] else NA
    desc[[length(desc) + 1]] <- data.frame(cohort = cn, approach = mode, n = nrow(z), events = sum(z$event), median_os_months = med, rmst12 = rmst_km(k, 12),
      immortal_person_time_pct = if (mode == "naive") 100 * sum(pmin(pmax(z$ngs_rel_m, 0), z$time)) / sum(z$time) else 0)
  }
  write_out(do.call(rbind, desc), "B_overall_survival_by_approach", cfg)
  # dependence between entry time and survival (quasi-independence check) among late-tested patients
  d <- add_covariates(cohorts$first_line_any_regimen(add_biomarker(r, t, cfg, "anytime")), t, cfg); z <- set_outcome(d, "os")
  qi <- tryCatch({ zz <- z[z$entry > 0, ]; zz$delay <- zz$entry; f <- coxph(Surv(entry, time, event) ~ delay, data = zz); s <- summary(f)$coefficients
    data.frame(term = "HR per month of delay in NGS report (entry > 0 only)", hr = s[1, 2], p = s[1, 5], n = sum(z$entry > 0)) }, error = function(e) data.frame(term = "NOT_ESTIMABLE"))
  write_out(qi, "B_quasi_independence", cfg)
  # treatment-class interaction from paper A under each approach
  for (mode in c("naive", "pretreatment", "anytime")) {
    tag <- paste0("interaction__", mode); x <- contrast_analysis(paper_A_data(t, cfg, mode), "biomarker", "tx", "os", cfg$taus_os, c(cov, "line_gt1"), cfg, tag)
    res[[tag]] <- x
  }
  write_paper(res, "B", cfg); invisible(res)
}

# ---- Paper C: KRAS / STK11-KEAP1 co-alteration classes on first-line platinum doublet chemotherapy (prognostic only) ----
paper_C <- function(t, cfg) {
  d <- add_covariates({ z <- add_biomarker(regimens(t, cfg), t, cfg, "anytime"); z[z$line == 1 & z$class == "platinum_doublet", ] }, t, cfg)
  d <- d[d$histology != "Squamous", ]; d$all <- "first_line_platinum_doublet_nonsquamous"; res <- list()
  d$co_class <- ifelse(d$KRAS == 1 & d$biomarker == 1, "KRAS+SK", ifelse(d$KRAS == 1, "KRAS_only", ifelse(d$biomarker == 1, "SK_only", "neither")))
  tb <- as.data.frame(table(co_class = d$co_class), responseName = "n", stringsAsFactors = FALSE); tb$deaths <- as.numeric(tapply(d$death, d$co_class, sum)[tb$co_class]); write_out(tb, "C_class_counts", cfg, c("n", "deaths"))
  for (outcome in c("os", "pfs")) {
    taus <- if (outcome == "os") cfg$taus_os else cfg$taus_pfs
    for (k in c("KRAS_only", "SK_only", "KRAS+SK")) {            # each class against 'neither'
      z <- d[d$co_class %in% c(k, "neither"), ]; z$g <- as.integer(z$co_class == k)
      tag <- paste0(outcome, "__", gsub("\\+", "_and_", k), "_vs_neither"); res[[tag]] <- contrast_analysis(z, "g", "all", outcome, taus, COVARS, cfg, tag)
    }
    d$kras_stratum <- ifelse(d$KRAS == 1, "2_KRAS_mutant", "1_KRAS_wildtype")   # is the STK11/KEAP1 effect KRAS-dependent?
    tag <- paste0(outcome, "__SK_by_KRAS_interaction"); res[[tag]] <- contrast_analysis(d, "biomarker", "kras_stratum", outcome, taus, COVARS, cfg, tag)
  }
  write_paper(res, "C", cfg); invisible(res)
}
