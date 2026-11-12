# Synthetic tables in the GENIE BPC file layout with KNOWN effects. For software validation only - never a result.
synth_tables <- function(n = 6000, seed = 1, bM = c(platinum_doublet = 0.3, ICI_mono = 0.7, ICI_chemo = 0.25),
                         late_frac = 0.5, cfg = list(days_per_month = 30.4375)) {
  set.seed(seed); dpm <- cfg$days_per_month
  id <- sprintf("GENIE-SYN-%06d", seq_len(n))
  site <- sample(c("MSK", "DFCI", "VICC", "UHN"), n, TRUE, c(.45, .35, .12, .08))
  age <- pmin(89, pmax(35, rnorm(n, 66, 9))); female <- rbinom(n, 1, .52); adeno <- rbinom(n, 1, .75)
  liver <- rbinom(n, 1, .15); brain <- rbinom(n, 1, .2)
  kras <- rbinom(n, 1, ifelse(adeno == 1, .35, .05))
  stk <- rbinom(n, 1, ifelse(kras == 1, .28, .10) * ifelse(adeno == 1, 1, .3) * ifelse(liver == 1, 1.5, 1))
  keap <- rbinom(n, 1, ifelse(stk == 1, .45, .10)); M <- as.integer(stk + keap > 0)
  pdl1 <- ifelse(runif(n) < .4, NA, round(100 * rbeta(n, .5, 1.2)))
  de_novo <- rbinom(n, 1, .7)
  cls <- sample(c("platinum_doublet", "ICI_mono", "ICI_chemo", "other"), n, TRUE, c(.5, .2, .15, .15))
  xb <- .5 * liver + .3 * brain + .02 * (age - 66)
  haz <- function(m, cl, x = xb) .05 * exp(ifelse(cl == "other", .3, bM[cl]) * m + ifelse(cl == "ICI_chemo", -.2, 0) + x)
  u1 <- rexp(n); u2 <- rexp(n); s2 <- runif(n, 4, 10); gets2 <- cls == "platinum_doublet" & runif(n) < .6
  death_m <- function(m) { T1 <- u1 / haz(m, cls); ifelse(gets2 & T1 > s2, s2 + u2 / haz(m, "ICI_mono"), T1) }
  d0 <- death_m(0); d1 <- death_m(1); death <- ifelse(M == 1, d1, d0)
  has2 <- gets2 & (u1 / haz(M, cls)) > s2
  cens <- runif(n, 8, 48)
  rep_rel <- ifelse(runif(n) < late_frac, runif(n, 15, 400), -runif(n, 1, 120))
  keep <- pmin(death, cens) * dpm > pmax(rep_rel, 0)         # only patients alive (and followed) at the NGS report enter GENIE
  truth <- list()
  for (tau in c(12, 24)) {
    i1 <- cls == "platinum_doublet"; truth[[sprintf("chemo%g", tau)]] <- mean(pmin(d1[i1], tau) - pmin(d0[i1], tau))
    ia <- cls == "ICI_mono"; t1 <- c(d1[ia], (u2 / haz(1, "ICI_mono"))[has2]); t0 <- c(d0[ia], (u2 / haz(0, "ICI_mono"))[has2])
    truth[[sprintf("ici%g", tau)]] <- mean(pmin(t1, tau) - pmin(t0, tau))
    truth[[sprintf("did%g", tau)]] <- truth[[sprintf("ici%g", tau)]] - truth[[sprintf("chemo%g", tau)]]
  }
  dmets <- ifelse(de_novo == 1, NA, sample(120:900, n, TRUE)); adv <- ifelse(de_novo == 1, 0, dmets); t0 <- adv + sample(14:60, n, TRUE)
  names1 <- c(platinum_doublet = "Carboplatin, Pemetrexed Disodium", ICI_mono = "Pembrolizumab",
              ICI_chemo = "Carboplatin, Pembrolizumab, Pemetrexed Disodium", other = "Erlotinib Hydrochloride")
  pfs <- function(os, h) { p <- rexp(length(os), 2.2 * h); pmin(os, p) }
  os1 <- pmin(death, cens); st1 <- as.integer(death <= cens)
  reg1 <- data.frame(record_id = id, ca_seq = 0, regimen_number = 1, regimen_drugs = unname(names1[cls]), drugs_ct_yn = "No",
    dx_reg_start_int = t0, os_g_status = st1, tt_os_g_days = round(os1 * dpm), stringsAsFactors = FALSE)
  p1 <- pfs(os1, haz(M, cls)); p1 <- ifelse(has2, pmin(p1, s2), p1)
  reg1$pfs_i_or_m_g_status <- as.integer(p1 < os1 | st1 == 1); reg1$tt_pfs_i_or_m_g_days <- round(p1 * dpm)
  r2 <- which(has2 & cens > s2)
  os2 <- pmin(death, cens)[r2] - s2[r2]; p2 <- pfs(os2, haz(M[r2], "ICI_mono", xb[r2]))
  reg2 <- data.frame(record_id = id[r2], ca_seq = 0, regimen_number = 2, regimen_drugs = "Nivolumab", drugs_ct_yn = "No",
    dx_reg_start_int = t0[r2] + round(s2[r2] * dpm), os_g_status = st1[r2], tt_os_g_days = round(os2 * dpm),
    pfs_i_or_m_g_status = as.integer(p2 < os2 | st1[r2] == 1), tt_pfs_i_or_m_g_days = round(p2 * dpm), stringsAsFactors = FALSE)
  k <- keep; kid <- id[k]
  smp <- paste0(id, "-S1")
  mut <- rbind(data.frame(Hugo_Symbol = "STK11", Variant_Classification = "Nonsense_Mutation", HGVSp_Short = "p.X1", Tumor_Sample_Barcode = smp[stk == 1 & k]),
               data.frame(Hugo_Symbol = "KEAP1", Variant_Classification = "Missense_Mutation", HGVSp_Short = "p.X2", Tumor_Sample_Barcode = smp[keap == 1 & k]),
               data.frame(Hugo_Symbol = "KRAS", Variant_Classification = "Missense_Mutation", HGVSp_Short = "p.G12C", Tumor_Sample_Barcode = smp[kras == 1 & k]),
               data.frame(Hugo_Symbol = "EGFR", Variant_Classification = "Missense_Mutation", HGVSp_Short = "p.L858R", Tumor_Sample_Barcode = smp[k & runif(n) < .04]))
  reg <- rbind(reg1, reg2); reg <- reg[reg$record_id %in% kid, ]
  list(tables = list(
    pt_char = data.frame(record_id = id, institution = site, birth_year = round(2016 - age), naaccr_sex_code = ifelse(female == 1, "Female", "Male"),
      naaccr_race_code_primary = sample(c("White", "Black", "Chinese", "Unknown"), n, TRUE, c(.75, .1, .08, .07)), stringsAsFactors = FALSE)[k, ],
    ca_dx_index = data.frame(record_id = id, ca_seq = 0, age_dx = age - t0 / 365.25, dob_ca_dx_days = round((age - t0 / 365.25) * 365.25),
      stage_dx_iv = ifelse(de_novo == 1, "Stage IV", "Stage I-III"), dmets_stage_i_iii = ifelse(de_novo == 1, NA, 1), dx_to_dmets_days = dmets,
      ca_hist_adeno_squamous = ifelse(adeno == 1, "Adenocarcinoma", "Squamous cell"),
      ca_lung_cigarette = sample(c("Never used", "Former user (quit >1 year)", "Current user"), n, TRUE, c(.15, .6, .25)),
      dmets_brain = brain, dx_to_dmets_brain_days = ifelse(brain == 1, adv, NA), dmets_liver = liver, dx_to_dmets_liver_days = ifelse(liver == 1, adv, NA),
      tt_os_dx_days = round(t0 + os1 * dpm), os_dx_status = st1, stringsAsFactors = FALSE)[k, ],
    ca_drugs = reg,
    prissmm_pathology = rbind(data.frame(record_id = id, dx_path_proc_days = t0 - sample(5:40, n, TRUE), pdl1_perc = pdl1)[k, ],
                              data.frame(record_id = id, dx_path_proc_days = t0 + 300, pdl1_perc = 99)[k & runif(n) < .2, ]),   # post-baseline: must be ignored
    cpt = data.frame(record_id = id, ca_seq = 0, cpt_genie_sample_id = smp,
      cpt_seq_assay_id = unname(c(MSK = "MSK-IMPACT468", DFCI = "DFCI-ONCOPANEL-3", VICC = "VICC-01-T7", UHN = "UHN-48-V1")[site]),
      dx_cpt_rep_days = round(t0 + rep_rel), dx_path_proc_cpt_days = round(t0 + rep_rel) - 20, stringsAsFactors = FALSE)[k, ],
    mutations_extended = mut), truth = truth)
}
write_synth <- function(tabs, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- c(pt_char = "patient_level_dataset.csv", ca_dx_index = "cancer_level_dataset_index.csv", ca_drugs = "regimen_cancer_level_dataset.csv",
         prissmm_pathology = "pathology_report_level_dataset.csv", cpt = "cancer_panel_test_level_dataset.csv", mutations_extended = "data_mutations_extended.txt")
  for (nm in names(tabs)) utils::write.table(tabs[[nm]], file.path(dir, f[[nm]]), sep = if (grepl("txt$", f[[nm]])) "\t" else ",", row.names = FALSE, na = "", qmethod = "double")
}
