# Source-to-derived-variable mapping. All day counts are relative to diagnosis of the index NSCLC.
# One row per REGIMEN START of interest (a patient can contribute more than one row; `record_id` clusters them).
PD1_PDL1 <- c("pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "cemiplimab", "avelumab")
PLATINUM <- c("carboplatin", "cisplatin")
PARTNER  <- c("pemetrexed", "paclitaxel", "docetaxel", "gemcitabine", "etoposide", "vinorelbine")
ADJUNCT  <- c("bevacizumab")
PROTEIN_ALTERING <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                      "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation")
LOSS_OF_FUNCTION <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site")
# Panels ASSUMED to tile STK11 and KEAP1. feasibility() prints mutation counts per panel so this can be checked.
PANELS_BOTH <- c("MSK-IMPACT", "DFCI-ONCOPANEL", "VICC-01-T5A", "VICC-01-T7")
num <- function(x) suppressWarnings(as.numeric(x))

classify_regimen <- function(regimen_drugs, trial = "No") {
  one <- function(s, tr) {
    if (is.na(s) || !nzchar(trimws(s))) return("other:missing")
    d <- tolower(trimws(strsplit(gsub("\\([^)]*\\)", "", s), ",")[[1]])); d <- d[nzchar(d)]
    if (identical(tolower(trimws(as.character(tr))), "yes") || any(grepl("investigational", d))) return("other:trial_or_investigational")
    has <- function(x, names, anywhere = FALSE) vapply(x, function(z) any(if (anywhere) vapply(names, grepl, logical(1), x = z, fixed = TRUE) else startsWith(z, names)), logical(1))
    ici <- has(d, PD1_PDL1); plat <- has(d, PLATINUM); part <- has(d, PARTNER, TRUE) & !plat; adj <- has(d, ADJUNCT)
    rest <- !(ici | plat | part | adj)
    if (!any(ici)) {
      if (any(plat) && any(part) && !any(rest)) return("platinum_doublet")
      return("other:no_pd1_pdl1")
    }
    if (sum(ici) > 1 || any(rest)) return("other:ici_plus_nonprotocol_agent")
    if (!any(plat) && !any(part) && !any(adj)) return("ICI_mono")
    if (any(plat)) return("ICI_chemo")
    "other:ici_nonplatinum_combo"
  }
  mapply(one, regimen_drugs, rep_len(trial, length(regimen_drugs)), USE.NAMES = FALSE)
}

index_dx <- function(t) {
  dx <- t$ca_dx_index[order(t$ca_dx_index$record_id, num(t$ca_dx_index$ca_seq)), ]
  dx <- dx[!duplicated(dx$record_id), ]
  dx$de_novo <- as.integer(!is.na(dx$stage_dx_iv) & dx$stage_dx_iv == "Stage IV")
  recur <- dx$de_novo == 0 & !is.na(num(dx$dmets_stage_i_iii)) & num(dx$dmets_stage_i_iii) == 1
  dx$adv_day <- ifelse(dx$de_novo == 1, 0, ifelse(recur, num(dx$dx_to_dmets_days), NA))
  dx
}

# All regimens of the index cancer started on/after the advanced-disease date, with line number and class.
regimens <- function(t, cfg) {
  dx <- index_dx(t); dx <- dx[!is.na(dx$adv_day), ]
  r <- merge(t$ca_drugs, dx[c("record_id", "ca_seq", "adv_day", "de_novo")], by = c("record_id", "ca_seq"))
  r$t0 <- num(r$dx_reg_start_int)
  r <- r[!is.na(r$t0), ]
  r$prior_systemic <- as.integer(r$record_id %in% r$record_id[r$t0 < r$adv_day])
  r <- r[r$t0 >= r$adv_day, ]
  r <- r[order(r$record_id, r$t0, num(r$regimen_number)), ]
  r$line <- stats::ave(seq_len(nrow(r)), r$record_id, FUN = seq_along)
  r$class <- classify_regimen(r$regimen_drugs, r$drugs_ct_yn)
  r$os_m <- num(r$tt_os_g_days) / cfg$days_per_month; r$death <- num(r$os_g_status)
  r$pfs_m <- num(r$tt_pfs_i_or_m_g_days) / cfg$days_per_month; r$prog <- num(r$pfs_i_or_m_g_status)
  r
}

sample_alterations <- function(t) {
  s <- t$cpt[c("record_id", "ca_seq", "cpt_genie_sample_id", "cpt_seq_assay_id", "dx_cpt_rep_days")]
  names(s)[3] <- "sample"; s$rep_day <- num(s$dx_cpt_rep_days)
  m <- t$mutations_extended; m <- m[m$Variant_Classification %in% PROTEIN_ALTERING, ]
  deep <- list()
  if (!is.null(t$cna) && "Hugo_Symbol" %in% names(t$cna)) for (g in c("STK11", "KEAP1")) {
    row <- t$cna[t$cna$Hugo_Symbol == g, setdiff(names(t$cna), "Hugo_Symbol"), drop = FALSE]
    if (nrow(row)) { v <- num(unlist(row[1, ])); deep[[g]] <- names(row)[!is.na(v) & v <= -2] }
  }
  for (g in c("STK11", "KEAP1", "KRAS", "TP53")) {
    mg <- m[m$Hugo_Symbol == g, ]
    s[[paste0(g, "_alt")]] <- s$sample %in% c(mg$Tumor_Sample_Barcode, deep[[g]])
    s[[paste0(g, "_lof")]] <- s$sample %in% c(mg$Tumor_Sample_Barcode[mg$Variant_Classification %in% LOSS_OF_FUNCTION], deep[[g]])
  }
  e <- m[m$Hugo_Symbol == "EGFR", ]
  sens <- e[grepl("L858R|G719|L861Q|S768I", e$HGVSp_Short) | e$Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins"), ]
  drv <- sens$Tumor_Sample_Barcode
  if (!is.null(t$fusions) && all(c("Hugo_Symbol", "Tumor_Sample_Barcode") %in% names(t$fusions)))
    drv <- c(drv, t$fusions$Tumor_Sample_Barcode[t$fusions$Hugo_Symbol %in% c("ALK", "ROS1")])
  s$driver <- s$sample %in% drv
  s$covers_both <- Reduce(`|`, lapply(PANELS_BOTH, function(p) startsWith(toupper(s$cpt_seq_assay_id), p)))
  s$covers_both[is.na(s$covers_both)] <- FALSE          # unknown assay: cannot assert coverage
  s
}

# Attach biomarker status to regimen rows.
#   ngs = "pretreatment": only reports dated <= t0 (no immortal time; entry = 0)
#   ngs = "anytime":      any evaluable report; entry_m = months from t0 to first evaluable report (delayed entry)
#   ngs = "naive":        any evaluable report, entry forced to 0 (the biased approach, used ONLY in paper B)
add_biomarker <- function(r, t, cfg, ngs = c("anytime", "pretreatment", "naive"), lof_only = FALSE) {
  ngs <- match.arg(ngs); s <- sample_alterations(t)
  x <- merge(r[c("record_id", "ca_seq", "regimen_number", "t0")], s, by = c("record_id", "ca_seq"))
  x$rel <- x$rep_day - x$t0
  drv <- unique(x$record_id[x$driver])
  x <- x[x$covers_both & !is.na(x$rel), ]
  if (ngs == "pretreatment") x <- x[x$rel <= 0, ]
  key <- paste(x$record_id, x$regimen_number)
  k <- if (lof_only) "_lof" else "_alt"
  agg <- data.frame(key = unique(key), stringsAsFactors = FALSE)
  for (g in c("STK11", "KEAP1", "KRAS", "TP53")) agg[[g]] <- as.integer(tapply(x[[paste0(g, if (g %in% c("STK11", "KEAP1")) k else "_alt")]], key, any)[agg$key])
  agg$first_rel <- as.numeric(tapply(x$rel, key, min)[agg$key])
  r$key <- paste(r$record_id, r$regimen_number)
  d <- merge(r, agg, by = "key")
  d <- d[!d$record_id %in% drv, ]
  d$biomarker <- as.integer(d$STK11 == 1 | d$KEAP1 == 1)
  d$ngs_rel_m <- d$first_rel / cfg$days_per_month
  d$entry_m <- if (ngs == "anytime") pmax(0, d$ngs_rel_m) else 0
  d$ngs_mode <- ngs
  d
}

add_covariates <- function(d, t, cfg) {
  dx <- index_dx(t)
  d <- merge(d, dx[c("record_id", "age_dx", "dob_ca_dx_days", "ca_hist_adeno_squamous", "ca_lung_cigarette", "dmets_brain",
                     "dx_to_dmets_brain_days", "dmets_liver", "dx_to_dmets_liver_days")], by = "record_id")
  p <- t$pt_char[c("record_id", "institution", "birth_year", "naaccr_sex_code", "naaccr_race_code_primary")]
  names(p)[2] <- "site"
  d$institution <- NULL
  d <- merge(d, p, by = "record_id")
  d$age <- num(d$age_dx) + d$t0 / 365.25
  d$female <- as.integer(tolower(d$naaccr_sex_code) == "female")
  d$race <- ifelse(d$naaccr_race_code_primary %in% "White", "White", ifelse(d$naaccr_race_code_primary %in% "Black", "Black", "Other_unknown"))
  h <- as.character(d$ca_hist_adeno_squamous)
  d$histology <- ifelse(startsWith(h, "Adeno") & !is.na(h), "Adenocarcinoma", ifelse(startsWith(h, "Squamous") & !is.na(h), "Squamous", "Other_NOS"))
  s <- as.character(d$ca_lung_cigarette)
  d$smoking <- ifelse(is.na(s), "Unknown", ifelse(startsWith(s, "Never"), "Never", ifelse(startsWith(s, "Current"), "Current", ifelse(startsWith(s, "Former"), "Former", "Unknown"))))
  for (site in c("brain", "liver")) {
    day <- num(d[[paste0("dx_to_dmets_", site, "_days")]])
    d[[paste0(site, "_mets")]] <- as.integer(!is.na(num(d[[paste0("dmets_", site)]])) & num(d[[paste0("dmets_", site)]]) == 1 & !is.na(day) & day <= d$t0)
  }
  d$year_t0 <- floor(num(d$birth_year) + (num(d$dob_ca_dx_days) + d$t0) / 365.25)   # approximate (birth year only)
  pl <- t$prissmm_pathology[c("record_id", "dx_path_proc_days", "pdl1_perc")]
  pl$pdl1_perc <- num(pl$pdl1_perc); pl$day <- num(pl$dx_path_proc_days); pl <- pl[!is.na(pl$pdl1_perc) & !is.na(pl$day), ]
  d$pdl1_cat <- "Unknown"; d$pdl1_day <- NA_real_
  if (nrow(pl)) {
    by_id <- split(pl, pl$record_id)
    for (i in which(d$record_id %in% names(by_id))) {
      z <- by_id[[d$record_id[i]]]; z <- z[z$day <= d$t0[i], ]
      if (nrow(z)) { j <- which.max(z$day); v <- z$pdl1_perc[j]; d$pdl1_day[i] <- z$day[j]
        d$pdl1_cat[i] <- if (v < 1) "<1%" else if (v < 50) "1-49%" else ">=50%" }
    }
  }
  d
}
COVARS <- c("age", "female", "race", "site", "histology", "smoking", "brain_mets", "liver_mets", "de_novo",
            "prior_systemic", "year_t0", "pdl1_cat", "TP53")
