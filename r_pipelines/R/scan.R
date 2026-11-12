# Feasibility scan over many candidate questions. Counts only: no outcome model is fitted, nothing is estimated.
# Output: F_cells.csv (group sizes and event counts, small cells blanked) and F_verdicts.csv (one line per question).
gene_samples <- function(t, gene, classes = PROTEIN_ALTERING, hgvs = NULL) {
  m <- t$mutations_extended; m <- m[m$Hugo_Symbol == gene & m$Variant_Classification %in% classes, ]
  if (!is.null(hgvs)) m <- m[grepl(hgvs, m$HGVSp_Short), ]
  unique(m$Tumor_Sample_Barcode)
}
# One row per patient: genotype from all sequenced samples of the index cancer (timing ignored here; counts only).
patient_genotypes <- function(t) {
  s <- sample_alterations(t); mu <- t$mutations_extended
  flag <- function(samples) tapply(s$sample %in% samples, s$record_id, any)
  g <- data.frame(record_id = names(flag(character(0))), stringsAsFactors = FALSE)
  put <- function(nm, samples) g[[nm]] <<- as.integer(flag(samples)[g$record_id])
  g$large_panel <- as.integer(tapply(s$covers_both, s$record_id, any)[g$record_id])
  for (x in c("STK11", "KEAP1", "KRAS", "TP53")) g[[x]] <- as.integer(tapply(s[[paste0(x, "_alt")]], s$record_id, any)[g$record_id])
  put("SMARCA4", gene_samples(t, "SMARCA4")); put("KRAS_G12C", gene_samples(t, "KRAS", hgvs = "G12C"))
  put("KRAS_G12D", gene_samples(t, "KRAS", hgvs = "G12D")); put("KRAS_G12V", gene_samples(t, "KRAS", hgvs = "G12V"))
  put("EGFR_ex19del", gene_samples(t, "EGFR", "In_Frame_Del")); put("EGFR_L858R", gene_samples(t, "EGFR", hgvs = "L858R"))
  put("EGFR_other_sens", gene_samples(t, "EGFR", hgvs = "G719|L861Q|S768I"))
  put("BRAF_V600E", gene_samples(t, "BRAF", hgvs = "V600E")); put("ERBB2_ins", gene_samples(t, "ERBB2", "In_Frame_Ins"))
  put("MET_splice", gene_samples(t, "MET", c("Splice_Site", "Splice_Region")))
  fus <- function(gene) if (!is.null(t$fusions) && all(c("Hugo_Symbol", "Tumor_Sample_Barcode") %in% names(t$fusions))) t$fusions$Tumor_Sample_Barcode[t$fusions$Hugo_Symbol == gene] else character(0)
  for (x in c("ALK", "ROS1", "RET")) put(paste0(x, "_fusion"), fus(x))
  g$EGFR_sens <- as.integer(g$EGFR_ex19del | g$EGFR_L858R | g$EGFR_other_sens)
  g$driver <- as.integer(g$EGFR_sens | g$ALK_fusion | g$ROS1_fusion)
  g$SK <- as.integer(g$STK11 | g$KEAP1)
  mc <- table(mu$Tumor_Sample_Barcode[mu$Variant_Classification %in% PROTEIN_ALTERING])
  s$mut_count <- as.numeric(mc[s$sample]); s$mut_count[is.na(s$mut_count)] <- 0; big <- s[s$covers_both, ]
  g$mut_count <- as.numeric(tapply(big$mut_count, big$record_id, max)[g$record_id])
  g$panel_family <- as.character(tapply(sub("-?[0-9.]+[A-Z]?$", "", big$cpt_seq_assay_id), big$record_id, function(x) x[1])[g$record_id])
  g
}
band <- function(x) as.character(cut(x, c(-Inf, 9, 19, 49, 99, Inf), c("<10", "10-19", "20-49", "50-99", "100+")))
verdict <- function(n, ev) if (min(n) >= 50 && min(ev) >= 30) "promising" else if (min(n) >= 20 && min(ev) >= 10) "marginal" else "insufficient"

scan <- function(t, cfg) {
  g <- patient_genotypes(t); r <- regimens(t, cfg); r <- merge(r, g, by = "record_id")
  lo <- tolower(r$regimen_drugs); single <- !grepl(",", gsub("\\([^)]*\\)", "", r$regimen_drugs))
  r$egfr_tki <- grepl("osimertinib|erlotinib|gefitinib|afatinib|dacomitinib", lo) & single
  r$osimertinib <- grepl("osimertinib", lo)
  r$docetaxel_based <- grepl("docetaxel", lo) & !grepl(paste(c(PD1_PDL1, PLATINUM), collapse = "|"), lo) & !grepl("investigational", lo)
  r$pem_doublet <- r$class == "platinum_doublet" & grepl("pemetrexed", lo); r$bev <- grepl("bevacizumab", lo)
  r$platinum_agent <- ifelse(grepl("cisplatin", lo), "cisplatin", ifelse(grepl("carboplatin", lo), "carboplatin", NA))
  l1 <- r[r$line == 1, ]; r$prior_platinum <- r$record_id %in% l1$record_id[grepl("carboplatin|cisplatin", tolower(l1$regimen_drugs))]
  r <- add_covariates(r, t, cfg)
  nd <- r[r$large_panel == 1 & r$driver == 0, ]                       # non-driver, genotype evaluable on a broad panel
  chemo1 <- nd[nd$line == 1 & nd$class == "platinum_doublet", ]
  ici <- first_ici(nd)
  cells <- list(); verd <- list()
  add <- function(id, title, design, novelty, d, groups, event = "death", key = NULL) {
    # groups: character vector of column names in d defining the comparison cells; `key` = cells that must be adequate
    if (!nrow(d)) { verd[[id]] <<- data.frame(id = id, question = title, design = design, cells = 0, smallest_cell_n = "0", smallest_cell_events = "0", verdict = "insufficient", novelty_prior_unverified = novelty); return(invisible()) }
    lab <- do.call(paste, c(lapply(groups, function(v) paste0(v, "=", d[[v]])), sep = " | "))
    a <- data.frame(id = id, group = sort(unique(lab)), stringsAsFactors = FALSE)
    a$n <- as.numeric(table(lab)[a$group]); a$patients <- as.numeric(tapply(d$record_id, lab, function(x) length(unique(x)))[a$group])
    a$events <- as.numeric(tapply(d[[event]], lab, function(x) sum(x == 1, na.rm = TRUE))[a$group]); a$event_type <- event
    cells[[id]] <<- a
    k <- if (is.null(key)) a else a[grepl(key, a$group), ]; if (!nrow(k)) k <- a
    verd[[id]] <<- data.frame(id = id, question = title, design = design, cells = nrow(a), smallest_cell_n = band(min(k$n)), smallest_cell_events = band(min(k$events)),
                              verdict = if (nrow(k) < 2) "insufficient" else verdict(k$n, k$events), novelty_prior_unverified = novelty)
  }
  two <- function(d1, d2, n1, n2) { d1$tx <- n1; d2$tx <- n2; rbind(d1, d2) }
  cvi <- two(chemo1, ici, "platinum_doublet_1L", "ICI_mono_first_use")
  p0 <- nd[nd$line == 1 & nd$class %in% c("ICI_mono", "ICI_chemo"), ]
  add("Q01", "STK11/KEAP1 x first-line ICI monotherapy vs ICI + platinum chemotherapy (original question)", "treatment x biomarker interaction", "moderate; Elkrief 2023 overlaps", p0, c("class", "SK"))
  add("Q02", "STK11/KEAP1 x anti-PD-(L)1 monotherapy vs first-line platinum doublet (paper A)", "treatment x biomarker interaction", "moderate; prior null in Flatiron", cvi, c("tx", "SK"))
  d <- nd[nd$line == 1 & !grepl("trial", nd$class), ]; r0 <- regimens(t, cfg); ab <- add_biomarker(r0[r0$line == 1, ], t, cfg, "anytime"); d$ngs_after_start <- as.integer(d$record_id %in% ab$record_id[ab$entry_m > 0])
  add("Q03", "Impact of NGS timing / delayed entry on biomarker survival estimates (paper B, methods)", "methods", "moderate; extends Brown 2022", d, c("ngs_after_start", "SK"))
  d <- chemo1[chemo1$histology != "Squamous", ]; d$co_class <- ifelse(d$KRAS & d$SK, "KRAS+SK", ifelse(d$KRAS == 1, "KRAS_only", ifelse(d$SK == 1, "SK_only", "neither")))
  add("Q04", "KRAS / STK11-KEAP1 co-alteration classes on first-line platinum doublet (paper C)", "prognostic", "low-moderate; cohort paper partly covers", d, "co_class")
  d <- cvi[cvi$KRAS == 1, ]; d$allele <- ifelse(d$KRAS_G12C == 1, "G12C", ifelse(d$KRAS_G12D == 1, "G12D", ifelse(d$KRAS_G12V == 1, "G12V", "other")))
  add("Q05", "KRAS allele (G12C / G12D / G12V / other) and outcomes on chemotherapy and immunotherapy", "prognostic by allele", "moderate", d, c("tx", "allele"), key = "G12C|G12D|G12V")
  for (q in list(c("Q06", "TP53", "TP53 x immunotherapy vs chemotherapy"), c("Q07", "SMARCA4", "SMARCA4 x immunotherapy vs chemotherapy"), c("Q08", "KRAS", "KRAS mutation x immunotherapy vs chemotherapy")))
    add(q[1], paste(q[3], "(OS; is the prognostic effect treatment-specific?)"), "treatment x biomarker interaction", if (q[2] == "SMARCA4") "moderate-high if cells allow" else "low-moderate", cvi, c("tx", q[2]))
  for (q in list(c("Q09", "liver_mets", "Liver metastases at treatment start"), c("Q10", "brain_mets", "Brain metastases at treatment start"))) add(q[1], paste(q[3], "x immunotherapy vs chemotherapy"), "treatment x clinical-factor interaction", "low-moderate", cvi, c("tx", q[2]))
  d <- cvi; d$never_smoker <- as.integer(d$smoking == "Never"); add("Q11", "Never-smoking status x immunotherapy vs chemotherapy (non-driver NSCLC)", "treatment x clinical-factor interaction", "moderate", d, c("tx", "never_smoker"))
  d <- cvi; d$age75 <- as.integer(d$age >= 75); add("Q12", "Age >= 75 x immunotherapy vs chemotherapy", "treatment x clinical-factor interaction", "low-moderate", d, c("tx", "age75"))
  add("Q13", "Sex x immunotherapy vs chemotherapy", "treatment x clinical-factor interaction", "low", cvi, c("tx", "female"))
  d <- ici; add("Q14", "PD-L1 level and outcomes on anti-PD-(L)1 monotherapy (how much PD-L1 is recorded?)", "prognostic / data completeness", "low", d, "pdl1_cat")
  d <- ici[!is.na(ici$mut_count), ]; if (nrow(d)) d$tmb_proxy_tertile <- as.integer(stats::ave(d$mut_count, d$panel_family, FUN = function(x) cut(rank(x, ties.method = "first"), 3, labels = FALSE)))
  add("Q15", "Panel mutation count (TMB proxy, tertiles within panel) and outcomes on anti-PD-(L)1 monotherapy", "prognostic", "low; well studied", d, "tmb_proxy_tertile")
  d <- nd[nd$line == 2 & nd$prior_platinum & (nd$class == "ICI_mono" | nd$docetaxel_based), ]; d$second_line <- ifelse(d$class == "ICI_mono", "ICI_mono", "docetaxel_based")
  add("Q16", "Second-line anti-PD-(L)1 vs docetaxel after platinum (target trial emulation; benchmark against RCTs)", "treatment comparison", "moderate as RCT-benchmarking methods paper", d, "second_line")
  add("Q17", "STK11/KEAP1 x second-line anti-PD-(L)1 vs docetaxel", "treatment x biomarker interaction", "moderate-high if cells allow", d, c("second_line", "SK"))
  d <- nd[nd$line == 1 & nd$pem_doublet & nd$histology != "Squamous", ]; add("Q18", "Bevacizumab added to first-line platinum-pemetrexed, non-squamous", "treatment comparison", "low-moderate", d, "bev")
  d <- nd[nd$line == 1 & nd$class == "platinum_doublet" & !is.na(nd$platinum_agent), ]; add("Q19", "Carboplatin vs cisplatin first-line doublets", "treatment comparison", "low", d, "platinum_agent")
  e <- r[r$EGFR_sens == 1 & r$egfr_tki, ]; e <- e[order(e$record_id, e$t0), ]; e <- e[!duplicated(e$record_id), ]
  add("Q20", "TP53 co-mutation and PFS on first EGFR TKI in EGFR-mutant NSCLC", "prognostic", "low-moderate; curated PFS is the added value", e, "TP53", event = "prog")
  e$egfr_type <- ifelse(e$EGFR_ex19del == 1, "ex19del", ifelse(e$EGFR_L858R == 1, "L858R", "other")); add("Q21", "EGFR exon 19 deletion vs L858R on first EGFR TKI", "prognostic", "low", e, "egfr_type", event = "prog", key = "ex19del|L858R")
  add("Q22", "Osimertinib vs earlier-generation EGFR TKI as first TKI", "treatment comparison", "low", e, "osimertinib", event = "prog")
  # ---- patient-level questions
  dx <- merge(index_dx(t), g, by = "record_id"); dx$death <- num(dx$os_dx_status)
  dx$dmets_event <- as.integer(!is.na(num(dx$dmets_stage_i_iii)) & num(dx$dmets_stage_i_iii) == 1)
  dx$brain_event <- as.integer(!is.na(num(dx$dmets_brain)) & num(dx$dmets_brain) == 1)
  early <- dx[dx$de_novo == 0 & dx$large_panel == 1, ]
  for (q in list(c("Q23", "KRAS"), c("Q24", "SK"), c("Q25", "TP53"), c("Q26", "EGFR_sens")))
    add(q[1], paste0("Stage I-III at diagnosis: ", q[2], " and time to distant metastasis (delayed entry needed)"), "prognostic, early stage", "moderate", early, q[2], event = "dmets_event")
  dx$geno <- ifelse(dx$EGFR_sens == 1, "EGFR", ifelse(dx$ALK_fusion | dx$ROS1_fusion | dx$RET_fusion, "fusion", ifelse(dx$KRAS == 1 & dx$SK == 1, "KRAS+SK", ifelse(dx$KRAS == 1, "KRAS", "other"))))
  add("Q27", "Genotype and brain metastasis (cumulative incidence by driver / co-mutation class)", "prognostic, organotropism", "moderate", dx, "geno", event = "brain_event", key = "EGFR|KRAS")
  add("Q28", "Genotype and brain metastasis in stage I-III at diagnosis", "prognostic, organotropism", "moderate", dx[dx$de_novo == 0, ], "geno", event = "brain_event", key = "EGFR|KRAS")
  rare <- do.call(rbind, lapply(c("BRAF_V600E", "ERBB2_ins", "MET_splice", "ALK_fusion", "ROS1_fusion", "RET_fusion", "KRAS_G12C", "SMARCA4"), function(v) { z <- dx[dx[[v]] == 1 & !is.na(dx[[v]]), ]; if (nrow(z)) { z$alteration <- v; z } }))
  if (!is.null(rare)) add("Q29", "Rare driver subsets: outcomes and treatment patterns (descriptive)", "descriptive", "low unless a subset is large", rare, "alteration")
  pt <- merge(dx, t$pt_char[c("record_id", "naaccr_race_code_primary", "naaccr_sex_code")], by = "record_id")
  pt$race <- ifelse(pt$naaccr_race_code_primary %in% "White", "White", ifelse(pt$naaccr_race_code_primary %in% "Black", "Black", ifelse(grepl("Chinese|Asian|Japanese|Korean|Vietnamese|Filipino", pt$naaccr_race_code_primary), "Asian", "Other_unknown")))
  add("Q30", "Race and somatic alteration prevalence / outcomes", "descriptive / disparities", "moderate if Black and Asian groups are large", pt[pt$large_panel == 1, ], "race", key = "Black|Asian")
  write_out(do.call(rbind, cells), "F_cells", cfg, c("n", "patients", "events"))
  v <- do.call(rbind, verd); v <- v[order(match(v$verdict, c("promising", "marginal", "insufficient")), v$id), ]
  utils::write.csv(v, file.path(cfg$out_dir, "F_verdicts.csv"), row.names = FALSE)
  invisible(v)
}
