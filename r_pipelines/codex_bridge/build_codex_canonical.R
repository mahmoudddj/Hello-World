#!/usr/bin/env Rscript
# Builds the five "canonical" input tables that the Codex pipeline (GENIE_BPC_R_Pipeline v1.0.0) expects,
# directly from the GENIE BPC NSCLC v2.0-public release files, so its unconfigured adapter is not needed.
#   Rscript codex_bridge/build_codex_canonical.R <GENIE_SECURE_ROOT>        (run from the r_pipelines folder)
# Reads GENIE_BPC_DATA_DIR; writes <root>/canonical/*.csv and <root>/source_mapping.csv (verified = NO until YOU check it).
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1); root <- args[1]
source("config.R"); for (f in c("load", "derive")) source(file.path("R", paste0(f, ".R")))
L <- load_tables(cfg$data_dir); stopifnot(schema_ok(L$schema)); t <- L$tables
dx <- index_dx(t); dx$cancer_id <- as.character(dx$ca_seq)
r <- regimens(t, cfg)
drug_cols <- grep("^drugs_drug_[0-9]+$", names(r), value = TRUE)
comp <- function(i) {            # component flags and earliest start day per component, from per-drug columns when present
  if (length(drug_cols)) { nm <- unlist(r[i, drug_cols]); dy <- num(unlist(r[i, sub("drugs_drug_", "dx_drug_start_int_", drug_cols)])) }
  else { nm <- strsplit(gsub("\\([^)]*\\)", "", r$regimen_drugs[i]), ",")[[1]]; dy <- rep(r$t0[i], length(nm)) }
  ok <- !is.na(nm) & nzchar(trimws(nm)); nm <- tolower(trimws(gsub("\\(.*$", "", nm[ok]))); dy <- dy[ok]
  is <- function(set, anywhere = FALSE) vapply(nm, function(z) any(if (anywhere) vapply(set, grepl, logical(1), x = z, fixed = TRUE) else startsWith(z, set)), logical(1))
  ici <- is(PD1_PDL1); pl <- is(PLATINUM); cy <- is(PARTNER, TRUE) & !pl; ot <- !(ici | pl | cy)
  day <- function(k) if (any(k) && any(!is.na(dy[k]))) min(dy[k], na.rm = TRUE) else NA
  c(ici = any(ici), platinum = any(pl), cytotoxic = any(cy), other_systemic = any(ot) || identical(tolower(as.character(r$drugs_ct_yn[i])), "yes"),
    ici_day = day(ici), platinum_day = day(pl), cytotoxic_day = day(cy))
}
cc <- as.data.frame(t(vapply(seq_len(nrow(r)), comp, numeric(7))))
regs <- data.frame(id = r$record_id, cancer_id = as.character(r$ca_seq), regimen_id = paste(r$record_id, r$regimen_number, sep = "_R"),
                   line_advanced = r$line, index_day = r$t0, cc, stringsAsFactors = FALSE)
first <- r[r$line == 1, c("record_id", "t0")]
p <- merge(dx, t$pt_char[c("record_id", "institution", "birth_year", "naaccr_sex_code")], by = "record_id", suffixes = c(".dx", ""))
p <- merge(p, first, by = "record_id", all.x = TRUE)
base <- ifelse(is.na(p$t0), 0, p$t0)
met <- function(site) as.integer(!is.na(num(p[[paste0("dmets_", site)]])) & num(p[[paste0("dmets_", site)]]) == 1 & !is.na(num(p[[paste0("dx_to_dmets_", site, "_days")]])) & num(p[[paste0("dx_to_dmets_", site, "_days")]]) <= base)
h <- as.character(p$ca_hist_adeno_squamous); s <- as.character(p$ca_lung_cigarette); dead <- num(p$os_dx_status)
patients <- data.frame(id = p$record_id, cancer_id = p$cancer_id, nsclc = 1, advanced_day = p$adv_day, age = num(p$age_dx) + base / 365.25,
  sex = ifelse(tolower(p$naaccr_sex_code) == "female", "female", ifelse(tolower(p$naaccr_sex_code) == "male", "male", "other")),
  ecog = NA,                                         # not available in the public release
  smoking = ifelse(is.na(s), NA, ifelse(startsWith(s, "Never"), "never", ifelse(startsWith(s, "Current"), "current", ifelse(startsWith(s, "Former"), "former", NA)))),
  histology = ifelse(is.na(h), NA, ifelse(startsWith(h, "Adeno"), "adenocarcinoma", ifelse(startsWith(h, "Squamous"), "squamous", "other_nsclc"))),
  brain_met = met("brain"), liver_met = met("liver"), year = floor(num(p$birth_year) + (num(p$dob_ca_dx_days) + base) / 365.25),
  site = p$institution, covariate_day = p$t0, death_day = ifelse(dead == 1, num(p$tt_os_dx_days), NA),
  last_alive_day = ifelse(dead == 1, NA, num(p$tt_os_dx_days)), dead = dead, stringsAsFactors = FALSE)
sa <- sample_alterations(t); sa <- merge(sa, dx[c("record_id", "ca_seq")], by = c("record_id", "ca_seq"))
sday <- if ("dx_path_proc_cpt_days" %in% names(t$cpt)) num(t$cpt$dx_path_proc_cpt_days[match(sa$sample, t$cpt$cpt_genie_sample_id)]) else NA
sday[!is.na(sday) & !is.na(sa$rep_day) & sday > sa$rep_day] <- NA
assays <- data.frame(id = sa$record_id, cancer_id = as.character(sa$ca_seq), assay_id = sa$sample, sample_day = sday, report_day = sa$rep_day,
  stk11_assessed = as.integer(sa$covers_both), keap1_assessed = as.integer(sa$covers_both), kras_assessed = 1,
  drivers_assessed = as.integer(sa$covers_both), driver_positive = as.integer(sa$driver), stringsAsFactors = FALSE)
assays <- assays[!duplicated(assays$assay_id), ]
variants <- do.call(rbind, lapply(c("STK11", "KEAP1", "KRAS"), function(g) { a <- sa$sample[sa[[paste0(g, "_alt")]]]
  data.frame(assay_id = a, gene = rep(g, length(a)), qualifying = rep(1, length(a)), stringsAsFactors = FALSE) }))
variants <- variants[variants$assay_id %in% assays$assay_id, ]
pl <- t$prissmm_pathology; pl <- pl[pl$record_id %in% patients$id & !is.na(num(pl$pdl1_perc)) & !is.na(num(pl$dx_path_proc_days)), ]
pdl1 <- data.frame(id = pl$record_id, cancer_id = patients$cancer_id[match(pl$record_id, patients$id)], test_id = paste0("PDL1_", seq_len(nrow(pl))),
  sample_day = num(pl$dx_path_proc_days), report_day = num(pl$dx_path_proc_days), tps = pmin(100, pmax(0, num(pl$pdl1_perc))), stringsAsFactors = FALSE)
out <- file.path(root, "canonical"); if (dir.exists(out)) stop("canonical/ already exists: move it aside before rebuilding", call. = FALSE)
dir.create(out, recursive = TRUE)
tabs <- list(patients = patients, regimens = regs[regs$id %in% patients$id, ], assays = assays[assays$id %in% patients$id, ], variants = variants, pdl1 = pdl1)
for (nm in names(tabs)) utils::write.csv(tabs[[nm]], file.path(out, paste0(nm, ".csv")), row.names = FALSE, na = "")
src <- list(patients = "cancer_level_dataset_index + patient_level_dataset", regimens = "regimen_cancer_level_dataset",
            assays = "cancer_panel_test_level_dataset + data_mutations_extended + data_fusions", variants = "data_mutations_extended + data_CNA", pdl1 = "pathology_report_level_dataset")
map <- do.call(rbind, lapply(names(tabs), function(nm) data.frame(canonical_table = nm, canonical_column = names(tabs[[nm]]), source_table = src[[nm]],
  source_column = "see codex_bridge/build_codex_canonical.R", derivation = "see codex_bridge/build_codex_canonical.R", time_origin = "days from index NSCLC diagnosis",
  release = "SET_EXACT_AUTHORIZED_RELEASE", verified = "NO")))
utils::write.csv(map, file.path(root, "source_mapping.csv"), row.names = FALSE)
cat("Canonical tables written. Patients:", nrow(patients), " first-line regimens:", sum(regs$line_advanced == 1), "\n",
    "Review source_mapping.csv, set release and verified = YES yourself, and remove 'ecog' from covariates in the Codex config.\n")
