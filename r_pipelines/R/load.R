# Read GENIE BPC NSCLC release files (v2.0-public layout). Never writes patient-level data.
TABLE_GLOBS <- list(
  pt_char = c("patient_level_dataset", "pt_char"),
  ca_dx_index = c("cancer_level_dataset_index", "ca_dx_index"),
  ca_drugs = c("regimen_cancer_level_dataset", "ca_drugs"),
  prissmm_pathology = c("pathology_report_level_dataset", "prissmm_pathology"),
  cpt = c("cancer_panel_test_level_dataset", "cpt"),
  mutations_extended = c("mutations_extended"),
  cna = c("data_CNA", "data_cna", "cna"),
  fusions = c("data_fusions", "fusions"))
OPTIONAL_TABLES <- c("cna", "fusions")
REQUIRED_COLUMNS <- list(
  pt_char = c("record_id", "institution", "birth_year", "naaccr_sex_code", "naaccr_race_code_primary"),
  ca_dx_index = c("record_id", "ca_seq", "age_dx", "dob_ca_dx_days", "stage_dx_iv", "dmets_stage_i_iii",
    "dx_to_dmets_days", "ca_hist_adeno_squamous", "ca_lung_cigarette", "dmets_brain", "dx_to_dmets_brain_days",
    "dmets_liver", "dx_to_dmets_liver_days", "tt_os_dx_days", "os_dx_status"),
  ca_drugs = c("record_id", "ca_seq", "regimen_number", "regimen_drugs", "drugs_ct_yn", "dx_reg_start_int",
    "os_g_status", "tt_os_g_days", "pfs_i_or_m_g_status", "tt_pfs_i_or_m_g_days"),
  prissmm_pathology = c("record_id", "dx_path_proc_days", "pdl1_perc"),
  cpt = c("record_id", "ca_seq", "cpt_genie_sample_id", "cpt_seq_assay_id", "dx_cpt_rep_days"),
  mutations_extended = c("Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Tumor_Sample_Barcode"))

find_file <- function(dir, keys) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  f <- f[grepl("\\.(csv|txt|tsv)$", f, ignore.case = TRUE) & !grepl("non_index", basename(f))]
  for (k in keys) { hit <- f[grepl(k, basename(f), fixed = TRUE)]; if (length(hit)) return(sort(hit)[1]) }
  NA_character_
}
read_table <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = if (sep == "\t") "#" else "",
                    check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"), fill = TRUE)
}
load_tables <- function(dir) {
  if (!nzchar(dir) || !dir.exists(dir)) stop("Set GENIE_BPC_DATA_DIR to the folder with the release files.", call. = FALSE)
  tabs <- list(); rep <- list()
  for (nm in names(TABLE_GLOBS)) {
    f <- find_file(dir, TABLE_GLOBS[[nm]])
    opt <- nm %in% OPTIONAL_TABLES
    if (is.na(f)) { rep[[length(rep) + 1]] <- data.frame(table = nm, file = "", column = "<table>", present = FALSE, optional = opt); next }
    tabs[[nm]] <- read_table(f)
    rep[[length(rep) + 1]] <- data.frame(table = nm, file = basename(f), column = "<table>", present = TRUE, optional = opt)
    for (cc in REQUIRED_COLUMNS[[nm]])
      rep[[length(rep) + 1]] <- data.frame(table = nm, file = basename(f), column = cc, present = cc %in% names(tabs[[nm]]), optional = FALSE)
  }
  list(tables = tabs, schema = do.call(rbind, rep))
}
schema_ok <- function(schema) all(schema$present[!schema$optional])
