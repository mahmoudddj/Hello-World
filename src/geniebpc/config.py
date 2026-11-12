"""Study constants. Column names follow the GENIE BPC NSCLC v2.0-public layout as exposed by the
public genieBPC R package test data (verified 2026-11-12; see docs/data_inventory.md). They must be
re-verified against the real release at gate G1 (src/geniebpc/feasibility.py does this)."""
from __future__ import annotations

DATA_ENV = "GENIE_BPC_DATA_DIR"
DAYS_PER_MONTH = 30.4375

# logical table -> filename glob candidates (first match wins)
TABLE_GLOBS = {
    "pt_char": ["*patient_level_dataset*.csv", "*pt_char*.csv"],
    "ca_dx_index": ["*cancer_level_dataset_index*.csv", "*ca_dx_index*.csv"],
    "ca_drugs": ["*regimen_cancer_level_dataset*.csv", "*ca_drugs*.csv"],
    "prissmm_pathology": ["*pathology_report_level_dataset*.csv", "*prissmm_pathology*.csv"],
    "cpt": ["*cancer_panel_test_level_dataset*.csv", "*cpt*.csv"],
    "mutations_extended": ["*mutations_extended*.txt", "*mutations_extended*.csv"],
    "cna": ["*CNA*.txt", "*cna*.txt", "*cna*.csv"],
    "fusions": ["*fusions*.txt", "*data_sv*.txt", "*fusions*.csv"],
}
OPTIONAL_TABLES = {"cna", "fusions"}

REQUIRED_COLUMNS = {
    "pt_char": ["record_id", "institution", "birth_year", "naaccr_sex_code", "naaccr_race_code_primary"],
    "ca_dx_index": ["record_id", "ca_seq", "age_dx", "stage_dx_iv", "dmets_stage_i_iii", "dx_to_dmets_days",
                    "ca_hist_adeno_squamous", "ca_lung_cigarette", "dmets_brain", "dx_to_dmets_brain_days",
                    "dmets_liver", "dx_to_dmets_liver_days", "dob_ca_dx_days"],
    "ca_drugs": ["record_id", "ca_seq", "regimen_number", "regimen_drugs", "drugs_ct_yn", "dx_reg_start_int",
                 "os_g_status", "tt_os_g_days"],
    "prissmm_pathology": ["record_id", "dx_path_proc_days", "pdl1_perc"],
    "cpt": ["record_id", "ca_seq", "cpt_genie_sample_id", "cpt_seq_assay_id", "dx_cpt_rep_days"],
    "mutations_extended": ["Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Tumor_Sample_Barcode"],
}

# --- treatment classification (matched on the leading drug name, case-insensitive) ---
PD1_PDL1 = ["pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "cemiplimab", "avelumab"]
PLATINUM = ["carboplatin", "cisplatin"]
PARTNER_CHEMO = ["pemetrexed", "paclitaxel", "nab-paclitaxel", "docetaxel", "gemcitabine", "etoposide",
                 "vinorelbine"]
ALLOWED_ADJUNCT = ["bevacizumab"]
ARM_MONO, ARM_CHEMO = "ICI_mono", "ICI_chemo"

# --- genomics ---
PROTEIN_ALTERING = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                    "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"}
LOSS_OF_FUNCTION = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                    "Translation_Start_Site"}
BIOMARKER_GENES = ["STK11", "KEAP1"]
# Panels assumed (UNVERIFIED until G1 checks genomic_information / panel gene lists) to tile both genes.
# Small hotspot panels (UHN-48/50, VICC-01-SOLIDTUMOR) are treated as NOT evaluable for KEAP1.
PANEL_PREFIX_BOTH_GENES = ["MSK-IMPACT", "DFCI-ONCOPANEL", "VICC-01-T5A", "VICC-01-T7", "UHN-555", "UHN-OCA"]
DRIVER_FUSION_GENES = ["ALK", "ROS1"]

# --- gate thresholds (pre-specified in protocol/SAP.md section 9) ---
MIN_N_PER_CELL = 20
MIN_EVENTS_PER_CELL = 10
TAUS_MONTHS = (12.0, 24.0)
