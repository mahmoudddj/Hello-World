# Local settings. Edit here or set the environment variables. Nothing in this file is a credential.
cfg <- list(
  data_dir   = Sys.getenv("GENIE_BPC_DATA_DIR", unset = ""),   # folder holding the NSCLC v2.0-public files
  out_dir    = Sys.getenv("GENIE_BPC_OUT_DIR", unset = "results"),
  min_cell   = 10L,        # counts 1..9 are suppressed in every exported table
  boot       = 1000L,      # patient-level bootstrap replicates
  seed       = 20261112L,
  taus_os    = c(12, 24),  # months
  taus_pfs   = c(6, 12),
  min_n_cell = 20L,        # analysis stop rule: patients per comparison cell
  min_ev_cell = 10L,       # analysis stop rule: events per comparison cell
  days_per_month = 30.4375
)
