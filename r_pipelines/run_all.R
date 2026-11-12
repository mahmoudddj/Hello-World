#!/usr/bin/env Rscript
# Usage (from this folder):  Rscript run_all.R test | feasibility | A | B | C | all
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "help"
source("config.R"); for (f in c("load", "derive", "estimators", "privacy", "papers", "scan", "synth")) source(file.path("R", paste0(f, ".R")))
if (mode == "help") { cat("Usage: Rscript run_all.R test | feasibility | A | B | C | all\n"); quit(status = 0) }
if (mode == "test") { source("tests/run_tests.R"); quit(status = 0) }
L <- load_tables(cfg$data_dir)
if (!schema_ok(L$schema)) { write_out(L$schema, "00_schema_check", cfg, character(0)); stop("Required tables/columns missing: see 00_schema_check.csv", call. = FALSE) }
t <- L$tables
if (mode %in% c("feasibility", "all")) { feasibility(t, L$schema, cfg); scan(t, cfg); cat("Feasibility tables written to", cfg$out_dir, "\n") }
if (mode %in% c("A", "all")) { paper_A(t, cfg); cat("Paper A done\n") }
if (mode %in% c("B", "all")) { paper_B(t, cfg); cat("Paper B done\n") }
if (mode %in% c("C", "all")) { paper_C(t, cfg); cat("Paper C done\n") }
writeLines(c(paste("Run finished:", format(Sys.time())), paste("R:", R.version.string), paste("survival:", utils::packageVersion("survival")),
             paste("bootstrap replicates:", cfg$boot), paste("minimum cell:", cfg$min_cell)), file.path(cfg$out_dir, "run_info.txt"))
cat("Outputs are aggregate tables/figures only. Review them before sharing.\n")
