# Aggregate-only export helpers. Automated suppression is a first pass, not a disclosure certificate: review before sharing.
suppress <- function(df, cols, k) {
  flag <- rep(FALSE, nrow(df))
  for (cc in intersect(cols, names(df))) { x <- df[[cc]]; bad <- !is.na(x) & x > 0 & x < k; flag <- flag | bad }
  # blank every count on a flagged row so the small value cannot be recovered from its neighbours on that row
  for (cc in intersect(cols, names(df))) df[[cc]][flag] <- NA
  df$suppressed <- ifelse(flag, paste0("<", k), "")
  df
}
write_out <- function(df, name, cfg, count_cols = c("n", "events", "patients", "count")) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(df)) return(invisible())
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols & !names(df) %in% count_cols] <- lapply(df[num_cols & !names(df) %in% count_cols], signif, 4)
  utils::write.csv(suppress(df, count_cols, cfg$min_cell), file.path(cfg$out_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
}
# One set of combined files per paper: <prefix>_status / _cells / _estimates / _cox / _balance / _curve_grid (+ figures/).
write_paper <- function(results, prefix, cfg) {
  bind <- function(field) { x <- lapply(results, function(r) if (is.null(r[[field]]) || !is.data.frame(r[[field]])) NULL else cbind(analysis = r$label, outcome = r$outcome, r[[field]]))
    x <- x[!vapply(x, is.null, logical(1))]; if (!length(x)) return(NULL)
    cols <- unique(unlist(lapply(x, names))); do.call(rbind, lapply(x, function(z) { z[setdiff(cols, names(z))] <- NA; z[cols] })) }
  status <- do.call(rbind, lapply(results, function(r) data.frame(analysis = r$label, outcome = r$outcome, status = r$status, n = r$n)))
  write_out(status, paste0(prefix, "_status"), cfg, "n")
  cl <- bind("cells"); cl$nonevents <- cl$n - cl$events
  write_out(cl, paste0(prefix, "_cells"), cfg, c("n", "events", "patients", "nonevents"))
  for (f in c("estimates", "cox", "balance", "curves")) write_out(bind(f), paste0(prefix, "_", if (f == "curves") "curve_grid" else f), cfg, character(0))
  dir.create(file.path(cfg$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  for (res in results) if (identical(res$status, "COMPLETED")) {
    grDevices::png(file.path(cfg$out_dir, "figures", paste0(prefix, "_", res$label, ".png")), width = 1400, height = 700, res = 150)
    cv <- res$curves; sl <- unique(cv$stratum); graphics::par(mfrow = c(1, length(sl)), mar = c(4, 4, 3, 1))
    for (s in sl) {
      graphics::plot(NA, xlim = c(0, 36), ylim = c(0, 1), xlab = "Months", ylab = paste("Adjusted", toupper(res$outcome)), main = s)
      for (g in 0:1) { z <- cv[cv$stratum == s & cv$group == g, ]; graphics::lines(z$month, z$survival, type = "b", col = c("#1f77b4", "#d62728")[g + 1], lwd = 2, pch = 16) }
      graphics::legend("topright", c("group 0 (wild-type / reference)", "group 1 (altered / class)"), col = c("#1f77b4", "#d62728"), lwd = 2, bty = "n")
    }
    grDevices::dev.off()
  }
}
