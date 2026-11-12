# Weighted Kaplan-Meier with delayed entry, RMST, inverse-probability weights, cluster bootstrap, Cox helpers.
suppressPackageStartupMessages(library(survival))

km <- function(time, event, w = NULL, entry = NULL) {
  if (is.null(w)) w <- rep(1, length(time))
  ev <- event == 1
  ut <- sort(unique(time[ev]))
  if (!length(ut)) return(list(time = numeric(0), surv = numeric(0), nrisk = numeric(0)))
  below <- function(x, ww) { o <- order(x); c(0, cumsum(ww[o]))[findInterval(ut, x[o], left.open = TRUE) + 1] }  # sum ww[x < t]
  risk <- sum(w) - below(time, w); nrisk <- length(time) - below(time, rep(1, length(time)))
  if (!is.null(entry)) { risk <- risk - (sum(w) - below(entry, w)); nrisk <- nrisk - (length(entry) - below(entry, rep(1, length(entry)))) }
  d <- as.numeric(tapply(w[ev], factor(time[ev], levels = ut), sum))
  f <- ifelse(risk > 0, 1 - d / risk, 1)
  list(time = ut, surv = cumprod(pmin(pmax(f, 0), 1)), nrisk = nrisk)
}
rmst_km <- function(k, tau) { i <- k$time < tau; sum(diff(c(0, k$time[i], tau)) * c(1, k$surv[i])) }
surv_km <- function(k, t) vapply(t, function(x) { j <- which(k$time <= x); if (length(j)) k$surv[max(j)] else 1 }, numeric(1))
rmst <- function(time, event, tau, w = NULL, entry = NULL) rmst_km(km(time, event, w, entry), tau)

lump <- function(x, min_n = 10) {       # merge sparse categories so that propensity models stay estimable
  x <- as.character(x); tb <- table(x); small <- names(tb)[tb < min_n]
  if (length(small) && length(tb) - length(small) >= 1) x[x %in% small] <- "Other_sparse"
  factor(x)
}
model_frame <- function(d, covars) {
  z <- d[covars]
  for (nm in covars) {
    if (is.character(z[[nm]]) || is.factor(z[[nm]])) z[[nm]] <- lump(z[[nm]])
    else if (anyNA(z[[nm]])) z[[nm]][is.na(z[[nm]])] <- stats::median(z[[nm]], na.rm = TRUE)
  }
  z[vapply(z, function(x) length(unique(x)) > 1, logical(1))]
}
# Stabilised weights for binary `group`, fitted separately within each level of `by`; truncated at the 99th percentile.
ipw <- function(d, group, covars, by = NULL, bound = 0.02) {
  w <- rep(1, nrow(d)); strata <- if (is.null(by)) rep(1, nrow(d)) else d[[by]]
  for (s in unique(strata)) {
    i <- which(strata == s); g <- d[[group]][i]
    if (length(unique(g)) < 2) next
    z <- model_frame(d[i, , drop = FALSE], covars); z$.g <- g
    fit <- suppressWarnings(stats::glm(.g ~ ., data = z, family = stats::binomial()))
    ps <- pmin(pmax(as.numeric(stats::fitted(fit)), bound), 1 - bound)
    ww <- ifelse(g == 1, mean(g) / ps, (1 - mean(g)) / (1 - ps))
    w[i] <- pmin(ww, stats::quantile(ww, 0.99, names = FALSE))
  }
  w
}
smd_table <- function(d, group, covars, w, by = NULL) {
  strata <- if (is.null(by)) rep("all", nrow(d)) else as.character(d[[by]]); out <- list()
  for (s in unique(strata)) {
    i <- strata == s; z <- model_frame(d[i, , drop = FALSE], covars); if (!ncol(z)) next
    x <- stats::model.matrix(~ ., z)[, -1, drop = FALSE]; g <- d[[group]][i] == 1; ww <- w[i]
    for (j in seq_len(ncol(x))) {
      den <- sqrt((stats::var(x[g, j]) + stats::var(x[!g, j])) / 2)
      f <- function(wt) if (is.na(den) || den == 0) 0 else (stats::weighted.mean(x[g, j], wt[g]) - stats::weighted.mean(x[!g, j], wt[!g])) / den
      out[[length(out) + 1]] <- data.frame(stratum = s, covariate = colnames(x)[j], smd_unweighted = f(rep(1, sum(i))), smd_weighted = f(ww))
    }
  }
  do.call(rbind, out)
}

# RMST by group (0/1) within stratum levels; difference (group 1 - group 0) per stratum; difference of differences.
cell_rmst <- function(d, group, stratum, taus, w) {
  out <- c(); lev <- sort(unique(d[[stratum]]))
  ent <- if (any(d$entry > 0)) d$entry else NULL
  for (tau in taus) {
    diffs <- c()
    for (s in lev) {
      r <- c()
      for (g in 0:1) {
        i <- d[[stratum]] == s & d[[group]] == g
        r[as.character(g)] <- if (sum(i)) rmst(d$time[i], d$event[i], tau, w[i], if (is.null(ent)) NULL else ent[i]) else NA
        out[sprintf("rmst%g|%s|g%d", tau, s, g)] <- r[as.character(g)]
      }
      diffs[as.character(s)] <- r["1"] - r["0"]; out[sprintf("diff%g|%s", tau, s)] <- diffs[as.character(s)]
    }
    if (length(lev) == 2) out[sprintf("did%g", tau)] <- diffs[2] - diffs[1]
  }
  out
}
boot_cluster <- function(d, fn, B, seed) {
  est <- fn(d); set.seed(seed)
  rows <- split(seq_len(nrow(d)), d$record_id); ids <- names(rows)
  reps <- matrix(NA_real_, B, length(est))
  for (b in seq_len(B)) {
    r <- tryCatch(fn(d[unlist(rows[sample(ids, length(ids), replace = TRUE)], use.names = FALSE), , drop = FALSE]), error = function(e) NULL)
    if (!is.null(r) && length(r) == length(est)) reps[b, ] <- r
  }
  ok <- stats::complete.cases(reps); R <- reps[ok, , drop = FALSE]
  se <- apply(R, 2, stats::sd); ci <- apply(R, 2, stats::quantile, c(.025, .975), names = FALSE)
  data.frame(parameter = names(est), estimate = as.numeric(est), se = se, ci_low = ci[1, ], ci_high = ci[2, ],
             p_wald = ifelse(grepl("^(diff|did)", names(est)) & se > 0, 2 * stats::pnorm(-abs(est / se)), NA), boot_ok = sum(ok), row.names = NULL)
}
cells <- function(d, group, stratum) {
  a <- stats::aggregate(cbind(n = 1, events = d$event, person_months = d$time - d$entry) ~ d[[stratum]] + d[[group]], FUN = sum)
  names(a)[1:2] <- c("stratum", "group"); a$patients <- as.numeric(tapply(d$record_id, list(d[[stratum]], d[[group]]), function(x) length(unique(x)))[cbind(as.character(a$stratum), as.character(a$group))])
  a
}
gate_cells <- function(cl, cfg, expected) nrow(cl) == expected && all(cl$n >= cfg$min_n_cell & cl$events >= cfg$min_ev_cell)

# Outcome set-up: `os` or `pfs`; drops rows not under observation after (delayed) entry.
set_outcome <- function(d, outcome) {
  d$time <- if (outcome == "os") d$os_m else d$pfs_m
  d$event <- if (outcome == "os") d$death else d$prog
  d$entry <- d$entry_m
  d[!is.na(d$time) & d$event %in% 0:1 & d$time > d$entry & d$time > 0, ]
}
# Full analysis of one contrast: counts, gate, IPW RMST contrasts with cluster bootstrap, Cox with interaction.
contrast_analysis <- function(d, group, stratum, outcome, taus, covars, cfg, label, adjust = TRUE) {
  d <- set_outcome(d, outcome); cl <- cells(d, group, stratum)
  nlev <- length(unique(d[[stratum]])); res <- list(label = label, outcome = outcome, cells = cl, n = nrow(d))
  if (!gate_cells(cl, cfg, 2 * nlev)) { res$status <- "STOPPED_SPARSE_CELLS"; return(res) }
  by <- if (nlev > 1) stratum else NULL
  fn <- function(z) cell_rmst(z, group, stratum, taus, if (adjust) ipw(z, group, covars, by) else rep(1, nrow(z)))
  res$estimates <- boot_cluster(d, fn, cfg$boot, cfg$seed)
  w <- if (adjust) ipw(d, group, covars, by) else rep(1, nrow(d))
  res$balance <- smd_table(d, group, covars, w, by)
  res$cox <- tryCatch(cox_contrast(d, group, stratum, covars), error = function(e) data.frame(term = "COX_FAILED", note = conditionMessage(e)))
  res$curves <- curve_grid(d, group, stratum, w, cfg)
  res$status <- "COMPLETED"; res
}
cox_contrast <- function(d, group, stratum, covars) {
  z <- model_frame(d, covars); z$.g <- d[[group]]; z$.s <- factor(d[[stratum]]); z$.id <- d$record_id
  z$entry <- d$entry; z$time <- d$time; z$event <- d$event
  two <- nlevels(z$.s) == 2
  rhs <- paste(c(".g", if (two) c(".g:.s", "strata(.s)"), setdiff(names(z), c(".g", ".s", ".id", "entry", "time", "event"))), collapse = " + ")
  fit <- coxph(stats::as.formula(paste("Surv(entry, time, event) ~", rhs, "+ cluster(.id)")), data = z)
  b <- stats::coef(fit); V <- stats::vcov(fit); k <- grep("^\\.g", names(b))
  row <- function(nm, est, se) data.frame(term = nm, hr = exp(est), ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se), p = 2 * stats::pnorm(-abs(est / se)))
  lev <- levels(z$.s); out <- row(paste0("HR group1 vs group0 | ", lev[1]), b[k[1]], sqrt(V[k[1], k[1]]))
  if (two && length(k) == 2 && !is.na(b[k[2]])) {
    out <- rbind(out, row(paste0("HR group1 vs group0 | ", lev[2]), sum(b[k]), sqrt(sum(V[k, k]))),
                 row("ratio of HRs (multiplicative interaction)", b[k[2]], sqrt(V[k[2], k[2]])))
  }
  rownames(out) <- NULL; out
}
# Survival on a coarse monthly grid (no individual event times), only while >= min_cell patients remain at risk.
curve_grid <- function(d, group, stratum, w, cfg, grid = seq(0, 36, 3)) {
  out <- list()
  for (s in sort(unique(d[[stratum]]))) for (g in 0:1) {
    i <- d[[stratum]] == s & d[[group]] == g; if (!sum(i)) next
    k <- km(d$time[i], d$event[i], w[i], if (any(d$entry > 0)) d$entry[i] else NULL)
    atrisk <- vapply(grid, function(x) sum(d$time[i] >= x & d$entry[i] <= x), numeric(1))
    out[[length(out) + 1]] <- data.frame(stratum = s, group = g, month = grid, survival = ifelse(atrisk >= cfg$min_cell, round(surv_km(k, grid), 3), NA))
  }
  do.call(rbind, out)
}
