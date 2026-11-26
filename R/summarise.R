## =============================================================================
##  summarise.R -- operating characteristics
##
##  The guarantee of framework.tex is per basket: P(reject basket c) <= alpha for every basket
##  with p_c <= p0. The platform-wide FWER (any null basket rejected) is NOT controlled at alpha
##  and is reported separately, never pooled with the per-basket error.
##  Open-ended designs: n_mean, N_mean and duration are truncated at the simulation horizon n_cap;
##  the columns open / open_any report how often that horizon was reached.
## =============================================================================

## Per basket: d has the columns of simulate_trials() plus the grouping columns `by`.
summarise_baskets <- function(d, by = c("recruitment", "scenario", "engine")) {
  g <- interaction(d[, c(by, "basket")], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(d, g), function(x) data.frame(
    x[1, c(by, "basket", "row", "col", "p_true", "null")],
    trials = nrow(x), reject = mean(x$reject), reject_se = sqrt(mean(x$reject) * (1 - mean(x$reject)) / nrow(x)),
    stopped = mean(x$stopped), open = if (is.null(x$open)) NA else mean(x$open),   # open: still recruiting at the simulation horizon
    n_mean = mean(x$n), time_mean = mean(x$time),
    post_bias = mean(x$post_mean - x$p_true), cover90 = mean(x$cover), pex = mean(x$pex),
    prior_mean = mean(x$prior_mean), log10E = mean(x$logE) / log(10), rhat_max = suppressWarnings(max(x$rhat, na.rm = TRUE)))))
  out$rhat_max[!is.finite(out$rhat_max)] <- NA
  rownames(out) <- NULL
  out
}

## Per platform: per-basket type I error (largest and mean over null baskets), FWER, power and size.
summarise_platform <- function(d, by = c("recruitment", "scenario", "engine")) {
  g <- interaction(d[, by], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(d, g), function(x) {
    tr <- split(x, x$trial); nul <- x[x$null, ]; alt <- x[!x$null, ]
    typeI <- if (nrow(nul)) tapply(nul$reject, nul$basket, mean) else NA
    data.frame(x[1, by], trials = length(tr), baskets_null = length(unique(nul$basket)),
      typeI_basket_max = max(typeI), typeI_basket_mean = mean(typeI),
      FWER = if (nrow(nul)) mean(vapply(tr, function(t) any(t$reject & t$null), TRUE)) else NA,
      power_basket_mean = if (nrow(alt)) mean(alt$reject) else NA,
      power_any = if (nrow(alt)) mean(vapply(tr, function(t) any(t$reject & !t$null), TRUE)) else NA,
      all_correct = mean(vapply(tr, function(t) all(t$reject == !t$null), TRUE)),
      open_any = if (is.null(x$open)) NA else mean(vapply(tr, function(t) any(t$open), TRUE)),
      N_mean = mean(vapply(tr, function(t) sum(t$n), 0)), duration = mean(vapply(tr, function(t) max(t$time), 0)),
      cover90 = mean(x$cover), mean_abs_error = mean(abs(x$post_mean - x$p_true)))
  }))
  rownames(out) <- NULL
  out
}
