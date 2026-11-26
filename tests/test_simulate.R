## The calendar-time simulator (R/simulate.R) and the summaries (R/summarise.R).
p0 <- 0.1; grid <- theta_grid(p0); lay <- grid_layout(2, 2)
design <- list(layout = lay, n_max = 15, rate = 1, entry = c(0, 0, 0, 5), look_every = 5, alpha = 0.05)
p <- c(0.1, 0.3, 0.1, 0.5)
engines <- list(none = engine_none(grid, 0.2), oneway = engine_oneway(grid, 0.2, n_max = 15))

for (nm in names(engines)) {
  d <- simulate_trials(p, design, engines[[nm]], n_sim = 6, seed = 5)
  check(sprintf("%s: one row per trial and basket", nm), nrow(d) == 24 && !anyNA(d[, c("reject", "n", "time", "logE")]))
  check(sprintf("%s: every basket decided (rejected, stopped or at n_max)", nm), all(d$reject | d$stopped | d$n == 15))
  check(sprintf("%s: rejected exactly when E >= 1/alpha", nm), all(d$reject == (d$logE >= log(20))))
  check(sprintf("%s: stopped only when 1/alpha is out of reach", nm),
        all(d$logE[d$stopped] + (15 - d$n[d$stopped]) * log(1 / p0) < log(20)))
  check(sprintf("%s: reproducible with the same seed", nm), identical(d, simulate_trials(p, design, engines[[nm]], n_sim = 6, seed = 5)))
}

## Predictability: outcomes that agree up to patient 10 of basket 1 give identical e-values up to
## the analysis at time 10, and the SAME prior for basket 1's block 11-15 (it is fixed at time 10).
set.seed(8); X1 <- array(runif(6 * 4 * 15) < rep(p, each = 6), c(6, 4, 15)); X2 <- X1; X2[, 1, 11:15] <- !X1[, 1, 11:15]
for (nm in names(engines)) {
  t1 <- attr(simulate_trials(p, design, engines[[nm]], 6, X = X1, trace = TRUE), "trace")
  t2 <- attr(simulate_trials(p, design, engines[[nm]], 6, X = X2, trace = TRUE), "trace")
  same_before <- identical(t1[[1]]$logE, t2[[1]]$logE) && identical(t1[[2]]$logE, t2[[2]]$logE)
  same_prior <- identical(t1[[3]]$prior_mean[, 1], t2[[3]]$prior_mean[, 1])          # prior used for block 11-15
  changed <- !identical(t1[[3]]$y[, 1], t2[[3]]$y[, 1])
  check(sprintf("%s: e-values before the change are identical", nm), same_before)
  check(sprintf("%s: prior of a block does not depend on the block's data", nm), same_prior && changed)
}

## Summaries: per-basket type I error and FWER are separate quantities, FWER >= each basket's rate
d <- simulate_trials(c(0.1, 0.1, 0.1, 0.1), modifyList(design, list(alpha = 0.3)), engines$none, n_sim = 40, seed = 2)
d$recruitment <- "balanced"; d$scenario <- "null"; d$engine <- "none"
b <- summarise_baskets(d); pl <- summarise_platform(d)
check("summaries: one row per basket", nrow(b) == 4)
check("summaries: FWER >= largest per-basket type I error", pl$FWER >= pl$typeI_basket_max)
check("summaries: largest per-basket type I error = max over baskets", abs(pl$typeI_basket_max - max(b$reject)) < 1e-12)
