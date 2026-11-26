## The simulation (additional_functions.R, Part 3) and the summaries (Part 4).
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

## Regression: expanding the initial K x G prior to S trials must copy rows, not recycle
## the matrix's column-major storage. Distinct basket priors detect both row and support errors.
initial_logpost <- t(log(norm_cells(grid, qlogis(c(0.18, 0.32, 0.55, 0.75)), c(0.3, 0.5, 0.7, 0.9))))
distinct_engine <- list(name = "distinct deterministic priors", grid = grid,
  posterior = function(Y, N, need = NULL, seed = NULL) {
    lp <- array(NA_real_, c(nrow(Y), ncol(Y), length(grid$theta)))
    for (k in seq_len(ncol(Y))) {
      lp[, k, ] <- loglik_grid(Y[, k], N[, k], grid) + rep(initial_logpost[k, ], each = nrow(Y))
    }
    list(logpost = lp, pex = matrix(NA_real_, nrow(Y), ncol(Y)))
  })
batch_design <- modifyList(design, list(n_max = 10, entry = 0))
batch_X <- array(FALSE, c(3, 4, 10))
for (s in seq_len(3)) for (k in seq_len(4)) {
  batch_X[s, k, ] <- (seq_len(10) + 2 * s + 3 * k) %% (k + 2) == 0
}
batch <- simulate_trials(p, batch_design, distinct_engine, 3, X = batch_X, trace = TRUE)
first <- attr(batch, "trace")[[1]]
initial_w <- restrict(initial_logpost, grid)
expected_logE <- matrix(NA_real_, 3, 4)
for (s in seq_len(3)) for (k in seq_len(4)) {
  successes <- sum(batch_X[s, k, 1:5])
  expected_logE[s, k] <- log(sum(initial_w[k, ] * (grid$support / p0)^successes *
    ((1 - grid$support) / (1 - p0))^(5 - successes)))
}
check("initial priors: batched first-block factors equal direct normalized mixtures",
      max(abs(first$logE - expected_logE)) < 1e-12)
check("initial priors: every trial uses its basket's unchanged prior mean",
      max(abs(first$prior_mean - matrix(as.vector(initial_w %*% grid$support), 3, 4, byrow = TRUE))) < 1e-12)

## A deterministic engine with fixed outcomes must give the same decisions, factors and
## posterior summaries when trials are simulated together or one at a time.
same_trials <- vapply(seq_len(3), function(s) {
  one <- simulate_trials(p, batch_design, distinct_engine, 1, X = batch_X[s, , , drop = FALSE])
  together <- batch[batch$trial == s, names(one), drop = FALSE]
  together$trial <- 1L
  attr(together, "trace") <- NULL
  rownames(together) <- rownames(one) <- NULL
  isTRUE(all.equal(together, one, tolerance = 1e-12))
}, TRUE)
check("simulation batching: fixed outcomes match separate one-trial runs", all(same_trials))
