## Validity does not depend on the borrowing engine: exact type I error of one null basket.
## Basket 1 has p = p0; the other baskets follow a fixed (adversarially chosen) path of responses. Blocks of 5 up to n = 30 (6 blocks); all 6^6 block outcomes
## of basket 1 are enumerated, so the type I error is exact (no simulation).
p0 <- 0.1; alpha <- 0.05; grid <- theta_grid(p0); block <- 5; B <- 6

exact_typeI <- function(engine, K, rate_other = 0.6, prior_form = "grid") {
  lf <- vector("list", B)
  for (k in seq_len(B)) {                                  # priors for block k, for every own count before it
    y1 <- 0:(block * (k - 1)); nk <- block * (k - 1)
    Y <- cbind(y1, matrix(round(rate_other * nk), length(y1), K - 1)); N <- matrix(nk, length(y1), K)
    need <- matrix(FALSE, length(y1), K); need[, 1] <- TRUE
    w <- restrict(matrix(engine$posterior(Y, N, need, seed = k)$logpost[, 1, ], length(y1)), grid)
    if (prior_form == "beta") w <- scaled_beta(w, grid)
    lf[[k]] <- sapply(0:block, function(yb) log_factor(w, rep(yb, nrow(w)), rep(block, nrow(w)), grid))
  }
  paths <- as.matrix(expand.grid(rep(list(0:block), B))); prev <- cbind(0, t(apply(paths, 1, cumsum))[, -B])
  logE <- 0; rej <- FALSE
  for (k in seq_len(B)) { logE <- logE + lf[[k]][cbind(prev[, k] + 1, paths[, k] + 1)]; rej <- rej | logE >= log(1 / alpha) }
  sum(exp(rowSums(dbinom(paths, block, p0, log = TRUE))) * rej)
}

## Platforms whose other baskets respond at rate 0.15 (prior pulled to just above p0, the hardest
## case for the null basket) or 0.6 (prior pulled far above p0).
tw <- engine_twoway(grid, grid_layout(2, 2), 0.2, n_max = 30, mcmc = list(burn = 300, keep = 500, chains = 2, chunk = 150))
ow <- engine_oneway(grid, 0.2, n_max = 30)
res <- c(none = exact_typeI(engine_none(grid, 0.2), 3),
         oneway_r15 = exact_typeI(ow, 3, 0.15), oneway_r60 = exact_typeI(ow, 3, 0.6),
         oneway_beta_r15 = exact_typeI(ow, 3, 0.15, prior_form = "beta"),
         twoway_r15 = exact_typeI(tw, 4, 0.15), twoway_r60 = exact_typeI(tw, 4, 0.6))
for (nm in names(res)) check(sprintf("exact type I error <= alpha: engine/platform '%s'", nm), res[[nm]] <= alpha,
                             sprintf("%.4f", res[[nm]]))
