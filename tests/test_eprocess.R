## The validity layer (R/eprocess.R): support, normalisation, block factor, (super)martingale identity.
p0 <- 0.1; grid <- theta_grid(p0)
lp0 <- qlogis(p0)
check("grid: logit(p0) is a cell boundary, so cells are entirely below or above p0",
      all(grid$upper[!grid$above] <= lp0 + 1e-12) && all(grid$lower[grid$above] >= lp0 - 1e-12))
check("grid: all support points in (p0, 1]", all(grid$support > p0 & grid$support <= 1))
for (q in c(0.05, 0.2, 0.37)) { g <- theta_grid(q)
  check(sprintf("grid aligned for p0 = %.2f", q), abs(g$lower[which(g$above)[1]] - qlogis(q)) < 1e-12) }

set.seed(11); G <- length(grid$theta)
lp <- matrix(rnorm(5 * G, sd = 3), 5)
w <- restrict(lp, grid)
check("restrict: rows sum to 1", max(abs(rowSums(w) - 1)) < 1e-12)
direct <- exp(lp[2, grid$above]) / sum(exp(lp[2, grid$above]))
check("restrict: weights proportional to the posterior on cells above p0", max(abs(w[2, ] - direct)) < 1e-12)

## block factor against a direct sum
worst <- 0
for (i in 1:200) {
  wi <- restrict(matrix(rnorm(G, sd = 4), 1), grid); n <- sample(0:15, 1); y <- sample(0:n, 1)
  direct <- sum(wi * (grid$support / p0)^y * ((1 - grid$support) / (1 - p0))^(n - y))
  worst <- max(worst, abs(exp(log_factor(wi, y, n, grid)) / direct - 1))
}
check("block factor = direct sum (200 random priors and blocks)", worst < 1e-12, paste("max rel. diff", fmt(worst)))
check("empty block has factor 1", all(abs(log_factor(w, rep(0, 5), rep(0, 5), grid)) < 1e-15))

## E_p[factor] = 1 at p = p0 and <= 1 below p0, for any prior on (p0, 1] (Proposition 1)
mart <- function(w, n, p) vapply(seq_len(nrow(w)), function(i)
  sum(dbinom(0:n, n, p) * exp(log_factor(w[rep(i, n + 1), , drop = FALSE], 0:n, rep(n, n + 1), grid))), 0)
m1 <- c(mart(w, 5, p0), mart(w, 12, p0)); m2 <- c(mart(w, 5, 0.05), mart(w, 12, 0.01))
check("E_p0[factor] = 1", max(abs(m1 - 1)) < 1e-12, paste("max |E - 1|", fmt(max(abs(m1 - 1)))))
check("E_p[factor] <= 1 for p < p0", all(m2 <= 1), paste("max", fmt(max(m2))))

## the checks of the validity layer reject invalid priors
check_error("unnormalised prior is rejected", check_priors(w * 1.01, grid))
check_error("negative weight is rejected", { v <- w; v[1, 1] <- -1e-3; v[1, 2] <- v[1, 2] + 1e-3; check_priors(v, grid) })
check_error("NA weight is rejected", { v <- w; v[1, 1] <- NA; check_priors(v, grid) })
check_error("support at p0 is rejected", { g2 <- grid; g2$support[1] <- p0; check_priors(w, g2) })
check_error("wrong number of support points is rejected", check_priors(w[, -1], grid))

## scaled Beta benchmark: a valid prior with (almost) the same restricted moments
b <- scaled_beta(w, grid)
check("scaled Beta: valid prior", isTRUE(check_priors(b, grid)))
dm <- abs(b %*% grid$support - w %*% grid$support)
check("scaled Beta: restricted mean within 0.005 of the exact one", max(dm) < 0.005, paste("max diff", fmt(max(dm))))
