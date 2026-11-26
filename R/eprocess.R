## =============================================================================
##  eprocess.R -- the validity layer of framework.tex. Nothing here knows about borrowing.
##
##  Basket c is analysed at times tau_0 < tau_1 < ...; block k holds its patients between
##  tau_{k-1} and tau_k. Before block k a prior pi_k on (p0, 1] is fixed from the platform
##  history and frozen for the whole block. The block factor is
##      e_k = int prod_{i in block} (p / p0)^Y_i ((1 - p) / (1 - p0))^(1 - Y_i) dpi_k(p),
##  E_k = e_1 * ... * e_k, and H0: p_c <= p0 is rejected at the first k with E_k >= 1/alpha.
##  Under the conditional null every factor has conditional mean <= 1 (framework.tex,
##  Proposition 1), so P(sup_k E_k >= 1/alpha) <= alpha for ANY rule that chooses the prior
##  from the past only and puts it on (p0, 1]. This file checks exactly those two things
##  (support and normalisation); predictability is guaranteed by the simulator (simulate.R).
##
##  Here a prior is a discrete measure on the grid points p_g > p0: one row of weights per
##  basket, rows summing to 1. A borrowing engine (borrow_*.R) returns the posterior of
##  theta = logit(p) on the same grid; restrict() turns it into the prior of eq. (8).
## =============================================================================

## Grid of theta cells of width h on `range`, with logit(p0) on a cell boundary, so every cell
## lies entirely below or entirely above p0. Cells above p0 form the support of the priors.
theta_grid <- function(p0, h = 0.05, range = c(-10, 6)) {
  g <- seq(ceiling((range[1] - qlogis(p0)) / h + 0.5), floor((range[2] - qlogis(p0)) / h + 0.5))
  theta <- qlogis(p0) + (g - 0.5) * h
  above <- g >= 1
  p <- plogis(theta)
  list(theta = theta, lower = theta - h / 2, upper = theta + h / 2, p = p, h = h, p0 = p0, range = range,
       above = above, support = p[above],
       l1 = log(p[above] / p0), l0 = log((1 - p[above]) / (1 - p0)))   # log LR per success / failure
}

## Probability of each grid cell under N(mean, sd^2), renormalised over the grid.
## Returns a (cells x length(mean)) matrix. Upper tails use lower.tail = FALSE for accuracy.
norm_cells <- function(grid, mean, sd) {
  G <- length(grid$theta); M <- max(length(mean), length(sd))
  mean <- rep(rep_len(mean, M), each = G); sd <- rep(rep_len(sd, M), each = G)
  a <- (rep(grid$lower, M) - mean) / sd; b <- (rep(grid$upper, M) - mean) / sd
  P <- matrix(ifelse(a > 0, pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE), pnorm(b) - pnorm(a)), G, M)
  sweep(P, 2, colSums(P), "/")
}

row_max <- function(x) x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]

## The prior of framework.tex eq. (8): the posterior restricted to (p0, 1] and renormalised.
## logpost: (baskets x cells) log posterior on the grid, up to a constant per row.
restrict <- function(logpost, grid) {
  lw <- logpost[, grid$above, drop = FALSE]
  w <- exp(lw - row_max(lw)); w <- w / rowSums(w)
  check_priors(w, grid)
  w
}

## The only requirements of the validity layer: support in (p0, 1], finite weights >= 0 summing to 1.
check_priors <- function(w, grid, tol = 1e-9) {
  if (!all(grid$support > grid$p0 & grid$support <= 1)) stop("prior support must lie in (p0, 1]")
  if (!is.matrix(w) || ncol(w) != length(grid$support)) stop("prior weights must be a matrix with one column per support point")
  if (any(!is.finite(w)) || any(w < 0)) stop("prior weights must be finite and >= 0")
  if (any(abs(rowSums(w) - 1) > tol)) stop("prior weights must sum to 1")
  invisible(TRUE)
}

## log block factor for each row of w: y responses among n patients (n = 0 gives 0).
log_factor <- function(w, y, n, grid) {
  lr <- outer(y, grid$l1) + outer(n - y, grid$l0)
  m <- row_max(lr)
  m + log(rowSums(w * exp(lr - m)))
}

## Benchmark representation (decision D4): the scaled Beta with the same restricted mean and
## second moment (p = p0 + (1 - p0) u, u ~ Beta(a, b)), put on the same cells. It is a
## different admissible prior on (p0, 1], so validity is unaffected (framework.tex, remark
## 'Computational approximation').
scaled_beta <- function(w, grid) {
  p0 <- grid$p0; m <- as.vector(w %*% grid$support); s2 <- as.vector(w %*% grid$support^2)
  u <- (m - p0) / (1 - p0); uv <- (s2 - m^2) / (1 - p0)^2
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6); uv <- pmin(pmax(uv, 1e-12), u * (1 - u) * (1 - 1e-6))
  sz <- u * (1 - u) / uv - 1; a <- u * sz; b <- (1 - u) * sz
  ub <- c((plogis(grid$lower[grid$above]) - p0) / (1 - p0), 1)              # cell boundaries on the u scale
  out <- t(vapply(seq_along(a), function(i) diff(pbeta(ub, a[i], b[i])), numeric(length(grid$support))))
  out <- pmax(out, 0); out <- out / rowSums(out)
  check_priors(out, grid)
  out
}
