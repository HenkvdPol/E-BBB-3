## =============================================================================
##  reference_twoway.R -- reference calculations for the two-way engine (validation only)
##
##  twoway_exact():  deterministic quadrature for ONE state, layouts with at most 3 columns.
##                   Column levels c_l = mu + b_l on a grid (mu integrated in closed form), and
##                   given c each row effect a_j in one dimension (rows are independent given c).
##                   Same tables and theta step as the engine (borrow_twoway.R), so differences
##                   come from the integration over (mu, a, b, sigmas) only.
##  twoway_jags():   the model written out independently in JAGS (continuous theta and sigmas).
##  oneway_equivalent_hyper(): an L x 1 two-way layout is a one-way EXNEX model with
##                   mu' = mu + b_1 and tau^2 = sigma_a^2 + sigma_e^2; the (mu', tau) grid and prior.
## =============================================================================

## log p(c | sigma_b) for column levels c (rows of a matrix): mu ~ TN(m0, s0^2, [lo, hi]) integrated out.
logp_columns <- function(cmat, sb, prior) {
  L <- ncol(cmat); m0 <- prior$mu_mean; s0 <- prior$mu_sd; lo <- prior$mu_range[1]; hi <- prior$mu_range[2]
  v <- 1 / (1 / s0^2 + L / sb^2); mh <- v * (m0 / s0^2 + rowSums(cmat) / sb^2)
  -rowSums(cmat^2) / (2 * sb^2) + mh^2 / (2 * v) - L * log(sb) + 0.5 * log(v) +
    log(pmax(pnorm((hi - mh) / sqrt(v)) - pnorm((lo - mh) / sqrt(v)), 1e-300))
}
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

twoway_exact <- function(y, n, layout, tab, prior, grid, c_half = 4, a_max = 5) {
  K <- length(y); J <- max(layout$row); L <- max(layout$col)
  if (L > 3) stop("twoway_exact: at most 3 columns (transpose the layout if needed)")
  rb <- seq(1, tab$Gm, by = round(tab$rb_step / tab$m_step))                 # RB lattice -> m lattice
  centre <- qlogis((sum(y) + 0.5) / (sum(n) + 1))
  ic <- which(abs(tab$rb_m - centre) <= c_half + 1e-9); nc <- length(ic)     # column levels on the lattice
  A <- round(a_max / tab$rb_step); off <- -A:A; av <- off * tab$rb_step; Na <- length(off)
  cg <- as.matrix(expand.grid(rep(list(seq_len(nc)), L))); ci <- matrix(ic[cg], ncol = L); cv <- matrix(tab$rb_m[ci], ncol = L)
  lpc <- apply(sapply(seq_along(prior$sig_b$value), function(i) logp_columns(cv, prior$sig_b$value[i], prior) + prior$sig_b$logprior[i]),
               1, lse)
  rho <- tab$row_of(y, n)
  ld <- function(k, ie, idx) tab$logD[cbind(rho[k], (ie - 1) * tab$Gm + rb[pmin(pmax(idx, 1), tab$Grb)])]
  rows <- lapply(seq_len(J), function(j) which(layout$row == j))
  comb <- expand.grid(ia = seq_along(prior$sig_a$value), ie = seq_along(prior$sig_e$value))
  la_of <- function(ia) { la <- dnorm(av, 0, prior$sig_a$value[ia], log = TRUE); la - lse(la) }
  vt <- function(j, t, ie, la) { v <- rep(la[t], nrow(cg)); for (k in rows[[j]]) v <- v + ld(k, ie, ci[, layout$col[k]] + off[t]); v }
  ## pass 1: log R_j(c) = log sum_a N(a) prod_{k in row j} D_k(c_l + a), and the joint weight of (c, sigmas)
  logR <- array(0, c(nrow(cg), J, nrow(comb))); logW <- matrix(0, nrow(cg), nrow(comb))
  for (h in seq_len(nrow(comb))) {
    la <- la_of(comb$ia[h])
    for (j in seq_len(J)) if (length(rows[[j]])) {
      V <- vapply(seq_len(Na), function(t) vt(j, t, comb$ie[h], la), numeric(nrow(cg)))
      mx <- row_max(V); logR[, j, h] <- mx + log(rowSums(exp(V - mx)))
    }
    logW[, h] <- lpc + rowSums(logR[, , h, drop = FALSE]) + prior$sig_a$logprior[comb$ia[h]] + prior$sig_e$logprior[comb$ie[h]]
  }
  logW <- logW - lse(logW)
  ## pass 2: distribution of (m_c, sigma_e) on the Rao-Blackwell lattice
  q <- matrix(0, K, tab$Gs * tab$Grb); edge_c <- sum(exp(logW)[rowSums(cg == 1 | cg == nc) > 0, ]); edge_a <- 0
  for (h in seq_len(nrow(comb))) {
    la <- la_of(comb$ia[h]); wc <- exp(logW[, h]); ie <- comb$ie[h]
    for (j in seq_len(J)) for (t in seq_len(Na)) {
      if (!length(rows[[j]])) next
      P <- wc * exp(vt(j, t, ie, la) - logR[, j, h])
      if (t %in% c(1, Na)) edge_a <- max(edge_a, sum(P))
      for (k in rows[[j]]) {
        s <- rowsum(P, cg[, layout$col[k]])[, 1]; idx <- pmin(pmax(ic + off[t], 1), tab$Grb)
        cols <- (ie - 1) * tab$Grb + idx
        q[k, ] <- q[k, ] + tabulate_w(cols, s, ncol(q))
      }
    }
  }
  th <- theta_from_q(Matrix::Matrix(q / rowSums(q), sparse = TRUE), y, n, tab, grid, prior$w_ex)
  list(logpost = th$logpost, pex = th$pex, edge_c = edge_c, edge_a = edge_a)
}
tabulate_w <- function(idx, w, n) { out <- numeric(n); s <- rowsum(w, idx); out[as.integer(rownames(s))] <- s[, 1]; out }

## JAGS: the same model with continuous theta (truncated to the grid range) and continuous sigmas
## (half-normal, truncated at sigma_max), or fixed sigmas (prior$fixed).
twoway_jags <- function(y, n, layout, prior, grid, n_iter = 20000, n_burn = 5000, n_chains = 4, seed = 1, sigma_max = 2.5) {
  if (!requireNamespace("rjags", quietly = TRUE)) stop("rjags is not installed")
  fx <- prior$fixed
  sig <- function(nm) if (is.null(fx)) sprintf("s%s ~ dnorm(0, 1 / %g^2) T(0, %g)", nm, prior$scale[[nm]], sigma_max) else ""
  model <- paste(
    "model {",
    "  mu ~ dnorm(m0, 1 / s0^2) T(lo, hi)",
    sig("a"), sig("b"), sig("e"),
    "  for (j in 1:J) { a[j] ~ dnorm(0, 1 / sa^2) }",
    "  for (l in 1:L) { b[l] ~ dnorm(0, 1 / sb^2) }",
    "  for (k in 1:K) {",
    "    d[k] ~ dbern(w)",
    "    mk[k] <- d[k] * (mu + a[row[k]] + b[col[k]]) + (1 - d[k]) * mN",
    "    pk[k] <- d[k] / se^2 + (1 - d[k]) / vN",
    "    theta[k] ~ dnorm(mk[k], pk[k]) T(tlo, thi)",
    "    y[k] ~ dbin(ilogit(theta[k]), n[k])",
    "  }",
    "}", sep = "\n")
  pb <- prior$pbar
  data <- list(y = y, n = n, row = layout$row, col = layout$col, K = length(y), J = max(layout$row), L = max(layout$col),
               m0 = prior$mu_mean, s0 = prior$mu_sd, lo = prior$mu_range[1], hi = prior$mu_range[2], w = prior$w_ex,
               mN = qlogis(pb), vN = 1 / pb + 1 / (1 - pb), tlo = grid$lower[1], thi = grid$upper[length(grid$upper)])
  if (!is.null(fx)) data <- c(data, list(sa = fx[["a"]], sb = fx[["b"]], se = fx[["e"]]))
  inits <- lapply(seq_len(n_chains), function(i) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + i))
  jm <- rjags::jags.model(textConnection(model), data, inits, n.chains = n_chains, quiet = TRUE)
  update(jm, n_burn, progress.bar = "none")
  sm <- rjags::coda.samples(jm, c("theta", "d"), n_iter, progress.bar = "none")
  x <- do.call(rbind, sm); th <- x[, sprintf("theta[%d]", seq_along(y)), drop = FALSE]
  list(theta = th, d = x[, sprintf("d[%d]", seq_along(y)), drop = FALSE],
       rhat = coda::gelman.diag(sm[, sprintf("theta[%d]", seq_along(y))], multivariate = FALSE)$psrf[, 1],
       ess = coda::effectiveSize(sm[, sprintf("theta[%d]", seq_along(y))]))
}

## Summaries compared across methods: P(p > p0), restricted mean and sd, P(EX), and the block
## factor of a block of 5 with 0, 2 or 5 responses under the restricted prior.
summ_grid <- function(logpost, pex, grid) {
  pr <- exp(logpost - row_max(logpost)); pr <- pr / rowSums(pr); w <- restrict(logpost, grid)
  rm <- as.vector(w %*% grid$support); rs <- sqrt(pmax(as.vector(w %*% grid$support^2) - rm^2, 0))
  f <- sapply(c(0, 2, 5), function(yb) exp(log_factor(w, rep(yb, nrow(w)), rep(5, nrow(w)), grid)))
  data.frame(p_above = rowSums(pr[, grid$above, drop = FALSE]), rmean = rm, rsd = rs, pex = pex,
             f0 = f[, 1], f2 = f[, 2], f5 = f[, 3])
}
summ_draws <- function(theta, d, grid) {
  p <- plogis(theta); ab <- theta > qlogis(grid$p0)
  rm <- colSums(p * ab) / colSums(ab); rs <- sqrt(pmax(colSums(p^2 * ab) / colSums(ab) - rm^2, 0))
  f <- sapply(c(0, 2, 5), function(yb) colSums(ab * (p / grid$p0)^yb * ((1 - p) / (1 - grid$p0))^(5 - yb)) / colSums(ab))
  data.frame(p_above = colMeans(ab), rmean = rm, rsd = rs, pex = colMeans(d), f0 = f[, 1], f2 = f[, 2], f5 = f[, 3])
}

## The one-way hyperprior equivalent to an L x 1 two-way layout (all baskets in one column).
oneway_equivalent_hyper <- function(prior, mu_step = 0.1, tau_digits = 2) {
  mu <- seq(prior$mu_range[1] - 1, prior$mu_range[2] + 1, by = mu_step)
  lmu <- apply(sapply(seq_along(prior$sig_b$value), function(i)
    logp_columns(matrix(mu), prior$sig_b$value[i], prior) + prior$sig_b$logprior[i]), 1, lse)
  tau <- round(sqrt(outer(prior$sig_a$value^2, prior$sig_e$value^2, "+")), tau_digits)
  pt <- tapply(exp(outer(prior$sig_a$logprior, prior$sig_e$logprior, "+")), tau, sum)
  h <- expand.grid(mu = mu, tau = as.numeric(names(pt)))
  h$logprior <- rep(lmu, length(pt)) + rep(log(as.vector(pt)), each = length(mu))
  h
}
