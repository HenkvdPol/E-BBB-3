## =============================================================================
##  borrow_twoway.R -- two-way EXNEX (framework.tex eq. 9) for a J x L basket grid
##
##  Basket c = (j, l): tumour type j (row), biomarker group l (column).
##    delta_c ~ Bernoulli(w), independent across baskets
##    EX  (delta_c = 1): theta_c ~ N(mu + a_j + b_l, sigma_e^2)
##    NEX (delta_c = 0): theta_c ~ N(logit pbar, 1/pbar + 1/(1 - pbar))
##    mu ~ N(logit p0, 10^2) on [-8, 4];  a_j ~ N(0, sigma_a^2);  b_l ~ N(0, sigma_b^2)
##    sigma_a, sigma_b, sigma_e ~ HalfNormal(1/sqrt(3)), on the grid 0.05, 0.15, ..., 2.45
##  Scales match the one-way benchmark's total expected EX variance before truncation/gridding:
##  E(sigma_a^2 + sigma_b^2 + sigma_e^2) = 1 = E(tau^2) for tau ~ HN(1). This does not match
##  the full variance distribution; finite grids also perturb the expectation slightly (decision D2).
##  The mu, w and NEX specifications follow the benchmark; its mu integration uses a grid.
##  Identifiability (D5): proper random effects, no sum-to-zero constraint; theta_c is identified.
##  theta lives on the cells of the e-process grid (eprocess.R); given the EX mean m and sigma_e
##  its prior is the renormalised cell probabilities of N(m, sigma_e^2).
##
##  Computation. Given m_c = mu + a_j + b_l and sigma_e, theta_c and delta_c depend only on the
##  basket's own data, so both are summed out on the theta grid at each table point:
##      D_c(m, sigma_e) = w * sum_theta L_c(theta) N(theta | m, sigma_e) + (1 - w) * sum_theta L_c(theta) V(theta).
##  What remains, (mu, a, b, sigmas), is integrated by
##    "mcmc"  (production): Metropolis-within-Gibbs, vectorised over states and chains, with exact
##            Gibbs moves along the three directions that leave every m_c unchanged;
##    quadrature for a state without any patients (prior predictive, no sampling);
##    reference_twoway.R: deterministic quadrature and JAGS, for validation only.
##  The posterior of theta_c is then Rao-Blackwellised over the draws of (m_c, sigma_e):
##      p(theta_c | y) = E[ L_c(theta) {w N(theta | m_c, sigma_e) + (1 - w) V(theta)} / D_c(m_c, sigma_e) ].
##  The finite sums at table points are direct evaluations. Interpolation/clamping of m, finite MCMC,
##  and deposition on the Rao-Blackwell lattice approximate the posterior weights. The e-process
##  sums the supplied restricted weights directly; validity requires their support and normalisation,
##  plus freezing before the next block, rather than exact posterior integration.
## =============================================================================

## Layout of the trial: one row per basket in the trial (missing cells are simply left out).
grid_layout <- function(J, L, missing = NULL) {
  lay <- expand.grid(row = seq_len(J), col = seq_len(L))
  lay <- lay[order(lay$row, lay$col), ]
  lay <- lay[!paste(lay$row, lay$col) %in% vapply(missing, paste, "", collapse = " "), ]
  lay$basket <- sprintf("T%d.B%d", lay$row, lay$col)
  rownames(lay) <- NULL
  lay
}

hn_grid <- function(scale, sigma = seq(0.05, 2.45, by = 0.1)) {
  lp <- dnorm(sigma, 0, scale, log = TRUE); list(value = sigma, logprior = lp - log(sum(exp(lp))))
}
twoway_prior <- function(p0, pbar, w_ex = 0.5, mu_sd = 10, mu_range = c(-8, 4), scale = 1 / sqrt(3),
                         sig_a = hn_grid(scale), sig_b = hn_grid(scale), sig_e = hn_grid(scale)) {
  list(w_ex = w_ex, mu_mean = qlogis(p0), mu_sd = mu_sd, mu_range = mu_range, scale = c(a = scale, b = scale, e = scale),
       sig_a = sig_a, sig_b = sig_b, sig_e = sig_e, pbar = pbar)
}
## Point-mass hyperpriors, for the reference comparisons with fixed standard deviations.
fix_sigmas <- function(prior, a, b, e) {
  pm <- function(v) list(value = v, logprior = 0)
  modifyList(prior, list(sig_a = pm(a), sig_b = pm(b), sig_e = pm(e), fixed = c(a = a, b = b, e = e)))
}

## Tables shared by all computational methods:
##   logD  [(n, y), (sigma_e, m)]  on the m lattice (step m_step), for the sampler;
##   DEXrb [(n, y), (sigma_e, m)]  EX marginal likelihood on the coarser Rao-Blackwell lattice;
##   KexT  [(sigma_e, m), cell]    EX prior of theta given (m, sigma_e) on the Rao-Blackwell lattice.
twoway_tables <- function(grid, prior, n_max, m_step = 0.05, rb_step = 0.1) {
  m <- seq(grid$range[1], grid$range[2], by = m_step); Gm <- length(m)
  rb <- seq(1, Gm, by = round(rb_step / m_step)); Grb <- length(rb)
  se <- prior$sig_e$value; Gs <- length(se)
  P <- do.call(cbind, lapply(se, function(s) norm_cells(grid, m, s)))       # cells x (Gs * Gm)
  V <- nex_prior(grid, prior$pbar)
  nn <- rep(0:n_max, times = 1:(n_max + 1)); yy <- sequence(1:(n_max + 1)) - 1
  Lik <- exp(loglik_grid(yy, nn, grid))
  DEX <- Lik %*% P; DN <- as.vector(Lik %*% V)
  rbcol <- as.vector(outer(rb, (seq_len(Gs) - 1) * Gm, "+"))                 # (sigma_e, m_rb) -> column of P
  list(m = m, m_step = m_step, Gm = Gm, rb_m = m[rb], rb_step = rb_step, Grb = Grb, Gs = Gs, n_max = n_max,
       logD = log(prior$w_ex * DEX + (1 - prior$w_ex) * DN), DEXrb = DEX[, rbcol, drop = FALSE], DN = DN,
       KexT = t(P[, rbcol, drop = FALSE]), V = V,
       row_of = function(y, n) { if (any(n > n_max)) stop("n exceeds n_max of the engine"); n * (n + 1) / 2 + y + 1 })
}

## log D_c(m, sigma_e) by linear interpolation in m (m clamped to the lattice).
lookup_logD <- function(tab, rho, s, m) {
  pos <- (pmin(pmax(m, tab$m[1]), tab$m[tab$Gm]) - tab$m[1]) / tab$m_step
  i0 <- pmin(floor(pos), tab$Gm - 2); f <- pos - i0; col <- (s - 1) * tab$Gm + i0 + 1
  rho <- as.vector(rho); col <- as.vector(col)
  (1 - f) * tab$logD[cbind(rho, col)] + f * tab$logD[cbind(rho, col + 1)]
}

## Posterior of theta_c on the grid from a distribution Q of (m_c, sigma_e) on the Rao-Blackwell
## lattice: Q is a sparse (pairs x (Gs * Grb)) matrix, rows summing to 1; pair i belongs to a
## basket with data (y[i], n[i]).
theta_from_q <- function(Q, y, n, tab, grid, w_ex) {
  tq <- triplets(Q); rho <- tab$row_of(y, n)[tq$i]
  dex <- tab$DEXrb[cbind(rho, tq$j)]; d <- w_ex * dex + (1 - w_ex) * tab$DN[rho]
  U <- Matrix::sparseMatrix(tq$i, tq$j, x = tq$x / d, dims = dim(Q))
  mix <- w_ex * as.matrix(U %*% tab$KexT) + (1 - w_ex) * outer(Matrix::rowSums(U), tab$V)
  pex <- as.vector(tapply(tq$x * w_ex * dex / d, factor(tq$i, seq_len(nrow(Q))), sum))
  list(logpost = loglik_grid(y, n, grid) + log(mix), pex = pex)
}

## Cloud-in-cell deposit of draws (m, sigma_e index) on the Rao-Blackwell lattice.
deposit <- function(pair, m, s, tab, n_pairs) {
  pos <- (pmin(pmax(m, tab$rb_m[1]), tab$rb_m[tab$Grb]) - tab$rb_m[1]) / tab$rb_step
  i0 <- pmin(floor(pos), tab$Grb - 2); f <- pos - i0; cell <- (s - 1) * tab$Grb + i0 + 1
  Q <- Matrix::sparseMatrix(c(pair, pair), c(cell, cell + 1), x = c(1 - f, f), dims = c(n_pairs, tab$Gs * tab$Grb))
  Matrix::Diagonal(x = 1 / Matrix::rowSums(Q)) %*% Q
}
triplets <- function(Q) { t <- methods::as(Q, "TsparseMatrix"); list(i = t@i + 1L, j = t@j + 1L, x = t@x) }

## Deterministic prior-predictive quadrature without data: m = mu + a + b, mu ~ TN,
## a + b ~ N(0, sa^2 + sb^2). Integrate mu on a 0.02 grid and put m-bin masses on the RB lattice.
empty_q <- function(tab, prior) {
  mu <- seq(prior$mu_range[1], prior$mu_range[2], by = 0.02)
  wmu <- dnorm(mu, prior$mu_mean, prior$mu_sd); wmu <- wmu / sum(wmu)
  v <- outer(prior$sig_a$value^2, prior$sig_b$value^2, "+"); pv <- exp(outer(prior$sig_a$logprior, prior$sig_b$logprior, "+"))
  pv <- tapply(pv, round(v, 10), sum); v <- as.numeric(names(pv))           # distinct variances of a + b
  edges <- c(-Inf, tab$rb_m[-1] - tab$rb_step / 2, Inf)
  qm <- numeric(tab$Grb)
  for (i in seq_along(v)) qm <- qm + pv[[i]] * as.vector(diff(pnorm(outer(edges, mu, "-") / sqrt(v[i]))) %*% wmu)
  as.vector(outer(qm, exp(prior$sig_e$logprior)))          # column order (sigma_e, m)
}

## Metropolis-within-Gibbs for (mu, a, b, sigma_a, sigma_b, sigma_e), vectorised over the rows
## of Y (states) and chains. Returns the draws of m_c (chains x baskets x keep) and of the
## sigma_e index (chains x keep); chain r belongs to state ((r - 1) %% S) + 1.
twoway_mcmc <- function(Y, N, layout, tab, prior, burn = 500, keep = 1000, chains = 2) {
  S <- nrow(Y); K <- ncol(Y); J <- max(layout$row); L <- max(layout$col); Rc <- S * chains
  st <- rep(seq_len(S), chains); ch <- rep(seq_len(chains), each = S)
  rho <- matrix(tab$row_of(Y, N), S, K)[st, , drop = FALSE]
  IR <- outer(layout$row, seq_len(J), "=="); IC <- outer(layout$col, seq_len(L), "==")
  sa <- prior$sig_a; sb <- prior$sig_b; se <- prior$sig_e; Gs <- length(se$value)
  lo <- prior$mu_range[1]; hi <- prior$mu_range[2]; m0 <- prior$mu_mean; s0 <- prior$mu_sd
  lp_mu <- function(x) ifelse(x >= lo & x <= hi, -(x - m0)^2 / (2 * s0^2), -Inf)
  gumbel <- function(lp) max.col(lp - log(-log(matrix(runif(length(lp)), nrow(lp)))))   # categorical draws
  start <- function(g, r) pmin(length(g$value), which.min(abs(g$value - 0.4)) + (r - 1) * 3)
  ## start: pooled rate (chain 1) and one unit higher (chain 2), random effects 0
  pooled <- qlogis((rowSums(Y) + 0.5) / (rowSums(N) + 1))[st]
  mu <- pmin(pmax(pooled + (ch - 1), lo + 0.01), hi - 0.01)
  a <- matrix(0, Rc, J); b <- matrix(0, Rc, L)
  ia <- start(sa, ch); ib <- start(sb, ch); ie <- start(se, ch)
  m <- mu + a[, layout$row, drop = FALSE] + b[, layout$col, drop = FALSE]
  lD <- matrix(lookup_logD(tab, rho, rep(ie, K), m), Rc, K)
  step <- matrix(0.3, Rc, 3); acc <- matrix(0, Rc, 3)
  M <- array(NA_real_, c(Rc, K, keep)); E <- matrix(NA_integer_, Rc, keep)
  mh_block <- function(cur, prop, inc, idx, sig) {            # a or b: one proposal per row/column
    mp <- m + (prop - cur)[, idx, drop = FALSE]
    lDp <- matrix(lookup_logD(tab, rho, rep(ie, K), mp), Rc, K)
    lr <- (lDp - lD) %*% inc + (cur^2 - prop^2) / (2 * sig^2)
    ok <- log(matrix(runif(length(lr)), Rc)) < lr
    okk <- ok[, idx, drop = FALSE]
    list(x = ifelse(ok, prop, cur), m = ifelse(okk, mp, m), lD = ifelse(okk, lDp, lD), rate = rowMeans(ok))
  }
  for (it in seq_len(burn + keep)) {
    ## 1. mu
    mup <- mu + step[, 1] * rnorm(Rc); mp <- m + (mup - mu)
    lDp <- matrix(lookup_logD(tab, rho, rep(ie, K), mp), Rc, K)
    ok <- log(runif(Rc)) < rowSums(lDp - lD) + lp_mu(mup) - lp_mu(mu)
    mu[ok] <- mup[ok]; m[ok, ] <- mp[ok, ]; lD[ok, ] <- lDp[ok, ]; acc[, 1] <- acc[, 1] + ok
    ## 2-3. row effects a_j and column effects b_l (rows, resp. columns, are independent given the rest)
    r <- mh_block(a, a + step[, 2] * matrix(rnorm(Rc * J), Rc), IR, layout$row, sa$value[ia])
    a <- r$x; m <- r$m; lD <- r$lD; acc[, 2] <- acc[, 2] + r$rate
    r <- mh_block(b, b + step[, 3] * matrix(rnorm(Rc * L), Rc), IC, layout$col, sb$value[ib])
    b <- r$x; m <- r$m; lD <- r$lD; acc[, 3] <- acc[, 3] + r$rate
    ## 4. exact Gibbs moves along (mu + k, a - k), (mu + k, b - k), (a + k, b - k): m is unchanged
    va <- sa$value[ia]^2; vb <- sb$value[ib]^2
    for (dir in 1:2) {
      n_e <- if (dir == 1) J else L; s_e <- if (dir == 1) rowSums(a) else rowSums(b); v_e <- if (dir == 1) va else vb
      P <- 1 / s0^2 + n_e / v_e; mk <- ((m0 - mu) / s0^2 + s_e / v_e) / P; sk <- 1 / sqrt(P)
      u <- runif(Rc, pnorm((lo - mu - mk) / sk), pnorm((hi - mu - mk) / sk))
      k <- mk + sk * qnorm(u); k[!is.finite(k)] <- 0
      mu <- mu + k; if (dir == 1) a <- a - k else b <- b - k
    }
    P <- J / va + L / vb; k <- (rowSums(b) / vb - rowSums(a) / va) / P + rnorm(Rc) / sqrt(P)
    a <- a + k; b <- b - k
    ## 5. sigma_a, sigma_b: exact draws on their grids
    ia <- gumbel(matrix(sa$logprior, Rc, length(sa$value), byrow = TRUE) - outer(rep(J, Rc), log(sa$value)) -
                 outer(rowSums(a^2), 1 / (2 * sa$value^2)))
    ib <- gumbel(matrix(sb$logprior, Rc, length(sb$value), byrow = TRUE) - outer(rep(L, Rc), log(sb$value)) -
                 outer(rowSums(b^2), 1 / (2 * sb$value^2)))
    ## 6. sigma_e: random walk on its grid
    if (Gs > 1) {
      iep <- ie + sample(c(-2L, -1L, 1L, 2L), Rc, TRUE); ok <- iep >= 1 & iep <= Gs; iep[!ok] <- ie[!ok]
      lDp <- matrix(lookup_logD(tab, rho, rep(iep, K), m), Rc, K)
      ok <- ok & log(runif(Rc)) < rowSums(lDp - lD) + se$logprior[iep] - se$logprior[ie]
      ie[ok] <- iep[ok]; lD[ok, ] <- lDp[ok, ]
    }
    ## adapt the random-walk steps during burn-in (target acceptance 0.35)
    if (it <= burn && it %% 50 == 0) { step <- pmin(pmax(step * exp(2 * (acc / 50 - 0.35)), 0.01), 5); acc[] <- 0 }
    if (it > burn) { M[, , it - burn] <- m; E[, it - burn] <- ie }
  }
  list(M = M, E = E, state = st, chain = ch, step = step)
}

## Potential scale reduction (Gelman-Rubin) of m_c per state and basket.
rhat_m <- function(dr, S) {
  K <- dim(dr$M)[2]; C <- max(dr$chain); out <- matrix(NA_real_, S, K)
  if (C < 2) return(out)
  for (k in seq_len(K)) {
    x <- dr$M[, k, ]; n <- ncol(x); mc <- rowMeans(x); vc <- rowSums((x - mc)^2) / (n - 1)
    W <- rowMeans(matrix(vc, S)); B <- n * apply(matrix(mc, S), 1, var)
    out[, k] <- sqrt(((n - 1) / n * W + B / n) / W)
  }
  out
}

engine_twoway <- function(grid, layout, pbar, n_max, prior = twoway_prior(grid$p0, pbar),
                          mcmc = list(burn = 500, keep = 1000, chains = 2, chunk = 150), tab = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("the two-way engine needs the Matrix package")
  if (is.null(tab)) tab <- twoway_tables(grid, prior, n_max)
  q0 <- empty_q(tab, prior)
  list(name = "two-way EXNEX", grid = grid, layout = layout, tab = tab, prior = prior,
       posterior = function(Y, N, need = NULL, seed = 1) {
    S <- nrow(Y); K <- ncol(Y); G <- length(grid$theta)
    if (K != nrow(layout)) stop("number of baskets differs from the layout")
    if (is.null(need)) need <- matrix(TRUE, S, K)
    lp <- array(NA_real_, c(S, K, G)); pex <- matrix(NA_real_, S, K); rhat <- matrix(NA_real_, S, K)
    empty <- rowSums(N) == 0
    for (i in which(empty & rowSums(need) > 0)) {                         # prior predictive, deterministic quadrature
      k <- which(need[i, ])
      Q <- Matrix::Matrix(matrix(q0, length(k), length(q0), byrow = TRUE), sparse = TRUE)
      th <- theta_from_q(Q, Y[i, k], N[i, k], tab, grid, prior$w_ex)
      lp[i, k, ] <- th$logpost; pex[i, k] <- th$pex
    }
    todo <- which(!empty & rowSums(need) > 0)
    if (length(todo)) set.seed(seed)
    for (ids in split(todo, ceiling(seq_along(todo) / mcmc$chunk))) {
      dr <- twoway_mcmc(Y[ids, , drop = FALSE], N[ids, , drop = FALSE], layout, tab, prior, mcmc$burn, mcmc$keep, mcmc$chains)
      rhat[ids, ] <- rhat_m(dr, length(ids))
      nd <- need[ids, , drop = FALSE]; pid <- matrix(0L, length(ids), K); pid[nd] <- seq_len(sum(nd))
      sel <- nd[dr$state, , drop = FALSE]                                  # chain rows x baskets that are needed
      kk <- col(sel)[sel]; rr <- row(sel)[sel]
      pair <- rep(pid[cbind(dr$state[rr], kk)], times = mcmc$keep)
      mm <- dr$M[cbind(rep(rr, times = mcmc$keep), rep(kk, times = mcmc$keep), rep(seq_len(mcmc$keep), each = length(rr)))]
      ss <- dr$E[cbind(rep(rr, times = mcmc$keep), rep(seq_len(mcmc$keep), each = length(rr)))]
      Q <- deposit(pair, mm, ss, tab, sum(nd))
      w <- which(nd, arr.ind = TRUE); w <- w[order(pid[w]), , drop = FALSE]   # pair id -> (state, basket)
      th <- theta_from_q(Q, Y[ids, , drop = FALSE][w], N[ids, , drop = FALSE][w], tab, grid, prior$w_ex)
      for (p in seq_len(nrow(w))) lp[ids[w[p, 1]], w[p, 2], ] <- th$logpost[p, ]
      pex[cbind(ids[w[, 1]], w[, 2])] <- th$pex
    }
    list(logpost = lp, pex = pex, rhat = rhat)
  })
}
