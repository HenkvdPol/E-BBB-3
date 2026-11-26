## =============================================================================
##  borrow_oneway.R -- borrowing engines without structure: none and one-way EXNEX
##
##  Every engine is a list with
##    $posterior(Y, N, need, seed)  Y, N: (states x baskets) responses and patients.
##        Returns  logpost: [states, baskets, cells] log posterior of theta_c given ALL data of
##                          the state (NA where need is FALSE), and
##                 pex:     P(basket c is exchangeable | data) (NA without borrowing).
##  The e-process only ever sees restrict(logpost) (eprocess.R).
##
##  One-way EXNEX (Neuenschwander et al. 2016), as in legacy/additional_functions.R:
##    EX  (prob w):     theta_c ~ N(mu, tau^2),  mu ~ N(logit p0, 10^2) on [-8, 4],  tau ~ HN(1)
##    NEX (prob 1 - w): theta_c ~ N(logit pbar, 1/pbar + 1/(1 - pbar)),  pbar = (p0 + p1) / 2
##  Given (mu, tau) the baskets are independent, so the posterior is a sum over a (mu, tau)
##  grid of one-dimensional integrals on the theta grid (deterministic, no sampling).
## =============================================================================

## Log likelihood of every state on the grid: (states x cells), no binomial coefficient.
loglik_grid <- function(y, n, grid) outer(y, log(grid$p)) + outer(n - y, log(1 - grid$p))

nex_prior <- function(grid, pbar) as.vector(norm_cells(grid, qlogis(pbar), sqrt(1 / pbar + 1 / (1 - pbar))))

## Table with one row per data state (n, y) of a basket. Rows for n <= n_max are built at once; rows
## for larger n are added the first time they are needed, so the engines put no limit on a basket's
## sample size (open-ended designs). compute(y, n) returns a named list of matrices (one row per
## state) or vectors (one element per state); they become fields of the table. tab$row_of(y, n)
## returns the row of each state, building missing rows first. Call row_of() BEFORE reading a field:
## in tab$LE[tab$row_of(y, n), ] R would fetch the old LE before the new rows are added.
state_table <- function(compute, n_max) {
  tab <- new.env()
  key <- function(y, n) n * (n + 1) / 2 + y + 1                                    # unique number of (n, y)
  nn <- rep(0:n_max, times = 1:(n_max + 1)); yy <- sequence(1:(n_max + 1)) - 1
  first <- compute(yy, nn); for (nm in names(first)) tab[[nm]] <- first[[nm]]
  tab$fields <- names(first); tab$pos <- seq_along(yy)                             # key -> row
  tab$row_of <- function(y, n) {
    if (any(y < 0 | y > n)) stop("impossible data state (y > n or y < 0)")
    k <- key(as.vector(y), as.vector(n))
    if (length(k) && max(k) > length(tab$pos)) tab$pos[max(k)] <- NA
    new <- which(is.na(tab$pos[k]) & !duplicated(k))
    if (length(new)) {
      add <- compute(as.vector(y)[new], as.vector(n)[new]); r0 <- NROW(tab[[tab$fields[1]]])
      for (nm in tab$fields) tab[[nm]] <- if (is.matrix(tab[[nm]])) rbind(tab[[nm]], add[[nm]]) else c(tab[[nm]], add[[nm]])
      tab$pos[k[new]] <- r0 + seq_along(new)
    }
    tab$pos[k]
  }
  tab
}

## No borrowing: theta_c ~ NEX prior, only the basket's own data.
engine_none <- function(grid, pbar) {
  logV <- log(nex_prior(grid, pbar))
  list(name = "none", grid = grid, posterior = function(Y, N, need = NULL, seed = NULL) {
    lp <- array(NA_real_, c(nrow(Y), ncol(Y), length(grid$theta)))
    for (k in seq_len(ncol(Y))) lp[, k, ] <- loglik_grid(Y[, k], N[, k], grid) + rep(logV, each = nrow(Y))
    list(logpost = lp, pex = matrix(NA_real_, nrow(Y), ncol(Y)))
  })
}

## The hyperparameter grid of the legacy benchmark (legacy/Setting.R, section 6).
oneway_hyper <- function(p0, mu_sd = 10, tau_scale = 1) {
  h <- expand.grid(mu = seq(-8, 4, by = 0.15), tau = seq(0.05, 3.95, by = 0.10))
  h$logprior <- dnorm(h$mu, qlogis(p0), mu_sd, log = TRUE) + dnorm(h$tau, 0, tau_scale, log = TRUE)
  h
}

## kernel = "cell": EX prior as cell probabilities (default); "density": normal density at the
## cell midpoints, renormalised -- the legacy choice, kept for the regression test only.
## n_max only sets how much of the (n, y) table is built at once; larger n are added when needed.
engine_oneway <- function(grid, pbar, n_max, w_ex = 0.5, hyper = oneway_hyper(grid$p0), kernel = "cell") {
  W <- if (kernel == "cell") norm_cells(grid, hyper$mu, hyper$tau) else {
    D <- dnorm(outer(grid$theta, hyper$mu, "-") / rep(hyper$tau, each = length(grid$theta))); sweep(D, 2, colSums(D), "/") }
  V <- if (kernel == "cell") nex_prior(grid, pbar) else {
    v <- dnorm(grid$theta, qlogis(pbar), sqrt(1 / pbar + 1 / (1 - pbar))); v / sum(v) }
  tab <- state_table(function(y, n) {                                                  # marginal likelihoods per (n, y)
    Lik <- exp(loglik_grid(y, n, grid)); list(LE = Lik %*% W, LN = as.vector(Lik %*% V)) }, n_max)
  list(name = "one-way EXNEX", grid = grid, W = W, V = V, hyper = hyper, tab = tab,
       posterior = function(Y, N, need = NULL, seed = NULL) {
    S <- nrow(Y); K <- ncol(Y); G <- length(grid$theta)
    if (is.null(need)) need <- matrix(TRUE, S, K)
    r <- matrix(tab$row_of(Y, N), S, K); LE <- tab$LE; LN <- tab$LN
    D <- lapply(seq_len(K), function(k) w_ex * LE[r[, k], , drop = FALSE] + (1 - w_ex) * LN[r[, k]])
    lw <- matrix(hyper$logprior, S, nrow(hyper), byrow = TRUE) + Reduce(`+`, lapply(D, log))
    w <- exp(lw - row_max(lw)); w <- w / rowSums(w)                                    # posterior of (mu, tau)
    lp <- array(NA_real_, c(S, K, G)); pex <- matrix(NA_real_, S, K)
    for (k in seq_len(K)) {
      i <- which(need[, k]); if (!length(i)) next
      u <- w[i, , drop = FALSE] / D[[k]][i, , drop = FALSE]                            # w_h / D_kh
      mix <- w_ex * (u %*% t(W)) + (1 - w_ex) * outer(rowSums(u), V)                  # prior of theta_k given the rest
      lp[i, k, ] <- loglik_grid(Y[i, k], N[i, k], grid) + log(mix)
      pex[i, k] <- w_ex * rowSums(u * LE[r[i, k], , drop = FALSE])
    }
    list(logpost = lp, pex = pex)
  })
}
