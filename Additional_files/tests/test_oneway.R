## One-way EXNEX and no borrowing (additional_functions.R), with the regression test against Additional_files/legacy.
p0 <- 0.1; pbar <- 0.2; grid <- theta_grid(p0)
set.seed(3); K <- 5; S <- 12
N <- matrix(sample(c(0, 5, 10, 20, 30), S * K, TRUE), S); Y <- matrix(rbinom(S * K, N, 0.25), S)
post_mean <- function(lp, g) { x <- exp(lp - max(lp)); sum(x * g$p) / sum(x) }

## 1. regression: the port reproduces legacy posterior_exnex() on the legacy grid and kernel
legacy <- file.path(root, "Additional_files", "legacy", "additional_functions.R")
if (!file.exists(legacy)) skip("one-way port equals legacy engine", "Additional_files/legacy not found") else {
  old <- new.env(); sys.source(legacy, envir = old)
  s <- list(p0 = p0, n_max = 30, n_ext = 0, futility = list(p_threshold = pbar),
            prior = list(mu_mean = qlogis(p0), mu_sd = 10, tau_scale = 1, nex_p = pbar),
            grid = list(theta_min = -10, theta_max = 6, theta_step = 0.02, mu_min = -8, mu_max = 4, mu_step = 0.15, tau_max = 4, tau_step = 0.10))
  e_old <- old$build_engine(s)
  ref <- old$posterior_exnex(e_old, Y, N, matrix(0.5, S, K), moments = TRUE)
  th <- seq(-10, 6, by = 0.02); lg <- list(theta = th, p = plogis(th), p0 = p0)
  e_new <- engine_oneway(lg, pbar, n_max = 30, kernel = "density")
  po <- e_new$posterior(Y, N)
  above <- pmin(pmax((th + 0.01 - qlogis(p0)) / 0.02, 0), 1)             # legacy's partial cell at p0
  mom <- function(i, k, f) { x <- exp(po$logpost[i, k, ] - max(po$logpost[i, k, ])); sum(x * above * f) / sum(x * above) }
  d1 <- max(abs(outer(1:S, 1:K, Vectorize(function(i, k) mom(i, k, lg$p))) - ref$rmean))
  d2 <- max(abs(outer(1:S, 1:K, Vectorize(function(i, k) mom(i, k, lg$p^2))) - ref$rsec))
  pe <- outer(1:S, 1:K, Vectorize(function(i, k) { x <- exp(po$logpost[i, k, ] - max(po$logpost[i, k, ])); sum(x * above) / sum(x) }))
  d3 <- max(abs(pe - ref$pp_eff))
  check("one-way port = legacy: restricted mean", d1 < 1e-10, paste("max diff", fmt(d1)))
  check("one-way port = legacy: restricted second moment", d2 < 1e-10, paste("max diff", fmt(d2)))
  check("one-way port = legacy: P(p > p0 | data)", d3 < 1e-10, paste("max diff", fmt(d3)))
  ## 2. the new default (aligned grid of width 0.05, cell probabilities) against legacy
  e_def <- engine_oneway(grid, pbar, n_max = 30); pd <- e_def$posterior(Y, N)
  rm_def <- outer(1:S, 1:K, Vectorize(function(i, k) { w <- restrict(matrix(pd$logpost[i, k, ], 1), grid); sum(w * grid$support) }))
  d4 <- max(abs(rm_def - ref$rmean))
  check("default grid vs legacy: restricted mean within 0.003", d4 < 0.003, paste("max diff", fmt(d4)))
}

## 3. no borrowing = one-way with EX weight 0
e1 <- engine_oneway(grid, pbar, n_max = 30); en <- engine_none(grid, pbar); ez <- engine_oneway(grid, pbar, n_max = 30, w_ex = 0)
a <- en$posterior(Y, N)$logpost; b <- ez$posterior(Y, N)$logpost
nz <- function(x) sweep(x, c(1, 2), apply(x, c(1, 2), max))     # log posterior with max 0 per basket
check("no borrowing = one-way EXNEX with w = 0", max(abs(nz(a) - nz(b))) < 1e-9)

## 4. exchangeability: permuting the baskets permutes the posteriors
pm <- c(3, 1, 5, 2, 4); p1 <- e1$posterior(Y, N); p2 <- e1$posterior(Y[, pm], N[, pm])
check("one-way: permutation of baskets permutes posteriors", max(abs(nz(p1$logpost)[, pm, ] - nz(p2$logpost))) < 1e-9)

## 5. no data: the same prior predictive for every basket; deterministic
p0d <- e1$posterior(matrix(0, 1, K), matrix(0, 1, K))
check("one-way: identical prior predictive without data", max(abs(nz(p0d$logpost)[1, 1, ] - nz(p0d$logpost)[1, K, ])) < 1e-12)
check("one-way: deterministic (two calls identical)", identical(e1$posterior(Y, N)$logpost, p1$logpost))

## 6. 'need': only the requested baskets are computed
nd <- matrix(FALSE, S, K); nd[2, 3] <- TRUE
pn <- e1$posterior(Y, N, nd)
check("one-way: need selects the baskets", all(is.na(pn$logpost[-2, , ])) && all(is.na(pn$logpost[2, -3, ])) &&
        max(abs(pn$logpost[2, 3, ] - p1$logpost[2, 3, ])) < 1e-12)

## 7. borrowing behaviour: a null basket among active ones is pulled up relative to no borrowing
Yb <- matrix(c(3, 10, 10, 10, 10), 1); Nb <- matrix(30, 1, 5)
check("one-way: borrowing pulls a null basket towards the others",
      post_mean(e1$posterior(Yb, Nb)$logpost[1, 1, ], grid) > post_mean(en$posterior(Yb, Nb)$logpost[1, 1, ], grid))
