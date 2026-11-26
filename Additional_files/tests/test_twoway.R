## Two-way EXNEX (additional_functions.R, Part 2): tables, indexing, missing/empty cells,
## reproducibility, and agreement with the deterministic reference (validation/reference_functions.R).
p0 <- 0.1; pbar <- 0.2; grid <- theta_grid(p0); prior <- twoway_prior(p0, pbar)
tab <- twoway_tables(grid, prior, 30)
pmean <- function(lp) { x <- exp(lp - max(lp)); sum(x * grid$p) / sum(x) }
means <- function(po) apply(po$logpost[1, , , drop = FALSE], 2, function(v) pmean(as.vector(v)))

## 1. tables
check("tables: no data gives log D = 0", max(abs(tab$logD[tab$row_of(0, 0), ])) < 1e-12)
set.seed(4); worst <- 0
for (i in 1:20) {
  n <- sample(1:30, 1); y <- sample(0:n, 1); s <- sample(tab$Gs, 1); j <- sample(tab$Gm, 1)
  L <- exp(loglik_grid(y, n, grid))[1, ]; P <- norm_cells(grid, tab$m[j], prior$sig_e$value[s])[, 1]
  direct <- log(prior$w_ex * sum(L * P) + (1 - prior$w_ex) * sum(L * tab$V))
  worst <- max(worst, abs(tab$logD[tab$row_of(y, n), (s - 1) * tab$Gm + j] - direct))
}
check("tables: log D = direct sum over the theta cells", worst < 1e-10, paste("max diff", fmt(worst)))
r <- tab$row_of(7, 20); mid <- (tab$m[100] + tab$m[101]) / 2
check("lookup: exact at lattice points and linear in between",
      abs(lookup_logD(tab, r, 3, tab$m[100]) - tab$logD[r, 2 * tab$Gm + 100]) < 1e-12 &&
      abs(lookup_logD(tab, r, 3, mid) - mean(tab$logD[r, 2 * tab$Gm + 100:101])) < 1e-12)

## 2. exact prior predictive of (m, sigma_e) without data, against Monte Carlo draws
q0 <- empty_q(tab, prior); set.seed(5); M <- 4e5
sa <- sample(prior$sig_a$value, M, TRUE, exp(prior$sig_a$logprior)); sb <- sample(prior$sig_b$value, M, TRUE, exp(prior$sig_b$logprior))
mu <- qnorm(runif(M, pnorm(-8, prior$mu_mean, 10), pnorm(4, prior$mu_mean, 10)), prior$mu_mean, 10)
m <- pmin(pmax(mu + rnorm(M, 0, sa) + rnorm(M, 0, sb), -10), 6)
hist_mc <- tabulate(round((m + 10) / tab$rb_step) + 1, tab$Grb) / M
tv <- 0.5 * sum(abs(rowSums(matrix(q0, tab$Grb)) - hist_mc))
check("empty platform: exact distribution of m sums to 1", abs(sum(q0) - 1) < 1e-10)
check("empty platform: exact distribution of m agrees with Monte Carlo", tv < 0.01, paste("TV", fmt(tv)))

## 3. a missing cell and an open cell without patients give the same posterior elsewhere
l9 <- grid_layout(3, 3); l8 <- grid_layout(3, 3, missing = list(c(3, 3)))
y9 <- c(9, 9, 9, 3, 3, 3, 3, 3, 0); n9 <- c(rep(30, 8), 0)
e9 <- engine_twoway(grid, l9, pbar, 30, prior = prior, tab = tab); e8 <- engine_twoway(grid, l8, pbar, 30, prior = prior, tab = tab)
p9 <- e9$posterior(matrix(y9, 1), matrix(n9, 1), seed = 3); p8 <- e8$posterior(matrix(y9[1:8], 1), matrix(n9[1:8], 1), seed = 3)
check("missing cell = empty cell for the other baskets", max(abs(p9$logpost[1, 1:8, ] - p8$logpost[1, , ])) < 1e-9)

## 4. reproducibility
p9b <- e9$posterior(matrix(y9, 1), matrix(n9, 1), seed = 3); p9c <- e9$posterior(matrix(y9, 1), matrix(n9, 1), seed = 4)
check("same seed gives identical posteriors", identical(p9$logpost, p9b$logpost))
d <- max(abs(means(p9) - means(p9c)))
check("another seed changes posterior means by < 0.01 (Monte Carlo error)", d < 0.01, paste("max diff", fmt(d)))

## 5. row/column indexing: transposing the layout and the data transposes the posterior
l23 <- grid_layout(2, 3); l32 <- grid_layout(3, 2)
y23 <- c(9, 3, 1, 6, 2, 4); n23 <- c(30, 20, 10, 30, 25, 15)            # (1,1) (1,2) (1,3) (2,1) (2,2) (2,3)
tr <- match(paste(l23$col, l23$row), paste(l32$row, l32$col))          # basket of l32 that is the transpose
e23 <- engine_twoway(grid, l23, pbar, 30, prior = prior, tab = tab); e32 <- engine_twoway(grid, l32, pbar, 30, prior = prior, tab = tab)
yt <- numeric(6); nt <- numeric(6); yt[tr] <- y23; nt[tr] <- n23
m23 <- means(e23$posterior(matrix(y23, 1), matrix(n23, 1), seed = 1)); m32 <- means(e32$posterior(matrix(yt, 1), matrix(nt, 1), seed = 1))
d <- max(abs(m23 - m32[tr]))
check("transposed layout and data give transposed posteriors (within Monte Carlo error)", d < 0.01, paste("max diff", fmt(d)))

## 6. structure: an active row is borrowed within the row; null rows are protected relative to one-way EXNEX
po2 <- means(p9); po1 <- engine_oneway(grid, pbar, 30)$posterior(matrix(y9, 1), matrix(n9, 1))
po1 <- apply(po1$logpost[1, , , drop = FALSE], 2, function(v) pmean(as.vector(v)))
check("active row: its baskets have the highest posterior means", min(po2[1:3]) > max(po2[4:8]), fmt(po2))
check("null baskets outside the active row: lower posterior mean than one-way EXNEX", all(po2[4:8] < po1[4:8]),
      paste("two-way", fmt(mean(po2[4:8])), "one-way", fmt(mean(po1[4:8]))))

## 7. priors handed to the e-process are valid; E_p0[factor] = 1
w <- restrict(matrix(p9$logpost[1, , ], 9), grid)
ef <- sapply(1:9, function(k) sum(dbinom(0:5, 5, p0) * exp(log_factor(w[rep(k, 6), ], 0:5, rep(5, 6), grid))))
check("two-way priors: valid and E_p0[factor] = 1", isTRUE(check_priors(w, grid)) && max(abs(ef - 1)) < 1e-12)

## 8. agreement with the deterministic reference (2 x 2, fixed sigmas: exact quadrature is fast)
pf <- fix_sigmas(prior, 0.45, 0.45, 0.45); tf <- twoway_tables(grid, pf, 30); l22 <- grid_layout(2, 2)
y4 <- c(9, 3, 4, 0); n4 <- c(30, 30, 20, 5)
ex <- twoway_exact(y4, n4, l22, tf, pf, grid, c_half = 3.5, a_max = 2.5)
mc <- engine_twoway(grid, l22, pbar, 30, prior = pf, tab = tf)$posterior(matrix(y4, 1), matrix(n4, 1), seed = 2)
a <- summ_grid(ex$logpost, ex$pex, grid); b <- summ_grid(matrix(mc$logpost[1, , ], 4), mc$pex[1, ], grid)
check("2 x 2, fixed sigmas: restricted mean within 0.003 of the quadrature", max(abs(a$rmean - b$rmean)) < 0.003, fmt(max(abs(a$rmean - b$rmean))))
check("2 x 2, fixed sigmas: block factor (5/5) within 5% of the quadrature", max(abs(log(a$f5 / b$f5))) < 0.05, fmt(max(abs(log(a$f5 / b$f5)))))
check("2 x 2, fixed sigmas: P(EX) within 0.05 of the quadrature (production sampler: Monte Carlo error)", max(abs(a$pex - b$pex)) < 0.05,
      fmt(max(abs(a$pex - b$pex))))

## 9. slow: L x 1 layout with full hyperpriors against the equivalent one-way model, and JAGS
if (!SLOW) skip("3 x 1 layout vs equivalent one-way model; JAGS", "slow test (use --slow)") else {
  l31 <- grid_layout(3, 1); y3 <- c(3, 9, 14); n3 <- rep(30, 3)
  ow <- engine_oneway(grid, pbar, 30, hyper = oneway_equivalent_hyper(prior))$posterior(matrix(y3, 1), matrix(n3, 1))
  mc <- engine_twoway(grid, l31, pbar, 30, prior = prior, tab = tab)$posterior(matrix(y3, 1), matrix(n3, 1), seed = 1)
  a <- summ_grid(matrix(ow$logpost[1, , ], 3), ow$pex[1, ], grid); b <- summ_grid(matrix(mc$logpost[1, , ], 3), mc$pex[1, ], grid)
  check("3 x 1: restricted mean within 0.003 of the one-way equivalent", max(abs(a$rmean - b$rmean)) < 0.003, fmt(max(abs(a$rmean - b$rmean))))
  check("3 x 1: block factor (5/5) within 5% of the one-way equivalent", max(abs(log(a$f5 / b$f5))) < 0.05, fmt(max(abs(log(a$f5 / b$f5)))))
  check("3 x 1: P(EX) within 0.05 of the one-way equivalent", max(abs(a$pex - b$pex)) < 0.05, fmt(max(abs(a$pex - b$pex))))
  if (!requireNamespace("rjags", quietly = TRUE)) skip("2 x 2 vs JAGS", "rjags not installed") else {
    jg <- twoway_jags(y4, n4, l22, pf, grid, n_iter = 20000, n_burn = 3000); j <- summ_draws(jg$theta, jg$d, grid)
    ref <- summ_grid(ex$logpost, ex$pex, grid); dj <- max(abs(ref$rmean - j$rmean))
    check("2 x 2, fixed sigmas: JAGS restricted mean within 0.003 of the quadrature", dj < 0.003, fmt(dj))
  }
}
