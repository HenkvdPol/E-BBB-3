## Open-ended baskets: futility layer, growing tables, telescoping, exact DP, simulator.
p0 <- 0.1; p1 <- 0.3; pbar <- 0.2; grid <- theta_grid(p0); fut <- futility_side(grid, p1)
lpV <- matrix(log(nex_prior(grid, pbar)), 1)
w0 <- as.vector(restrict(lpV, grid)); w0f <- as.vector(restrict_below(lpV, fut))

## 1. futility layer: support below p1, E_p1[factor] = 1, <= 1 above p1, invalid priors rejected
set.seed(21); lp <- matrix(rnorm(4 * length(grid$theta), sd = 3), 4); wf <- restrict_below(lp, fut)
check("futility support lies in [0, p1)", all(fut$support < p1) && length(fut$support) == sum(grid$p < p1))
mart_f <- function(w, n, p) vapply(seq_len(nrow(w)), function(i)
  sum(dbinom(0:n, n, p) * exp(log_factor(w[rep(i, n + 1), , drop = FALSE], 0:n, rep(n, n + 1), fut))), 0)
m1 <- mart_f(wf, 7, p1); m2 <- c(mart_f(wf, 7, 0.4), mart_f(wf, 7, 0.6))
check("E_p1[futility factor] = 1", max(abs(m1 - 1)) < 1e-12, fmt(max(abs(m1 - 1))))
check("E_p[futility factor] <= 1 for p > p1", all(m2 <= 1), fmt(max(m2)))
check_error("futility prior with mass at p1 is rejected", { f2 <- fut; f2$support[length(f2$support)] <- p1; check_priors_below(wf, f2) })
check_error("unnormalised futility prior is rejected", check_priors_below(wf * 1.01, fut))

## 2. tables grow beyond n_max and agree with tables built at once
small <- engine_oneway(grid, pbar, n_max = 10); big <- engine_oneway(grid, pbar, n_max = 40)
Y <- matrix(c(3, 12, 0, 25, 7, 1), 2); N <- matrix(c(10, 40, 5, 40, 30, 3), 2)
a <- small$posterior(Y, N)$logpost; b <- big$posterior(Y, N)$logpost
check("one-way: rows beyond n_max are added on demand and equal a full table", max(abs(a - b)) < 1e-10)
tw_small <- twoway_tables(grid, twoway_prior(p0, pbar), 5); tw_big <- twoway_tables(grid, twoway_prior(p0, pbar), 12)
rs <- tw_small$row_of(c(4, 9, 0), c(12, 11, 7)); rb <- tw_big$row_of(c(4, 9, 0), c(12, 11, 7))   # rows first, then read
check("two-way: rows beyond n_max are added on demand and equal a full table",
      max(abs(tw_small$logD[rs, ] - tw_big$logD[rb, ])) < 1e-12 && max(abs(tw_small$DEXrb[rs, ] - tw_big$DEXrb[rb, ])) < 1e-12)
check_error("impossible data state is rejected", small$tab$row_of(5, 3))

## 3. telescoping: without borrowing the simulated e-values equal the fixed mixtures at (y, n)
design <- list(layout = grid_layout(2, 2), n_max = Inf, n_cap = 60, rate = 1, entry = 0, look_every = 5, alpha = 0.05,
               futility = "evalue", p1 = p1, beta = 0.2)
d <- simulate_trials(c(0.1, 0.2, 0.3, 0.45), design, engine_none(grid, pbar), 6, seed = 3)
tE <- vapply(seq_len(nrow(d)), function(i) mixture_logE(w0, grid, d$y[i], d$n[i]), 0)
tF <- vapply(seq_len(nrow(d)), function(i) mixture_logE(w0f, fut, d$y[i], d$n[i]), 0)
check("telescoping: simulated log E = mixture log E at (y, n)", max(abs(d$logE - tE)) < 1e-9, fmt(max(abs(d$logE - tE))))
check("telescoping: simulated log E^f = mixture log E^f at (y, n)", max(abs(d$logEf - tF)) < 1e-9, fmt(max(abs(d$logEf - tF))))

## 4. the exact DP equals brute-force enumeration of the block-by-block recursion (3 blocks of 5)
e1 <- engine_none(grid, pbar)
enumerate <- function(p) {
  paths <- as.matrix(expand.grid(0:5, 0:5, 0:5)); pr <- 0; pf <- 0
  for (i in seq_len(nrow(paths))) {
    y <- 0; n <- 0; lE <- 0; lF <- 0; done <- FALSE
    for (k in 1:3) {
      lp <- matrix(e1$posterior(matrix(y, 1), matrix(n, 1))$logpost[1, 1, ], 1)        # prior frozen before block k
      lE <- lE + log_factor(restrict(lp, grid), paths[i, k], 5, grid); lF <- lF + log_factor(restrict_below(lp, fut), paths[i, k], 5, fut)
      y <- y + paths[i, k]; n <- n + 5
      if (lE >= log(20)) { pr <- pr + prod(dbinom(paths[i, ], 5, p)); done <- TRUE; break }
      if (lF >= log(5)) { pf <- pf + prod(dbinom(paths[i, ], 5, p)); done <- TRUE; break }
    }
  }
  c(pr, pf)
}
for (p in c(0.1, 0.3)) {
  dp <- exact_basket(p, w0, grid, 0.05, block = 5, n_max = Inf, horizon = 15, futility = "evalue", w0f = w0f, fut = fut, beta = 0.2)$summary
  en <- enumerate(p)
  check(sprintf("exact DP = path enumeration of the recursion (p = %.1f)", p),
        abs(dp$reject - en[1]) < 1e-12 && abs(dp$futility - en[2]) < 1e-12, fmt(c(dp$reject, en[1], dp$futility, en[2])))
}

## 5. exact guarantees over a long horizon (no borrowing)
o0 <- exact_basket(p0, w0, grid, 0.05, block = 5, horizon = 1000)$summary
f0 <- exact_basket(p0, w0, grid, 0.05, block = 5, horizon = 1000, futility = "evalue", w0f = w0f, fut = fut, beta = 0.2)$summary
f1 <- exact_basket(p1, w0, grid, 0.05, block = 5, horizon = 1000, futility = "evalue", w0f = w0f, fut = fut, beta = 0.2)$summary
check("open-ended, no futility: P(reject by n = 1000 | p0) <= alpha", o0$reject <= 0.05, fmt(o0$reject))
check("open-ended, no futility: a null basket usually never stops", o0$open > 0.9, fmt(o0$open))
check("open-ended with e-value futility: P(futility | p1) <= beta", f1$futility <= 0.2, fmt(f1$futility))
check("open-ended with e-value futility: every basket stops (p0 and p1)", f0$open < 1e-12 && f1$open < 1e-12, fmt(c(f0$last_n, f1$last_n)))

## 6. simulator: open-ended statuses, the simulation horizon, and the required n_cap
d2 <- simulate_trials(c(0.1, 0.1, 0.1, 0.1), modifyList(design, list(futility = "none", n_cap = 40)), e1, 4, seed = 9)
check("simulator: without futility, null baskets are open at the horizon or rejected", all(d2$reject | d2$open) && any(d2$open))
check("simulator: open baskets have exactly n_cap patients", all(d2$n[d2$open] == 40) && all(d2$n <= 40))
check("simulator: every basket rejected, stopped or open (e-value futility)", all(d$reject | d$stopped | d$open))
check_error("simulator: n_max = Inf without n_cap is refused", simulate_trials(rep(0.1, 4), modifyList(design, list(n_cap = NULL)), e1, 2))
d3 <- simulate_trials(c(0.1, 0.2, 0.3, 0.45), design, engine_oneway(grid, pbar, 10), 3, seed = 4)
check("simulator: open-ended with one-way borrowing runs beyond the engine's initial table", max(d3$n) > 10 && all(d3$reject | d3$stopped | d3$open))
