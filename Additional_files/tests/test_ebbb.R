## The priors of eBBB (engine_ebbb in additional_functions.R, Part 2) and the settings of the study.
## The closed forms below are written out here from the definitions (they are not eBBB's code).
p0 <- 0.1; p1 <- 0.3; grid <- theta_grid(p0); lay <- grid_layout(4, 3); K <- nrow(lay)

## Closed-form block factor of a Beta(a, b) prior restricted to (p0, 1]
beta_factor <- function(a, b, y, n) {
  exp(lbeta(a + y, b + n - y) - lbeta(a, b) +
      pbeta(p0, a + y, b + n - y, lower.tail = FALSE, log.p = TRUE) - pbeta(p0, a, b, lower.tail = FALSE, log.p = TRUE) -
      y * log(p0) - (n - y) * log(1 - p0))
}
grid_factor <- function(lp, y, n, g = grid) exp(log_factor(restrict(matrix(lp, 1), g), y, n, g))

## Beta on the grid
a <- c(1, 1, 3, 10, 0.5 + 15); b <- c(1, 20, 8, 30, 15)
check("Beta on the grid: probabilities sum to 1", max(abs(rowSums(exp(beta_cells(grid, a, b))) - 1)) < 1e-3,
      fmt(max(abs(rowSums(exp(beta_cells(grid, a, b))) - 1))))
rel <- max(sapply(seq_along(a), function(i) sapply(list(c(0, 5), c(2, 5), c(5, 5)), function(bl)
  abs(grid_factor(beta_cells(grid, a[i], b[i]), bl[1], bl[2]) / beta_factor(a[i], b[i], bl[1], bl[2]) - 1))))
check("Beta on the grid: block factors within 1e-3 of the closed form", rel < 1e-3, fmt(rel))

## Data on the 4 x 3 grid
set.seed(4); S <- 30
N <- matrix(sample(c(0, 5, 10, 20, 30), S * K, TRUE), S); Y <- matrix(rbinom(S * K, N, 0.25), S)
same_row <- outer(lay$row, lay$row, "==") & !diag(K); same_col <- outer(lay$col, lay$col, "==") & !diag(K)

## No borrowing and power priors: the Beta parameters of their definitions
nb <- engine_ebbb(grid, lay, "no borrowing")$posterior(Y, N)
check("no borrowing: Beta(1 + y, 1 + n - y)", isTRUE(all.equal(nb$logpost[, 3, ], beta_cells(grid, 1 + Y[, 3], 1 + N[, 3] - Y[, 3]))))
pe <- engine_ebbb(grid, lay, "power", w_row = 0.5, w_col = 0.5)$posterior(Y, N)
k <- 5; rr <- same_row[k, ]; cc <- same_col[k, ]
ya <- as.vector(1 + Y[, k] + 0.5 * Y %*% rr + 0.5 * Y %*% cc)
yb <- as.vector(1 + N[, k] - Y[, k] + 0.5 * (N - Y) %*% rr + 0.5 * (N - Y) %*% cc)
check("power, equal: own data + half of the same row and column", isTRUE(all.equal(pe$logpost[, k, ], beta_cells(grid, ya, yb))))
check("power with weights 0 = no borrowing",
      isTRUE(all.equal(engine_ebbb(grid, lay, "power")$posterior(Y, N)$logpost, nb$logpost)))

## EX/NEX: posterior probability of EX by the beta-binomial likelihood of the basket's own data
ex <- engine_ebbb(grid, lay, "EX/NEX", w_col = 1, p_ex = 0.5)$posterior(Y, N)
aE <- as.vector(1 + Y %*% cc); bE <- as.vector(1 + (N - Y) %*% cc)
mlE <- lbeta(aE + Y[, k], bE + N[, k] - Y[, k]) - lbeta(aE, bE); mlN <- lbeta(1 + Y[, k], 1 + N[, k] - Y[, k])
check("EX/NEX: P(EX) from the beta-binomial marginal likelihoods",
      max(abs(ex$pex[, k] - 1 / (1 + exp(mlN - mlE)))) < 1e-12)
mix <- log(ex$pex[, k] * exp(beta_cells(grid, aE + Y[, k], bE + N[, k] - Y[, k])) +
           (1 - ex$pex[, k]) * exp(beta_cells(grid, 1 + Y[, k], 1 + N[, k] - Y[, k])))
check("EX/NEX: mixture of the two updated components", max(abs(exp(ex$logpost[, k, ]) - exp(mix))) < 1e-12)
check("EX/NEX with p_ex = 1 = power, column", isTRUE(all.equal(
  engine_ebbb(grid, lay, "EX/NEX", w_col = 1, p_ex = 1)$posterior(Y, N)$logpost,
  engine_ebbb(grid, lay, "power", w_col = 1)$posterior(Y, N)$logpost)))
check("EX/NEX with p_ex = 0 = no borrowing", isTRUE(all.equal(
  engine_ebbb(grid, lay, "EX/NEX", w_col = 1, p_ex = 0)$posterior(Y, N)$logpost, nb$logpost)))

## Similarity: Jensen-Shannon divergence and weights
check("JSD of equal Betas is 0", max(jsd_beta(c(2, 5), c(9, 3), c(2, 5), c(9, 3))) < 1e-12)
check("JSD is symmetric", abs(jsd_beta(2, 9, 7, 4) - jsd_beta(7, 4, 2, 9)) < 1e-12)
check("JSD of very different Betas is close to 1", jsd_beta(1, 60, 60, 1) > 0.99, fmt(jsd_beta(1, 60, 60, 1)))
Ys <- matrix(c(3, 3, rep(c(0, 27), 5)), 1); Ns <- matrix(rep(30, 12), 1)   # baskets 1, 2 equal; others far away
si <- engine_ebbb(grid, grid_layout(4, 3), "similarity", epsilon = 2, tau = 0.3)$posterior(Ys, Ns)
w12 <- (1 - jsd_beta(4, 28, 4, 28))^2
check("similarity: an identical basket is pooled, dissimilar ones are ignored",
      isTRUE(all.equal(si$logpost[1, 1, ], as.vector(beta_cells(grid, 4 + w12 * 3, 28 + w12 * 27)))))
check("similarity: empty baskets are ignored",
      isTRUE(all.equal(engine_ebbb(grid, lay, "similarity")$posterior(matrix(c(2, rep(0, 11)), 1), matrix(c(10, rep(0, 11)), 1))$logpost[1, 1, ],
                       as.vector(beta_cells(grid, 3, 9)))))

## Layout: transposing the grid (rows <-> columns) transposes the power prior
tl <- grid_layout(3, 4); perm <- order(lay$col, lay$row)       # basket (j, l) of lay is (l, j) of tl
pc_t <- engine_ebbb(grid, tl, "power", w_row = 1)$posterior(Y[, perm], N[, perm])
pc <- engine_ebbb(grid, lay, "power", w_col = 1)$posterior(Y, N)
check("power: column borrowing on J x L = row borrowing on the transposed L x J grid",
      isTRUE(all.equal(pc_t$logpost, pc$logpost[, perm, , drop = FALSE])))

## Oracle: exact likelihood-ratio factors and the true rate as summary
q <- list(true = c(0.1, 0.3, 0.45, 0.05), eff = c(0.3, 0.3, 0.45, 0.3), fut = c(0.1, 0.1, 0.1, 0.05))
or <- engine_ebbb(grid, grid_layout(2, 2), "oracle", q = q)
po <- or$posterior(matrix(0, 1, 4), matrix(0, 1, 4))
og <- or$grid; fut <- futility_side(og, p1)
f_eff <- sapply(1:4, function(k) grid_factor(po$logpost[1, k, ], 2, 5, og))
f_fut <- sapply(1:4, function(k) exp(log_factor(restrict_below(matrix(po$logpost[1, k, ], 1), fut), 2, 5, fut)))
check("oracle: efficacy factor = (q/p0)^y ((1-q)/(1-p0))^(n-y)",
      max(abs(f_eff / ((q$eff / p0)^2 * ((1 - q$eff) / (1 - p0))^3) - 1)) < 1e-12)
check("oracle: futility factor = (q/p1)^y ((1-q)/(1-p1))^(n-y)",
      max(abs(f_fut / ((q$fut / p1)^2 * ((1 - q$fut) / (1 - p1))^3) - 1)) < 1e-12)
pm <- sapply(1:4, function(k) { pr <- exp(po$logpost[1, k, ]); sum(pr * og$p) / sum(pr) })
check("oracle: posterior summaries are the true rate", max(abs(pm - q$true)) < 1e-9)

## Validity: exact type I error of a null basket, every eBBB prior (as test_validity.R)
block <- 5; B <- 6; alpha <- 0.05
exact_typeI <- function(engine, K, rate_other) {
  g <- engine$grid; lf <- vector("list", B)
  for (k in seq_len(B)) {
    y1 <- 0:(block * (k - 1)); nk <- block * (k - 1)
    Yk <- cbind(y1, matrix(round(rate_other * nk), length(y1), K - 1)); Nk <- matrix(nk, length(y1), K)
    need <- matrix(FALSE, length(y1), K); need[, 1] <- TRUE
    w <- restrict(matrix(engine$posterior(Yk, Nk, need)$logpost[, 1, ], length(y1)), g)
    lf[[k]] <- sapply(0:block, function(yb) log_factor(w, rep(yb, nrow(w)), rep(block, nrow(w)), g))
  }
  paths <- as.matrix(expand.grid(rep(list(0:block), B))); prev <- cbind(0, t(apply(paths, 1, cumsum))[, -B])
  logE <- 0; rej <- FALSE
  for (k in seq_len(B)) { logE <- logE + lf[[k]][cbind(prev[, k] + 1, paths[, k] + 1)]; rej <- rej | logE >= log(1 / alpha) }
  sum(exp(rowSums(dbinom(paths, block, p0, log = TRUE))) * rej)
}
l22 <- grid_layout(2, 2)
eng <- list("no borrowing" = engine_ebbb(grid, l22, "no borrowing"),
            "power, equal" = engine_ebbb(grid, l22, "power", w_row = 0.5, w_col = 0.5),
            "power, column" = engine_ebbb(grid, l22, "power", w_col = 1),
            "EX/NEX" = engine_ebbb(grid, l22, "EX/NEX", w_col = 1, p_ex = 0.5),
            "similarity" = engine_ebbb(grid, l22, "similarity"),
            "oracle" = engine_ebbb(grid, l22, "oracle", q = list(true = rep(p0, 4), eff = rep(p1, 4), fut = rep(p0, 4))))
for (nm in names(eng)) for (r in c(0.15, 0.6)) {
  t1 <- exact_typeI(eng[[nm]], 4, r)
  check(sprintf("exact type I error <= alpha: %s, other baskets at %.2f", nm, r), t1 <= alpha, sprintf("%.4f", t1))
}

## The study settings: every model of setting.R builds, on the unbalanced 4 x 3 layout
st <- new.env(); sys.source(file.path(root, "setting.R"), envir = st)
me <- make_engine; environment(me) <- st; md <- make_design; environment(md) <- st
rec <- st$Recruitment$unbalanced; inn <- !is.na(rec$rate)
for (m in st$Models) {
  d <- md(rec, m, "fixed")
  e <- me(m, d$layout, st$Scenario[["strong heterogeneity"]][inn])
  pz <- e$posterior(matrix(0, 1, sum(inn)), matrix(0, 1, sum(inn)), seed = 1)
  ok <- identical(dim(pz$logpost), c(1L, sum(inn), length(e$grid$theta))) && d$futility == "none" && nrow(d$layout) == 11
  check(sprintf("setting.R: '%s' builds on the unbalanced 4 x 3 layout (11 baskets), fixed design without futility", m), ok)
}
check("setting.R: every scenario has 12 rates in (0, 1)", all(sapply(st$Scenario, function(s) length(s) == 12 && all(s > 0 & s < 1))))
check("setting.R: eBBB scenarios present", all(c("global null", "biomarker 1 active", "interaction", "all active") %in% names(st$Scenario)))

## Fixed design (efficacy or n_max): every basket ends rejected or at n_max
d <- md(st$Recruitment$balanced, "no borrowing", "fixed"); e <- me("no borrowing", d$layout, rep(p1, 12))
sim <- simulate_trials(rep(c(p0, p1), 6), d, e, n_sim = 4, seed = 3)
check("fixed design: every basket rejected or at n_max, none stopped for futility",
      all(sim$reject | sim$n == st$n_max) && !any(sim$stopped))
