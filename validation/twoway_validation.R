## =============================================================================
##  twoway_validation.R -- the production two-way engine against reference calculations
##  Run:    Rscript validation/twoway_validation.R        (about 20-30 min on one core)
##  Output: validation/results/twoway_validation.csv      (every method, state and basket)
##          validation/results/twoway_validation_summary.csv
##  No trial simulation: posterior calculations for fixed data states only.
##
##  A. 3 x 1 layout, full hyperpriors: deterministic quadrature (twoway_exact) and the equivalent
##     one-way EXNEX model (one-way grid engine) as references.
##  B. 3 x 3 layout, sigma_a = sigma_b = sigma_e = 0.45 fixed: quadrature as reference.
##  C. 3 x 3 layout, full hyperpriors: JAGS (4 x 25000 draws) and a long run of the sampler
##     (4 chains x 20000) as references; the production sampler with 5 seeds.
## =============================================================================
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "validation", "..")) })
for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
out_dir <- file.path(root, "validation", "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p0 <- 0.1; pbar <- 0.2; grid <- theta_grid(p0); prior <- twoway_prior(p0, pbar); tab <- twoway_tables(grid, prior, 30)
prod <- list(burn = 500, keep = 1000, chains = 2, chunk = 150)          # production settings
long <- list(burn = 2000, keep = 20000, chains = 4, chunk = 150)
l33 <- grid_layout(3, 3); l31 <- grid_layout(3, 1)
v <- function(...) c(...)
states33 <- list(                                     # baskets in the order (1,1) (1,2) (1,3) (2,1) ... (3,3)
  "homogeneous null"       = list(y = rep(3, 9), n = rep(30, 9)),
  "one active basket"      = list(y = v(10, rep(3, 8)), n = rep(30, 9)),
  "active row"             = list(y = v(9, 9, 9, rep(3, 6)), n = rep(30, 9)),
  "active column"          = list(y = v(9, 3, 3, 9, 3, 3, 9, 3, 3), n = rep(30, 9)),
  "interaction (diagonal)" = list(y = v(9, 3, 3, 3, 9, 3, 3, 3, 9), n = rep(30, 9)),
  "unbalanced, empty (3,2)" = list(y = v(3, 6, 1, 7, 1, 4, 9, 0, 2), n = v(30, 20, 10, 25, 5, 15, 30, 0, 10)),
  "missing (3,3), empty (2,2)" = list(y = v(9, 9, 9, 3, 0, 3, 3, 3), n = v(30, 30, 30, 30, 0, 30, 30, 30),
                                      layout = grid_layout(3, 3, missing = list(c(3, 3)))),
  "strong heterogeneity"   = list(y = v(1, 3, 14, 3, 9, 3, 18, 3, 1), n = rep(30, 9)),
  "early (5 per basket)"   = list(y = v(0, 1, 2, 1, 0, 0, 2, 1, 3), n = rep(5, 9)))
states31 <- list("heterogeneous" = list(y = v(3, 9, 14), n = rep(30, 3)), "homogeneous null" = list(y = rep(3, 3), n = rep(30, 3)),
                 "unbalanced" = list(y = v(0, 2, 5), n = v(5, 10, 30)), "empty basket" = list(y = v(1, 0, 2), n = v(5, 0, 5)))

rows <- list(); timing <- list()
add <- function(section, state, method, s, secs) {
  s$section <- section; s$state <- state; s$method <- method; s$basket <- seq_len(nrow(s))
  rows[[length(rows) + 1]] <<- s; timing[[length(timing) + 1]] <<- data.frame(section, state, method, seconds = secs)
}
timed <- function(expr) { t0 <- Sys.time(); x <- expr; list(x = x, s = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
grid_summ <- function(po) summ_grid(matrix(po$logpost[1, , ], dim(po$logpost)[2]), po$pex[1, ], grid)
mcmc_run <- function(st, lay, pr, tb, set, seed) {
  e <- engine_twoway(grid, lay, pbar, 30, prior = pr, mcmc = set, tab = tb); e$posterior(matrix(st$y, 1), matrix(st$n, 1), seed = seed) }

## A. 3 x 1 layout, full hyperpriors
ow <- engine_oneway(grid, pbar, 30, hyper = oneway_equivalent_hyper(prior))
for (nm in names(states31)) { st <- states31[[nm]]; cat("A:", nm, "\n")
  r <- timed(twoway_exact(st$y, st$n, l31, tab, prior, grid, tail_step = 1)); add("A", nm, "exact", summ_grid(r$x$logpost, r$x$pex, grid), r$s)
  r <- timed(ow$posterior(matrix(st$y, 1), matrix(st$n, 1))); add("A", nm, "one-way equivalent", grid_summ(r$x), r$s)
  for (sd in 1:5) { r <- timed(mcmc_run(st, l31, prior, tab, prod, sd)); add("A", nm, sprintf("mcmc seed %d", sd), grid_summ(r$x), r$s) }
  r <- timed(twoway_jags(st$y, st$n, l31, prior, grid, n_iter = 25000, n_burn = 5000)); add("A", nm, "jags", summ_draws(r$x$theta, r$x$d, grid), r$s)
}
## B. 3 x 3, fixed sigmas
pf <- fix_sigmas(prior, 0.45, 0.45, 0.45); tf <- twoway_tables(grid, pf, 30)
for (nm in names(states33)) { st <- states33[[nm]]; lay <- if (is.null(st$layout)) l33 else st$layout; cat("B:", nm, "\n")
  r <- timed(twoway_exact(st$y, st$n, lay, tf, pf, grid, c_half = 3.5, a_max = 2.5, tail_step = 3))
  add("B", nm, "exact", summ_grid(r$x$logpost, r$x$pex, grid), r$s)
  for (sd in 1:5) { r <- timed(mcmc_run(st, lay, pf, tf, prod, sd)); add("B", nm, sprintf("mcmc seed %d", sd), grid_summ(r$x), r$s) }
  r <- timed(twoway_jags(st$y, st$n, lay, pf, grid, n_iter = 25000, n_burn = 5000)); add("B", nm, "jags", summ_draws(r$x$theta, r$x$d, grid), r$s)
}
## C. 3 x 3, full hyperpriors
for (nm in names(states33)) { st <- states33[[nm]]; lay <- if (is.null(st$layout)) l33 else st$layout; cat("C:", nm, "\n")
  r <- timed(twoway_jags(st$y, st$n, lay, prior, grid, n_iter = 25000, n_burn = 5000)); add("C", nm, "jags", summ_draws(r$x$theta, r$x$d, grid), r$s)
  r <- timed(mcmc_run(st, lay, prior, tab, long, 99)); add("C", nm, "mcmc long", grid_summ(r$x), r$s)
  for (sd in 1:5) { r <- timed(mcmc_run(st, lay, prior, tab, prod, sd)); add("C", nm, sprintf("mcmc seed %d", sd), grid_summ(r$x), r$s) }
}
res <- do.call(rbind, rows); tim <- do.call(rbind, timing)
write.csv(res, file.path(out_dir, "twoway_validation.csv"), row.names = FALSE)

## Largest differences against the reference of each section
ref_of <- c(A = "exact", B = "exact", C = "jags")
cmp <- do.call(rbind, lapply(split(res, list(res$section, res$state), drop = TRUE), function(d) {
  ref <- d[d$method == ref_of[[d$section[1]]], ]; ref <- ref[order(ref$basket), ]
  do.call(rbind, lapply(setdiff(unique(d$method), ref$method[1]), function(m) { x <- d[d$method == m, ]; x <- x[order(x$basket), ]
    data.frame(section = d$section[1], state = d$state[1], method = sub(" seed [0-9]+", "", m),
               rmean = max(abs(x$rmean - ref$rmean)), rsd = max(abs(x$rsd - ref$rsd)), p_above = max(abs(x$p_above - ref$p_above)),
               pex = max(abs(x$pex - ref$pex)), logf0 = max(abs(log(x$f0 / ref$f0))), logf5 = max(abs(log(x$f5 / ref$f5)))) }))
}))
summ <- aggregate(cbind(rmean, rsd, p_above, pex, logf0, logf5) ~ section + method, cmp, max)
tim$method <- sub(" seed [0-9]+", "", tim$method)
summ <- merge(summ, aggregate(seconds ~ section + method, tim, mean), sort = FALSE)
## Monte Carlo spread of the production sampler over the 5 seeds (largest sd over states and baskets)
mc <- res[grepl("^mcmc seed", res$method), ]; mc$logf5 <- log(mc$f5)
spread <- aggregate(cbind(rmean, pex, logf5) ~ section + state + basket, mc, sd)
spread <- aggregate(cbind(rmean, pex, logf5) ~ section, spread, max)
names(spread) <- c("section", "sd_rmean", "sd_pex", "sd_logf5")
write.csv(summ, file.path(out_dir, "twoway_validation_summary.csv"), row.names = FALSE)
write.csv(spread, file.path(out_dir, "twoway_validation_mcspread.csv"), row.names = FALSE)
print(format(summ, digits = 3), row.names = FALSE); print(format(spread, digits = 3), row.names = FALSE)
sha <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA)
writeLines(c(paste("git", sha), capture.output(sessionInfo())), file.path(out_dir, "twoway_validation_session.txt"))
