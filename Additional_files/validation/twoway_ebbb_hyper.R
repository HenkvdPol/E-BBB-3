## =============================================================================
##  twoway_ebbb_hyper.R -- the two-way engine with the hyperparameters of eBBB, against JAGS
##  Run from the repository folder:
##      Rscript Additional_files/validation/twoway_ebbb_hyper.R      (about 10 min on one core)
##  Output: Additional_files/validation/results/twoway_ebbb_hyper.csv (every method, state, basket)
##          Additional_files/validation/results/twoway_ebbb_hyper_summary.csv, _session.txt
##  No trial simulation: posterior calculations for fixed data states only.
##
##  4 x 3 layout (as in eBBB and setting.R), prior of setting.R / make_engine():
##  mu ~ N(logit p0, 2^2), sigma_a, sigma_b ~ HN(1), sigma_eps ~ HN(0.5), NEX N(logit p0, 2^2),
##  p_ex = 0.5. References: JAGS (4 x 25000 draws) and a long run of the sampler
##  (4 chains x 20000). The production sampler (setting.R: mcmc) with 5 seeds.
## =============================================================================
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "Additional_files/validation", "..", "..")) })
source(file.path(root, "setting.R"))
source(file.path(root, "additional_functions.R"))
for (f in list.files(file.path(root, "Additional_files", "validation"), pattern = "_functions\\.R$", full.names = TRUE)) source(f)
out_dir <- file.path(root, "Additional_files", "validation", "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- theta_grid(p0); pbar <- (p0 + p1) / 2
lay <- grid_layout(J, L); lay_miss <- grid_layout(J, L, missing = list(c(3, 3)))
eng <- make_engine("two-way EX/NEX", lay, rep(p0, nrow(lay)))
prior <- environment(eng$posterior)$prior
tab <- twoway_tables(grid, prior, 30)
long <- list(burn = 2000, keep = 20000, chains = 4, chunk = 150)
v <- function(...) c(...)
## baskets in the order (1,1) (1,2) (1,3) (2,1) ... (4,3)
states <- list(
  "global null"          = list(y = rep(2, 12), n = rep(20, 12)),
  "biomarker 1 active"   = list(y = v(6, 2, 2, 6, 2, 2, 6, 2, 2, 6, 2, 2), n = rep(20, 12)),
  "interaction"          = list(y = v(6, 2, 2, 6, 2, 2, 6, 2, 2, 2, 2, 2), n = rep(20, 12)),
  "all active"           = list(y = rep(6, 12), n = rep(20, 12)),
  "unbalanced, missing (3,3)" = list(y = v(3, 6, 1, 7, 1, 0, 2, 1, 0, 1, 0), n = v(30, 20, 10, 25, 5, 0, 10, 10, 5, 5, 5),
                                     layout = lay_miss),
  "early (5 per basket)" = list(y = v(0, 1, 2, 1, 0, 0, 2, 1, 3, 0, 1, 0), n = rep(5, 12)))

rows <- list()
add <- function(state, method, s, secs) {
  s$state <- state; s$method <- method; s$basket <- seq_len(nrow(s)); s$seconds <- secs
  rows[[length(rows) + 1]] <<- s
}
timed <- function(expr) { t0 <- Sys.time(); x <- expr; list(x = x, s = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
grid_summ <- function(po) summ_grid(matrix(po$logpost[1, , ], dim(po$logpost)[2]), po$pex[1, ], grid)
run <- function(st, set, seed) {
  l <- if (is.null(st$layout)) lay else st$layout
  e <- engine_twoway(grid, l, pbar, 30, prior = prior, mcmc = set, tab = tab)
  e$posterior(matrix(st$y, 1), matrix(st$n, 1), seed = seed)
}
for (nm in names(states)) { st <- states[[nm]]; cat(nm, "\n")
  l <- if (is.null(st$layout)) lay else st$layout
  r <- timed(twoway_jags(st$y, st$n, l, prior, grid, n_iter = 25000, n_burn = 5000)); add(nm, "jags", summ_draws(r$x$theta, r$x$d, grid), r$s)
  r <- timed(run(st, long, 99)); add(nm, "long mcmc", grid_summ(r$x), r$s)
  for (sd in 1:5) { r <- timed(run(st, mcmc, sd)); add(nm, sprintf("mcmc seed %d", sd), grid_summ(r$x), r$s) }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "twoway_ebbb_hyper.csv"), row.names = FALSE)

## Largest differences per state: production sampler (5 seeds) vs JAGS and vs the long run
ref <- function(m) res[res$method == m, ]
cmp <- function(a, b, what) {
  k <- c("p_above", "rmean", "pex")
  d <- merge(a, b, by = c("state", "basket"))
  out <- aggregate(cbind(p_above = abs(d$p_above.x - d$p_above.y), rmean = abs(d$rmean.x - d$rmean.y),
                         pex = abs(d$pex.x - d$pex.y), f5_rel = abs(d$f5.x / d$f5.y - 1)) ~ d$state, FUN = max)
  names(out)[1] <- "state"; out$comparison <- what; out
}
prod <- res[grepl("^mcmc seed", res$method), ]
summ <- rbind(cmp(prod, ref("jags"), "production vs JAGS"), cmp(prod, ref("long mcmc"), "production vs long run"),
              cmp(ref("long mcmc"), ref("jags"), "long run vs JAGS"))
summ$seconds_production <- NA
for (i in seq_len(nrow(summ))) summ$seconds_production[i] <- mean(prod$seconds[prod$state == summ$state[i]])
write.csv(summ, file.path(out_dir, "twoway_ebbb_hyper_summary.csv"), row.names = FALSE)
print(summ, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "twoway_ebbb_hyper_session.txt"))
