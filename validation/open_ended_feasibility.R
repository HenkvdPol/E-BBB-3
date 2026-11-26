## =============================================================================
##  open_ended_feasibility.R -- exact operating characteristics of open-ended baskets
##  Run:    Rscript validation/open_ended_feasibility.R          (a few minutes; no simulation)
##  Output: validation/results/open_ended_feasibility.csv, open_ended_reject_by_n.csv, _session.txt
##
##  One basket with p0 = 0.10, p1 = 0.30, alpha = 0.05, beta = 0.20, followed up to 3000 patients.
##  Exact (dynamic programme, R/open_ended.R) for starting priors that do not change with other data:
##    none                       NEX prior, no borrowing
##    one-way / two-way EXNEX    borrowed predictive given 8 other baskets whose data are complete
##  Designs:
##    fixed n_max 30             the current study (stop when 1/alpha is out of reach)
##    open, no futility          no maximum sample size, efficacy only
##    open, e-value futility     no maximum sample size, anytime-valid futility at 1/beta
## =============================================================================
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "validation", "..")) })
for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
out_dir <- file.path(root, "validation", "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p0 <- 0.10; p1 <- 0.30; pbar <- (p0 + p1) / 2; alpha <- 0.05; beta <- 0.20; horizon <- 3000
grid <- theta_grid(p0); fut <- futility_side(grid, p1)
p_grid <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40)

## ---- starting priors: basket 1 has no data yet; the other baskets' data are complete -------------
start <- list()
add_start <- function(name, lp) start[[name]] <<- list(w0 = as.vector(restrict(lp, grid)), w0f = as.vector(restrict_below(lp, fut)))
add_start("none", matrix(log(nex_prior(grid, pbar)), 1))
ow <- engine_oneway(grid, pbar, 30)
others <- list("8 null baskets (3/30)" = rep(3, 8), "8 active baskets (9/30)" = rep(9, 8), "4 null + 4 active" = c(3, 3, 3, 3, 9, 9, 9, 9))
for (nm in names(others))
  add_start(paste("one-way EXNEX,", nm), matrix(ow$posterior(matrix(c(0, others[[nm]]), 1), matrix(c(0, rep(30, 8)), 1))$logpost[1, 1, ], 1))
tw <- engine_twoway(grid, grid_layout(3, 3), pbar, 30, mcmc = list(burn = 2000, keep = 20000, chains = 4, chunk = 150))
cells <- function(m) as.vector(t(m))
others2 <- list(                                       # basket (1,1) is the one followed; y of the others out of 30
  "own row active (9/30), rest null"    = cells(matrix(c(0, 9, 9,  3, 3, 3,  3, 3, 3), 3, byrow = TRUE)),
  "own row and column null, rest active" = cells(matrix(c(0, 3, 3,  3, 9, 9,  3, 9, 9), 3, byrow = TRUE)),
  "all other baskets null"              = cells(matrix(c(0, 3, 3,  3, 3, 3,  3, 3, 3), 3, byrow = TRUE)))
for (nm in names(others2))
  add_start(paste("two-way EXNEX,", nm),
            matrix(tw$posterior(matrix(others2[[nm]], 1), matrix(c(0, rep(30, 8)), 1), seed = 7)$logpost[1, 1, ], 1))

## ---- designs ------------------------------------------------------------------------------------
designs <- list(
  "fixed n_max 30, blocks of 5"          = list(block = 5, n_max = 30, futility = "unreachable"),
  "open, no futility, blocks of 5"       = list(block = 5, n_max = Inf, futility = "none"),
  "open, e-value futility, blocks of 5"  = list(block = 5, n_max = Inf, futility = "evalue"),
  "open, e-value futility, every patient" = list(block = 1, n_max = Inf, futility = "evalue"))

rows <- list(); by_n <- list(); t0 <- Sys.time()
for (sn in names(start)) for (dn in names(designs)) for (p in p_grid) {
  d <- designs[[dn]]
  r <- exact_basket(p, start[[sn]]$w0, grid, alpha, block = d$block, n_max = d$n_max, horizon = horizon,
                    futility = d$futility, w0f = start[[sn]]$w0f, fut = fut, beta = beta)
  rows[[length(rows) + 1]] <- cbind(prior = sn, design = dn, r$summary)
  by_n[[length(by_n) + 1]] <- data.frame(prior = sn, design = dn, p = p, n = as.numeric(names(r$reject_by)), reject_by_n = r$reject_by)
}
res <- do.call(rbind, rows); byn <- do.call(rbind, by_n)
write.csv(res, file.path(out_dir, "open_ended_feasibility.csv"), row.names = FALSE)
write.csv(byn, file.path(out_dir, "open_ended_reject_by_n.csv"), row.names = FALSE)
cat(sprintf("computed in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## ---- checks of the guarantees (exact up to the horizon) -----------------------------------------
null_rows <- res[res$p <= p0, ]; alt_rows <- res[res$p >= p1 & res$design != "fixed n_max 30, blocks of 5", ]
efut <- grepl("e-value futility", res$design)
cat(sprintf("largest P(reject) for p <= p0 over all priors and designs: %.4f (alpha = %.2f)\n", max(null_rows$reject), alpha))
cat(sprintf("largest P(futility stop) for p >= p1, e-value futility: %.4f (beta = %.2f)\n", max(res$futility[efut & res$p >= p1]), beta))
cat(sprintf("largest P(still open at n = %d), e-value futility: %.2e\n", horizon, max(res$open[efut])))
fmt_tab <- res[res$p %in% c(0.10, 0.20, 0.30), c("prior", "design", "p", "reject", "futility", "open", "EN", "N_median", "N_q90")]
print(format(fmt_tab, digits = 3), row.names = FALSE)
sha <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA)
writeLines(c(paste("git", sha), capture.output(sessionInfo())), file.path(out_dir, "open_ended_feasibility_session.txt"))
