## Run from anywhere:  Rscript checks/regression_mixture.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
old <- new.env(); sys.source(file.path(check_dir, "old_versions", "additional_functions_v3.R"), envir = old)
runs <- complete_runs(settings)
handmade <- function(S, K, N) array(((outer(outer(1:S * 0.7548777, 1:K * 0.5698403, "+"), 1:N * 0.3819660, "+")) %% 1), c(S, K, N))
checked <- diffs <- 0
for (rn in c("grid_p10-30_K5_n30", "simon_p70-90_K5_n27")) {
  s <- complete_settings(settings, runs[runs$name == rn, ]); eng <- build_engine(s)
  U <- handmade(40, s$K, max(s$n_max + s$n_ext))
  for (su in names(s$setups)) for (dn in c("e-mixture (framework)", "e-mixture (no borrowing)")) {
    if (!(dn %in% s$setups[[su]]$designs)) next; d <- s$designs[[dn]]
    for (sn in names(s$schedules)) for (scn in names(s$scenarios)) {
      C <- make_data(U, s$scenarios[[scn]]); en <- s$setups[[su]]$entry
      a <- run_trials_mixture(d, C, s$schedules[[sn]], en, s, eng, new_cache())
      b <- old$run_trials_mixture(d, C, s$schedules[[sn]], en, s, eng, old$new_cache())
      f <- c(status = identical(a$status, b$status), n = identical(as.integer(a$n), as.integer(b$n)), e = isTRUE(all.equal(a$e_value, b$e_value)), time = isTRUE(all.equal(a$time, b$time))); same <- all(f); if (!same && diffs < 3) print(f)
      checked <- checked + 1; if (!same) { diffs <- diffs + 1; cat("DIFFERENT:", rn, su, dn, sn, scn, "\n") } }
  }
}
cat("framework e-value, new vs previous engine:", checked, "combinations,", diffs, "different\n")
## predictive probabilities of an extension block sum to 1
s <- settings$for_hypothesis(settings, 0.1, 0.3); set.seed(3); worst <- 0
for (i in 1:50) { a <- runif(1, 0.3, 20); b <- runif(1, 0.3, 40); ne <- sample(3:25, 1)
  tot <- sum(vapply(0:ne, function(y) dbinom(y, ne, 0.1) * mixture_factor(a, b, y, ne, 0.1), 0)); worst <- max(worst, abs(tot - 1)) }
cat("prior-predictive probabilities of an extension block sum to 1: largest error", signif(worst, 3), "\n")
