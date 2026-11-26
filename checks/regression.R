## Run from anywhere:  Rscript checks/regression.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
## Old vs new engine on hand-made data (40 fixed trials): every existing design must give identical decisions
source(file.path(check_dir, "load_settings.R"))
old <- new.env(); sys.source(file.path(check_dir, "old_versions", "additional_functions_old.R"), envir = old)
runs <- complete_runs(settings)
handmade <- function(S, K, N) array(((outer(outer(1:S * 0.7548777, 1:K * 0.5698403, "+"), 1:N * 0.3819660, "+")) %% 1), c(S, K, N))
diffs <- 0; checked <- 0
for (rn in c("grid_p10-30_K5_n30", "simon_p70-90_K5_n27", "grid_p10-30_K3_n20")) {
  s <- complete_settings(settings, runs[runs$name == rn, ]); so <- old$complete_settings(settings, runs[runs$name == rn, ])
  eng <- build_engine(s); eng_o <- old$build_engine(so); banks <- get_banks(s)
  U <- handmade(40, s$K, max(s$n_max + s$n_ext))
  for (su in names(s$setups)) for (dn in s$setups[[su]]$designs) {
    d <- s$designs[[dn]]; if (isTRUE(d$accelerate) || isTRUE(d$extend) || d$type == "mixture") next
    for (sn in names(s$schedules)) for (scn in names(s$scenarios)) {
      C <- make_data(U, s$scenarios[[scn]]); block <- s$schedules[[sn]]; en <- s$setups[[su]]$entry
      if (d$type == "bayes") {
        a <- run_trials_bayes(d, C, block, en, s, eng, new_cache(), 0.97)
        b <- old$run_trials_bayes(d, C[, , 1:(max(s$n_max) + 1)], block, en, so, eng_o, old$new_cache(), 0.97)
      } else {
        bk <- banks_for(d, block, s$setups[[su]], s, banks)
        a <- run_trials_evalue(d, C, block, en, s, eng, new_cache(), bk)
        b <- old$run_trials_evalue(d, C[, , 1:(max(s$n_max) + 1)], block, en, so, eng_o, old$new_cache(), bk)
      }
      same <- identical(a$status, b$status) && identical(as.integer(a$n), as.integer(b$n)) &&
              (is.null(b$e_value) || isTRUE(all.equal(a$e_value, b$e_value)))
      checked <- checked + 1; if (!same) { diffs <- diffs + 1; cat("DIFFERENT:", rn, su, dn, sn, scn, "\n") }
    }
  }
}
cat("design x set-up x schedule x scenario combinations compared:", checked, " different:", diffs, "\n")
