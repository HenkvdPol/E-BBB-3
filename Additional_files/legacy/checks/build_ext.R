## Run from anywhere:  Rscript checks/build_ext.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings); t0 <- Sys.time(); before <- length(list.files(settings$sim$bank_dir))
for (i in seq_len(nrow(runs))) { s <- complete_settings(settings, runs[i, ]); invisible(get_banks(s)) }
cat("new design files:", length(list.files(settings$sim$bank_dir)) - before, " seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
s <- complete_settings(settings, runs[runs$name == "simon_p30-50_K5_n46", ]); b <- get_banks(s)
ex <- ext_banks_for(s$designs[["e-based + mEXNEX (flexible)"]], s, b)[[1]]
cat("extension bank n_ext =", dim(ex$bet_idx)[2], " theta1 values:", length(ex$theta1), " cp matrix:", dim(ex$cp), "\n")
## conditional power from the bank = 1 - cost of the library's P-max dynamic program at t = 0
d <- design_pmax(s$n_ext[1], s$p0, 0.5, s$alpha, analysis_times = seq_len(s$n_ext[1]))
j <- which(abs(ex$theta1 - 0.5) < 1e-9)
cat("max |cp(bank) - (1 - cost DP)| over all e-values at theta1 = 0.5:", signif(max(abs(ex$cp[, j] - (1 - d$cost[, 1]))), 3), "\n")
i <- sapply(c(0.5, 1, 2, 5, 10), function(m) findInterval(m, ex$m_grid))
print(data.frame(e_value = c(0.5, 1, 2, 5, 10), cond_power_theta1_0.5 = round(ex$cp[i, j], 3), cond_power_theta1_0.35 = round(ex$cp[i, 21 - 20 + 0], 3)))
