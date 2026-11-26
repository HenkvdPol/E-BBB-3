## Run from anywhere:  Rscript checks/grow_check.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings); worst <- 0; n_checked <- 0
for (rn in c("grid_p10-30_K5_n30", "simon_p70-90_K5_n27", "simon_p05-25_K5_n17")) {
  s <- complete_settings(settings, runs[runs$name == rn, ]); b <- get_banks(s)
  for (sn in names(s$schedules)) { bk <- banks_for(s$designs[["e-GROW"]], s$schedules[[sn]], s$setups$basket, s, b)[[1]]
    i <- as.integer(bk$bet_idx[-c(1, length(bk$m_grid)), , 1]); bets <- bk$b_val0[i[i > 0] + 1]
    kelly <- (bk$theta1 - s$p0) / (1 - s$p0); worst <- max(worst, abs(bets - kelly)); n_checked <- n_checked + length(bets) } }
cat("e-GROW: largest |bet - (theta1 - p0)/(1 - p0)| over", n_checked, "state-patient entries:", signif(worst, 3), "\n")
