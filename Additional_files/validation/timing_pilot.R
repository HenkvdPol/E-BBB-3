## =============================================================================
##  timing_pilot.R -- run time of the study, from a pilot with 20 trials per cell
##  Run from the repository folder:
##      Rscript Additional_files/validation/timing_pilot.R        (about 45 min on 2 cores)
##  Output (Additional_files/validation/results/):
##      timing_pilot_cells.csv   per pilot cell: set-up, engine and other seconds, number of fits
##      timing_pilot_calls.csv   per fit of the borrowing model: number of data states and seconds
##      timing_pilot_twoway.csv  two-way fit time against the number of data states (calibration)
##      timing_pilot_estimate.csv  predicted seconds per cell of the study at R = 1000
##      timing_pilot_session.txt
##  The pilot results (rejection rates etc.) are NOT interpreted: 20 trials per cell.
##
##  A. Pilot: every model x both recruitment designs x 3 scenarios (global null, all active, strong
##     heterogeneity) x both horizons, R = 20, the seeds of simulate.R. Every fit of the borrowing
##     model is timed with its number of data states.
##  B. Calibration of the two-way fit: t(S) for S = 20, 150, 300 and 600 data states on both layouts.
##  C. Prediction for R = 1000: per fit, the number of data states x 50 (an upper bound: identical
##     states are fitted once); the two-way fits by t(S) of B, the other models linearly; the rest of
##     the simulation linearly in R. Unmeasured scenarios: the largest measured value of the same
##     model, recruitment and horizon.
## =============================================================================
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "Additional_files/validation", "..", "..")) })
setwd(root)
source("setting.R")
source("additional_functions.R")
out_dir <- file.path(root, "Additional_files", "validation", "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
R_pilot <- 20
R_study <- 1000
n_cores <- min(2, parallel::detectCores())

## The cells of the study, numbered as in simulate.R (the fit seeds depend on the number)
sim <- expand.grid(Recruitment = names(Recruitment), Scenario = names(Scenario), Model = Models, stringsAsFactors = FALSE)
sim$id <- seq_len(nrow(sim))
pilot_scen <- c("global null", "all active", "strong heterogeneity")
cells <- rbind(cbind(sim[sim$Scenario %in% pilot_scen, ], horizon = "fixed"),
               cbind(sim[sim$Scenario %in% pilot_scen, ], horizon = "open"))
cells <- cells[order(cells$Model != "two-way EX/NEX", cells$Recruitment != "unbalanced"), ]   # slowest first

################################################################################
# A. Pilot
################################################################################
run_cell <- function(i){
  cl <- cells[i, ]
  t0 <- proc.time()[["elapsed"]]
  p_true <- Scenario[[cl$Scenario]][!is.na(Recruitment[[cl$Recruitment]]$rate)]
  design <- make_design(Recruitment[[cl$Recruitment]], cl$Model, cl$horizon)
  engine <- make_engine(cl$Model, design$layout, p_true)
  setup <- proc.time()[["elapsed"]] - t0

  # Time every fit of the borrowing model
  calls <- NULL
  fit <- engine$posterior
  engine$posterior <- function(Y, N, need = NULL, seed = NULL){
    s0 <- proc.time()[["elapsed"]]
    out <- fit(Y, N, need, seed)
    calls <<- rbind(calls, data.frame(S = nrow(Y), n_max = max(N), seconds = proc.time()[["elapsed"]] - s0))
    return(out)
  }
  seed_outcome <- startSeed0 + 1000 * match(cl$Scenario, names(Scenario)) + match(cl$Recruitment, names(Recruitment))
  seed_fit <- startSeed0 + 10^6 + 10^4 * cl$id
  t1 <- proc.time()[["elapsed"]]
  d <- simulate_trials(p_true, design, engine, R_pilot, seed = seed_outcome, fit_seed = seed_fit)
  total <- proc.time()[["elapsed"]] - t1
  calls$call <- seq_len(nrow(calls))
  calls <- cbind(cl[rep(1, nrow(calls)), c("id", "Recruitment", "Scenario", "Model", "horizon")], calls)
  return(list(cell = data.frame(cl[, c("id", "Recruitment", "Scenario", "Model", "horizon")], R = R_pilot,
                                setup = setup, engine = sum(calls$seconds), other = total - sum(calls$seconds),
                                fits = nrow(calls), states = sum(calls$S), duration = max(d$time), n_max_seen = max(d$n)),
              calls = calls))
}
cat("A. pilot:", nrow(cells), "cells, R =", R_pilot, "\n")
tA <- proc.time()[["elapsed"]]
res <- parallel::mclapply(seq_len(nrow(cells)), run_cell, mc.cores = n_cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, TRUE, "try-error")
if (any(bad)) stop("pilot cell failed: ", conditionMessage(attr(res[[which(bad)[1]]], "condition")))
cell_res <- do.call(rbind, lapply(res, `[[`, "cell"))
call_res <- do.call(rbind, lapply(res, `[[`, "calls"))
write.csv(cell_res, file.path(out_dir, "timing_pilot_cells.csv"), row.names = FALSE)
write.csv(call_res, file.path(out_dir, "timing_pilot_calls.csv"), row.names = FALSE)
cat(sprintf("   pilot wall time %.1f min\n", (proc.time()[["elapsed"]] - tA) / 60))

################################################################################
# B. Calibration of the two-way fit time (one process, nothing else running)
################################################################################
cat("B. two-way calibration\n")
cal <- NULL
for (rec in names(Recruitment)) {
  design <- make_design(Recruitment[[rec]], "two-way EX/NEX", "fixed")
  K <- nrow(design$layout)
  engine <- make_engine("two-way EX/NEX", design$layout, rep(p0, K))
  engine$posterior(matrix(0, 1, K), matrix(0, 1, K), seed = 1)                 # tables, as after the first fit
  for (S in c(20, 150, 300, 600)) {
    set.seed(S + K)
    N <- matrix(sample(seq(0, 30, 5), S * K, TRUE), S)                          # states as in the fixed design
    Y <- matrix(rbinom(S * K, N, sample(c(p0, p1), S * K, TRUE)), S)
    s0 <- proc.time()[["elapsed"]]
    engine$posterior(Y, N, seed = 2)
    cal <- rbind(cal, data.frame(Recruitment = rec, baskets = K, S = S, seconds = proc.time()[["elapsed"]] - s0))
    cat(sprintf("   %s, S = %d: %.1f s\n", rec, S, tail(cal$seconds, 1)))
  }
}
write.csv(cal, file.path(out_dir, "timing_pilot_twoway.csv"), row.names = FALSE)
coef_tw <- lapply(split(cal, cal$Recruitment), function(x) coef(lm(seconds ~ S, x)))

################################################################################
# C. Prediction for R = 1000
################################################################################
scale <- R_study / R_pilot
call_res$S_study <- ifelse(call_res$call == 1, call_res$S, call_res$S * scale)     # the empty-platform fit stays 1 state
call_res$sec_study <- ifelse(call_res$Model == "two-way EX/NEX",
  mapply(function(rec, S) sum(coef_tw[[rec]] * c(1, S)), call_res$Recruitment, call_res$S_study),
  call_res$seconds * call_res$S_study / call_res$S)
eng_study <- aggregate(sec_study ~ id + horizon, call_res, sum)
est <- merge(cell_res, eng_study, by = c("id", "horizon"))
est$seconds_study <- est$setup + est$sec_study + est$other * scale
## All 154 cells per horizon: unmeasured scenarios get the largest measured value of the same model,
## recruitment and horizon
worst <- aggregate(seconds_study ~ Model + Recruitment + horizon, est, max)
full <- merge(rbind(cbind(sim, horizon = "fixed"), cbind(sim, horizon = "open")), worst, by = c("Model", "Recruitment", "horizon"))
write.csv(est, file.path(out_dir, "timing_pilot_estimate.csv"), row.names = FALSE)

cat("\nPredicted CPU hours at R = 1000, all 154 cells (upper estimate):\n")
tab <- aggregate(seconds_study ~ Model + horizon, full, function(x) sum(x) / 3600)
print(reshape(tab, idvar = "Model", timevar = "horizon", direction = "wide"), digits = 3, row.names = FALSE)
tot <- aggregate(seconds_study ~ horizon, full, function(x) c(cpu_hours = sum(x) / 3600, longest_cell_hours = max(x) / 3600))
print(tot, digits = 3)
cat("\nWall time with 8 cores (as sim.slurm): about max(CPU hours / 8, longest cell)\n")
for (h in unique(full$horizon)) {
  x <- full$seconds_study[full$horizon == h] / 3600
  cat(sprintf("   %s: %.1f h (CPU %.1f h, longest cell %.1f h)\n", h, max(sum(x) / 8, max(x)), sum(x), max(x)))
}
cat("\nMonte Carlo standard error of a rejection rate:\n")
se <- outer(c(500, 1000, 2000), c(0.05, 0.2, 0.5), function(R, p) sqrt(p * (1 - p) / R))
dimnames(se) <- list(R = c(500, 1000, 2000), rate = c(0.05, 0.2, 0.5)); print(round(se, 4))
writeLines(c(capture.output(sessionInfo()), paste("cores:", n_cores)), file.path(out_dir, "timing_pilot_session.txt"))
