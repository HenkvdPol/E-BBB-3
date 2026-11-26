## =============================================================================
##  Setting.R  --  all settings AND the complete run
##
##  Run:    Rscript Setting.R              (e.g. as the only line in a slurm job)
##  Needs:  additional_functions.R and design_optimal_evalues.R in this folder.
##  Output: Result/<run>/   one folder per run (tables and figures)
##          Result/         all runs together (table and overview figures)
##
##  Evaluates the test of framework.tex: per basket an anytime-valid e-value whose
##  alternative is a prior from a Bayesian borrowing model (EXNEX) fitted to all
##  baskets, frozen at every analysis and restricted to (p0, 1]. Type I error <= alpha
##  per basket at every sample size, under optional stopping and continuation.
##  Primary design: "e-mixture (framework)". Alternatives within the framework:
##  design-optimal e-values (Baas et al.) that bet on a borrowed point alternative.
##  Comparators outside it: Bayesian posterior-probability designs, Simon, binomial test.
## =============================================================================

start_time <- Sys.time()
this_dir <- local({                       # folder of this file, wherever the job starts
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})

settings <- list()

## ---- 1. Error rates (per cohort) -----------------------------------------------------------
settings$alpha <- 0.05          # one-sided; e-value threshold 1/alpha = 20
settings$beta  <- 0.20          # target type II error (e-based design, Simon design)

## ---- 2. The runs: hypotheses H0: p <= p0 vs H1: p >= p1, number of cohorts K, cohort size n_max
##  series "grid" : p0 = 0.10, p1 = 0.30; K = 3, 5, 10 x n_max = 20, 30, 40
##  series "simon": the 8 hypotheses of Table 1 of Simon (1989) (p1 - p0 = 0.20), K = 5,
##                  n_max = n of Simon's optimal design for (p0, p1, alpha, beta) (n_max = NA)
simon_p0 <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
settings$runs <- rbind(
  data.frame(series = "grid", p0 = 0.10, p1 = 0.30, expand.grid(K = c(3, 5, 10), n_max = c(20, 30, 40))),
  data.frame(series = "simon", p0 = simon_p0, p1 = round(simon_p0 + 0.20, 2), K = 5, n_max = NA)
  ## Table 2 of Simon (p1 - p0 = 0.15; n up to 84, so slow):
  ## , data.frame(series = "simon15", p0 = c(simon_p0, 0.80), p1 = round(c(simon_p0, 0.80) + 0.15, 2), K = 5, n_max = NA)
)
rm(simon_p0)

## Everything that depends on the hypotheses of a run
settings$for_hypothesis <- function(s, p0, p1) {
  s$p0 <- p0; s$p1 <- p1
  s$borrow$theta1_nex    <- p1                       # alternative without borrowing
  s$borrow$theta1_min    <- p0 + 0.05                # range of theta1 with borrowing
  s$borrow$theta1_max    <- min(0.95, p1 + 0.30)
  s$futility$p_threshold <- (p0 + p1) / 2            # Bayesian futility: P(p > midpoint | data) < gamma
  s$prior$mu_mean        <- log(p0 / (1 - p0))       # EX: mu ~ N(logit(p0), mu_sd^2)
  s$prior$nex_p          <- (p0 + p1) / 2            # NEX: theta ~ N(logit(nex_p), 1/nex_p + 1/(1 - nex_p))
  s
}
## Analysis schedules: number of new patients per cohort between two analyses
settings$make_schedules <- function(n_max) {
  setNames(list(1, round(n_max / 3), round(n_max / 2)),
           c("Fully sequential", sprintf("Blocked (%d)", round(n_max / 3)), "Two-stage"))
}
## Scenarios (d = p1 - p0): 0, 1, about half, all but one and all cohorts active
## (p = p1, the others p0); 'mixed': rates spread from p0 - d/4 to p1 + d/2;
## 'strong': all at p1 + d/2 (for p0 = 0.10, p1 = 0.30: 0.05 to 0.40 and 0.40).
settings$make_scenarios <- function(K, p0, p1) {
  d <- p1 - p0; lo <- max(p0 - d / 4, p0 / 2); hi <- min(p1 + d / 2, (1 + p1) / 2)
  act <- sort(unique(c(0, 1, floor(K / 2), K - 1, K)))
  sc <- lapply(act, function(a) c(rep(p1, a), rep(p0, K - a)))
  names(sc) <- ifelse(act == 0, "0 active (global null)",
                      ifelse(act == K, sprintf("%d active (global alt)", K), sprintf("%d active", act)))
  c(sc, list("mixed" = seq(lo, hi, length.out = K), "strong" = rep(hi, K)))
}
## Umbrella set-up: two cohorts open at a time, every round(n_max / 3) patients
settings$make_entry <- function(K, n_max) floor((seq_len(K) - 1) / 2) * round(n_max / 3)

## ---- 3. Design-optimal e-values ----------------------------------------------------------------
settings$evalue_library <- file.path(this_dir, "design_optimal_evalues.R")
settings$evalue <- list(
  futility_hz = TRUE,   # stop when 1/alpha cannot be reached even by betting everything
  power_cap   = 0.95    # e-based target power = min(1 - beta, power_cap * P-max power at theta1)
)

## ---- 4. Borrowing: how the other cohorts choose theta1 -------------------------------------
##  theta1_k = pi_k * s_k + (1 - pi_k) * theta1_nex   (theta1 range: section 2)
##  s_k = mean (or quantile) of the EX predictive of p_k given the OTHER cohorts' data;
##  theta1_k = theta1_nex as long as the other cohorts have no patients.
##  The borrowing prior only chooses the alternative, so type I error <= alpha always.
##  summary, quantile and redesign are defaults; a design can set its own (section 8).
settings$borrow <- list(
  summary     = "mean",              # "mean" or "quantile"
  quantile    = 0.25,
  redesign    = "analysis",          # new theta1 at every analysis, or "entry" (once, when a cohort opens)
  theta1_step = 0.01
)

## ---- 4b. Flexible inclusion (designs with accelerate = TRUE and/or extend = TRUE, section 8) ------
##  The e-value stays valid for any data-dependent inclusion (README, 'Flexible inclusion').
##  Calendar time: 1 time unit = 1 new patient per cohort at normal speed.
settings$flex <- list(
  speed_up = 2,                          # accelerated cohort: speed_up patients per time unit
  m_speed  = sqrt(1 / settings$alpha),   # accelerate while the e-value >= m_speed (halfway to 1/alpha on the log scale)
  ext_frac = 0.5,                        # extension: n_ext = round(ext_frac * n_max) extra patients
  ext_at   = 0.75,                       # decision point: last planned analysis at or before ext_at * n_max
  cp_min   = 0.30                        # extend if conditional power planned < 1 - beta and extended >= cp_min
)

## ---- 5. Bayesian comparators (reject if P(p > p0 | data) > Delta) --------------------------
settings$futility    <- list(on = TRUE, gamma = 0.10, min_n = 10)   # threshold: section 2
settings$calibration <- list(tol = 1e-5)   # Delta calibrated under the global null (target per set-up)

## ---- 6. EXNEX priors (means: section 2) and integration grid -------------------------------------
settings$prior <- list(mu_sd = 10, tau_scale = 1)         # tau ~ HalfNormal(tau_scale)
settings$grid  <- list(theta_min = -10, theta_max = 6, theta_step = 0.02,
                       mu_min = -8, mu_max = 4, mu_step = 0.15, tau_max = 4, tau_step = 0.10)

## ---- 7. EXNEX weights pi_k --------------------------------------------------------------------
##  rule            pi_k                                                    parameters
##  none            0 (no borrowing)                                        -
##  full            1 (BHM)                                                 -
##  fixed           constant (EXNEX)                                        pi_fixed
##  hellinger_mean  mean over j != k of 1 - H(Beta_k, Beta_j), with         c_e
##                  Beta_k = Beta(y_k + 1, n_k - y_k + 1); 0 if the observed
##                  rate differs > c_e from all other cohorts (mEXNEX)
##  hellinger_max   max over j != k of 1 - H(Beta_k, Beta_j) (mEXNEX_min)   c_e (optional)
##  delayed         0 while n_k < n_burnin, then base_rule                  n_burnin, base_rule

## ---- 8. Designs ---------------------------------------------------------------------------------
##  type "mixture": the e-value of framework.tex (PRIMARY): per block the likelihood ratio mixed
##                  over the posterior of p_c given all platform data (EXNEX, weight pi_fixed, own
##                  data included), restricted to (p0, 1] and frozen at the analysis before the block
##  type "evalue":  design-optimal e-values (Baas et al.): bets from the library's dynamic program
##                  for a point alternative theta1, borrowed from the other cohorts (section 4).
##                  Within the framework: a point-mass prior on (p0, 1] (see README)
##  type "bayes":   comparators outside the framework (posterior probability > Delta)
##  Flexible inclusion (section 4b): accelerate = TRUE and/or extend = TRUE; base = the design
##  it is a variant of (same colour in the figures).
fw  <- list(type = "mixture", pi_rule = "fixed", pi_fixed = 0.5)
mex <- list(type = "evalue", evalue = "ebased", borrow = TRUE, pi_rule = "hellinger_mean", c_e = 0.10)
del <- modifyList(mex, list(pi_rule = "delayed", n_burnin = 10, base_rule = "hellinger_mean"))
ebs <- list(type = "evalue", evalue = "ebased", borrow = FALSE)
settings$designs <- list(
  ## framework.tex (primary)
  "e-mixture (framework)"           = fw,
  "e-mixture (framework, flexible)" = c(fw, accelerate = TRUE, extend = TRUE, base = "e-mixture (framework)"),
  "e-mixture (no borrowing)"        = list(type = "mixture", pi_rule = "none"),
  ## design-optimal e-values with a borrowed point alternative
  "e-GROW"                     = list(type = "evalue", evalue = "grow",   borrow = FALSE),
  "e-Pmax"                     = list(type = "evalue", evalue = "pmax",   borrow = FALSE),
  "e-ESSmin"                   = list(type = "evalue", evalue = "essmin", borrow = FALSE),
  "e-based"                    = list(type = "evalue", evalue = "ebased", borrow = FALSE),
  "e-based + mEXNEX"           = mex,
  "e-based + mEXNEX (c_e 0.2)" = modifyList(mex, list(c_e = 0.20)),
  "e-based + mEXNEX (25% q.)"  = modifyList(mex, list(summary = "quantile", quantile = 0.25)),
  "e-based + delayed"          = del,
  ## Experimental only: not included in active setups until true entry-time events are implemented.
  "e-based + EXNEX at entry (disabled)" = list(type = "evalue", evalue = "ebased", borrow = TRUE, pi_rule = "fixed",
                                                pi_fixed = 0.5, redesign = "entry"),
  "e-based (accelerated)"      = c(ebs, accelerate = TRUE, base = "e-based"),
  "e-based (extended)"         = c(ebs, extend = TRUE, base = "e-based"),
  "e-based (flexible)"         = c(ebs, accelerate = TRUE, extend = TRUE, base = "e-based"),
  "e-based + mEXNEX (flexible)"  = c(mex, accelerate = TRUE, extend = TRUE, base = "e-based + mEXNEX"),
  "e-based + delayed (flexible)" = c(del, accelerate = TRUE, extend = TRUE, base = "e-based + delayed"),
  ## comparators outside the framework
  "Bayes Umbrella"             = list(type = "bayes", pi_rule = "none"),
  "Bayes mEXNEX"               = list(type = "bayes", pi_rule = "hellinger_mean", c_e = 0.10)
)
settings$primary <- "e-mixture (framework)"     # reference for the design effects (design_effects.csv)
rm(fw, mex, del, ebs)

## ---- 9. Set-ups (entry times: section 2) ----------------------------------------------------------
settings$setups <- list(
  basket   = list(designs = c("e-mixture (framework)", "e-mixture (framework, flexible)", "e-mixture (no borrowing)",
                              "e-GROW", "e-Pmax", "e-ESSmin", "e-based", "e-based + mEXNEX", "e-based + mEXNEX (c_e 0.2)",
                              "e-based + mEXNEX (25% q.)", "e-based + delayed",
                              "e-based (accelerated)", "e-based (extended)", "e-based (flexible)",
                              "e-based + mEXNEX (flexible)", "Bayes Umbrella", "Bayes mEXNEX"),
                  calibration = "mean"),
  umbrella = list(designs = c("e-mixture (framework)", "e-mixture (framework, flexible)",
                              "e-GROW", "e-Pmax", "e-ESSmin", "e-based", "e-based + delayed",
                              "e-based (flexible)", "e-based + delayed (flexible)", "Bayes Umbrella"),
                  calibration = "per_group")   # one Delta per group of cohorts with the same analysis times
)

## ---- 10. Simulation ---------------------------------------------------------------------------------
##  n_sim can also be set from the slurm file:  export N_SIM=1000
##  Cores: the number slurm allocated (--cpus-per-task), otherwise all cores of the machine.
settings$sim <- list(
  n_sim     = as.integer(Sys.getenv("N_SIM", 1000)),  # trials per scenario (1 = check that it runs)
  n_cal     = as.integer(Sys.getenv("N_SIM", 1000)),  # trials to calibrate Delta of the Bayesian designs
  seed      = 20261124,
  n_cores   = if (.Platform$OS.type == "windows") 1L else
                as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores())),
  out_dir   = file.path(this_dir, "Result"),
  bank_dir  = file.path(this_dir, "design_banks")       # cache of the e-value designs, one file per bank
)


## =============================================================================
##  RUN: every run of the table, then an overview per series and set-up
## =============================================================================
source(settings$evalue_library)                      # design_optimal_evalues.R
source(file.path(this_dir, "additional_functions.R"))
dir.create(settings$sim$out_dir, showWarnings = FALSE)
runs <- complete_runs(settings)                      # fills in Simon's n_max and the folder names
message(sprintf("%d runs, %d trials per scenario, %d cores, output in %s",
                nrow(runs), settings$sim$n_sim, settings$sim$n_cores, settings$sim$out_dir))

all <- list()
for (i in seq_len(nrow(runs))) {
  message(sprintf("run %d of %d: %s (%s)", i, nrow(runs), runs$name[i], format(Sys.time(), "%H:%M")))
  all[[i]] <- run_one_setting(complete_settings(settings, runs[i, ]))
}
all <- do.call(rbind, all)
write.csv(all, file.path(settings$sim$out_dir, "all_settings_results.csv"), row.names = FALSE)
saveRDS(list(settings = settings, runs = runs, results = all), file.path(settings$sim$out_dir, "all_results.rds"))

for (sr in unique(runs$series)) for (su in names(settings$setups)) {
  save_plot(plot_overview(all, sr, su, settings),
            file.path(settings$sim$out_dir, sprintf("overview_%s_%s", sr, su)), width = 11, height = 15)
}
message(sprintf("done in %.1f minutes; results in %s",
                as.numeric(difftime(Sys.time(), start_time, units = "mins")), settings$sim$out_dir))
