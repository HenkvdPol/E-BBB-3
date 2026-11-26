## =============================================================================
##  run_framework_study.R -- simulation study of framework.tex: one e-process, four priors
##
##  Run:    N_SIM=1000 Rscript run_framework_study.R      (e.g. the only line of a slurm job)
##  Smoke:  N_SIM=2 Rscript run_framework_study.R
##  Open-ended variant (no maximum sample size, anytime-valid futility; docs/open-ended-design.md):
##          OPEN_ENDED=1 N_SIM=2 Rscript run_framework_study.R   -> Result/framework_study_open/
##  Output: Result/framework_study/  baskets.csv, platform.csv, fig_*.pdf/png, settings.rds, session.txt
##
##  Each basket is tested with the framework e-process (R/eprocess.R). At every analysis the
##  borrowing model is fitted to all platform data, its posterior of p_c is restricted to (p0, 1] and
##  frozen as the prior of the next block (R/simulate.R). Borrowing models compared:
##    none                         NEX prior, the basket's own data only
##    one-way EXNEX                the earlier benchmark: all baskets one exchangeable list
##    one-way EXNEX, scaled Beta   the same, with the earlier prototype's moment-matched prior
##    two-way EXNEX                framework.tex eq. (9): tumour type x biomarker group
##  Settings and decisions: docs/implementation-plan.md. Per-basket type I error and the
##  platform-wide FWER are reported separately; only the former is controlled at alpha.
## =============================================================================
start_time <- Sys.time()
this_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
for (f in list.files(file.path(this_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

## ---- 1. Hypotheses, error rate and trial ------------------------------------------------------
settings <- list(
  p0 = 0.10, p1 = 0.30,          # H0: p_c <= p0 per basket; p1 only centres the NEX prior and the scenarios
  alpha = 0.05,                  # per basket; e-value threshold 1/alpha = 20
  J = 3, L = 3,                  # tumour types (rows) x biomarker groups (columns)
  n_max = 30, look_every = 5,    # time units; balanced recruitment: 1 patient per basket per unit
  n_sim = as.integer(Sys.getenv("N_SIM", 1000)), seed = 20261126,
  n_cores = if (.Platform$OS.type == "windows") 1L else as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores())),
  mcmc = list(burn = 500, keep = 1000, chains = 2, chunk = 150),   # two-way sampler (validation/twoway_validation.R)
  out_dir = file.path(this_dir, "Result", "framework_study"),
  ## open-ended variant (OPEN_ENDED=1): no maximum sample size per basket; a basket stops when it is
  ## rejected (E >= 1/alpha) or for futility (E^f >= 1/beta, testing p >= p1). n_cap only ends the
  ## computation; baskets still open at n_cap are reported as open (column 'open').
  open_ended = Sys.getenv("OPEN_ENDED") == "1", beta = 0.20, n_cap = 150)
if (settings$open_ended) settings$out_dir <- file.path(this_dir, "Result", "framework_study_open")
pbar <- (settings$p0 + settings$p1) / 2

## ---- 2. Recruitment: balanced, and unbalanced with a missing cell and a late basket -----------
cells <- function(m) as.vector(t(m))                    # 3 x 3 matrix by rows -> basket order (1,1) (1,2) ... (3,3)
recruitment <- list(
  balanced   = list(layout = grid_layout(3, 3), rate = 1, entry = 0),
  unbalanced = list(layout = grid_layout(3, 3, missing = list(c(3, 3))),          # cell (3,3) not in the trial
                    rate = cells(matrix(c(1, 1, 0.5, 1, 0.5, 0.5, 0.5, 0.5, NA), 3, byrow = TRUE))[-9],
                    entry = cells(matrix(c(0, 0, 0, 0, 0, 20, 0, 0, NA), 3, byrow = TRUE))[-9]))  # (2,3) opens at t = 20

## ---- 3. Scenarios: true response rates per basket (rows: tumour types, columns: biomarkers) ---
s0 <- settings$p0; s1 <- settings$p1
scen <- function(...) cells(matrix(c(...), 3, byrow = TRUE))
scenarios <- list(
  "global null"                 = scen(s0, s0, s0,  s0, s0, s0,  s0, s0, s0),
  "global alternative"          = scen(s1, s1, s1,  s1, s1, s1,  s1, s1, s1),
  "active tumour type (row 1)"  = scen(s1, s1, s1,  s0, s0, s0,  s0, s0, s0),
  "active biomarker (column 1)" = scen(s1, s0, s0,  s1, s0, s0,  s1, s0, s0),
  "one active basket"           = scen(s1, s0, s0,  s0, s0, s0,  s0, s0, s0),
  "row and column active"       = scen(s1, s1, s1,  s1, s0, s0,  s1, s0, s0),
  "additive row + column"       = cells(plogis(qlogis(s0) + outer(c(0, 0.5, 1), c(0, 0.5, 1), "+"))),   # 0.10 ... 0.45
  "interaction (diagonal)"      = scen(s1, s0, s0,  s0, s1, s0,  s0, s0, s1),
  "mixed active and null"       = scen(s1, s1, s0,  s1, s0, s0,  s0, s0, s1),
  "strong heterogeneity"        = scen(0.05, s0, 0.45,  s0, s1, s0,  0.60, s0, 0.05))

## ---- 4. Borrowing models ---------------------------------------------------------------------
make_designs <- function(rec) {
  grid <- theta_grid(settings$p0)
  base <- list(layout = rec$layout, n_max = settings$n_max, rate = rec$rate, entry = rec$entry,
               look_every = settings$look_every, alpha = settings$alpha)
  if (settings$open_ended) base <- modifyList(base, list(n_max = Inf, n_cap = settings$n_cap, futility = "evalue",
                                                          p1 = settings$p1, beta = settings$beta))
  ow <- engine_oneway(grid, pbar, settings$n_max)
  list("none"                       = list(engine = engine_none(grid, pbar), design = base),
       "one-way EXNEX"              = list(engine = ow, design = base),
       "one-way EXNEX, scaled Beta" = list(engine = ow, design = c(base, prior_form = "beta")),
       "two-way EXNEX"              = list(engine = engine_twoway(grid, rec$layout, pbar, settings$n_max, mcmc = settings$mcmc),
                                           design = base))
}

## =============================================================================
##  RUN
## =============================================================================
dir.create(settings$out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("framework study%s: %d trials per cell, %d cores, output %s\n", if (settings$open_ended) " (open-ended)" else "",
            settings$n_sim, settings$n_cores, settings$out_dir))
chunks <- parallel::splitIndices(settings$n_sim, min(settings$n_cores, settings$n_sim))
all <- list(); timing <- list()
for (ri in seq_along(recruitment)) {
  rec <- recruitment[[ri]]; designs <- make_designs(rec)
  for (si in seq_along(scenarios)) {
    p <- scenarios[[si]]; if (length(p) > nrow(rec$layout)) p <- p[-9]          # unbalanced: no cell (3,3)
    for (dn in names(designs)) {
      t0 <- Sys.time()
      ## same seed for every borrowing model: the models see the same patients (paired comparison)
      res <- parallel::mclapply(seq_along(chunks), function(ch) {
        d <- simulate_trials(p, designs[[dn]]$design, designs[[dn]]$engine, length(chunks[[ch]]),
                             seed = settings$seed + 1000 * (100 * ri + si) + ch)
        d$trial <- chunks[[ch]][d$trial]; d }, mc.cores = settings$n_cores, mc.preschedule = FALSE)
      bad <- vapply(res, inherits, TRUE, "try-error"); if (any(bad)) stop(res[[which(bad)[1]]])
      d <- do.call(rbind, res); d$recruitment <- names(recruitment)[ri]; d$scenario <- names(scenarios)[si]; d$engine <- dn
      all[[length(all) + 1]] <- d
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      timing[[length(timing) + 1]] <- data.frame(recruitment = names(recruitment)[ri], scenario = names(scenarios)[si], engine = dn, seconds = secs)
      cat(sprintf("  %-10s | %-28s | %-27s %7.1f s\n", names(recruitment)[ri], names(scenarios)[si], dn, secs))
    }
  }
}
d <- do.call(rbind, all)
lev <- function(x, l) factor(x, levels = l)
d$scenario <- lev(d$scenario, names(scenarios)); d$engine <- lev(d$engine, names(designs))
b <- summarise_baskets(d); pl <- summarise_platform(d)
write.csv(b, file.path(settings$out_dir, "baskets.csv"), row.names = FALSE)
write.csv(pl, file.path(settings$out_dir, "platform.csv"), row.names = FALSE)
write.csv(do.call(rbind, timing), file.path(settings$out_dir, "timing.csv"), row.names = FALSE)
saveRDS(list(settings = settings, scenarios = scenarios, recruitment = recruitment, trials = d), file.path(settings$out_dir, "results.rds"))
for (r in names(recruitment)) {
  save_plot(plot_baskets(b, r, settings$alpha), file.path(settings$out_dir, paste0("fig_baskets_", r)), 11, 14)
  save_plot(plot_tradeoff(pl, r), file.path(settings$out_dir, paste0("fig_tradeoff_", r)), 11, 8)
}
sha <- tryCatch(system2("git", c("-C", this_dir, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA)
writeLines(c(paste("git", sha), sprintf("N_SIM = %d, cores = %d, seed = %d", settings$n_sim, settings$n_cores, settings$seed),
             capture.output(sessionInfo())), file.path(settings$out_dir, "session.txt"))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
