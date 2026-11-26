## Run from anywhere:  Rscript checks/exact_checks.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings)
one_cohort <- function(p0, p1, n) {        # settings for a single cohort with horizon n
  s <- settings$for_hypothesis(settings, p0, p1); s$K <- 1; s$n_max <- n; s$sim$n_cores <- 2; s }
times_for <- function(block, n) cohort_times(block, n, 0)[[1]]
mean_n_reject <- function(ev) { d <- diff(c(0, ev$by_time$cum_reject)); sum(ev$by_time$t * d) / sum(d) }

## A. Enumerate all 2^n response sequences of one cohort and run the trial code on them
cat("A. Trial code on all 2^n sequences vs exact library result (evaluate_design)\n")
enum <- function(p0, p1, n) {
  s <- one_cohort(p0, p1, n)
  Y <- sapply(seq_len(n), function(t) bitwAnd(bitwShiftR(0:(2^n - 1), t - 1), 1L))
  C <- array(0L, c(nrow(Y), 1, n + 1)); C[, 1, -1] <- t(apply(Y, 1, cumsum))
  ny <- rowSums(Y)
  out <- list()
  for (sched in c(1, round(n / 3))) for (type in c("grow", "pmax", "essmin", "ebased")) {
    at <- times_for(sched, n)
    bank <- compact_bank(build_bank(type, n, at, FALSE, s))
    res <- run_trials_evalue(list(type = "evalue", evalue = type, borrow = FALSE), C, sched, 0, s, NULL, NULL, list(bank))
    d <- make_evalue_design(type, p1, n, at, s)
    for (p in c(p0, p1)) {
      w <- p^ny * (1 - p)^(n - ny); ev <- evaluate_design(d, p)
      out[[length(out) + 1]] <- data.frame(p0 = p0, n = n, analyses = length(at), design = type, p = p,
        reject_code = sum(w * (res$status == 1)), reject_exact = ev$reject,
        ESS_code = sum(w * res$n), ESS_exact = ev$ess)
    }
  }
  do.call(rbind, out)
}
A <- rbind(enum(0.05, 0.25, 17), enum(0.70, 0.90, 20))
print(format(A, digits = 4), row.names = FALSE)
cat("type I of the trial code <= alpha everywhere:", all(A$reject_code[A$p == A$p0] <= settings$alpha + 1e-12), "\n")
cat("largest |difference| code vs exact: reject", signif(max(abs(A$reject_code - A$reject_exact)), 2),
    " ESS", signif(max(abs(A$ESS_code - A$ESS_exact)), 2), "\n\n")

## B. Exact single-cohort results at Simon's hypotheses and n (no borrowing)
cat("B. Single cohort, exact: e-value designs vs Simon at Simon's optimal n\n")
B <- list()
for (i in which(runs$series == "simon")) {
  p0 <- runs$p0[i]; p1 <- runs$p1[i]; n <- runs$n_max[i]; s <- one_cohort(p0, p1, n)
  sd <- simon_design(p0, p1, settings$alpha, settings$beta, n = n)
  B[[length(B) + 1]] <- data.frame(hyp = sprintf("%.2f/%.2f", p0, p1), n = n, design = "Simon", schedule = "two-stage",
    typeI = sd$typeI, power = sd$power, early = 0, n_rej = n, ESS0 = sd$EN0,
    ESS1 = sd$n1 + (n - sd$n1) * (1 - pbinom(sd$r1, sd$n1, p1)))
  for (sched in c(1, round(n / 3), round(n / 2))) for (type in c("grow", "pmax", "essmin", "ebased")) {
    d <- make_evalue_design(type, p1, n, times_for(sched, n), s); e0 <- evaluate_design(d, p0); e1 <- evaluate_design(d, p1)
    B[[length(B) + 1]] <- data.frame(hyp = sprintf("%.2f/%.2f", p0, p1), n = n, design = paste0("e-", type),
      schedule = c("sequential", "blocked", "two-stage")[match(sched, c(1, round(n / 3), round(n / 2)))],
      typeI = e0$reject, power = e1$reject, early = e1$by_time$cum_reject[n], n_rej = mean_n_reject(e1),
      ESS0 = e0$ess, ESS1 = e1$ess)
  }
}
B <- do.call(rbind, B); saveRDS(B, file.path(out_dir, "exact_simon.rds"))
print(format(B, digits = 3), row.names = FALSE)
