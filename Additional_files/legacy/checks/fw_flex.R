## Run from anywhere:  Rscript checks/fw_flex.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
s <- settings$for_hypothesis(settings, 0.10, 0.30); s$K <- 2; s$n_max <- c(10, 10); s$n_ext <- c(5L, 5L)
s$schedules <- list("Fully sequential" = 1, "Blocked (3)" = 3); eng <- build_engine(s)
N <- 15; Y1 <- sapply(seq_len(N), function(t) bitwAnd(bitwShiftR(0:(2^N - 1), t - 1), 1L)); S <- nrow(Y1); ny <- rowSums(Y1)
lag <- function(Y) cbind(0L, Y[, -N])
paths <- list("all responses" = matrix(1L, S, N), "no responses" = matrix(0L, S, N), "copies cohort 1 (lag 1)" = lag(Y1),
              "opposite of cohort 1 (lag 1)" = 1L - lag(Y1))
dz <- s$designs[c("e-mixture (framework)", "e-mixture (framework, flexible)")]
out <- list()
for (pn in names(paths)) {
  C <- array(0L, c(S, 2, N + 1)); C[, 1, -1] <- t(apply(Y1, 1, cumsum)); C[, 2, -1] <- t(apply(paths[[pn]], 1, cumsum))
  for (sn in names(s$schedules)) for (dn in names(dz)) {
    r <- run_trials_mixture(dz[[dn]], C, s$schedules[[sn]], c(0, 0), s, eng, new_cache()); rej <- r$status[, 1] == 1
    w <- function(p) p^ny * (1 - p)^(N - ny)
    out[[length(out) + 1]] <- data.frame(other_cohort = pn, schedule = sn, design = dn, typeI_p0 = sum(w(0.1) * rej),
      typeI_0.05 = sum(w(0.05) * rej), power_0.30 = sum(w(0.3) * rej), P_extended_0.30 = sum(w(0.3) * r$extended[, 1]),
      n_max_used = max(r$n[, 1]), time_to_reject_0.30 = sum(w(0.3) * rej * r$time[, 1]) / sum(w(0.3) * rej))
  }
}
res <- do.call(rbind, out); print(format(res, digits = 3), row.names = FALSE)
cat("largest type I error over all paths, schedules and designs:", signif(max(res$typeI_p0, res$typeI_0.05), 4), "\n\n")
## single cohort, drifting response probability <= p0
s1 <- s; s1$K <- 1; s1$n_max <- 10; s1$n_ext <- 5L
C1 <- array(0L, c(S, 1, N + 1)); C1[, 1, -1] <- t(apply(Y1, 1, cumsum)); prev <- lag(Y1)
drifts <- list("constant p0" = matrix(0.10, S, N), "0.10 after a response, else 0.05" = ifelse(prev == 1, 0.10, 0.05),
               "rising 0.02 -> 0.10" = matrix(seq(0.02, 0.10, length.out = N), S, N, byrow = TRUE),
               "0.10 once 2+ responses, else 0.03" = ifelse(t(apply(prev, 1, cumsum)) >= 2, 0.10, 0.03))
d1 <- list("framework, flexible" = dz[[2]], "no borrowing, flexible" = list(type = "mixture", pi_rule = "none", accelerate = TRUE, extend = TRUE))
dr <- list()
for (sn in names(s$schedules)) for (dn in names(d1)) { r <- run_trials_mixture(d1[[dn]], C1, s$schedules[[sn]], 0, s1, eng, new_cache())
  dr[[length(dr) + 1]] <- data.frame(schedule = sn, design = dn, t(vapply(drifts, function(pm) sum(apply(ifelse(Y1 == 1, pm, 1 - pm), 1, prod) * (r$status[, 1] == 1)), 0))) }
dr <- do.call(rbind, dr); print(format(dr, digits = 3), row.names = FALSE)
cat("largest type I error under drift:", signif(max(dr[, -(1:2)]), 4), "\n")
saveRDS(list(adv = res, drift = dr), file.path(out_dir, "fw_flex.rds"))
