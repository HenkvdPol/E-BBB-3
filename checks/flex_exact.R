## Run from anywhere:  Rscript checks/flex_exact.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
## Single cohort, every possible response sequence (2^(n_max + n_ext)): exact type I error and power of
## the flexible designs, and an independent scalar implementation of the same rules on every 8th sequence.
source(file.path(check_dir, "load_settings.R"))
one <- function(p0, p1, n, ne) { s <- settings$for_hypothesis(settings, p0, p1); s$K <- 1; s$n_max <- n; s$n_ext <- ne; s$sim$n_cores <- 2; s }
scalar <- function(y, bk, ek, block, acc, ext, s, cpp, cpe) {   # straightforward re-implementation, one sequence
  n <- s$n_max; ne <- if (ext) s$n_ext else 0; p0 <- s$p0; thr <- 1 / s$alpha; g <- bk$m_grid; top <- length(g)
  M <- 1; i <- which(g == 1); hit <- FALSE; t <- 0; cap <- n; hzn <- n + ne; extended <- FALSE; decided <- FALSE; speed <- 1
  td <- if (ext) ek$td else NA
  planned <- sort(unique(c(seq(block, n, by = block), n))); tau_prev <- 0; l <- 0
  repeat {
    l <- l + 1; tau <- if (l <= length(planned)) planned[l] else tau_prev + block
    target <- min(cap, t + speed * (tau - tau_prev)); tnow <- if (target >= cap) tau_prev + (cap - t) / speed else tau
    tau_prev <- tau; if (target <= t) next
    while (t < target) { t <- t + 1
      b <- if (!extended) bk$b_val0[as.integer(bk$bet_idx[i, t, 1]) + 1] else ek$b_val0[as.integer(ek$bet_idx[i, t - td, 1]) + 1]
      if (hit) b <- 0
      f <- if (y[t] == 1) 1 + b * (1 / p0 - 1) else 1 - b
      M <- M * f; if (i != top) i <- findInterval(g[i] * f, g); if (M >= thr) hit <- TRUE }
    if (hit) return(c(1, t, tnow, extended))
    if (ext && !decided && t >= td) {
      plan <- if (t < n) cpp[[t + 1]][i] else 0; more <- cpe[[t - td + 1]][i]; decided <- TRUE
      if (plan < 1 - s$beta && more >= s$flex$cp_min) { extended <- TRUE; cap <- n + ne } else hzn <- n }
    if (t < cap) { stp <- if (!extended) as.integer(bk$stop[i, t + 1, 1]) == 1 else as.integer(ek$stop[i, t - td + 1, 1]) == 1
                   if (stp || g[i] < p0^(hzn - t) / s$alpha) return(c(2, t, tnow, extended)) }
    if (t >= cap) return(c(3, t, tnow, extended))
    if (acc) speed <- if (M >= s$flex$m_speed) s$flex$speed_up else 1
  }
}
out <- list(); mism <- 0
for (h in list(c(0.10, 0.30), c(0.30, 0.50))) {
  n <- 12; ne <- 6; s <- one(h[1], h[2], n, ne); N <- n + ne
  Y <- sapply(seq_len(N), function(t) bitwAnd(bitwShiftR(0:(2^N - 1), t - 1), 1L))
  C <- array(0L, c(nrow(Y), 1, N + 1)); C[, 1, -1] <- t(apply(Y, 1, cumsum)); ny <- rowSums(Y)
  for (block in c(1, 4)) {
    at <- cohort_times(block, n, 0)[[1]]; td <- ext_point(at, n, s); r <- n - td + ne
    bk <- compact_bank(build_bank("ebased", n, at, FALSE, s))
    ek <- c(compact_bank(build_bank("pmax", r, seq_len(r), FALSE, s)), list(td = td))
    cpp <- lapply(0:n, function(t) cond_power(bk, s$p0, t)[, 1]); cpe <- lapply(0:r, function(t) cond_power(ek, s$p0, t)[, 1])
    for (v in list(c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE))) {
      d <- list(type = "evalue", evalue = "ebased", borrow = FALSE, accelerate = v[1], extend = v[2])
      r <- run_trials_evalue(d, C, block, 0, s, NULL, NULL, list(bk), list(ek))
      sub <- seq(1, nrow(Y), by = 8)
      sc <- t(vapply(sub, function(i) scalar(Y[i, ], bk, ek, block, v[1], v[2], s, cpp, cpe), numeric(4)))
      mism <- mism + sum(sc[, 1] != r$status[sub] | sc[, 2] != r$n[sub] | abs(sc[, 3] - r$time[sub]) > 1e-9 | sc[, 4] != r$extended[sub])
      for (p in c(h[1] - 0.05, h[1], h[2])) { w <- p^ny * (1 - p)^(N - ny); rej <- r$status == 1
        out[[length(out) + 1]] <- data.frame(p0 = h[1], schedule = c("sequential", "blocked (4)")[match(block, c(1, 4))],
          design = c("e-based", "accelerated", "extended", "flexible")[1 + v[1] + 2 * v[2]], p = p,
          reject = sum(w * rej), E_n = sum(w * r$n), P_extended = sum(w * r$extended), P_reject_after_n_max = sum(w * (rej & r$n > n)),
          time_to_reject = sum(w * rej * r$time) / sum(w * rej), max_n = max(r$n)) }
    }
  }
}
res <- do.call(rbind, out)
print(format(res, digits = 3), row.names = FALSE)
cat("\nscalar re-implementation vs engine, mismatching sequences:", mism, "\n")
cat("type I error (p <= p0) of every design <= alpha:", all(res$reject[res$p <= res$p0 + 1e-9] <= settings$alpha), 
    "; largest:", signif(max(res$reject[res$p <= res$p0 + 1e-9]), 4), "\n")
saveRDS(res, file.path(out_dir, "flex_exact.rds"))
