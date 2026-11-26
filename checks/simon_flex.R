## Run from anywhere:  Rscript checks/simon_flex.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings); t0 <- Sys.time(); before <- length(list.files(settings$sim$bank_dir))
for (i in seq_len(nrow(runs))) invisible(get_banks(complete_settings(settings, runs[i, ])))
cat("extension designs built:", length(list.files(settings$sim$bank_dir)) - before, "in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n\n")

## Exact single-cohort recursion on the e-value grid (planned analyses, normal speed), with extension
exact_flex <- function(bk, ek, at, p, s, extend, block) {
  g <- bk$m_grid; top <- length(g); n <- s$n_max[1]; ne <- if (extend) s$n_ext[1] else 0; p0 <- s$p0
  td <- if (extend) ek$td else Inf; A_all <- sort(unique(c(at, if (extend) c(seq(n + block, n + ne, by = block), n + ne))))
  step <- function(m, b) { up <- findInterval(g * (1 + b * (1 / p0 - 1)), g); dn <- findInterval(g * (1 - b), g); up[top] <- dn[top] <- top
                           .accumulate(c(up, dn), c(p * m, (1 - p) * m), top) }
  mA <- numeric(top); mA[which(g == 1)] <- 1; mB <- numeric(top); hzA <- n + ne; decided <- FALSE
  rej <- ess <- ext_p <- rej_after <- 0
  for (t in seq_len(n + ne)) {
    if (t <= n) { ess <- ess + sum(mA); mA <- step(mA, bk$b_val0[as.integer(bk$bet_idx[, t, 1]) + 1]) }
    if (t > td) { ess <- ess + sum(mB); mB <- step(mB, ek$b_val0[as.integer(ek$bet_idx[, t - td, 1]) + 1]) }
    if (!(t %in% A_all)) next
    rej <- rej + mA[top] + mB[top]; if (t > n) rej_after <- rej_after + mB[top]; mA[top] <- mB[top] <- 0
    if (extend && !decided && t >= td) {
      plan <- if (t < n) cond_power(bk, p0, t)[, 1] else 0; more <- cond_power(ek, p0, t - td)[, 1]
      go <- plan < 1 - s$beta & more >= s$flex$cp_min; ext_p <- sum(mA[go]); mB <- mB + mA * go; mA <- mA * !go
      decided <- TRUE; hzA <- n }
    if (t < n) mA[as.integer(bk$stop[, t + 1, 1]) == 1L | g < p0^(hzA - t) / s$alpha] <- 0 else mA[] <- 0
    if (t > td && t < n + ne) mB[g < p0^(n + ne - t) / s$alpha] <- 0
  }
  c(reject = rej, ESS = ess, extended = ext_p, reject_after = rej_after)
}
tab <- list()
for (i in which(runs$series == "simon")) {
  s <- complete_settings(settings, runs[i, ]); b <- get_banks(s); s1 <- s; s1$K <- 1; s1$n_max <- s$n_max[1]; s1$n_ext <- s$n_ext[1]
  for (sn in names(s$schedules)) {
    block <- s$schedules[[sn]]; su <- s$setups$basket
    bk <- banks_for(s$designs[["e-based"]], block, su, s, b)[[1]]
    ek <- ext_banks_for(s$designs[["e-based (extended)"]], block, su, s, b)[[1]]
    at <- cohort_times(block, s$n_max, su$entry)[[1]]
    f0 <- exact_flex(bk, ek, at, s$p0, s1, FALSE, block); f1 <- exact_flex(bk, ek, at, s$p1, s1, FALSE, block)
    e0 <- exact_flex(bk, ek, at, s$p0, s1, TRUE, block); e1 <- exact_flex(bk, ek, at, s$p1, s1, TRUE, block)
    tab[[length(tab) + 1]] <- data.frame(hyp = sprintf("%.2f/%.2f", s$p0, s$p1), n_max = s$n_max[1], n_ext = s$n_ext[1], td = ek$td,
      schedule = c("sequential", "blocked", "two-stage")[match(sn, names(s$schedules))],
      typeI = f0[["reject"]], typeI_ext = e0[["reject"]], power = f1[["reject"]], power_ext = e1[["reject"]],
      ESS0 = f0[["ESS"]], ESS0_ext = e0[["ESS"]], ESS1 = f1[["ESS"]], ESS1_ext = e1[["ESS"]],
      P_ext_H0 = e0[["extended"]], P_ext_H1 = e1[["extended"]])
  }
}
tab <- do.call(rbind, tab); saveRDS(tab, file.path(out_dir, "simon_flex.rds"))
print(format(tab, digits = 3), row.names = FALSE)
cat("type I <= alpha everywhere (grid):", all(c(tab$typeI, tab$typeI_ext) <= settings$alpha), "\n\n")

## Naive extension of the exact binomial test: test at n_max; if not significant but within 2 responses of
## the critical value, add n_ext patients and test all patients at the critical value of a test with n_max + n_ext
cat("Naive extension of the binomial test (same rule, no correction):\n")
nb <- do.call(rbind, lapply(which(runs$series == "simon"), function(i) {
  p0 <- runs$p0[i]; n <- runs$n_max[i]; ne <- round(settings$flex$ext_frac * n); a <- settings$alpha
  crit <- function(N) min(which(pbinom(0:N - 1, N, p0, lower.tail = FALSE) <= a)) - 1
  r <- crit(n); r2 <- crit(n + ne); x <- max(0, r - 2):(r - 1)
  extra <- sum(dbinom(x, n, p0) * pbinom(r2 - x - 1, ne, p0, lower.tail = FALSE))
  data.frame(hyp = sprintf("%.2f/%.2f", p0, runs$p1[i]), n_max = n, n_ext = ne, typeI_fixed = pbinom(r - 1, n, p0, lower.tail = FALSE),
             typeI_with_extension = pbinom(r - 1, n, p0, lower.tail = FALSE) + extra) }))
print(format(nb, digits = 3), row.names = FALSE); saveRDS(nb, file.path(out_dir, "naive_binomial.rds"))
