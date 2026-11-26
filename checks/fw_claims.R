## Run from anywhere:  Rscript checks/fw_claims.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
p0 <- 0.10; a <- 0.05; thr <- 1 / a; N <- 20
Y <- sapply(seq_len(N), function(t) bitwAnd(bitwShiftR(0:(2^N - 1), t - 1), 1L)); S <- nrow(Y)
cs <- t(apply(Y, 1, cumsum))
## single basket, analysis after every patient; prior for patient t chosen from data of patients < t (or <= t)
run_seq <- function(q_of) {                     # q_of(t, y_before, Y_t): predictive mean used for patient t
  E <- rep(1, S); hit <- rep(FALSE, S)
  for (t in seq_len(N)) { yb <- if (t == 1) rep(0L, S) else cs[, t - 1]; q <- q_of(t, yb, Y[, t])
    E <- E * ifelse(Y[, t] == 1, q / p0, (1 - q) / (1 - p0)); hit <- hit | E >= thr }
  hit
}
P <- function(hit, p) sum(p^cs[, N] * (1 - p)^(N - cs[, N]) * hit)
## posterior Beta(1 + y, 1 + m - y) of p given m patients with y responses; mean restricted to (p0, 1] or not
pm <- function(al, be, restrict) if (restrict) {
  num <- exp(lbeta(al + 1, be) - lbeta(al, be)) * pbeta(p0, al + 1, be, lower.tail = FALSE); num / pbeta(p0, al, be, lower.tail = FALSE)
} else al / (al + be)
## the mixture over a frozen prior = product of its predictive means, so these are framework block factors
rows <- list()
add <- function(what, hit) rows[[length(rows) + 1]] <<- data.frame(construction = what, typeI_at_p0 = P(hit, p0), typeI_at_0.02 = P(hit, 0.02), power_at_0.30 = P(hit, 0.30))
add("framework: posterior given earlier data, restricted to (p0, 1]", run_seq(function(t, yb, yt) pm(1 + yb, t - yb, TRUE)))
add("NOT restricted: full posterior given earlier data",             run_seq(function(t, yb, yt) pm(1 + yb, t - yb, FALSE)))
add("NOT frozen: restricted posterior incl. the patient's own outcome", run_seq(function(t, yb, yt) pm(1 + yb + yt, t - yb + 1 - yt, TRUE)))
mc <- function(seed, H, own = FALSE) function(t, yb, yt) {      # Monte Carlo average over draws > p0 (eq. mc); 1 if none (H' = 0)
  y_use <- yb + if (own) yt else 0; m_use <- t - 1 + own         # draws depend only on (t, y): one set per state
  tab <- vapply(0:m_use, function(y) { set.seed(seed * 1e6 + t * 1e3 + y); d <- rbeta(H, 1 + y, 1 + m_use - y); d <- d[d > p0]
    if (!length(d)) p0 else mean(d) }, 0)
  tab[y_use + 1]
}
## Monte Carlo: the factor averages the block likelihood ratio over the draws; for one patient that is the LR at the mean draw
for (sd in 1:3) add(sprintf("Monte Carlo, 200 frozen draws (seed %d)", sd), run_seq(mc(sd, 200)))
add("Monte Carlo, draws re-made after seeing the outcome", run_seq(mc(1, 200, own = TRUE)))
res <- do.call(rbind, rows); print(format(res, digits = 3), row.names = FALSE)

## H' = 0: probability that none of H posterior draws exceeds p0, for null cohorts in the EXNEX model
runs <- complete_runs(settings); s <- complete_settings(settings, runs[runs$name == "grid_p10-30_K5_n30", ]); eng <- build_engine(s)
st <- rbind(c(3, 3, 3, 3, 3), c(1, 1, 1, 1, 1), c(0, 1, 0, 1, 0), c(0, 0, 0, 0, 0))
nn <- matrix(30, 4, 5); po <- posterior_exnex(eng, st, nn, matrix(0.5, 4, 5))
cat("\nH' = 0 (no posterior draw above p0 = 0.10), five null cohorts of 30 patients, EXNEX (pi = 0.5):\n")
print(data.frame(responses_per_cohort = apply(st, 1, paste, collapse = ","), P_p_above_p0 = signif(po$pp_eff[, 1], 3),
                 P_Hprime0_H1000 = signif((1 - po$pp_eff[, 1])^1000, 3), P_Hprime0_H10000 = signif((1 - po$pp_eff[, 1])^10000, 3)), row.names = FALSE)

## Design-optimal bets inside a block: a frozen prior raises the bet after a response and lowers it after a non-response
s2 <- complete_settings(settings, runs[runs$name == "simon_p10-30_K5_n29", ]); b2 <- get_banks(s2)
for (sn in names(s2$schedules)) {
  bk <- banks_for(s2$designs[["e-based"]], s2$schedules[[sn]], s2$setups$basket, s2, b2)[[1]]; at <- cohort_times(s2$schedules[[sn]], s2$n_max, 0)[[1]]
  g <- bk$m_grid; top <- length(g); viol <- tot <- 0
  for (t in setdiff(seq_len(29 - 1), at)) {      # t and t + 1 in the same block (no analysis after t)
    b <- bk$b_val0[as.integer(bk$bet_idx[, t, 1]) + 1]; ok <- which(b > 0 & seq_along(g) < top)
    up <- findInterval(g[ok] * (1 + b[ok] * (1 / p0 - 1)), g); b_up <- bk$b_val0[as.integer(bk$bet_idx[cbind(up, t + 1, 1)]) + 1]
    keep <- up < top; viol <- viol + sum(b_up[keep] < b[ok][keep] - 1e-12); tot <- tot + sum(keep) }
  cat(sprintf("%s: within-block states where e-based LOWERS its bet after a response: %d of %d\n", sn, viol, tot))
}
saveRDS(res, file.path(out_dir, "fw_claims.rds"))
