## Run from anywhere:  Rscript checks/fw_claims2.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
p0 <- 0.10; thr <- 20
## 1. Restriction to (p0, 1]: the unrestricted posterior-predictive product is the Bayes factor B(1+y,1+n-y)/(p0^y (1-p0)^(n-y)),
##    a function of (n, y): exact first-passage probability by dynamic programming over (n, y), horizon up to 150
first_pass <- function(p, nmax, restricted) {
  mass <- c(1); hit <- 0
  ## restricted product depends on the path; use the unrestricted Bayes factor (path-free) and, for the restricted case,
  ## the (exact, path-free) mixture over the uniform prior restricted to (p0, 1]: B-ratio x tail ratio
  E <- function(n, y) { bf <- exp(lbeta(1 + y, 1 + n - y) - y * log(p0) - (n - y) * log(1 - p0))
    if (restricted) bf * pbeta(p0, 1 + y, 1 + n - y, lower.tail = FALSE) / pbeta(p0, 1, 1, lower.tail = FALSE) else bf }
  out <- numeric(nmax)
  for (n in 1:nmax) { new <- numeric(n + 1); new[1:n] <- mass * (1 - p); new[2:(n + 1)] <- new[2:(n + 1)] + mass * p
    cross <- E(n, 0:n) >= thr; hit <- hit + sum(new[cross]); new[cross] <- 0; mass <- new; out[n] <- hit }
  out
}
cat("1. Prior on (p0, 1] or not: P(evidence >= 20 by n) for the uniform-prior mixture (fully sequential)\n")
tab <- do.call(rbind, lapply(c(0.10, 0.05, 0.02), function(p) data.frame(true_p = p,
  restricted_n50 = first_pass(p, 50, TRUE)[50], restricted_n150 = first_pass(p, 150, TRUE)[150],
  unrestricted_n50 = first_pass(p, 50, FALSE)[50], unrestricted_n150 = first_pass(p, 150, FALSE)[150])))
print(format(tab, digits = 3), row.names = FALSE)

## 2. Conditional null of framework.tex: response probability may drift / depend on the history, as long as <= p0
s <- settings$for_hypothesis(settings, 0.10, 0.30); s$K <- 1; s$n_max <- 10; s$n_ext <- 5L
s$sim$bank_dir <- file.path(out_dir, "banks_k2"); s$sim$n_cores <- 2; s$setups$basket$entry <- 0
N <- 15; Y <- sapply(seq_len(N), function(t) bitwAnd(bitwShiftR(0:(2^N - 1), t - 1), 1L)); S <- nrow(Y)
C <- array(0L, c(S, 1, N + 1)); C[, 1, -1] <- t(apply(Y, 1, cumsum)); eng <- build_engine(s)
w_of <- function(pmat) apply(ifelse(Y == 1, pmat, 1 - pmat), 1, prod)
prev <- cbind(0L, Y[, -N])
drifts <- list("constant p0" = matrix(0.10, S, N),
               "0.10 after a response, else 0.05" = ifelse(prev == 1, 0.10, 0.05),
               "rising from 0.02 to 0.10" = matrix(seq(0.02, 0.10, length.out = N), S, N, byrow = TRUE),
               "0.10 while 2+ responses so far, else 0.03" = ifelse(t(apply(prev, 1, cumsum)) >= 2, 0.10, 0.03))
ws <- lapply(drifts, w_of)
dz <- list("e-mixture (no borrowing)" = list(type = "mixture", pi_rule = "none"),
           "e-based" = list(type = "evalue", evalue = "ebased", borrow = FALSE),
           "e-based (flexible)" = list(type = "evalue", evalue = "ebased", borrow = FALSE, accelerate = TRUE, extend = TRUE))
out <- list()
for (block in c(1, 3)) {
  at <- cohort_times(block, 10, 0)[[1]]; td <- ext_point(at, 10, s)
  bk <- compact_bank(build_bank("ebased", 10, at, FALSE, s))
  r_ext <- 10 - td + 5; ek <- c(compact_bank(build_bank("pmax", r_ext, seq_len(r_ext), FALSE, s)), list(td = td))
  for (dn in names(dz)) { d <- dz[[dn]]
    r <- if (d$type == "mixture") run_trials_mixture(d, C, block, 0, s, eng, new_cache()) else run_trials_evalue(d, C, block, 0, s, eng, NULL, list(bk), list(ek))
    for (nm in names(ws)) out[[length(out) + 1]] <- data.frame(schedule = c("sequential", "blocked (3)")[match(block, c(1, 3))], design = dn,
      response_probability = nm, typeI = sum(ws[[nm]] * (r$status[, 1] == 1))) }
}
res <- do.call(rbind, out); w <- reshape(res, idvar = c("schedule", "design"), timevar = "response_probability", direction = "wide")
names(w) <- sub("typeI.", "", names(w)); print(format(w, digits = 3), row.names = FALSE)
cat("all <= alpha:", all(res$typeI <= 0.05), "; largest:", signif(max(res$typeI), 4), "\n")
saveRDS(list(restrict = tab, drift = w), file.path(out_dir, "fw_claims2.rds"))
