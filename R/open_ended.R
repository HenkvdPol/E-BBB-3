## =============================================================================
##  open_ended.R -- exact operating characteristics of one basket, with or without a maximum
##                  sample size (feasibility of open-ended designs; docs/open-ended-design.md)
##
##  Telescoping. Suppose the prior of every block is the posterior from ONE fixed starting
##  distribution pi_0, restricted to (p0, 1]. That holds without borrowing, and with borrowing from
##  baskets whose data no longer change (pi_0 = the borrowed predictive given their data). Then
##      prod_k e_k = sum_g w0_g (p_g / p0)^y ((1 - p_g) / (1 - p0))^(n - y),   w0 = pi_0 restricted,
##  because each factor is a ratio of consecutive restricted normalising constants. The same holds
##  for the futility e-process with pi_0 restricted to [0, p1). E and E^f then depend on (y, n) only,
##  so a dynamic programme over y gives the exact probabilities of rejection, futility stop and of
##  still being open after any number of patients. No simulation is involved.
##  (With borrowing from baskets that are still recruiting the factors do not telescope; those
##  designs need simulation: R/simulate.R.)
## =============================================================================

## log E (or log E^f) at y responses out of n, for a vector y; side = grid (efficacy) or futility_side().
mixture_logE <- function(w0, side, y, n) {
  lr <- outer(y, side$l1) + outer(n - y, side$l0) + rep(log(w0), each = length(y))
  m <- row_max(lr); m + log(rowSums(exp(lr - m)))
}

## One basket with true response probability p, analysed after every `block` patients.
##  futility: "none", "evalue" (E^f >= 1/beta, needs w0f and fut) or "unreachable" (needs finite n_max).
##  Patients up to min(n_max, horizon) are followed; mass still open at the horizon is reported.
exact_basket <- function(p, w0, grid, alpha, block = 1, n_max = Inf, horizon = 3000,
                         futility = "none", w0f = NULL, fut = NULL, beta = NULL, report_at = c(30, 100, 300, 1000, 3000)) {
  thr <- log(1 / alpha); end <- min(n_max, horizon); step <- dbinom(0:block, block, p)
  prob <- 1; n <- 0                                   # prob[y + 1] = P(open, y responses, n patients)
  stop_n <- c(); rej <- c(); fail <- c(); by <- setNames(rep(NA_real_, length(report_at)), report_at); cum_rej <- 0
  while (n < end && sum(prob) > 1e-14) {
    b <- min(block, end - n); st <- if (b == block) step else dbinom(0:b, b, p)
    new <- numeric(length(prob) + b)
    for (j in 0:b) new[j + seq_along(prob)] <- new[j + seq_along(prob)] + prob * st[j + 1]
    prob <- new; n <- n + b
    live <- which(prob > 0); y <- live - 1
    r <- f <- 0
    if (length(live)) {
      hit <- mixture_logE(w0, grid, y, n) >= thr
      r <- sum(prob[live[hit]]); prob[live[hit]] <- 0
      if (futility == "evalue") {
        live <- live[!hit]; y <- live - 1
        if (length(live)) { hf <- mixture_logE(w0f, fut, y, n) >= log(1 / beta); f <- sum(prob[live[hf]]); prob[live[hf]] <- 0 }
      } else if (futility == "unreachable") {
        live <- live[!hit]; y <- live - 1
        if (length(live)) { hf <- mixture_logE(w0, grid, y, n) + (n_max - n) * log(1 / grid$p0) < thr; f <- sum(prob[live[hf]]); prob[live[hf]] <- 0 }
      }
    }
    stop_n <- c(stop_n, n); rej <- c(rej, r); fail <- c(fail, f); cum_rej <- cum_rej + r
    by[as.character(report_at[report_at == n])] <- cum_rej
  }
  open <- sum(prob); stopped <- rej + fail
  by[as.numeric(names(by)) > n] <- cum_rej             # nothing more can happen after the end
  q <- function(a) { cs <- cumsum(stopped); if (max(cs) < a) NA else stop_n[which(cs >= a)[1]] }
  list(summary = data.frame(p = p, reject = sum(rej), futility = sum(fail), complete = if (is.finite(n_max) && n >= n_max) open else 0,
                            open = if (is.finite(n_max) && n >= n_max) 0 else open,
                            EN = sum(stop_n * stopped) + n * open, N_median = q(0.5), N_q90 = q(0.9), last_n = n),
       reject_by = by)
}
