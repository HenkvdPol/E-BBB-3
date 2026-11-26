## =============================================================================
##  simulate.R -- the framework e-process in calendar time for a J x L basket trial
##
##  Basket c opens at entry[c] and recruits rate[c] patients per time unit up to n_max[c]
##  (n_max = Inf: open-ended, no maximum sample size). Global analyses at look_every, 2 look_every, ...
##  At every analysis tau:
##    1. each recruiting basket multiplies E_c (and, with e-value futility, E^f_c) by the block factor
##       of the patients that arrived since the previous analysis, with the priors frozen then;
##    2. decisions, in this order: reject if E_c >= 1/alpha; futility (below); complete at a finite
##       n_max; still open when the simulation horizon n_cap is reached (open-ended designs only);
##    3. the borrowing engine is fitted to ALL platform data at tau (stopped baskets included) and
##       gives the priors for the next block of every basket that still recruits.
##  So the prior of a block is a function of the data before the block only (predictable), and it
##  does not change while the block accrues (frozen). An empty block has factor 1.
##
##  design: layout (data.frame: basket, row, col), n_max (may be Inf), rate, entry (per basket or one
##          value), look_every, alpha, prior_form ("grid": restricted posterior; "beta": scaled Beta),
##          futility: "unreachable" (default) stop when 1/alpha cannot be reached by n_max even if all
##                                remaining patients respond (factor <= 1/p0 per patient); never
##                                stops an open-ended basket;
##                    "evalue"    anytime-valid futility: stop when E^f_c >= 1/beta, where E^f tests
##                                p >= p1 with the posterior restricted to [0, p1) (eprocess.R);
##                                needs p1 and beta; P(futility stop | p >= p1) <= beta;
##                    "none";
##          n_cap: simulation horizon per basket, required when n_max is infinite. It only stops the
##                 computation: baskets still open at n_cap are reported as open, not as decided.
##  Every futility rule and n_cap are optional stopping, so the type I error bound is unaffected.
## =============================================================================

planned_n <- function(design, t) pmin(design$n_sim, floor(design$rate * pmax(0, t - design$entry) + 1e-9))

complete_design <- function(design) {
  K <- nrow(design$layout)
  for (v in c("n_max", "rate", "entry")) design[[v]] <- rep_len(design[[v]], K)
  if (is.null(design$prior_form)) design$prior_form <- "grid"
  if (is.null(design$futility)) design$futility <- "unreachable"
  if (!design$futility %in% c("unreachable", "evalue", "none")) stop("futility must be 'unreachable', 'evalue' or 'none'")
  if (design$futility == "evalue" && (is.null(design$p1) || is.null(design$beta))) stop("e-value futility needs p1 and beta")
  if (any(is.infinite(design$n_max)) && is.null(design$n_cap)) stop("open-ended design (n_max = Inf): set the simulation horizon n_cap")
  design$n_sim <- pmin(design$n_max, if (is.null(design$n_cap)) Inf else design$n_cap)   # patients simulated per basket
  design
}

## One posterior call for the distinct data states among the trials (rows) that need one.
priors_and_summaries <- function(engine, y, n, need, seed) {
  rows <- which(rowSums(need) > 0)
  key <- do.call(paste, c(as.data.frame(cbind(y, n)[rows, , drop = FALSE]), sep = "."))
  u <- unique(key); iu <- match(key, u); first <- rows[match(u, key)]
  needu <- rowsum(need[rows, , drop = FALSE] * 1, iu) > 0
  po <- engine$posterior(y[first, , drop = FALSE], n[first, , drop = FALSE], needu, seed)
  state <- integer(nrow(y)); state[rows] <- iu
  list(po = po, state = state)
}

## logpost rows for the pairs (trial s, basket k) given by linear index `ix` into the S x K matrices.
pairs_logpost <- function(fit, ix, S) {
  s <- (ix - 1) %% S + 1; k <- (ix - 1) %/% S + 1; G <- dim(fit$po$logpost)[3]; i <- fit$state[s]
  matrix(fit$po$logpost[cbind(rep(i, G), rep(k, G), rep(seq_len(G), each = length(ix)))], length(ix), G)
}

simulate_trials <- function(p_true, design, engine, n_sim, seed = 1, X = NULL, trace = FALSE) {
  design <- complete_design(design); grid <- engine$grid
  K <- nrow(design$layout); S <- n_sim; nmx <- matrix(design$n_max, S, K, byrow = TRUE)
  ncap <- matrix(design$n_sim, S, K, byrow = TRUE)
  if (is.null(X)) {                                     # outcomes of patient i of basket k in trial s
    set.seed(seed)
    X <- array(runif(S * K * max(design$n_sim)) < rep(p_true, each = S), c(S, K, max(design$n_sim)))
  }
  CY <- array(0L, c(S, K, dim(X)[3] + 1))
  for (i in seq_len(dim(X)[3])) CY[, , i + 1] <- CY[, , i] + X[, , i]
  thr <- log(1 / design$alpha); lmax <- log(1 / grid$p0)
  efut <- design$futility == "evalue"
  if (efut) { fut <- futility_side(grid, design$p1); thr_f <- log(1 / design$beta) }
  times <- seq(design$look_every, max(design$entry + design$n_sim / design$rate) + design$look_every, by = design$look_every)
  n <- y <- matrix(0L, S, K); logE <- logEf <- matrix(0, S, K)
  status <- matrix(0L, S, K)          # 0 recruiting, 1 rejected, 2 stopped for futility, 3 complete at n_max, 4 open at n_cap
  res <- list(time = matrix(NA_real_, S, K), post_mean = matrix(NA_real_, S, K), p_above = matrix(NA_real_, S, K),
              ci_lo = matrix(NA_real_, S, K), ci_hi = matrix(NA_real_, S, K), pex = matrix(NA_real_, S, K),
              prior_mean = matrix(NA_real_, S, K), rhat = matrix(NA_real_, S, K))
  make_prior <- function(lp) { w <- restrict(lp, grid); if (design$prior_form == "beta") scaled_beta(w, grid) else w }
  ## priors for the first block: fitted to the empty platform (one state for all trials)
  fit <- priors_and_summaries(engine, n, n, matrix(TRUE, S, K), seed)
  first_priors <- make_prior(pairs_logpost(fit, (seq_len(K) - 1L) * S + 1L, S))  # K x support points
  W <- first_priors[rep(seq_len(K), each = S), , drop = FALSE]  # row s + (k - 1) S
  check_priors(W, grid)
  if (efut) {                                                   # futility priors on [0, p1), same rows
    Wf <- restrict_below(pairs_logpost(fit, (seq_len(K) - 1L) * S + 1L, S), fut)[rep(seq_len(K), each = S), , drop = FALSE]
    check_priors_below(Wf, fut)
  }
  tr <- list()
  for (a in seq_along(times)) {
    act <- status == 0
    n_new <- n; n_new[act] <- matrix(planned_n(design, times[a]), S, K, byrow = TRUE)[act]
    y_new <- matrix(CY[cbind(rep(seq_len(S), K), rep(seq_len(K), each = S), as.vector(n_new) + 1)], S, K)
    upd <- which(act & n_new > n)
    if (length(upd)) logE[upd] <- logE[upd] + log_factor(W[upd, , drop = FALSE], (y_new - y)[upd], (n_new - n)[upd], grid)
    if (efut && length(upd)) logEf[upd] <- logEf[upd] + log_factor(Wf[upd, , drop = FALSE], (y_new - y)[upd], (n_new - n)[upd], fut)
    n <- n_new; y <- y_new
    status[act & logE >= thr] <- 1L
    if (efut) status[act & status == 0 & logEf >= thr_f] <- 2L
    status[act & status == 0 & n >= nmx] <- 3L
    if (design$futility == "unreachable") status[act & status == 0 & logE + (nmx - n) * lmax < thr] <- 2L
    status[act & status == 0 & n >= ncap] <- 4L
    now <- act & status != 0
    res$time[now] <- times[a]
    res$prior_mean[now] <- (W %*% grid$support)[now]
    need <- status == 0 | now
    if (trace) tr[[a]] <- list(time = times[a], n = n, y = y, logE = logE, logEf = logEf, prior_mean = matrix(W %*% grid$support, S, K))
    if (!any(need)) break
    fit <- priors_and_summaries(engine, y, n, need, seed + a)
    ix <- which(status == 0)                                  # new frozen priors for the next block
    if (length(ix)) { lp <- pairs_logpost(fit, ix, S); W[ix, ] <- make_prior(lp); if (efut) Wf[ix, ] <- restrict_below(lp, fut) }
    ix <- which(now)                                          # summaries of the posterior at the decision
    if (length(ix)) {
      lp <- pairs_logpost(fit, ix, S); pr <- exp(lp - row_max(lp)); pr <- pr / rowSums(pr); cp <- t(apply(pr, 1, cumsum))
      res$post_mean[ix] <- pr %*% grid$p; res$p_above[ix] <- rowSums(pr[, grid$above, drop = FALSE])
      res$ci_lo[ix] <- grid$p[max.col(cp >= 0.05, "first")]; res$ci_hi[ix] <- grid$p[max.col(cp >= 0.95, "first")]
      s <- (ix - 1) %% S + 1; k <- (ix - 1) %/% S + 1
      res$pex[ix] <- fit$po$pex[cbind(fit$state[s], k)]
      if (!is.null(fit$po$rhat)) res$rhat[ix] <- fit$po$rhat[cbind(fit$state[s], k)]
    }
  }
  out <- data.frame(trial = rep(seq_len(S), K), basket = rep(design$layout$basket, each = S),
                    row = rep(design$layout$row, each = S), col = rep(design$layout$col, each = S),
                    p_true = rep(p_true, each = S), null = rep(p_true <= grid$p0, each = S),
                    reject = as.vector(status == 1), stopped = as.vector(status == 2), open = as.vector(status == 4),
                    n = as.vector(n), y = as.vector(y), logE = as.vector(logE), logEf = as.vector(logEf), lapply(res, as.vector))
  out$cover <- out$ci_lo <= out$p_true & out$p_true <= out$ci_hi
  if (trace) attr(out, "trace") <- tr
  out
}
