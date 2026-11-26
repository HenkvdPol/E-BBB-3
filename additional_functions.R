## =============================================================================
##  additional_functions.R  --  all functions (nothing is run here)
##
##  1. EXNEX model: posterior and borrowing prior of each cohort
##  2. pi-rules (EXNEX weights)
##  3. Borrowing prior -> alternative theta1
##  4. Design-optimal e-value designs (uses design_optimal_evalues.R)
##  5. Trial simulation in calendar time: e-value designs (incl. acceleration and
##     extension beyond n_max) and Bayesian designs
##  6. Summaries, references (Simon, binomial) and figures
##  7. One run: complete the settings, run everything, overview rows
## =============================================================================

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))
row_max <- function(x) if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMaxs(x) else apply(x, 1, max)
`%||%` <- function(a, b) if (is.null(a)) b else a

## -----------------------------------------------------------------------------
## 1. EXNEX model
## -----------------------------------------------------------------------------
##  y_k ~ Bin(n_k, p_k), theta_k = logit(p_k)
##  EX  (prob pi_k):  theta_k ~ N(mu, tau^2),  mu ~ N(m, s^2),  tau ~ HalfNormal(s_tau)
##  NEX (1 - pi_k) :  theta_k ~ N(logit(nex_p), 1/nex_p + 1/(1 - nex_p))
##  Given (mu, tau) the cohorts are independent, so the posterior is a sum over
##  a (mu, tau) grid of pre-computed binomial-normal integrals (no MCMC).

## Pre-computes the integrals for every (n, y) with n <= max(n_max).
build_engine <- function(settings) {
  g <- settings$grid; pr <- settings$prior
  theta <- seq(g$theta_min, g$theta_max, by = g$theta_step)
  above <- function(t) pmin(pmax((theta + g$theta_step / 2 - t) / g$theta_step, 0), 1)
  p <- expit(theta)
  hyp <- expand.grid(mu = seq(g$mu_min, g$mu_max, by = g$mu_step),
                     tau = seq(g$tau_step / 2, g$tau_max, by = g$tau_step))
  W <- dnorm(sweep(outer(theta, hyp$mu, "-"), 2, hyp$tau, "/"))        # EX: theta | mu, tau
  W <- sweep(W, 2, colSums(W), "/")
  V <- dnorm(theta, logit(pr$nex_p), sqrt(1 / pr$nex_p + 1 / (1 - pr$nex_p)))
  V <- V / sum(V)                                                       # NEX
  Nmax <- max(settings$n_max + settings$n_ext)               # cohorts can be extended beyond n_max
  nn <- rep(0:Nmax, times = 1:(Nmax + 1)); yy <- sequence(1:(Nmax + 1)) - 1
  Lik <- t(vapply(seq_along(nn), function(i) dbinom(yy[i], nn[i], p), numeric(length(theta))))
  times <- function(v) Lik * rep(v, each = nrow(Lik))
  q_grid <- seq(0.005, 0.995, by = 0.005)
  list(
    Nmax = Nmax, H = nrow(hyp),
    logprior = dnorm(hyp$mu, pr$mu_mean, pr$mu_sd, log = TRUE) + dnorm(hyp$tau, 0, pr$tau_scale, log = TRUE),
    LE = Lik %*% W, TE = times(above(logit(settings$p0))) %*% W,            # EX: likelihood, P(p > p0)
    FE = times(above(logit(settings$futility$p_threshold))) %*% W,          # EX: P(p > p_fut)
    LN = as.vector(Lik %*% V), TN = as.vector(times(above(logit(settings$p0))) %*% V),
    FN = as.vector(times(above(logit(settings$futility$p_threshold))) %*% V),
    ME = times(above(logit(settings$p0)) * p) %*% W,                        # EX: E(p 1{p > p0}), E(p^2 1{p > p0})
    SE = times(above(logit(settings$p0)) * p^2) %*% W,
    MN = as.vector(times(above(logit(settings$p0)) * p) %*% V), SN = as.vector(times(above(logit(settings$p0)) * p^2) %*% V),
    mEX = as.vector(p %*% W),                                               # E(p | mu, tau)
    CDF = sapply(q_grid, function(q) colSums(W * (1 - above(logit(q))))),  # P(p <= q | mu, tau)
    q_grid = q_grid)
}

## Posterior for many data sets at once (rows = trials, columns = cohorts):
##   pp_eff = P(p_k > p0 | y), pp_fut = P(p_k > p_fut | y),
##   loo    = mean or `level`-quantile of the EX predictive of p_k given the OTHER cohorts,
##   moments = TRUE: rmean, rsec = E(p_k | y, p_k > p0), E(p_k^2 | y, p_k > p0): the posterior
##            of p_k given ALL data restricted to (p0, 1] (the prior of framework.tex).
posterior_exnex <- function(eng, y, n, pi, summary = "mean", level = 0.25, moments = FALSE) {
  S <- nrow(y); K <- ncol(y); H <- eng$H
  r <- n * (n + 1) / 2 + y + 1                                          # row of (n, y)
  out <- list(pp_eff = eng$TN[r] / eng$LN[r], pp_fut = eng$FN[r] / eng$LN[r], loo = rep(NA_real_, S * K))
  dim(out$pp_eff) <- dim(out$pp_fut) <- dim(out$loo) <- c(S, K)
  if (moments) { out$rmean <- matrix(eng$MN[r] / eng$TN[r], S, K); out$rsec <- matrix(eng$SN[r] / eng$TN[r], S, K) }
  summary_of <- function(w) {                                           # w: weights over (mu, tau)
    if (summary == "mean") return(as.vector(w %*% eng$mEX))
    eng$q_grid[max.col(w %*% eng$CDF >= level, ties.method = "first")]
  }
  todo <- which(rowSums(pi > 0) > 0)                                    # rows with borrowing
  chunk <- max(1, floor(3e6 / (K * H)))
  for (st in seq_len(ceiling(length(todo) / chunk))) {
    ids <- todo[((st - 1) * chunk + 1):min(st * chunk, length(todo))]
    s <- length(ids)
    D <- lapply(seq_len(K), function(k) eng$LE[r[ids, k], , drop = FALSE] * pi[ids, k] +
                                        (1 - pi[ids, k]) * eng$LN[r[ids, k]])
    logw_all <- matrix(eng$logprior, s, H, byrow = TRUE) + Reduce(`+`, lapply(D, log))
    w <- exp(logw_all - row_max(logw_all)); w <- w / rowSums(w)
    for (k in seq_len(K)) {
      rk <- r[ids, k]; pk <- pi[ids, k]
      wD <- w / D[[k]]; swD <- rowSums(wD)
      out$pp_eff[ids, k] <- pk * rowSums(wD * eng$TE[rk, , drop = FALSE]) + (1 - pk) * eng$TN[rk] * swD
      out$pp_fut[ids, k] <- pk * rowSums(wD * eng$FE[rk, , drop = FALSE]) + (1 - pk) * eng$FN[rk] * swD
      if (moments) {
        out$rmean[ids, k] <- (pk * rowSums(wD * eng$ME[rk, , drop = FALSE]) + (1 - pk) * eng$MN[rk] * swD) / out$pp_eff[ids, k]
        out$rsec[ids, k]  <- (pk * rowSums(wD * eng$SE[rk, , drop = FALSE]) + (1 - pk) * eng$SN[rk] * swD) / out$pp_eff[ids, k]
      }
      lw <- logw_all - log(D[[k]])                                       # leave cohort k out
      wk <- exp(lw - row_max(lw)); out$loo[ids, k] <- summary_of(wk / rowSums(wk))
    }
  }
  out
}

## Same, with a cache: data states repeat a lot across trials and analyses.
## States are sorted first (the posterior is symmetric in the cohorts).
new_cache <- function() new.env(hash = TRUE)
posterior_cached <- function(eng, design, y, n, cache) {
  S <- nrow(y); K <- ncol(y); mom <- design$type == "mixture"
  code <- n * (eng$Nmax + 1) + y
  ord <- matrix(((order(rep(seq_len(S), K), as.vector(code)) - 1) %/% S) + 1, S, K, byrow = TRUE)
  sc <- matrix(code[cbind(rep(seq_len(S), K), as.vector(ord))], S, K)
  key <- do.call(paste, c(as.data.frame(sc), sep = "."))
  miss <- unique(key[!vapply(mget(key, envir = cache, ifnotfound = list(NULL)), length, 1L)])
  if (length(miss)) {
    m <- sc[match(miss, key), , drop = FALSE]
    nm <- m %/% (eng$Nmax + 1); ym <- m %% (eng$Nmax + 1)
    pi <- compute_pi(design, ym, nm)
    po <- posterior_exnex(eng, ym, nm, pi, design$summary %||% "mean", design$quantile %||% 0.25, mom)
    val <- cbind(po$pp_eff, po$pp_fut, pi, po$loo, if (mom) cbind(po$rmean, po$rsec))
    list2env(setNames(split(val, row(val)), miss), envir = cache)
  }
  v <- matrix(unlist(mget(key, envir = cache), use.names = FALSE), S, (4 + 2 * mom) * K, byrow = TRUE)
  back <- function(m) { o <- matrix(NA_real_, S, K); o[cbind(rep(seq_len(S), K), as.vector(ord))] <- m; o }
  out <- list(pp_eff = back(v[, 1:K]), pp_fut = back(v[, K + 1:K]), pi = back(v[, 2 * K + 1:K]), loo = back(v[, 3 * K + 1:K]))
  if (mom) { out$rmean <- back(v[, 4 * K + 1:K]); out$rsec <- back(v[, 5 * K + 1:K]) }
  out
}

## -----------------------------------------------------------------------------
## 2. pi-rules (overview in Setting.R)
## -----------------------------------------------------------------------------
hellinger_beta <- function(a1, b1, a2, b2) {
  sqrt(pmax(0, 1 - exp(lbeta((a1 + a2) / 2, (b1 + b2) / 2) - (lbeta(a1, b1) + lbeta(a2, b2)) / 2)))
}

## mEXNEX weights. type "mean" (mEXNEX) or "max" (mEXNEX_min); elig = possible partners.
pi_hellinger <- function(y, n, type = "mean", c_e = NULL, elig = n > 0) {
  S <- nrow(y); K <- ncol(y); a <- y + 1; b <- n - y + 1
  keep <- elig
  if (!is.null(c_e)) {                       # 'extreme' cohort: rate differs > c_e from all others
    est <- y / pmax(n, 1)
    mind <- sapply(seq_len(K), function(k) {
      d <- abs(est[, k] - est[, -k, drop = FALSE]); d[!elig[, -k, drop = FALSE]] <- Inf
      apply(d, 1, min)
    })
    keep <- elig & !(matrix(mind, S, K) > c_e + 1e-5)
  }
  pi <- sapply(seq_len(K), function(k) {
    s <- sapply(setdiff(seq_len(K), k), function(j) 1 - hellinger_beta(a[, k], b[, k], a[, j], b[, j]))
    s <- matrix(s, S); m <- keep[, -k, drop = FALSE]
    if (type == "mean") ifelse(rowSums(m) > 0, rowSums(s * m) / pmax(rowSums(m), 1), 0)
    else { s[!m] <- 0; apply(s, 1, max) }
  })
  pi <- matrix(pi, S, K); pi[!keep] <- 0
  pi
}

compute_pi <- function(design, y, n) {
  S <- nrow(y); K <- ncol(y)
  switch(design$pi_rule,
    none  = matrix(0, S, K),
    full  = matrix(1, S, K),
    fixed = matrix(design$pi_fixed, S, K),
    hellinger_mean = pi_hellinger(y, n, "mean", design$c_e),
    hellinger_max  = pi_hellinger(y, n, "max",  design$c_e),
    delayed = {
      elig <- n >= design$n_burnin
      pi <- pi_hellinger(y, n, sub("hellinger_", "", design$base_rule), design$c_e, elig = elig)
      pi[!elig] <- 0
      pi
    })
}

## -----------------------------------------------------------------------------
## 3. Borrowing prior -> theta1
## -----------------------------------------------------------------------------
theta1_grid <- function(settings) {
  b <- settings$borrow
  round(seq(b$theta1_min, b$theta1_max, by = b$theta1_step), 6)
}
## theta1_k = pi_k * s_k + (1 - pi_k) * theta1_nex, clipped and rounded; returns the
## index in theta1_grid(). theta1_nex while the other cohorts have no patients
## (s_k would only be the EX prior). n: patients per cohort.
theta1_index <- function(po, n, settings) {
  b <- settings$borrow
  th <- po$pi * po$loo + (1 - po$pi) * b$theta1_nex
  th[po$pi <= 0 | is.na(th) | rowSums(n) - n == 0] <- b$theta1_nex
  th <- pmin(pmax(th, b$theta1_min), b$theta1_max)
  matrix(as.integer(round((th - b$theta1_min) / b$theta1_step)) + 1L, nrow(po$pi))
}

## -----------------------------------------------------------------------------
## 4. Design-optimal e-value designs (design_optimal_evalues.R)
## -----------------------------------------------------------------------------
## One design for alternative theta1, horizon n, analysis times `at`.
make_evalue_design <- function(type, theta1, n, at, settings) {
  p0 <- settings$p0; a <- settings$alpha
  switch(type,
    grow = design_grow(n, p0, theta1, a, analysis_times = at),
    pmax = design_pmax(n, p0, theta1, a, analysis_times = at),
    essmin = design_essmin(n, p0, theta1, a, analysis_times = at),
    ebased = {                     # power target capped at what is attainable (ebased_capped())
      pm <- design_pmax(n, p0, theta1, a, analysis_times = at)
      target <- min(1 - settings$beta, settings$evalue$power_cap * evaluate_design(pm, theta1)$reject)
      suppressWarnings(design_ebased(n, p0, theta1, a, beta = 1 - target, analysis_times = at))
    })
}

## Name of a bank (also its file name in the cache): everything the designs depend on.
bank_key <- function(type, n, at, borrow, settings) {
  th <- if (borrow) with(settings$borrow, sprintf("%g-%g-%g", theta1_min, theta1_max, theta1_step))
        else sprintf("%g", settings$borrow$theta1_nex)
  sprintf("%s_p0=%g_a=%g%s_th=%s_n=%d_at=%s", type, settings$p0, settings$alpha,
          if (type == "ebased") sprintf("_b=%g_cap=%g", settings$beta, settings$evalue$power_cap) else "",
          th, n, if (identical(as.integer(at), seq_len(n))) "all" else paste(at, collapse = ","))
}

## A bank holds the betting policies of one design type for every theta1 on the
## grid (borrowing designs) or for theta1_nex only:
##   bet_idx[i, t, j] = index in b_val of the bet for patient t at grid e-value
##   m_grid[i] when betting on theta1 number j;  stop[i, t, j] = futility stop.
build_bank <- function(type, n, at, borrow, settings) {
  th <- if (borrow) theta1_grid(settings) else settings$borrow$theta1_nex
  ds <- parallel::mclapply(th, make_evalue_design, type = type, n = n, at = at,
                           settings = settings, mc.cores = settings$sim$n_cores)
  b_val <- sort(unique(round(na.omit(unlist(lapply(ds, function(d) if (type == "grow") d$bet else d$b_grid))), 12)))
  stopifnot(length(b_val) < 255)                    # bets are stored as raw bytes
  nm <- length(ds[[1]]$m_grid)
  bet_idx <- array(NA_integer_, c(nm, n, length(th))); stp <- array(FALSE, c(nm, n, length(th)))
  for (j in seq_along(th)) {
    bet_idx[, , j] <- match(round(ds[[j]]$bet, 12), round(b_val, 12))
    stp[, , j] <- ds[[j]]$stop
  }
  list(theta1 = th, m_grid = ds[[1]]$m_grid, b_val = b_val, bet_idx = bet_idx, stop = stp)
}

## Compact copy of a bank for the simulation: bet indices and stop flags as raw
## bytes (index 0 = no bet), 4x less memory than integers.
compact_bank <- function(bk) {
  idx <- bk$bet_idx; idx[is.na(idx)] <- 0L
  list(theta1 = bk$theta1, m_grid = bk$m_grid, b_val0 = c(0, bk$b_val),
       bet_idx = array(as.raw(idx), dim(idx)), stop = array(as.raw(bk$stop), dim(bk$stop)))
}

## All banks needed by one run. Each bank is cached as its own file in
## settings$sim$bank_dir (file name = bank_key), so a run only loads what it
## needs and missing banks are built once. Designs that can be extended also need
## the extension design: P-max for n_ext patients, analysed after every patient.
needed_banks <- function(settings) {
  out <- list()
  for (su in settings$setups) for (dn in su$designs) {
    d <- settings$designs[[dn]]; if (d$type != "evalue") next
    for (block in settings$schedules) {
      at <- cohort_times(block, settings$n_max, su$entry)
      for (k in seq_len(settings$K)) {
        key <- bank_key(d$evalue, settings$n_max[k], at[[k]], d$borrow, settings)
        out[[key]] <- list(type = d$evalue, n = settings$n_max[k], at = at[[k]], borrow = d$borrow)
      }
    }
    if (isTRUE(d$extend)) for (block in settings$schedules) {
      at <- cohort_times(block, settings$n_max, su$entry)
      for (k in seq_len(settings$K)) {
        r <- settings$n_max[k] - ext_point(at[[k]], settings$n_max[k], settings) + settings$n_ext[k]
        out[[bank_key("pmax", r, seq_len(r), d$borrow, settings)]] <- list(type = "pmax", n = r, at = seq_len(r), borrow = d$borrow)
      }
    }
  }
  out
}
get_banks <- function(settings) {
  dir <- settings$sim$bank_dir; dir.create(dir, showWarnings = FALSE)
  need <- needed_banks(settings)
  for (key in names(need)) {
    f <- file.path(dir, paste0(key, ".rds"))
    if (!file.exists(f)) {
      message("  building e-value designs: ", key)
      x <- need[[key]]
      saveRDS(compact_bank(build_bank(x$type, x$n, x$at, x$borrow, settings)), f)
    }
  }
  lapply(setNames(names(need), names(need)), function(key) readRDS(file.path(dir, paste0(key, ".rds"))))
}

## The banks of one design under one schedule, one per cohort.
banks_for <- function(design, block, setup, settings, banks) {
  at <- cohort_times(block, settings$n_max, setup$entry)
  lapply(seq_len(settings$K), function(k) banks[[bank_key(design$evalue, settings$n_max[k], at[[k]], design$borrow, settings)]])
}

## Conditional power of a bank's designs: P(grid e-value reaches 1/alpha by the end
## of the bank's horizon | grid e-value m_grid[i] after `from` patients, p = theta1[j]),
## for all i and j: backward recursion with the bank's own bets and futility stops.
cond_power <- function(bk, p0, from = 0) {
  g <- bk$m_grid; top <- length(g); n <- dim(bk$bet_idx)[2]
  sapply(seq_along(bk$theta1), function(j) {
    th <- bk$theta1[j]; V <- as.numeric(seq_along(g) == top)
    for (t in rev(seq_len(n))[seq_len(max(n - from, 0))]) {      # t = n, ..., from + 1
      b <- bk$b_val0[as.integer(bk$bet_idx[, t, j]) + 1L]
      up <- findInterval(g * (1 + b * (1 / p0 - 1)), g); dn <- findInterval(g * (1 - b), g)
      up[top] <- dn[top] <- top
      V <- th * V[up] + (1 - th) * V[dn]
      V[as.integer(bk$stop[, t, j]) == 1L] <- 0                  # the design stops before patient t
    }
    V
  })
}

## Decision point for extension: the last planned analysis of the cohort at or
## before ext_at * n_max (interim of a two-stage design, 3/4 of n_max if sequential).
ext_point <- function(at, n_max, settings) {
  a <- at[at <= settings$flex$ext_at * n_max]; if (length(a)) max(a) else min(at)
}

## The extension designs of one design under one schedule, one per cohort: P-max
## for the patients after the decision point td up to n_max + n_ext, analysed after
## every patient, so it can be entered later (accelerated cohorts) and is optimal from
## every e-value (the dynamic program is solved for all states).
ext_banks_for <- function(design, block, setup, settings, banks) {
  at <- cohort_times(block, settings$n_max, setup$entry)
  lapply(seq_len(settings$K), function(k) {
    td <- ext_point(at[[k]], settings$n_max[k], settings); r <- settings$n_max[k] - td + settings$n_ext[k]
    c(banks[[bank_key("pmax", r, seq_len(r), design$borrow, settings)]], list(td = td))
  })
}

## -----------------------------------------------------------------------------
## 5. Trial simulation
## -----------------------------------------------------------------------------
## Planned cohort sizes at each analysis (rows) for calendar time in patients per cohort.
look_plan <- function(block, n_max, entry) {
  times <- sort(unique(c(seq(block, max(entry + n_max), by = block), entry + n_max)))
  sapply(seq_along(n_max), function(k) pmin(pmax(times - entry[k], 0), n_max[k]))
}
cohort_times <- function(block, n_max, entry) {
  lp <- matrix(look_plan(block, n_max, entry), ncol = length(n_max))
  lapply(seq_along(n_max), function(k) sort(unique(lp[lp[, k] > 0, k])))
}

## Simulated data: C[s, k, i + 1] = responders among the first i patients of cohort k.
## U holds the uniforms, so all designs and scenarios use the same random numbers.
make_data <- function(U, p) {
  C <- array(0L, dim(U) + c(0, 0, 1))
  for (i in seq_len(dim(U)[3])) C[, , i + 1] <- C[, , i] + (U[, , i] < rep(p, each = dim(U)[1]))
  C
}

## e-value design in calendar time. One time unit = one new patient per cohort at
## normal speed (cohort k starts at entry[k]). Analyses at the planned times (every
## `block` time units, and when a cohort reaches n_max at normal speed), after that
## every `block` time units; a cohort that reaches its maximum is analysed at once.
## Per cohort the real e-value M (it decides) and the grid (shadow) e-value (it looks
## up the bets), as in simulate_design() of design_optimal_evalues.R.
## Flexible inclusion (design$accelerate, design$extend; settings$flex):
##  accelerate: while its e-value is >= m_speed, a cohort includes speed_up patients
##              per time unit (decided at every analysis);
##  extend:     at its decision point td (first analysis with >= td patients) a cohort
##              without a decision is extended by n_ext patients if the planned design
##              will probably fall short (conditional power < 1 - beta) and the design
##              re-optimised from the current e-value for the remaining patients + n_ext
##              (P-max) is promising (conditional power >= cp_min). An extended cohort
##              bets with the re-optimised design from then on.
## status: 1 = H0 rejected, 2 = futility stop, 3 = maximum reached without rejection.
## time: calendar time from the cohort's entry to its decision.
run_trials_evalue <- function(design, C, block, entry, settings, eng, cache, banks, ext = NULL) {
  if (identical(design$redesign, "entry")) {
    stop("redesign = 'entry' is temporarily disabled: cohort-entry-time updates are not yet represented as separate events. See docs/known-limitations.md.")
  }
  S <- dim(C)[1]; K <- dim(C)[2]; p0 <- settings$p0; thr <- 1 / settings$alpha; fl <- settings$flex
  acc <- isTRUE(design$accelerate); extend <- isTRUE(design$extend)
  NM  <- matrix(settings$n_max, S, K, byrow = TRUE)
  NE  <- if (extend) matrix(settings$n_ext, S, K, byrow = TRUE) else 0L * NM
  cap <- NM; HZ <- NM + NE                    # current maximum; horizon of the hopeless zone
  ENT <- matrix(entry, S, K, byrow = TRUE)
  g <- banks[[1]]$m_grid; top <- length(g)
  if (extend) { stopifnot(identical(ext[[1]]$m_grid, g)); TD <- vapply(ext, `[[`, 0, "td")
                cp_memo <- new.env(); cp_at <- function(k, which, t) {       # conditional power, cached
                  key <- paste(k, which, t)
                  if (is.null(cp_memo[[key]])) cp_memo[[key]] <- if (which == "plan") cond_power(banks[[k]], p0, t)
                                                                else cond_power(ext[[k]], p0, t - TD[k])
                  cp_memo[[key]] } }
  status <- ncur <- ycur <- matrix(0L, S, K); speed <- matrix(1L, S, K)
  M <- matrix(1, S, K); idx <- matrix(which(g == 1), S, K); hit <- extended <- decided_ext <- matrix(FALSE, S, K)
  jth <- matrix(1L, S, K); time <- matrix(NA_real_, S, K); th_sum <- th_cnt <- numeric(K)
  new_theta1 <- function(rows) {
    n <- ncur[rows, , drop = FALSE]
    theta1_index(posterior_cached(eng, design, ycur[rows, , drop = FALSE], n, cache), n, settings)
  }
  if (design$borrow) jth[] <- new_theta1(seq_len(S))
  look <- function(bank, field, cells, t) as.integer(bank[[field]][cbind(idx[cells], t, jth[cells])])

  planned <- sort(unique(c(seq(block, max(entry + settings$n_max), by = block), entry + settings$n_max)))
  tau_prev <- 0; l <- 0L
  while (any(status == 0L)) {
    l <- l + 1L; tau <- if (l <= length(planned)) planned[l] else tau_prev + block
    start <- pmax(tau_prev, ENT); act <- status == 0L
    target <- ncur; target[act] <- as.integer(pmin(cap[act], ncur[act] + speed[act] * pmax(tau - start[act], 0)))
    grow <- target > ncur
    tau_prev <- tau
    if (!any(grow)) next
    tnow <- ifelse(target >= cap, start + (cap - ncur) / speed, tau)   # time of this analysis per cohort
    for (k in seq_len(K)) {                                   # patients one by one
      sk <- which(grow[, k]); if (!length(sk)) next
      bk <- banks[[k]]; nk <- ncur[sk, k]; tk <- target[sk, k]
      for (st in seq_len(max(tk - nk))) {
        r <- sk[nk + st <= tk]; t <- ncur[r, k] + st; cells <- cbind(r, k)
        y <- C[cbind(r, k, t + 1L)] - C[cbind(r, k, t)]
        main <- !extended[cells]; b <- numeric(length(r))
        if (any(main))  b[main]  <- bk$b_val0[look(bk, "bet_idx", cells[main, , drop = FALSE], t[main]) + 1L]
        if (any(!main)) b[!main] <- ext[[k]]$b_val0[look(ext[[k]], "bet_idx", cells[!main, , drop = FALSE],
                                                         t[!main] - TD[k]) + 1L]
        b[hit[cells]] <- 0
        f <- ifelse(y == 1L, 1 + b * (1 / p0 - 1), 1 - b)
        M[cells] <- M[cells] * f
        i <- idx[cells]; mv <- i != top; i[mv] <- findInterval(g[i[mv]] * f[mv], g); idx[cells] <- i
        hit[cells] <- hit[cells] | M[cells] >= thr
        th_sum[k] <- th_sum[k] + sum(bk$theta1[jth[cells]]); th_cnt[k] <- th_cnt[k] + length(r)
      }
      ncur[sk, k] <- tk; ycur[sk, k] <- C[cbind(sk, k, tk + 1L)]
    }
    ## analysis: efficacy -> new theta1 (borrowing) -> extension -> futility -> maximum
    decide <- function(cells, code) { status[cells] <<- code; time[cells] <<- tnow[cells] - ENT[cells] }
    decide(which(grow & hit & status == 0L), 1L)
    if (design$borrow && any(status == 0L)) {
      act <- which(rowSums(status == 0L) > 0)
      upd <- status[act, , drop = FALSE] == 0L
      jn <- new_theta1(act); jo <- jth[act, , drop = FALSE]; jo[upd] <- jn[upd]; jth[act, ] <- jo
    }
    if (extend) {                             # extension decision at the decision point
      cand <- which(grow & status == 0L & !decided_ext & ncur >= matrix(TD, S, K, byrow = TRUE), arr.ind = TRUE)
      if (nrow(cand)) {
        cp_plan <- cp_ext <- numeric(nrow(cand))
        for (i in seq_len(nrow(cand))) {
          k <- cand[i, 2]; ci <- cand[i, , drop = FALSE]; tt <- ncur[ci]
          if (tt < settings$n_max[k]) cp_plan[i] <- cp_at(k, "plan", tt)[idx[ci], jth[ci]]
          cp_ext[i] <- cp_at(k, "ext", tt)[idx[ci], jth[ci]]
        }
        go <- cp_plan < 1 - settings$beta & cp_ext >= fl$cp_min
        decided_ext[cand] <- TRUE
        extended[cand[go, , drop = FALSE]] <- TRUE; cap[cand[go, , drop = FALSE]] <- (NM + NE)[cand[go, , drop = FALSE]]
        HZ[cand[!go, , drop = FALSE]] <- NM[cand[!go, , drop = FALSE]]
      }
    }
    cand <- which(grow & status == 0L & ncur < cap, arr.ind = TRUE)
    if (nrow(cand)) {
      stp <- logical(nrow(cand))
      for (k in unique(cand[, 2])) {
        w <- which(cand[, 2] == k); cw <- cand[w, , drop = FALSE]; tt <- ncur[cw]; main <- !extended[cw]
        if (any(main))  stp[w[main]]  <- look(banks[[k]], "stop", cw[main, , drop = FALSE], tt[main] + 1L) == 1L
        if (any(!main)) stp[w[!main]] <- look(ext[[k]], "stop", cw[!main, , drop = FALSE], tt[!main] - TD[k] + 1L) == 1L
      }
      hz <- settings$evalue$futility_hz & g[idx[cand]] < p0^(HZ[cand] - ncur[cand]) / settings$alpha
      decide(cand[stp | hz, , drop = FALSE], 2L)
    }
    decide(which(grow & status == 0L & ncur >= cap), 3L)
    if (acc) { a <- status == 0L; speed[a] <- ifelse(M[a] >= fl$m_speed, as.integer(fl$speed_up), 1L) }
  }
  list(status = status, n = ncur, e_value = M, theta1 = th_sum / pmax(th_cnt, 1),
       time = time, extended = extended, entry = entry)
}

## Bayesian design: reject if P(p_k > p0 | data) > Delta; futility if
## P(p_k > p_threshold | data) < gamma (from min_n patients on).
run_trials_bayes <- function(design, C, block, entry, settings, eng, cache, Delta) {
  S <- dim(C)[1]; K <- dim(C)[2]; fu <- settings$futility
  plan <- matrix(look_plan(block, settings$n_max, entry), ncol = K)
  times <- sort(unique(c(seq(block, max(entry + settings$n_max), by = block), entry + settings$n_max)))
  NM <- matrix(settings$n_max, S, K, byrow = TRUE)
  status <- ncur <- matrix(0L, S, K); time <- matrix(NA_real_, S, K)
  Dm <- matrix(rep_len(Delta, K), S, K, byrow = TRUE)
  for (l in seq_len(nrow(plan))) {
    target <- matrix(plan[l, ], S, K, byrow = TRUE)
    grow <- status == 0L & target > ncur
    if (!any(grow)) next
    ncur[grow] <- target[grow]
    y <- matrix(C[cbind(rep(seq_len(S), K), rep(seq_len(K), each = S), as.vector(ncur) + 1L)], S, K)
    rows <- which(rowSums(grow) > 0)
    po <- posterior_cached(eng, design, y[rows, , drop = FALSE], ncur[rows, , drop = FALSE], cache)
    g <- grow[rows, , drop = FALSE]; nr <- ncur[rows, , drop = FALSE]; fin <- nr >= NM[rows, , drop = FALSE]
    eff <- g & po$pp_eff > Dm[rows, , drop = FALSE]
    fut <- g & !eff & !fin & fu$on & nr >= fu$min_n & po$pp_fut < fu$gamma
    st <- status[rows, , drop = FALSE]; st[eff] <- 1L; st[fut] <- 2L; st[g & !eff & !fut & fin] <- 3L
    status[rows, ] <- st
    new <- status != 0L & is.na(time); time[new] <- times[l] - matrix(entry, S, K, byrow = TRUE)[new]
  }
  list(status = status, n = ncur, time = time, entry = entry)
}

## The e-value of framework.tex. At every analysis the prior of cohort c for its next
## block is frozen: the posterior of p_c given ALL platform data so far (EXNEX, own data
## included), restricted to (p0, 1]. It is represented by the distribution on (p0, 1]
## with the same mean and variance (p = p0 + (1 - p0) u, u ~ Beta(a, b)); like the
## Monte Carlo average of framework.tex this replaces the exact restricted posterior by
## another frozen prior on (p0, 1], which keeps the type I error guarantee. The block
## factor is the likelihood ratio of the new block mixed over this prior (closed form);
## E = product of the block factors; H0 rejected at the first analysis with E >= 1/alpha.
## Planned analyses at normal speed; stop when even all responses cannot reach 1/alpha.
mixture_prior <- function(m, s2, p0) {           # (a, b) of the scaled Beta with mean m and 2nd moment s2
  u <- (m - p0) / (1 - p0); uv <- (s2 - m^2) / (1 - p0)^2
  bad <- !is.finite(u) | !is.finite(uv); u[bad] <- 1e-6; uv[bad] <- 1e-12
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6); uv <- pmin(pmax(uv, 1e-12), u * (1 - u) * (1 - 1e-6))
  sz <- u * (1 - u) / uv - 1
  list(a = u * sz, b = (1 - u) * sz)
}
## Block factor: integral of (p/p0)^y ((1-p)/(1-p0))^(n-y) over the scaled Beta(a, b) prior
## = sum_j choose(y, j) ((1-p0)/p0)^j B(a + j, b + n - y) / B(a, b).
mixture_factor <- function(a, b, y, n, p0) {
  J <- max(c(y, 0)); r <- log((1 - p0) / p0)
  lt <- matrix(vapply(0:J, function(j) ifelse(j <= y, lchoose(pmax(y, j), j) + j * r + lbeta(a + j, b + n - y) - lbeta(a, b), -Inf),
                      numeric(length(y))), length(y))
  mx <- row_max(lt); exp(mx + log(rowSums(exp(lt - mx))))
}
## Predictive probability that a block of n_ext patients brings E to 1/alpha, under the
## frozen prior: P(Y = y) = dbinom(y, n_ext, p0) * factor(y) (the factor is the ratio of
## the prior-predictive to the null probability of the block).
extension_power <- function(a, b, E, n_ext, p0, thr) {
  pw <- numeric(length(a))
  for (y in 0:max(n_ext)) {
    ok <- y <= n_ext; if (!any(ok)) next
    f <- mixture_factor(a[ok], b[ok], rep(y, sum(ok)), n_ext[ok], p0)
    pw[ok] <- pw[ok] + dbinom(y, n_ext[ok], p0) * f * (E[ok] * f >= thr)
  }
  pw
}
## The framework e-value in calendar time (same time line, analyses and flexible
## inclusion as run_trials_evalue()):
##  accelerate: speed_up patients per time unit while E >= m_speed;
##  extend:     the e-value has no horizon, so a cohort without a decision at n_max is
##              extended by n_ext patients when the predictive probability that the
##              extension brings E to 1/alpha is >= cp_min (extension_power()).
run_trials_mixture <- function(design, C, block, entry, settings, eng, cache) {
  S <- dim(C)[1]; K <- dim(C)[2]; p0 <- settings$p0; thr <- 1 / settings$alpha; fl <- settings$flex
  acc <- isTRUE(design$accelerate); extend <- isTRUE(design$extend)
  NM  <- matrix(settings$n_max, S, K, byrow = TRUE)
  NE  <- if (extend) matrix(settings$n_ext, S, K, byrow = TRUE) else 0L * NM
  cap <- NM; HZ <- NM + NE; ENT <- matrix(entry, S, K, byrow = TRUE)
  status <- ncur <- ycur <- matrix(0L, S, K); speed <- matrix(1L, S, K)
  E <- matrix(1, S, K); time <- matrix(NA_real_, S, K); extended <- decided_ext <- matrix(FALSE, S, K)
  q_sum <- q_cnt <- numeric(K)
  prior <- function(rows) {                     # frozen prior for the next block, from the history now
    po <- posterior_cached(eng, design, ycur[rows, , drop = FALSE], ncur[rows, , drop = FALSE], cache)
    mixture_prior(po$rmean, po$rsec, p0)
  }
  pr <- prior(seq_len(S)); A <- pr$a; B <- pr$b
  planned <- sort(unique(c(seq(block, max(entry + settings$n_max), by = block), entry + settings$n_max)))
  tau_prev <- 0; l <- 0L
  while (any(status == 0L)) {
    l <- l + 1L; tau <- if (l <= length(planned)) planned[l] else tau_prev + block
    start <- pmax(tau_prev, ENT); act <- status == 0L
    target <- ncur; target[act] <- as.integer(pmin(cap[act], ncur[act] + speed[act] * pmax(tau - start[act], 0)))
    grow <- target > ncur
    tau_prev <- tau
    if (!any(grow)) next
    tnow <- ifelse(target >= cap, start + (cap - ncur) / speed, tau)   # time of this analysis per cohort
    ynew <- matrix(C[cbind(rep(seq_len(S), K), rep(seq_len(K), each = S), as.vector(target) + 1L)], S, K)
    nb <- target - ncur
    E[grow] <- E[grow] * mixture_factor(A[grow], B[grow], (ynew - ycur)[grow], nb[grow], p0)
    qm <- p0 + (1 - p0) * A / (A + B)                                    # prior mean, reported as theta1
    q_sum <- q_sum + colSums(qm * nb * grow); q_cnt <- q_cnt + colSums(nb * grow)
    ncur[grow] <- target[grow]; ycur[grow] <- ynew[grow]
    ## analysis: efficacy -> new frozen priors -> extension -> maximum -> hopeless
    decide <- function(cells, code) { status[cells] <<- code; time[cells] <<- tnow[cells] - ENT[cells] }
    decide(which(grow & E >= thr), 1L)
    rows <- which(rowSums(status == 0L) > 0)
    if (length(rows)) {
      pr <- prior(rows); upd <- status[rows, , drop = FALSE] == 0L
      a <- A[rows, , drop = FALSE]; a[upd] <- pr$a[upd]; A[rows, ] <- a
      b <- B[rows, , drop = FALSE]; b[upd] <- pr$b[upd]; B[rows, ] <- b
    }
    if (extend) {
      cand <- which(grow & status == 0L & !decided_ext & ncur >= NM)
      if (length(cand)) {
        go <- extension_power(A[cand], B[cand], E[cand], NE[cand], p0, thr) >= fl$cp_min
        decided_ext[cand] <- TRUE; extended[cand[go]] <- TRUE
        cap[cand[go]] <- (NM + NE)[cand[go]]; HZ[cand[!go]] <- NM[cand[!go]]
      }
    }
    decide(which(grow & status == 0L & ncur >= cap), 3L)
    if (settings$evalue$futility_hz) decide(which(grow & status == 0L & E * (1 / p0)^(HZ - ncur) < thr), 2L)
    if (acc) { a <- status == 0L; speed[a] <- ifelse(E[a] >= fl$m_speed, as.integer(fl$speed_up), 1L) }
  }
  list(status = status, n = ncur, e_value = E, theta1 = q_sum / pmax(q_cnt, 1), time = time,
       extended = extended, entry = entry)
}

## Calibrate Delta under the global null: smallest Delta with type I error <= alpha.
## "mean": one Delta (mean over cohorts); "per_group": one Delta per group of
## cohorts with the same analysis times.
calibrate_delta <- function(design, C_null, block, setup, settings, eng, cache) {
  rej <- function(D) colMeans(run_trials_bayes(design, C_null, block, setup$entry, settings, eng, cache, D)$status == 1L)
  bisect <- function(err) {
    lo <- 0.5; hi <- 1
    if (err(lo) <= settings$alpha) return(lo)
    while (hi - lo > settings$calibration$tol) { mid <- (lo + hi) / 2; if (err(mid) <= settings$alpha) hi <- mid else lo <- mid }
    hi
  }
  sig <- vapply(cohort_times(block, settings$n_max, setup$entry), paste, "", collapse = ",")
  grp <- if (setup$calibration == "per_group") match(sig, unique(sig)) else rep(1L, length(sig))
  D <- rep(1, length(grp))
  for (sweep in seq_len(if (max(grp) == 1) 1 else 2)) for (gi in unique(grp)) {
    D[grp == gi] <- bisect(function(d) { Dk <- D; Dk[grp == gi] <- d; mean(rej(Dk)[grp == gi]) })
  }
  D
}

## One design: all schedules and scenarios.
run_design <- function(dname, setup, settings, eng, banks, data) {
  design <- settings$designs[[dname]]; cache <- new_cache(); out <- list()
  for (sname in names(settings$schedules)) {
    block <- settings$schedules[[sname]]
    if (design$type == "bayes") Delta <- calibrate_delta(design, data$cal, block, setup, settings, eng, cache)
    if (design$type == "evalue") bk <- banks_for(design, block, setup, settings, banks)
    ext <- if (design$type == "evalue" && isTRUE(design$extend)) ext_banks_for(design, block, setup, settings, banks)
    for (scn in names(settings$scenarios)) {
      res <- switch(design$type,
        bayes   = run_trials_bayes(design, data$eval[[scn]], block, setup$entry, settings, eng, cache, Delta),
        mixture = run_trials_mixture(design, data$eval[[scn]], block, setup$entry, settings, eng, cache),
        evalue  = run_trials_evalue(design, data$eval[[scn]], block, setup$entry, settings, eng, cache, bk, ext))
      C <- data$eval[[scn]]; S <- dim(C)[1]; K <- dim(C)[2]
      res$y <- matrix(C[cbind(rep(seq_len(S), K), rep(seq_len(K), each = S), as.vector(res$n) + 1L)], S, K)
      out[[length(out) + 1]] <- summarise_run(res, settings$scenarios[[scn]], settings,
                                              list(design = dname, schedule = sname, scenario = scn))
    }
  }
  out
}

## -----------------------------------------------------------------------------
## 6. Summaries, references and figures
## -----------------------------------------------------------------------------
## Operating characteristics of one design / schedule / scenario.
##  sample size to reject: n_reject_H1 = mean n at rejection of active cohorts;
##    early_reject_H1 = P(active cohort rejected before n_max); n_early_reject_H1 = mean n of those
##  calendar time (from a cohort's entry; 1 = one patient at normal speed): time_H0 / time_H1 =
##    mean time to the decision; time_reject_H1 = mean time to rejection; duration = mean time
##    until every cohort has a decision
##  beyond n_max: extend_H0 / extend_H1 = P(cohort extended); reject_ext_H1 = P(active cohort
##    rejected after n_max)
summarise_run <- function(res, p, settings, id) {
  st <- res$status; n <- res$n; S <- nrow(st); K <- length(p)
  null <- p <= settings$p0 + 1e-9; alt <- p >= settings$p1 - 1e-9; rej <- st == 1L
  NM <- matrix(settings$n_max, S, K, byrow = TRUE); early <- rej & n < NM
  tm <- res$time; ext <- if (is.null(res$extended)) matrix(FALSE, S, K) else res$extended
  m <- function(x) if (length(x)) mean(x) else NA
  q <- function(x, pr) if (length(x) && !all(is.na(x))) unname(quantile(x, pr, type = 1, na.rm = TRUE)) else NA
  TP <- rowSums(rej[, alt, drop = FALSE]); FP <- rowSums(rej[, null, drop = FALSE])
  P <- matrix(p, S, K, byrow = TRUE); ph <- res$y / pmax(n, 1)                      # naive estimate at the decision
  lE <- if (is.null(res$e_value)) matrix(NA_real_, S, K) else log10(pmax(res$e_value, 1e-6))
  trial <- data.frame(id,
    typeI      = m(colMeans(rej)[null]),
    typeI_max  = if (any(null)) max(colMeans(rej)[null]) else NA,
    power      = m(colMeans(rej)[alt]),
    FWER       = if (any(null)) mean(rowSums(rej[, null, drop = FALSE]) > 0) else NA,
    ESS_H0     = m(colMeans(n)[null]),
    ESS_H1     = m(colMeans(n)[alt]),
    n_reject_H1 = m(n[, alt, drop = FALSE][rej[, alt, drop = FALSE]]),
    early_reject_H1   = m(colMeans(early)[alt]),
    n_early_reject_H1 = m(n[, alt, drop = FALSE][early[, alt, drop = FALSE]]),
    futility_H0 = m(colMeans(st == 2L)[null]),
    futility_H1 = m(colMeans(st == 2L)[alt]),
    theta1_H0  = if (is.null(res$theta1)) NA else m(res$theta1[null]),
    theta1_H1  = if (is.null(res$theta1)) NA else m(res$theta1[alt]),
    time_H0    = m(colMeans(tm)[null]),
    time_H1    = m(colMeans(tm)[alt]),
    time_reject_H1 = m(tm[, alt, drop = FALSE][rej[, alt, drop = FALSE]]),
    duration   = mean(apply(tm + matrix(res$entry, S, K, byrow = TRUE), 1, max)),
    extend_H0  = m(colMeans(ext)[null]),
    extend_H1  = m(colMeans(ext)[alt]),
    reject_ext_H1 = m(colMeans(rej & n > NM)[alt]),
    ## free statistics: computed from the outputs above, no extra simulation
    n_q50_H0   = q(n[, null], 0.5), n_q90_H0 = q(n[, null], 0.9),
    n_q50_H1   = q(n[, alt], 0.5),  n_q90_H1 = q(n[, alt], 0.9),
    trial_n    = mean(rowSums(n)), trial_n_q90 = q(rowSums(n), 0.9),
    true_rejections  = if (any(alt)) mean(TP) else NA,
    false_rejections = if (any(null)) mean(FP) else NA,
    FDR        = if (any(null)) mean(FP / pmax(TP + FP, 1)) else NA,
    power_any  = if (any(alt)) mean(TP > 0) else NA,
    power_all  = if (any(alt)) mean(TP == sum(alt)) else NA,
    all_correct = mean(TP == sum(alt) & FP == 0),
    n_per_discovery = if (any(alt) && mean(TP) > 0) mean(rowSums(n)) / mean(TP) else NA,
    log10E_H0  = q(lE[, null], 0.5), log10E_H1 = q(lE[, alt], 0.5),
    bias_H0    = m((ph - P)[, null]), bias_H1 = m((ph - P)[, alt]),
    bias_reject_H1 = m((ph - P)[, alt, drop = FALSE][rej[, alt, drop = FALSE]]),
    rmse_H1    = if (any(alt)) sqrt(mean(((ph - P)[, alt])^2)) else NA)
  ## cumulative curves: by patient (rejection, futility) and by calendar time (rejection)
  Tn <- max(settings$n_max + settings$n_ext) + 1
  cum <- function(grp, what, by_time = FALSE) {
    if (!any(grp)) return(NULL)
    x <- if (by_time) ceiling(tm - 1e-9) else n
    v <- x[, grp, drop = FALSE][(st == what)[, grp, drop = FALSE]] + 1L
    cnt <- tabulate(v, max(Tn, v))
    data.frame(id, group = if (identical(grp, null)) "H0" else "H1",
               measure = if (by_time) "reject_time" else if (what == 1L) "reject" else "futility",
               t = seq_along(cnt) - 1, value = cumsum(cnt) / (S * sum(grp)))
  }
  list(trial = trial, curves = rbind(cum(null, 1L), cum(null, 2L), cum(alt, 1L), cum(alt, 2L),
                                     cum(null, 1L, TRUE), cum(alt, 1L, TRUE)))
}

## Simon (1989) two-stage design: stop after n1 patients if <= r1 responses; reject
## H0 if > r responses out of n. Simon has no early stop for efficacy.
simon_reject <- function(n1, r1, n, r, p) {
  s1 <- (r1 + 1):n1
  sum(dbinom(s1, n1, p) * pbinom(r - s1, n - n1, p, lower.tail = FALSE))
}
## All designs with total n: for each (n1, r1) the smallest r with type I <= alpha.
simon_candidates <- function(n, p0, p1, alpha) {
  r <- 0:(n - 1); out <- list()
  for (n1 in seq_len(n - 1)) {
    x1 <- 0:n1
    rej <- function(p) {                         # rows r1 = 0..n1-1, columns r
      a <- dbinom(x1, n1, p) * outer(x1, r, function(x, rr) pbinom(rr - x, n - n1, p, lower.tail = FALSE))
      matrix(apply(a, 2, function(v) rev(cumsum(rev(v)))), n1 + 1)[-1, , drop = FALSE]
    }
    R0 <- rej(p0); R0[outer(0:(n1 - 1), r, ">")] <- Inf
    ok <- R0 <= alpha + 1e-12; has <- which(rowSums(ok) > 0)
    if (!length(has)) next
    rs <- max.col(ok[has, , drop = FALSE], "first"); r1 <- has - 1
    out[[n1]] <- data.frame(n1 = n1, r1 = r1, n = n, r = r[rs], typeI = R0[cbind(has, rs)],
                            power = rej(p1)[cbind(has, rs)], EN0 = n1 + (n - n1) * (1 - pbinom(r1, n1, p0)))
  }
  do.call(rbind, out)
}
## n given: the design with the smallest EN(p0) among those with power >= 1 - beta
## (or the most powerful one if there is none). n NULL: Simon's 'optimal' (smallest
## EN(p0)) or 'minimax' (smallest n, then smallest EN(p0)) design.
simon_design <- function(p0, p1, alpha, beta, n = NULL, type = "optimal", n_limit = 100) {
  best <- NULL
  for (nn in if (is.null(n)) 2:n_limit else n) {
    cand <- simon_candidates(nn, p0, p1, alpha)
    if (is.null(cand)) next
    ok <- cand[cand$power >= 1 - beta - 1e-12, ]
    if (nrow(ok)) {
      b <- ok[which.min(ok$EN0), ]
      if (is.null(best) || b$EN0 < best$EN0 - 1e-9) best <- b
      if (type == "minimax") break
    } else if (!is.null(n)) best <- cand[order(-cand$power, cand$EN0)[1], ]
    if (!is.null(best) && is.null(n) && nn >= 1.5 * best$n) break   # longer designs never had a smaller EN(p0) in Simon's tables
  }
  best
}
## Simon's design with n = n_max and the exact binomial test with n = n_max, at
## p0, p1 and the 'strong' rate. n_reject: mean n at rejection (no early rejection).
references <- function(settings) {
  N <- max(settings$n_max); sd <- simon_design(settings$p0, settings$p1, settings$alpha, settings$beta, n = N)
  r_b <- min(which(pbinom(0:N - 1, N, settings$p0, lower.tail = FALSE) <= settings$alpha)) - 1
  do.call(rbind, lapply(c(settings$p0, settings$p1, settings$scenarios$strong[1]), function(p) rbind(
    data.frame(design = sprintf("Simon (n1 = %d, r1 = %d, r = %d)", sd$n1, sd$r1, sd$r), p = p,
               reject = simon_reject(sd$n1, sd$r1, N, sd$r, p),
               ESS = sd$n1 + (N - sd$n1) * (1 - pbinom(sd$r1, sd$n1, p)), n1 = sd$n1,
               futility = pbinom(sd$r1, sd$n1, p), n_reject = N),
    data.frame(design = sprintf("Binomial test (reject if >= %d)", r_b), p = p,
               reject = pbinom(r_b - 1, N, p, lower.tail = FALSE), ESS = N, n1 = NA, futility = 0, n_reject = N))))
}

## Colours and line types: solid = Bayesian, dashed = e-value, long-dashed = e-value with borrowing.
## Line type = family: framework e-values solid, design-optimal e-values dashed, Bayesian
## comparators dotted. Flexible variants have the colour of their base design and line
## type dot-dash (accelerated), two-dash (extended) or long-dash (both).
design_style <- function(settings, dl) {
  pal <- c("#2a78d6", "#1baf7a", "#6250d6", "#eda100", "#e34948", "#7a1f5c", "#eb6834", "#e87ba4",
           "#52514e", "#9c6b30", "#17becf", "#8c8c00")
  dl <- intersect(dl, names(settings$designs))
  bases <- names(settings$designs)[vapply(settings$designs, function(d) is.null(d$base), TRUE)]
  base <- vapply(dl, function(d) settings$designs[[d]]$base %||% d, "")
  cols <- setNames(pal[(match(base, bases) - 1) %% length(pal) + 1], dl)
  ltys <- setNames(vapply(dl, function(d) { x <- settings$designs[[d]]; a <- isTRUE(x$accelerate); e <- isTRUE(x$extend)
    if (a && e) "longdash" else if (a) "dotdash" else if (e) "twodash"
    else if (x$type == "mixture") "solid" else if (x$type == "bayes") "13" else "22" }, ""), dl)
  list(cols = c(cols, "Simon two-stage" = "#0b0b0b"), ltys = c(ltys, "Simon two-stage" = "dotted"))
}
theme_oc <- function() {
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom", legend.title = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"),
    plot.background = ggplot2::element_rect(fill = "#fcfcfb", colour = NA))
}

## Figure 1 of the e-value paper: cumulative type I error, power and futility
## stopping by patient (rows) for each analysis schedule (columns).
plot_cumulative <- function(curves, setup_name, settings, refs) {
  library(ggplot2)
  dl <- settings$setups[[setup_name]]$designs; sty <- design_style(settings, dl)
  rows <- c("Type I error\n(global null)", "Power\n(global alternative)",
            "Futility stop\n(global null)", "Futility stop\n(global alternative)")
  pick <- function(scn, grp, meas, i) { d <- curves[curves$scenario == scn & curves$group == grp & curves$measure == meas, ]
                                        if (nrow(d)) d$row <- rows[i]; d }
  d <- rbind(pick(settings$scenario_null, "H0", "reject", 1), pick(settings$scenario_alt, "H1", "reject", 2),
             pick(settings$scenario_null, "H0", "futility", 3), pick(settings$scenario_alt, "H1", "futility", 4))
  N <- max(settings$n_max); t <- 0:N; sim <- refs[grepl("Simon", refs$design), ]
  s0 <- sim[sim$p == settings$p0, ]; s1 <- sim[sim$p == settings$p1, ]
  d <- rbind(d[, c("design", "schedule", "row", "t", "value")], data.frame(
    design = "Simon two-stage", schedule = names(settings$schedules)[3], row = rep(rows, each = N + 1), t = t,
    value = c((t >= N) * s0$reject, (t >= N) * s1$reject, (t >= s0$n1) * s0$futility, (t >= s1$n1) * s1$futility)))
  d$row <- factor(d$row, levels = rows); d$schedule <- factor(d$schedule, levels = names(settings$schedules))
  d$design <- factor(d$design, levels = c(dl, "Simon two-stage"))
  binom_power <- refs$reject[grepl("Binomial", refs$design) & refs$p == settings$p1]
  ggplot(d, aes(t, value, colour = design, linetype = design)) +
    geom_hline(data = data.frame(row = factor(rows[1:2], levels = rows), y = c(settings$alpha, binom_power)),
               aes(yintercept = y), colour = "grey40", linetype = "dotdash", linewidth = 0.4) +
    geom_step(linewidth = 0.6) + facet_grid(row ~ schedule, scales = "free_y") +
    scale_colour_manual(values = sty$cols) + scale_linetype_manual(values = sty$ltys) +
    labs(x = "Patients in cohort", y = "Cumulative probability",
         title = sprintf("%s trial (%s): cumulative operating characteristics per cohort",
                         tools::toTitleCase(setup_name), settings$label),
         subtitle = sprintf("H0: p <= %.2f vs H1: p >= %.2f, alpha = %.2f. Dash-dot: alpha and power of the exact binomial test.",
                            settings$p0, settings$p1, settings$alpha)) +
    theme_oc() + guides(colour = guide_legend(nrow = 4))
}

## Type I error, power and expected sample size by scenario.
plot_scenarios <- function(trial, setup_name, settings, meas = NULL, what = "results by scenario") {
  library(ggplot2)
  dl <- settings$setups[[setup_name]]$designs; sty <- design_style(settings, dl)
  meas <- meas %||% c(typeI = "Type I error", power = "Power", early_reject_H1 = "Rejected before\nn_max (active)",
            ESS_H0 = "Expected n,\nnull cohorts", ESS_H1 = "Expected n,\nactive cohorts",
            n_reject_H1 = "Mean n at\nrejection (active)", time_reject_H1 = "Mean time to\nrejection (active)")
  d <- do.call(rbind, lapply(names(meas), function(v) data.frame(trial[, c("design", "schedule", "scenario")], row = meas[[v]], value = trial[[v]])))
  d <- d[!is.na(d$value), ]
  d$row <- factor(d$row, levels = meas); d$schedule <- factor(d$schedule, levels = names(settings$schedules))
  d$design <- factor(d$design, levels = dl)
  d$scen <- factor(sub(" .*", "", d$scenario), levels = unique(sub(" .*", "", names(settings$scenarios))))
  alpha_line <- if ("typeI" %in% names(meas))
    geom_hline(data = data.frame(row = factor(meas[["typeI"]], levels = meas), y = settings$alpha),
               aes(yintercept = y), colour = "grey40", linetype = "dotdash", linewidth = 0.4)
  ggplot(d, aes(scen, value, colour = design, linetype = design, group = design)) + alpha_line +
    geom_line(linewidth = 0.5) + geom_point(size = 1.5) + facet_grid(row ~ schedule, scales = "free_y") +
    scale_colour_manual(values = sty$cols) + scale_linetype_manual(values = sty$ltys) +
    labs(x = sprintf("Number of active cohorts (p = %.2f, the others p = %.2f)", settings$p1, settings$p0), y = NULL,
         title = sprintf("%s trial (%s): %s", tools::toTitleCase(setup_name), settings$label, what),
         caption = sprintf("mixed: response rates spread from %.3f to %.3f; strong: all at %.3f",
                           min(settings$scenarios$mixed), max(settings$scenarios$mixed), settings$scenarios$strong[1])) +
    theme_oc() + guides(colour = guide_legend(nrow = 4))
}

## Cumulative probability of rejection by calendar time (time from the cohort's
## entry; 1 = one patient at normal speed): type I error under the global null,
## power under the global alternative and in the 'strong' scenario.
plot_time <- function(curves, setup_name, settings) {
  library(ggplot2)
  dl <- settings$setups[[setup_name]]$designs; sty <- design_style(settings, dl)
  rows <- c("Type I error\n(global null)", "Power\n(global alternative)", "Power\n(strong)")
  pick <- function(scn, grp, i) { d <- curves[curves$scenario == scn & curves$group == grp & curves$measure == "reject_time", ]
                                  if (nrow(d)) d$row <- rows[i]; d }
  d <- rbind(pick(settings$scenario_null, "H0", 1), pick(settings$scenario_alt, "H1", 2), pick("strong", "H1", 3))
  d$row <- factor(d$row, levels = rows); d$schedule <- factor(d$schedule, levels = names(settings$schedules))
  d$design <- factor(d$design, levels = dl)
  ggplot(d, aes(t, value, colour = design, linetype = design)) +
    geom_hline(data = data.frame(row = factor(rows[1], levels = rows), y = settings$alpha),
               aes(yintercept = y), colour = "grey40", linetype = "dotdash", linewidth = 0.4) +
    geom_vline(xintercept = max(settings$n_max), colour = "grey60", linewidth = 0.3) +
    geom_step(linewidth = 0.6) + facet_grid(row ~ schedule, scales = "free_y") +
    scale_colour_manual(values = sty$cols) + scale_linetype_manual(values = sty$ltys) +
    labs(x = "Calendar time from the cohort's entry (1 = one patient at normal speed); grey line: n_max",
         y = "Cumulative probability",
         title = sprintf("%s trial (%s): rejection by calendar time", tools::toTitleCase(setup_name), settings$label)) +
    theme_oc() + guides(colour = guide_legend(nrow = 4))
}

## Measures of fig5 (by scenario): decisions per trial, trial size, estimation, evidence.
meas_decisions <- c(true_rejections = "True rejections\nper trial", false_rejections = "False rejections\nper trial",
                    FDR = "False discovery\nrate", all_correct = "P(all decisions\ncorrect)",
                    trial_n = "Expected total n\nof the trial", bias_reject_H1 = "Bias of y/n,\nrejected active cohorts",
                    log10E_H1 = "Median log10 e-value\nat decision (active)")

## Inclusion of a design (shape in fig4): flexible variants share the colour of their base design.
design_inclusion <- function(d, settings) vapply(d, function(x) {
  if (x == "Simon two-stage") return("Simon")
  z <- settings$designs[[x]]; a <- isTRUE(z$accelerate); e <- isTRUE(z$extend)
  if (a && e) "accelerated + extended" else if (a) "accelerated" else if (e) "extended" else "fixed inclusion" }, "")

## fig4: what each design buys. Power under the global alternative against the expected n of
## a null cohort (global null) and against the mean time to rejection; Simon as reference.
plot_tradeoff <- function(trial, setup_name, settings, refs) {
  library(ggplot2)
  dl <- settings$setups[[setup_name]]$designs; sty <- design_style(settings, dl)
  pan <- c("Expected n of a null cohort (global null)", "Mean time to rejection (global alternative)")
  nul <- trial[trial$scenario == settings$scenario_null, c("design", "schedule", "ESS_H0")]
  alt <- trial[trial$scenario == settings$scenario_alt, c("design", "schedule", "power", "time_reject_H1")]
  d <- merge(nul, alt, by = c("design", "schedule"))
  dd <- rbind(data.frame(d[, c("design", "schedule", "power")], panel = pan[1], x = d$ESS_H0),
              data.frame(d[, c("design", "schedule", "power")], panel = pan[2], x = d$time_reject_H1))
  s <- refs[grepl("^Simon", refs$design), ]; at <- function(p) s[abs(s$p - p) < 1e-9, ]
  dd <- rbind(dd, do.call(rbind, lapply(names(settings$schedules), function(sn) data.frame(
    design = "Simon two-stage", schedule = sn, power = at(settings$p1)$reject, panel = pan,
    x = c(at(settings$p0)$ESS, max(settings$n_max))))))
  dd <- dd[!is.na(dd$x) & !is.na(dd$power), ]
  dd$schedule <- factor(dd$schedule, levels = names(settings$schedules)); dd$panel <- factor(dd$panel, levels = pan)
  dd$inclusion <- factor(design_inclusion(as.character(dd$design), settings),
                         levels = c("fixed inclusion", "accelerated", "extended", "accelerated + extended", "Simon"))
  dd$design <- factor(dd$design, levels = c(dl, "Simon two-stage"))
  ggplot(dd, aes(x, power, colour = design, shape = inclusion)) + geom_point(size = 2.6) +
    facet_grid(schedule ~ panel, scales = "free_x") +
    scale_colour_manual(values = sty$cols) +
    scale_shape_manual(values = c("fixed inclusion" = 16, "accelerated" = 17, "extended" = 15, "accelerated + extended" = 18, "Simon" = 4)) +
    labs(x = NULL, y = "Power (global alternative)",
         title = sprintf("%s trial (%s): power against sample size and time", tools::toTitleCase(setup_name), settings$label),
         caption = "Up and to the left is better. Simon: no early rejection, so its time to rejection is n_max.") +
    theme_oc() + theme(legend.box = "vertical") +
    guides(colour = guide_legend(ncol = 4, override.aes = list(shape = 16)), shape = guide_legend(nrow = 1))
}

## Effect of each design element, as differences in the results table (same schedule
## and scenario): borrowing (vs the same design without borrowing), flexible inclusion
## (vs its base design) and every design vs the primary design (settings$primary).
design_effects <- function(trial, settings, dl) {
  metrics <- c("typeI_max", "power", "ESS_H0", "ESS_H1", "trial_n", "time_reject_H1", "duration",
               "FDR", "all_correct", "n_per_discovery", "bias_reject_H1")
  counterpart <- function(d) {
    x <- settings$designs[[d]]
    if (!is.null(x$base)) return(c(x$base, "flexible inclusion"))
    if (x$type == "evalue" && isTRUE(x$borrow)) return(c("e-based", "borrowing"))
    if (x$type == "mixture" && x$pi_rule != "none") return(c("e-mixture (no borrowing)", "borrowing"))
    if (x$type == "bayes" && x$pi_rule != "none") return(c("Bayes Umbrella", "borrowing"))
    NULL
  }
  key <- c("schedule", "scenario"); out <- list()
  for (d in dl) {
    cmp <- list(counterpart(d)); if (d != settings$primary) cmp <- c(cmp, list(c(settings$primary, "vs primary design")))
    for (cc in cmp) {
      if (is.null(cc) || !(cc[1] %in% dl)) next
      m <- merge(trial[trial$design == d, c(key, metrics)], trial[trial$design == cc[1], c(key, metrics)],
                 by = key, suffixes = c("", ".ref"))
      dd <- m[, key]; for (v in metrics) dd[[paste0("d_", v)]] <- m[[v]] - m[[paste0(v, ".ref")]]
      out[[length(out) + 1]] <- data.frame(design = d, compared_with = cc[1], effect = cc[2], dd)
    }
  }
  do.call(rbind, out)
}

## Overview of one series and one set-up, with Simon's design as reference.
##  series "grid":  x = K, columns = n_max, blocked schedule
##  other series:   x = p0 (hypotheses), columns = the three schedules
## Rows: worst type I error, power, P(rejected before n_max), mean n at rejection,
## expected n of null and active cohorts (global null / global alternative).
plot_overview <- function(all, series, setup_name, settings) {
  library(ggplot2)
  sched <- c("Fully sequential", "Blocked (n_max / 3)", "Two-stage (n_max / 2)")
  a <- all[all$series == series & all$setup %in% c(setup_name, "reference"), ]
  is_ref <- a$setup == "reference"
  if (series == "grid") {
    a <- a[is_ref | a$schedule_pos == 2, ]
    a$x <- factor(a$K, levels = sort(unique(a$K)))
    a$col <- factor(sprintf("n_max = %d", a$n_max), levels = sprintf("n_max = %d", sort(unique(a$n_max))))
    xlab <- "Number of cohorts K"; ttl <- sprintf("p0 = %.2f, p1 = %.2f, blocked analyses", a$p0[1], a$p1[1])
  } else {
    ref <- a[is_ref, ]; a <- a[!is_ref, ]                   # Simon in every schedule column
    a <- rbind(a, do.call(rbind, lapply(1:3, function(j) { r <- ref; r$schedule_pos <- j; r })))
    lab <- sprintf("%.2f\n%.2f\n%d", a$p0, a$p1, a$n_max)       # evenly spaced hypotheses
    a$x <- factor(lab, levels = unique(lab[order(a$p0)])); a$col <- factor(sched[a$schedule_pos], levels = sched)
    xlab <- "p0 (top), p1 (middle) and n_max = Simon's optimal n (bottom)"
    ttl <- sprintf("K = %d, the hypotheses of Simon (1989)", a$K[1])
  }
  worst <- aggregate(typeI_max ~ design + x + col, a, max)
  alt <- a[a$is_alt, ]; nul <- a[a$is_null, ]
  meas <- c("Worst type I error\n(any scenario)", "Power\n(global alternative)", "Rejected before n_max\n(global alternative)",
            "Mean n at rejection\n(global alternative)", "Mean time to rejection\n(global alternative)",
            "Expected n, null cohort\n(global null)", "Expected n, active cohort\n(global alternative)")
  part <- function(x, i, v) data.frame(design = x$design, x = x$x, col = x$col, row = meas[i], value = v)
  d <- rbind(part(worst, 1, worst$typeI_max), part(alt, 2, alt$power), part(alt, 3, alt$early_reject_H1),
             part(alt, 4, alt$n_reject_H1), part(alt, 5, alt$time_reject_H1), part(nul, 6, nul$ESS_H0),
             part(alt, 7, alt$ESS_H1))
  d <- d[!is.na(d$value), ]
  dl <- c(settings$setups[[setup_name]]$designs, "Simon two-stage"); sty <- design_style(settings, dl)
  d$row <- factor(d$row, levels = meas); d$design <- factor(d$design, levels = dl)
  ggplot(d, aes(x, value, colour = design, linetype = design, group = design)) +
    geom_hline(data = data.frame(row = factor(meas[1], levels = meas), y = settings$alpha),
               aes(yintercept = y), colour = "grey40", linetype = "dotdash", linewidth = 0.4) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.5) + facet_grid(row ~ col, scales = "free_y") +
    scale_colour_manual(values = sty$cols) + scale_linetype_manual(values = sty$ltys) +
    labs(x = xlab, y = NULL, title = sprintf("%s trial, series '%s': %s", tools::toTitleCase(setup_name), series, ttl)) +
    theme_oc() + guides(colour = guide_legend(nrow = 4))
}

## Save a figure as PDF (works on every machine, also compute nodes without a
## display) and as PNG when the machine supports it.
save_plot <- function(p, file, width, height) {
  ggplot2::ggsave(paste0(file, ".pdf"), p, width = width, height = height)
  args <- list(paste0(file, ".png"), p, width = width, height = height, dpi = 150)
  if (capabilities("cairo")) args$type <- "cairo"
  tryCatch(do.call(ggplot2::ggsave, args),
           error = function(e) message("PNG not written (", conditionMessage(e), "); the PDF is there"))
}

## -----------------------------------------------------------------------------
## 7. One run (one row of the table of runs in Setting.R)
## -----------------------------------------------------------------------------
## Fill in everything that depends on the run: hypotheses (p0, p1), number of
## cohorts K and cohort size n_max (run = one row of complete_runs()).
complete_settings <- function(base, run) {
  K <- run$K; n_max <- run$n_max
  s <- base$for_hypothesis(base, run$p0, run$p1)
  s$series <- run$series; s$K <- K; s$n_max <- rep(n_max, K)
  s$n_ext <- rep(as.integer(round(base$flex$ext_frac * n_max)), K)    # extra patients when extended
  s$schedules <- base$make_schedules(n_max)
  s$scenarios <- base$make_scenarios(K, run$p0, run$p1)
  s$scenario_null <- names(s$scenarios)[1]
  s$scenario_alt  <- names(s$scenarios)[grepl("global alt", names(s$scenarios))]
  s$setups$basket$entry   <- rep(0, K)
  s$setups$umbrella$entry <- base$make_entry(K, n_max)
  for (d in names(s$designs)) if (isTRUE(s$designs[[d]]$borrow)) {           # borrowing defaults
    s$designs[[d]]$summary  <- s$designs[[d]]$summary  %||% base$borrow$summary
    s$designs[[d]]$quantile <- s$designs[[d]]$quantile %||% base$borrow$quantile
    s$designs[[d]]$redesign <- s$designs[[d]]$redesign %||% base$borrow$redesign
  }
  s$label <- sprintf("p0 = %.2f, p1 = %.2f, K = %d, n_max = %d", run$p0, run$p1, K, n_max)
  s$sim$out_dir <- file.path(base$sim$out_dir, run$name)
  s
}

## Run every set-up, design, schedule and scenario of one setting; save its
## tables and figures in its own folder and return all rows of the results table.
run_one_setting <- function(settings) {
  t0 <- Sys.time()
  dir.create(settings$sim$out_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(settings$sim$seed)
  eng   <- build_engine(settings)
  banks <- get_banks(settings)
  refs  <- references(settings)
  K <- settings$K; N <- max(settings$n_max)
  U     <- array(runif(settings$sim$n_sim * K * N), c(settings$sim$n_sim, K, N))
  U_cal <- array(runif(settings$sim$n_cal * K * N), c(settings$sim$n_cal, K, N))
  NE    <- max(settings$n_ext)                    # patients beyond n_max (extended cohorts), drawn last
  U     <- array(c(U, runif(settings$sim$n_sim * K * NE)), c(settings$sim$n_sim, K, N + NE))
  data  <- list(eval = lapply(settings$scenarios, function(p) make_data(U, p)),
                cal  = make_data(U_cal, settings$scenarios[[settings$scenario_null]]))
  out <- list()
  for (su in names(settings$setups)) {
    setup <- settings$setups[[su]]
    res <- parallel::mclapply(setup$designs, run_design, setup = setup, settings = settings,
                              eng = eng, banks = banks, data = data, mc.cores = settings$sim$n_cores)
    failed <- vapply(res, inherits, TRUE, "try-error")
    if (any(failed)) stop("design failed: ", paste(setup$designs[failed], collapse = ", "), "\n", res[failed][[1]])
    res <- unlist(res, recursive = FALSE)
    trial  <- do.call(rbind, lapply(res, `[[`, "trial"))
    curves <- do.call(rbind, lapply(res, `[[`, "curves"))
    write.csv(trial,  file.path(settings$sim$out_dir, paste0(su, "_results.csv")), row.names = FALSE)
    write.csv(curves, file.path(settings$sim$out_dir, paste0(su, "_curves.csv")),  row.names = FALSE)
    save_plot(plot_cumulative(curves, su, settings, refs),
              file.path(settings$sim$out_dir, paste0(su, "_fig1_cumulative")), width = 11, height = 10.5)
    save_plot(plot_scenarios(trial, su, settings),
              file.path(settings$sim$out_dir, paste0(su, "_fig2_scenarios")), width = 11, height = 13)
    save_plot(plot_time(curves, su, settings),
              file.path(settings$sim$out_dir, paste0(su, "_fig3_time")), width = 11, height = 8)
    save_plot(plot_tradeoff(trial, su, settings, refs),
              file.path(settings$sim$out_dir, paste0(su, "_fig4_tradeoff")), width = 11, height = 10)
    save_plot(plot_scenarios(trial, su, settings, meas_decisions, "decisions, trial size and estimation"),
              file.path(settings$sim$out_dir, paste0(su, "_fig5_decisions")), width = 11, height = 15)
    write.csv(design_effects(trial, settings, setup$designs),
              file.path(settings$sim$out_dir, paste0(su, "_design_effects.csv")), row.names = FALSE)
    out[[su]] <- cbind(setup = su, series = settings$series, p0 = settings$p0, p1 = settings$p1, K = K, n_max = N,
                       schedule_pos = match(trial$schedule, names(settings$schedules)),
                       is_null = trial$scenario == settings$scenario_null,
                       is_alt = trial$scenario == settings$scenario_alt, trial)
  }
  write.csv(refs, file.path(settings$sim$out_dir, "references.csv"), row.names = FALSE)
  out$reference <- simon_rows(out[[1]][1:2, ], refs, settings)
  message(sprintf("  %s done in %.1f min", settings$label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  do.call(rbind, out)
}

## Simon's design as two rows of the results table (global null, global alternative),
## for the overview figures. Simon cannot reject early: n at rejection = n_max.
simon_rows <- function(template, refs, settings) {
  s <- refs[grepl("^Simon", refs$design), ]; at <- function(p) s[abs(s$p - p) < 1e-9, ]
  r <- template; r[] <- NA
  r$setup <- "reference"; r$series <- settings$series; r$p0 <- settings$p0; r$p1 <- settings$p1
  r$K <- settings$K; r$n_max <- max(settings$n_max); r$schedule_pos <- 0; r$schedule <- s$design[1]
  r$design <- "Simon two-stage"; r$scenario <- c(settings$scenario_null, settings$scenario_alt)
  r$is_null <- c(TRUE, FALSE); r$is_alt <- c(FALSE, TRUE)
  r$typeI[1] <- r$typeI_max[1] <- at(settings$p0)$reject; r$ESS_H0[1] <- at(settings$p0)$ESS
  r$futility_H0[1] <- at(settings$p0)$futility
  r$power[2] <- at(settings$p1)$reject; r$ESS_H1[2] <- at(settings$p1)$ESS; r$futility_H1[2] <- at(settings$p1)$futility
  r$n_reject_H1[2] <- max(settings$n_max); r$early_reject_H1[2] <- 0
  r$time_H0[1] <- r$ESS_H0[1]; r$time_H1[2] <- r$ESS_H1[2]; r$time_reject_H1[2] <- max(settings$n_max)  # normal speed
  r$extend_H0[1] <- r$extend_H1[2] <- r$reject_ext_H1[2] <- 0
  r
}

## The table of runs with Simon's optimal n filled in where n_max = NA, and a
## folder name per run.
complete_runs <- function(settings) {
  r <- settings$runs
  for (i in which(is.na(r$n_max))) r$n_max[i] <- simon_design(r$p0[i], r$p1[i], settings$alpha, settings$beta)$n
  r$name <- sprintf("%s_p%02.0f-%02.0f_K%d_n%d", r$series, 100 * r$p0, 100 * r$p1, r$K, r$n_max)
  r
}
