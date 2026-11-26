# =============================================================================
# Design-optimal e-values for single-arm trials with binary outcomes
#
# R re-implementation of the dynamic programme in
#   Baas, ter Schure & van Rosmalen (2026). Adaptive clinical trials based on
#   design-optimal e-values with automatic curtailment: an application to
#   single-arm trials with binary data. arXiv:2605.28653.
# Equation numbers in the comments refer to that paper.
#
# Setting (Section 2). Y_1, ..., Y_n iid Bernoulli(theta), H0: theta <= theta0.
# A bettor starts with capital M_0 = 1 and bets a fraction B_t of M_{t-1} on
# the event {Y_t = 1}:
#     M_t = M_{t-1} * (1 + B_t * (Y_t / theta0 - 1))                    (2)
# H0 is rejected once M_t >= 1 / alpha. Any predictable bets give a test
# supermartingale, so type I error control never depends on how B_t is chosen.
# In the "ticket" form prod_t (1 + lambda_t * (Y_t - theta0)) the same e-process
# has lambda_t = B_t / theta0.
#
# Optimal bets (Section 3 and Appendix A). The state is X_t = (t, M_t). The
# e-value is discretised to a grid M-hat containing 0, 1 and 1/alpha, and every
# transition is rounded *down* onto the grid:
#     (t, M_t) -> (t + 1, floor_grid(M_t * (1 + B_{t+1} (Y_{t+1}/theta0 - 1))))  (15)
# with 0 and 1/alpha absorbing. Bets come from a grid B-hat. The Bellman
# equations (7)-(8) are solved backwards in t under theta = theta1; ties go to
# the smallest bet. Three objectives are implemented:
#   design_pmax()    maximise power at the horizon n                  (3.1.1)
#   design_essmin()  minimise E[min(n + 1, tau_{1/alpha})]            (3.1.2)
#   design_ebased()  minimise expected sample size with optimal futility
#                    stopping, subject to power >= 1 - beta           (3.2, A.2)
# and, for comparison, design_grow() (constant Kelly bet, Section 2.3).
#
# Using a design on data. The grid process M-hat ("shadow" e-value) only
# selects the bets: the bet for participant t + 1 is looked up at M-hat_t. The
# reported e-value M-tilde applies the same bets to the data without rounding,
# so M-tilde_t >= M-hat_t on every path (Appendix A.1). Power computed on the
# grid is therefore a lower bound for the real e-process, and expected sample
# size an upper bound. e_process() and simulate_design() track both processes.
#
# Blocked designs (Section 3.3). `analysis_times` lists the participant counts
# at which data are analysed (default: after every participant). Bets are still
# placed participant by participant in arrival order, but stopping (efficacy,
# futility) and recruitment costs only happen at analysis times.
#
# Base R only.
# =============================================================================


# ---- Grids and small helpers -------------------------------------------------

#' e-value grid of Appendix A.1: 0, `n_log` log-spaced values from `m_min` to
#' 1 - 2 * machine epsilon, and `n_lin` equidistant values from 1 to 1/alpha.
baas_m_grid <- function(alpha, n_log = 1000, n_lin = 1000, m_min = 1e-5) {
  eps <- .Machine$double.eps
  low <- exp(seq(log(m_min), log(1 - 2 * eps), length.out = n_log))
  low[c(1, n_log)] <- c(m_min, 1 - 2 * eps)
  high <- seq(1, 1 / alpha, length.out = n_lin)
  high[n_lin] <- 1 / alpha
  c(0, low, high)
}

#' Bet grid of Appendix A.1: {0, 1e-4, 1e-3, 0.01, 0.02, ..., 0.99, 0.999, 0.9999, 1}.
baas_b_grid <- function() {
  c(0, 1e-4, 1e-3, (1:99) / 100, 0.999, 0.9999, 1)
}

#' Kelly / GROW bet (eq. 4).
kelly_bet <- function(theta0, theta1) {
  (theta1 - theta0) / (1 - theta0)
}

#' Hopeless zone (Remark 3): from e-value m at time t, 1/alpha cannot be
#' reached by participant n even when betting everything on every remaining
#' participant.
hopeless <- function(m, t, n, theta0, alpha) {
  m < theta0^(n - t) / alpha
}

# Sum weights `w` into bins `idx` (1..nbins).
.accumulate <- function(idx, w, nbins) {
  out <- numeric(nbins)
  keep <- w > 0
  if (any(keep)) {
    s <- rowsum(w[keep], idx[keep])
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}


# ---- Problem set-up ----------------------------------------------------------

.setup <- function(n, theta0, theta1, alpha, analysis_times, m_grid, b_grid) {
  stopifnot(
    length(n) == 1, n >= 1, n == round(n),
    theta0 > 0, theta0 < 1, theta1 > 0, theta1 < 1,
    alpha > 0, alpha < 1,
    !is.unsorted(m_grid, strictly = TRUE),
    m_grid[1] == 0, any(m_grid == 1), m_grid[length(m_grid)] == 1 / alpha,
    !is.unsorted(b_grid, strictly = TRUE), b_grid[1] == 0, max(b_grid) <= 1
  )
  n <- as.integer(n)
  analysis_times <- sort(unique(c(as.integer(analysis_times), n)))
  if (any(analysis_times < 1 | analysis_times > n)) {
    stop("`analysis_times` must lie in 1..n")
  }

  # Decision epochs: t = 0 and every interim analysis. At an epoch the next
  # block of participants is recruited (or, for design_ebased(), the trial may
  # stop for futility). block_after[t + 1] is that block's size (0 if t is not
  # an epoch).
  epochs <- c(0L, analysis_times[analysis_times < n])
  block_after <- integer(n)
  block_after[epochs + 1L] <- diff(c(0L, analysis_times))

  nm <- length(m_grid)
  nb <- length(b_grid)
  top <- nm
  # Transition indices for every (state, bet): floor onto the grid; values
  # >= 1/alpha map to the top state. 1/alpha is absorbing (Section 3.1).
  up <- matrix(findInterval(outer(m_grid, 1 + b_grid * (1 / theta0 - 1)), m_grid), nm, nb)
  down <- matrix(findInterval(outer(m_grid, 1 - b_grid), m_grid), nm, nb)
  up[top, ] <- top
  down[top, ] <- top

  list(
    n = n, theta0 = theta0, theta1 = theta1, alpha = alpha,
    analysis_times = analysis_times, is_epoch = block_after > 0L,
    block_after = block_after, m_grid = m_grid, b_grid = b_grid,
    up = up, down = down
  )
}


# ---- Backward induction ------------------------------------------------------

# Solve the Bellman equations (7)-(8) for the discretised MDP (15), written as
# cost minimisation. C_t(m) is the expected cost-to-go from state (t, m):
#   "pmax":   C_n(m) = I(m < 1/alpha); no stage costs. C_t = 1 - P(reject).
#   "essmin": stage cost n_{t+1} * I(m < 1/alpha) at each decision epoch
#             (eq. 12 and Section 3.3) plus the final cost I(m_n < 1/alpha).
#             C_0 = E[min(n + 1, tau_{1/alpha})] when fully sequential.
#   "ebased": stage cost n_{t+1} for each recruited block, a futility stop
#             at decision epochs that costs `lambda`, and terminal cost
#             lambda * I(m_n < 1/alpha). This is the Lagrangian of (14):
#             E[sample size] + lambda * P(no rejection).
# Ties (up to a relative tolerance of 1e-10) go to the smallest bet; for
# "ebased", a tie between stopping and continuing goes to stopping.
.solve_dp <- function(s, objective = c("pmax", "essmin", "ebased"), lambda = NULL) {
  objective <- match.arg(objective)
  if (objective == "ebased" && (is.null(lambda) || lambda < 0)) {
    stop("objective 'ebased' needs lambda >= 0")
  }
  n <- s$n
  nm <- length(s$m_grid)
  nb <- length(s$b_grid)
  top <- nm
  below <- c(rep(1, nm - 1), 0) # I(m < 1/alpha)
  rows <- seq_len(nm)
  fail_cost <- if (objective == "ebased") lambda else 1

  cost <- matrix(NA_real_, nm, n + 1)
  bet_idx <- matrix(NA_integer_, nm, n)
  stop_mat <- matrix(FALSE, nm, n)

  C <- fail_cost * below
  cost[, n + 1] <- C
  for (t in (n - 1):0) {
    Q <- s$theta1 * C[s$up] + (1 - s$theta1) * C[s$down]
    dim(Q) <- c(nm, nb)
    tol <- 1e-10 * max(1, max(C))
    qmin <- Q[cbind(rows, max.col(-Q, ties.method = "first"))]
    j <- max.col((Q <= qmin + tol) + 0, ties.method = "first")
    Cn <- Q[cbind(rows, j)]

    if (s$is_epoch[t + 1] && objective != "pmax") {
      Cn <- Cn + s$block_after[t + 1] * below
    }
    if (s$is_epoch[t + 1] && objective == "ebased") {
      st <- below == 1 & lambda <= Cn + tol
      Cn[st] <- lambda
      stop_mat[st, t + 1] <- TRUE
      j[st] <- NA_integer_
    }
    j[c(1, top)] <- NA_integer_ # absorbing states 0 and 1/alpha: no bet
    bet_idx[, t + 1] <- j
    cost[, t + 1] <- Cn
    C <- Cn
  }
  list(bet_idx = bet_idx, stop = stop_mat, cost = cost)
}

.new_design <- function(s, type, bet, stop_mat, cost, extra = list()) {
  dimnames(bet) <- dimnames(stop_mat) <- list(NULL, participant = seq_len(s$n))
  if (!is.null(cost)) dimnames(cost) <- list(NULL, t = 0:s$n)
  d <- c(
    list(type = type),
    s[c("n", "theta0", "theta1", "alpha", "analysis_times", "is_epoch",
        "block_after", "m_grid", "b_grid")],
    list(bet = bet, stop = stop_mat, cost = cost),
    extra
  )
  structure(d, class = "evalue_design")
}

.design_from_dp <- function(s, type, sol, extra = list()) {
  bet <- matrix(s$b_grid[sol$bet_idx], nrow(sol$bet_idx), ncol(sol$bet_idx))
  .new_design(s, type, bet, sol$stop, sol$cost, extra)
}


# ---- Public design constructors ----------------------------------------------
#
# All return an object of class "evalue_design" with
#   bet[i, k]   bet B_k for participant k given grid e-value m_grid[i] after
#               k - 1 participants (NA in the absorbing states 0 and 1/alpha,
#               and where design_ebased() stops)
#   stop[i, k]  TRUE if design_ebased() stops for futility in that state
#   cost[i, t+1] optimal cost-to-go from (t, m_grid[i]) (see .solve_dp())

#' Power-maximising design-optimal e-value (Section 3.1.1).
design_pmax <- function(n, theta0, theta1, alpha = 0.05,
                        analysis_times = seq_len(n),
                        m_grid = baas_m_grid(alpha), b_grid = baas_b_grid()) {
  s <- .setup(n, theta0, theta1, alpha, analysis_times, m_grid, b_grid)
  .design_from_dp(s, "pmax", .solve_dp(s, "pmax"))
}

#' Expected-sample-size-minimising design-optimal e-value (Section 3.1.2).
design_essmin <- function(n, theta0, theta1, alpha = 0.05,
                          analysis_times = seq_len(n),
                          m_grid = baas_m_grid(alpha), b_grid = baas_b_grid()) {
  s <- .setup(n, theta0, theta1, alpha, analysis_times, m_grid, b_grid)
  .design_from_dp(s, "essmin", .solve_dp(s, "essmin"))
}

#' e-value-based design with optimal futility stopping and a power constraint
#' (Section 3.2, Appendix A.2).
#'
#' Solves min E[sample size] s.t. power >= 1 - beta under theta1 through the
#' Lagrangian E[sample size] + lambda * P(no rejection). Power is nondecreasing
#' in lambda, so a bisection finds the smallest lambda whose policy meets the
#' constraint; that policy has the smallest expected sample size among the
#' feasible Lagrangian policies. (The paper uses a discrete Newton method and
#' accepts any lambda with type II error in [beta - eps, beta].) A warning is
#' given if the attained type II error is more than `eps` below beta.
#'
#' Power is evaluated from starting e-value `m0`.
design_ebased <- function(n, theta0, theta1, alpha = 0.05, beta = 0.2,
                          analysis_times = seq_len(n),
                          m_grid = baas_m_grid(alpha), b_grid = baas_b_grid(),
                          eps = 0.01, m0 = 1, lambda_rel_tol = 1e-4,
                          verbose = FALSE) {
  s <- .setup(n, theta0, theta1, alpha, analysis_times, m_grid, b_grid)
  fit <- function(lambda) {
    d <- .design_from_dp(s, "ebased", .solve_dp(s, "ebased", lambda),
                         list(lambda = lambda, beta = beta))
    type2 <- 1 - .forward(d, theta1, m0)$reject
    if (verbose) cat(sprintf("lambda = %12.4f  type II error = %.5f\n", lambda, type2))
    list(design = d, type2 = type2)
  }

  lo <- 0
  hi <- fit(n)
  while (hi$type2 > beta) {
    if (hi$design$lambda > 1e9) {
      stop(sprintf(
        "Power %.4f cannot be reached; the largest attainable power on this grid is about %.4f (see design_pmax()).",
        1 - beta, 1 - hi$type2
      ))
    }
    lo <- hi$design$lambda
    hi <- fit(2 * lo)
  }
  while (hi$design$lambda - lo > lambda_rel_tol * hi$design$lambda) {
    mid <- fit((lo + hi$design$lambda) / 2)
    if (mid$type2 <= beta) hi <- mid else lo <- mid$design$lambda
  }
  if (beta - hi$type2 > eps) {
    warning(sprintf(
      "Type II error %.4f is more than eps = %g below beta = %g; power jumps in lambda here (a randomised policy would be needed to make the constraint tight).",
      hi$type2, eps, beta
    ))
  }
  hi$design
}

#' GROW e-value: constant Kelly bet (Section 2.3), on the same grid so it can
#' be evaluated with the same tools. `bet` overrides the Kelly bet.
design_grow <- function(n, theta0, theta1, alpha = 0.05,
                        analysis_times = seq_len(n),
                        m_grid = baas_m_grid(alpha),
                        bet = kelly_bet(theta0, theta1)) {
  s <- .setup(n, theta0, theta1, alpha, analysis_times, m_grid, b_grid = 0)
  nm <- length(m_grid)
  bets <- matrix(bet, nm, s$n)
  bets[c(1, nm), ] <- NA_real_
  .new_design(s, "grow", bets, matrix(FALSE, nm, s$n), cost = NULL)
}


# ---- Operating characteristics (forward recursion) ----------------------------

# Propagate the distribution of the grid e-value M-hat under success
# probability `theta`. At t = 0 and every analysis time, mass at 1/alpha is
# counted as a rejection; before n, mass is counted as a futility stop if the
# design stops there (design_ebased()) or, with `futility_hz`, if the state is
# in the hopeless zone (as in Section 4 of the paper).
.forward <- function(d, theta, m0 = 1, futility_hz = TRUE, keep_occupancy = FALSE) {
  g <- d$m_grid
  nm <- length(g)
  top <- nm
  n <- d$n
  hz_bound <- d$theta0^(n - (0:n)) / d$alpha

  p <- numeric(nm)
  p[findInterval(m0, g)] <- 1
  cum_rej <- cum_fut <- p_active <- numeric(n + 1)
  rej <- fut <- ess <- 0
  occ <- if (keep_occupancy) matrix(0, nm, n + 1) else NULL

  for (t in 0:n) {
    if (keep_occupancy) occ[, t + 1] <- p
    if (t == n || d$is_epoch[t + 1]) {
      rej <- rej + p[top]
      p[top] <- 0
      if (t < n) {
        st <- d$stop[, t + 1]
        if (futility_hz) st <- st | g < hz_bound[t + 1]
        fut <- fut + sum(p[st])
        p[st] <- 0
        ess <- ess + sum(p) * d$block_after[t + 1]
      }
    }
    cum_rej[t + 1] <- rej
    cum_fut[t + 1] <- fut
    p_active[t + 1] <- sum(p)
    if (t < n) {
      b <- d$bet[, t + 1]
      b[is.na(b)] <- 0
      up <- findInterval(g * (1 + b * (1 / d$theta0 - 1)), g)
      down <- findInterval(g * (1 - b), g)
      up[top] <- top
      down[top] <- top
      p <- .accumulate(c(up, down), c(theta * p, (1 - theta) * p), nm)
    }
  }
  list(
    reject = rej, futility = fut, ess = ess,
    by_time = data.frame(t = 0:n, cum_reject = cum_rej, cum_futility = cum_fut,
                         p_active = p_active),
    occupancy = occ
  )
}

#' Exact operating characteristics of a design on the e-value grid.
#'
#' Returns P(reject H0) (power, or type I error if theta <= theta0), P(futility
#' stop before n), expected sample size, and cumulative curves by participant
#' (as in Figure 1). Rejection and futility only happen at analysis times.
#' Futility stops are design_ebased() stop decisions and, if `futility_hz`,
#' entries into the hopeless zone. `occupancy[i, t + 1]` = P(M-hat_t = m_grid[i],
#' trial still running).
evaluate_design <- function(design, theta, m0 = 1, futility_hz = TRUE) {
  res <- .forward(design, theta, m0, futility_hz, keep_occupancy = TRUE)
  c(list(theta = theta), res)
}


# ---- Using a design on data ---------------------------------------------------

#' Bet for participant t + 1 when the grid (shadow) e-value after t
#' participants is m_hat. NA means: absorbing state or stop.
next_bet <- function(design, t, m_hat) {
  stopifnot(t >= 0, t < design$n)
  design$bet[findInterval(m_hat, design$m_grid), t + 1]
}

#' Run a design on observed outcomes y (0/1, in arrival order).
#'
#' Columns: bet B_t and lambda_t = B_t / theta0 for participant t; `e_value` is
#' the real e-value M-tilde_t; `e_value_grid` is the shadow M-hat_t used to look
#' up the next bet; `rejected` = M-tilde has reached 1/alpha at or before t (in a
#' blocked design, H0 is rejected at the next analysis time); `hopeless` =
#' M-tilde_t in the hopeless zone; `stop_advised` = design_ebased() advises not
#' recruiting participant t (non-binding; the bet is then 0).
e_process <- function(design, y, m0 = 1) {
  y <- as.integer(y)
  if (!all(y %in% 0:1)) stop("`y` must be 0/1")
  if (length(y) > design$n) stop("more outcomes than the design horizon n")
  g <- design$m_grid
  top <- length(g)
  thr <- 1 / design$alpha
  i <- findInterval(m0, g)
  m <- m0
  rejected <- m >= thr
  k <- length(y)
  out <- data.frame(t = seq_len(k), y = y, bet = NA_real_, lambda = NA_real_,
                    e_value = NA_real_, e_value_grid = NA_real_,
                    rejected = NA, hopeless = NA, stop_advised = NA)
  for (t in seq_len(k)) {
    stop_adv <- design$stop[i, t]
    b <- if (rejected) 0 else design$bet[i, t]
    if (is.na(b)) b <- 0
    f <- if (y[t] == 1) 1 + b * (1 / design$theta0 - 1) else 1 - b
    m <- m * f
    i <- if (i == top) top else findInterval(g[i] * f, g)
    rejected <- rejected || m >= thr
    out$bet[t] <- b
    out$e_value[t] <- m
    out$e_value_grid[t] <- g[i]
    out$rejected[t] <- rejected
    out$hopeless[t] <- !rejected && hopeless(m, t, design$n, design$theta0, design$alpha)
    out$stop_advised[t] <- stop_adv
  }
  out$lambda <- out$bet / design$theta0
  out
}

#' Monte Carlo simulation of the real e-process M-tilde (bets looked up at the
#' shadow M-hat). H0 is rejected at the first analysis time at which M-tilde
#' has reached 1/alpha (betting stops once it has). Futility stops follow the
#' same rules as evaluate_design(), applied to M-hat. `min_gap` is the
#' smallest M-tilde - M-hat over all paths and times (it is >= 0, Appendix A.1).
simulate_design <- function(design, theta, nsim = 1e4, m0 = 1, futility_hz = TRUE) {
  g <- design$m_grid
  top <- length(g)
  n <- design$n
  thr <- 1 / design$alpha
  idx <- rep(findInterval(m0, g), nsim)
  m <- rep(m0, nsim)
  hit <- m >= thr
  active <- rep(TRUE, nsim)
  rejected <- futile <- rep(FALSE, nsim)
  size <- numeric(nsim)
  min_gap <- Inf
  for (t in 0:n) {
    if (t == n || design$is_epoch[t + 1]) {
      r <- active & hit
      rejected[r] <- TRUE
      active[r] <- FALSE
      if (t < n) {
        st <- active & design$stop[cbind(idx, t + 1)]
        if (futility_hz) st <- st | (active & hopeless(g[idx], t, n, design$theta0, design$alpha))
        futile[st] <- TRUE
        active[st] <- FALSE
        size[active] <- size[active] + design$block_after[t + 1]
      }
    }
    if (t < n) {
      y <- stats::rbinom(nsim, 1, theta)
      b <- design$bet[cbind(idx, t + 1)]
      b[is.na(b) | hit | !active] <- 0
      f <- ifelse(y == 1, 1 + b * (1 / design$theta0 - 1), 1 - b)
      m <- m * f
      moved <- idx != top
      idx[moved] <- findInterval(g[idx[moved]] * f[moved], g)
      hit <- hit | m >= thr
      min_gap <- min(min_gap, m - g[idx])
    }
  }
  list(theta = theta, nsim = nsim, reject = mean(rejected),
       futility = mean(futile), ess = mean(size), min_gap = min_gap)
}


# ---- Printing and plotting -----------------------------------------------------

print.evalue_design <- function(x, ...) {
  labels <- c(pmax = "power-maximising (P-max) e-value",
              essmin = "expected-sample-size-minimising (ESS-min) e-value",
              ebased = "e-value-based design with futility stopping",
              grow = "GROW e-value (constant Kelly bet)")
  cat("Design-optimal e-value design:", labels[[x$type]], "\n")
  cat(sprintf("  n = %d, theta0 = %g, theta1 = %g, alpha = %g\n",
              x$n, x$theta0, x$theta1, x$alpha))
  if (length(x$analysis_times) < x$n) {
    cat("  analysis times:", x$analysis_times, "\n")
  }
  if (x$type == "ebased") {
    cat(sprintf("  beta = %g, Lagrange multiplier lambda = %.4f\n", x$beta, x$lambda))
  }
  if (x$type == "grow") {
    cat(sprintf("  grid: %d e-values\n", length(x$m_grid)))
  } else {
    cat(sprintf("  grid: %d e-values, %d bets\n", length(x$m_grid), length(x$b_grid)))
  }
  b1 <- next_bet(x, 0, 1)
  cat(sprintf("  first bet B_1 = %s (Kelly bet = %.4f)\n",
              format(b1, digits = 4), kelly_bet(x$theta0, x$theta1)))
  oc1 <- .forward(x, x$theta1)
  oc0 <- .forward(x, x$theta0)
  cat(sprintf("  under theta1: power = %.4f, ESS = %.2f, P(futility) = %.4f\n",
              oc1$reject, oc1$ess, oc1$futility))
  cat(sprintf("  under theta0: type I error = %.4f, ESS = %.2f, P(futility) = %.4f\n",
              oc0$reject, oc0$ess, oc0$futility))
  invisible(x)
}

.bet_palette <- function(kelly, n_col = 200) {
  breaks <- seq(0, 1, length.out = n_col + 1)
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  low <- grDevices::colorRamp(c("#08306b", "#4292c6", "#ffffff"))
  high <- grDevices::colorRamp(c("#ffffff", "#ef3b2c", "#99000d", "#000000"))
  col <- character(n_col)
  lo <- mid <= kelly
  col[lo] <- grDevices::rgb(low(mid[lo] / kelly), maxColorValue = 255)
  col[!lo] <- grDevices::rgb(high((mid[!lo] - kelly) / (1 - kelly)), maxColorValue = 255)
  list(breaks = breaks, col = col)
}

#' Heat map of a betting strategy, as in Figures 2-3 of the paper. The y-axis
#' is the grid index (log scale below 1, linear from 1 to 1/alpha). Colour:
#' blue = bet below the Kelly bet, white = Kelly bet, red to black = up to 1.
#' Grey: hopeless zone (Remark 3) and futility stops. Faded: states that cannot
#' be reached from e-value m0.
plot_policy <- function(design, main = NULL, m_min = 1e-3, m0 = 1,
                        kelly = kelly_bet(design$theta0, design$theta1)) {
  g <- design$m_grid
  n <- design$n
  rows <- which(g >= m_min & g < 1 / design$alpha)
  y <- seq_along(rows)
  pal <- .bet_palette(kelly)
  graphics::image(seq_len(n), y, t(design$bet[rows, , drop = FALSE]),
                  col = pal$col, breaks = pal$breaks, axes = FALSE,
                  xlab = "Participant index", ylab = "e-value",
                  main = if (is.null(main)) design$type else main)

  hz <- outer(g[rows], 0:(n - 1), hopeless, n = n, theta0 = design$theta0,
              alpha = design$alpha)
  grey <- hz | design$stop[rows, , drop = FALSE]
  graphics::image(seq_len(n), y, t(ifelse(grey, 1, NA)), col = "grey65", add = TRUE)

  occ <- .forward(design, 0.5, m0, futility_hz = FALSE, keep_occupancy = TRUE)$occupancy
  unreachable <- occ[rows, seq_len(n), drop = FALSE] == 0 & !grey
  graphics::image(seq_len(n), y, t(ifelse(unreachable, 1, NA)),
                  col = grDevices::adjustcolor("white", alpha.f = 0.6), add = TRUE)

  graphics::axis(1)
  ticks <- c(10^(floor(log10(m_min)):0), 1 / design$alpha)
  at <- vapply(ticks, function(v) which.min(abs(g[rows] - v)), 1L)
  graphics::axis(2, at = at, labels = format(ticks, scientific = FALSE, drop0trailing = TRUE),
                 las = 1)
  graphics::box()
  invisible(NULL)
}
