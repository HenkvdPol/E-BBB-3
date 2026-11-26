# File with additional functions
#
# The e-value test of framework.tex (Additional_files/framework.pdf) with Bayesian borrowing.
#   Part 1. The e-value: priors on (p0, 1], block factors, futility.
#   Part 2. Borrowing models: none, one-way EXNEX and two-way EXNEX.
#   Part 3. Simulation of the trial in calendar time.
#   Part 4. Summaries and figures.
#
# The type I error guarantee only needs Part 1: every prior must lie on (p0, 1], sum to 1, be
# chosen from the data before the block and stay fixed during the block. Part 2 only decides
# which prior that is, so it changes the power, never the validity.

# Necessary libraries (ggplot2 is loaded in the figure functions only)
library(Matrix)


################################################################################
################################################################################
# Part 1. The e-value
################################################################################
################################################################################

################################################################################
# Grid of theta = logit(p) on which all priors and posteriors live
################################################################################
theta_grid <- function(p0, h = 0.05, range = c(-10, 6)){
  # p0 = null response rate. logit(p0) is put on a cell boundary, so that every
  #      cell lies entirely below or entirely above p0.
  # h = width of the cells on the logit scale
  # range = smallest and largest theta

  # Cell numbers g: the cell (logit(p0) + (g - 1) h, logit(p0) + g h)
  g <- seq(ceiling((range[1] - qlogis(p0)) / h + 0.5),
           floor((range[2] - qlogis(p0)) / h + 0.5))
  theta <- qlogis(p0) + (g - 0.5) * h
  p <- plogis(theta)

  # The cells above p0 form the support of the priors of the e-value.
  above <- g >= 1

  grid <- list(theta = theta, lower = theta - h / 2, upper = theta + h / 2,
               p = p, h = h, p0 = p0, range = range, above = above,
               support = p[above])

  # Log likelihood ratio against p0 of one response (l1) and one non-response (l0)
  grid$l1 <- log(p[above] / p0)
  grid$l0 <- log((1 - p[above]) / (1 - p0))

  return(grid)
}

################################################################################
# Probability of every grid cell under a normal distribution
################################################################################
norm_cells <- function(grid, mean, sd){
  # grid = theta grid (theta_grid)
  # mean, sd = mean(s) and standard deviation(s) of the normal distribution(s)
  # Returns a (cells x distributions) matrix; every column sums to 1 over the grid.
  G <- length(grid$theta)
  M <- max(length(mean), length(sd))
  mean <- rep(rep_len(mean, M), each = G)
  sd <- rep(rep_len(sd, M), each = G)

  # Standardised cell boundaries
  a <- (rep(grid$lower, M) - mean) / sd
  b <- (rep(grid$upper, M) - mean) / sd

  # For cells in the upper tail we use upper-tail probabilities (more accurate).
  P <- ifelse(a > 0,
              pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE),
              pnorm(b) - pnorm(a))
  P <- matrix(P, G, M)

  return(sweep(P, 2, colSums(P), "/"))
}

################################################################################
# Largest value in every row of a matrix
################################################################################
row_max <- function(x){
  return(x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))])
}

################################################################################
# Check that priors meet the conditions of the e-value
################################################################################
check_priors <- function(w, grid, tol = 1e-9){
  # w = matrix of prior weights, one row per basket, one column per support point
  # grid = theta grid
  # The only requirements of the type I error guarantee: support in (p0, 1] and
  # finite weights >= 0 that sum to 1.
  if(!all(grid$support > grid$p0 & grid$support <= 1)){
    stop("prior support must lie in (p0, 1]")
  }
  if(!is.matrix(w) || ncol(w) != length(grid$support)){
    stop("prior weights must be a matrix with one column per support point")
  }
  if(any(!is.finite(w)) || any(w < 0)){
    stop("prior weights must be finite and >= 0")
  }
  if(any(abs(rowSums(w) - 1) > tol)){
    stop("prior weights must sum to 1")
  }
  invisible(TRUE)
}

################################################################################
# The prior of framework.tex eq. (8): the posterior restricted to (p0, 1]
################################################################################
restrict <- function(logpost, grid){
  # logpost = (baskets x cells) log posterior on the grid, up to a constant per row
  # grid = theta grid
  lw <- logpost[, grid$above, drop = FALSE]

  # Normalise the weights on the cells above p0
  w <- exp(lw - row_max(lw))
  w <- w / rowSums(w)

  check_priors(w, grid)
  return(w)
}

################################################################################
# Log block factor of the e-value
################################################################################
log_factor <- function(w, y, n, grid){
  # w = prior weights, one row per block
  # y, n = responses and patients in the block (n = 0 gives factor 1)
  # grid = theta grid, or the futility side (futility_side) for the futility e-value
  # The factor is sum_g w_g (p_g / p0)^y ((1 - p_g) / (1 - p0))^(n - y).
  lr <- outer(y, grid$l1) + outer(n - y, grid$l0)
  m <- row_max(lr)
  return(m + log(rowSums(w * exp(lr - m))))
}

################################################################################
# Scaled Beta with the same mean and variance as the restricted posterior
################################################################################
scaled_beta <- function(w, grid){
  # w = restricted posterior weights (one row per basket)
  # grid = theta grid
  # Benchmark (decision D4): p = p0 + (1 - p0) u with u ~ Beta(a, b), put on the same
  # cells. It is another prior on (p0, 1], so the validity is unaffected.
  p0 <- grid$p0

  # Mean and variance of u = (p - p0) / (1 - p0)
  m <- as.vector(w %*% grid$support)
  s2 <- as.vector(w %*% grid$support^2)
  u <- (m - p0) / (1 - p0)
  uv <- (s2 - m^2) / (1 - p0)^2
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  uv <- pmin(pmax(uv, 1e-12), u * (1 - u) * (1 - 1e-6))

  # Beta parameters with this mean and variance
  sz <- u * (1 - u) / uv - 1
  a <- u * sz
  b <- (1 - u) * sz

  # Probability of every cell: differences of the Beta distribution function
  ub <- c((plogis(grid$lower[grid$above]) - p0) / (1 - p0), 1)
  out <- t(vapply(seq_along(a), function(i){diff(pbeta(ub, a[i], b[i]))},
                  numeric(length(grid$support))))
  out <- pmax(out, 0)
  out <- out / rowSums(out)

  check_priors(out, grid)
  return(out)
}

################################################################################
# Futility: support of the futility e-value (testing p >= p1)
################################################################################
futility_side <- function(grid, p1){
  # grid = theta grid
  # p1 = reference alternative. The futility e-value tests H1: P(Y = 1 | past) >= p1
  #      with priors on [0, p1); stop for futility when it reaches 1/beta. Then
  #      P(futility stop) <= beta whenever H1 holds (framework.tex, Proposition 2).
  # The grid is aligned at logit(p0), so the support is the cell midpoints below p1.
  if(p1 <= grid$p0){
    stop("futility needs p1 > p0")
  }
  keep <- grid$p < p1
  q <- grid$p[keep]

  return(list(p1 = p1, keep = keep, support = q,
              l1 = log(q / p1), l0 = log((1 - q) / (1 - p1))))
}

################################################################################
# Futility: check the futility priors
################################################################################
check_priors_below <- function(w, fut, tol = 1e-9){
  # w = futility prior weights; fut = futility side (futility_side)
  if(!all(fut$support >= 0 & fut$support < fut$p1)){
    stop("futility prior support must lie in [0, p1)")
  }
  if(!is.matrix(w) || ncol(w) != length(fut$support)){
    stop("futility prior weights must have one column per support point")
  }
  if(any(!is.finite(w)) || any(w < 0)){
    stop("futility prior weights must be finite and >= 0")
  }
  if(any(abs(rowSums(w) - 1) > tol)){
    stop("futility prior weights must sum to 1")
  }
  invisible(TRUE)
}

################################################################################
# Futility: the posterior restricted to [0, p1)
################################################################################
restrict_below <- function(logpost, fut){
  # logpost = (baskets x cells) log posterior on the grid
  # fut = futility side. The futility block factor is log_factor(w, y, n, fut).
  lw <- logpost[, fut$keep, drop = FALSE]
  w <- exp(lw - row_max(lw))
  w <- w / rowSums(w)

  check_priors_below(w, fut)
  return(w)
}


################################################################################
################################################################################
# Part 2. Borrowing models
#
# Every model ("engine") is a list with a function posterior(Y, N, need, seed):
#   Y, N = (data states x baskets) responses and patients;
#   need = (data states x baskets) TRUE where the posterior is needed.
# It returns logpost[state, basket, cell], the log posterior of theta given ALL data
# of the state (NA where need is FALSE), and pex = P(basket exchangeable | data).
################################################################################
################################################################################

################################################################################
# Log likelihood of data states on the grid
################################################################################
loglik_grid <- function(y, n, grid){
  # y, n = responses and patients (vectors); returns (states x cells), without
  # the binomial coefficient.
  return(outer(y, log(grid$p)) + outer(n - y, log(1 - grid$p)))
}

################################################################################
# NEX prior: theta ~ N(logit pbar, 1/pbar + 1/(1 - pbar)) on the grid
################################################################################
nex_prior <- function(grid, pbar){
  # pbar = (p0 + p1) / 2: unit-information prior centred between p0 and p1
  return(as.vector(norm_cells(grid, qlogis(pbar), sqrt(1 / pbar + 1 / (1 - pbar)))))
}

################################################################################
# Table with one row per data state (n, y) of a basket, growing when needed
################################################################################
state_table <- function(compute, n_max){
  # compute = function(y, n) returning a named list of matrices (one row per state)
  #           or vectors (one element per state); they become fields of the table.
  # n_max = rows for all n <= n_max are built at once. Rows for larger n are added
  #         the first time they are needed, so there is no maximum sample size.
  # tab$row_of(y, n) returns the row of every state. Call it BEFORE reading a field:
  # in tab$LE[tab$row_of(y, n), ] R reads the old LE before the new rows are added.
  tab <- new.env()

  # Unique number of the state (n, y)
  key <- function(y, n){
    return(n * (n + 1) / 2 + y + 1)
  }

  # All states with n <= n_max
  nn <- rep(0:n_max, times = 1:(n_max + 1))
  yy <- sequence(1:(n_max + 1)) - 1
  first <- compute(yy, nn)
  for(nm in names(first)){
    tab[[nm]] <- first[[nm]]
  }
  tab$fields <- names(first)
  tab$pos <- seq_along(yy)   # key -> row

  tab$row_of <- function(y, n){
    if(any(y < 0 | y > n)){
      stop("impossible data state (y > n or y < 0)")
    }
    k <- key(as.vector(y), as.vector(n))
    if(length(k) && max(k) > length(tab$pos)){
      tab$pos[max(k)] <- NA
    }

    # States that are not in the table yet: compute and add them.
    new <- which(is.na(tab$pos[k]) & !duplicated(k))
    if(length(new)){
      add <- compute(as.vector(y)[new], as.vector(n)[new])
      r0 <- NROW(tab[[tab$fields[1]]])
      for(nm in tab$fields){
        if(is.matrix(tab[[nm]])){
          tab[[nm]] <- rbind(tab[[nm]], add[[nm]])
        }else{
          tab[[nm]] <- c(tab[[nm]], add[[nm]])
        }
      }
      tab$pos[k[new]] <- r0 + seq_along(new)
    }
    return(tab$pos[k])
  }

  return(tab)
}

################################################################################
# Borrowing model 'none': the NEX prior and the basket's own data only
################################################################################
engine_none <- function(grid, pbar){
  # grid = theta grid; pbar = centre of the NEX prior
  logV <- log(nex_prior(grid, pbar))

  posterior <- function(Y, N, need = NULL, seed = NULL){
    lp <- array(NA_real_, c(nrow(Y), ncol(Y), length(grid$theta)))
    for(k in seq_len(ncol(Y))){
      lp[, k, ] <- loglik_grid(Y[, k], N[, k], grid) + rep(logV, each = nrow(Y))
    }
    return(list(logpost = lp, pex = matrix(NA_real_, nrow(Y), ncol(Y))))
  }

  return(list(name = "none", grid = grid, posterior = posterior))
}

################################################################################
# One-way EXNEX: grid of the hyperparameters (mu, tau) and their prior
################################################################################
oneway_hyper <- function(p0, mu_sd = 10, tau_scale = 1){
  # The benchmark of Additional_files/legacy: mu ~ N(logit p0, 10^2) on [-8, 4], tau ~ HN(1)
  h <- expand.grid(mu = seq(-8, 4, by = 0.15), tau = seq(0.05, 3.95, by = 0.10))
  h$logprior <- dnorm(h$mu, qlogis(p0), mu_sd, log = TRUE) +
    dnorm(h$tau, 0, tau_scale, log = TRUE)
  return(h)
}

################################################################################
# Borrowing model 'one-way EXNEX' (Neuenschwander et al. 2016)
################################################################################
engine_oneway <- function(grid, pbar, n_max, w_ex = 0.5, hyper = oneway_hyper(grid$p0),
                          kernel = "cell"){
  # grid = theta grid; pbar = centre of the NEX prior
  # n_max = the (n, y) table is built up to n_max at once and grows beyond it
  # w_ex = prior probability that a basket is exchangeable (EX)
  # hyper = grid of (mu, tau) with log prior
  # kernel = "cell": EX prior as cell probabilities (default); "density": normal
  #          density at the cell midpoints (the legacy choice; regression test only)
  # Model: EX theta_c ~ N(mu, tau^2), NEX theta_c ~ N(logit pbar, 1/pbar + 1/(1 - pbar)).
  # Given (mu, tau) the baskets are independent, so the posterior is a sum over the
  # (mu, tau) grid of one-dimensional sums on the theta grid (no sampling).

  # EX prior of theta for every (mu, tau): (cells x hyperparameters)
  if(kernel == "cell"){
    W <- norm_cells(grid, hyper$mu, hyper$tau)
    V <- nex_prior(grid, pbar)
  }else{
    D <- dnorm(outer(grid$theta, hyper$mu, "-") / rep(hyper$tau, each = length(grid$theta)))
    W <- sweep(D, 2, colSums(D), "/")
    v <- dnorm(grid$theta, qlogis(pbar), sqrt(1 / pbar + 1 / (1 - pbar)))
    V <- v / sum(v)
  }

  # Marginal likelihood of every data state under EX (per (mu, tau)) and under NEX
  tab <- state_table(function(y, n){
    Lik <- exp(loglik_grid(y, n, grid))
    return(list(LE = Lik %*% W, LN = as.vector(Lik %*% V)))
  }, n_max)

  posterior <- function(Y, N, need = NULL, seed = NULL){
    S <- nrow(Y)
    K <- ncol(Y)
    G <- length(grid$theta)
    if(is.null(need)){
      need <- matrix(TRUE, S, K)
    }
    r <- matrix(tab$row_of(Y, N), S, K)
    LE <- tab$LE
    LN <- tab$LN

    # Likelihood of basket k given (mu, tau): D_kh = w LE + (1 - w) LN
    D <- lapply(seq_len(K), function(k){
      w_ex * LE[r[, k], , drop = FALSE] + (1 - w_ex) * LN[r[, k]]
    })

    # Posterior of (mu, tau) given all baskets
    lw <- matrix(hyper$logprior, S, nrow(hyper), byrow = TRUE) + Reduce(`+`, lapply(D, log))
    w <- exp(lw - row_max(lw))
    w <- w / rowSums(w)

    lp <- array(NA_real_, c(S, K, G))
    pex <- matrix(NA_real_, S, K)
    for(k in seq_len(K)){
      i <- which(need[, k])
      if(!length(i)){
        next
      }
      # Prior of theta_k given the other baskets: (mu, tau) weights w_h / D_kh
      u <- w[i, , drop = FALSE] / D[[k]][i, , drop = FALSE]
      mix <- w_ex * (u %*% t(W)) + (1 - w_ex) * outer(rowSums(u), V)

      # Posterior = own likelihood x that prior
      lp[i, k, ] <- loglik_grid(Y[i, k], N[i, k], grid) + log(mix)
      pex[i, k] <- w_ex * rowSums(u * LE[r[i, k], , drop = FALSE])
    }
    return(list(logpost = lp, pex = pex))
  }

  return(list(name = "one-way EXNEX", grid = grid, W = W, V = V, hyper = hyper, tab = tab,
              posterior = posterior))
}

################################################################################
# Two-way EXNEX: layout of the baskets in the J x L grid
################################################################################
grid_layout <- function(J, L, missing = NULL){
  # J = number of tumour types (rows); L = number of biomarker groups (columns)
  # missing = list of cells c(j, l) that are not in the trial
  lay <- expand.grid(row = seq_len(J), col = seq_len(L))
  lay <- lay[order(lay$row, lay$col), ]
  lay <- lay[!paste(lay$row, lay$col) %in% vapply(missing, paste, "", collapse = " "), ]
  lay$basket <- sprintf("T%d.B%d", lay$row, lay$col)
  rownames(lay) <- NULL
  return(lay)
}

################################################################################
# Two-way EXNEX: half-normal prior of a standard deviation on a grid
################################################################################
hn_grid <- function(scale, sigma = seq(0.05, 2.45, by = 0.1)){
  lp <- dnorm(sigma, 0, scale, log = TRUE)
  return(list(value = sigma, logprior = lp - log(sum(exp(lp)))))
}

################################################################################
# Two-way EXNEX: the prior (framework.tex eq. 9)
################################################################################
twoway_prior <- function(p0, pbar, w_ex = 0.5, mu_sd = 10, mu_range = c(-8, 4), scale = 1 / sqrt(3),
                         sig_a = hn_grid(scale), sig_b = hn_grid(scale), sig_e = hn_grid(scale)){
  # Basket c = (j, l), tumour type j and biomarker group l:
  #   delta_c ~ Bernoulli(w_ex), independent across baskets
  #   EX  (delta_c = 1): theta_c ~ N(mu + a_j + b_l, sigma_e^2)
  #   NEX (delta_c = 0): theta_c ~ N(logit pbar, 1/pbar + 1/(1 - pbar))
  #   mu ~ N(logit p0, mu_sd^2) on mu_range; a_j ~ N(0, sigma_a^2); b_l ~ N(0, sigma_b^2)
  #   sigma_a, sigma_b, sigma_e ~ HN(scale) on the grid 0.05, 0.15, ..., 2.45
  # scale = 1/sqrt(3): E(sigma_a^2 + sigma_b^2 + sigma_e^2) = 1 = E(tau^2) of the
  # one-way benchmark with tau ~ HN(1) (decision D2). No sum-to-zero constraint (D5).
  return(list(w_ex = w_ex, mu_mean = qlogis(p0), mu_sd = mu_sd, mu_range = mu_range,
              scale = c(a = scale, b = scale, e = scale),
              sig_a = sig_a, sig_b = sig_b, sig_e = sig_e, pbar = pbar))
}

################################################################################
# Two-way EXNEX: fixed standard deviations (for the reference comparisons)
################################################################################
fix_sigmas <- function(prior, a, b, e){
  point <- function(v){
    return(list(value = v, logprior = 0))
  }
  return(modifyList(prior, list(sig_a = point(a), sig_b = point(b), sig_e = point(e),
                                fixed = c(a = a, b = b, e = e))))
}

################################################################################
# Two-way EXNEX: tables of the one-dimensional integrals
################################################################################
twoway_tables <- function(grid, prior, n_max, m_step = 0.05, rb_step = 0.1){
  # Given the EX mean m = mu + a_j + b_l and sigma_e, theta_c and delta_c only depend
  # on the basket's own data, so both are summed out on the theta grid:
  #   D_c(m, sigma_e) = w sum_theta L_c N(theta | m, sigma_e) + (1 - w) sum_theta L_c V.
  # Tables (one row per data state (n, y), growing when needed):
  #   logD  = log D on the m lattice (step m_step), for the sampler;
  #   DEXrb = EX part on the coarser Rao-Blackwell lattice (step rb_step);
  # and KexT = EX prior of theta given (m, sigma_e) on the Rao-Blackwell lattice.
  m <- seq(grid$range[1], grid$range[2], by = m_step)
  Gm <- length(m)
  rb <- seq(1, Gm, by = round(rb_step / m_step))
  Grb <- length(rb)
  se <- prior$sig_e$value
  Gs <- length(se)

  # EX prior of theta for every (sigma_e, m): cells x (Gs * Gm)
  P <- do.call(cbind, lapply(se, function(s){norm_cells(grid, m, s)}))
  V <- nex_prior(grid, prior$pbar)

  # Columns of P that belong to the Rao-Blackwell lattice, in the order (sigma_e, m)
  rbcol <- as.vector(outer(rb, (seq_len(Gs) - 1) * Gm, "+"))

  tab <- state_table(function(y, n){
    Lik <- exp(loglik_grid(y, n, grid))
    DEX <- Lik %*% P
    DN <- as.vector(Lik %*% V)
    return(list(logD = log(prior$w_ex * DEX + (1 - prior$w_ex) * DN),
                DEXrb = DEX[, rbcol, drop = FALSE], DN = DN))
  }, n_max)

  fixed <- list(m = m, m_step = m_step, Gm = Gm, rb_m = m[rb], rb_step = rb_step, Grb = Grb,
                Gs = Gs, n_max = n_max, KexT = t(P[, rbcol, drop = FALSE]), V = V)
  for(nm in names(fixed)){
    tab[[nm]] <- fixed[[nm]]
  }
  return(tab)
}

################################################################################
# Two-way EXNEX: log D_c(m, sigma_e) by linear interpolation in m
################################################################################
lookup_logD <- function(tab, rho, s, m){
  # rho = table rows (data states); s = index of sigma_e; m = EX means (clamped to the lattice)
  pos <- (pmin(pmax(m, tab$m[1]), tab$m[tab$Gm]) - tab$m[1]) / tab$m_step
  i0 <- pmin(floor(pos), tab$Gm - 2)
  f <- pos - i0
  col <- (s - 1) * tab$Gm + i0 + 1
  rho <- as.vector(rho)
  col <- as.vector(col)
  return((1 - f) * tab$logD[cbind(rho, col)] + f * tab$logD[cbind(rho, col + 1)])
}

################################################################################
# Two-way EXNEX: entries of a sparse matrix
################################################################################
triplets <- function(Q){
  t <- methods::as(Q, "TsparseMatrix")
  return(list(i = t@i + 1L, j = t@j + 1L, x = t@x))
}

################################################################################
# Two-way EXNEX: posterior of theta from the distribution of (m, sigma_e)
################################################################################
theta_from_q <- function(Q, y, n, tab, grid, w_ex){
  # Q = sparse (pairs x (Gs * Grb)) matrix: distribution of (m_c, sigma_e) of every
  #     pair (state, basket) on the Rao-Blackwell lattice; rows sum to 1
  # y, n = data of the basket of every pair
  # Rao-Blackwellisation: p(theta | y) = E[L {w N(theta | m, s) + (1 - w) V} / D(m, s)].
  tq <- triplets(Q)
  rho <- tab$row_of(y, n)[tq$i]
  dex <- tab$DEXrb[cbind(rho, tq$j)]
  d <- w_ex * dex + (1 - w_ex) * tab$DN[rho]

  U <- Matrix::sparseMatrix(tq$i, tq$j, x = tq$x / d, dims = dim(Q))
  mix <- w_ex * as.matrix(U %*% tab$KexT) + (1 - w_ex) * outer(Matrix::rowSums(U), tab$V)
  pex <- as.vector(tapply(tq$x * w_ex * dex / d, factor(tq$i, seq_len(nrow(Q))), sum))

  return(list(logpost = loglik_grid(y, n, grid) + log(mix), pex = pex))
}

################################################################################
# Two-way EXNEX: put draws of (m, sigma_e) on the Rao-Blackwell lattice
################################################################################
deposit <- function(pair, m, s, tab, n_pairs){
  # pair = pair number of every draw; m = draws of m; s = draws of the sigma_e index
  # Every draw is split over the two nearest lattice points (cloud-in-cell).
  pos <- (pmin(pmax(m, tab$rb_m[1]), tab$rb_m[tab$Grb]) - tab$rb_m[1]) / tab$rb_step
  i0 <- pmin(floor(pos), tab$Grb - 2)
  f <- pos - i0
  cell <- (s - 1) * tab$Grb + i0 + 1
  Q <- Matrix::sparseMatrix(c(pair, pair), c(cell, cell + 1), x = c(1 - f, f),
                            dims = c(n_pairs, tab$Gs * tab$Grb))
  return(Matrix::Diagonal(x = 1 / Matrix::rowSums(Q)) %*% Q)
}

################################################################################
# Two-way EXNEX: distribution of (m, sigma_e) without any data
################################################################################
empty_q <- function(tab, prior){
  # Prior predictive by deterministic quadrature: m = mu + a + b with mu ~ TN and
  # a + b ~ N(0, sigma_a^2 + sigma_b^2); mu on a 0.02 grid, m in bins of the lattice.
  mu <- seq(prior$mu_range[1], prior$mu_range[2], by = 0.02)
  wmu <- dnorm(mu, prior$mu_mean, prior$mu_sd)
  wmu <- wmu / sum(wmu)

  # Distinct variances of a + b and their prior probabilities
  v <- outer(prior$sig_a$value^2, prior$sig_b$value^2, "+")
  pv <- exp(outer(prior$sig_a$logprior, prior$sig_b$logprior, "+"))
  pv <- tapply(pv, round(v, 10), sum)
  v <- as.numeric(names(pv))

  # Probability of every bin of the lattice
  edges <- c(-Inf, tab$rb_m[-1] - tab$rb_step / 2, Inf)
  qm <- numeric(tab$Grb)
  for(i in seq_along(v)){
    qm <- qm + pv[[i]] * as.vector(diff(pnorm(outer(edges, mu, "-") / sqrt(v[i]))) %*% wmu)
  }

  # Column order (sigma_e, m)
  return(as.vector(outer(qm, exp(prior$sig_e$logprior))))
}

################################################################################
# Two-way EXNEX: Metropolis-within-Gibbs sampler of (mu, a, b, sigmas)
################################################################################
twoway_mcmc <- function(Y, N, layout, tab, prior, burn = 500, keep = 1000, chains = 2){
  # Y, N = (data states x baskets); layout = grid_layout; tab = twoway_tables
  # burn, keep, chains = burn-in, kept draws and chains per data state
  # All states and chains are sampled at once (vectorised). Chain r belongs to
  # state ((r - 1) %% S) + 1. Returns the draws of m_c (chains x baskets x keep)
  # and of the sigma_e index (chains x keep).
  S <- nrow(Y)
  K <- ncol(Y)
  J <- max(layout$row)
  L <- max(layout$col)
  Rc <- S * chains
  st <- rep(seq_len(S), chains)
  ch <- rep(seq_len(chains), each = S)
  rho <- matrix(tab$row_of(Y, N), S, K)[st, , drop = FALSE]

  # Which basket is in which row and column
  IR <- outer(layout$row, seq_len(J), "==")
  IC <- outer(layout$col, seq_len(L), "==")

  sa <- prior$sig_a
  sb <- prior$sig_b
  se <- prior$sig_e
  Gs <- length(se$value)
  lo <- prior$mu_range[1]
  hi <- prior$mu_range[2]
  m0 <- prior$mu_mean
  s0 <- prior$mu_sd

  # Log prior of mu (truncated normal)
  lp_mu <- function(x){
    return(ifelse(x >= lo & x <= hi, -(x - m0)^2 / (2 * s0^2), -Inf))
  }
  # Categorical draws from log probabilities (one per row)
  gumbel <- function(lp){
    return(max.col(lp - log(-log(matrix(runif(length(lp)), nrow(lp))))))
  }
  # Starting index of a standard deviation: near 0.4, three steps higher in chain 2
  start <- function(g, r){
    return(pmin(length(g$value), which.min(abs(g$value - 0.4)) + (r - 1) * 3))
  }

  # Starting values: pooled rate (chain 1) and one unit higher (chain 2), random effects 0
  pooled <- qlogis((rowSums(Y) + 0.5) / (rowSums(N) + 1))[st]
  mu <- pmin(pmax(pooled + (ch - 1), lo + 0.01), hi - 0.01)
  a <- matrix(0, Rc, J)
  b <- matrix(0, Rc, L)
  ia <- start(sa, ch)
  ib <- start(sb, ch)
  ie <- start(se, ch)
  m <- mu + a[, layout$row, drop = FALSE] + b[, layout$col, drop = FALSE]
  lD <- matrix(lookup_logD(tab, rho, rep(ie, K), m), Rc, K)

  step <- matrix(0.3, Rc, 3)   # random-walk steps of mu, a and b
  acc <- matrix(0, Rc, 3)      # acceptance counts
  M <- array(NA_real_, c(Rc, K, keep))
  E <- matrix(NA_integer_, Rc, keep)

  # Random-walk update of the row effects a or the column effects b: one proposal
  # per row (column); rows (columns) are independent given the rest.
  mh_block <- function(cur, prop, inc, idx, sig){
    mp <- m + (prop - cur)[, idx, drop = FALSE]
    lDp <- matrix(lookup_logD(tab, rho, rep(ie, K), mp), Rc, K)
    lr <- (lDp - lD) %*% inc + (cur^2 - prop^2) / (2 * sig^2)
    ok <- log(matrix(runif(length(lr)), Rc)) < lr
    okk <- ok[, idx, drop = FALSE]
    return(list(x = ifelse(ok, prop, cur), m = ifelse(okk, mp, m), lD = ifelse(okk, lDp, lD),
                rate = rowMeans(ok)))
  }

  for(it in seq_len(burn + keep)){
    # 1. mu
    mup <- mu + step[, 1] * rnorm(Rc)
    mp <- m + (mup - mu)
    lDp <- matrix(lookup_logD(tab, rho, rep(ie, K), mp), Rc, K)
    ok <- log(runif(Rc)) < rowSums(lDp - lD) + lp_mu(mup) - lp_mu(mu)
    mu[ok] <- mup[ok]
    m[ok, ] <- mp[ok, ]
    lD[ok, ] <- lDp[ok, ]
    acc[, 1] <- acc[, 1] + ok

    # 2. row effects a_j
    prop <- a + step[, 2] * matrix(rnorm(Rc * J), Rc)
    upd <- mh_block(a, prop, IR, layout$row, sa$value[ia])
    a <- upd$x
    m <- upd$m
    lD <- upd$lD
    acc[, 2] <- acc[, 2] + upd$rate

    # 3. column effects b_l
    prop <- b + step[, 3] * matrix(rnorm(Rc * L), Rc)
    upd <- mh_block(b, prop, IC, layout$col, sb$value[ib])
    b <- upd$x
    m <- upd$m
    lD <- upd$lD
    acc[, 3] <- acc[, 3] + upd$rate

    # 4. Exact Gibbs moves along (mu + k, a - k), (mu + k, b - k) and (a + k, b - k).
    #    They leave every m_c unchanged, so only the priors matter.
    va <- sa$value[ia]^2
    vb <- sb$value[ib]^2
    for(dir in 1:2){
      if(dir == 1){
        n_e <- J
        s_e <- rowSums(a)
        v_e <- va
      }else{
        n_e <- L
        s_e <- rowSums(b)
        v_e <- vb
      }
      P <- 1 / s0^2 + n_e / v_e
      mk <- ((m0 - mu) / s0^2 + s_e / v_e) / P
      sk <- 1 / sqrt(P)
      u <- runif(Rc, pnorm((lo - mu - mk) / sk), pnorm((hi - mu - mk) / sk))
      k <- mk + sk * qnorm(u)
      k[!is.finite(k)] <- 0
      mu <- mu + k
      if(dir == 1){
        a <- a - k
      }else{
        b <- b - k
      }
    }
    P <- J / va + L / vb
    k <- (rowSums(b) / vb - rowSums(a) / va) / P + rnorm(Rc) / sqrt(P)
    a <- a + k
    b <- b - k

    # 5. sigma_a and sigma_b: exact draws on their grids
    ia <- gumbel(matrix(sa$logprior, Rc, length(sa$value), byrow = TRUE) -
                   outer(rep(J, Rc), log(sa$value)) - outer(rowSums(a^2), 1 / (2 * sa$value^2)))
    ib <- gumbel(matrix(sb$logprior, Rc, length(sb$value), byrow = TRUE) -
                   outer(rep(L, Rc), log(sb$value)) - outer(rowSums(b^2), 1 / (2 * sb$value^2)))

    # 6. sigma_e: random walk on its grid
    if(Gs > 1){
      iep <- ie + sample(c(-2L, -1L, 1L, 2L), Rc, TRUE)
      ok <- iep >= 1 & iep <= Gs
      iep[!ok] <- ie[!ok]
      lDp <- matrix(lookup_logD(tab, rho, rep(iep, K), m), Rc, K)
      ok <- ok & log(runif(Rc)) < rowSums(lDp - lD) + se$logprior[iep] - se$logprior[ie]
      ie[ok] <- iep[ok]
      lD[ok, ] <- lDp[ok, ]
    }

    # Adapt the random-walk steps during burn-in (target acceptance 0.35)
    if(it <= burn && it %% 50 == 0){
      step <- pmin(pmax(step * exp(2 * (acc / 50 - 0.35)), 0.01), 5)
      acc[] <- 0
    }
    # Keep the draws after burn-in
    if(it > burn){
      M[, , it - burn] <- m
      E[, it - burn] <- ie
    }
  }

  return(list(M = M, E = E, state = st, chain = ch, step = step))
}

################################################################################
# Two-way EXNEX: convergence diagnostic (Gelman-Rubin R-hat of m_c)
################################################################################
rhat_m <- function(dr, S){
  # dr = output of twoway_mcmc; S = number of data states
  K <- dim(dr$M)[2]
  C <- max(dr$chain)
  out <- matrix(NA_real_, S, K)
  if(C < 2){
    return(out)
  }
  for(k in seq_len(K)){
    x <- dr$M[, k, ]
    n <- ncol(x)
    mc <- rowMeans(x)
    vc <- rowSums((x - mc)^2) / (n - 1)
    W <- rowMeans(matrix(vc, S))
    B <- n * apply(matrix(mc, S), 1, var)
    out[, k] <- sqrt(((n - 1) / n * W + B / n) / W)
  }
  return(out)
}

################################################################################
# Borrowing model 'two-way EXNEX' (framework.tex eq. 9)
################################################################################
engine_twoway <- function(grid, layout, pbar, n_max, prior = twoway_prior(grid$p0, pbar),
                          mcmc = list(burn = 500, keep = 1000, chains = 2, chunk = 150), tab = NULL){
  # grid = theta grid; layout = grid_layout; pbar = centre of the NEX prior
  # n_max = initial size of the (n, y) tables (they grow when needed)
  # mcmc = sampler settings; chunk = data states sampled together
  # Without any data the posterior is the prior predictive (quadrature, no sampling).
  # Otherwise (mu, a, b, sigmas) are sampled and theta is Rao-Blackwellised.
  if(is.null(tab)){
    tab <- twoway_tables(grid, prior, n_max)
  }
  q0 <- empty_q(tab, prior)

  posterior <- function(Y, N, need = NULL, seed = 1){
    S <- nrow(Y)
    K <- ncol(Y)
    G <- length(grid$theta)
    if(K != nrow(layout)){
      stop("number of baskets differs from the layout")
    }
    if(is.null(need)){
      need <- matrix(TRUE, S, K)
    }
    lp <- array(NA_real_, c(S, K, G))
    pex <- matrix(NA_real_, S, K)
    rhat <- matrix(NA_real_, S, K)

    # States without any patients: the prior predictive
    empty <- rowSums(N) == 0
    for(i in which(empty & rowSums(need) > 0)){
      k <- which(need[i, ])
      Q <- Matrix::Matrix(matrix(q0, length(k), length(q0), byrow = TRUE), sparse = TRUE)
      th <- theta_from_q(Q, Y[i, k], N[i, k], tab, grid, prior$w_ex)
      lp[i, k, ] <- th$logpost
      pex[i, k] <- th$pex
    }

    # States with patients: sample, in chunks of states
    todo <- which(!empty & rowSums(need) > 0)
    if(length(todo)){
      set.seed(seed)
    }
    for(ids in split(todo, ceiling(seq_along(todo) / mcmc$chunk))){
      dr <- twoway_mcmc(Y[ids, , drop = FALSE], N[ids, , drop = FALSE], layout, tab, prior,
                        mcmc$burn, mcmc$keep, mcmc$chains)
      rhat[ids, ] <- rhat_m(dr, length(ids))

      # Number the needed pairs (state, basket)
      nd <- need[ids, , drop = FALSE]
      pid <- matrix(0L, length(ids), K)
      pid[nd] <- seq_len(sum(nd))

      # Draws of (m, sigma_e) of every needed pair, from all chains
      sel <- nd[dr$state, , drop = FALSE]
      kk <- col(sel)[sel]
      rr <- row(sel)[sel]
      pair <- rep(pid[cbind(dr$state[rr], kk)], times = mcmc$keep)
      mm <- dr$M[cbind(rep(rr, times = mcmc$keep), rep(kk, times = mcmc$keep),
                       rep(seq_len(mcmc$keep), each = length(rr)))]
      ss <- dr$E[cbind(rep(rr, times = mcmc$keep), rep(seq_len(mcmc$keep), each = length(rr)))]
      Q <- deposit(pair, mm, ss, tab, sum(nd))

      # Posterior of theta of every needed pair
      w <- which(nd, arr.ind = TRUE)
      w <- w[order(pid[w]), , drop = FALSE]
      th <- theta_from_q(Q, Y[ids, , drop = FALSE][w], N[ids, , drop = FALSE][w], tab, grid, prior$w_ex)
      for(p in seq_len(nrow(w))){
        lp[ids[w[p, 1]], w[p, 2], ] <- th$logpost[p, ]
      }
      pex[cbind(ids[w[, 1]], w[, 2])] <- th$pex
    }
    return(list(logpost = lp, pex = pex, rhat = rhat))
  }

  return(list(name = "two-way EXNEX", grid = grid, layout = layout, tab = tab, prior = prior,
              posterior = posterior))
}


################################################################################
################################################################################
# Part 3. Simulation of the trial in calendar time
################################################################################
################################################################################

################################################################################
# Number of patients of every basket at time t
################################################################################
planned_n <- function(design, t){
  return(pmin(design$n_sim, floor(design$rate * pmax(0, t - design$entry) + 1e-9)))
}

################################################################################
# Complete and check a design
################################################################################
complete_design <- function(design){
  # design = list with
  #   layout = grid_layout; n_max = maximum per basket (Inf: open-ended, no maximum)
  #   rate, entry = patients per time unit and opening time per basket
  #   look_every = time between analyses; alpha = type I error per basket
  #   prior_form = "grid" (restricted posterior) or "beta" (scaled Beta benchmark)
  #   futility = "unreachable" (default): stop when 1/alpha cannot be reached by n_max
  #                                       even if all remaining patients respond;
  #              "evalue": stop when the futility e-value reaches 1/beta (needs p1, beta);
  #              "none"
  #   n_cap = simulation horizon, required when n_max = Inf. It only ends the
  #           computation: baskets still recruiting there are reported as 'open'.
  K <- nrow(design$layout)
  for(v in c("n_max", "rate", "entry")){
    design[[v]] <- rep_len(design[[v]], K)
  }
  if(is.null(design$prior_form)){
    design$prior_form <- "grid"
  }
  if(is.null(design$futility)){
    design$futility <- "unreachable"
  }
  if(!design$futility %in% c("unreachable", "evalue", "none")){
    stop("futility must be 'unreachable', 'evalue' or 'none'")
  }
  if(design$futility == "evalue" && (is.null(design$p1) || is.null(design$beta))){
    stop("e-value futility needs p1 and beta")
  }
  if(any(is.infinite(design$n_max)) && is.null(design$n_cap)){
    stop("open-ended design (n_max = Inf): set the simulation horizon n_cap")
  }

  # Patients simulated per basket
  if(is.null(design$n_cap)){
    design$n_sim <- design$n_max
  }else{
    design$n_sim <- pmin(design$n_max, design$n_cap)
  }
  return(design)
}

################################################################################
# Fit the borrowing model once for every distinct data state
################################################################################
priors_and_summaries <- function(engine, y, n, need, seed){
  # y, n = (trials x baskets) data; need = where a posterior is needed
  # Trials with the same data share one fit.
  rows <- which(rowSums(need) > 0)
  key <- do.call(paste, c(as.data.frame(cbind(y, n)[rows, , drop = FALSE]), sep = "."))
  u <- unique(key)
  iu <- match(key, u)
  first <- rows[match(u, key)]
  needu <- rowsum(need[rows, , drop = FALSE] * 1, iu) > 0

  po <- engine$posterior(y[first, , drop = FALSE], n[first, , drop = FALSE], needu, seed)
  state <- integer(nrow(y))
  state[rows] <- iu
  return(list(po = po, state = state))
}

################################################################################
# Log posterior of the pairs (trial, basket)
################################################################################
pairs_logpost <- function(fit, ix, S){
  # ix = linear index into the (trials x baskets) matrices; S = number of trials
  s <- (ix - 1) %% S + 1
  k <- (ix - 1) %/% S + 1
  G <- dim(fit$po$logpost)[3]
  i <- fit$state[s]
  return(matrix(fit$po$logpost[cbind(rep(i, G), rep(k, G), rep(seq_len(G), each = length(ix)))],
                length(ix), G))
}

################################################################################
# Simulate trials with the e-value of framework.tex
################################################################################
simulate_trials <- function(p_true, design, engine, n_sim, seed = 1, fit_seed = NULL, X = NULL,
                            trace = FALSE){
  # p_true = true response rate per basket
  # design = see complete_design; engine = borrowing model
  # n_sim = number of trials (simulated together)
  # seed = seed of the patient outcomes
  # fit_seed = seed of the borrowing model (default seed + 10^6); a separate stream,
  #            so the model's sampling never reuses the outcome stream
  # X = optional (trials x baskets x patients) outcomes, instead of simulating them
  # trace = keep the e-values and prior means at every analysis
  #
  # Global analyses at look_every, 2 look_every, ... At every analysis:
  #   1. every recruiting basket multiplies E (and E^f) by the factor of the patients
  #      since the previous analysis, with the prior frozen at the previous analysis;
  #   2. decisions: reject if E >= 1/alpha; futility; complete at a finite n_max;
  #      'open' at the simulation horizon n_cap;
  #   3. the borrowing model is fitted to ALL platform data (stopped baskets included)
  #      and gives the priors of the next block.
  # So every prior depends on the data before its block only, and is frozen during it.
  design <- complete_design(design)
  if(is.null(fit_seed)){
    fit_seed <- seed + 10^6
  }
  grid <- engine$grid
  K <- nrow(design$layout)
  S <- n_sim
  nmx <- matrix(design$n_max, S, K, byrow = TRUE)
  ncap <- matrix(design$n_sim, S, K, byrow = TRUE)

  # Outcomes of patient i of basket k in trial s, and their cumulative sums
  if(is.null(X)){
    set.seed(seed)
    X <- array(runif(S * K * max(design$n_sim)) < rep(p_true, each = S),
               c(S, K, max(design$n_sim)))
  }
  CY <- array(0L, c(S, K, dim(X)[3] + 1))
  for(i in seq_len(dim(X)[3])){
    CY[, , i + 1] <- CY[, , i] + X[, , i]
  }

  # Thresholds and analysis times
  thr <- log(1 / design$alpha)
  lmax <- log(1 / grid$p0)   # largest log factor per patient (all respond)
  efut <- design$futility == "evalue"
  if(efut){
    fut <- futility_side(grid, design$p1)
    thr_f <- log(1 / design$beta)
  }
  times <- seq(design$look_every,
               max(design$entry + design$n_sim / design$rate) + design$look_every,
               by = design$look_every)

  # State of every basket. status: 0 recruiting, 1 rejected, 2 stopped for futility,
  # 3 complete at n_max, 4 open at the simulation horizon
  n <- matrix(0L, S, K)
  y <- matrix(0L, S, K)
  logE <- matrix(0, S, K)
  logEf <- matrix(0, S, K)
  status <- matrix(0L, S, K)
  res <- list(time = matrix(NA_real_, S, K), post_mean = matrix(NA_real_, S, K),
              p_above = matrix(NA_real_, S, K), ci_lo = matrix(NA_real_, S, K),
              ci_hi = matrix(NA_real_, S, K), pex = matrix(NA_real_, S, K),
              prior_mean = matrix(NA_real_, S, K), rhat = matrix(NA_real_, S, K))

  make_prior <- function(lp){
    w <- restrict(lp, grid)
    if(design$prior_form == "beta"){
      w <- scaled_beta(w, grid)
    }
    return(w)
  }

  # Priors of the first block: fitted to the empty platform (the same for all trials).
  # W has one row per (trial, basket): row s + (k - 1) S.
  fit <- priors_and_summaries(engine, n, n, matrix(TRUE, S, K), fit_seed)
  lp_first <- pairs_logpost(fit, (seq_len(K) - 1L) * S + 1L, S)
  W <- make_prior(lp_first)[rep(seq_len(K), each = S), , drop = FALSE]
  check_priors(W, grid)
  if(efut){
    Wf <- restrict_below(lp_first, fut)[rep(seq_len(K), each = S), , drop = FALSE]
    check_priors_below(Wf, fut)
  }

  tr <- list()
  for(a in seq_along(times)){
    act <- status == 0

    # New patients since the previous analysis
    n_new <- n
    n_new[act] <- matrix(planned_n(design, times[a]), S, K, byrow = TRUE)[act]
    y_new <- matrix(CY[cbind(rep(seq_len(S), K), rep(seq_len(K), each = S), as.vector(n_new) + 1)], S, K)

    # 1. Block factors with the frozen priors
    upd <- which(act & n_new > n)
    if(length(upd)){
      dy <- (y_new - y)[upd]
      dn <- (n_new - n)[upd]
      logE[upd] <- logE[upd] + log_factor(W[upd, , drop = FALSE], dy, dn, grid)
      if(efut){
        logEf[upd] <- logEf[upd] + log_factor(Wf[upd, , drop = FALSE], dy, dn, fut)
      }
    }
    n <- n_new
    y <- y_new

    # 2. Decisions
    status[act & logE >= thr] <- 1L
    if(efut){
      status[act & status == 0 & logEf >= thr_f] <- 2L
    }
    status[act & status == 0 & n >= nmx] <- 3L
    if(design$futility == "unreachable"){
      status[act & status == 0 & logE + (nmx - n) * lmax < thr] <- 2L
    }
    status[act & status == 0 & n >= ncap] <- 4L
    now <- act & status != 0
    res$time[now] <- times[a]
    res$prior_mean[now] <- (W %*% grid$support)[now]

    need <- status == 0 | now
    if(trace){
      tr[[a]] <- list(time = times[a], n = n, y = y, logE = logE, logEf = logEf,
                      prior_mean = matrix(W %*% grid$support, S, K))
    }
    if(!any(need)){
      break
    }

    # 3. Fit the borrowing model to all data now
    fit <- priors_and_summaries(engine, y, n, need, fit_seed + a)

    # New frozen priors for the next block of every recruiting basket
    ix <- which(status == 0)
    if(length(ix)){
      lp <- pairs_logpost(fit, ix, S)
      W[ix, ] <- make_prior(lp)
      if(efut){
        Wf[ix, ] <- restrict_below(lp, fut)
      }
    }

    # Posterior summaries of the baskets decided now
    ix <- which(now)
    if(length(ix)){
      lp <- pairs_logpost(fit, ix, S)
      pr <- exp(lp - row_max(lp))
      pr <- pr / rowSums(pr)
      cp <- t(apply(pr, 1, cumsum))
      res$post_mean[ix] <- pr %*% grid$p
      res$p_above[ix] <- rowSums(pr[, grid$above, drop = FALSE])
      res$ci_lo[ix] <- grid$p[max.col(cp >= 0.05, "first")]
      res$ci_hi[ix] <- grid$p[max.col(cp >= 0.95, "first")]
      s <- (ix - 1) %% S + 1
      k <- (ix - 1) %/% S + 1
      res$pex[ix] <- fit$po$pex[cbind(fit$state[s], k)]
      if(!is.null(fit$po$rhat)){
        res$rhat[ix] <- fit$po$rhat[cbind(fit$state[s], k)]
      }
    }
  }

  # One row per (trial, basket)
  out <- data.frame(trial = rep(seq_len(S), K), basket = rep(design$layout$basket, each = S),
                    row = rep(design$layout$row, each = S), col = rep(design$layout$col, each = S),
                    p_true = rep(p_true, each = S), null = rep(p_true <= grid$p0, each = S),
                    reject = as.vector(status == 1), stopped = as.vector(status == 2),
                    open = as.vector(status == 4), n = as.vector(n), y = as.vector(y),
                    logE = as.vector(logE), logEf = as.vector(logEf), lapply(res, as.vector))
  out$cover <- out$ci_lo <= out$p_true & out$p_true <= out$ci_hi
  if(trace){
    attr(out, "trace") <- tr
  }
  return(out)
}

################################################################################
# Design of one simulation cell (from the settings)
################################################################################
make_design <- function(recruitment, model, horizon){
  # recruitment = element of Recruitment (setting.R); model = one of Models
  # horizon = 'fixed' (at most n_max per basket) or 'open' (no maximum, e-value futility)
  # Uses the settings p0, p1, alpha, beta, J, L, n_max, n_cap and look_every.
  in_trial <- !is.na(recruitment$rate)
  all_cells <- grid_layout(J, L)
  missing <- lapply(which(!in_trial), function(i){c(all_cells$row[i], all_cells$col[i])})

  design <- list(layout = grid_layout(J, L, missing = missing), n_max = n_max,
                 rate = recruitment$rate[in_trial], entry = recruitment$entry[in_trial],
                 look_every = look_every, alpha = alpha)
  if(model == "one-way EXNEX, scaled Beta"){
    design$prior_form <- "beta"
  }
  if(horizon == "open"){
    design$n_max <- Inf
    design$n_cap <- n_cap
    design$futility <- "evalue"
    design$p1 <- p1
    design$beta <- beta
  }
  return(design)
}

################################################################################
# Borrowing model of one simulation cell (from the settings)
################################################################################
make_engine <- function(model, layout){
  # model = one of Models (setting.R); layout = layout of the trial
  grid <- theta_grid(p0)
  pbar <- (p0 + p1) / 2
  if(model == "none"){
    return(engine_none(grid, pbar))
  }
  if(model %in% c("one-way EXNEX", "one-way EXNEX, scaled Beta")){
    return(engine_oneway(grid, pbar, n_max = 30))
  }
  if(model == "two-way EXNEX"){
    return(engine_twoway(grid, layout, pbar, n_max = 30, mcmc = mcmc))
  }
  stop("unknown model")
}


################################################################################
################################################################################
# Part 4. Summaries and figures
#
# The guarantee is per basket: P(reject) <= alpha for every basket with p <= p0.
# The platform-wide FWER (any null basket rejected) is not controlled at alpha and is
# reported separately. Open-ended designs: sample sizes are truncated at n_cap; the
# columns open / open_any report how often n_cap was reached.
################################################################################
################################################################################

################################################################################
# Results per basket
################################################################################
summarise_baskets <- function(d, by = c("recruitment", "scenario", "engine")){
  # d = output of simulate_trials with the grouping columns 'by'
  g <- interaction(d[, c(by, "basket")], drop = TRUE, lex.order = TRUE)

  one_basket <- function(x){
    rej <- mean(x$reject)
    out <- data.frame(x[1, c(by, "basket", "row", "col", "p_true", "null")],
                      trials = nrow(x),
                      reject = rej,
                      reject_se = sqrt(rej * (1 - rej) / nrow(x)),
                      stopped = mean(x$stopped),
                      open = if(is.null(x$open)) NA else mean(x$open),
                      n_mean = mean(x$n),
                      time_mean = mean(x$time),
                      post_bias = mean(x$post_mean - x$p_true),
                      cover90 = mean(x$cover),
                      pex = mean(x$pex),
                      prior_mean = mean(x$prior_mean),
                      log10E = mean(x$logE) / log(10),
                      rhat_max = suppressWarnings(max(x$rhat, na.rm = TRUE)))
    return(out)
  }

  out <- do.call(rbind, lapply(split(d, g), one_basket))
  out$rhat_max[!is.finite(out$rhat_max)] <- NA
  rownames(out) <- NULL
  return(out)
}

################################################################################
# Results per platform: per-basket type I error, FWER, power and sample size
################################################################################
summarise_platform <- function(d, by = c("recruitment", "scenario", "engine")){
  # d = output of simulate_trials with the grouping columns 'by'
  g <- interaction(d[, by], drop = TRUE, lex.order = TRUE)

  one_platform <- function(x){
    tr <- split(x, x$trial)
    nul <- x[x$null, ]
    alt <- x[!x$null, ]
    if(nrow(nul)){
      typeI <- tapply(nul$reject, nul$basket, mean)
    }else{
      typeI <- NA
    }
    out <- data.frame(x[1, by],
      trials = length(tr),
      baskets_null = length(unique(nul$basket)),
      typeI_basket_max = max(typeI),
      typeI_basket_mean = mean(typeI),
      FWER = if(nrow(nul)) mean(vapply(tr, function(t){any(t$reject & t$null)}, TRUE)) else NA,
      power_basket_mean = if(nrow(alt)) mean(alt$reject) else NA,
      power_any = if(nrow(alt)) mean(vapply(tr, function(t){any(t$reject & !t$null)}, TRUE)) else NA,
      all_correct = mean(vapply(tr, function(t){all(t$reject == !t$null)}, TRUE)),
      open_any = if(is.null(x$open)) NA else mean(vapply(tr, function(t){any(t$open)}, TRUE)),
      N_mean = mean(vapply(tr, function(t){sum(t$n)}, 0)),
      duration = mean(vapply(tr, function(t){max(t$time)}, 0)),
      cover90 = mean(x$cover),
      mean_abs_error = mean(abs(x$post_mean - x$p_true)))
    return(out)
  }

  out <- do.call(rbind, lapply(split(d, g), one_platform))
  rownames(out) <- NULL
  return(out)
}

# Colours of the borrowing models in the figures
engine_colours <- c("none" = "#7a7a7a", "one-way EXNEX" = "#0072B2",
                    "one-way EXNEX, scaled Beta" = "#56B4E9", "two-way EXNEX" = "#D55E00")

################################################################################
# Figure: rejection rate per basket
################################################################################
plot_baskets <- function(b, recruitment, alpha){
  # b = summarise_baskets output; null baskets (open symbols) against alpha,
  # active baskets (filled symbols) show the power
  library(ggplot2)
  b <- b[b$recruitment == recruitment, ]
  b$type <- ifelse(b$null, "null basket (type I error)", "active basket (power)")
  p <- ggplot(b, aes(basket, reject, colour = engine, shape = type)) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = pmax(reject - 2 * reject_se, 0), ymax = pmin(reject + 2 * reject_se, 1)),
                   position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_shape_manual(values = c("null basket (type I error)" = 1, "active basket (power)" = 16)) +
    scale_colour_manual(values = engine_colours) +
    facet_wrap(~ scenario, ncol = 2) +
    labs(x = "basket (tumour type . biomarker group)", y = "rejection rate", colour = "borrowing", shape = NULL,
         title = sprintf("Rejection per basket (%s recruitment); dashed line: alpha per basket", recruitment)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.box = "vertical",
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

################################################################################
# Figure: power against expected total sample size
################################################################################
plot_tradeoff <- function(pl, recruitment){
  # pl = summarise_platform output; one panel per scenario with active baskets
  library(ggplot2)
  pl <- pl[pl$recruitment == recruitment & !is.na(pl$power_basket_mean), ]
  p <- ggplot(pl, aes(N_mean, power_basket_mean, colour = engine)) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = engine_colours) +
    facet_wrap(~ scenario, ncol = 3, scales = "free_x") +
    labs(x = "expected total sample size", y = "mean power over active baskets", colour = "borrowing",
         title = sprintf("Power against sample size (%s recruitment)", recruitment)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  return(p)
}

################################################################################
# Save a figure as pdf and png
################################################################################
save_plot <- function(p, file, width, height){
  ggplot2::ggsave(paste0(file, ".pdf"), p, width = width, height = height)
  ggplot2::ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 150)
}
