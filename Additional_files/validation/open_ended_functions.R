# File with the exact calculation of one basket without a maximum sample size
# (feasibility of the open-ended design; Additional_files/docs/open-ended-design.md).
#
# Telescoping. If the prior of every block is the posterior from ONE fixed starting
# distribution pi_0, restricted to (p0, 1], the product of the block factors is
#   E = sum_g w0_g (p_g / p0)^y ((1 - p_g) / (1 - p0))^(n - y),  w0 = pi_0 restricted,
# because every factor is a ratio of consecutive normalising constants. This holds
# without borrowing, and with borrowing from baskets whose data no longer change
# (pi_0 = the borrowed predictive given their data). The same holds for the futility
# e-value with pi_0 restricted to [0, p1). E and E^f then only depend on (y, n), so a
# dynamic programme over y gives the exact probabilities of rejection, futility stop
# and of still being open after any number of patients. With borrowing from baskets
# that are still recruiting the factors do not telescope; that needs simulate.R.
# Needs the functions of additional_functions.R.


################################################################################
# log E (or log E^f) after y responses out of n patients
################################################################################
mixture_logE <- function(w0, side, y, n){
  # w0 = starting prior weights; side = theta grid (efficacy) or futility_side (futility)
  # y = vector of numbers of responses; n = number of patients
  lr <- outer(y, side$l1) + outer(n - y, side$l0) + rep(log(w0), each = length(y))
  m <- row_max(lr)
  return(m + log(rowSums(exp(lr - m))))
}

################################################################################
# Exact operating characteristics of one basket
################################################################################
exact_basket <- function(p, w0, grid, alpha, block = 1, n_max = Inf, horizon = 3000,
                         futility = "none", w0f = NULL, fut = NULL, beta = NULL,
                         report_at = c(30, 100, 300, 1000, 3000)){
  # p = true response rate; w0 = starting prior on (p0, 1]; grid = theta grid
  # alpha = type I error; block = patients between analyses
  # n_max = maximum per basket (Inf = open-ended); horizon = patients followed at most
  # futility = "none", "evalue" (E^f >= 1/beta; needs w0f, fut and beta) or
  #            "unreachable" (1/alpha out of reach by n_max; needs a finite n_max)
  # report_at = sample sizes at which P(rejected so far) is reported
  thr <- log(1 / alpha)
  end <- min(n_max, horizon)
  step <- dbinom(0:block, block, p)

  # prob[y + 1] = P(still open, y responses, n patients)
  prob <- 1
  n <- 0
  stop_n <- c()
  rej <- c()
  fail <- c()
  by <- setNames(rep(NA_real_, length(report_at)), report_at)
  cum_rej <- 0

  while(n < end && sum(prob) > 1e-14){
    # The next block of b patients
    b <- min(block, end - n)
    if(b == block){
      st <- step
    }else{
      st <- dbinom(0:b, b, p)
    }
    new <- numeric(length(prob) + b)
    for(j in 0:b){
      new[j + seq_along(prob)] <- new[j + seq_along(prob)] + prob * st[j + 1]
    }
    prob <- new
    n <- n + b

    # Decisions at this analysis, for every number of responses still open
    live <- which(prob > 0)
    y <- live - 1
    r <- 0
    f <- 0
    if(length(live)){
      hit <- mixture_logE(w0, grid, y, n) >= thr
      r <- sum(prob[live[hit]])
      prob[live[hit]] <- 0
      if(futility == "evalue"){
        live <- live[!hit]
        y <- live - 1
        if(length(live)){
          hf <- mixture_logE(w0f, fut, y, n) >= log(1 / beta)
          f <- sum(prob[live[hf]])
          prob[live[hf]] <- 0
        }
      }else if(futility == "unreachable"){
        live <- live[!hit]
        y <- live - 1
        if(length(live)){
          hf <- mixture_logE(w0, grid, y, n) + (n_max - n) * log(1 / grid$p0) < thr
          f <- sum(prob[live[hf]])
          prob[live[hf]] <- 0
        }
      }
    }
    stop_n <- c(stop_n, n)
    rej <- c(rej, r)
    fail <- c(fail, f)
    cum_rej <- cum_rej + r
    by[as.character(report_at[report_at == n])] <- cum_rej
  }

  open <- sum(prob)
  stopped <- rej + fail
  by[as.numeric(names(by)) > n] <- cum_rej   # nothing more can happen after the end

  # Sample size by which a fraction a of the baskets has stopped
  q <- function(a){
    cs <- cumsum(stopped)
    if(max(cs) < a){
      return(NA)
    }
    return(stop_n[which(cs >= a)[1]])
  }

  ended <- is.finite(n_max) && n >= n_max
  summary <- data.frame(p = p, reject = sum(rej), futility = sum(fail),
                        complete = if(ended) open else 0,
                        open = if(ended) 0 else open,
                        EN = sum(stop_n * stopped) + n * open,
                        N_median = q(0.5), N_q90 = q(0.9), last_n = n)
  return(list(summary = summary, reject_by = by))
}
