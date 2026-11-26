## Run from anywhere:  Rscript checks/fw_engine.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings); s <- complete_settings(settings, runs[runs$name == "grid_p10-30_K5_n30", ]); eng <- build_engine(s); p0 <- s$p0
## 1. closed-form block factor vs numerical integration over the scaled Beta prior
set.seed(1); worst <- 0
for (i in 1:200) { a <- runif(1, 0.2, 30); b <- runif(1, 0.2, 60); n <- sample(1:20, 1); y <- sample(0:n, 1)
  f <- function(u) { p <- p0 + (1 - p0) * u; (p / p0)^y * ((1 - p) / (1 - p0))^(n - y) * dbeta(u, a, b) }
  num <- integrate(f, 0, 1, rel.tol = 1e-10, subdivisions = 2000)$value
  worst <- max(worst, abs(mixture_factor(a, b, y, n, p0) / num - 1)) }
cat("1. block factor, closed form vs numerical integration: largest relative difference", signif(worst, 3), "(200 random cases)\n")
## 2. restricted posterior moments from the engine vs direct computation on the theta grid
g <- s$grid; theta <- seq(g$theta_min, g$theta_max, by = g$theta_step); pth <- expit(theta)
hyp <- expand.grid(mu = seq(g$mu_min, g$mu_max, by = g$mu_step), tau = seq(g$tau_step / 2, g$tau_max, by = g$tau_step))
W <- dnorm(sweep(outer(theta, hyp$mu, "-"), 2, hyp$tau, "/")); W <- sweep(W, 2, colSums(W), "/")
V <- dnorm(theta, logit(s$prior$nex_p), sqrt(1 / s$prior$nex_p + 1 / (1 - s$prior$nex_p))); V <- V / sum(V)
above <- pmin(pmax((theta + g$theta_step / 2 - logit(p0)) / g$theta_step, 0), 1)
direct <- function(y, n, pi, k) {             # posterior of theta_k given all cohorts, on the grid
  D <- sapply(seq_along(y), function(j) { L <- dbinom(y[j], n[j], pth); pi[j] * colSums(L * W) + (1 - pi[j]) * sum(L * V) })
  lw <- eng$logprior + rowSums(log(D)) - log(D[, k]); w <- exp(lw - max(lw)); w <- w / sum(w)       # leave k out
  post <- dbinom(y[k], n[k], pth) * (pi[k] * as.vector(W %*% w) + (1 - pi[k]) * V); post <- post / sum(post)
  list(post = post, rm = sum(post * above * pth) / sum(post * above), rs = sum(post * above * pth^2) / sum(post * above))
}
states <- list(c(3, 1, 0, 5, 2), c(9, 8, 2, 3, 4), c(0, 0, 1, 0, 0), c(6, 6, 6, 6, 6))
nn <- list(c(10, 10, 10, 10, 10), c(30, 30, 12, 20, 15), c(20, 20, 20, 20, 20), c(20, 20, 20, 20, 20))
worst_m <- 0; tab <- list()
for (i in seq_along(states)) {
  y <- matrix(states[[i]], 1); n <- matrix(nn[[i]], 1); pi <- matrix(0.5, 1, 5)
  po <- posterior_exnex(eng, y, n, pi, moments = TRUE)
  for (k in 1:5) { d <- direct(states[[i]], nn[[i]], rep(0.5, 5), k)
    worst_m <- max(worst_m, abs(po$rmean[k] - d$rm), abs(po$rsec[k] - d$rs))
    ## 3. block factor with the matched prior vs with the exact restricted posterior (block of 5 with 0..5 responses)
    pr <- mixture_prior(po$rmean[k], po$rsec[k], p0); keep <- above > 0; wr <- d$post * above / sum(d$post * above)
    for (yb in c(0, 2, 5)) { ex <- sum(wr * (pth / p0)^yb * ((1 - pth) / (1 - p0))^(5 - yb))
      tab[[length(tab) + 1]] <- data.frame(state = i, cohort = k, y_n = paste0(states[[i]][k], "/", nn[[i]][k]), block = paste0(yb, "/5"),
        exact = ex, matched = mixture_factor(pr$a, pr$b, yb, 5, p0)) } }
}
cat("2. restricted posterior mean and second moment, engine vs direct: largest difference", signif(worst_m, 3), "\n")
tab <- do.call(rbind, tab); tab$ratio <- tab$matched / tab$exact
cat("3. block factor with the matched prior / with the exact restricted posterior: range", round(range(tab$ratio), 3),
    "; median", round(median(tab$ratio), 3), "\n")
print(format(tab[tab$cohort == 1 | tab$state == 2, ], digits = 3), row.names = FALSE)
