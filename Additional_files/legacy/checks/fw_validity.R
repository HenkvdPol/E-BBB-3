## Run from anywhere:  Rscript checks/fw_validity.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
s <- settings$for_hypothesis(settings, 0.10, 0.30)
s$K <- 2; s$n_max <- c(10, 10); s$n_ext <- c(5L, 5L); s$series <- "check"
s$schedules <- list("Fully sequential" = 1, "Blocked (3)" = 3); s$scenarios <- list(x = c(0.1, 0.1))
s$setups$basket$entry <- c(0, 0); s$sim$bank_dir <- file.path(out_dir, "banks_k2"); s$sim$n_cores <- 2
for (d in names(s$designs)) if (isTRUE(s$designs[[d]]$borrow)) s$designs[[d]] <- modifyList(list(summary = "mean", quantile = 0.25, redesign = "analysis"), s$designs[[d]])
designs <- c("e-mixture (framework)", "e-mixture (no borrowing)", "e-based + mEXNEX", "e-based + mEXNEX (flexible)")
s$setups$basket$designs <- designs; s$setups$umbrella <- NULL
eng <- build_engine(s); banks <- get_banks(s)
N <- 15; Y1 <- sapply(seq_len(N), function(t) bitwAnd(bitwShiftR(0:(2^N - 1), t - 1), 1L)); S <- nrow(Y1); ny <- rowSums(Y1)
lag <- function(Y) cbind(0L, Y[, -N])
paths <- list("all responses" = matrix(1L, S, N), "no responses" = matrix(0L, S, N),
              "responses, then none" = matrix(rep(c(rep(1L, 7), rep(0L, 8)), each = S), S),
              "copies cohort 1 (lag 1)" = lag(Y1), "opposite of cohort 1 (lag 1)" = 1L - lag(Y1))
out <- list()
for (pn in names(paths)) {
  C <- array(0L, c(S, 2, N + 1)); C[, 1, -1] <- t(apply(Y1, 1, cumsum)); C[, 2, -1] <- t(apply(paths[[pn]], 1, cumsum))
  for (sn in names(s$schedules)) for (dn in designs) { d <- s$designs[[dn]]; block <- s$schedules[[sn]]
    r <- switch(d$type,
      mixture = run_trials_mixture(d, C, block, c(0, 0), s, eng, new_cache()),
      evalue  = run_trials_evalue(d, C, block, c(0, 0), s, eng, new_cache(), banks_for(d, block, s$setups$basket, s, banks),
                                  if (isTRUE(d$extend)) ext_banks_for(d, block, s$setups$basket, s, banks)))
    rej <- r$status[, 1] == 1
    out[[length(out) + 1]] <- data.frame(other_cohort = pn, schedule = sn, design = dn,
      typeI_at_p0 = sum(0.10^ny * 0.90^(N - ny) * rej), typeI_at_0.05 = sum(0.05^ny * 0.95^(N - ny) * rej),
      power_at_0.30 = sum(0.30^ny * 0.70^(N - ny) * rej))
  }
}
res <- do.call(rbind, out); saveRDS(res, file.path(out_dir, "fw_validity.rds"))
print(format(res, digits = 3), row.names = FALSE)
cat("\nlargest type I error of cohort 1 over all other-cohort paths, schedules and designs:", signif(max(res$typeI_at_p0), 4),
    " (all <= alpha:", all(res$typeI_at_p0 <= s$alpha & res$typeI_at_0.05 <= s$alpha), ")\n")
