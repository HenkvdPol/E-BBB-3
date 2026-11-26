## Run from anywhere:  Rscript checks/bank_test.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
s <- settings$for_hypothesis(settings, 0.3, 0.5); s$sim$n_cores <- 2
for (type in c("grow", "ebased")) { b <- compact_bank(build_bank(type, 8, c(4, 8), TRUE, s))
  cat(type, "borrow bank: theta1 values", length(b$theta1), " distinct bets", length(b$b_val0) - 1,
      " bets in [0,1]:", all(b$b_val0 >= 0 & b$b_val0 <= 1), " dims", dim(b$bet_idx), "\n") }
## GROW bets must differ by theta1 (Kelly bet), check two theta1 values
b <- compact_bank(build_bank("grow", 8, c(4, 8), TRUE, s)); i1 <- which(b$m_grid == 1)
cat("GROW bet at m = 1 for theta1 =", b$theta1[c(1, 21)], ":", b$b_val0[as.integer(b$bet_idx[i1, 1, c(1, 21)]) + 1],
    " Kelly:", kelly_bet(0.3, b$theta1[c(1, 21)]), "\n")
