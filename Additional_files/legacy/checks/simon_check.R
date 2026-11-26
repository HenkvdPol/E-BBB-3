## Run from anywhere:  Rscript checks/simon_check.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(dirname(check_dir), "additional_functions.R"))
pub <- read.table(header = TRUE, text = "
p0 p1 a b o_r1 o_n1 o_r o_n m_r1 m_n1 m_r m_n
0.05 0.25 0.10 0.10 0 9 2 24 0 13 2 20
0.05 0.25 0.05 0.20 0 9 2 17 0 12 2 16
0.05 0.25 0.05 0.10 0 9 3 30 0 15 3 25
0.10 0.30 0.10 0.10 1 12 5 35 1 16 4 25
0.10 0.30 0.05 0.20 1 10 5 29 1 15 5 25
0.10 0.30 0.05 0.10 2 18 6 35 2 22 6 33
0.20 0.40 0.10 0.10 3 17 10 37 3 19 10 36
0.20 0.40 0.05 0.20 3 13 12 43 4 18 10 33
0.20 0.40 0.05 0.10 4 19 15 54 5 24 13 45
0.30 0.50 0.10 0.10 7 22 17 46 7 28 15 39
0.30 0.50 0.05 0.20 5 15 18 46 6 19 16 39
0.30 0.50 0.05 0.10 8 24 24 63 7 24 21 53
0.40 0.60 0.10 0.10 7 18 22 46 11 28 20 41
0.40 0.60 0.05 0.20 7 16 23 46 17 34 20 39
0.40 0.60 0.05 0.10 11 25 32 66 12 29 27 54
0.50 0.70 0.10 0.10 11 21 26 45 11 23 23 39
0.50 0.70 0.05 0.20 8 15 26 43 12 23 23 37
0.50 0.70 0.05 0.10 13 24 36 61 14 27 32 53
0.60 0.80 0.10 0.10 6 11 26 38 18 27 24 35
0.60 0.80 0.05 0.20 7 11 30 43 8 13 25 35
0.60 0.80 0.05 0.10 12 19 37 53 15 26 32 45
0.70 0.90 0.10 0.10 6 9 22 28 11 16 20 25
0.70 0.90 0.05 0.20 4 6 22 27 19 23 21 26
0.70 0.90 0.05 0.10 11 15 29 36 13 18 26 32")
t0 <- Sys.time()
res <- do.call(rbind, lapply(seq_len(nrow(pub)), function(i) { x <- pub[i, ]
  o <- simon_design(x$p0, x$p1, x$a, x$b, type = "optimal"); m <- simon_design(x$p0, x$p1, x$a, x$b, type = "minimax")
  data.frame(x[1:4], opt = sprintf("%d/%d %d/%d", o$r1, o$n1, o$r, o$n), opt_pub = sprintf("%d/%d %d/%d", x$o_r1, x$o_n1, x$o_r, x$o_n),
             EN0 = round(o$EN0, 1), mm = sprintf("%d/%d %d/%d", m$r1, m$n1, m$r, m$n), mm_pub = sprintf("%d/%d %d/%d", x$m_r1, x$m_n1, x$m_r, x$m_n)) }))
res$same <- res$opt == res$opt_pub & res$mm == res$mm_pub
print(res, row.names = FALSE)
cat("all identical:", all(res$same), " seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
## fixed n = Simon's optimal n gives the optimal design again (used by references())
f <- simon_design(0.3, 0.5, 0.05, 0.2, n = 46); cat("fixed n = 46:", f$r1, f$n1, f$r, f$n, "\n")
