## =============================================================================
##  ebbb_compare.R -- the eBBB priors of additional_functions.R against the eBBB code itself
##  Run from the repository folder, with the path to eBBB's lr_evalue.R (not part of this repository):
##      Rscript Additional_files/validation/ebbb_compare.R <path>/eBBB/code/LR_evalue_multiple/lr_evalue.R
##  Output: Additional_files/validation/results/ebbb_compare.csv (largest relative difference per prior and block)
##  The eBBB file is only sourced into its own environment; nothing of it is copied or changed.
##
##  For 40 random data states on the 4 x 3 grid (n per basket 0, 5, 10, 20 or 30) and every basket:
##  the block factor of a block of 5 with 0, 2 or 5 responses under the restricted prior, from
##  engine_ebbb() on the grid and from eBBB's closed forms (block_factor), and the EX/NEX weights.
## =============================================================================
args <- commandArgs(TRUE)
if (length(args) != 1 || !file.exists(args[1])) stop("give the path to eBBB's lr_evalue.R")
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "Additional_files/validation", "..", "..")) })
ours <- new.env(); sys.source(file.path(root, "additional_functions.R"), envir = ours)
his <- new.env(); sys.source(args[1], envir = his)
out_dir <- file.path(root, "Additional_files", "validation", "results"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p0 <- 0.10; J <- 4; L <- 3; K <- J * L
grid <- ours$theta_grid(p0); lay <- ours$grid_layout(J, L)
idx_his <- lay$row + (lay$col - 1) * J                      # our basket k (row-major) -> eBBB index (column-major)
same_row <- his$grid_same_row(J, L); same_col <- his$grid_same_col(J, L)
set.seed(7); S <- 40
N <- matrix(sample(c(0, 5, 10, 20, 30), S * K, TRUE), S); Y <- matrix(rbinom(S * K, N, 0.25), S)
blocks <- list(c(0, 5), c(2, 5), c(5, 5))
methods <- list(
  "no borrowing"  = list(m = "no borrowing", w = c(0, 0),     his = function(c, y, n) his$borrowing_prior(c, y, n, same_row, same_col, 0, 0, 1, 1)),
  "power, equal"  = list(m = "power",        w = c(0.5, 0.5), his = function(c, y, n) his$borrowing_prior(c, y, n, same_row, same_col, 0.5, 0.5, 1, 1)),
  "power, column" = list(m = "power",        w = c(0, 1),     his = function(c, y, n) his$borrowing_prior(c, y, n, same_row, same_col, 0, 1, 1, 1)),
  "EX/NEX"        = list(m = "EX/NEX",       w = c(0, 1),     his = function(c, y, n) his$exnex_prior(c, y, n, same_row, same_col, 0, 1, 1, 1, 0.5)),
  "similarity"    = list(m = "similarity",   w = c(0, 0),     his = function(c, y, n) his$similarity_prior(c, y, n, 1, 1, 2, 0.3)))
out <- NULL
for (nm in names(methods)) {
  md <- methods[[nm]]
  po <- ours$engine_ebbb(grid, lay, md$m, 1, 1, w_row = md$w[1], w_col = md$w[2], p_ex = 0.5, epsilon = 2, tau = 0.3)$posterior(Y, N)
  for (s in 1:S) for (k in 1:K) {
    yh <- numeric(K); nh <- numeric(K); yh[idx_his] <- Y[s, ]; nh[idx_his] <- N[s, ]
    pr <- md$his(idx_his[k], yh, nh)
    w <- ours$restrict(matrix(po$logpost[s, k, ], 1), grid)
    for (bl in blocks) {
      out <- rbind(out, data.frame(prior = nm, block = sprintf("%d/%d", bl[1], bl[2]),
                                   rel_diff = exp(ours$log_factor(w, bl[1], bl[2], grid)) / his$block_factor(pr, bl[1], bl[2], p0) - 1,
                                   ex_diff = if (nm == "EX/NEX") po$pex[s, k] - pr$w[["ex"]] else NA))
    }
  }
}
res <- aggregate(cbind(max_abs_rel_diff = abs(rel_diff)) ~ prior + block, out, max)

## Oracle: point mass at q, exact on the grid (add_points)
q <- list(true = c(0.10, 0.30, 0.45, 0.05), eff = c(0.30, 0.30, 0.45, 0.30), fut = c(0.10, 0.10, 0.10, 0.05))
or <- ours$engine_ebbb(grid, ours$grid_layout(2, 2), "oracle", q = q)
po <- or$posterior(matrix(0, 1, 4), matrix(0, 1, 4))
for (bl in blocks) {
  d <- sapply(1:4, function(k) exp(ours$log_factor(ours$restrict(matrix(po$logpost[1, k, ], 1), or$grid), bl[1], bl[2], or$grid)) /
                  his$block_factor(list(type = "point", q = q$eff[k]), bl[1], bl[2], p0) - 1)
  res <- rbind(res, data.frame(prior = "oracle", block = sprintf("%d/%d", bl[1], bl[2]), max_abs_rel_diff = max(abs(d))))
}
res$max_abs_ex_weight_diff <- ifelse(res$prior == "EX/NEX", max(abs(out$ex_diff), na.rm = TRUE), NA)
write.csv(res, file.path(out_dir, "ebbb_compare.csv"), row.names = FALSE)
print(res, digits = 3)
