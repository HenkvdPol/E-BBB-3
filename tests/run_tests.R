## =============================================================================
##  run_tests.R -- runs every tests/test_*.R (base R only; no test package needed)
##
##  Rscript tests/run_tests.R            fast tests (CI)
##  Rscript tests/run_tests.R --slow     also the slow reference comparisons
##  Rscript tests/run_tests.R oneway     only files whose name contains 'oneway'
##  Each result is PASS, FAIL or SKIP (with the reason); the exit status is 1 if anything FAILs.
##  Details: tests/output/test_results.csv
## =============================================================================
root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(if (length(f)) dirname(normalizePath(f[1])) else "tests", ".."))
})
for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
args <- commandArgs(TRUE); SLOW <- "--slow" %in% args; pattern <- setdiff(args, "--slow")

.results <- data.frame(file = character(), test = character(), status = character(), info = character())
.file <- ""
record <- function(test, status, info = "") {
  .results[nrow(.results) + 1, ] <<- list(.file, test, status, info)
  cat(sprintf("  %-4s %s%s\n", status, test, if (nzchar(info)) paste0("  [", info, "]") else ""))
}
check <- function(test, ok, info = "") record(test, if (isTRUE(ok)) "PASS" else "FAIL", info)
check_error <- function(test, expr) record(test, if (inherits(try(expr, silent = TRUE), "try-error")) "PASS" else "FAIL", "error expected")
skip <- function(test, reason) record(test, "SKIP", reason)
fmt <- function(x, d = 3) paste(signif(x, d), collapse = ", ")

files <- list.files(file.path(root, "tests"), pattern = "^test_.*\\.R$", full.names = TRUE)
if (length(pattern)) files <- files[grepl(paste(pattern, collapse = "|"), basename(files))]
t_all <- Sys.time()
for (f in files) {
  .file <- basename(f); cat(sprintf("\n== %s\n", .file)); t0 <- Sys.time()
  r <- try(sys.source(f, envir = new.env()), silent = TRUE)
  if (inherits(r, "try-error")) record("file ran without error", "FAIL", conditionMessage(attr(r, "condition")))
  cat(sprintf("  (%.1f s)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
dir.create(file.path(root, "tests", "output"), showWarnings = FALSE)
write.csv(.results, file.path(root, "tests", "output", "test_results.csv"), row.names = FALSE)
tab <- table(factor(.results$status, c("PASS", "FAIL", "SKIP")))
cat(sprintf("\n%d PASS, %d FAIL, %d SKIP (%.1f min%s)\n", tab[["PASS"]], tab[["FAIL"]], tab[["SKIP"]],
            as.numeric(difftime(Sys.time(), t_all, units = "mins")), if (SLOW) ", with slow tests" else ""))
if (tab[["FAIL"]] > 0) quit(status = 1)
