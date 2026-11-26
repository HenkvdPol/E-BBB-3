# Parse every tracked R source file without executing simulations.
files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE,
                    all.files = TRUE)
files <- files[!grepl("(^|/)(Result|design_banks|output)/", files)]
if (!length(files)) stop("No R files found")
for (f in files) {
  tryCatch(parse(file = f), error = function(e) {
    stop(sprintf("Parse failed for %s: %s", f, conditionMessage(e)), call. = FALSE)
  })
}
cat(sprintf("Parsed %d R files successfully.\n", length(files)))
