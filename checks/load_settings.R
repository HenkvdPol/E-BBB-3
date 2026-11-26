## Helper, sourced by the check scripts: evaluates Setting.R up to the RUN part (no simulation).
## Evaluate Setting.R up to the RUN part (no simulation), with paths in basket_multi
dir <- dirname(check_dir)          # the package folder (parent of checks/)
txt <- readLines(file.path(dir, "Setting.R"))
eval(parse(text = txt[1:(grep("^##  RUN:", txt) - 2)]))
settings$evalue_library <- file.path(dir, "design_optimal_evalues.R")
settings$sim$out_dir <- file.path(dir, "Result"); settings$sim$bank_dir <- file.path(dir, "design_banks")
source(settings$evalue_library); source(file.path(dir, "additional_functions.R"))
