## Run from anywhere:  Rscript checks/check_results.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
## Checks of all output of the 1-trial run
source(file.path(check_dir, "load_settings.R"))
R <- settings$sim$out_dir; runs <- complete_runs(settings)
x <- readRDS(file.path(R, "all_results.rds")); all <- x$results
csv <- read.csv(file.path(R, "all_settings_results.csv"))
problems <- character(); bad <- function(cond, msg) if (any(cond, na.rm = TRUE)) problems <<- c(problems, sprintf("%s (%d rows)", msg, sum(cond, na.rm = TRUE)))
cat("rows in all_settings_results.csv:", nrow(csv), " in all_results.rds:", nrow(all), "\n")

## 1. files
top <- c("all_settings_results.csv", "all_results.rds", as.vector(outer(c("overview_grid_basket", "overview_grid_umbrella",
         "overview_simon_basket", "overview_simon_umbrella"), c(".pdf", ".png"), paste0)))
per <- c("basket_results.csv", "umbrella_results.csv", "basket_curves.csv", "umbrella_curves.csv", "references.csv",
         as.vector(outer(c("basket_fig1_cumulative", "basket_fig2_scenarios", "umbrella_fig1_cumulative", "umbrella_fig2_scenarios"), c(".pdf", ".png"), paste0)))
miss <- c(top[!file.exists(file.path(R, top))], unlist(lapply(runs$name, function(r) file.path(r, per)[!file.exists(file.path(R, r, per))])))
cat("folders:", sum(dir.exists(file.path(R, runs$name))), "of", nrow(runs), "; missing files:", length(miss), if (length(miss)) paste(miss, collapse = ", "), "\n")

## 2. rows per run: designs x schedules x scenarios (+ 2 Simon rows)
for (i in seq_len(nrow(runs))) {
  s <- complete_settings(settings, runs[i, ]); a <- all[all$series == runs$series[i] & all$p0 == runs$p0[i] & all$K == runs$K[i] & all$n_max == runs$n_max[i], ]
  want <- 3 * length(s$scenarios) * (length(s$setups$basket$designs) + length(s$setups$umbrella$designs)) + 2
  bad(nrow(a) != want, sprintf("%s: %d rows instead of %d", runs$name[i], nrow(a), want))
  bad(!setequal(unique(a$scenario[a$setup != "reference"]), names(s$scenarios)), paste(runs$name[i], "scenarios differ"))
  ## Simon reference = exact Simon design with n = n_max (in the simon series: Simon's published optimal design)
  ref <- read.csv(file.path(R, runs$name[i], "references.csv")); sd <- simon_design(s$p0, s$p1, s$alpha, s$beta, n = runs$n_max[i])
  bad(abs(ref$reject[grepl("Simon", ref$design) & abs(ref$p - s$p0) < 1e-9] - simon_reject(sd$n1, sd$r1, sd$n, sd$r, s$p0)) > 1e-12, paste(runs$name[i], "Simon type I"))
  bad(abs(a$power[a$setup == "reference" & a$is_alt] - sd$power) > 1e-12, paste(runs$name[i], "Simon power in overview rows"))
  ## theta1 within the borrowing range, = p1 without borrowing
  ev <- a[a$setup != "reference", ]; bor <- vapply(ev$design, function(d) isTRUE(s$designs[[d]]$borrow), TRUE)
  th <- c(ev$theta1_H0[bor], ev$theta1_H1[bor])
  bad(any(th < s$borrow$theta1_min - 1e-9 | th > s$borrow$theta1_max + 1e-9, na.rm = TRUE), paste(runs$name[i], "theta1 outside range"))
  isev <- vapply(ev$design, function(d) s$designs[[d]]$type == "evalue", TRUE)
  bad(abs(c(ev$theta1_H0[isev & !bor], ev$theta1_H1[isev & !bor]) - s$p1) > 1e-9, paste(runs$name[i], "theta1 != p1 without borrowing"))
  bad(!is.na(c(ev$theta1_H0[!isev], ev$theta1_H1[!isev])), paste(runs$name[i], "theta1 for Bayesian design"))
  ## fully sequential: designs without borrowing give the same results in the basket and umbrella set-up
  for (d in intersect(s$setups$basket$designs, s$setups$umbrella$designs)) if (isev[match(d, ev$design)] && !isTRUE(s$designs[[d]]$borrow)) {
    b <- ev[ev$setup == "basket" & ev$design == d & ev$schedule_pos == 1, ]; u <- ev[ev$setup == "umbrella" & ev$design == d & ev$schedule_pos == 1, ]
    cols <- c("typeI", "power", "ESS_H0", "ESS_H1", "n_reject_H1", "early_reject_H1")
    bad(!isTRUE(all.equal(b[, cols], u[, cols], check.attributes = FALSE)), paste(runs$name[i], d, "basket != umbrella (fully sequential)"))
  }
  ## curves: in [0, 1] and non-decreasing
  for (su in c("basket", "umbrella")) { cu <- read.csv(file.path(R, runs$name[i], paste0(su, "_curves.csv")))
    bad(cu$value < -1e-12 | cu$value > 1 + 1e-12, paste(runs$name[i], su, "curve outside [0, 1]"))
    dec <- unlist(tapply(cu$value, interaction(cu$design, cu$schedule, cu$scenario, cu$group, cu$measure, drop = TRUE), function(v) any(diff(v) < -1e-12)))
    bad(dec, paste(runs$name[i], su, "curve decreasing")) }
}

## 3. values in every row
pr <- c("typeI", "typeI_max", "power", "FWER", "early_reject_H1", "futility_H0", "futility_H1")
for (v in pr) bad(all[[v]] < 0 | all[[v]] > 1, paste(v, "outside [0, 1]"))
bad(all$early_reject_H1 > all$power + 1e-12, "early rejection > power")
bad(all$typeI > all$typeI_max + 1e-12, "mean type I > worst type I")
for (v in c("ESS_H0", "ESS_H1", "n_reject_H1", "n_early_reject_H1")) bad(all[[v]] < 1 | all[[v]] > all$n_max, paste(v, "outside [1, n_max]"))
bad(all$n_early_reject_H1 >= all$n_max, "early rejection at n_max")
bad(!is.na(all$n_early_reject_H1) & is.na(all$early_reject_H1), "n of early rejection without early rejection")
nul0 <- grepl("^0 active", all$scenario); altK <- grepl("global alt|strong", all$scenario)
bad(!is.na(all$power[nul0]), "power in global null"); bad(!is.na(all$typeI[altK]), "type I without null cohorts")
bad(is.na(all$typeI[nul0]), "no type I in global null"); bad(is.na(all$power[altK]), "no power in global alternative")
cat("\nproblems found:", length(problems), "\n"); if (length(problems)) writeLines(paste(" -", problems))

## 4. summary of the 1-trial results (meaningless statistically, shows that every design produces numbers)
ev <- all[all$setup != "reference", ]
cat("\nNA per column (rows with null and active cohorts expected to have NA in some columns):\n")
print(colSums(is.na(ev[, c("typeI", "power", "ESS_H0", "ESS_H1", "n_reject_H1", "early_reject_H1", "n_early_reject_H1", "theta1_H1")])))
cat("\nmean over all rows by design (1 trial!):\n")
print(aggregate(cbind(typeI, power, ESS_H0, ESS_H1, n_reject_H1) ~ design, ev, mean, na.action = na.pass), digits = 2)
