## Run from anywhere:  Rscript checks/unit_flex.R   (paths are relative to this file)
check_dir <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
out_dir <- file.path(check_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
source(file.path(check_dir, "load_settings.R"))
runs <- complete_runs(settings)
handmade <- function(S, K, N) array(((outer(outer(1:S * 0.7548777, 1:K * 0.5698403, "+"), 1:N * 0.3819660, "+")) %% 1), c(S, K, N))
bad <- character(); chk <- function(ok, msg) if (!isTRUE(all(ok))) bad <<- c(bad, msg)
for (rn in c("grid_p10-30_K5_n30", "simon_p05-25_K5_n17", "simon_p70-90_K5_n27")) {
  s <- complete_settings(settings, runs[runs$name == rn, ]); eng <- build_engine(s); banks <- get_banks(s)
  U <- handmade(30, s$K, max(s$n_max + s$n_ext))
  for (su in names(s$setups)) { setup <- s$setups[[su]]
    for (sn in names(s$schedules)) { block <- s$schedules[[sn]]
      for (scn in names(s$scenarios)) { C <- make_data(U, s$scenarios[[scn]]); R <- list()
        for (dn in setup$designs) { d <- s$designs[[dn]]
          if (d$type == "mixture") { r <- run_trials_mixture(d, C, block, setup$entry, s, eng, new_cache()); tag <- paste(rn, su, sn, scn, dn)
            ne <- if (isTRUE(d$extend)) s$n_ext[1] else 0
            chk(all(r$status %in% 1:3) && all(r$n <= s$n_max[1] + ne) && !any(r$n > s$n_max[1] & !r$extended), paste(tag, "mixture: n rules"))
            chk(if (isTRUE(d$accelerate)) all(r$time <= r$n + 1e-9) else all(abs(r$time - r$n) < 1e-9), paste(tag, "mixture: time rules"))
            chk(all((r$status == 1) == (r$e_value >= 20)), paste(tag, "mixture: rejected iff E >= 1/alpha"))
            chk(all(r$theta1 > s$p0), paste(tag, "mixture prior mean not above p0")); next }
          if (d$type != "evalue") next
          R[[dn]] <- run_trials_evalue(d, C, block, setup$entry, s, eng, new_cache(), banks_for(d, block, setup, s, banks),
                                       if (isTRUE(d$extend)) ext_banks_for(d, block, setup, s, banks))
          r <- R[[dn]]; tag <- paste(rn, su, sn, scn, dn)
          chk(all(r$status %in% 1:3) && !anyNA(r$time), paste(tag, "undecided cohort"))
          chk(r$n <= s$n_max[1] + if (isTRUE(d$extend)) s$n_ext[1] else 0, paste(tag, "n above maximum"))
          chk(!(r$n > s$n_max[1] & !r$extended), paste(tag, "beyond n_max without extension"))
          chk(!r$extended | isTRUE(d$extend), paste(tag, "extended without extend"))
          if (!isTRUE(d$accelerate)) chk(abs(r$time - r$n) < 1e-9, paste(tag, "time != n at normal speed"))
          else chk(r$time <= r$n + 1e-9, paste(tag, "accelerated slower than normal"))
        }
        if (su == "basket") {                                 # same bets: speeding up never loses a rejection
          a <- R[["e-based (accelerated)"]]; f <- R[["e-based"]]; both <- a$status == 1 & f$status == 1
          chk(a$status[f$status == 1] == 1, paste(rn, sn, scn, "accelerated lost a rejection"))
          chk(a$time[both] <= f$time[both] + 1e-9, paste(rn, sn, scn, "accelerated rejected later"))
        }
      }
    }
  }
  ## the whole pipeline of one run on hand-made data: every design, summaries and all figures
  refs <- references(s); data <- list(eval = lapply(s$scenarios, function(p) make_data(U, p)), cal = make_data(U[, , 1:s$n_max[1]], s$scenarios[[1]]))
  out <- list()
  for (su in names(s$setups)) {
    res <- unlist(lapply(s$setups[[su]]$designs, run_design, setup = s$setups[[su]], settings = s, eng = eng, banks = banks, data = data), recursive = FALSE)
    trial <- do.call(rbind, lapply(res, `[[`, "trial")); curves <- do.call(rbind, lapply(res, `[[`, "curves"))
    for (v in c("typeI", "power", "early_reject_H1", "extend_H0", "extend_H1", "reject_ext_H1")) chk(is.na(trial[[v]]) | (trial[[v]] >= 0 & trial[[v]] <= 1), paste(rn, su, v, "outside [0,1]"))
    chk(is.na(trial$reject_ext_H1) | trial$reject_ext_H1 <= trial$power + 1e-12, paste(rn, su, "rejected after n_max > power"))
    fixed <- !grepl("flexible|extended", trial$design); chk(trial$extend_H1[fixed] %in% c(0, NA), paste(rn, su, "fixed design extended"))
    p <- list(plot_cumulative(curves, su, s, refs), plot_scenarios(trial, su, s), plot_time(curves, su, s),
              plot_tradeoff(trial, su, s, refs), plot_scenarios(trial, su, s, meas_decisions, "decisions"))
    na <- function(x) is.na(x)
    chk(na(trial$FDR) | (trial$FDR >= 0 & trial$FDR <= 1), paste(rn, su, "FDR outside [0,1]"))
    chk(na(trial$power_all) | trial$power_all <= trial$power + 1e-12, paste(rn, su, "all-active power > power"))
    chk(na(trial$power_any) | trial$power_any >= trial$power - 1e-12, paste(rn, su, "any-active power < power"))
    chk(na(trial$false_rejections) | na(trial$FWER) | trial$false_rejections >= trial$FWER - 1e-12, paste(rn, su, "false rejections < FWER"))
    chk(na(trial$n_q90_H0) | na(trial$n_q50_H0) | trial$n_q90_H0 >= trial$n_q50_H0, paste(rn, su, "quantiles"))
    chk(is.finite(trial$trial_n) & trial$trial_n > 0, paste(rn, su, "trial n"))
    chk(is.na(trial$log10E_H1) == (trial$design %in% c("Bayes Umbrella", "Bayes mEXNEX") | is.na(trial$power)), paste(rn, su, "e-value at decision"))
    eff <- design_effects(trial, s, s$setups[[su]]$designs)
    chk(all(c("borrowing", "flexible inclusion", "vs primary design") %in% eff$effect), paste(rn, su, "effects missing"))
    e1 <- eff[eff$effect == "vs primary design", ][1, ]; ref <- trial[trial$design == s$primary & trial$schedule == e1$schedule & trial$scenario == e1$scenario, ]
    me <- trial[trial$design == e1$design & trial$schedule == e1$schedule & trial$scenario == e1$scenario, ]
    chk(isTRUE(all.equal(e1$d_power, me$power - ref$power)), paste(rn, su, "effect arithmetic"))
    if (rn == "grid_p10-30_K5_n30" && su == "basket") { save_plot(p[[4]], file.path(out_dir, "fig4_basket"), 11, 10); save_plot(p[[5]], file.path(out_dir, "fig5_basket"), 11, 15)
      print(format(trial[trial$schedule == names(s$schedules)[2] & trial$scenario == "mixed",
        c("design", "true_rejections", "false_rejections", "FDR", "all_correct", "trial_n", "bias_reject_H1", "log10E_H1")], digits = 2), row.names = FALSE)
      print(eff[eff$scenario == "mixed" & eff$schedule == names(s$schedules)[2], c("design", "compared_with", "effect", "d_power", "d_ESS_H0", "d_time_reject_H1")], row.names = FALSE) }
    for (q in p) ggplot2::ggsave(file.path(out_dir, "q.png"), q, width = 11, height = 12)
    out[[su]] <- cbind(setup = su, series = s$series, p0 = s$p0, p1 = s$p1, K = s$K, n_max = max(s$n_max),
                       schedule_pos = match(trial$schedule, names(s$schedules)), is_null = trial$scenario == s$scenario_null,
                       is_alt = trial$scenario == s$scenario_alt, trial)
    if (rn == "grid_p10-30_K5_n30") save_plot(plot_time(curves, su, s), paste0(file.path(out_dir, "fig3_"), su), 11, 8)
  }
  out$reference <- simon_rows(out[[1]][1:2, ], refs, s); all <- do.call(rbind, out)
  for (su in names(s$setups)) save_plot(plot_overview(all, all$series[1], su, s), paste0(file.path(out_dir, "ov_"), rn, "_", su), 11, 15)
  cat(rn, ": pipeline ok, rows", nrow(all), "\n")
  if (rn == "grid_p10-30_K5_n30") print(format(all[all$setup == "basket" & all$schedule_pos == 2 & all$is_alt,
     c("design", "power", "time_reject_H1", "n_reject_H1", "extend_H1", "reject_ext_H1", "duration")], digits = 3), row.names = FALSE)
}
cat("\nproblems:", length(bad), "\n"); if (length(bad)) writeLines(head(bad, 20))
