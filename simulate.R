################################################################################
# Simulation of the e-value test with Bayesian borrowing (framework.tex)
################################################################################

# Run from this folder:
#   Rscript simulate.R                         test setting: 2 trials per cell
#   Rscript simulate.R fixed 1000 20261126 TRUE  (horizon, R, start seed, use_cluster)
# On the cluster: Additional_files/slurm/sim.slurm (one line of parameters.txt per task).

source('setting.R')

if(length(commandArgs(trailingOnly = TRUE)) == 0){
  # Test setting when in R and don't want to run the full simulation.
  args <- c('fixed', 2, startSeed0, TRUE)
}else{
  args <- commandArgs(trailingOnly = TRUE)
}

cat(args)

horizon <- args[1]                 # 'fixed': at most n_max per basket; 'open': no maximum
R <- as.integer(args[2])           # Monte Carlo repeats (trials per simulation cell)
startSeed <- as.integer(args[3])   # seed for reproducibility
use_cluster <- as.logical(args[4]) # parallel computing or not

cat('\n', horizon, R, startSeed, use_cluster)

# All simulation cells: recruitment design x scenario x borrowing model
sim <- expand.grid('Recruitment' = names(Recruitment),
                   'Scenario' = names(Scenario),
                   'Model' = Models,
                   stringsAsFactors = FALSE)
sim$id <- 1:nrow(sim)

################################################################################
# Libraries and additional functions to create simulation.
################################################################################

# Necessary libraries
library(Matrix)
library(ggplot2)
library(foreach)
library(doParallel)

# Additional functions in use.
source('additional_functions.R')

dir.create(file.path('Result', horizon), showWarnings = FALSE, recursive = TRUE)


################################################################################
# Simulation function
################################################################################

genSim <- function(R, recruitment, scenario, model, horizon, startSeed, id){
  # R = number of Monte Carlo repeats (trials, simulated together)
  # recruitment = name of the recruitment design (Recruitment in setting.R)
  # scenario = name of the scenario (Scenario in setting.R)
  # model = borrowing model (Models in setting.R)
  # horizon = 'fixed' or 'open'
  # startSeed = seed for reproducibility
  # id = the simulation number
  cat('\n id =', id, ':', recruitment, '|', scenario, '|', model)
  start <- proc.time()

  # True response rates of the baskets in the trial
  p_true <- Scenario[[scenario]][!is.na(Recruitment[[recruitment]]$rate)]

  # The design and the borrowing model of this cell (the oracle needs the true rates)
  design <- make_design(Recruitment[[recruitment]], model, horizon)
  engine <- make_engine(model, design$layout, p_true)

  # Seeds. The patient outcomes only depend on the recruitment design and the scenario,
  # so all borrowing models see the same patients. The borrowing model (the two-way
  # sampler) has its own seeds, which never overlap with those of the outcomes.
  seed_outcome <- startSeed + 1000 * match(scenario, names(Scenario)) +
    match(recruitment, names(Recruitment))
  seed_fit <- startSeed + 10^6 + 10^4 * id

  Result <- simulate_trials(p_true, design, engine, R, seed = seed_outcome, fit_seed = seed_fit)

  # We also add the parameters here to store the information.
  Result$recruitment <- recruitment
  Result$scenario <- scenario
  Result$engine <- model
  Result$horizon <- horizon
  Result$Simulation <- id
  Result$seconds <- (proc.time() - start)[['elapsed']]

  # We save the results for every simulation cell.
  save(Result, file = paste0('Result/', horizon, '/sim-', R, '-', startSeed, '-', id, '.Rdata'))

  # Return the id.
  return(id)
}

if(use_cluster){
  message("SLURM_CPUS_PER_TASK: ", Sys.getenv("SLURM_CPUS_PER_TASK"))
  message("SLURM_JOB_ID: ", Sys.getenv("SLURM_JOB_ID"))
  message("detectCores: ", parallel::detectCores())

  # Want to know the time
  start <- proc.time()
  # Workers: the cores of the slurm job, otherwise all cores of this computer
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
  registerDoParallel(cores = n_cores)

  message("Registered parallel backend: ", getDoParName())
  message("Number of workers: ", getDoParWorkers())

  # The cells differ a lot in run time (two-way EXNEX is slowest), so every worker
  # takes the next cell when it is done (no prescheduling).
  foreach(i = 1:nrow(sim),
          .options.multicore = list(preschedule = FALSE)) %dopar% {genSim(R,
                                                                          sim$Recruitment[i],
                                                                          sim$Scenario[i],
                                                                          sim$Model[i],
                                                                          horizon,
                                                                          startSeed,
                                                                          i)}

  # Stop clusters.
  stopImplicitCluster()
  print(proc.time() - start)
}

if(!use_cluster){
  start <- proc.time()

  for(i in 1:nrow(sim)){
    genSim(R,
           sim$Recruitment[i],
           sim$Scenario[i],
           sim$Model[i],
           horizon,
           startSeed,
           i)
  }
  print(proc.time() - start)
}


################################################################################
# Combine all cells: summaries and figures
################################################################################

files <- paste0('Result/', horizon, '/sim-', R, '-', startSeed, '-', sim$id, '.Rdata')

if(all(file.exists(files))){
  all_results <- list()
  for(i in seq_along(files)){
    load(files[i])
    all_results[[i]] <- Result
  }
  d <- do.call(rbind, all_results)
  d$scenario <- factor(d$scenario, levels = names(Scenario))
  d$engine <- factor(d$engine, levels = Models)

  # Per basket and per platform. Per-basket type I error and FWER are separate columns.
  baskets <- summarise_baskets(d)
  platform <- summarise_platform(d)
  write.csv(baskets, file.path('Result', horizon, 'baskets.csv'), row.names = FALSE)
  write.csv(platform, file.path('Result', horizon, 'platform.csv'), row.names = FALSE)

  # Run time per cell
  timing <- unique(d[, c('Simulation', 'recruitment', 'scenario', 'engine', 'seconds')])
  write.csv(timing, file.path('Result', horizon, 'timing.csv'), row.names = FALSE)

  # Figures per recruitment design
  for(r in names(Recruitment)){
    save_plot(plot_baskets(baskets, r, alpha), file.path('Result', horizon, paste0('fig_baskets_', r)), 11, 14)
    save_plot(plot_tradeoff(platform, r), file.path('Result', horizon, paste0('fig_tradeoff_', r)), 11, 8)
  }

  # What was run, for reproducibility
  sha <- tryCatch(system2('git', c('rev-parse', 'HEAD'), stdout = TRUE), error = function(e){NA})
  writeLines(c(paste('git', sha), paste('arguments:', paste(args, collapse = ' ')),
               capture.output(sessionInfo())),
             file.path('Result', horizon, 'session.txt'))
  cat('\n Summaries and figures in', file.path('Result', horizon), '\n')
}
