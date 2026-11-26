################################################################################
# Simulation settings
################################################################################

# Hypotheses per basket: H0: p <= p0. The reference alternative p1 is the response
# rate of an active basket in the scenarios, the oracle's bet for inactive baskets
# and, in the open-ended design, the futility hypothesis p >= p1.
p0 <- 0.10
p1 <- 0.30

# Type I error per basket: a basket is rejected when its e-value reaches 1/alpha = 20.
alpha <- 0.05

# Type II error bound of the futility e-value (open-ended design only):
# a basket stops for futility when its futility e-value reaches 1/beta = 5.
beta <- 0.20

# Maximum number of patients per basket in the fixed design. The fixed design stops
# a basket at efficacy (E >= 1/alpha) or at n_max, as in eBBB.
n_max <- 30

# Open-ended design (no maximum): the simulation stops a basket at n_cap patients
# only to end the computation; it is then reported as 'open'.
n_cap <- 150

# Time between two analyses. One time unit = one new patient per basket at rate 1.
look_every <- 5

# Grid of baskets, as in eBBB: J tumour types (rows) x L biomarker groups (columns).
J <- 4
L <- 3

# Borrowing models (priors of the e-value) that are compared: the priors of eBBB
# (LR_evalue_multiple) and our two-way EX/NEX in place of its Stan fit.
Models <- c('no borrowing', 'power, equal', 'power, column', 'EX/NEX', 'similarity',
            'two-way EX/NEX', 'oracle')

# Settings of the eBBB priors.
# Initial Beta prior (before restriction to (p0, 1]).
a0 <- 1
b0 <- 1
# Power prior weights of the baskets in the same row (tumour type) and column (biomarker).
w_equal <- c(row = 0.5, col = 0.5)   # 'power, equal'
w_col   <- c(row = 0,   col = 1)     # 'power, column' and the EX component of 'EX/NEX'
# Prior probability that a basket is exchangeable ('EX/NEX' and 'two-way EX/NEX').
p_ex <- 0.5
# Similarity prior: weight (1 - JSD)^epsilon, set to 0 below tau.
epsilon <- 2
tau <- 0.3

# Two-way EX/NEX hyperparameters, as in eBBB (exnex.stan):
# mu ~ N(logit p0, mu_sd^2); sigma_a ~ HN(scale_a); sigma_b ~ HN(scale_b);
# sigma_eps ~ HN(scale_e); NEX: theta ~ N(logit p0, nex_sd^2).
twoway_hyper <- list(mu_sd = 2, scale_a = 1, scale_b = 1, scale_e = 0.5, nex_sd = 2)

# Recruitment designs. Per basket in the order (1,1), (1,2), (1,3), (2,1), ..., (4,3)
# (row = tumour type, column = biomarker group):
# rate = patients per time unit; entry = opening time; NA = cell not in the trial.
# 'unbalanced': accrual falls with the tumour type and the biomarker group, as in eBBB,
# cell (3,3) is not in the trial and basket (2,3) opens at time 20.
Recruitment <- list(
  'balanced'   = list(rate  = c(1, 1, 1,      1, 1, 1,      1, 1, 1,      1, 1, 1),
                      entry = c(0, 0, 0,      0, 0, 0,      0, 0, 0,      0, 0, 0)),
  'unbalanced' = list(rate  = c(1, 1, 0.5,    1, 1, 0.5,    0.5, 0.5, NA, 0.5, 0.5, 0.25),
                      entry = c(0, 0, 0,      0, 0, 20,     0, 0, NA,     0, 0, 0)))

# Scenarios: true response rate per basket, in the same order (4 rows of 3).
# The first four are the scenarios of eBBB.
Scenario <- list(
  'global null'          = c(p0, p0, p0,     p0, p0, p0,     p0, p0, p0,     p0, p0, p0),
  'biomarker 1 active'   = c(p1, p0, p0,     p1, p0, p0,     p1, p0, p0,     p1, p0, p0),
  'interaction'          = c(p1, p0, p0,     p1, p0, p0,     p1, p0, p0,     p0, p0, p0),  # biomarker 1, not tumour type 4
  'all active'           = c(p1, p1, p1,     p1, p1, p1,     p1, p1, p1,     p1, p1, p1),
  'tumour type 1 active' = c(p1, p1, p1,     p0, p0, p0,     p0, p0, p0,     p0, p0, p0),
  'one active basket'    = c(p1, p0, p0,     p0, p0, p0,     p0, p0, p0,     p0, p0, p0),
  'tumour type 1 and biomarker 1 active'
                         = c(p1, p1, p1,     p1, p0, p0,     p1, p0, p0,     p1, p0, p0),
  # logit p = logit p0 + tumour effect (0, 1/3, 2/3, 1) + biomarker effect (0, 0.5, 1)
  'additive tumour + biomarker'
                         = as.vector(t(plogis(qlogis(p0) + outer(c(0, 1, 2, 3) / 3, c(0, 0.5, 1), '+')))),
  'diagonal interaction' = c(p1, p0, p0,     p0, p1, p0,     p0, p0, p1,     p0, p0, p0),
  'mixed active and null'= c(p1, p1, p0,     p1, p0, p0,     p0, p0, p1,     p0, p1, p0),
  'strong heterogeneity' = c(0.05, p0, 0.45, p0, p1, p0,     0.60, p0, 0.05, p1, 0.05, p0))

# Keep the intended null basket exactly at p0: plogis(qlogis(p0)) can round above it.
Scenario[['additive tumour + biomarker']][1] <- p0

# Two-way EX/NEX sampler per data state: burn-in, kept draws, chains, and the number of
# data states sampled together (accuracy: Additional_files/docs/two-way-validation.md).
mcmc <- list(burn = 500, keep = 1000, chains = 2, chunk = 150)

# Number of Monte Carlo repeats (trials per simulation cell).
R <- 1000

# Use parallel computing or not.
use_cluster <- TRUE

# Starting seed for reproducibility.
startSeed0 <- 20261126
