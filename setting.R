################################################################################
# Simulation settings
################################################################################

# Hypotheses per basket: H0: p <= p0. The reference alternative p1 centres the NEX
# prior at (p0 + p1) / 2, defines the scenarios and, in the open-ended design, the
# futility hypothesis p >= p1.
p0 <- 0.10
p1 <- 0.30

# Type I error per basket: a basket is rejected when its e-value reaches 1/alpha = 20.
alpha <- 0.05

# Type II error bound of the futility e-value (open-ended design only):
# a basket stops for futility when its futility e-value reaches 1/beta = 5.
beta <- 0.20

# Maximum number of patients per basket in the fixed design.
n_max <- 30

# Open-ended design (no maximum): the simulation stops a basket at n_cap patients
# only to end the computation; it is then reported as 'open'.
n_cap <- 150

# Time between two analyses. One time unit = one new patient per basket at rate 1.
look_every <- 5

# Grid of baskets: J tumour types (rows) x L biomarker groups (columns).
J <- 3
L <- 3

# Borrowing models (priors of the e-value) that are compared.
Models <- c('none', 'one-way EXNEX', 'one-way EXNEX, scaled Beta', 'two-way EXNEX')

# Recruitment designs. Per basket in the order (1,1), (1,2), (1,3), (2,1), ..., (3,3):
# rate = patients per time unit; entry = opening time; NA = cell not in the trial.
Recruitment <- list(
  'balanced'   = list(rate  = c(1, 1, 1,     1, 1, 1,     1, 1, 1),
                      entry = c(0, 0, 0,     0, 0, 0,     0, 0, 0)),
  'unbalanced' = list(rate  = c(1, 1, 0.5,   1, 0.5, 0.5, 0.5, 0.5, NA),
                      entry = c(0, 0, 0,     0, 0, 20,    0, 0, NA)))

# Scenarios: true response rate per basket, in the same order (rows = tumour types).
Scenario <- list(
  'global null'                 = c(p0, p0, p0,    p0, p0, p0,    p0, p0, p0),
  'global alternative'          = c(p1, p1, p1,    p1, p1, p1,    p1, p1, p1),
  'active tumour type (row 1)'  = c(p1, p1, p1,    p0, p0, p0,    p0, p0, p0),
  'active biomarker (column 1)' = c(p1, p0, p0,    p1, p0, p0,    p1, p0, p0),
  'one active basket'           = c(p1, p0, p0,    p0, p0, p0,    p0, p0, p0),
  'row and column active'       = c(p1, p1, p1,    p1, p0, p0,    p1, p0, p0),
  # logit p = logit p0 + row effect + column effect, effects 0, 0.5 and 1 (0.10 to 0.45)
  'additive row + column'       = as.vector(t(plogis(qlogis(p0) + outer(c(0, 0.5, 1), c(0, 0.5, 1), '+')))),
  'interaction (diagonal)'      = c(p1, p0, p0,    p0, p1, p0,    p0, p0, p1),
  'mixed active and null'       = c(p1, p1, p0,    p1, p0, p0,    p0, p0, p1),
  'strong heterogeneity'        = c(0.05, p0, 0.45, p0, p1, p0,   0.60, p0, 0.05))

# Two-way EXNEX sampler per data state: burn-in, kept draws, chains, and the number of
# data states sampled together (accuracy: Additional_files/docs/two-way-validation.md).
mcmc <- list(burn = 500, keep = 1000, chains = 2, chunk = 150)

# Number of Monte Carlo repeats (trials per simulation cell).
R <- 1000

# Use parallel computing or not.
use_cluster <- TRUE

# Starting seed for reproducibility.
startSeed0 <- 20261126
