E-BBB-3: anytime-valid e-values with Bayesian borrowing for basket trials
=========================================================================

The method is described in Additional_files/framework.pdf (the lead document) and
summarised in setup.pdf. Every basket is tested with an e-value: at each analysis a
Bayesian borrowing model is fitted to all baskets, its posterior is restricted to
(p0, 1] and frozen as the prior of the next block of patients. The type I error is at
most alpha per basket, whether or not the borrowing model is right, at any time, and
also without a maximum sample size.


Files
-----
setting.R               All settings: p0, p1, alpha, beta, n_max, the grid of baskets,
                        recruitment designs, scenarios, borrowing models, R, seeds.
simulate.R              Runs the simulation. Every simulation cell (recruitment design x
                        scenario x borrowing model) is one call of genSim(), saved as
                        Result/<horizon>/sim-<R>-<seed>-<id>.Rdata. When all cells are
                        there, it writes baskets.csv, platform.csv, timing.csv, the
                        figures and session.txt to Result/<horizon>/.
additional_functions.R  All functions used by simulate.R:
                          Part 1  the e-value (priors on (p0, 1], block factor, futility);
                          Part 2  borrowing models: none, one-way EXNEX, two-way EXNEX;
                          Part 3  the simulation of a trial in calendar time;
                          Part 4  summaries and figures.
setup.pdf               Three-page set-up of the study.
Additional_files/       Everything else:
  framework.tex/.pdf      the method and its proofs (lead document);
  setup.tex               source of setup.pdf;
  docs/                   design notes: open-ended design, two-way model and its
                          validation, known limitations, validation plan, review,
                          entry-time proposal, implementation plan;
  tests/                  tests of additional_functions.R (run_tests.R, parse_all.R);
  validation/             reference calculations (quadrature, JAGS) and the exact
                          open-ended feasibility analysis, with their results;
  slurm/                  sim.slurm and parameters.txt for the cluster;
  legacy/                 the earlier code (design-optimal e-values, Simon, Bayesian
                          comparators), unchanged;
  PROJECT_HISTORY.md      decisions and results over time;
  CLAUDE_NEXT_STEP.md     instructions for the next development step.


How to run (from this folder)
-----------------------------
Rscript simulate.R                            test setting: 2 trials per cell (fixed design)
Rscript simulate.R open 2 20261126 TRUE       test setting of the open-ended design
Rscript simulate.R fixed 1000 20261126 TRUE   the study (arguments: horizon, R, seed, use_cluster)

horizon 'fixed': at most n_max = 30 patients per basket.
horizon 'open' : no maximum sample size; a basket stops when it is rejected or when the
                 futility e-value (testing p >= p1) reaches 1/beta. n_cap = 150 only ends
                 the computation; baskets still recruiting there are reported as 'open'.

On the cluster: mkdir -p logs; sbatch Additional_files/slurm/sim.slurm
(two array tasks: the fixed and the open-ended design, 1000 trials per cell, 8 cores).
The full study has not been run yet: approve the compute first.

Tests and validation:
Rscript Additional_files/tests/run_tests.R          all tests (about 1.5 min)
Rscript Additional_files/tests/run_tests.R --slow   plus the slow reference comparisons
Rscript Additional_files/validation/open_ended_feasibility.R   exact open-ended results (1 min)
Rscript Additional_files/validation/twoway_validation.R        two-way engine vs references (30 min)

Packages: Matrix, foreach, doParallel, ggplot2; rjags (with JAGS) only for the JAGS reference.


Results so far (exact calculations and tests; no full simulation)
------------------------------------------------------------------
- The type I error is at most alpha for every borrowing model (exact enumeration).
- The open-ended design with e-value futility stops every basket (by 125 patients in the
  exact analysis), with E[N] 20-26 under p0 and power 0.90-0.95 at p1, against power
  0.49-0.62 for the fixed design with n_max = 30 (Additional_files/docs/open-ended-design.md).
- The two-way EXNEX posterior agrees with deterministic quadrature and JAGS within Monte
  Carlo error (Additional_files/docs/two-way-validation.md).
- Open points: Additional_files/docs/known-limitations.md.
