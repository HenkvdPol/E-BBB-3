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
setting.R               All settings: p0, p1, alpha, beta, n_max, the 4 x 3 grid of
                        baskets, recruitment designs, scenarios, borrowing models (priors
                        of eBBB and the two-way EX/NEX) and their parameters, R, seeds.
simulate.R              Runs the simulation. Every simulation cell (recruitment design x
                        scenario x borrowing model) is one call of genSim(), saved as
                        Result/<horizon>/sim-<R>-<seed>-<id>.Rdata. When all cells are
                        there, it writes baskets.csv, platform.csv, timing.csv, the
                        figures and session.txt to Result/<horizon>/.
additional_functions.R  All functions used by simulate.R:
                          Part 1  the e-value (priors on (p0, 1], block factor, futility);
                          Part 2  borrowing models: the priors of eBBB (no borrowing,
                                  power, EX/NEX, similarity, oracle) and the two-way
                                  EX/NEX; the earlier one-way EXNEX as an option;
                          Part 3  the simulation of a trial in calendar time;
                          Part 4  summaries and figures.
setup.pdf               Three-page set-up of the study.
Additional_files/       Everything else:
  framework.tex/.pdf      the method and its proofs (lead document);
  setup.tex               source of setup.pdf;
  docs/                   design notes: alignment with eBBB, open-ended design,
                          two-way model and its validation, known limitations,
                          validation plan, review, entry-time proposal,
                          implementation plan;
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

horizon 'fixed': a basket stops when it is rejected (E >= 1/alpha) or at n_max = 30
                 patients, as in eBBB.
horizon 'open' : no maximum sample size; a basket stops when it is rejected or when the
                 futility e-value (testing p >= p1) reaches 1/beta. n_cap = 150 only ends
                 the computation; baskets still recruiting there are reported as 'open'.

On the cluster: mkdir -p logs; sbatch Additional_files/slurm/sim.slurm
(two array tasks: the fixed and the open-ended design, 1000 trials per cell, 8 cores).
The full study has not been run yet: approve the compute first. Timing pilot: about 4.4 (fixed)
and 5.8 (open) CPU hours at R = 1000; under 1 h per horizon on 8 cores is provisional
(Additional_files/docs/timing-pilot.md).

Tests and validation:
Rscript Additional_files/tests/run_tests.R          all tests (about 2 min)
Rscript Additional_files/tests/run_tests.R --slow   plus the slow reference comparisons
Rscript Additional_files/validation/open_ended_feasibility.R   exact open-ended results (1 min)
Rscript Additional_files/validation/twoway_validation.R        two-way engine vs references (30 min)
Rscript Additional_files/validation/twoway_ebbb_hyper.R        two-way, eBBB hyperparameters, vs JAGS
Rscript Additional_files/validation/timing_pilot.R             run time of the study (20 trials per cell)

Packages: Matrix, foreach, doParallel, ggplot2; rjags (with JAGS) only for the JAGS reference.


Results so far (exact calculations and tests; no full simulation)
------------------------------------------------------------------
- The per-basket type I error bound follows from the framework; exact enumeration
  checks it in selected finite configurations for every borrowing model.
- The priors of eBBB give the same block factors as eBBB's closed forms within 3e-4
  (relative) in the tested data states; the EX/NEX weights and the oracle are exact
  (Additional_files/docs/ebbb-alignment.md).
- In the earlier single-basket feasibility calculations with fixed starting priors,
  e-value futility stops every basket by 125 patients, with E[N] 20-26 under p0 and
  power 0.90-0.95 at p1, against 0.49-0.62 for n_max = 30. These are results for those
  configurations; the current 4 x 3 platform study remains to be run
  (Additional_files/docs/open-ended-design.md).
- The two-way engine has been compared with quadrature and JAGS. For the current
  hyperparameters in six tested states, production versus JAGS restricted-mean
  differences reached 0.0077 and 5/5 block-factor relative differences reached 17%.
  Numerical refinement remains part of validation (Additional_files/docs/two-way-validation.md).
- Open points: Additional_files/docs/known-limitations.md.
