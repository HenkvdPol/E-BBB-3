# E-BBB-3: anytime-valid e-values with Bayesian borrowing for basket and platform trials

Research code for basket and platform trials with a binary response. [`framework.tex`](framework.tex) is the lead
document: it specifies the test and proves its validity. This repository implements that test and evaluates it by
simulation.

## The test (framework.tex)

* For each basket c, test H0,c: p_c <= p0 with a one-sided alpha per basket.
* **At each analysis:**
  * fit a Bayesian borrowing model to all platform data;
  * restrict its posterior of p_c to (p0, 1];
  * freeze it as the prior of the next block of patients.
* **Block factor:** the likelihood ratio of the new block, mixed over that prior.
* **E-value:** the product of the block factors. Reject at the first analysis with E_c >= 1/alpha.
* **Guarantee:** under the conditional null, P(reject basket c) <= alpha at every sample size. This holds under
  optional stopping and continuation, whatever the other baskets show, and whether or not the borrowing model is
  right.
* **What it is not:** per basket, not a platform-wide FWER. Allocation, outcome order and delays matter
  (framework.tex, 'Scope of the guarantee').

## Repository layout

| Path | Content |
|---|---|
| `framework.tex`, `setup.tex` | the statistical proposal, and the research and software set-up |
| `R/eprocess.R` | **validity layer**: grid, restriction to (p0, 1], prior checks, block factor; knows nothing about borrowing |
| `R/borrow_oneway.R` | borrowing engines: none, and one-way EXNEX (the earlier benchmark; regression-tested against `legacy/`) |
| `R/borrow_twoway.R` | two-way EXNEX for a tumour type x biomarker grid (framework.tex eq. 9) |
| `R/reference_twoway.R` | reference calculations for validation: deterministic quadrature, JAGS, equivalent one-way model |
| `R/simulate.R` | the e-process in calendar time; priors fitted at an analysis and frozen for the next block |
| `R/summarise.R`, `R/plots.R` | operating characteristics (per-basket type I error and FWER kept apart) and figures |
| `run_framework_study.R` | **the framework-only simulation study**: the single file to run (settings + run) |
| `tests/` | `run_tests.R` and `test_*.R`: unit, regression and exact validity tests (base R) |
| `validation/` | two-way engine against reference calculations; results in `validation/results/` |
| `legacy/` | the earlier full code base, unchanged: design-optimal e-values (GROW, P-max, ESS-min, e-based), flexible inclusion, Bayesian comparators, Simon, and its checks |
| `docs/` | implementation plan, known limitations, validation plan, two-way specification and validation, entry-time proposal |
| `PROJECT_HISTORY.md` | decisions and results over time |

Generated output (`Result/`, `design_banks/`, `*/output/`) is not version controlled.

## Borrowing models in the study

Every model gives a valid test. They differ in power and borrowing, not in type I error control.

| Name | Model |
|---|---|
| none | NEX prior N(logit pbar, 1/pbar + 1/(1 - pbar)), pbar = (p0 + p1)/2; the basket's own data only |
| one-way EXNEX | all baskets one exchangeable list (Neuenschwander et al. 2016). EX: N(mu, tau²), mu ~ N(logit p0, 10²) on [-8, 4], tau ~ HN(1); w = 0.5 |
| one-way EXNEX, scaled Beta | as above, but the restricted posterior is replaced by the scaled Beta with the same mean and variance (the earlier prototype's prior) |
| two-way EXNEX | EX: theta_jl ~ N(mu + a_j + b_l, sigma_e²), with a_j, b_l row and column effects; sigma_a, sigma_b, sigma_e ~ HN(1/sqrt 3), so the total EX variance matches the one-way benchmark; w = 0.5 per basket; same NEX |

Two-way details, decisions and validation: [`docs/two-way-exnex-spec.md`](docs/two-way-exnex-spec.md) and
[`docs/two-way-validation.md`](docs/two-way-validation.md).

## Running

Needs R >= 4.3. Packages:
* base R and `Matrix` (a recommended package) for the computations;
* `ggplot2` for the figures;
* `rjags` (with JAGS) only for the JAGS reference.

```sh
Rscript tests/run_tests.R                   # all fast tests (~1.5 min)
Rscript tests/run_tests.R --slow            # plus slow reference comparisons
Rscript validation/twoway_validation.R      # two-way engine vs references (~30 min, no trial simulation)
N_SIM=2 Rscript run_framework_study.R       # smoke run of the study (~9 min on 2 cores)
N_SIM=1000 Rscript run_framework_study.R    # the study: only after the compute cost is approved
```

`run_framework_study.R` uses `SLURM_CPUS_PER_TASK` cores when set, otherwise all cores. It writes to
`Result/framework_study/`:
* `baskets.csv`: per basket — rejection rate, mean n and time, posterior bias, 90% coverage, P(EX), restricted
  prior mean, log10 E;
* `platform.csv`: per scenario — largest and mean per-basket type I error, FWER, power, all correct, expected N,
  duration;
* figures, `timing.csv`, `results.rds` and `session.txt` (git SHA, seed, `sessionInfo()`).

**Run-time estimate for `N_SIM=1000`.** 2 recruitment designs x 10 scenarios x 4 models.
* Engine cost: about 0.035 s per data state for two-way EXNEX and 0.007 s for one-way EXNEX (measured on 2 cores,
  reference BLAS); no borrowing is negligible.
* Upper bound, assuming all data states distinct: about 3 CPU-hours. Two-way takes about 2 h of it, mostly the
  unbalanced design with its 16 analyses.
* Wall time: about 25–35 min on 8 cores, or 1.5 h on 2 cores.

The legacy study is run as before: `Rscript legacy/Setting.R` (see `legacy/README.md`).

## Status (2026-11-26)

Tests run on this branch (R 4.3.3) and their results are listed in the pull request and in `PROJECT_HISTORY.md`.
No full simulation has been run. Open points: [`docs/known-limitations.md`](docs/known-limitations.md); in short:
* partial exchangeability (row-only or column-only borrowing);
* hyperprior sensitivity runs;
* entry-time events ([`docs/entry-time-redesign.md`](docs/entry-time-redesign.md));
* outcome delays and drift scenarios;
* re-adding the design-optimal comparators and flexible inclusion on the new code.
