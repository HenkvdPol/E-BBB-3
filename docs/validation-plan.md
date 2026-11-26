# Validation and reproducibility plan

## A. Mathematical invariants

For every basket and every block:
- The prior is measurable from the platform history available at block start.
- The prior is a normalized probability measure supported on `(p0, 1]`.
- It is held fixed for all outcomes in the block.
- The block likelihood ratio uses only outcomes in that new block.
- Basket allocation is known before the corresponding outcome.
- The null is conditional on the pre-outcome platform history.

These conditions, not correctness of EXNEX, are the basis of the type-I error proof.

## B. Automated checks

Run from repository root:

```sh
Rscript scripts/parse_all.R            # every R file parses
Rscript tests/run_tests.R              # framework code in R/: unit, regression and exact validity tests (~1.5 min)
Rscript tests/run_tests.R --slow       # plus the slow reference comparisons (3 x 1 vs one-way equivalent, JAGS)
Rscript legacy/checks/fw_engine.R      # legacy framework engine
Rscript legacy/checks/fw_validity.R    # legacy exact validity, and further legacy checks in legacy/checks/
Rscript validation/twoway_validation.R # two-way engine vs references (~30 min; writes validation/results/)
```

`tests/run_tests.R` reports PASS, FAIL or SKIP per test and stops with exit status 1 on any failure. Details are
written to `tests/output/test_results.csv`. Main tests:

| File | What it tests |
|---|---|
| `test_eprocess.R` | grid alignment; restriction; block factor against a direct sum; E_p0[factor] = 1 and E_p[factor] <= 1 below p0; invalid priors are rejected; scaled-Beta benchmark |
| `test_oneway.R` | regression: the port reproduces legacy `posterior_exnex()` to < 1e-15; no borrowing = EX weight 0; permutation symmetry; determinism |
| `test_twoway.R` | tables against a direct sum; exact empty-platform prior against Monte Carlo; missing cell = empty cell; reproducibility; transposition; borrowing structure; agreement with the quadrature reference |
| `test_validity.R` | exact type I error of a null basket (all 6^6 block paths enumerated) <= alpha for every engine, including scaled Beta, with adversarial other baskets |
| `test_simulate.R` | decisions; reproducibility; predictability (outcomes that differ only inside a block leave that block's prior unchanged); summaries |

* **Pull requests:** the GitHub workflow runs parsing, `tests/run_tests.R` and `legacy/checks/fw_engine.R`.
* **Manual dispatch:** the workflow runs the slow tests and the legacy exact, validity, regression and unit checks.

A green CI run is evidence that those particular scripts passed. It is not a proof for all possible data-generating
processes.

## C. Required simulation scenarios

For the one-way prototype and future two-way model, include:

1. Global null.
2. One active basket with all other baskets null.
3. Active row, null column and vice versa.
4. A single active basket with strong interaction.
5. Row-only heterogeneity.
6. Column-only heterogeneity.
7. Basket-specific interaction heterogeneity.
8. All active and mixed response rates.
9. Temporal drift and non-constant response probabilities.
10. Staggered entry and delayed outcome observation.
11. Data-dependent allocation, with the allocation rule explicit.
12. Sparse and unbalanced two-way grids.

Report per-basket type-I error, platform FWER, power, expected sample size,
calendar time, extension frequency, acceleration frequency and borrowing weights.
Distinguish per-basket type-I error from platform FWER in every table and figure.

**Status (2026-11-26).** `run_framework_study.R` defines items 1–8 and 12 for the 3 x 3 grid:

| Plan item | Study scenario or design |
|---|---|
| 1 | global null |
| 2 | one active basket |
| 3 | active tumour type (row 1); active biomarker (column 1); row and column active |
| 4, 7 | interaction (diagonal); mixed active and null |
| 5, 6 | additive row + column |
| 8 | global alternative; strong heterogeneity |
| 10 (staggered entry) and 12 | the unbalanced recruitment design: missing cell, rates 0.5–1, a basket opening at t = 20 |

Not yet defined: item 9 (drift), delayed outcomes (item 10) and item 11 (data-dependent allocation). Extension and
acceleration are legacy-only. Only smoke runs (`N_SIM = 2`) have been executed.

## D. Approximation assessment

Compare:
- The exact restricted posterior on the deterministic theta grid.
- The current moment-matched scaled-Beta prior.
- An empirical posterior-draw prior, if implemented with draws generated before
  the next block and then frozen.

For each, compare block factors, predictive distributions, power, expected sample
size and borrowing summaries. Evaluate rare-response and extreme-response states,
where the restricted posterior can be skewed or concentrated near `p0`.

**Status (2026-11-26).**
* The new code uses the exact restricted grid posterior.
* The framework study runs the scaled-Beta representation as a separate design, so power and sample size can be
  compared directly.
* For the two-way engine, the production sampler is compared with deterministic quadrature, JAGS and the equivalent
  one-way model in `validation/twoway_validation.R`; results are in `docs/two-way-validation.md`.
* A frozen posterior-draw prior is not implemented.

## E. Reproducibility record

For every full simulation, save:
- Git commit SHA.
- R version and `sessionInfo()`.
- OS and CPU/parallel settings.
- Random seed.
- `N_SIM` and calibration simulation count.
- Complete settings object.
- Whether design banks were rebuilt or reused.
- Output checksums and timestamp.

`run_framework_study.R` writes the git SHA, `N_SIM`, cores, seed and `sessionInfo()` to
`Result/framework_study/session.txt`, and the full settings to `results.rds`.

Do not commit large `Result/`, `design_banks/` or `checks/output/` directories.
Do not launch the full simulation without researcher approval.
