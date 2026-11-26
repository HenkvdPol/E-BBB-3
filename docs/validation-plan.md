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
Rscript validation/open_ended_feasibility.R  # exact open-ended single-basket results (~1 min)
```

`tests/run_tests.R` reports PASS, FAIL or SKIP per test and stops with exit status 1 on any failure. Details are
written to `tests/output/test_results.csv`. Main tests:

| File | What it tests |
|---|---|
| `test_eprocess.R` | grid alignment; restriction; block factor against a direct sum; E_p0[factor] = 1 and E_p[factor] <= 1 below p0; invalid priors are rejected; scaled-Beta benchmark |
| `test_oneway.R` | regression: the port reproduces legacy `posterior_exnex()` to < 1e-15; no borrowing = EX weight 0; permutation symmetry; determinism |
| `test_twoway.R` | tables against a direct sum; deterministic empty-platform quadrature against Monte Carlo; missing cell = empty cell; reproducibility; transposition; borrowing structure; agreement with the quadrature reference |
| `test_validity.R` | exact type I error of a null basket (all 6^6 block paths enumerated) <= alpha for every engine, including scaled Beta, with adversarial other baskets |
| `test_open_ended.R` | futility layer (E_p1[factor] = 1, <= 1 above p1); tables growing beyond n_max equal full tables; telescoping in the simulator; exact DP = path enumeration of the recursion; exact alpha and beta bounds over 1000 patients; open-ended simulator statuses |
| `test_simulate.R` | decisions; reproducibility; predictability (outcomes that differ only inside a block leave that block's prior unchanged); summaries |

* **Pull requests:** the GitHub workflow runs parsing, `tests/run_tests.R` and `legacy/checks/fw_engine.R`.
* **Manual dispatch:** the workflow runs the slow tests and the legacy exact, validity, regression and unit checks.

A green CI run is evidence that those particular scripts passed. It is not a proof for all possible data-generating
processes. Engine-level exact enumeration also does not establish that the simulator passes the correct weights
to the factor calculation; that requires separate integration tests.

## C. Required simulation scenarios

The broader research plan includes the following scenarios. The current framework-only study uses the subset
mapped below; drift, delayed outcomes and adaptive allocation remain later work.

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
calendar time and borrowing weights. Extension and acceleration frequencies belong to the later flexible-design study.
Distinguish per-basket type-I error from platform FWER in every table and figure.

**Status (2026-11-26).** `run_framework_study.R` covers the static heterogeneity themes in items 1–8 and 12
for the 3 x 3 grid. Some scenarios address more than one theme; a single active cell and a diagonal interaction
are distinct departures from an additive model.

| Plan item | Study scenario or design |
|---|---|
| 1 | global null |
| 2 | one active basket |
| 3 | active tumour type (row 1); active biomarker (column 1); row and column active |
| 4 | one active basket |
| 5, 6 | active tumour type (row 1); active biomarker (column 1) |
| 7 | interaction (diagonal); mixed active and null |
| 8 | global alternative; strong heterogeneity |
| 10 (staggered entry) and 12 | the unbalanced recruitment design: missing cell, rates 0.5–1, a basket opening at t = 20 |

The additive row + column scenario additionally checks graded effects. Not yet defined: item 9 (drift), delayed
outcomes (item 10) and item 11 (data-dependent allocation). Extension and acceleration are legacy-only.
The implementation report records only smoke runs (`N_SIM = 2`); these do not establish operating characteristics
and need repeating after simulator corrections.

## D. Approximation assessment

Compare:
- The restricted posterior represented by normalized theta-grid weights.
- The current moment-matched scaled-Beta prior.
- An empirical posterior-draw prior, if implemented with draws generated before
  the next block and then frozen.

For each, compare block factors, predictive distributions, power, expected sample
size and borrowing summaries. Evaluate rare-response and extreme-response states,
where the restricted posterior can be skewed or concentrated near `p0`.

**Status (2026-11-26).**
* The new code restricts and normalizes the computed grid posterior, then evaluates the block-factor finite sum
  directly. The two-way weights retain MCMC, interpolation, truncation and grid approximation; they are not an
  exact continuous posterior.
* The framework study runs the scaled-Beta representation as a separate design, so power and sample size can be
  compared directly.
* For the two-way engine, the production sampler is compared with deterministic quadrature, JAGS and the one-way
  reduction in `validation/twoway_validation.R`; results and reference limitations are in `docs/two-way-validation.md`.
  The quadrature shares production tables; JAGS uses continuous parameter versions; the one-way reduction is exact
  before conditional theta truncation. None alone establishes convergence of every production approximation.
* Before a production power study, compare a few difficult early/empty/heterogeneous states under finer theta/m/
  Rao-Blackwell grids and longer chains. Check restricted means, tail-sensitive block factors and convergence
  diagnostics. Finite-grid summation is direct; accuracy of its estimated weights remains a separate question.
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
