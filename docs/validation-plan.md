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
Rscript checks/fw_engine.R
Rscript checks/fw_validity.R
Rscript checks/exact_checks.R
Rscript checks/regression.R
Rscript checks/unit_flex.R
```

The GitHub pull-request workflow runs syntax parsing and `fw_engine.R`. The manual
workflow runs the more expensive exact, validity, regression and unit checks.
A green CI run is evidence that those particular scripts passed; it is not a proof
for all possible data-generating processes.

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

## D. Approximation assessment

Compare:
- The exact restricted posterior on the deterministic theta grid.
- The current moment-matched scaled-Beta prior.
- An empirical posterior-draw prior, if implemented with draws generated before
  the next block and then frozen.

For each, compare block factors, predictive distributions, power, expected sample
size and borrowing summaries. Evaluate rare-response and extreme-response states,
where the restricted posterior can be skewed or concentrated near `p0`.

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

Do not commit large `Result/`, `design_banks/` or `checks/output/` directories.
Do not launch the full simulation without researcher approval.
