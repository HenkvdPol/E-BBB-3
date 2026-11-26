# Alignment with eBBB (T. Jacobs)

Date: 2026-11-26. Reference: the R and Stan code of eBBB (`tijn-jacobs/eBBB`, folder `code/`), received as a zip file.
That code was only read and run locally for the comparison below; it is not part of this repository, and the eBBB
repository was not changed.

## Decisions (researcher, 2026-11-26)

* **Adopted from eBBB:**
  * the priors: no borrowing, power (equal and column weights), conjugate EX/NEX, similarity and the oracle;
  * the stopping rule of the fixed design: efficacy (`E >= 1/alpha`) or `n_max`, with no futility stop;
  * the two-way EX/NEX hyperparameters;
  * the grid of 4 tumour types x 3 biomarker groups;
  * eBBB's four scenarios.
* **Two-way model:** our validated grid sampler replaces eBBB's Stan fit.
* **Kept from E-BBB-3:** calendar-time recruitment with global analyses every 5 time units, `p0 = 0.10`,
  `p1 = 0.30`, `n_max = 30`, and our scenario family (rebuilt for 4 x 3).
* **Kept as an option:** the open-ended design with e-value futility (`simulate.R open`).

## What eBBB does and where it is here

| eBBB (`LR_evalue_multiple/lr_evalue.R`, `simulate_platform.R`, `exnex.stan`) | E-BBB-3 | |
|---|---|---|
| `borrowing_prior`, weights `w_row = w_col = 0` | `engine_ebbb(..., "no borrowing")` | same prior |
| `borrowing_prior`, `w_equal = (0.5, 0.5)`, `w_col = (0, 1)` | `engine_ebbb(..., "power", w_row, w_col)` | same prior |
| `exnex_prior`, `p_ex = 0.5`, EX borrows from the same column | `engine_ebbb(..., "EX/NEX")` | same prior; `pex` = posterior EX weight |
| `similarity_prior`, `epsilon = 2`, `tau = 0.3`, `jsd_beta` (200-point midpoint rule, base 2) | `engine_ebbb(..., "similarity")`, `jsd_beta()` | same prior |
| oracle: point mass at `p_true` if `p_true > p0`, else `p1` | `engine_ebbb(..., "oracle")` | same prior; exact points on the grid (`add_points`) |
| `exnex.stan`: `mu ~ N(logit p0, 2^2)`, `sigma_a, sigma_b ~ HN(1)`, `sigma_eps ~ HN(0.5)`, NEX `N(logit p0, 2^2)`, `p_ex = 0.5` | `twoway_prior()` via `make_engine("two-way EX/NEX")` | same model and hyperparameters; grid sampler instead of Stan; `mu` truncated to `[-10, 6]` (mass outside < 1e-4) |
| `a0 = b0 = 1` | `a0`, `b0` in `setting.R` | same |
| J = 4, L = 3 | `J`, `L` in `setting.R` | same |
| scenarios: global null, biomarker 1 active, interaction (biomarker 1 except tumour type 4), all active | first four entries of `Scenario` | same patterns, at `p0 = 0.10` / `p1 = 0.30` |
| stop at efficacy or `n_max` | `make_design(..., "fixed")`: `futility = "none"` | same rule |
| prior restricted to `(p0, 1]`, frozen per block, one e-value per basket, `alpha = 0.05` | Part 1 and `simulate_trials()` | same method |

The oracle needs a bet for the futility e-value of the open-ended design. It bets on the true rate if that rate is
below `p1`, otherwise on `p0`. This mirrors eBBB's efficacy choice; eBBB itself has no futility.

## Numerical agreement with eBBB's closed forms

`Additional_files/validation/ebbb_compare.R` sources eBBB's `lr_evalue.R` from a path given on the command line. It
compares block factors for 40 random data states on the 4 x 3 grid, every basket, and blocks of 5 with 0, 2 or 5
responses. Results are in `results/ebbb_compare.csv`.

| Prior | Largest relative difference of the block factor |
|---|---|
| no borrowing | 1.5e-4 |
| power, equal | 2.1e-4 |
| power, column | 2.0e-4 |
| EX/NEX | 2.2e-4 (EX weights: 1e-16) |
| similarity | 2.3e-4 |
| oracle | 3e-15 (exact) |

* **Remaining difference.** It comes from putting the continuous Beta on the theta grid: the density at the cell
  midpoints times the cell widths.
* **Tail mass.** The Beta mass beyond the grid (p > 0.9975 or below the lowest cell) goes into the end cells. Without
  this, the flat `Beta(1, 1)` prior of an empty basket lost 0.25% of its mass near p = 1. That changed a 5/5 block
  factor by 1.2%.
* **Rejected alternative.** Exact cell probabilities (`pbeta` differences) were worse than the midpoint rule (up to
  2e-3), because the grid puts each cell's mass on its midpoint.

## Two-way EX/NEX with the eBBB hyperparameters

See the addendum in `two-way-validation.md` (`validation/twoway_ebbb_hyper.R`: production sampler vs JAGS and a long
run, 4 x 3 layout, six data states).

## Differences that remain by design

| | eBBB | E-BBB-3 |
|---|---|---|
| Recruitment | one stream of `N_total = 800` outcomes, basket shares by row and column | calendar time: rate and opening time per basket (`Recruitment` in `setting.R`) |
| Analyses | per basket, every 5 own outcomes and at `n_max` | global analyses every `look_every = 5` time units; block sizes differ per basket |
| Refits of the platform model | Stan every 5 platform outcomes | at every global analysis (all models) |
| Hypotheses | `p0 = 0.30`, `p1 = 0.50` | `p0 = 0.10`, `p1 = 0.30` |
| Maximum per basket | `n_max = 100` | `n_max = 30` (fixed design); none in the open-ended design |
| Prior representation | closed-form Beta, point mass or posterior draws | normalized weights on the theta grid (agreement above) |
| Scenarios | 4 | 11 (eBBB's 4 and 7 of ours) |

## Tests

* **`Additional_files/tests/test_ebbb.R`** (40 checks). The closed forms are written out in the test, not copied from
  eBBB:
  * Beta on the grid against the closed-form factor;
  * the Beta parameters of no borrowing and the power priors;
  * the EX/NEX weights and mixture, and the limits `p_ex = 0, 1`;
  * the JSD and similarity weights;
  * transposition of the grid;
  * the oracle's factors and summaries;
  * the exact type I error of every eBBB prior (6^6 paths);
  * every model of `setting.R` on the unbalanced 4 x 3 layout;
  * the fixed design stopping only at efficacy or `n_max`.
* **Two-way engine, old defaults.** It gives bit-identical posteriors to `main` before this change (checked locally
  with the same seed).
