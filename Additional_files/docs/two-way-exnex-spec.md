# Two-way EXNEX model: implementation specification

## Purpose

Replace the current one-way cohort exchangeability model with a model that
represents the two-way basket structure in `framework.tex`:

- row `j`: tumour type;
- column `l`: biomarker group;
- basket `c = (j,l)`.

The current one-way model must remain available as a benchmark and fallback.

## Proposed EX component

For basket `c=(j,l)`, let

```math
theta_jl = mu + a_j + b_l + epsilon_jl
```

with

```math
a_j ~ N(0, sigma_a^2)
b_l ~ N(0, sigma_b^2)
epsilon_jl ~ N(0, sigma_epsilon^2)
```

and proper hyperpriors on `mu`, `sigma_a`, `sigma_b`, and `sigma_epsilon`.
The existing framework uses half-normal priors for the standard deviations;
the final scales and centering constraints must be justified and prespecified.

Introduce basket-specific exchangeability indicators:

```math
delta_jl ~ Bernoulli(w_jl)
```

and

```math
theta_jl | delta_jl = 1 ~ N(mu + a_j + b_l, sigma_epsilon^2)
theta_jl | delta_jl = 0 ~ N(m_jl, v_jl).
```

The NEX component protects against basket-specific departures, including
treatment effects that differ sharply at a particular tumour/biomarker
intersection. Consider whether partial exchangeability (row-only or
column-only borrowing) is needed; do not silently assume a single EX/NEX
indicator is sufficient for every scientific setting.

## Required design decisions

Before coding, specify:

1. Whether `delta_jl` is independent across baskets or has a structured prior.
2. Whether row and column effects are centered/sum-to-zero for identifiability.
3. The hyperprior scales and sensitivity ranges.
4. Whether a basket can borrow through its row, its column, both, or neither.
5. How empty rows, empty columns and missing cells are handled.
6. Whether external data `D0` enter the model and how they are discounted.
7. Whether the posterior target is conditional on all platform data, including
   the target basket's own previous data.
8. How posterior uncertainty is computed and validated.

## Computational strategy

Start with a small 3 x 3 grid and compare at least two approaches:
- deterministic quadrature / grid integration on a deliberately small model;
- a validated Bayesian computation approach such as Stan or another MCMC/INLA
  implementation, subject to model support and diagnostics.

Do not commit to INLA unless the chosen mixture and latent structure are
supported by the intended likelihood and prior. Compare posterior marginals and
restricted posterior moments against a trusted reference on small examples.

Keep the interface separate from the e-process layer. The borrowing module should
return a proper posterior or a normalized restricted prior for each basket. The
e-process module should only require a normalized, predictable prior on `(p0,1]`.
This separation is essential: changing the borrowing model must not change the
validity proof.

## Validation scenarios

At minimum, simulate:
- no row, column or interaction heterogeneity;
- row-only heterogeneity;
- column-only heterogeneity;
- one basket-specific interaction;
- one active tumour type across biomarker groups;
- one active biomarker across tumour types;
- one active intersection only;
- active and null baskets mixed in the same row and column;
- extreme non-exchangeable baskets;
- sparse/unbalanced grids and late-opening baskets.

Evaluate posterior calibration, false borrowing, borrowing weights, power,
per-basket type-I error, FWER, expected sample size and sensitivity to prior
scales. Verify the e-process validity independently of the fitted model.

## Acceptance criteria

- The old one-way engine still reproduces its existing reference checks.
- The new model passes simulation-based parameter recovery and posterior
  predictive checks on small known examples.
- Restricted priors are normalized and supported on `(p0,1]`.
- The restricted prior is computed using only information available before the
  next block and is frozen over that block.
- Exact or trusted-reference comparisons are documented.
- The e-process validity layer is unchanged and independently tested.
- No full production simulation is run until the researcher approves it.

---

## Implemented model and decisions (2026-11-26, branch `feature/lean-framework-twoway-exnex`)

Code: `additional_functions.R` (engine), `Additional_files/validation/reference_functions.R` (references). Validation: `Additional_files/docs/two-way-validation.md`.
Decisions D1–D4 were confirmed by the researcher; the others are documented defaults
(`Additional_files/docs/implementation-plan.md`, section 7).

### Answers to the required design decisions

1. **delta_jl.** Independent across baskets, `w_jl = w = 0.5`, the value of the one-way framework design (D3:
   binary EX/NEX first).
2. **Centering.** No hard sum-to-zero constraint (D5). Row and column effects are proper random effects around `mu`,
   so `theta_jl` is identified and is the only quantity the e-process needs. The sampler moves exactly along the
   three non-identified directions (mu+k, a-k), (mu+k, b-k) and (a+k, b-k). These moves help exploration of the
   posterior ridge; convergence and Monte Carlo precision still need checking.
3. **Hyperpriors** (D2: match the one-way benchmark's total expected EX variance before discretization).

   | Parameter | Prior | Reason |
   |---|---|---|
   | `mu` | N(logit p0, 10²) on [-8, 4] | same continuous specification as the one-way benchmark, which integrates mu on a grid |
   | `sigma_a`, `sigma_b`, `sigma_e` | HN(1/sqrt 3) each, on the grid 0.05, 0.15, …, 2.45 | Before truncation/gridding, E(sigma_a² + sigma_b² + sigma_e²) = 1 = E(tau²) for tau ~ HN(1); each continuous two-way SD has about 2e-5 mass above 2.45 |
   | NEX | N(logit pbar, 1/pbar + 1/(1 - pbar)), pbar = (p0 + p1)/2 | as in the benchmark |

   This matches a prior expectation, not the full distribution: the untruncated variance sums are
   chi-squared(3)/3 and chi-squared(1), respectively. The finite grids also perturb the expectations slightly.
   Sensitivity ranges to run: each scale × 0.5 and × 2, and `w` in {0.25, 0.75}. These are not yet run.
4. **Borrowing route.** Through both row and column, or not at all (NEX). Row-only and column-only components
   (partial exchangeability) are a planned sensitivity model. They need a decision on the four component weights
   and are not implemented.
5. **Empty rows, columns and cells.**
   * A cell that is not in the trial is left out of the layout.
   * An open basket without patients stays in the model and gets its prior predictive through its row and column.
   * Both give identical posteriors for the other baskets; this is tested.
   * Empty rows and columns contribute their prior only.
6. **External data D0.** None (D8).
7. **Posterior target.** The posterior given all platform data, including the basket's own data (framework.tex,
   eq. 8). It is restricted to (p0, 1] and frozen for the next block.
8. **Uncertainty computation and validation.** See below and `Additional_files/docs/two-way-validation.md`.

### Computation

* **Finite sums over theta and delta.** Given `m_c = mu + a_j + b_l` and `sigma_e`, the posterior of
  `theta_c` and `delta_c` depends only on the basket's own data. The 1-D tables of `D_c(m, sigma_e)` sum these
  variables out directly on the theta grid (cell width 0.05, aligned at logit p0). The sampler linearly
  interpolates `log D` on an m lattice of step 0.05 and clamps m to the lattice limits.
* **Production sampler** for `(mu, a, b, sigma_a, sigma_b, sigma_e)`: Metropolis-within-Gibbs, vectorised over
  data states and chains.
  * Moves: random walk for `mu`, `a` and `b`; exact Gibbs draws for `sigma_a` and `sigma_b` on their grids; a random
    walk on the `sigma_e` grid; exact moves along the non-identified directions.
  * Steps are adapted during burn-in only.
  * Default: 2 chains × (500 burn-in + 1000 kept).
* **Posterior of theta_c on the grid.** Rao-Blackwellised over the draws of `(m_c, sigma_e)`, deposited by
  cloud-in-cell on a lattice of step 0.1. Finite MCMC, interpolation, clamping and the finite grids all
  approximate the posterior weights. The e-process evaluates the finite sum using the supplied restricted
  weights directly; this is not an exact continuous posterior calculation.
* **Platform without patients.** The prior predictive uses deterministic quadrature, with mu step 0.02 and
  m-bin probabilities placed on the Rao-Blackwell lattice. It has no sampling error but retains quadrature error.
* **Validity under sampling.** The sampler's randomness is independent of future outcomes, and its output is
  fixed at the analysis. Provided its actual weights are normalized and supported on `(p0, 1]`, the prior is
  predictable and frozen, and validity does not depend on Monte Carlo accuracy.
* **References.**
  * Deterministic quadrature over the column levels `c_l = mu + b_l` (L ≤ 3), with `mu` integrated in closed form
    and each row effect in 1-D. It shares the production likelihood tables, so it checks hyperparameter
    integration without independently establishing convergence of those tables' grids.
  * JAGS, with continuous theta and sigmas: an independent implementation with its own Monte Carlo error
    and a different numerical representation.
  * For L × 1 layouts, the one-way reduction on its own grid engine. The Gaussian reduction is exact before
    conditional theta truncation; the implemented finite grids make it a numerical cross-check.
