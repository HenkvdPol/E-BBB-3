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
