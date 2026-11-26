# Implementation plan: lean framework study and two-way EXNEX

Branch: `feature/lean-framework-twoway-exnex`.

**Status:** D1–D4 were confirmed by the user on 2026-11-26:
* D1: new `R/` code, with the current code in `legacy/`;
* D2: hyperpriors matched to the one-way benchmark;
* D3: binary EX/NEX first;
* D4: the exact restricted grid posterior.

D5–D9 are the proposals below; they are implemented as defaults and remain open to revision.
`framework.tex` (repository version, with the proof of Proposition 1) is the lead document.

## 0. Order of work

1. **Simplify.** A short, readable code base for the framework method only.
2. **Framework-only simulation study.** One e-process, three borrowing engines (none, one-way EXNEX, two-way EXNEX) on
   a J × L basket grid. Smoke tests only; pilot and full runs only after approval.
3. **Two-way EXNEX module.** A 3 × 3 grid first, validated against trusted reference calculations, plugged into the
   same restricted-prior interface.
4. **Operating characteristics.** Define all scenarios of the validation plan; run nothing larger than a smoke test
   without approval.
5. **Separate:** entry-time redesign, as a design proposal only. The guard stays.
6. **Later ('extras'):** design-optimal comparators, flexible inclusion, Bayesian comparators and Simon. These stay
   available in `legacy/` and are re-integrated afterwards.

## 1. Code structure (proposal)

```
R/eprocess.R        validity layer: prior objects on (p0, 1], checks, block factor, update, rejection
R/borrow_none.R     no borrowing: NEX posterior from the basket's own data, restricted to (p0, 1]
R/borrow_oneway.R   one-way EXNEX (ported from additional_functions.R; regression-tested against it)
R/borrow_twoway.R   two-way EXNEX (new)
R/simulate.R        calendar-time trial loop for a J x L grid; priors frozen per block
R/summarise.R       per-basket type I error, FWER, power, ESS, time, borrowing, calibration
R/plots.R
run_framework_study.R   the only file to run for the framework study (settings + run)
tests/run_tests.R       base-R test runner (no package needed in CI); tests/test_*.R
legacy/             current Setting.R, additional_functions.R, design_optimal_evalues.R and checks/, unchanged
```

**Interface between borrowing and the e-process.**
* A borrowing engine is a function `prior_fn(y, n, layout, settings)`. It gets the counts per basket at an analysis
  and returns one prior per basket: `list(p = support points, w = weights)`, with all points in (p0, 1] and weights
  ≥ 0 summing to 1.
* The e-process sees nothing else. Changing the engine cannot change the validity layer.

**Validity layer** (`R/eprocess.R`, mathematics unchanged from framework.tex):
* block factor e = Σ_g w_g (p_g/p0)^y ((1 − p_g)/(1 − p0))^(n − y);
* E ← E · e;
* reject at the first analysis with E ≥ 1/α;
* explicit assertions on support, normalisation and finiteness.

Predictability and freezing are guaranteed by the structure of the simulator: the prior is computed from the counts at
the previous analysis and stored, and the block factor uses only the new block's counts.

## 2. Framework-only simulation study

* **Defaults (to confirm):**
  * p0 = 0.10, p1 = 0.30, α = 0.05 per basket;
  * a 3 × 3 grid (rows = tumour types, columns = biomarker groups);
  * n_max = 30 per basket;
  * an analysis every 5 time units (calendar time);
  * accrual per basket balanced (1 per time unit) or unbalanced (rates 0.25–1);
  * staggered entry and empty baskets as options.
* **Engines compared:** none, one-way EXNEX (benchmark; cohorts treated as an exchangeable list), two-way EXNEX.
* **Stopping:** stop a basket when 1/α can no longer be reached. This is validity-neutral and only shortens trials.
* **Measures:**
  * per-basket type I error, each null basket separately, never pooled with FWER;
  * FWER, power per active basket, expected n, time to decision;
  * borrowing: prior mean, posterior EX probability and false borrowing (prior mass far above the basket's own
    estimate);
  * posterior calibration: coverage of 90% credible intervals for p_c at the final analysis.

## 3. Two-way EXNEX

**Model** (framework.tex eq. 9). For basket c = (j, l):
* δ_c ~ Bernoulli(w);
* EX: θ_c = μ + a_j + b_l + ε_c;
* NEX: θ_c ~ N(m, v);
* a_j ~ N(0, σ_a²), b_l ~ N(0, σ_b²), ε_c ~ N(0, σ_ε²).

**Proposed defaults, matching the total expected EX variance of the one-way benchmark before numerical
truncation and discretization (not its complete prior distribution)** (decision D2):

| Component | Proposal | Justification |
|---|---|---|
| μ | N(logit p0, 10²) truncated to [−8, 4] | identical to the one-way benchmark |
| σ_a, σ_b, σ_ε | half-normal with scale 1/√3 each | E[σ_a² + σ_b² + σ_ε²] = 1 = E[τ²] under the benchmark's τ ~ HN(1) |
| w | 0.5, independent across baskets | the value of the one-way framework design |
| NEX | N(logit p̄, 1/p̄ + 1/(1 − p̄)), p̄ = (p0 + p1)/2 | unit-information prior, as in the one-way benchmark |
| Sensitivity | each scale × 0.5 and × 2; w ∈ {0.25, 0.75} | – |

**Other choices:**
* **Identifiability:** proper random effects without hard sum-to-zero constraints. θ_c is what the prior for the
  e-process needs, and it is identified. Row and column effects are reported centred (μ + ā + b̄, a_j − ā, b_l − b̄).
  The alternative is hard sum-to-zero constraints, which change the implied prior on θ_c.
* **Missing cells and empty baskets:** cells not in the trial are left out of the model. An open basket without
  patients stays in: it gets the prior predictive through its row and column.
* **Other settings:** no external data D0; own data included in the posterior (framework.tex); δ independent across
  baskets.

**Partial exchangeability** (decision D3). An optional four-component version: both row and column, row only, column
only, or NEX. The engine handles any number of components per basket; their weights need a decision.

**Computation.**

Reference (trusted), for validation only:
* **R1 — deterministic quadrature on reduced models.** Two set-ups: a 2 × 2 grid with all parameters; and a 3 × 3 grid
  with hyperparameters fixed and hard sum-to-zero constraints. In both, δ and ε are integrated per basket in 1-D on the
  θ grid and the rest on a tensor grid.
* **R2 — JAGS (installed).** The full 3 × 3 model with all hyperpriors; long chains with R-hat and effective-sample-size
  diagnostics.

Production candidates for simulation, which needs about 20 ms or less per 3 × 3 fit:
* **P1 — deterministic nested integration.** Hyperparameters (σ_a, σ_b, σ_ε) on a small grid. For each: a Gaussian
  (Laplace) approximation of (μ, a, b), with the EX/NEX mixture integrated exactly per basket from 1-D tables. Then the
  posterior of θ_c on a grid. Reproducible, no random draws.
* **P2 — MCMC vectorised over all simulated trials** (Metropolis-within-Gibbs on (μ, a, b, log σ)), with a
  Rao-Blackwellised grid posterior for θ_c. Draws are made at the analysis and then frozen.

The production method is chosen by agreement with R1 and R2 (restricted mean and second moment, total variation of the
grid prior, block factors) and by run time; the comparison is documented.

**Tests:**
* normalisation, support, finite moments;
* missing cells, row/column indexing;
* deterministic reproducibility;
* the e-process is invariant to the engine: exact enumeration of all sequences of one basket at small n with priors
  from each engine gives type I error ≤ α;
* the one-way port equals the legacy engine.

## 4. Operating characteristics (definitions now, runs after approval)

* **Scenarios:**
  * global null; row-only and column-only heterogeneity; a basket-specific interaction;
  * one active basket among nulls; an active tumour type across biomarkers; an active biomarker across tumour types;
  * mixed active and null baskets in the same row and column; strongly heterogeneous effects;
  * sparse and unbalanced grids; empty and late-opening baskets.
* **Misspecification:** data generated outside the model, for example a strong interaction.
* **Measures:** as in section 2. Per-basket type I error and FWER are always reported separately.

## 5. Entry-time redesign (separate)

A design proposal in `docs/entry-time-redesign.md`.
* **Event-driven simulator:** cohort entry as its own event. Before the first patient of a cohort, the borrowed
  alternative is computed from the history at entry. In framework terms this is an analysis of that basket with an
  empty block.
* **The guard on `redesign = "entry"` stays** until this is implemented and tested.

## 6. Commits and pull request

Small commits, in this order:
1. plan;
2. move to `legacy/` with path fixes, plus a check that nothing changed;
3. validity layer and tests;
4. one-way port and regression test;
5. simulator and summaries, with a smoke test;
6. two-way reference implementations (R1, R2);
7. production method and validation report;
8. scenario definitions;
9. docs (README, known limitations, validation plan) and dated entries appended to `PROJECT_HISTORY.md`.

The pull request lists the tests actually run (with configuration and run time), the tests not run, unresolved choices
and computational limits.

## 7. Decisions needing confirmation

| | Decision | Proposal |
|---|---|---|
| D1 | Code structure | New lean `R/` code for the framework study; current code moved unchanged to `legacy/` |
| D2 | Two-way hyperpriors | Match the one-way benchmark (table in section 3), with sensitivity analyses |
| D3 | Partial exchangeability | Binary EX/NEX first (framework.tex); row-only and column-only components as a sensitivity model once weights are agreed |
| D4 | Prior handed to the e-process | The exact restricted posterior on a grid of p values (framework.tex eq. 8 up to quadrature); the scaled Beta stays as benchmark |
| D5 | Identifiability | Random effects without hard constraints; centred summaries |
| D6 | Study defaults | Section 2 |
| D7 | Production computation | **Implemented: P2** (vectorised MCMC, Rao-Blackwellised grid prior). P1 (deterministic nested Laplace) was not pursued: the EX/NEX mixture makes the hyper-level posterior heavy-tailed and possibly multimodal, so a Laplace step would itself need validation. P2 is checked against R1 and R2 in `docs/two-way-validation.md` |
| D8 | External data D0 | None |
| D9 | Pilot runs (N_SIM ≤ 200) | Only after approval; until then smoke tests (N_SIM ≤ 2) and exact checks only |

## 8. As implemented (2026-11-26)

The independent review and corrections are recorded in `docs/review-2026-11-26.md`. The original first-prior
replication across trials was incorrect; it is corrected in `d0b8d5b`, with three simulator regressions and passing
CI. The three-page `setup.tex`/`setup.pdf` now describe this scoped study. Proposal text above is retained as the
record of decisions; rates actually used are 0.5 and 1, the dedicated false-borrowing statistic is not implemented,
and the posterior grid weights are numerical approximations. The finite e-factor sum uses these weights directly.


* **theta grid.** Cells of width 0.05 (legacy: 0.02), aligned at logit(p0), so the restriction to (p0, 1] is exact
  on cells.
  * The one-way port on the legacy grid reproduces legacy to < 1e-15.
  * On the new grid it agrees to 3e-5 in the restricted mean.
* **Two-way production engine.** Vectorised Metropolis-within-Gibbs with exact moves along the non-identified
  directions; theta and delta integrated exactly; Rao-Blackwellised grid prior. Defaults: 2 chains, 500 burn-in,
  1000 kept.
* **R1: deterministic quadrature.** Over the column levels c_l = mu + b_l (mu in closed form) and each row effect
  in 1-D; layouts with up to 3 columns.
  * Lesson learned: the grid for c must cover the whole prior range of mu. A window around the data alone cuts the
    configurations in which a whole column is NEX, which biased P(EX) upwards by up to 0.045.
* **R2: JAGS.** Continuous theta and sigmas.
* **R3: one-way model equivalent to an L x 1 layout.** Agrees with R1 to 4 digits.
* **Tests.** Base-R runner, `tests/run_tests.R`; testthat is not needed.
* **Partial exchangeability (D3).** Not implemented; it needs the weights of the row-only and column-only
  components.
