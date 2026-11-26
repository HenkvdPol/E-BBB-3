# Known limitations and implementation decisions

This document separates verified properties of the code from planned work. Update it when an item is resolved.
Last update: 2026-11-26 (branch `feature/lean-framework-twoway-exnex`: new framework code in `R/`, earlier code in
`legacy/`).

## Statistical scope

1. **Per-basket error control only.** The framework targets the one-sided type I error `alpha` for each basket
   separately. It does not control the platform-wide FWER or FDR. The study reports per-basket error and FWER separately
   (`R/summarise.R`); it does not currently calculate FDR.
2. **Conditional null.** Validity assumes `P(Y_i = 1 | G_i) <= p0` for every patient assigned to the tested
   basket. A marginal-average null is insufficient in general.
3. **Outcome process.** Allocation must be known before the outcome. The protocol must specify outcome ordering,
   ties, pending outcomes and potentially informative observation delays. The new simulator has immediate outcomes
   and deterministic accrual only.
4. **Predictable prior.** The prior must be a normalized probability measure on `(p0, 1]`, selected using history
   available before the next block and frozen during that block.
   * `R/eprocess.R::check_priors()` checks support and normalization at every analysis, including the assembled
     initial priors after the batching correction described in `docs/review-2026-11-26.md`.
   * `R/simulate.R` guarantees predictability.
   * `tests/test_simulate.R` checks it directly.
5. **Two-way model: binary EX/NEX only.**
   * *Resolved:* the two-way model exists (`R/borrow_twoway.R`).
   * *Open:* each basket borrows through its row and column together or not at all. Row-only and column-only
     components (partial exchangeability) are not implemented; their weights need a decision.
6. **Restricted posterior representation.**
   * *Implemented in `R/`:* the new code supplies normalized, restricted grid weights and evaluates their finite
     block-factor sum directly. Posterior weights are approximate (quadrature; MCMC and numerical grids for two-way
     fitting), rather than an exact continuous restricted posterior.
   * The moment-matched scaled Beta of the earlier prototype remains as a benchmark design,
     "one-way EXNEX, scaled Beta".
   * Both are valid priors.
7. **Design-optimality scope.** Baas et al.'s finite-horizon optimality results are for the stated IID Bernoulli
   point-alternative problem. A borrowed alternative updated over time is an adaptive extension; do not claim the
   original optimality result for it without a separate proof. These designs are in `legacy/` only.
8. **Monte Carlo in the two-way prior.**
   * The two-way engine integrates `(mu, a, b, sigmas)` by MCMC. The prior is Rao-Blackwellised on the grid and
     fixed at the analysis, so validity is unaffected.
   * Its accuracy is documented in `docs/two-way-validation.md`. With production settings, P(EX) carries a Monte
     Carlo error of a few hundredths; the restricted moments and block factors are much more precise.
   * The no-borrowing and one-way engines are deterministic.
9. **Grids.**
   * theta is restricted to `[-10, 6]` (p from 4.5e-5 to 0.9975) in cells of width 0.05.
   * The standard deviations of the two-way model live on the grid 0.05, 0.15, …, 2.45.
   * The mean `m` of the EX component is clamped to `[-10, 6]` in the tables.
   * The theta grid is shared. SD grids differ: the one-way benchmark uses 0.05, 0.15, ..., 3.95; the two-way grid
     above and EX-mean clamping apply to the two-way approximation. JAGS uses continuous theta and SDs.
   * HN(1/sqrt(3)) for the three two-way SDs matches total expected EX variance to HN(1) before discretization,
     not the full prior distribution. The finite-grid expected variances are only approximately equal.
10. **Common p0.** The new code assumes one `p0` for all baskets: the theta grid is aligned at `logit(p0)`.
    Basket-specific `p0` needs one grid per distinct `p0`.

11. **Open-ended baskets (added 2026-11-26).** No maximum sample size is needed: Proposition 1 holds over an
    unbounded horizon (`docs/open-ended-design.md`).
    * Without a futility rule, a null basket runs forever with probability at least 1 - alpha.
    * The optional e-value futility (framework.tex, Proposition 2) bounds P(futility | p >= p1) by beta and stops
      every basket with a constant response rate eventually.
    * The exact results cover one basket without borrowing, or borrowing from baskets whose data are complete.
      Borrowing from still-recruiting baskets needs simulation; only a smoke run was done.
    * The theta grid (width 0.05) becomes coarse for n in the thousands; validity is unaffected.

## Software limitations

1. **Entry-time redesign disabled.** The legacy option `redesign = "entry"` still stops with an error. The new
   simulator supports delayed recruitment through `entry`, but posterior refits occur only at scheduled global
   analyses; there is no dedicated arbitrary-entry refit event. Design proposal: `docs/entry-time-redesign.md`. Do not re-enable it before the
   event-driven simulator and its tests exist.
2. **ESS-min.** It was added to the legacy wrapper on 2026-11-26. The legacy checks that cover it (`exact_checks.R`)
   have not been rerun since.
3. **No completed full simulation.** Neither the legacy study nor the new framework study (`run_framework_study.R`)
   has been run with `N_SIM = 1000`. Only smoke runs (`N_SIM = 2`) were executed. A full run needs the researcher's
   approval.
4. **Run time of the two-way engine.**
   * About 0.035 s per data state (2 cores, reference BLAS), against 0.007 s for one-way EXNEX.
   * Each call also has a fixed cost of about 1 s.
   * About 3 CPU-hours is a provisional throughput extrapolation assuming distinct states, not an upper bound:
     it omits fixed fit, engine setup and parallel overhead. The 25–35 min on 8 cores remains an unverified
     estimate to check in an approved timing pilot.
5. **Quadrature reference limited to at most 3 columns.** Larger layouts can be checked against JAGS only (or after
   transposing).
6. **Reproducibility.** Package versions are not pinned (no `renv.lock`). `run_framework_study.R` and
   `validation/twoway_validation.R` record the git SHA and `sessionInfo()`.
7. **RNG streams before production.** Outcome seeds are `base + chunk`, while fitting seeds add the analysis
   index. This can reuse seeds across worker chunks and analyses. Use distinct outcome/fitting streams and
   validate paired comparisons before relying on ordinary independent-replication Monte Carlo uncertainty.
8. **Posterior summaries.** Empirical 90% interval coverage at stopping is a diagnostic, not a frequentist
   guarantee. Platform `mean_abs_error` is mean absolute estimation error; per-basket `post_bias` is signed bias.
   The proposed dedicated false-borrowing measure is not implemented.
9. **Later extras not yet on the new code.** Design-optimal comparators, flexible inclusion (acceleration and
   extension), Bayesian comparators and Simon designs are only in `legacy/`.

## Immediate safeguards

- Do not remove the `redesign = "entry"` guard until event-time semantics are implemented and tested.
- Do not interpret per-basket `alpha` as platform-wide alpha.
- Do not interpret moment-matching or Monte Carlo accuracy as proof of validity. Validity follows from
  predictability, normalization and support; `tests/test_validity.R` enumerates selected finite null configurations using each engine, and the simulator
  has separate batching/predictability regressions. These tests do not replace the general validity proof.
- Do not launch full simulations without explicit approval.
