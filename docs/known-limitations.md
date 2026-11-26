# Known limitations and implementation decisions

This document separates verified properties of the current code from planned work.
It should be updated when an item is resolved.

## Statistical scope

1. **Per-basket error control only.** The framework targets one-sided type-I error
   `alpha` for each basket separately. It does not control platform-wide FWER or FDR.
2. **Conditional null.** Validity assumes
   `P(Y_i = 1 | G_i) <= p0` for every patient assigned to the tested basket.
   A marginal-average null is insufficient in general.
3. **Outcome process.** Allocation must be known before the outcome. The protocol
   must specify outcome ordering, ties, pending outcomes and potentially informative
   observation delays.
4. **Predictable prior.** The prior must be a normalized probability measure on
   `(p0, 1]`, selected using history available before the next block and frozen
   during that block.
5. **Two-way model absent.** The R prototype currently fits one-way EXNEX across
   cohorts. Tumour-type row effects, biomarker column effects and basket-specific
   interactions have not been implemented.
6. **Restricted posterior approximation.** The primary e-mixture matches the first
   two moments of the restricted posterior with a scaled Beta distribution. This
   is a valid alternative prior if properly normalized and frozen, but it is not
   the exact posterior restriction in `framework.tex`.
7. **Design-optimality scope.** Baas et al.'s finite-horizon optimality results are
   for the stated IID Bernoulli point-alternative problem. A borrowed alternative
   that is updated over time is an adaptive extension; do not claim the original
   global optimality result for that extension without a separate proof.
8. **Monte Carlo vs deterministic approximation.** The framework text describes
   posterior draws as one possible numerical approximation. The current R code
   uses deterministic quadrature and moment matching, not retained Monte Carlo
   draws for the primary e-mixture.

## Software limitations

1. **Entry-time redesign disabled.** The simulator's current calendar-time loop
   does not represent cohort opening as an independent event. The old
   `redesign = "entry"` condition could miss the intended update when a cohort
   enters between global analyses and accrues patients before the next analysis.
   This option is therefore excluded from active set-ups and now fails explicitly
   if selected. It must not be silently re-enabled.
2. **ESS-min added to wrapper and comparisons.** The DP constructor existed but
   was previously not exposed by the main wrapper or simulation settings. It is
   now exposed as `e-ESSmin` and included in exact-check enumeration. The checks
   still need to be run in an R environment after this change.
3. **No completed final full simulation.** The project history records that the
   final full run is outstanding. Historical simulation results are not a
   substitute for rerunning the final commit.
4. **No two-way engine.** The two-way EXNEX model remains a research development
   task; do not add a partially specified model to the production simulation.
5. **Reproducibility.** Package versions are not yet pinned with `renv.lock`.
   Record `sessionInfo()` with final simulation results and consider adding `renv`
   after the required packages and supported R version are agreed.

## Immediate safeguards

- Do not remove the `redesign = "entry"` guard until event-time semantics are
  implemented and tested.
- Do not interpret per-basket `alpha` as platform-wide alpha.
- Do not interpret moment-matching accuracy as proof of validity; validity follows
  from predictability, normalization and support.
- Do not launch full simulations without explicit approval.
