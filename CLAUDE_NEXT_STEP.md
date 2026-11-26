# Claude task: implement and validate the two-way EXNEX borrowing model

You are working in the GitHub repository `HenkvdPol/E-BBB-3`. Read these files
before changing code:

1. `framework.tex`
2. `PROJECT_HISTORY.md`
3. `README.md`
4. `docs/known-limitations.md`
5. `docs/two-way-exnex-spec.md`
6. `docs/validation-plan.md`
7. `Setting.R`, `additional_functions.R`, `design_optimal_evalues.R`
8. The relevant scripts under `checks/`

## Objective

Implement the next research step: a two-way EXNEX borrowing model for baskets
indexed by tumour type `j` and biomarker group `l`, while leaving the existing
one-way engine available as a benchmark. The primary e-process and its validity
conditions must remain unchanged.

## Non-negotiable statistical constraints

- The tested null remains the conditional null in `framework.tex`.
- For each basket, the alternative prior supplied to the e-process must be a
  normalized probability measure supported on `(p0, 1]`.
- The prior must be measurable from information available before the next block
  and frozen throughout that block.
- The borrowing model is allowed to be misspecified; it must not enter the
  type-I error proof.
- Do not claim platform-wide FWER control from per-basket alpha.
- Do not claim the finite-horizon optimality result of Baas et al. for an
  adaptively updated borrowed alternative without a separate proof.
- Keep the one-way model and all its reference checks intact.

## Model

Start from the model specified in `framework.tex` and
`docs/two-way-exnex-spec.md`:

`theta_jl = mu + a_j + b_l + epsilon_jl` in the EX component, with row,
column and basket-specific effects; basket-specific EX/NEX indicators; and a
proper NEX distribution. Make identifiability constraints and hyperprior scales
explicit. Consider whether row-only and column-only partial borrowing are needed.
Do not invent prior scales without documenting sensitivity analyses.

## Work sequence

1. First write a short implementation plan and list the model choices that need
   confirmation. Do not start a full simulation.
2. Implement a small 3 x 3 reference model with a transparent data structure
   mapping each basket to `(j,l)`.
3. Compare the new posterior with a trusted reference on small datasets. Include
   empty cells, unbalanced recruitment and basket-specific interaction cases.
4. Add a new module (for example `two_way_exnex.R`) without changing the
   primary e-process functions. Its public interface should return posterior
   summaries or normalized restricted priors for each basket.
5. Add unit tests for posterior normalization, support, finite moments,
   missing cells, row/column indexing and deterministic reproducibility.
6. Add tests showing that the e-process validity layer is invariant to which
   borrowing engine supplies the prior, provided predictability, normalization
   and support hold.
7. Add operating-characteristic scenarios covering row-only, column-only,
   interaction-only and mixed active/null baskets.
8. Document computational cost and diagnostics. Compare at least a small
   deterministic reference calculation with the selected production method.
9. Update README, known limitations and validation plan. Preserve the project
   history; append dated decisions rather than rewriting historical entries.

## Separate issue: entry-time redesign

Do not re-enable `redesign = "entry"` as part of this task. The current simulator
does not model cohort entry as a separate calendar-time event. Create a separate
design proposal for an event-driven simulator that updates the borrowed
alternative at the actual entry time before the first patient in that cohort is
observed. The current explicit guard must remain until this is implemented and
tested. Do not approximate entry-time updates by updating at the next global
analysis after patients have already accrued.

## Validation and delivery

- Run syntax parsing and the focused new unit tests.
- Run `checks/fw_engine.R` and the existing regression checks to ensure the
  one-way implementation is unchanged.
- Run exact/validity checks only when computationally feasible and report
  runtime and configuration.
- Do not launch `N_SIM=1000` or any full production simulation without the
  researcher's explicit approval.
- Work on a feature branch, make small commits, and open a pull request.
- In the PR description, distinguish completed implementation, tests actually
  run, tests not run, unresolved modelling choices and computational limitations.
- Never report a test as passing unless it was executed and passed.

The primary deliverable is a validated two-way borrowing module that plugs into
the existing restricted-prior interface without altering the e-process validity
layer.
