# E-BBB-3: Anytime-valid e-values for basket and platform trials

Research prototype for basket/platform clinical trials with binary response outcomes.
The primary method is specified in [`framework.tex`](framework.tex). The design-optimal
single-arm e-value comparators follow Baas, ter Schure & van Rosmalen (2026),
[arXiv:2605.28653](https://arxiv.org/abs/2605.28653).

## Method hierarchy

1. **Primary:** basket-specific e-processes. At each interim analysis, select a
   predictable probability measure supported on `(p0, 1]`, freeze it for the next
   block, mix the block likelihood ratio, and multiply block factors. Reject basket
   `c` when `E_c >= 1/alpha`.
2. **Design-optimal e-value comparators:** GROW, P-max, ESS-min and e-based
   betting policies. These are single-arm finite-horizon designs; the current
   borrowed-alternative variants are adaptive extensions and should not be described
   as retaining the paper's original global optimality result.
3. **External comparators:** Bayesian EXNEX posterior-probability rules, Simon
   two-stage designs and exact binomial tests.

The guarantee is **per basket**. It is not automatically family-wise error control
across all baskets. The conditional null in `framework.tex` is stronger than a
marginal-average response-rate null; outcome ordering, allocation and delays matter.

## Repository layout

The R scripts currently remain at the repository root because `Setting.R` and the
check scripts resolve paths relative to their own locations. Do not move them into
subdirectories without first refactoring and testing those paths.

- `framework.tex`: statistical proposal and validity argument.
- `setup.tex` / `setup.pdf`: concise research/software set-up.
- `Setting.R`: settings and simulation entry point.
- `additional_functions.R`: EXNEX quadrature, borrowing, e-processes and simulation.
- `design_optimal_evalues.R`: GROW, P-max, ESS-min and e-based dynamic programs.
- `checks/`: exact, validity, regression and unit checks.
- `PROJECT_HISTORY.md`: historical decisions, results and open items.
- `docs/`: limitations, validation plan and two-way-model specification.
- `CLAUDE_NEXT_STEP.md`: scoped instructions for the next Claude development task.

Generated `Result/`, `design_banks/` and `checks/output/` are reproducible outputs
and are excluded from version control.

## Current implementation versus research target

The current borrowing engine is a **one-way EXNEX model over cohorts**. The
proposed tumour-type × biomarker-group model with row, column and basket-specific
effects is not yet implemented. The current framework e-mixture computes the
restricted posterior's first two moments and fits a scaled Beta distribution on
`(p0, 1]`; it does not evaluate the exact restricted posterior mixture. This
approximation preserves the e-process validity conditions when normalized,
predictable and frozen, but its operating-characteristic impact needs assessment.

The `redesign = "entry"` point-alternative option is deliberately disabled in active
set-ups and guarded in the simulator. Correct entry-time updates need a separate
calendar-time event model; updating at the next global analysis after a cohort has
already accrued patients is not equivalent to updating at cohort entry.

## Run checks

Use an R installation (R >= 4.3 recommended). Base R is sufficient for the core
engine and design library. `ggplot2` is needed for plots; `matrixStats` is optional
and used for faster row maxima.

From the repository root:

```sh
Rscript checks/fw_engine.R
Rscript checks/fw_validity.R
Rscript checks/exact_checks.R
Rscript checks/regression.R
Rscript checks/unit_flex.R
```

The pull-request workflow parses all R files and runs `fw_engine.R`. The full
validation suite is available through manual GitHub Actions dispatch because exact
enumeration and bank construction can be computationally intensive.

A smoke run is:

```sh
N_SIM=1 Rscript Setting.R
```

A full run (currently configured for 1,000 trials per scenario) should only be
started after the researcher approves the compute cost:

```sh
N_SIM=1000 Rscript Setting.R
```

Record the R version, package versions, seed, Git commit SHA, simulation count and
compute resources with every full run. See [`docs/validation-plan.md`](docs/validation-plan.md).

## Development status

See [`docs/known-limitations.md`](docs/known-limitations.md) and
[`docs/two-way-exnex-spec.md`](docs/two-way-exnex-spec.md). The next substantive
research task is to implement and validate the two-way borrowing model without
changing the e-process validity layer.
