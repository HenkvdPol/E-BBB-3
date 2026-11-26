# E-BBB-3: Anytime-valid e-values for basket and platform trials

Research code for a clinical-trial design project using e-values with Bayesian
borrowing across baskets. The primary method is specified in
`framework.tex`; the single-arm design-optimal e-value comparators are based
on Baas, ter Schure & van Rosmalen (2026),
[arXiv:2605.28653](https://arxiv.org/abs/2605.28653).

## Method hierarchy

1. **Primary:** framework e-value: at each analysis, freeze a prior supported
   on `(p0, 1]`, mix the likelihood ratio for the next block, multiply block
   factors, and reject basket `c` when `E_c >= 1/alpha`.
2. **Comparators within the e-value family:** design-optimal P-max, e-based
   and GROW policies, optionally with a borrowed point alternative.
3. **External comparators:** Bayesian EXNEX posterior-probability rules,
   Simon two-stage designs and exact binomial tests.

The type-I error guarantee is per basket. It is not automatically a
family-wise error guarantee across the entire platform.

## Repository status

The repository currently contains the research README, `setup.tex` and
`.gitignore`. The complete source-only project archive is supplied separately
as `E-BBB-3-source-only.zip`; import its contents into this clone to add the
R implementation and checks. The archive excludes generated `design_banks/`
and `Result/` files because these are reproducible caches and outputs.

## Project files in the source archive

- `framework.tex`: research framework and validity argument.
- `setup.tex` / `setup.pdf`: concise project and reproducibility set-up.
- `Setting.R`: settings and full simulation entry point.
- `additional_functions.R`: borrowing, e-value construction, simulation and summaries.
- `design_optimal_evalues.R`: design-optimal e-value library.
- `checks/`: exact, validity, regression and unit checks.
- `PROJECT_HISTORY.md`: development decisions and open research items.

## Current scope and known limitation

The uploaded implementation currently uses a one-way EXNEX model over a list
of cohorts. The proposed two-way tumour-type by biomarker-group model is not
implemented yet. That is the main methodological development still required.
The DRUP motivating example and final literature references also remain open.

## Reproducibility

Run checks before any full simulation. From the project root:

```sh
Rscript checks/fw_engine.R
Rscript checks/fw_validity.R
Rscript checks/exact_checks.R
Rscript checks/regression.R
N_SIM=1 Rscript Setting.R       # smoke test
N_SIM=1000 Rscript Setting.R    # full run, only after approval
```

The project history says not to launch full simulations without the
researcher's approval. Missing design-bank caches are rebuilt automatically;
they can require substantial compute time. Record R/package versions, seed,
commit SHA and compute resources with every run.
