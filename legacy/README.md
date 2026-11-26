# legacy/: the full comparison code (unchanged)

This folder holds the earlier code base. It was moved here unchanged on 2026-11-26 (branch
`feature/lean-framework-twoway-exnex`), so the new framework code in `R/` can stay short.

| File | Content |
|---|---|
| `Setting.R` | All settings and the complete run of the earlier study. Run with `Rscript legacy/Setting.R`. |
| `additional_functions.R` | One-way EXNEX engine; framework design with a scaled-Beta prior; design-optimal designs; Bayesian comparators; flexible inclusion; figures. |
| `design_optimal_evalues.R` | GROW, P-max, ESS-min and e-based designs (Baas et al.). |
| `checks/` | Exact and regression checks. Run each one with `Rscript legacy/checks/<name>.R`; see `checks/README_checks.md`. |

**What it is used for now:**
* **Benchmark.** `tests/test_oneway.R` checks that the one-way EXNEX engine in `R/` reproduces the restricted
  posterior moments of `posterior_exnex()` here.
* **The later extras.** Design-optimal comparators, flexible inclusion, Bayesian comparators and Simon will be built on
  the new code later.

All paths inside are relative to the file's own folder, so the move changed no results.
`Rscript legacy/checks/fw_engine.R` gives the same output before and after the move (verified 2026-11-26).

The guard on `redesign = "entry"` in `run_trials_evalue()` stays (see `docs/entry-time-redesign.md`).
* `README_legacy.txt`: the detailed README of this code base (the former top-level `README.txt`).
