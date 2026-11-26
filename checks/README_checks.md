# Check scripts (individual and exact checks, no simulation)

These are the checks behind section 9 of the main README, kept so they can be rerun after a change.

**How to run.**
* Run any script with `Rscript checks/<name>.R` from any directory. Paths are relative to the package.
* Output (.rds, figures, small design caches) goes to `checks/output/`.
* `load_settings.R` is a helper: it reads `Setting.R` up to the run part, without running a simulation.

Times are rough, for 2 cores, with `design_banks/` present.

| Script | What it checks | Time | Result when last run |
|---|---|---|---|
| `simon_check.R` | The Simon search reproduces all 48 designs of Simon (1989), Table 1 (optimal and minimax) | ~2 min | all identical |
| `exact_checks.R` | Trial engine on all 2^n response sequences vs the library's exact `evaluate_design` for GROW, P-max, ESS-min and e-based; exact single-cohort results at Simon's hypotheses | ~3 min | historical run: equal; type I ≤ α |
| `flex_exact.R` | Accelerated and extended design-optimal e-values, all 2^18 sequences; an independent one-sequence implementation compared with the engine | ~5 min | type I ≤ α (max 0.0487); 0 mismatches |
| `simon_flex.R` | Extension at Simon's hypotheses (exact grid recursion); naive extension of the binomial test | ~3 min | power 0.80 → 0.82–0.85; binomial up to 0.054 |
| `fw_engine.R` | Framework e-value: closed-form block factor vs integration; restricted-posterior moments vs direct grid; matched prior vs exact restricted posterior | ~1 min | 5e-12; < 1e-15; −1.4% to +8% |
| `fw_validity.R` | Two cohorts, all 2^15 sequences of cohort 1, adversarial paths of cohort 2; framework and borrowing designs | ~5 min | type I ≤ α (max 0.0484) |
| `fw_claims.R` | framework.tex claims: frozen prior, Monte Carlo with frozen and re-made draws, H' = 0, design-optimal bets inside blocks | ~3 min | not frozen 0.058; re-made draws 0.051; H' = 0 in 8–18% |
| `fw_claims2.R` | Restriction to (p0, 1] (horizon 150); drifting and history-dependent null rates | ~2 min | unrestricted 0.745 at p = 0.02; drift ≤ α |
| `fw_flex.R` | Flexible framework e-value: adversarial other cohort and drift | ~3 min | type I ≤ α (max 0.017) |
| `grow_check.R` | e-GROW bets = Kelly fraction, so e-GROW = framework factor with a point-mass prior | < 1 min | 4.7e-13 |
| `bank_test.R` | Design banks with borrowing for GROW and e-based (bets per θ1) | < 1 min | correct |
| `regression.R` | New engine vs the original engine (`old_versions/additional_functions_old.R`) for all design-optimal and Bayesian designs, hand-made data | ~2 min | 840/840 identical |
| `regression_mixture.R` | Framework engine vs its previous version (`old_versions/additional_functions_v3.R`); predictive probabilities of an extension sum to 1 | ~3 min | 126/126 identical; 8e-15 |
| `unit_flex.R` | Every design × schedule × scenario on hand-made data: decision rules, n and time rules, new statistics, design effects, all figures | ~10 min | 0 problems |
| `check_results.R` | **After a run of `Setting.R`:** checks every row and file in `Result/` (counts, ranges, NA patterns, Simon values, curves) | < 1 min | 0 problems (1-trial run) |
| `build_ext.R` | Builds all design banks the 17 runs need (normally not needed: `design_banks/` is complete) | minutes | – |
| `clean_cache.R` | **Deletes** files in `design_banks/` that no current run needs. Only use after changing settings on purpose. | < 1 min | – |

`old_versions/` holds the two earlier versions of `additional_functions.R` that the regression checks compare with.


## Changes requiring a fresh run

The `essmin` design type has now been added to `make_evalue_design()` and to
`Setting.R`; `exact_checks.R` now enumerates it alongside the other design
families. The historical timings and results above predate this addition.
Run `Rscript checks/exact_checks.R` and the full validation workflow before
treating the ESS-min comparison as validated.

The `redesign = "entry"` variant is no longer in active set-ups and the
simulator now stops explicitly if it is selected. A correct implementation
requires event-driven cohort-entry handling; see `docs/known-limitations.md`.
