# Project history: anytime-valid e-values with Bayesian borrowing for basket and platform trials

A self-contained record of the work done with Claude (11–25 November 2026), so it can be resumed later without the
conversation. The current method and results are in `README.md`; this file records how we got there, what was decided
and what is still open.

## How to resume

* **Code:** `basket_multi.zip`. It holds everything: `Setting.R`, `additional_functions.R`, your library
  `design_optimal_evalues.R`, the design cache `design_banks/`, `checks/`, `README.md` / `README.txt` and this file.
* **Project docs** (claude.ai project 'E-value n-of-1 cohorts', folder `claude/evalue-step3-grid/`): `Setting.R`,
  `additional_functions.R`, `README.md`, `PROJECT_HISTORY.md`.
  * The library and the design cache are only in the zip.
  * Without the cache, missing designs are rebuilt automatically (about 1.5 h on 2 cores, much less on 8).
* **In a new chat in the project:** ask Claude to read `PROJECT_HISTORY.md` and `README.md` first. Claude can also
  search the earlier chats of the project (this conversation included) for details.
* **Working agreements** (from the conversation):
  * `Setting.R` is the only file to run (slurm job); all output goes to `Result/`; the README also as `.txt`.
  * **Do not run simulations without asking.** Claude gives run-time estimates and does individual or exact checks;
    the user runs the full simulations (500–1,000 trials) on the cluster.
  * **`framework.tex` is primary:** it motivates and plans the work; other designs are alternatives or comparators.
  * Files are saved to the project folder; earlier steps are kept in `claude/mexnex-step1/` and `claude/evalue-step2/`.

## Timeline

1. **mEXNEX (11 Nov), `claude/mexnex-step1/`.**
   * Reimplemented Daniells et al. (github.com/LibbyDaniells/mEXNEX): settings, model and simulation files.
   * Overview of π-rules; basket and umbrella set-up; type I error, power and sample size; Figure 1 of the e-value
     paper.
   * The repository's JAGS models use an undefined `sigmae`, so τ was taken as the EX standard deviation with prior
     HN(1).
   * A deterministic quadrature engine over (μ, τ) replaced JAGS; it was checked against JAGS and importance sampling.
2. **Design-optimal e-values (24 Nov), `claude/evalue-step2/`.**
   * Integrated your library (P-max, e-based, GROW).
   * Borrowing via θ1: the point alternative comes from the other cohorts (leave-one-out EX predictive and mEXNEX π).
   * Comparators: Bayesian designs with Δ calibrated under the global null; Simon and the binomial test.
   * Simplified to three files; slurm-ready `Setting.R` with `Result/`.
3. **Grid of settings (24 Nov).**
   * K ∈ {3, 5, 10} × n_max ∈ {20, 30, 40}, with borrowing variants.
   * Design cache with one file per bank: peak memory 3.5 → 1.8 GB.
   * Measured: K = 5, n_max = 30, 1,000 trials = 10.7 min on 2 cores.
4. **Simon's hypotheses (24 Nov).**
   * The 8 hypotheses of Simon's Table 1, with n_max = Simon's optimal n. The Simon search reproduces all 48
     published designs.
   * New measures: sample size to reject and early rejection.
   * Review fixes: θ1 = p1 while the other cohorts have no data; cache keyed by hypotheses; GROW bank with
     borrowing; 'e-Pmax + mEXNEX' removed.
   * 1-trial run of all 17 runs checked (67 min, mostly building 186 designs).
5. **Flexible inclusion (24 Nov).**
   * Acceleration for cohorts that do well, and extension beyond n_max.
   * Finding: for design-optimal e-values an extension decided at n_max is useless, because their bets are
     all-or-nothing near the horizon.
   * Fix: decide at ≤ 3/4 n_max with exact conditional powers, and continue with a re-optimised P-max design.
   * A naive extension of the binomial test inflates the type I error.
6. **framework.tex (25 Nov).**
   * Point-by-point accordance analysis.
   * The framework e-value ('e-mixture') implemented with exact restricted-posterior moments and a closed-form block
     factor.
   * Exact checks of the framework's claims; suggested additions to framework.tex (README section 7).
7. **Framework primary (25 Nov).**
   * Code and README restructured around framework.tex.
   * Flexible framework e-value: extension at n_max by predictive probability.
   * Statistics that cost no computing time (FDR, decisions per trial, bias, evidence), a design-effects table, and the
     trade-off and decision figures.
   * Check scripts kept in `checks/`.

## Main findings so far (exact unless stated)

* **Validity:** every e-value design keeps type I error ≤ α per basket. This holds with borrowing, extension,
  acceleration, adversarial other baskets and drifting null rates (enumeration of all response sequences).
* **The motivation holds:** Bayesian mEXNEX calibrated under the global null reaches type I error 0.15 for null baskets
  next to active ones (10,000 trials, K = 5, n_max = 30).
* **Fixed inclusion, single cohort:** e-based matches Simon's power with Simon's n. Analysed after every patient it
  needs fewer patients under H0 and rejects at about half of Simon's n.
* **Borrowing through θ1** has bought little so far (10,000 trials, K = 5, n_max = 30).
* **Design-optimal e-values:**
  * Extension: power 0.80 → 0.82–0.85 for 1–5 more patients per cohort.
  * Acceleration: time to rejection 14–21% shorter.
* **The framework e-value is conservative at small n:** at n = 10 its type I error is 0.003–0.017 and its power about
  half that of e-based + mEXNEX. Flexible inclusion helps (blocked: 0.17 → 0.25).
* **framework.tex claims:**
  * The restriction to (p0, 1] is essential: unrestricted, P(E ≥ 20 by n = 150) = 0.745 at p = 0.02.
  * Freezing is essential: 0.058 if not frozen.
  * Monte Carlo is valid with frozen draws but needs a rule for H' = 0, which happens in 8–18% of analyses with
    1,000 draws.

## Open items

1. **Full run** on the cluster (`N_SIM=1000`, about 4–4.5 h on 8 cores); then read fig4, fig5 and
   `*_design_effects.csv`.
2. **Two-way grid model** (tumour type × biomarker group) for the framework prior: needs a 2-D basket layout and an
   MCMC, INLA or approximate engine.
3. **DRUP example** in framework.tex.
4. **Edits to framework.tex** suggested in README section 7:
   * a general betting remark;
   * Monte Carlo rules (frozen draws, H' = 0);
   * ties;
   * why (p0, 1];
   * optional continuation;
   * the orange 'True?'.
5. **Possible:**
   * the framework e-value with mEXNEX weights;
   * tuning `cp_min`, `ext_frac` and `m_speed`;
   * Simon's Table 2 (commented in `Setting.R`).

## Files delivered along the way (superseded by `basket_multi.zip`)

* `basket_exnex` (step 1);
* `basket_evalue` and `basket_simple` (step 2);
* `basket_slurm` (single-setting slurm version);
* earlier versions of `basket_multi.zip`.


---

## 2026-11-26 — Research review and repository integration

- Added an explicit supermartingale proof and clarified the conditional-null,
  per-basket multiplicity, allocation and outcome-observation assumptions in
  `framework.tex`.
- Documented that the current framework e-mixture uses a moment-matched
  scaled-Beta prior on `(p0, 1]`, not the exact restricted EXNEX posterior or
  a retained Monte Carlo mixture.
- Added ESS-min to the design-optimal wrapper, active simulation settings and
  exact sequence checks.
- Disabled the experimental `redesign = "entry"` design in active set-ups.
  The current calendar-time simulator does not represent cohort entry as an
  independent event; selecting this option now stops with an explanatory error.
- Added known-limitations, two-way EXNEX specification and validation-plan
  documents, plus a scoped Claude task prompt.
- Added GitHub Actions focused validation (R parsing and `fw_engine.R`) and a
  manually dispatched full check suite. No R checks or full simulations were
  executed as part of this review because R is unavailable in the authoring
  environment.

---

## 2026-11-26 — Lean framework code, framework-only study and two-way EXNEX (branch `feature/lean-framework-twoway-exnex`)

**Decisions** (confirmed by the researcher: D1–D4; documented defaults: D5–D9; see `docs/implementation-plan.md`):
- **D1.** New, short framework code in `R/`. The earlier code base moved unchanged to `legacy/`: same output of
  `legacy/checks/fw_engine.R` before and after the move.
- **D2.** Two-way hyperpriors match the one-way benchmark: sigma_a, sigma_b, sigma_e ~ HN(1/sqrt 3), so the total EX
  variance has the same expectation as tau² with tau ~ HN(1). Same mu prior, w = 0.5, same NEX.
- **D3.** Binary EX/NEX per basket (framework.tex eq. 9). Row-only and column-only partial exchangeability later,
  once the component weights are agreed.
- **D4.** The e-process gets the restricted posterior itself on a theta grid aligned at logit(p0) (framework.tex
  eq. 8 up to quadrature). The scaled Beta stays as a benchmark design.
- **D5–D9.**
  - Random effects without a sum-to-zero constraint.
  - No external data D0.
  - Two-way production computation by MCMC, validated against deterministic quadrature and JAGS.
  - Only smoke runs until approval.

**Code.**
- `R/eprocess.R`: validity layer.
- `R/borrow_oneway.R`: none and one-way EXNEX.
- `R/borrow_twoway.R`: two-way EXNEX.
- `R/reference_twoway.R`: references.
- `R/simulate.R`, `R/summarise.R`, `R/plots.R`.
- `run_framework_study.R`: single entry point; 3 x 3 grid, 10 scenarios, balanced and unbalanced recruitment, four
  priors.
- `tests/`: base-R runner.
- `validation/twoway_validation.R`.

**Verified** (executed on R 4.3.3):
- One-way port = legacy `posterior_exnex()`: restricted mean, second moment and P(p > p0) to < 1e-15. New 0.05 grid
  vs legacy: 3e-5.
- Exact type I error of a null basket, with all 6^6 block paths and adversarial other baskets, at alpha = 0.05:

  | Prior | Exact type I error |
  |---|---|
  | none | 0.0073 |
  | one-way EXNEX | 0.0060 / 0.0067 |
  | one-way EXNEX, scaled Beta | 0.0059 |
  | two-way EXNEX | 0.0059 / 0.0067 |

  All are ≤ alpha. The framework e-value remains conservative at n ≤ 30.
- Predictability, checked directly: outcomes that differ only inside a block leave that block's prior unchanged.
- Two-way posterior accuracy (`docs/two-way-validation.md`): the quadrature and the equivalent one-way model agree
  to 2e-4.
  - Production sampler vs references: typically 0.0004–0.0007 in the restricted mean and 1–2% in the 5/5 block
    factor for baskets with n ≥ 15.
  - Up to 0.009 and 12–15% for empty or 5-patient baskets (Monte Carlo error, no bias).
- Smoke run of the study (`N_SIM=2`): all 80 cells run (9 min on 2 cores).

**Findings.**
- Borrowing from strongly responding baskets *lowers* the exact type I error of a null basket: the prior moves away
  from p0. The hardest platform is one with response rates just above p0.
- A quadrature window around the data underestimates the all-NEX configurations under the wide mu prior, and so
  overstates P(EX). References must cover the whole prior range.

**Open.**
- Full study run (estimate: 3 CPU-hours upper bound, 25–35 min on 8 cores): **needs approval**.
- Partial exchangeability; hyperprior sensitivity.
- Entry-time events (proposal: `docs/entry-time-redesign.md`; the guard stays).
- Drift, delays and allocation scenarios.
- Moving the design-optimal comparators and flexible inclusion onto the new code.

---

## 2026-11-26 — Independent review of the simplified framework study

The entries above are preserved as the historical record. This entry qualifies the current implementation and
evidence after reviewing PR #6 at `a8a915f`.

- Confirmed that the active study is limited to the framework e-process with no borrowing, one-way EX/NEX,
  two-way EX/NEX and the scaled-Beta representation benchmark. Earlier executable code is unchanged in `legacy/`.
- Found and corrected a first-prior batching defect in `R/simulate.R`: reshaping before row replication moved
  weights between support points and could break normalization when multiple trials were simulated per worker.
  The corrected code repeats prior rows directly and checks the assembled matrix.
- Added three meaningful regressions: direct initial-factor agreement, correct initial prior means, and
  fixed-outcome agreement between batched and separate one-trial execution.
- Fixed one-trial/one-worker study chunking with `parallel::splitIndices()`.
- The correction commit `d0b8d5bb4a085bd7e84dfa7b727029ca1a604af6` passed GitHub Actions
  [run 45805989231](https://github.com/HenkvdPol/E-BBB-3/actions/runs/45805989231):
  39 R files parsed; 69 PASS, 0 FAIL, 1 skipped slow-reference group. The original reviewed-head CI had
  66 PASS, 0 FAIL, 1 SKIP. The previously reported 70/70 slow-suite was a separate author-run result.
- Renamed the platform `abs_bias` output to `mean_abs_error` to describe its calculation accurately;
  per-basket signed `post_bias` is unchanged.
- Rewrote `setup.tex` and `setup.pdf` as a maximum-three-page description of the current study. Synchronized
  the computational remark in `framework.tex`, README files and numerical/validation descriptions.
- Clarified that finite grid-factor summation is direct, while posterior fitting is approximate; the two-way
  hyperpriors match expected EX variance before discretization, not the full one-way prior law.
- Removed unsupported claims of unbiased finite MCMC and of a guaranteed three-CPU-hour upper bound.
- Independently recomputed the saved validation summaries from 1,216 rows/22 data states. The reported
  comparison maxima and Monte Carlo spreads were reproduced to rounding. These fixed-state checks are
  distinct from an end-to-end simulation study.
- R was unavailable in the local review environment. Executed R evidence was read from actual CI logs;
  the batching reproduction and saved-summary calculations used Python.
- Before production, separate outcome/fitting RNG streams and check sparse-state numerical sensitivity.
  Keep the study scope fixed; no new comparator families or event-time extensions are needed for that step.
- No pilot or full operating-characteristic simulation was run. Approval for larger runs remains required
  under the existing D9 agreement.

Full review: `docs/review-2026-11-26.md`. The next scoped task is in `CLAUDE_NEXT_STEP.md`.

---

## 2026-11-26 (afternoon) — Open-ended baskets: no maximum sample size

**Request.** Check whether a trial that potentially never stops (an infinite sample size) is still part of the
design, and add it if not.

**Finding.**
* framework.tex already allows it: Proposition 1 bounds the probability of rejecting at *any* analysis over an
  unbounded patient sequence.
* The new code did not: the simulator stopped baskets at a finite `n_max`, and the engines' tables were built only
  up to `n_max`.

**Added** (the fixed 80-cell study is unchanged; ChatGPT's fixes `d0b8d5b` and `075b954` are kept):
* `n_max = Inf` in the simulator, with a simulation horizon `n_cap` that only ends the computation. Baskets still
  open there are reported as `open`.
* Engine tables that grow on demand (`state_table()`).
* Optional anytime-valid futility.
  * E^f tests p >= p1 with the posterior restricted to [0, p1).
  * Stop at E^f >= 1/beta; P(futility | p >= p1) <= beta (framework.tex, new Proposition 2).
* The remarks 'No maximum sample size' and 'Open-ended design' in framework.tex. Proposition 1 and its proof are
  unchanged.
* `R/open_ended.R`: the exact dynamic programme, based on the telescoping of the factors for a fixed starting
  prior.
* `validation/open_ended_feasibility.R` and its results.
* `OPEN_ENDED=1` option in `run_framework_study.R`.
* `docs/open-ended-design.md`.
* A paragraph in `setup.tex`, which stays at 3 pages.

**Verified** (exact, p0 = 0.10, p1 = 0.30, alpha = 0.05, beta = 0.20, up to 3000 patients; 7 starting priors
including one-way and two-way borrowing from completed baskets):
* largest P(reject | p <= p0) = 0.023;
* largest P(futility | p >= p1) = 0.160;
* with e-value futility every basket stops, by n = 125 at the latest.

**Compared with the fixed design (n_max = 30), open-ended with e-value futility:**

| Measure | Fixed n_max = 30 | Open-ended, e-value futility |
|---|---|---|
| E[N] at p0 | 30 | 20–26 |
| Power at p1 | 0.49–0.62 | 0.90–0.95 |
| E[N] at p1 | about 23 | 27–30 |

Without futility, null baskets never stop with probability about 0.98.

**Tests.** `tests/run_tests.R --slow`: 94 PASS, 0 FAIL.

**Open.** Open-ended operating characteristics with borrowing from still-recruiting baskets need a simulation
pilot, which requires approval. ChatGPT's validation items from `CLAUDE_NEXT_STEP.md` also remain open (RNG
streams, grid and MCMC sensitivity, pilot specification).

---

## 2026-11-26 — Restructure in the researcher's code design

**Request.** One file with the functions, one with the simulation, one with the settings; everything else in a
single folder. Follow the code design of HenkvdPol/Class_Imbalance_Prediction (`Setting/setting.R`, `simulate.R`,
`additional_functions.R`).

**New layout.** At the top level only `setting.R`, `simulate.R`, `additional_functions.R`, `README.txt`, `setup.pdf`
and `Additional_files/`. The `.github/` and `.gitignore` files stay at the top level, as GitHub requires.

| Old path | New path |
|---|---|
| `R/eprocess.R`, `R/borrow_oneway.R`, `R/borrow_twoway.R`, `R/simulate.R`, `R/summarise.R`, `R/plots.R` | `additional_functions.R` (Parts 1–4) |
| `run_framework_study.R` | `setting.R` (plain settings) + `simulate.R` (`genSim()` per cell, `foreach`/`doParallel`, command-line arguments, as in Class_Imbalance_Prediction) |
| `R/reference_twoway.R`, `R/open_ended.R` | `Additional_files/validation/reference_functions.R`, `open_ended_functions.R` |
| `framework.tex/.pdf`, `setup.tex`, `docs/`, `tests/`, `validation/`, `legacy/`, `scripts/parse_all.R`, `PROJECT_HISTORY.md`, `CLAUDE_NEXT_STEP.md` | `Additional_files/…` |
| `README.md` | `README.txt` |
| — | new `Additional_files/slurm/sim.slurm` and `parameters.txt` (job array: fixed and open-ended design) |

**Code style.** Every function has a `####` banner with a title and argument comments, and one step per line.

**Unchanged computations.** All 20 outputs of a regression battery are bit-identical to the old code, including the
two-way sampler with a fixed seed: posteriors of all engines, simulated trials (fixed, open-ended, futility
variants), the exact dynamic programme, the quadrature and the summaries.

**Changed: seeds.** `simulate.R` runs every cell's trials in one process, with outcome seeds and fitting seeds from
disjoint ranges. This resolves the RNG item of the review, and results no longer depend on the number of cores.

**Tests.** `Additional_files/tests/run_tests.R --slow`: all pass (see the commit message for counts). Smoke run of
`simulate.R`, both horizons, with 2 trials per cell.

---

## 2026-11-26 — Aligned with eBBB (T. Jacobs): priors, stopping rule, hyperparameters, grid

**Reference.** The eBBB code (`tijn-jacobs/eBBB`, `code/`), received as a zip file. It was read only and run
locally for comparison; it is not copied into this repository, and the eBBB repository was not changed.

**Researcher's decisions.**

| Item | Decision |
|---|---|
| Priors | Adopt eBBB's set (no borrowing, power equal and column, conjugate EX/NEX, similarity, oracle); our two-way EX/NEX replaces its Stan fit |
| Stopping rule, fixed design | Efficacy or `n_max`, as in eBBB; no futility stop |
| Open-ended design | Kept as an option (e-value futility) |
| Two-way hyperparameters | As in eBBB: `mu ~ N(logit p0, 2^2)`, `sigma_a, sigma_b ~ HN(1)`, `sigma_eps ~ HN(0.5)`, NEX `N(logit p0, 2^2)`, `p_ex = 0.5` |
| Grid | 4 tumour types x 3 biomarker groups |
| Kept from E-BBB-3 | Calendar-time recruitment, global analyses every 5 time units, `p0 = 0.10`, `p1 = 0.30`, `n_max = 30`, our scenarios (rebuilt for 4 x 3, with eBBB's four added) |

**Code.**
* `additional_functions.R`, new functions:
  * `beta_cells()`: a Beta prior on the theta grid, with the tail mass in the end cells;
  * `jsd_beta()`;
  * `add_points()`: exact oracle points;
  * `engine_ebbb()`: the five eBBB priors.
* `additional_functions.R`, changed functions:
  * `twoway_prior()`: separate SD scales and a general NEX prior. The defaults give bit-identical posteriors.
  * `hn_grid()`: reaches 4 x the scale.
  * `make_design()`: the fixed horizon has `futility = "none"`.
  * `make_engine(model, layout, p_true)`: the seven study models; the earlier models remain available.
  * Basket names are `T<j>-B<l>`.
* `setting.R`: J = 4, L = 3, `Models`, the eBBB prior settings, `twoway_hyper`, recruitment designs and 11 scenarios
  on 4 x 3. The study has 154 cells.
* `simulate.R`: passes the true rates to `make_engine()` for the oracle.

**Validation.**
* **eBBB priors.** Block factors match eBBB's closed forms within 2.3e-4 (relative); the EX weights and the oracle
  are exact (`validation/ebbb_compare.R`).
* **Two-way with the eBBB hyperparameters.** On 4 x 3 it agrees with JAGS within the Monte Carlo error of the
  production sampler (`validation/twoway_ebbb_hyper.R`).
* Details: `docs/ebbb-alignment.md` and the addendum to `docs/two-way-validation.md`.
* **New tests.** `tests/test_ebbb.R` has 40 tests, including the exact type I error of every eBBB prior.

**Docs.** Updated: `setup.tex`/`setup.pdf` (sections 3–4, still three pages), `README.txt`,
`docs/known-limitations.md` (item 12: differences that remain), `docs/open-ended-design.md`.


## 2026-11-26: sanity review after the three-file simplification

Reviewed main `ec57509383cad1c8007bffde114e681aafbd3673`, including the eBBB alignment.
The requested three-file structure and the familiar settings / genSim / per-cell output
pattern are in place. The earlier initial-prior replication correction and separate
outcome/fitting seed ranges remain intact.

**Necessary correction.** The first additive-scenario basket is intended to be exactly
`p0`. In R 4.4.2, its logit round trip is slightly above `p0`, so it was classified as
active and the oracle selected the wrong bet. `setting.R` now explicitly anchors that
one basket to `p0`; two regression checks preserve exact equality and one null basket
in each recruitment design. No global comparison tolerance was added.

**Executed validation.** Code commit `1f7dd9031db438f075d3042b1ed7d7596805ebd9` passed
[GitHub Actions run 45855905668](https://github.com/HenkvdPol/E-BBB-3/actions/runs/45855905668):
41 R files parsed; **132 PASS, 0 FAIL, 1 SKIP**. The log directly records
`p0 = 0.10000000000000001` versus an unanchored round trip of
`0.10000000000000002`. The legacy numerical checks also passed.
The slow reference group was skipped; the top-level study driver was parsed, not
executed, in this CI run. No new pilot or full simulation was launched.

**Documentation.** Corrected the oracle description to a truth-informed reference
(it is not a universal finite-horizon power upper bound), qualified grid precision as
observed in tested states, and separated earlier fixed-prior open-ended feasibility
results from the current 4 x 3 study. Updated the short setup's validation evidence
and obsolete runtime estimate. The rebuilt `setup.pdf` has three pages and was
visually checked. The source remains in `Additional_files/setup.tex`.
