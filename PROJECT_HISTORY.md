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
