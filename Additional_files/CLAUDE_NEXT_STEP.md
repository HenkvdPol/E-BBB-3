# Claude next step: validate the corrected framework study

Work in `HenkvdPol/E-BBB-3`, continuing PR #6 on
`feature/lean-framework-twoway-exnex`. Read the current branch before making changes.

## Read first

1. `framework.tex` and the revised three-page `setup.tex`.
2. `docs/review-2026-11-26.md`.
3. `README.md`, `docs/known-limitations.md` and the latest appended `PROJECT_HISTORY.md` entry.
4. `run_framework_study.R`, `R/simulate.R`, `R/eprocess.R`, `R/borrow_twoway.R`.
5. `tests/test_simulate.R`, `docs/two-way-validation.md` and `docs/validation-plan.md`.

The two-way engine is already implemented. The independent review corrected initial-prior row replication
across simulated trials and single-worker chunking. Do not undo those fixes. At correction commit `d0b8d5b`,
CI reported 69 passes, zero failures and one skipped slow-reference group.

## Objective

Complete a focused validation pass before an approved operating-characteristic pilot.
Keep the existing framework study fixed: the same e-process, three borrowing engines, scaled-Beta benchmark,
10 response patterns and two recruitment designs. Improve only the remaining reproducibility and numerical
validation needed to interpret that study.

## Work

1. **Separate simulation RNG streams.** The current data seed is base + chunk; posterior-fitting seeds add
   the analysis index. This can reuse streams across worker chunks and analyses. Use a small, explicit
   independent-stream scheme for outcomes and fitting. Preserve common potential patient outcomes across
   methods. Document whether reproducibility is defined for a fixed core count or across core counts.
   Add focused tests for the actual invariants; do not require independent MCMC runs to be numerically identical.

2. **Preserve the simulator safeguards.** Check normalized support weights after initial row replication.
   Keep the new batch-versus-single-trial regression and the predictability tests. Include a one-worker
   smoke configuration so more than one trial is processed within a worker.

3. **Check numerical sensitivity on fixed states.** Reuse selected sparse, early, heterogeneous and
   interaction states from `validation/twoway_validation.R`. Compare production defaults with more retained
   MCMC draws and a finer theta/EX-mean reconstruction grid. Report restricted moments, log block-factor
   differences and convergence diagnostics. Keep this a bounded fixed-state check, not a full trial study.
   Do not equate finite-MCMC agreement with unbiasedness or exact continuous-posterior computation.

4. **Prepare the pilot specification.** State which scenarios/engines, replication count, independent
   RNG handling, target Monte Carlo precision and anticipated compute cost will be used. Treat the old
   three-CPU-hour figure as a provisional throughput estimate that excludes fixed overhead.
   Prepare the command/configuration so any approval concerns a concrete run.

5. **Run allowed focused checks.** Syntax parsing and the fast test suite are allowed.
   Only smoke simulations with N_SIM <= 2 are authorized by the existing D9 agreement; pilot and production
   runs require the researcher's explicit approval. Never launch the default N_SIM=1000 study.

## Constraints

- Priors used by the e-process must be normalized probability measures on (p0,1], selected from past
  platform information and frozen throughout the next block.
- Keep the conditional-null and per-basket error guarantee unchanged; no platform-wide FWER/FDR claim.
- Keep `legacy/` unchanged. Do not add design-optimal, Simon or Bayesian decision comparators.
- Do not add partial exchangeability, arbitrary-entry refits, drift, delayed outcomes or adaptive allocation.
- Keep the legacy `redesign = "entry"` guard.
- Keep `setup.tex` and `setup.pdf` synchronized and at most three pages.
- Use clear names and small functions. Avoid making the code shorter by compressing multiple scientific
  operations into hard-to-read array expressions.

## Deliver

A small commit on the existing branch, actual test outputs and a concise validation note. Distinguish
tests executed from checks only proposed, numerical approximation from validity, and measured runtime
from extrapolation. Update the latest status without rewriting historical entries. Present the pilot
specification for approval only after the authorized validation work is complete.

---

**Status update 2026-11-26 (afternoon).** At the researcher's request, the optional **open-ended design** was added:
no maximum sample size, anytime-valid e-value futility, `docs/open-ended-design.md`, framework.tex Proposition 2.
* The fixed 80-cell study is unchanged.
* The validation items above (RNG streams, numerical sensitivity, pilot specification) are still open. They now
  apply to the open-ended variant (`OPEN_ENDED=1`) as well.

---

**Layout change 2026-11-26.**
* The code is now `setting.R`, `simulate.R` and `additional_functions.R` at the top level.
* All other files are in `Additional_files/`; see the mapping in `PROJECT_HISTORY.md`.
* Paths above that name `R/*.R`, `run_framework_study.R`, `tests/`, `validation/`, `docs/` or `legacy/` refer to the
  old layout.
* Follow the researcher's code design: banners per function, argument comments, plain settings in `setting.R`.
* The RNG item is resolved: disjoint outcome and fitting seeds per cell.

---

**eBBB alignment 2026-11-26.**
* The study now uses:
  * the priors of eBBB (T. Jacobs) and our two-way EX/NEX with eBBB's hyperparameters;
  * the 4 x 3 grid;
  * stopping at efficacy or `n_max`;
  * 11 scenarios, including eBBB's four;
  * 154 cells.
* See `docs/ebbb-alignment.md`.
* Do not copy eBBB code into this repository. Compare with `validation/ebbb_compare.R <path to lr_evalue.R>`.
* Still open:
  * approval of the full run;
  * a timing pilot. The unbalanced design has a basket at rate 0.25, so trials last up to 120 time units and the
    two-way cells are the slowest;
  * hyperprior sensitivity.
