# Open-ended baskets: no maximum sample size

Status 2026-11-26. Optional design added on branch `feature/lean-framework-twoway-exnex`. The fixed 80-cell study
(`simulate.R` (settings in `setting.R`) with defaults) is unchanged.

## 1. Question and short answer

**Question.** Can a basket be designed without a maximum sample size, so that the trial may never stop?

**Methodology: yes, already.**
* framework.tex defines the patients i = 1, 2, … and the analysis times τ_{c,1} < τ_{c,2} < … without an upper
  limit.
* Proposition 1 bounds the probability that E_c reaches 1/α at *any* analysis (Ville's inequality).
* A basket may therefore run forever, and any maximum, budget or calendar end is optional stopping.

**Code: it was not implemented before this change.**
* The new simulator stopped every basket at a finite `n_max`, and the engines' tables were built only up to `n_max`.
* The legacy code allowed finite extensions (`n_ext`) only.

**Now (2026-11-26).**
* `n_max = Inf` is allowed.
* An anytime-valid futility rule has been added (framework.tex, Proposition 2).
* The engines' tables grow with the sample size.
* An exact feasibility analysis is in `Additional_files/validation/open_ended_feasibility.R`.

## 2. Method

**Efficacy (unchanged).** Reject H0,c: p_c ≤ p0 at the first analysis with E_c ≥ 1/α.

**Futility (new, optional).** This is the mirror image of the efficacy test.
* E^f_c tests H1,c: P(Y_i = 1 | G_i) ≥ p1 with block factors
  ∫ ∏ (q/p1)^Y ((1 − q)/(1 − p1))^(1 − Y) dπ^f(q).
* The prior π^f is the borrowing posterior restricted to [0, p1), fixed at the analysis and frozen for the block.
* Stop for futility at the first analysis with E^f_c ≥ 1/β.
* Under H1,c, P(futility stop) ≤ β at any time (framework.tex, Proposition 2). The proof is that of Proposition 1
  with the inequality reversed.

**The combined open-ended design.** Stop at whichever comes first: rejection or futility. There is no maximum.

**Guarantees.**
* Each rule is optional stopping for the other process. So, whether or not the borrowing model is correct:
  * P(reject) ≤ α under H0,c;
  * P(futility) ≤ β under H1,c.
* For a constant p_c and priors that concentrate, p0 < p1 means every p_c lies above p0 or below p1. At least one
  process then grows, and the basket stops after finitely many patients with probability one.
* Then P(reject) ≥ 1 − β for p_c ≥ p1: **a power guarantee without a maximum sample size**.
* Rates between p0 and p1 may end in either decision.

**Without a futility rule** (efficacy only): a test of power one.
* For p_c > p0 the basket is rejected eventually with probability one.
* Under H0 it runs forever with probability at least 1 − α. In practice an open-ended basket needs a futility rule,
  a budget or a calendar end.

**Telescoping.** If the prior of every block is the restricted posterior from one fixed starting distribution π0,
then ∏ e_k is the mixture likelihood ratio ∫ (p/p0)^y ((1 − p)/(1 − p0))^(n − y) dπ0,R(p), where π0,R is π0
restricted to (p0, 1].
* This holds without borrowing, and with borrowing from baskets whose data no longer change.
* E (and E^f) then depend on (y, n) only, so the operating characteristics are computed exactly by a dynamic
  programme over y (`Additional_files/validation/open_ended_functions.R`).
* `Additional_files/tests/test_open_ended.R` checks the identity in the simulator (to 9e-16). It also checks the dynamic programme
  against a brute-force enumeration of the block-by-block recursion (equal to 1e-12).
* With borrowing from baskets that are still recruiting, the factors do not telescope, and simulation is needed.

## 3. Feasibility: exact results

`Additional_files/validation/open_ended_feasibility.R` (20 s; results in `Additional_files/validation/results/open_ended_*.csv`).
* One basket; p0 = 0.10, p1 = 0.30, α = 0.05, β = 0.20; followed up to 3000 patients.
* Seven starting priors:
  * none;
  * one-way EXNEX borrowing from 8 completed baskets: all null, all active, or mixed;
  * two-way EXNEX on a 3 x 3 grid: own row active; own row and column null; all other baskets null.
* Analyses after every 5 patients ("blocks of 5") or after every patient.

**Guarantees, checked exactly over every prior, design and p:**

| Check | Result | Bound |
|---|---|---|
| largest P(reject) for p ≤ p0 | 0.023 | ≤ α = 0.05 |
| largest P(futility stop) for p ≥ p1, e-value futility | 0.160 | ≤ β = 0.20 |
| largest P(still open after 3000 patients), e-value futility | 0 | every basket has stopped by n = 125 |

**No borrowing, blocks of 5.**

| Design | P(reject), p = 0.10 | E[N], p = 0.10 | P(reject), p = 0.20 | E[N], p = 0.20 | P(reject), p = 0.30 | E[N], p = 0.30 |
|---|---|---|---|---|---|---|
| fixed n_max = 30 (current study) | 0.006 | 29.9 | 0.18 | 28.2 | 0.62 | 22.8 |
| open-ended, no futility | 0.015 (98.5% never stop) | ≥ 2959 | 1.00 | 107 | 1.00 | 32 |
| open-ended, e-value futility | 0.009 | 23.6 | 0.47 | 45 | 0.93 | 28.9 |
| open-ended, e-value futility, every patient | 0.022 | 15.5 | 0.41 | 29.5 | 0.84 | 20.5 |

**With borrowing (open-ended, e-value futility, blocks of 5).**

| Starting prior | P(reject), p = 0.10 | E[N], p = 0.10 | P(reject), p = 0.30 | E[N], p = 0.30 |
|---|---|---|---|---|
| one-way, 8 null baskets | 0.007 | 19.7 | 0.90 | 29.8 |
| one-way, 8 active baskets | 0.012 | 25.7 | 0.95 | 27.9 |
| two-way, own row active | 0.010 | 19.8 | 0.90 | 26.9 |
| two-way, own row and column null, rest active | 0.009 | 20.0 | 0.90 | 27.3 |

The 90% quantile of N is at most 85 patients with blocks of 5, and at most 71 with analyses after every patient.

**P(reject by n) without futility, p = p0 = 0.10, no borrowing.** The bound α = 0.05 is never approached.

| n = 30 | n = 100 | n = 300 | n = 1000 | n = 3000 |
|---|---|---|---|---|
| 0.0062 | 0.0094 | 0.0118 | 0.0138 | 0.0150 |

**Reading.**
1. **Feasible and valid.** Validity holds over the whole horizon. With e-value futility, every basket stops: by
   125 patients at the latest in these settings, with a typical 90% quantile of 40–85.
2. **Efficient.** Against the current fixed design (n_max = 30), the open-ended design with e-value futility needs
   fewer patients at p0 (E[N] 20–26 against 30) and has much more power at p1 (0.90–0.95 against 0.49–0.62), at a
   similar E[N] under p1. A fixed n_max of 30 is too short for 80% power at p1 = 0.30 with this e-value; the
   open-ended design lets the data decide.
3. **Borrowing matters more without a maximum.** With n ≤ 30 the integer rejection boundaries of different priors
   often coincide (e.g. none and "8 active baskets": reject from y = 4, 5, 6, 7, 8, 9 at n = 5, …, 30), so borrowing
   changes little. Beyond 30 the boundaries differ.
4. **Analysing after every patient** brings the type I and II errors closer to their bounds (0.022 and 0.16) and
   lowers E[N]. Blocks of 5 are more conservative.
5. **Two-way starting priors.** Their predictive was computed with a long MCMC run (4 x 20,000), so the dynamic
   programme is exact for that (approximate) starting prior. Without borrowing and for one-way EXNEX it is exact.

## 4. Implementation

| Where | What |
|---|---|
| `additional_functions.R` (Part 1) | `futility_side()`, `restrict_below()`, `check_priors_below()`. The futility factor reuses `log_factor()`. The efficacy functions are unchanged. |
| `additional_functions.R` (Part 2) | `state_table()`: (n, y) tables built up to `n_max` at once and extended on demand. Results for n ≤ `n_max` are identical to before. |
| `additional_functions.R` (Part 3) | `n_max = Inf` with a required simulation horizon `n_cap` (only to end the computation; baskets still recruiting there are reported as `open`, not as decided); `futility = "unreachable"` (default, as before), `"evalue"` (with `p1`, `beta`) or `"none"`; new output columns `open` and `logEf` |
| `Additional_files/validation/open_ended_functions.R` | `exact_basket()`: the exact dynamic programme; `mixture_logE()` |
| `additional_functions.R` (Part 4) | `open` per basket; `open_any` per platform. E[N] is truncated at `n_cap`. |
| `simulate.R`, `setting.R` | `Rscript simulate.R open 1000 20261126 TRUE` runs the same 80 cells with `n_max = Inf`, e-value futility (β = 0.20) and `n_cap = 150`, written to `Result/open/`. The default run is unchanged. |
| `Additional_files/tests/test_open_ended.R` | 21 tests: futility layer, growing tables, telescoping, the dynamic programme against path enumeration, the exact guarantees, simulator statuses |

## 5. Practical points and limitations

* **A futility rule or another stop is needed.** Without one, null baskets run forever with probability ≥ 1 − α.
  Alternatives to e-value futility (a budget, a calendar end, a posterior rule) are all valid for α. Only the
  e-value rule carries the β guarantee.
* **Choice of p1 and β.** p1 now plays a real role: it is the boundary of the futility hypothesis. β = 0.20 follows
  the earlier designs. Rates between p0 and p1 may end either way.
* **Long trials.** The conditional null allows drift, but over a long period the clinical question itself (p0, the
  standard of care) may change. This is a protocol question, not a validity question.
* **Platform-wide error.** Per-basket α holds for every basket however long it runs. Platform-wide FWER is not
  controlled, and in an open platform that keeps adding baskets it is even less meaningful.
* **Numerics.**
  * The theta grid has cells of width 0.05. At n ≈ 1000 the posterior sd of theta is about 0.07, so the grid
    becomes coarse. Validity is unaffected; a finer grid (`theta_grid(p0, h = 0.02)`) costs about 2.5 times more.
  * Tables grow with the data states seen (memory about 80 KB per (n, y) row for the two-way engine).
* **Simulation.** With borrowing from baskets that are still recruiting, operating characteristics need simulation.
  Only a smoke run (`Rscript simulate.R open 2 20261126 TRUE`) has been executed; a pilot needs approval. Open-ended trials have more
  analyses, so the two-way engine (about 1 s fixed cost per analysis per worker) dominates the run time.
