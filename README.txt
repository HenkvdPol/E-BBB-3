===================================================================================
ANYTIME-VALID TESTING WITH BAYESIAN BORROWING FOR PLATFORM TRIALS: SIMULATION STUDY
===================================================================================

This code implements and evaluates the test of framework.tex ('Anytime-valid testing with Bayesian
borrowing for platform trials'). The framework sets the plan; everything here follows from it.

The test (framework.tex):
- Each basket c tests H0,c: p_c ≤ p0 against p_c ≥ p1 (binary endpoint, one-sided α = 0.05 per
  basket, no multiplicity correction).
- The evidence is an e-value. At every analysis the prior for the next block of patients is frozen:
  the posterior of p_c from a Bayesian borrowing model fitted to all baskets (EXNEX), restricted to
  (p0, 1].
- The block factor is the likelihood ratio of the new block mixed over that prior. The e-value is
  the product of the block factors, and H0,c is rejected at the first analysis with e-value ≥ 1/α =
  20.
- The type I error is ≤ α at every sample size, however many analyses are done and whatever the
  other baskets show. This holds whether or not the borrowing model is correct, and under optional
  stopping and continuation.

What the simulation compares:

A. The framework e-value
  Designs:  e-mixture (framework), e-mixture (framework, flexible), e-mixture (no borrowing)
  Role:     Primary. The test of framework.tex, with optional continuation, and without borrowing to
            measure what borrowing adds.

B. Design-optimal e-values
  Designs:  e-GROW, e-Pmax, e-based (Baas, ter Schure & van Rosmalen, arXiv:2605.28653
            (https://arxiv.org/abs/2605.28653)), with a borrowed point alternative θ1 (mEXNEX), and
            flexible variants
  Role:     Alternatives within the framework: a point-mass prior on (p0, 1] with bets optimised for
            a maximum sample size.

C. Comparators
  Designs:  Bayes Umbrella, Bayes mEXNEX (posterior probability > Δ), Simon's two-stage design,
            exact binomial test
  Role:     Outside the framework: no error guarantee under borrowing or flexible inclusion (the
            motivation of framework.tex).

1. The plan of framework.tex and where the code stands
======================================================

A test per basket with its own α; type I error ≤ α at any sample size, whatever the other baskets do
  In the code:  All e-value designs; Bayesian designs and Simon as comparators
  Status:       Done. Checked exactly, also with adversarial other baskets and drifting response
                rates (section 9).

Interim analyses at any time; optional stopping and continuation
  In the code:  Three schedules (sequential, blocked, two-stage); futility stopping; extension
                beyond n_max; faster inclusion for cohorts that do well. For the framework e-value:
                e-mixture (framework, flexible)
  Status:       Done

e-value construction: frozen prior on (p0, 1], block mixture likelihood ratio, product, reject at the first crossing
  In the code:  e-mixture (framework), function run_trials_mixture()
  Status:       Done. Exact closed form instead of Monte Carlo (section 2A).

Prior from a Bayesian borrowing model fitted to all baskets, updated before every analysis: EX/NEX, posterior restricted to (p0, 1]
  In the code:  One-way EXNEX on a (μ, τ) grid, weight π = 0.5, restricted posterior moments
                computed exactly
  Status:       Done for one dimension

Baskets on a two-way grid (tumour type × biomarker group) with additive row and column effects
  In the code:  Not yet: the simulated trials are a list of cohorts
  Status:       Open. With one row the framework's model is the one-way model used here. The two-way
                model needs a 2-D basket layout and an MCMC/INLA or approximate engine (section 7a).

Platform trials: baskets open over time
  In the code:  Umbrella set-up with staggered entry (framework e-value included)
  Status:       Done

Motivation: Bayesian borrowing has no frequentist error guarantee
  In the code:  Bayes mEXNEX, calibrated under the global null, reaches type I error 0.15 in mixed
                scenarios
  Status:       Shown (section 10)

Monte Carlo block factor
  In the code:  Closed form; Monte Carlo rules checked exactly
  Status:       Done. Two additions for the text (section 7).

The DRUP study as concrete example
  In the code:  Not in the code (text of framework.tex)
  Status:       Open

Double check: was the framework primary? Not fully. What changed
----------------------------------------------------------------

- Settings and README: both presented the design-optimal e-values as the method and the framework as
  a compliance appendix. Now Setting.R and this README start from framework.tex.
- Design order: the framework designs came last (positions 18–19), so they were last in every legend
  and had left-over colours. Now they come first, with the first colour and solid lines.
- Flexibility: optional continuation and faster inclusion, the central benefit claimed in
  framework.tex, existed only for the design-optimal e-values. The framework e-value now has them
  too.
- Comparisons: designs were compared only with Simon. A new table, <setup>_design_effects.csv, gives
  every design relative to the framework e-value, plus the effect of borrowing and of flexible
  inclusion.
- Still open: the two-way grid model and the DRUP example (table above).

2. Designs (section 8 of Setting.R)
===================================

A. The framework e-value (primary)
----------------------------------

e-mixture (framework)
  How it works:  At every analysis the prior for the next block is frozen: the posterior of p_c
                 given all platform data (EXNEX with π = 0.5, own data included), restricted to (p0,
                 1]. Block factor = likelihood ratio of the new block mixed over that prior; e-value
                 = product of block factors; reject at the first analysis with e-value ≥ 1/α; stop
                 when 1/α can no longer be reached.
  Set-up:        both

e-mixture (framework, flexible)
  How it works:  The same with optional continuation and faster inclusion: an extension beyond n_max
                 and acceleration (section 8).
  Set-up:        both

e-mixture (no borrowing)
  How it works:  The same with the NEX prior only (own data only): the effect of borrowing.
  Set-up:        basket

How it is computed:
- Prior. The engine computes the mean and variance of the restricted posterior exactly on the grid.
  The prior for the next block is the scaled Beta distribution on (p0, 1] with that mean and
  variance.
- Validity. Like the Monte Carlo average of framework.tex, this replaces the exact restricted
  posterior by another frozen prior on (p0, 1], so validity is unchanged.
- Block factor. In closed form: e = Σ_j choose(y, j) ((1 − p0)/p0)^j B(a + j, b + n − y) / B(a, b).
- Extension decision. The e-value of the framework has no horizon, so a cohort without a decision at
  n_max can be extended there.
  - It is extended by n_ext = n_max/2 patients if the predictive probability that the extension
    brings the e-value to 1/α is ≥ 0.30.
  - That probability is exact: under the prior, P(Y = y) = dbinom(y, n_ext, p0) × factor(y).
- Acceleration. While its e-value is ≥ √(1/α) ≈ 4.5, a cohort includes 2 patients per time unit
  instead of 1.

B. Design-optimal e-values (alternatives within the framework)
--------------------------------------------------------------

- What they are. They bet b ∈ [0, 1] per patient, from the dynamic program of
  design_optimal_evalues.R, for a point alternative θ1.
- Why they fit the framework. A bet b is the likelihood ratio at the point-mass prior p0 + b(1 −
  p0), which lies in (p0, 1].
- Where they differ. Inside a block their bets also depend on earlier outcomes of that block, which
  goes beyond the framework's block factor but not beyond its proof (section 7c).
- Borrowing. θ1 = π · s + (1 − π) · p1, with s from the other cohorts (leave-one-out EX predictive)
  and π from mEXNEX.

e-GROW, e-Pmax, e-based
  How it works:  θ1 = p1 (no borrowing). e-based: e-value-based design with futility stopping,
                 target power = min(0.80, 0.95 × P-max power at θ1). e-GROW = the framework's factor
                 with a point-mass prior at θ1.
  Set-up:        GROW, e-based: both; e-Pmax: basket

e-based + mEXNEX
  How it works:  θ1 re-chosen at every analysis; θ1 = p1 as long as the other cohorts have no
                 patients.
  Set-up:        basket

e-based + mEXNEX (c_e 0.2)
  How it works:  The same with c_e = 0.20: more borrowing.
  Set-up:        basket

e-based + mEXNEX (25% q.)
  How it works:  s = the 25% quantile instead of the mean: a conservative θ1.
  Set-up:        basket

e-based + delayed
  How it works:  Borrowing only once a cohort has ≥ 10 patients.
  Set-up:        both

e-based + EXNEX at entry
  How it works:  π = 0.5; θ1 chosen once, when a cohort opens.
  Set-up:        umbrella

e-based (accelerated / extended / flexible)
  How it works:  e-based with acceleration, extension (decided at about 3/4 of n_max, re-optimised
                 design), or both (section 8).
  Set-up:        basket (flexible: both)

e-based + mEXNEX (flexible), e-based + delayed (flexible)
  How it works:  With acceleration and extension.
  Set-up:        basket / umbrella

C. Comparators outside the framework
------------------------------------

Bayes Umbrella, Bayes mEXNEX
  How it works:  Reject if P(p > p0 | data) > Δ, with Δ calibrated under the global null; futility
                 stop if P(p > (p0 + p1)/2 | data) < 0.10 from 10 patients on. Fixed inclusion.
  Set-up:        basket (umbrella: Bayes Umbrella)

Simon two-stage, exact binomial test
  How it works:  Exact references with n = n_max. In the simon series this is Simon's optimal design
                 itself.
  Set-up:        both

Figures:
- Line type = family: framework e-values solid, design-optimal e-values dashed, Bayesian comparators
  dotted, Simon black dotted.
- Flexible variants have the colour of their base design and line type dot-dash (accelerated), two-
  dash (extended) or long-dash (both). In the trade-off figure, shape shows the inclusion.

EXNEX weights π (Setting.R, section 7):

none
  π_k:  0
  use:  no borrowing

full
  π_k:  1
  use:  BHM

fixed
  π_k:  constant
  use:  EXNEX (the framework e-value: 0.5)

hellinger_mean
  π_k:  mean over other cohorts of 1 − Hellinger distance between Beta(y+1, n−y+1) posteriors; 0 if
        the observed rate differs > c_e from all others
  use:  mEXNEX

hellinger_max
  π_k:  the same, maximum instead of mean
  use:  mEXNEX_min

delayed
  π_k:  0 until n_k ≥ n_burnin, then the base rule
  use:  no initial borrowing

3. What is evaluated
====================

Every results table (one row per design × schedule × scenario) has the columns below.

The new statistics (marked new) come from outputs the simulation already produces: per cohort the
decision, n, time, e-value and the number of responses at the decision. They add no simulation and
no noticeable computing time.

Error control (the framework's criterion)
  Columns:  typeI (mean over null cohorts), typeI_max (worst null cohort), FWER

Power and sample size
  Columns:  power; ESS_H0, ESS_H1 (expected n of null and active cohorts); new: n_q50_*, n_q90_*
            (median and 90th percentile of n)

Sample size to reject
  Columns:  n_reject_H1, early_reject_H1 (P(rejected before n_max)), n_early_reject_H1

Calendar time (from the cohort's entry; 1 = one patient at normal speed)
  Columns:  time_H0, time_H1, time_reject_H1, duration (until every cohort has a decision)

Optional continuation
  Columns:  extend_H0, extend_H1 (P(extended)), reject_ext_H1 (P(rejected after n_max))

Borrowing
  Columns:  theta1_H0, theta1_H1 (mean alternative bet on; framework e-value: mean of the frozen
            prior)

Decisions over the whole trial (new)
  Columns:  true_rejections, false_rejections (per trial), FDR (false discovery rate), power_any (≥
            1 active cohort rejected), power_all (all active rejected), all_correct (all null kept
            and all active rejected), trial_n (expected total n), trial_n_q90, n_per_discovery
            (patients per true rejection)

Estimation (new)
  Columns:  bias_H0, bias_H1, bias_reject_H1 (bias of the naive estimate y/n at the decision; after
            early rejection it is biased upwards), rmse_H1

Strength of evidence (new)
  Columns:  log10E_H0, log10E_H1 (median log10 e-value at the decision; e-value designs)

Effects of each design element (new, <setup>_design_effects.csv). These are differences in the
results above, for the same schedule and scenario:
- borrowing: the design minus the same design without borrowing (framework vs no borrowing, e-based
  + mEXNEX vs e-based, Bayes mEXNEX vs Bayes Umbrella);
- flexible inclusion: the flexible variant minus its base design;
- vs primary design: every design minus e-mixture (framework).

4. The runs (section 2 of Setting.R)
====================================

grid
  Hypotheses:         p0 = 0.10, p1 = 0.30
  Cohorts K:          3, 5, 10
  Cohort size n_max:  20, 30, 40
  Runs:               9

simon
  Hypotheses:         the 8 hypotheses of Table 1 of Simon (1989): p0 = 0.05, 0.10, …, 0.70 with p1
                      = p0 + 0.20
  Cohorts K:          5
  Cohort size n_max:  n of Simon's optimal design: 17, 29, 43, 46, 46, 43, 43, 27
  Runs:               8

Simon's Table 2 (p1 − p0 = 0.15, n up to 84) is a commented line in Setting.R.

What changes per run:
- Hypotheses (settings$for_hypothesis):
  - θ1 without borrowing = p1; with borrowing θ1 ranges from p0 + 0.05 to min(0.95, p1 + 0.30).
  - Bayesian futility threshold (p0 + p1)/2.
  - EX prior mean logit(p0); NEX prior centred at (p0 + p1)/2.
- Schedules: fully sequential; blocked every round(n_max/3); two-stage (round(n_max/2) + the rest).
  Analyses are in calendar time.
- Scenarios (d = p1 − p0):
  - 0, 1, ⌊K/2⌋, K − 1 and K cohorts active (p = p1, the others p0).
  - mixed: rates from p0 − d/4 to p1 + d/2.
  - strong: all at p1 + d/2.
- Umbrella entry: two cohorts open together, the next two round(n_max/3) time units later, and so
  on.
- Extension size: n_ext = round(n_max/2).

5. Files and how to run
=======================

Setting.R
  Content:  The only file to run. Top: all settings (table of runs, rules per hypothesis, flexible
            inclusion, designs and set-ups, with the framework designs first). Bottom: the run.

additional_functions.R
  Content:  All functions, nothing is run: EXNEX model (including the restricted posterior),
            π-rules, the framework e-value, design-optimal e-values and conditional power, trial
            simulation in calendar time, summaries and design effects, Simon designs, figures,
            running one row of the table.

design_optimal_evalues.R
  Content:  Your single-cohort e-value library, unchanged (used by family B).

design_banks/
  Content:  Cache of the design-optimal designs (350 files, 85 MB) for all 17 runs, so the run only
            simulates. The framework e-value needs none.

Result/
  Content:  Created by the run.

checks/
  Content:  The individual and exact checks of section 9, rerunnable (checks/README_checks.md).

PROJECT_HISTORY.md
  Content:  How the work developed, decisions, working agreements and open items: start here when
            resuming.

    Rscript Setting.R

In a slurm file:

    #SBATCH --cpus-per-task=8          # the run uses exactly the cores slurm allocates
    #SBATCH --mem=8G
    #SBATCH --time=09:00:00            # ~4-4.5 h needed on 8 cores (see below)
    export N_SIM=1000                  # optional: trials per scenario (default 1000)
    Rscript /path/to/Setting.R

- Quick check: N_SIM=1 runs everything in a few minutes.
- Fewer runs: delete rows from settings$runs.
- Faster: remove secondary designs from settings$setups (for example the c_e 0.2 and 25% q.
  variants).

Run time (estimates, 1,000 trials per scenario):
- Basis: K = 5, n_max = 30 took 10.7 min on 2 cores with the earlier design list (measured). The
  other runs are extrapolated with K² × n_max.
- This version: 16 basket and 9 umbrella designs. The flexible framework design adds about 10%. The
  new statistics add nothing noticeable.

                       2 cores  8 cores
---------------------  -------  --------
grid series (9 runs)   ~4.5 h   ~2.5–3 h
simon series (8 runs)  ~3 h     ~1.5 h
Total                  ~7–8 h   ~4–4.5 h

Memory: about 2.5 GB on 2 cores; use --mem=4G for 2 cores and --mem=8G for 8–10 cores.

6. Output (Result/)
===================

Per run, in Result/<series>_p<p0>-<p1>_K<K>_n<n_max>/:

- basket_results.csv, umbrella_results.csv: all statistics of section 3.
- *_design_effects.csv: effects of borrowing, of flexible inclusion, and of every design vs the
  framework e-value.
- *_curves.csv: cumulative rejection and futility by patient, and rejection by calendar time.
- *_fig1_cumulative: Figure 1 of the e-value paper (rejection and futility by patient, per schedule;
  Simon in the two-stage column).
- *_fig2_scenarios: by scenario: type I error, power, early rejection, expected n, n and time at
  rejection.
- *_fig3_time: rejection by calendar time (type I error, power, power in 'strong'); shows
  acceleration.
- *_fig4_tradeoff (new): what each design buys. Power against the expected n of a null cohort and
  against the time to rejection, per schedule, with Simon. Colour = design, shape = inclusion
  (fixed, accelerated, extended, both).
- *_fig5_decisions (new): by scenario: true and false rejections per trial, FDR, P(all decisions
  correct), total trial size, bias of y/n in rejected active cohorts, median log10 e-value at the
  decision.
- references.csv: Simon's design and the exact binomial test.

All runs together, in Result/:
- all_settings_results.csv: all rows, with series, p0, p1, K, n_max.
- overview_<series>_<setup>: worst type I error, power, early rejection, n and time at rejection,
  expected n; grid series by K and n_max, simon series by hypothesis and schedule; Simon as black
  dotted line.
- all_results.rds: everything in one R object.

7. Accordance with framework.tex in detail
==========================================

Two-way grid c = (j, l)
  Implementation:  A list of K cohorts
  In accordance?:  Partly (a)

Patients ordered by outcome time; basket known before the outcome (G_t)
  Implementation:  Outcomes immediate: inclusion order = outcome order
  In accordance?:  Yes (no outcome delay)

H0,c conditional on the history: drift allowed, nothing about other baskets
  Implementation:  Constant rates in the simulation; exact checks with drift and adversarial other
                   baskets
  In accordance?:  Yes

A separate α per basket
  Implementation:  α per cohort; FWER and FDR only reported
  In accordance?:  Yes

Analysis times are data-dependent stopping times
  Implementation:  Planned calendar times; acceleration and extension make them data-dependent,
                   using only e-values, counts and θ1
  In accordance?:  Yes

Prior on (p0, 1], determined by F_τ, frozen, may use own data
  Implementation:  Framework e-values: restricted posterior incl. own data. Design-optimal: point
                   mass at θ1 ∈ [p0 + 0.05, 0.95] from the other cohorts
  In accordance?:  Yes (b)

Block factor = mixture likelihood ratio over the frozen prior
  Implementation:  Framework e-values: yes. e-GROW: yes (point mass). Sequential e-Pmax / e-based:
                   yes (point mass per patient). Blocked e-Pmax / e-based: bets inside a block
                   depend on that block
  In accordance?:  Partly, still valid (c)

E = product; reject at the first analysis with E ≥ 1/α
  Implementation:  Framework e-values: yes. Design-optimal: betting stops once E ≥ 1/α inside a
                   block
  In accordance?:  Partly, still valid (d)

No futility rule
  Implementation:  Hopeless-zone stop (all designs), e-based futility rule
  In accordance?:  Yes (stopping is always allowed)

No maximum sample size
  Implementation:  Framework e-value: none needed. Design-optimal: n_max, extension by re-
                   optimisation
  In accordance?:  Yes

EX/NEX: θ_c = μ + a_j + b_l + ε_c, NEX N(m_c, v_c), half-normal SDs
  Implementation:  One-way EXNEX: θ ~ N(μ, τ²), τ ~ HN(1), μ ~ N(logit p0, 10²) on [−8, 4]; π = 0.5
                   (framework) or mEXNEX (design-optimal)
  In accordance?:  Partly (a)

Prior = posterior restricted to (p0, 1]
  Implementation:  Mean and variance exact, scaled Beta with those moments
  In accordance?:  Yes (e)

Block factor by Monte Carlo
  Implementation:  Closed form
  In accordance?:  Same purpose (f)

(a) One dimension instead of two.
- With one row, a_1 is absorbed in μ and b_l + ε_c ~ N(0, σ_b² + σ_ε²): the framework's model is
  then the one-way model used here, with τ² = σ_b² + σ_ε². Only the prior on τ differs.
- The two-way model needs a two-dimensional basket layout and an MCMC, INLA or approximate engine
  (the quadrature engine integrates over (μ, τ) only).
- Validity does not depend on the model.

(b) Own data.
- The framework e-value uses the full posterior including the basket's own data, as framework.tex
  proposes.
- The design-optimal e-values take θ1 from the other cohorts only (allowed: 'may use').

(c) Design-optimal bets inside a block.
- With a frozen prior the implied bet rises after a response and falls after a non-response.
- Blocked and two-stage e-based designs lower their bet after a response in 54–55% of the states
  inside a block, so they are not of the framework's form.
- They are valid, because the framework's proof only needs every bet fixed before its outcome.
- Suggested remark for framework.tex: below.

(d) Betting stops at 1/α.
- The design-optimal e-values stop betting once E ≥ 1/α inside a block, so E stays ≥ 1/α until the
  analysis.
- That is a bet fixed in advance: valid, and no rejection is lost.

(e) The restricted posterior. Checked:
- the closed form equals numerical integration (largest relative difference 5·10⁻¹²);
- the engine's restricted moments equal a direct grid computation (< 10⁻¹⁵);
- block factors are within 0.5% of those with the exact restricted posterior for ordinary blocks,
  and within −1.4% to +8% for extreme blocks.

(f) Monte Carlo needs two additions.
- The draws must be frozen. Draws made after seeing the outcome raise the type I error from 0.022 to
  0.051; frozen draws give 0.021–0.022.
- The formula is undefined when no draw exceeds p0 (H' = 0). That happens in 8–18% of the analyses
  of null cohorts with 0–1 responses of 30 (EXNEX, 1,000 draws). Set the factor to 1 then.

Checks of the framework's claims (exact, no simulation):

1. Restriction to (p0, 1] is essential for H0: p ≤ p0. Probability of E ≥ 20 by n = 150, uniform-
   prior mixture, fully sequential:

   | true p | restricted | not restricted | |---|---|---| | 0.10 | 0.027 | 0.028 | | 0.05 | 0.004 |
0.123 | | 0.02 | 0.0005 | 0.745 |

2. Freezing is essential. A prior that includes the patient's own outcome gives type I error 0.058
   (frozen: 0.022).
3. Whatever the other basket does.
   - Set-up: two cohorts, all 2^15 sequences of cohort 1, and five paths for cohort 2 (all
     responses, none, responses then none, a copy and the opposite of cohort 1).
   - Designs: both framework e-values with and without borrowing, the flexible framework e-value,
     e-based + mEXNEX with and without flexibility; sequential and blocked.
   - Type I error of cohort 1 ≤ α everywhere (largest 0.048).
4. Drifting and history-dependent response rates ≤ p0 (four patterns), framework and design-optimal
   e-values, fixed and flexible: type I error ≤ α everywhere (largest 0.048).
5. e-GROW is exactly the framework's factor with a point-mass prior at θ1: its bet equals (θ1 −
   p0)/(1 − p0) in all 443,778 state-patient entries (difference < 5·10⁻¹³).

Suggested additions to framework.tex:
- General form. After Proposition 1: > Remark. More generally, let $e_{c,k} = \prod_{i \in B_{c,k}}
  \{1 + b_i (Y_i/p_0 - 1)\}$ with bets > $b_i \in [0,1]$ that are $\G_i$-measurable. Then the
  conclusion of Proposition 1 still holds, since each factor > has conditional mean $1 + b_i\{\P(Y_i
  = 1 \mid \G_i)/p_0 - 1\} \le 1$. The mixture factor > (\ref{eq:factor}) is the special case $b_i =
  (q_i - p_0)/(1-p_0)$, with $q_i$ the mean of $\pi_{c,k-1}$ updated > with the earlier outcomes of
  block $k$. Design-optimal e-values are another special case.
- Monte Carlo: draws generated from $\F_{\tau_{c,k-1}}$ only; if $H' = 0$, set $e_{c,k} = 1$.
- Ties: fix the order of outcomes that arrive at the same time (for example by inclusion).
- Why (p0, 1]: the restriction makes E an e-value for the composite null p ≤ p0 (check 1).
- Optional continuation: extending a basket and speeding up inclusion are optional continuation in
  the framework's sense. The framework e-value needs no maximum sample size, so a basket can be
  extended at the planned end (section 8).
- The orange 'True?' (rows tumour types, columns biomarker groups) is still open.

8. Optional continuation in detail: extension beyond n_max and faster inclusion
===============================================================================

Validity. Every e-value here is a running product M_t = ∏ (1 + b_i (Y_i/p0 − 1)) with b_i ∈ [0, 1]
fixed before Y_i is known.
- Under p ≤ p0 each factor has conditional mean ≤ 1, so P(M_t ≥ 1/α for some t) ≤ α (Ville).
- This holds for any rule that decides, from data already seen, how many patients a basket includes,
  how fast, when it is analysed, whether it is extended, and which prior it uses.

Four rules keep it so:
1. Fixed order: outcomes enter in a fixed order; each bet uses only earlier outcomes.
2. Bets in [0, 1].
3. Reject at the first analysis with M ≥ 1/α; the rejection is final.
4. Never restart: extending means continuing the same product.

Without e-values it fails. A naive extension of the exact binomial test inflates the type I error:
- test at n_max;
- if within two responses of significance, add n_ext = n_max/2 patients;
- test at the critical value for n_max + n_ext.

p0 / p1      n_max + n_ext  type I without extension  type I with extension
-----------  -------------  ------------------------  ---------------------
0.05 / 0.25  17 + 8         0.0088                    0.0318
0.10 / 0.30  29 + 14        0.0216                    0.0325
0.20 / 0.40  43 + 22        0.0362                    0.0544
0.30 / 0.50  46 + 23        0.0363                    0.0496
0.40 / 0.60  46 + 23        0.0343                    0.0502
0.50 / 0.70  43 + 22        0.0330                    0.0477
0.60 / 0.80  43 + 22        0.0354                    0.0539
0.70 / 0.90  27 + 14        0.0202                    0.0463

The framework e-value (family A): extension at n_max.
- It has no horizon and does not bet all-or-nothing near n_max. So the decision can be taken at
  n_max: extend by n_ext if the exact predictive probability that the extension reaches 1/α is ≥
  cp_min = 0.30.
- The hopeless zone is measured against n_max + n_ext until the decision.

Design-optimal e-values (family B): extension needs re-optimisation.
- The paper's designs are optimal for a fixed n_max (Sections 1, 2.1). Their power guarantee
  survives other analysis times only while n_max stays the same (Sections 3.3, 5). Past n_max, power
  is not analysed; the remedy for a change is to re-solve the dynamic program (Section 3.4).
- Near n_max their bets become all-or-nothing, so a cohort still undecided at n_max has an e-value
  near 0. With the decision at n_max, no cohort was ever extended (exact, all 2^18 sequences).
- Fix:
  - decide at the last planned analysis ≤ 3/4 · n_max;
  - extend if the planned design will probably fall short (conditional power < 1 − β) and the design
    re-optimised for the remaining patients + n_ext (P-max, optimal from every e-value) has
    conditional power ≥ 0.30;
  - an extended cohort bets with that design.
  - Both conditional powers are exact (backward recursion, equal to the library's dynamic program).

Faster inclusion (both families).
- While its e-value is ≥ √(1/α), a cohort includes 2 patients per time unit. The bets per patient do
  not change; analyses stay in calendar time, so blocks become larger.
- For a design-optimal e-value with the same bets, acceleration never loses a rejection and never
  rejects later.

Settings: settings$flex (section 4b of Setting.R): speed_up = 2, m_speed = √(1/α), ext_frac = 0.5,
ext_at = 0.75, cp_min = 0.30.

Protocol checklist:
1. Pre-specify p0, α, the e-value (prior model, analyses), the futility, extension and acceleration
   rules. This is needed for power, sample size and duration, not for the type I error.
2. At each analysis: update the e-value in the fixed order, report it, reject if ≥ 1/α, apply the
   rules.
3. Unplanned changes keep the type I error as long as rules 1–4 hold. For good power, update the
   prior (framework e-value) or re-optimise the bets (design-optimal) from the current e-value.
4. Report per basket: e-value, n, calendar time of the decision, accelerated or extended.

9. Checks of this version (individual checks, no simulation)
============================================================

- Framework engine unchanged where it should be. The new calendar-time framework engine reproduces
  the previous one exactly without flexibility: 126 combinations (two designs, two runs, both set-
  ups, all schedules and scenarios), with identical decisions, n, e-values and times.
- All design-optimal and Bayesian designs unchanged: 840 combinations identical to the original
  engine.
- Extension rule. The prior-predictive probabilities of an extension block sum to 1 (largest error
  8·10⁻¹⁵).
- The flexible framework e-value keeps type I error ≤ α, exact:
  - against four adversarial other-cohort paths: largest 0.0076 (fixed framework: 0.017);
  - under four drift patterns: largest 0.0065.
- Rule checks on hand-made data for every design, schedule and scenario of three runs:
  - every cohort gets a decision;
  - n ≤ n_max, or ≤ n_max + n_ext for extended cohorts only;
  - time = n at normal speed and ≤ n when accelerated;
  - framework e-values reject exactly when E ≥ 1/α, with the prior mean always above p0;
  - an accelerated design-optimal cohort is rejected whenever its fixed version is, and never later.
- New statistics consistent:
  - FDR in [0, 1];
  - P(all active rejected) ≤ power ≤ P(any active rejected);
  - false rejections per trial ≥ FWER;
  - quantiles ordered;
  - e-value at the decision present exactly for the e-value designs;
  - design effects equal the differences of the table entries.
- The whole pipeline ran for three runs with all 25 design × set-up combinations, all tables and
  figures (fig1–fig5), and the figures were looked at.
- Two bugs found and fixed:
  - the flexible framework design looked for design-optimal extension designs;
  - the decisions figure looked for a type I row it does not have.
- Still valid from earlier versions:
  - the exact checks of the framework's claims (section 7);
  - the flexible design-optimal e-values (2^18 sequences, independent re-implementation);
  - Simon's 48 designs;
  - the 1-trial run of all 17 runs.

10. Results so far (exact or from earlier runs; the full run is next)
=====================================================================

The framework e-value vs the design-optimal e-value with borrowing (two cohorts, n_max = 10 + n_ext
= 5, all 2^15 sequences of cohort 1, four or five paths of cohort 2; ranges over the paths):

e-mixture (framework), sequential
  type I error:             0.014–0.017
  power at p = 0.30:        0.24–0.26
  P(extended) at p = 0.30:  –

e-mixture (framework), blocked
  type I error:             0.003
  power at p = 0.30:        0.17
  P(extended) at p = 0.30:  –

e-mixture (framework, flexible), sequential
  type I error:             0.005–0.008
  power at p = 0.30:        0.27–0.28
  P(extended) at p = 0.30:  0.15–0.17

e-mixture (framework, flexible), blocked
  type I error:             0.004
  power at p = 0.30:        0.25–0.26
  P(extended) at p = 0.30:  0.19

e-based + mEXNEX, sequential / blocked
  type I error:             0.041–0.048
  power at p = 0.30:        0.44–0.51
  P(extended) at p = 0.30:  –

e-based + mEXNEX (flexible)
  type I error:             0.041–0.047
  power at p = 0.30:        0.46–0.57
  P(extended) at p = 0.30:  –

- At this small size the framework e-value is conservative. It uses a fraction of its α and has
  about half the power of the design-optimal e-value. It spreads its bets over the whole prior
  instead of aiming at 1/α by n_max.
- Optional continuation helps the framework e-value most where it is weakest: blocked power 0.17 →
  0.25.
- The full run shows how this develops at realistic sizes.

Extension of the design-optimal e-value (one cohort, exact, blocked, Simon's hypotheses and n):
power 0.80 → 0.82–0.85 with unchanged type I error, at a cost of 1–5 patients per cohort.

p0 / p1      n_max + n_ext (decision at)  type I: fixed, extended  power: fixed, extended  E[n] under H0: fixed, extended  E[n] under H1: fixed, extended  P(extended): H0, H1
-----------  ---------------------------  -----------------------  ----------------------  ------------------------------  ------------------------------  -------------------
0.05 / 0.25  17 + 8 (12)                  0.048, 0.047             0.78, 0.82              12.3, 13.0                      11.1, 11.9                      0.09, 0.20
0.10 / 0.30  29 + 14 (20)                 0.048, 0.048             0.80, 0.83              15.5, 16.4                      16.7, 17.9                      0.07, 0.11
0.20 / 0.40  43 + 22 (28)                 0.048, 0.048             0.80, 0.83              19.4, 21.4                      23.6, 25.6                      0.09, 0.12
0.30 / 0.50  46 + 23 (30)                 0.048, 0.048             0.80, 0.84              22.8, 24.9                      27.5, 30.4                      0.10, 0.15
0.40 / 0.60  46 + 23 (30)                 0.048, 0.047             0.80, 0.84              24.4, 26.7                      29.3, 32.7                      0.12, 0.17
0.50 / 0.70  43 + 22 (28)                 0.047, 0.047             0.80, 0.85              24.0, 26.5                      28.9, 32.9                      0.14, 0.20
0.60 / 0.80  43 + 22 (28)                 0.048, 0.048             0.80, 0.83              19.5, 20.7                      25.5, 27.5                      0.06, 0.10
0.70 / 0.90  27 + 14 (18)                 0.048, 0.048             0.79, 0.82              13.5, 14.5                      20.5, 22.3                      0.09, 0.13

Acceleration of the design-optimal e-value (one cohort, exact, p0 = 0.30, n_max = 12): same power
and type I error, time to rejection 14–21% shorter, a few tenths of a patient more.

Single cohort without borrowing, fixed inclusion, at Simon's hypotheses (exact):
- e-based matches Simon's power with Simon's n.
- Analysed after every patient it needs fewer patients under H0 and rejects at about half of Simon's
  n.
- e-Pmax has more power (0.82–0.89) but rejects close to n_max; e-GROW has less power (0.62–0.78).

0.05 / 0.25
  n:                                                                                17
  Simon: power, E[n] under H0:                                                      0.81, 12.0
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.78, 10.5,
                                                                                    0.78, 8.1
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.78, 12.3,
                                                                                    0.64, 10.4
  e-Pmax, blocked: power, mean n at rejection:                                      0.82, 12.4

0.10 / 0.30
  n:                                                                                29
  Simon: power, E[n] under H0:                                                      0.81, 15.0
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 13.2,
                                                                                    0.79, 12.2
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 15.5,
                                                                                    0.73, 16.3
  e-Pmax, blocked: power, mean n at rejection:                                      0.87, 23.0

0.20 / 0.40
  n:                                                                                43
  Simon: power, E[n] under H0:                                                      0.80, 20.6
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 16.3,
                                                                                    0.80, 18.7
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 19.4,
                                                                                    0.79, 23.7
  e-Pmax, blocked: power, mean n at rejection:                                      0.89, 37.4

0.30 / 0.50
  n:                                                                                46
  Simon: power, E[n] under H0:                                                      0.80, 23.6
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 19.0,
                                                                                    0.79, 22.8
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 22.8,
                                                                                    0.79, 27.3
  e-Pmax, blocked: power, mean n at rejection:                                      0.87, 42.1

0.40 / 0.60
  n:                                                                                46
  Simon: power, E[n] under H0:                                                      0.80, 24.5
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 20.0,
                                                                                    0.78, 24.9
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 24.4,
                                                                                    0.78, 28.8
  e-Pmax, blocked: power, mean n at rejection:                                      0.85, 43.4

0.50 / 0.70
  n:                                                                                43
  Simon: power, E[n] under H0:                                                      0.80, 23.5
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 19.3,
                                                                                    0.77, 24.8
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 24.0,
                                                                                    0.76, 28.1
  e-Pmax, blocked: power, mean n at rejection:                                      0.84, 41.4

0.60 / 0.80
  n:                                                                                43
  Simon: power, E[n] under H0:                                                      0.80, 20.5
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.80, 14.4,
                                                                                    0.78, 22.5
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.80, 19.5,
                                                                                    0.78, 26.0
  e-Pmax, blocked: power, mean n at rejection:                                      0.89, 41.8

0.70 / 0.90
  n:                                                                                27
  Simon: power, E[n] under H0:                                                      0.80, 14.8
  e-based, sequential: power, E[n] under H0, P(reject early), mean n at rejection:  0.79, 11.2,
                                                                                    0.76, 18.5
  e-based, blocked: power, E[n] under H0, P(reject early), mean n at rejection:     0.79, 13.5,
                                                                                    0.46, 21.8
  e-Pmax, blocked: power, mean n at rejection:                                      0.83, 27.0

The motivation of framework.tex (K = 5, n_max = 30, blocked, 10,000 trials, earlier version):

Design            Type I, all null  Worst type I, mixed scenarios  Power, all at 0.30  Expected n, null cohort
----------------  ----------------  -----------------------------  ------------------  -----------------------
Bayes Umbrella    0.035             0.037                          0.837               19.2
Bayes mEXNEX      0.049             0.152                          0.923               15.4
e-based           0.048             0.052                          0.810               15.6
e-based + mEXNEX  0.049             0.052                          0.806               16.0

- Bayesian borrowing calibrated under the global null reaches type I error 0.15 for null baskets
  next to active ones.
- The e-value designs stay at α for every basket in every scenario.
