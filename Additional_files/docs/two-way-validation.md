# Two-way EXNEX: validation of the posterior calculations

Run on 2026-11-26 with `Additional_files/validation/twoway_validation.R` (R 4.3.3, JAGS 4, 2 cores; about 30 min).
* The computations ran at commit `845fee1`. Until the reviewed PR head `a8a915f`, the engine code changed only in its
  figure functions. The code was later moved into `additional_functions.R` with bit-identical results (2026-11-26
  restructure). These are archived results from that engine version, not a rerun of the current branch.
* The summary step first failed: the `aggregate()` call was missing `FUN =` and `data =`. After the fix, the
  summaries were recomputed from the saved results with `SUMMARY_ONLY=1`.
* Raw results, per method, state and basket: `Additional_files/validation/results/twoway_validation.csv`.

These are checks of **posterior calculations for fixed data states**, not trial simulations. Accuracy here matters
for power and borrowing only. The type I error guarantee requires a predictable, frozen probability measure on
(p0, 1], including normalized weights as actually used by the simulator. It does not require an exact posterior.
`Additional_files/tests/test_validity.R` enumerates rejection probabilities for specified small examples with every engine;
it does not test the simulator's handling of those weights.

## Design

p0 = 0.10, pbar = 0.20, theta cells of width 0.05, n up to 30 per basket.

**Methods.**

| | Method | Integration over (mu, a, b, sigmas) |
|---|---|---|
| **production** | the sampler of `additional_functions.R` with default settings: 2 chains x (500 burn-in + 1000 kept); 5 seeds | Monte Carlo (Metropolis-within-Gibbs) |
| quadrature | `twoway_exact()`: deterministic numerical reference on the 0.1 lattice | grid over the column levels (mu in closed form) and 1-D per row effect |
| one-way equivalent | L x 1 Gaussian reduction: mu' = mu + b_1 and tau² = sigma_a² + sigma_e²; a numerical cross-check after conditional theta truncation | deterministic grid over (mu', tau) |
| JAGS | the model written out independently: continuous theta and sigmas, 4 x 25,000 draws | MCMC |
| long run | the production sampler with 4 chains x 20,000 draws | Monte Carlo |

The quadrature and production sampler share the likelihood tables, so their agreement checks integration over
hyperparameters without establishing convergence of the underlying theta and m grids. The one-way Gaussian
reduction is exact before conditional theta truncation; integrating already normalized truncated conditionals
differs slightly from truncating the reduced Gaussian. JAGS uses continuous theta and sigmas and has its own
Monte Carlo error. These comparisons therefore assess several numerical approximations, not exact agreement
with a common continuous posterior oracle.

**Sections.**

| Section | Layout and priors | Reference |
|---|---|---|
| A | 3 x 1, full hyperpriors | quadrature (also compared: one-way equivalent, JAGS, production) |
| B | 3 x 3 with sigma_a = sigma_b = sigma_e = 0.45 fixed | quadrature (also: JAGS, production) |
| C | 3 x 3, full hyperpriors | JAGS (also: long run, production) |

**Data states** (y/n per basket, rows = tumour types). They cover the cases asked for: empty baskets, unbalanced
recruitment, strongly heterogeneous effects and interactions.
* **3 x 3:**
  * homogeneous null (3/30 everywhere);
  * one active basket (10/30);
  * active row; active column; diagonal interaction (9/30 against 3/30);
  * unbalanced with an empty basket, n = 0 to 30;
  * missing cell (3,3) plus an empty basket (2,2);
  * strong heterogeneity (1/30 to 18/30);
  * early (5 per basket).
* **3 x 1:** heterogeneous (3, 9, 14 of 30); homogeneous null; unbalanced (n = 5, 10, 30); with an empty basket.

**Measures.** Compared per basket:
* P(p > p0 | data);
* the restricted mean and sd: the mean and sd of the prior handed to the e-process;
* P(EX | data);
* the block factor of the next 5 patients with 0, 2 or 5 responses under the restricted prior (f0, f2, f5).

## Results

**Largest absolute difference from the section's reference**, over all states, baskets and seeds. For the factors,
the difference is |log ratio|; 0.1 is about 10%.

| Section | Method | restricted mean | restricted sd | P(p > p0) | P(EX) | log f0 | log f5 |
|---|---|---|---|---|---|---|---|
| A | one-way equivalent | 0.0002 | 0.0001 | 0.0004 | 0.0005 | 0.0004 | 0.002 |
| A | JAGS | 0.0044 | 0.0025 | 0.009 | 0.013 | 0.009 | 0.096 |
| A | **production** | 0.0094 | 0.0094 | 0.021 | 0.047 | 0.023 | 0.119 |
| B | JAGS | 0.0020 | 0.0019 | 0.006 | 0.010 | 0.011 | 0.083 |
| B | **production** | 0.0067 | 0.0034 | 0.026 | 0.037 | 0.039 | 0.080 |
| C | long run | 0.0026 | 0.0015 | 0.012 | 0.022 | 0.013 | 0.066 |
| C | **production** | 0.0074 | 0.0039 | 0.031 | 0.071 | 0.032 | 0.154 |

**Typical difference of the production sampler**, as the median / 90th percentile over seeds and baskets, by the
basket's own sample size.

| Section | Basket data | Cases | restricted mean | P(EX) | log f0 | log f5 |
|---|---|---|---|---|---|---|
| A | n = 0 | 5 | 0.0060 / 0.0094 | 0 / 0 | 0.010 / 0.019 | 0.055 / 0.088 |
| A | n = 5–10 | 20 | 0.0012 / 0.0037 | 0.007 / 0.021 | 0.009 / 0.019 | 0.017 / 0.071 |
| A | n ≥ 15 | 35 | 0.0006 / 0.0014 | 0.013 / 0.033 | 0.003 / 0.009 | 0.016 / 0.045 |
| B | n = 0 | 10 | 0.0018 / 0.0029 | 0 / 0 | 0.005 / 0.014 | 0.009 / 0.020 |
| B | n = 5–10 | 60 | 0.0014 / 0.0038 | 0.007 / 0.016 | 0.006 / 0.017 | 0.022 / 0.049 |
| B | n ≥ 15 | 330 | 0.0004 / 0.0011 | 0.004 / 0.013 | 0.002 / 0.007 | 0.013 / 0.032 |
| C | n = 0 | 10 | 0.0018 / 0.0064 | 0.003 / 0.005 | 0.013 / 0.025 | 0.019 / 0.030 |
| C | n = 5–10 | 60 | 0.0014 / 0.0034 | 0.005 / 0.013 | 0.008 / 0.019 | 0.027 / 0.066 |
| C | n ≥ 15 | 330 | 0.0007 / 0.0017 | 0.007 / 0.021 | 0.004 / 0.010 | 0.023 / 0.066 |

In section C the reference (JAGS) has Monte Carlo error of its own: the long run and JAGS differ by up to 0.0026 in
the restricted mean.

**Monte Carlo spread of the production sampler.** The largest standard deviation over the 5 seeds, for any state and
basket:

| Section | restricted mean | P(EX) | log f5 |
|---|---|---|---|
| A | 0.0074 | 0.027 | 0.069 |
| B | 0.0049 | 0.022 | 0.042 |
| C | 0.0032 | 0.033 | 0.075 |

**Run time for one 3 x 3 state.**

| Method | Run time |
|---|---|
| production | 3.3 s, of which about 2.5 s is engine set-up; 150 states take 7.5 s, about 0.035 s per state |
| long run | 16 s |
| JAGS | 11 s |
| quadrature with fixed sigmas | 62 s |
| one-way EXNEX | 0.9 s, including engine set-up |

## Conclusions

1. **The deterministic references agree closely in the tested states.** For the 3 x 1 layout, the largest
   differences between the two-way quadrature and one-way cross-check are about 0.000223 in the restricted mean,
   0.000372 in P(p > p0), and 0.000524 in P(EX). The largest absolute log factor ratio is 0.002196
   (about 0.22%). Their different numerical representations and truncation order should be kept in mind.
2. **JAGS agrees with the quadrature** to 0.002–0.004 in the restricted mean. The largest factor differences
   (about 10% for 5/5 responses) occur in baskets with little data, where the factor depends on the far tail of the
   prior.
3. **The production sampler generally agrees with the references within the observed variation.** Its seeds
   scatter on both sides of the reference in the empty-basket example (`twoway_validation.csv`). Five seeds
   do not establish unbiasedness or exclude finite-chain and discretization bias.
   * Baskets with n ≥ 15: typical errors are 0.0004–0.0007 in the restricted mean and about 1–2% in the 5/5
     factor.
   * The largest errors are for empty and 5-patient baskets: up to 0.009 in the restricted mean and 12–15% in the
     5/5 factor.
   * On the scale that decides rejection (log E against log 20 ≈ 3.0), a 10% error in one factor is 0.1.
4. **P(EX) is the least precise output** (Monte Carlo sd up to 0.03). It is a borrowing diagnostic only and does
   not enter the e-process.
5. **A finding about the references.** A quadrature window around the data alone cuts the configurations in which a
   whole column is NEX, which the wide prior of mu (sd 10) supports. That biased P(EX) upwards by up to 0.045. The
   quadrature now covers the whole lattice, with a coarser grid in the tails.

**Default settings.** The production defaults (2 x 1000 draws) are kept for the study. If convergence and the
autocorrelation structure are stable, doubling `keep` would approximately halve Monte Carlo variance, while
increasing the sampling cost. It would not remove discretization or finite-chain bias. Focused draw-count and
grid-refinement checks should precede interpretation of a production power study; see `Additional_files/docs/known-limitations.md`.

**Not yet validated:**
* layouts larger than 3 x 3;
* the hyperprior sensitivity variants;
* states with n > 30 (the engine supports them up to its `n_max`).
* systematic refinement of theta/m/Rao-Blackwell grids and draw counts for the planned study.

---

## Addendum 2026-11-26: the eBBB hyperparameters on the 4 x 3 grid

**Script.** `Additional_files/validation/twoway_ebbb_hyper.R`, about 25 minutes on one core. Results in
`results/twoway_ebbb_hyper.csv` and `_summary.csv`.

**Prior.** The prior of the study (`make_engine("two-way EX/NEX")`), as in eBBB's `exnex.stan`:
* `mu ~ N(logit p0, 2^2)` on `[-10, 6]`;
* `sigma_a, sigma_b ~ HN(1)` on the grid 0.05, …, 3.95;
* `sigma_e ~ HN(0.5)`;
* NEX `N(logit p0, 2^2)`;
* `P(EX) = 0.5`.

**Design.** Six data states on the 4 x 3 layout:
* global null, biomarker 1 active, interaction and all active, each with 20 patients per basket;
* unbalanced, with cell (3,3) missing and an empty basket;
* 5 patients per basket.

The production sampler (`mcmc` of `setting.R`, 5 seeds) was compared with JAGS (4 x 25000 draws, SDs truncated at
the grid ends) and with a long run of the sampler (4 x 20000).

| Comparison | Restricted mean, \|diff\| | 5/5 block factor, relative \|diff\| | P(p > p0), \|diff\| |
|---|---|---|---|
| production vs JAGS | median 0.0012, 95% 0.0046, max 0.0077 | median 3.3%, 95% 9.6%, max 17% | median 0.004, max 0.032 |
| long run vs JAGS | median 0.0006, max 0.0045 | median 1.5%, max 9.5% | max 0.014 |
| production, sd over 5 seeds | median 0.0013 | median 3.6% | — |

**Findings.**
* **Monte Carlo error dominates.** The production sampler's differences from JAGS are of the size of its
  seed-to-seed spread.
* **No systematic bias.** The long run agrees with JAGS to the level reported above for the 3 x 3 layout with the
  earlier hyperparameters. The mean signed difference in P(p > p0) is -0.002, and in P(EX) -0.0007.
* **Largest differences.** They are in the 5-patient state and in baskets with little data, as before.
* **Run time.** One data state took 6–8 s with the production settings, including engine set-up; JAGS took
  13–17 s.
* **Conclusion.** The validation conclusions above carry over to the eBBB hyperparameters on 4 x 3.
