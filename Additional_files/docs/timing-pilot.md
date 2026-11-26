# Timing pilot (2026-11-26)

**Script.** `Additional_files/validation/timing_pilot.R`. Results are in `validation/results/timing_pilot_*.csv`.
The run took 4.8 min for the pilot and 2 min for the calibration, on 2 cores.

**Run approved by the researcher.** 20 trials per cell, with the seeds of `simulate.R`. The pilot results
(rejection rates and so on) are **not** interpreted.

## What was measured

* **A. Pilot.** 84 of the 2 x 154 cells:
  * all 7 models;
  * both recruitment designs;
  * 3 scenarios (global null, all active, strong heterogeneity);
  * both horizons.

  Every fit of the borrowing model was timed, together with its number of data states.
* **B. Two-way fit time.** It was measured against the number of data states S (20, 150, 300, 600; random states
  with n ≤ 30), in one process. The time is linear in S:

  | Layout | Fixed cost per fit | Per data state |
  |---|---|---|
  | balanced, 12 baskets | 1.2 s | 0.043 s |
  | unbalanced, 11 baskets | 0.5 s | 0.041 s |

* **C. Prediction for R = 1000.**
  * **Data states per fit:** the pilot count x 50, assuming the pilot's frequency of distinct states is representative.
    Identical states are fitted once; scaling a finite pilot does not give an upper bound for another run.
  * **Fit time:** two-way fits by the line of B; the other models linearly.
  * **Rest of the simulation:** linear in R.
  * **Scenarios:** all 11 use the largest measured value of the same model, recruitment design and horizon.

## Result

**Predicted CPU hours, all 154 cells per horizon, R = 1000:**

| Model | fixed | open |
|---|---|---|
| no borrowing | 0.05 | 0.07 |
| power, equal | 0.04 | 0.07 |
| power, column | 0.04 | 0.07 |
| EX/NEX | 0.06 | 0.08 |
| similarity | 0.27 | 0.40 |
| two-way EX/NEX | 3.9 | 5.1 |
| oracle | 0.04 | 0.06 |
| **total** | **4.4** | **5.8** |

**Two-way cells.** These are the slowest:
* The balanced fixed cells take about 4.6 min each (7 fits).
* The unbalanced cells take 13–18 min (25 fits in the fixed design, up to 61 in the open design). There, basket
  (4,3) at rate 0.25 makes the trial last 120 time units, or up to 300 in the open design.
* The longest cell is about 0.3 h.

**Wall time with `sim.slurm`** (8 cores per horizon, the two horizons as two array tasks running at the same time):
The calculation gives ideal scheduling lower bounds of about 0.6 h (fixed) and 0.7 h (open).
Actual wall time includes worker imbalance and overhead; memory use was not profiled by this pilot.

**Uncertainty.**
* **Small pilot fits.** At R = 20, the unbalanced two-way fits took up to 2.9 x the line of B. The cause is the
  overhead per fit with few states, and growing tables when n > 30 in the open design.
* **At R = 1000** the cost per state is expected to dominate. Twice the prediction (about 9 and 12 CPU hours)
  is a provisional planning allowance, not an upper bound. More trials may need additional late fits;
  the open design reaches n = 150, beyond the n <= 30 calibration, and only three scenarios were measured.
* **Hardware.** The cluster's speed per core may differ from this machine's.

## Monte Carlo precision

The standard error of a rejection rate is sqrt(r(1 - r)/R):

| R | r = 0.05 | r = 0.2 | r = 0.5 |
|---|---|---|---|
| 500 | 0.0097 | 0.0179 | 0.0224 |
| 1000 | 0.0069 | 0.0126 | 0.0158 |
| 2000 | 0.0049 | 0.0089 | 0.0112 |

With R = 1000, a type I error of 0.05 is estimated to about ±0.014 (2 SE). R = 2000 halves the variance at about
twice the compute: an estimated 20 CPU hours in total, a few hours of wall time on 8 cores.
