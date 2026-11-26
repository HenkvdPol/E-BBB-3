# Entry-time redesign: design proposal (not implemented)

Status: proposal only, 2026-11-26.
* The guard in `Additional_files/legacy/additional_functions.R::run_trials_evalue()` stays: `redesign = "entry"` stops with an error.
* The new framework simulator (`additional_functions.R`) has no entry option yet.
* Nothing here is re-enabled until the implementation below exists and passes the tests in section 5.

## 1. What the option is meant to do

A cohort (basket) that opens during the trial should begin with a borrowed alternative (or prior) fitted to the
platform data **at the moment it opens**. Its first patient is then tested against an alternative that uses
everything known at entry. The disabled implementation instead set the alternative at the next global analysis.
By then the cohort may already have accrued patients, whose block was played with a stale alternative, and the
update arrives in the middle of what should have been the first block.

## 2. Where it sits in framework.tex

framework.tex allows the analysis times of basket c, τ_{c,0} < τ_{c,1} < …, to be any stopping times of the platform
filtration. Entry can therefore be its own analysis of basket c with an **empty block**:
* At the entry time e_c, set τ_{c,0} = e_c.
* Fit the borrowing model to F_{e_c}.
* Restrict the posterior to (p0, 1], freeze it, and use it for the block (e_c, τ_{c,1}].

No new theory is needed. The block factor of the empty block before entry is 1, and the prior of the first block is
F_{e_c}-measurable. What is missing is the **simulator**: the analysis times must be events in calendar time, not a
fixed global grid.

## 3. Event-driven simulator

* **Events, processed in time order:**
  * cohort entry (`entry`);
  * outcome observed (`outcome`), with an optional delay after enrolment;
  * scheduled analysis (`analysis`, global or per basket);
  * cohort closure.
* **State:** for every basket, the counts (y, n), the e-value E_c, the frozen prior π_c, status, and the time of its
  last analysis.

**Rules.**
1. At an `entry` event of basket c: fit the engine to the current data and set π_c. No factor is computed (empty
   block).
2. At an `analysis` event involving basket c:
   * multiply E_c by the factor of the outcomes observed since its last analysis, under π_c;
   * decide;
   * refit and set a new π_c.
3. At an `outcome` event: store the outcome only. It never changes π_c.
4. Ties at equal times are processed in a fixed, pre-specified order: entry, then outcomes, then analyses.
   Simultaneous outcomes are taken in enrolment order. framework.tex requires a fixed order; the order itself is a
   design choice.

**Vectorisation.**
* With deterministic accrual, the event times are the same in all trials. The event list can therefore be built once,
  and the vectorised simulator of `additional_functions.R` stays usable: its `times` become the merged list of analysis and entry
  times. At an entry time only the entering baskets are refitted.
* Random accrual needs per-trial event lists, which means per-trial loops or grouping by identical lists.

**Engines.** No change: every engine already returns the posterior for any data state. An entering basket has n = 0,
and its posterior is the borrowed prior predictive (`additional_functions.R` computes it through the other baskets' data).

## 4. Design-optimal designs (legacy)

* For the point-alternative designs, entry fixes θ1 at e_c; the dynamic-programming bets are then built for that θ1.
* Planned extension: `redesign = "entry"` means θ1 never changes after entry.
* This needs the same event list and belongs to the later 'extras' (`Additional_files/legacy/`).

## 5. Tests before re-enabling

1. **Predictability.** Perturb the outcomes of patients enrolled after e_c. π_c at entry is unchanged.
2. **Entry between analyses.** A cohort that opens between two global analyses and accrues patients before the next
   one has its first-block prior from F_{e_c}, not from the next analysis. Check this on a trace.
3. **Exact type I error.** Enumerate one basket entering at a non-analysis time, with adversarial other baskets, as in
   `Additional_files/tests/test_validity.R`. The error must be ≤ α.
4. **Regression.** With all entries at analysis times, results equal the current simulator.
5. **Ties.** Simultaneous entry and analysis follow the documented order.

## 6. Open choices for the researcher

* Whether the entry fit uses only other baskets' data (n_c = 0 anyway) or also external data D0. No D0 is used now
  (decision D8).
* Whether a cohort's own analyses are per basket (after every m patients of that basket) or global (calendar). The
  current simulator uses global calendar analyses.
* Outcome delays: the same event mechanism handles them. The conditional null must then hold given the history at
  enrolment (framework.tex, 'Scope of the guarantee').
