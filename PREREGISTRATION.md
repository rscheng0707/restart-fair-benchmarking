# Methodology Freeze — Phase 2 Pre-Commitment
**Project:** Time-Fair Benchmarking of Metaheuristics · ICMLSC 2027
**Frozen on:** 16 August 2026
**Status:** Phase 1 (exploratory pilots) closed. All decisions below are fixed before any main-matrix execution.

---

## Purpose of this document

The paper's contribution rests on testing claims that were asserted without evidence. That position is only defensible if our own analysis choices were fixed before we saw our results. This document records those choices, dated, so that the distinction between exploratory and confirmatory work is auditable rather than asserted.

**Any deviation after this date is reported as a deviation, with reason.** Deviations are not disqualifying; concealed deviations are.

> ### ⚠ AMENDED 16 August 2026 (same day, post-calibration)
> Calibration runs and a pre-execution review produced four design changes, all recorded in §11. In summary:
> **(D1)** Gallagher101 (f21) replaced by LunacekBiRastrigin (f24), same BBOB group.
> **(D2)** SAPSO/PSO cost ratio promoted from a fixed parameter to a **controlled experimental factor**.
> **(D3)** Cost-ratio levels set to 3×, 6×, 12× (originally 2×, 5×, 10×) — the 2× level is structurally unreachable.
> **(D4)** Multi-budget FE regime reduced to {0.1×, 1×}; the 10× level was 71.6% of total compute.
> **(D5)** Thread pinning and SMT demoted from experimental factors to fixed, reported execution conditions.
> All were made *before* any main-matrix execution, driven by calibration measurement and compute accounting rather than by results. Sections below are updated; original text is preserved in §11.

---

## 1. What Phase 1 established

Three pilot rounds (hand-written functions → scaled workload → BBOB suite). Findings that shaped the frozen design:

| Finding | Consequence |
|---|---|
| Hand-written test functions gave 20–140 ms runs; OpenMP overhead dominated, no thermal drift detectable in 90/90 configurations, 29/90 τ estimates unstable | Moved to the BBOB suite via IOHexperimenter; per-evaluation cost rose 4–27× |
| BBOB single-thread times span 0.136 s (Sphere) to 11.1 s (Gallagher101) | Workload now in range; Sphere requires a raised FE budget |
| Parallel efficiency at 16 threads varies 30.3% (Sphere) to 68.6% (Gallagher101) — a 2.26× spread | Motivating evidence that efficiency is problem-dependent. **Measured on the pilot-3 portfolio, which included Gallagher101 (since excluded, D1); the upper bound came from that function.** The main matrix re-measures efficiency on the final portfolio; RQ3 conclusions rest on those values, not these. |
| Thread pinning gives no benefit at 1–4 threads, +12–25% at 6–16 threads | Placement matters only past the intra-CCD boundary; pinning policy is a reported factor |
| SMT pass (24/32 threads) shows poor efficiency and unstable τ | SMT reported as a factor level, not a headline |
| Thermal drift significant in 4/130 configurations, concentrated in the longest runs | Warm-up/steady-state rule required, with empirical support |
| τ instability fell from 29/90 to 14/130 | τ fixable from pilots for most configurations; exceptions listed in §4 |
| SAPSO/PSO per-iteration cost ratio varies 1.05× to 6.65× by problem, reproduced across two independent runs | Cost ratio promoted to a controlled factor (D2); fixed parameters cannot hold the treatment constant |
| Surrogate overhead is near-constant (0.45–1.65 s) while PSO base cost spans 86× (0.29 s – 25.24 s) | Mechanistic explanation for the ratio collapse: a constant addend against an 86× range cannot yield a constant ratio |
| Under Rule B at fixed SAPSO parameters, k = 1 on four of five problems | The inversion mechanism could not operate on 80% of the portfolio — the reason for D2 |
| Gallagher101 measured at 25.24 s single-thread, 87% of total portfolio compute | Excluded (D1); its retention would have made the cost-ratio factor unaffordable |
| Gallagher101 cost differs ~20× between MSVC/Ryzen and GCC/x86 container while four other functions agree closely; Gallagher21 (21 peaks) costs 5.5 s vs Gallagher101 (101 peaks) at 25.2 s, ratio 4.6× ≈ peak-count ratio 4.8× | Cost scales with number of local optima, consistent with memory-bound behaviour. Recorded as a **substantive observation** (§12), not merely a nuisance |

Pilot and calibration data is exploratory throughout and is reported as such.

---

## 2. Algorithm portfolio — FIXED

| Slot | Algorithm | Role |
|---|---|---|
| Low-overhead baseline | PSO (`pso.hpp`), canonical spec with inertia weight | Restart beneficiary under Proposition 1; the protocol's own worked example |
| High per-iteration cost | Surrogate-assisted PSO (`sapso.hpp`), RBF multiquadric | Anchors coupling Rule B |
| RQ3 treatment | PSO, OpenMP implementation | Parallelism isolated by holding the algorithm fixed |
| RQ3 control | PSO, strictly sequential implementation | Verified bit-identical to OpenMP at 1 thread |

**Verification gate (already passed):** serial and OpenMP PSO produce bit-identical results under matched seeds at 1/2/4/8 threads, on both hand-written and BBOB problems. Re-run `test_equiv` after any code change. If it fails, the RQ3 comparison is confounded and no result may be reported.

**SAPSO parameters (AMENDED, D2).** `n_cand=3, refit_every=5, rbf_eps=1.0` fixed. **`archive_cap` is calibrated per (problem × target cost ratio)** to hold the treatment constant across problems — see §4.4. Calibration is performed *before* the main matrix, from calibration runs only, and the resulting parameter table is reported in full.

Rationale for the change from fixed parameters: surrogate overhead is near-constant while base evaluation cost spans 86× across the portfolio, so fixed parameters produce a cost ratio that varies 1.05×–6.65× by problem. That makes the *treatment* uncontrolled — the mechanism under study fires only where the ratio happens to be large. Calibrating the parameter to fix the ratio is standard controlled design, and is the opposite of tuning toward a preferred outcome.

SAPSO is not presented as a novel algorithm.

**Implementation provenance:** PSO implemented from the canonical published specification; SAPSO assembled from published surrogate-assisted components. Neither claimed as a contribution. Both released.

---

## 3. Problems and budgets — FIXED

BBOB suite via IOHexperimenter, instance 1, dimension 30.

| fid | Function | BBOB group | Measured PSO τ (s) | FE budget |
|---|---|---|---|---|
| 1 | Sphere | 1 — separable | 0.292 @10,000×D | **30,000 × D** (D6) |
| 8 | Rosenbrock | 2 — low/moderate conditioning | 0.779 | 10,000 × D |
| 10 | EllipsoidRotated | 3 — high conditioning, unimodal | 1.143 | 10,000 × D |
| 15 | RastriginRotated | 4 — multimodal, adequate global structure | 1.602 | 10,000 × D |
| 24 | LunacekBiRastrigin | 5 — multimodal, weak global structure | 0.842 | 10,000 × D |

All five BBOB groups are covered. τ values measured single-threaded, 15 runs, all with 95% CI half-width below 5%.

**Selection criterion:** coverage of the five official BBOB groups, fixed before results were seen.

### 3.0 Portfolio amendment (D1)
Gallagher101 (f21) was the original group-5 representative and is **excluded**. Measured at 25.24 s single-thread — 87% of total portfolio compute, and 5.4× the cost of the other four problems combined. Retaining it would have made the cost-ratio factor (D2) unaffordable at ~86 h; excluding it, together with D4 and D5, brings the full design to ≈6.5 h.

Replacement selected by measuring all four remaining group-5 functions on the target hardware:

| fid | Function | PSO τ (s) | τ CI | Note |
|---|---|---|---|---|
| 20 | Schwefel | 1.057 | ok | |
| 22 | Gallagher21 | 5.518 | SAPSO unstable (8.7%) | expensive; same family as excluded f21 |
| 23 | Katsuura | 1.734 | ok | |
| **24** | **LunacekBiRastrigin** | **0.842** | **ok (0.85%)** | **selected — cheapest, most stable** |

The exclusion is on measured compute grounds and is stated in the paper. Landscape coverage is unaffected: f24 occupies the same BBOB group as f21.

### 3.1 Sphere exception (D6) — FIXED
Sphere at the conventional 10,000×D gives a single-run time of **0.292 s** (measured, 15 runs). That is short enough to degrade τ precision: calibration observed SAPSO CV of 8.4% and 12.0% at the two lowest `archive_cap` values, against below 1% on most other configurations.

**FE multiplier set to 30,000 × D**, giving ≈0.88 s — inside the 0.5–1.0 s target band. Sphere's share of total portfolio compute rises from 6.3% to 16.7%; total base cost rises from 4.66 s to 5.24 s, which is affordable.

*(The earlier figure of 0.136 s in this document came from pilot 3, before the portfolio and swarm size were finalised; 0.292 s is the calibrated value and supersedes it.)*

Set once, from calibration data, before the main matrix. Not adjusted after results are seen. Reported in the paper as a stated departure from the uniform FE convention, with this reason.

**✅ VERIFICATION COMPLETE — AND IT FAILED, prompting D8.** The `archive_cap` values in §4.5b were calibrated at femult = 10,000. The cost ratio is a quotient of two times that both scale with FE count, so it *should* be approximately invariant to the multiplier — but this is an assumption, not a measurement. Re-run the Sphere calibration at femult = 30,000 and confirm the achieved ratios remain within ±10% of 3×/6×/12×:

```
calibrate.exe --fids 1 --caps 23,100,217 --femult 30000 --runs 5 > v1b.csv 2> v1b.log
python pick_caps.py v1b.csv --targets 3,6,12
```

Result: ratios shifted 15–17% upward, outside the band. Sphere's caps were recalibrated at femult = 30,000 (D8) and now read 17 / 83 / 199.

**Consequence worth carrying into the analysis.** Because the ratio is multiplier-dependent, the **multi-budget FE regime at 0.1× will not exhibit the nominal 3×/6×/12× ratios either** — it runs at one-tenth the calibrated budget. This is a consequence of the manipulation rather than a defect: the regime exists precisely to vary the budget. It must nevertheless be stated, and the paper should report the achieved ratios at 0.1× rather than the nominal labels. Measuring the 0.1× drift once per problem is cheap and preferable to a caveat.

Sphere is retained rather than dropped: it is the cheapest problem, and the per-evaluation cost spread is the mechanism under study.

### 3.2 FE convention
Base budget 10,000 × D follows CEC convention. **This convention is itself arbitrary** — a point made by Piotrowski et al. (2024) — which is why the multi-budget FE regime is included as a factor rather than treating one budget as authoritative.

---

## 4. Budget coupling and τ — FIXED

### 4.1 Coupling rules
All three reported. **Rule B is primary.**

- **Rule B (primary).** T = time for SAPSO to exhaust its FE allocation. Reproduces the protocol's own construction, testing the claim on its author's terms.
- **Rule A (robustness).** T = time for PSO to exhaust its FE allocation.
- **Rule C (control).** T = mean of the two algorithms' single-run durations **for that problem**: T = ½(τ_PSO + τ_SAPSO). Defined per problem rather than across the portfolio, because a portfolio mean would impose the same T on problems whose base cost differs by 5.5×, starving the expensive ones and inflating restart counts on the cheap ones. This is what `mainmatrix.cpp` implements.

Cut order if compute binds: C first, then the multi-budget FE regime. Rule B is never cut.

### 4.2 τ estimation
- Estimated from **separate pilot runs**, never from runs under analysis.
- Estimated **per configuration**: problem × implementation × thread count × toolchain × pinning policy.
- Fixed at freeze; held constant through the main matrix.
- **Source: `taupilot.cpp`**, run once under the fixed pinning policy (D5), producing `tau.csv` for every (fid × ratio × threads × impl) cell the main matrix visits. The pilot-3 CSVs do **not** serve this purpose: they predate the cost-ratio factor and were measured on a portfolio containing Gallagher101, since excluded.
- `mainmatrix.cpp` reads τ from file and **exits with a fatal error if any configuration is missing**. It never estimates τ itself — doing so would fit the restart count to the data it evaluates.
- τ is tied to the pinning policy. **If the policy changes, τ must be re-measured.**

**Unstable configurations.** Any cell whose 95% CI half-width exceeds 5% of τ is re-run with more pilot repetitions, **or** excluded with the exclusion recorded. No configuration proceeds with an imprecise τ. `taupilot` flags these on stderr.

### 4.3 Restart accounting
Both accountings implemented and both reported:
1. Post-hoc count of runs completed within T (Algorithm 1 as written)
2. Pre-computed k = ⌊T/τ⌋ (Definition 2)

The protocol is ambiguous about which is authoritative; reporting both is the honest resolution.

**Restart admission policy:** a new run starts only if elapsed + τ ≤ T. The protocol does not specify a policy; this one is stated, parameterised (`admit_only_if_fits`), and reported, since the choice measurably affects results.

---

### 4.4 Cost ratio as a controlled factor (D2, levels amended D3)

**Levels: 3×, 6×, 12×** (SAPSO per-iteration cost relative to PSO), calibrated per problem by setting `archive_cap`.

*Original levels were 2×, 5×, 10×; amended to 3×, 6×, 12× — see D3.* The 2× level proved **structurally unreachable**: SAPSO carries a floor cost independent of archive size (candidate generation and screening at `n_cand=3`), and on Sphere that floor alone gives a ratio of 2.64× at `archive_cap=16`, the smallest value tested. Reaching 2× would require varying `n_cand` as well, abandoning the single-knob mechanism that keeps the cost driver interpretable.

Measured achievable band across the portfolio: **2.64× – 22.25×**, bounded below by Sphere's floor and above by LunacekBiRastrigin's ceiling within the swept range. The chosen levels sit inside that band on every problem and preserve approximately geometric spacing.

*Why this replaced fixed parameters.* Measured at `archive_cap=100`, the ratio varied 1.05×–6.65× across the portfolio, giving k = ⌊T/τ⌋ = 1 on four of five problems under the primary coupling rule. Where k = 1 the baseline receives no restarts, the protocol's mechanism cannot operate, and no inversion is possible *by construction*. The study would have tested its central question on one problem out of five.

*Why the ratio varies.* Surrogate overhead is near-constant (0.45–1.65 s measured) while PSO base cost spans 86× (0.29–25.24 s). A constant addend against an 86× range cannot produce a constant ratio. This is a property of surrogate assistance, not of our implementation, and is reported as such.

*What the factor buys.* RQ4 changes from "does the inversion occur?" — unanswerable at k = 1 — to **"at what per-iteration cost ratio does time-fairness begin to invert FE-based verdicts?"** A quantitative threshold, answerable in either direction, and consistent with the treatment of coupling rules and restart accountings elsewhere in this design.

*Discipline.* `archive_cap` values are calibrated from calibration runs only, fixed before the main matrix, and reported as a table (problem × target ratio × archive_cap × achieved ratio). Calibrating a parameter to hold the treatment constant is controlled design; adjusting it after seeing results would not be, and is prohibited.

### 4.5 Measured calibration values — FIXED

Single-threaded, 15 runs, dim 30, swarm 120, FE = 10,000×D, `archive_cap=100`. All τ CI half-widths below 5%.

| fid | Problem | PSO τ (s) | SAPSO τ (s) | ratio | overhead (s) |
|---|---|---|---|---|---|
| 1 | Sphere | 0.292 | 1.942 | 6.65 | 1.65 |
| 8 | Rosenbrock | 0.779 | 1.234 | 1.58 | 0.45 |
| 10 | EllipsoidRotated | 1.143 | 1.597 | 1.40 | 0.45 |
| 15 | RastriginRotated | 1.602 | 2.069 | 1.29 | 0.47 |
| 24 | LunacekBiRastrigin | 0.842 | 1.286 | 1.53 | 0.44 |
| ~~21~~ | ~~Gallagher101~~ | ~~25.240~~ | ~~26.428~~ | ~~1.05~~ | ~~1.19~~ |

Reproducibility check: an independent 3-run repetition agreed to within 0.07% on Gallagher101 (25.240 vs 25.223) and reproduced all ratios. The 11.1 s figure for Gallagher101 recorded during pilot 3 is superseded — pilot 3 executed configurations interleaved, giving different cache and thermal state.

**τ table for the main matrix** (per problem × implementation × thread count × toolchain × pinning) is populated from the pilot 3 CSVs and recorded in the freeze appendix before first execution.

### 4.5b `archive_cap` calibration — VERIFIED AND FIXED

Verified 16 August 2026. Single-threaded, dim 30, swarm 120, FE = 10,000×D, 5 runs per point. **Acceptance criterion (achieved ratio within ±10% of target) met at all fifteen points.**

| fid | Problem | cap @3× | achieved | cap @6× | achieved | cap @12× | achieved |
|---|---|---|---|---|---|---|---|
| 1 | Sphere | **17** | **2.96** | **83** | **5.94** | **199** | **12.36** |
| 8 | Rosenbrock | 430 | 3.00 | 812 | 5.99 | 1285 | 12.03 |
| 10 | EllipsoidRotated | 580 | 2.96 | 1027 | 5.93 | 1604 | 11.98 |
| 15 | RastriginRotated | 705 | 2.94 | 1242 | 6.01 | 1908 | 12.09 |
| 24 | LunacekBiRastrigin | 473 | 3.06 | 893 | 6.21 | 1401 | 12.48 |

Worst deviation 4.04% (f24 @12×); best 0.06% (f8 @3×). **Achieved ratios, not targets, are the reported quantity.**

**Sphere recalibrated at femult = 30,000 (D8).** Its row above shows the recalibrated values. The originals (23/100/217, calibrated at femult = 10,000) produced 3.26×/6.97×/14.06× once the multiplier was raised — a 15–17% upward shift, outside the acceptance band. Recalibrated caps achieve 2.96×/5.94×/12.36× with CV 0.2–0.8%, the tightest measurements in the calibration set.

**Implied restart count under Rule B** (k = ⌊achieved ratio⌋):

| Problem | 3× | 6× | 12× |
|---|---|---|---|
| Sphere | 2 | 5 | 12 |
| Rosenbrock | 3 | 5 | 12 |
| EllipsoidRotated | 2 | 5 | 11 |
| RastriginRotated | 2 | 6 | 12 |
| LunacekBiRastrigin | 3 | 6 | 12 |

k now spans 2–12 across the portfolio, against k = 1 on four of five problems at fixed `archive_cap=100`. The restart mechanism can operate on every problem at every level — which is what D2 was introduced to achieve.

**Recorded exception:** Sphere showed SAPSO CV of 8.4% (cap 23) and 12.0% (cap 100), above the 5% guideline. Cause is its short absolute runtime (0.292 s base), which inflates relative variance. Both points nonetheless landed within 6.1% of target. Either re-measure with additional runs before Session 1, or report as a stated exception with the CV values disclosed. All other thirteen points are below 6% CV, eleven below 1%.

**Calibration parameters used:** `n_cand=3, refit_every=5, rbf_eps=1.0` throughout; `archive_cap` is the sole calibration knob, preserving a single interpretable cost mechanism.

**Note for the methods section.** Sphere requires `archive_cap=217` to reach 12× where RastriginRotated requires 1908 — an 8.8× spread driven by a 5.5× spread in base evaluation cost. This is the direct evidence that a single fixed parameter cannot hold the cost ratio constant across a benchmark suite, and it justifies per-problem calibration in one sentence.

### 4.6 Compute budget — FIXED

Full design: 3 budget regimes × 3 coupling rules × 3 cost-ratio levels × 7 thread levels (+ serial reference at 1 thread) × 2 algorithms × 15 seeds × 5 problems.

**Original estimate: ≈6.5 h total, ≈1.3 h per session.** *Superseded by measurement.*

**Measured (Session 1): 258.8 min ≈ 4.3 h per session, ≈21.6 h total.** The estimate assumed each time-fair condition costs T. It does not: under the D9 overrun policy, SAPSO runs to completion whenever τ_SAPSO > T, so Rules A and C cost τ_SAPSO rather than T. SAPSO accounts for 192.5 of the 258.8 min; the time-fair regime accounts for 204.7 min.

Remaining four sessions: ≈17 h. Sessions are independent and seed-blocked, so they may be run on separate days without confounding (§6.0).

Enumerated condition count: 4,320 per session, 21,600 total.

Cut order if compute binds further: coupling Rule C first, then the multi-budget FE regime entirely, then thread levels reduced to {1, 4, 8, 16}. **The cost-ratio factor is never cut** — it is what makes RQ4 answerable.
## 5. Experimental factors — FIXED

| Factor | Levels |
|---|---|
| Budget regime | FE-matched · time-fair (restart-fair) · multi-budget FE (0.1×, 1× — see D4) |
| Coupling rule | B (primary) · A · C |
| **Cost ratio** | **3× · 6× · 12×** (calibrated per problem, §4.4–4.5b) |
| Algorithm | PSO · SAPSO |
| Implementation | serial · OpenMP |
| Threads | 1 · 2 · 4 · 6 · 8 · 12 · 16 |
| ~~SMT level~~ | *fixed condition, not a factor — see D5* |
| ~~Pinning~~ | *fixed condition, not a factor — see D5* |
| Problem | fids 1, 8, 10, 15, 24 |
| Seeds | 15 independent runs, executed as 5 sessions × 3 seeds (§6.0) |

**Toolchain held fixed** at the MSVC configuration used across pilots 1–3, and reported. A second toolchain pass is optional and, if run, is labelled distinctly — compiler is a factor, not an implementation detail. `-march=native` / `/arch:` beyond the recorded baseline is prohibited unless recorded as a factor level.

**Fixed execution conditions (D5), reported but not varied:**
- **Thread pinning:** `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, applied to `taupilot` and every session. Justified by the pilot-3 affinity comparison (three passes): no effect at 1–4 threads, 12–25% improvement at 6–16 threads, 4/5 problems responding, with thread optima largely unchanged. Pinning reduces variance without moving rankings, so it is controlled rather than sampled — the same treatment given to active cooling.
- **SMT:** excluded. Thread levels stop at physical core count. Pilot 3 found poor parallel efficiency and unstable τ at 24/32 threads under `OMP_PLACES=cores`; the pinned-SMT pass is reported as pilot evidence, not carried into the main matrix.

Varying either as a factor would require executing the full matrix once per policy, with a separate τ file each — 13 h rather than 6.5 h — to re-derive a result the pilots already establish.

**Estimated compute:** ≈6.5 h total, ≈1.3 h per session (§4.6), executed as five seed-blocked sessions (§6.0). *Splitting by coupling rule or by any other factor is prohibited — it would confound session with that factor (§6.0).*

---

## 6. Execution protocol — FIXED

### 6.0 Session blocking

The main matrix (~14 h) is executed in **five sessions of three seeds each**:

| Session | Seeds | Est. duration |
|---|---|---|
| 1 | 1–3 | ~2.8 h |
| 2 | 4–6 | ~2.8 h |
| 3 | 7–9 | ~2.8 h |
| 4 | 10–12 | ~2.8 h |
| 5 | 13–15 | ~2.8 h |

**Split by seed, not by factor.** Each session executes the *complete* condition set — every problem, algorithm, implementation, thread level, coupling rule, and budget regime — at three seeds. Session is therefore a balanced block: every condition appears in every session, and any session-level disturbance (ambient temperature, background load, machine state) affects all conditions equally rather than landing on one factor level.

Splitting by coupling rule, problem, or thread count was considered and rejected: it would confound session with that factor, and no post-hoc analysis could separate them.

**Session index is recorded in every output row**, permitting an explicit test for session effects. Absence of a session effect is reported as a robustness check; presence is reported and controlled for.

**Warm-up is repeated at the start of every session** (§6.2). A cold machine on one evening is not equivalent to a warm machine on another; the warm-up rule applies per session, not once per experiment.

**Partial completion is analysable.** Any completed session constitutes a valid experiment at reduced replication. If execution is interrupted after session *k*, the dataset supports analysis at 3*k* seeds with correspondingly reduced statistical power, stated explicitly. This is the reason for the seed-wise split rather than a single 14-hour block.

**The split is fixed at freeze.** Session boundaries are not adjusted after any results are seen.

### 6.1 Standing requirements
1. **Active cooling** in place for every timed run. Mandatory.
2. **Warm-up:** a warm-up load precedes each measurement batch until CPU temperature and clock stabilise within a stated tolerance. Pre-steady-state runs are discarded. Warm-up duration and stabilisation criterion are recorded. Interleaving distributes thermal bias but does not remove it.
3. **Interleaved and randomised** (problem × configuration) order, ordering seed recorded.
4. **Logging:** temperature, clock, and throttling events throughout; ambient temperature per batch.
5. **Environment:** mains power; sleep and automatic restart disabled; background services minimised; machine otherwise unused.
6. **Checkpointing:** results flushed per configuration so an interrupted run costs one cell, not the batch.
7. **Reported per §3.8:** CPU model, core count, cache configuration, SIMD level, BIOS settings, BLAS, compiler and flags, OS and runtime versions.

---

## 7. Analysis plan — FIXED

**RQ1/RQ4 (ranking, inversion).** Algorithm orderings under each budget regime; reversals identified per problem and per coupling rule.

**Unit of analysis (D12) — applies to every inference below.** The inversion verdict is constant across conditions within a (problem × seed) block in the large majority of cases; thread count, implementation, and coupling rule alter timing without altering which algorithm wins. The independent unit is therefore the **(problem × seed × cost-ratio) block**, and:
- Block-level rates are computed first; condition counts are never used as a sample size.
- The ordered cost-ratio factor is tested with **Friedman's test** across blocks, which preserves the pairing (all three levels are measured on the same blocks). An unpaired χ² on pooled conditions discards that pairing and is not used.
- Interval estimates use block counts. With five sessions the design yields ≈61 blocks — the real statistical power of the study, against a superficially reassuring ~1,470 conditions.
- Any per-problem comparison rests on ~12 blocks per problem and is reported as exploratory unless it survives a block-level test.

**Tie rule (D11), applied before any inversion is counted.** A paired comparison is excluded when the two algorithms' `best_fitness` values agree within `rel_tol=1e-9, abs_tol=1e-12`, on either the FE-matched or the time-fair side. A tie expresses no verdict and therefore cannot be reversed; including it would inflate the inversion count. Excluded comparisons are **reported as a count alongside every inversion rate**, since their number is itself informative — see the Sphere saturation note below.

**Sphere saturation (observed, sessions 1–2).** All 120 excluded comparisons fell on Sphere (fid 1); the other four problems produced zero ties. Both algorithms reach the same objective value on Sphere at the D6 budget of 30,000×D, so it discriminates between them not at all. This is reported as a benchmark-selection finding: raising Sphere's FE budget for timing precision (D6) pushed it into saturation for outcome purposes. Sphere is retained in the design — its timing data remains valid and it is the cheapest problem, which the cost-spread argument requires — but it contributes no inversion evidence and this is stated rather than absorbed silently.

**Supporting evidence for RQ2 (three strands, of differing strength).**
1. *Controlled:* pilot-3 affinity study — three passes, same seeds, pinning policy varied within pass.
2. *Physical mechanism:* inter-CCD thermal asymmetry of 10.2 °C measured during Session 1 (§11b), giving a reason why placement affects achievable clock.
3. *Controlled within-seed:* the Session 2 pinned-vs-unpinned comparison (§11c) — identical seeds and conditions, policy varied — in which PSO and SAPSO move in opposite directions from 6 threads upward. Small magnitudes (a few percent), direction stable across estimators, and no effect on verdicts.

**RQ2 (stability).** Rank correlation computed *within* each problem between condition pairs, then pooled across problems — effective n is (problems × condition-pairs), not portfolio size. Reported with the distribution of per-problem τ, not only its mean. Supplemented by a pairwise reversal count with binomial confidence intervals, which does not degrade at small portfolio size. Bootstrap intervals over problems and seeds.

**Statistical power is stated explicitly, at the block level (D12).** The design yields ≈61 independent blocks, not ~1,470 conditions; power is reported against the former. Where the design cannot separate a small ranking shift from noise, that limit is reported as a bound on what the study concludes.

**Budget overrun (D9).** The realised wall-clock time of every time-fair condition is reported as a ratio to its nominal budget T, broken down by coupling rule and algorithm. Where an algorithm's single-run duration exceeds T, the protocol as specified provides no result and our policy admits one run regardless; the resulting overrun is reported as a measured quantity, not smoothed away.

**Consequence to state explicitly:** under Rule A, where T is set by the low-overhead baseline, the budget does not constrain the expensive algorithm at all — it overran by a mean factor of 33.7× in Session 1. Rule A therefore functions as a demonstration that a fixed-time budget is unenforceable against an algorithm whose minimum granularity exceeds it, rather than as a working comparison. This is a finding about the protocol's specification, and it is reported as such.

**RQ3 (parallelism).** Parallel efficiency per problem; evaluations obtained within T as a function of thread count; the resulting spread in effective budget across problems.

**RQ5 (Proposition 1).** Whether the low-overhead algorithm's restart count under time-fairness corresponds to improved outcomes. Framed as *testing a stated proposition*, not debunking.

**Session effects.** Wall-clock measurements are tested for a session effect (§6.0) before pooling. Because session is a balanced block — every condition appears in every session — a session effect inflates variance but does not bias any factor comparison. If present, it is reported and controlled for; if absent, that is reported as a robustness check.

**Precision discipline.** All analysis at full precision; rounding for display only. Apparent threshold effects verified unrounded before being claimed.

**Saturation.** If orderings are identical across all conditions, that is reported as the finding, not buried as a limitation.

---

## 8. Pre-registered interpretation

Fixed before results exist so that framing cannot be retrofitted.

| Outcome | Reading |
|---|---|
| Inversion under all coupling rules | Robust support for the protocol's central claim |
| Inversion under B only | Holds on its author's terms but sensitive to budget anchoring — a boundary condition, more informative than bare confirmation |
| Inversion under none | Claim is not general; conditions of failure are the contribution |
| Rankings stable across conditions | Reproducibility objection is weaker than assumed; time-fair benchmarking is safe to adopt |
| Rankings unstable | Objection substantiated and quantified; boundary conditions specified |
| Parallelism dominates | Definition 1's permission of internal parallelism needs constraint; concrete amendment proposed |
| Proposition 1 holds / fails | Support for, or boundary condition on, a stated but untested proposition |

Every cell is publishable. That property was designed in.

**Title selected at Gate 3**, once results exist, and must be a sentence the data carries without qualification.

---

## 9. Scope exclusions — FIXED

**GPU/CUDA implementation is out of scope.** Definition 1's parallelism clause extends to accelerators, and a GPU implementation would sharpen RQ3 considerably. It is excluded because (a) BBOB functions would require device-side reimplementation with bit-exact agreement against IOHexperimenter, a substantial correctness risk, and (b) the timeline does not permit it without endangering the rest. Stated in the discussion as the sharper instance left for future work, with the measured CPU efficiency spread as motivation.

**Second physical platform is out of scope.** Only one workstation is available. RQ2 is reframed as within-platform heterogeneity — thread count, toolchain, pinning policy, thermal state, software stack — varied one at a time on fully-characterised hardware. Justified on measurement control, not cost: shared cloud environments have non-deterministic hardware allocation, unobservable co-tenant variance, vCPU limits precluding the thread factor, and unverifiable BIOS/cache/turbo state. A reduced supplementary pass on a hosted notebook is reported separately as a low-control observation.

---

## 10. Pre-flight checklist

Before the first main-matrix run:

- [x] Citing check re-run on arXiv:2509.08986 — **still v1, unrevised, no peer-reviewed version, no empirical instantiation found (16 Aug 2026)**
- [x] Cost ratio measured on target hardware, 15 runs — §4.5
- [x] All calibration τ within the 5% threshold
- [x] Portfolio finalised (D1); group-5 replacement measured and selected
- [x] Cost-ratio factor levels set (D2)
- [ ] `make check` / equivalence gate re-run on current code
- [x] `archive_cap` calibration sweep — 45 points across two passes; achievable band measured at 2.64×–22.25×
- [x] Cost-ratio levels revised to 3×/6×/12× (D3)
- [x] **Verification run complete (§4.5b)** — all 15 points within ±10% of target; calibration table FROZEN
- [x] Sphere FE multiplier set to 30,000×D (D6); drivers updated
- [x] Pinning and SMT demoted to fixed conditions (D5)
- [x] Rule C definition reconciled with `mainmatrix.cpp` (per-problem mean)
- [x] Main-matrix driver written (`mainmatrix.cpp`) and τ pilot written (`taupilot.cpp`)
- [x] **Sphere re-calibration at femult=30,000** — check failed; caps recalibrated to 17/83/199 (D8)
- [ ] Optional: measure achieved cost ratios at the 0.1× multi-budget level, so the paper reports actual rather than nominal ratios
- [ ] **Run `taupilot.exe` under the fixed pinning policy** → `tau.csv`; resolve any flagged unstable cells
- [ ] `--dry-run` on session 1 to confirm the condition count (expect 3,120 on 16 cores)
- [x] Session 1 thermal record captured and analysed (§11b) — no throttling, warm-up validated
- [x] Session 1 complete (3,120 conditions, 258.8 min)
- [x] **Session 2 re-run under the fixed pinning policy** (D10) — 3,120 conditions, 253.6 min, `OMP_PROC_BIND=close OMP_PLACES=cores` confirmed in log
- [x] Tie-handling rule specified and locked (D11) before Session 3
- [x] Session 3 (seeds 7–9) — 250.1 min, pinning confirmed
- [x] Session 4 (seeds 10–12) — 277.2 min, pinning confirmed
- [x] Unit-of-analysis rule specified and locked (D12) before Session 5
- [ ] Session 5 (seeds 13–15) — final session
- [ ] Set `OMP_PROC_BIND` / `OMP_PLACES` **inside** the run script; verify the echoed values in the log before leaving any session unattended
- [ ] Sessions 2–5: log thermal data for the **full session duration**, not the first 70 min
- [x] Ambient conditions recorded: air-conditioned room at a 26 °C setpoint, constant across all sessions (§11b)
- [ ] Confirm `OMP_PROC_BIND=close` and `OMP_PLACES=cores` are set for **both** `taupilot` and every session — τ is tied to the pinning policy

- [ ] Sphere FE multiplier set and recorded (τ = 0.292 s is short; target 0.5–1.0 s)
- [ ] Warm-up duration and stabilisation criterion set and recorded
- [ ] Active cooling in place
- [ ] Sleep and automatic restart disabled
- [ ] Checkpointing verified by interrupting a short run
- [ ] Environment metadata captured per §6.7
- [ ] Main-matrix driver written (does not yet exist; needs `--seeds-from`, `--seeds-to`, `--session`, coupling-rule and cost-ratio sweeps)
- [ ] Citing check re-run once more immediately before submission

## 11. Post-freeze deviations log

| # | Date | Deviation | Reason | Made before results? |
|---|---|---|---|---|
| D1 | 16 Aug 2026 | Gallagher101 (f21) replaced by LunacekBiRastrigin (f24) | f21 measured at 25.24 s single-thread — 87% of total portfolio compute, 5.4× the other four problems combined. Retention would have made the cost-ratio factor unaffordable (~86 h vs 13.9 h). Replacement selected by measuring all four remaining group-5 functions; f24 cheapest (0.842 s) with the most stable τ. BBOB group coverage unchanged. | Yes — calibration only |
| D2 | 16 Aug 2026 | SAPSO/PSO cost ratio promoted from fixed parameter to controlled factor (2×, 5×, 10×) | At fixed `archive_cap=100`, ratio varied 1.05×–6.65× and k = 1 on four of five problems, so the inversion mechanism could not operate on 80% of the portfolio. Fixed parameters leave the *treatment* uncontrolled. Calibrating `archive_cap` per problem holds the treatment constant and converts RQ4 from a binary question into a threshold. | Yes — calibration only |

| D3 | 16 Aug 2026 | Cost-ratio levels changed from 2×/5×/10× to **3×/6×/12×** | The 2× level is unreachable by construction: SAPSO's floor cost (candidate generation and screening at `n_cand=3`, independent of `archive_cap`) gives a minimum ratio of 2.64× on Sphere at `archive_cap=16`. Measured achievable band across the portfolio is 2.64×–22.25×. Reaching 2× would require varying `n_cand` as a second knob, losing the single-mechanism interpretation. Revised levels sit inside the achievable band on all five problems with approximately geometric spacing. | Yes — calibration only |

| D4 | 16 Aug 2026 | Multi-budget FE regime reduced from {0.1×, 1×, 10×} to **{0.1×, 1×}** | Full enumeration costs 22.3 h, of which the 10× level is 71.6%. Removing it gives 6.5 h total (1.3 h per session), preserving the session-blocking design (§6.0) at evening-sized blocks. All research questions are retained; the multi-budget regime still functions as the rival to time-fairness (Piotrowski et al., 2024). **Stated limitation:** the regime now spans one order of magnitude rather than the three orders Piotrowski et al. advocate, so this study tests budget sensitivity but cannot claim their range. | Yes — before execution |

| D5 | 16 Aug 2026 | Thread pinning and SMT demoted from experimental factors to **fixed, reported execution conditions** | Varying pinning as a factor requires running the full matrix once per policy with a separate τ file each (13 h vs 6.5 h) to re-derive what pilot 3 already established across three passes: no effect at 1–4 threads, 12–25% at 6–16 threads, thread optima largely unchanged. Pinning moves timing magnitude, not rankings — and RQ2 concerns rankings. Fixing it also removes a variance source, as with active cooling. SMT excluded on the same basis: pilot 3 found poor efficiency and unstable τ at 24/32 threads. | Yes — before execution |
| D6 | 16 Aug 2026 | Sphere FE multiplier raised from 10,000×D to **30,000×D** | At 10,000×D, single-run time is 0.292 s and τ precision degrades (SAPSO CV 8.4% and 12.0% at the two lowest caps, against <1% elsewhere). 30,000×D gives ≈0.88 s, inside the target band. Raises Sphere's share of compute from 6.3% to 16.7%. Anticipated in the original freeze (§3.1) as an exception to be quantified; this fixes the value. **Conditional on the §3.1 verification check.** | Yes — before execution |

| D7 | 16 Aug 2026 | PSO conditions under the FE-matched and multi-budget regimes emitted **once**, not once per cost-ratio level | PSO does not use `archive_cap`, so under budget-only regimes its run is identical across cost-ratio levels — same parameters, same FE budget, same seed. Emitting per level produced 2,400 duplicate rows across the experiment. Under the **time-fair** regime PSO does differ by level, because T derives from SAPSO's runtime at that level and sets a different restart count k; those remain enumerated per level. Condition count per session falls 4,320 → 3,120. Found in pre-execution code review. | Yes — before execution |

| D8 | 18 Aug 2026 | Sphere `archive_cap` recalibrated at femult = 30,000: **17 / 83 / 199** (was 23 / 100 / 217) | The §3.1 verification check required by D6 **failed**: at the raised multiplier the original caps gave 3.26×/6.97×/14.06×, a 15–17% shift outside the ±10% band. **The cost ratio is not invariant to the FE multiplier** — surrogate overhead does not scale with FE count in the same proportion as base evaluation cost. Recalibration at the operative multiplier achieves 2.96×/5.94×/12.36% with CV 0.2–0.8%. D6 itself is vindicated: Sphere's CV fell from 8.4%/12.0% to below 1%. | Yes — before execution |

| D9 | 18 Aug 2026 | **Budget-overrun policy fixed as "complete the first run regardless of T"**, retaining the behaviour already exercised in Session 1 | The protocol does not specify what happens when an algorithm cannot complete a single run within the fixed-time budget T. Definition 1 requires aggregate execution to stay within T; Definition 3 sets ERT = ∞ and reports the success rate when no run attains the target — but neither states whether a run that *cannot* fit should be attempted at all. Our implementation admits the first run unconditionally and admits subsequent runs only if `elapsed + τ ≤ T`. Consequence, measured in Session 1: SAPSO overruns T by a mean factor of **33.7× under Rule A** (max 133×), **2.11× under Rule C**, and **1.16× under Rule B**. The alternative reading — treat an unfinishable run as did-not-finish — was considered and rejected on the grounds that Session 1 was already complete under the current policy and the primary rule (B) is barely affected. **This is the first amendment made after execution began**; it fixes existing behaviour rather than changing it, and no data is invalidated. | **No — after Session 1** |

| D10 | 19 Aug 2026 | **Session 2 executed without the fixed pinning policy** — `OMP_PROC_BIND` and `OMP_PLACES` were unset. **Operator error, not a design decision.** | The session log records `OMP_PROC_BIND=(unset) OMP_PLACES=(unset)` against Session 1's `close`/`cores`. The variables were set in the shell for Session 1 and not re-exported for Session 2; the batch script did not set them internally. **Consequence:** pinning policy is perfectly confounded with seed set (seeds 1–3 pinned, seeds 4–6 unpinned), which is exactly the assignment §6.0 prohibits — a factor level mapped onto a session block rather than varied within it. τ in `tau.csv` was measured pinned, so T and k_predicted are unaffected, but realised restarts fell short of k_predicted by 5.4% in Session 2 against 2.8% in Session 1. **Resolution: Session 2 is re-run under the correct policy.** The unpinned data is retained and reported as an unplanned pinning observation (§11c), not discarded. | **No — error during execution** |

| D11 | 19 Aug 2026 | **Tie-handling rule specified: a comparison is excluded when the two algorithms' `best_fitness` values agree within tolerance, on either side of the pairing.** Tolerance is `math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)`. | The inversion measure compares which algorithm wins under FE-matched budgets against which wins under time-fair budgets. No rule was specified for the case where the two algorithms return the **same** objective value. Treating a tie as "PSO does not win" — the implicit behaviour of a bare `<` comparison — counts a tied FE verdict paired with a time-fair PSO win as an *inversion*, but a tie is not a verdict and there is nothing to reverse. **A tie is therefore excluded from the denominator.** Tolerance rather than exact equality, because floating-point results that are algorithmically identical need not be bitwise identical. **Impact when applied to sessions 1–2:** 120 of 720 comparisons excluded, all on Sphere (fid 1); pooled inversion rate 32.1% → 29.7%; and the apparent cost-ratio trend (χ²=6.00, p=0.0497) **fails** under the corrected rule (χ²=2.36, p=0.31). The trend claim is withdrawn. | **No — after Sessions 1–2 were analysed** |

| D12 | 20 Aug 2026 | **Unit of analysis for the inversion measure fixed as the (problem × seed × cost-ratio) block, not the individual condition.** Significance tested with paired methods across blocks (Friedman for the ordered cost-ratio factor); binomial intervals computed on block counts. | Inversion verdicts are **not independent across conditions within a block**: measured over sessions 1–4, the verdict is constant across all 24 conditions in 42 of 49 (problem × seed) blocks, and constant across all 8 conditions in 136 of 147 (problem × seed × ratio) blocks. Thread count and implementation change timing but almost never change which algorithm wins. Treating each of the 1,176 conditions as an independent observation is pseudoreplication: it understates confidence intervals by roughly √(24) and invalidates the χ² tests reported for sessions 1–3. **Effective sample is 49 blocks, not 1,176 comparisons.** **Consequence:** the cost-ratio trend, reported as non-significant under the unpaired test (χ²=2.93, p=0.23), is **significant under the correct paired test (Friedman χ²=13.04, p=0.0015)** — the unpaired test discarded the pairing that the three ratio levels are measured on the same blocks. Conversely, the apparent problem effect (χ²=30.30, p=4×10⁻⁶) is **not supported** once clustering is accounted for; per-problem rates range 0–76% across sessions on ~12 blocks each. | **No — after Sessions 1–4 were analysed** |

| D13 | 6 Nov 2026 | **Three tables in the submitted manuscript were computed from `session2.csv` — the voided unpinned run — instead of `session2_pinned.csv`.** Corrected values recorded below; to be applied at camera-ready. | Discovered during post-acceptance verification of the submitted PDF against the raw session data.

**What is verified.** The affected tables reproduce exactly, to four decimal places, from the four unaffected sessions plus the *unpinned* session 2. Three alternative hypotheses were tested and excluded: pooling both session-2 files gives 1.2316 s / 0.1216 s / 10.13× and rule-B rates of 95.3 / 32.6; omitting session 2 entirely gives 1.2377 s / 0.1215 s / 10.19× and 95.6 / 27.8; the pinned file gives 1.2430 s / 0.1200 s / 10.35× and 94.9 / 29.6. Only the unpinned substitution matches the printed values. So the dataset used is established beyond doubt.

**What is not verified: why.** The reason is not recoverable from the artefacts, and the author does not recall making a deliberate choice. The most likely mechanism, stated here as a hypothesis and not a finding, is that the timing tables were computed before the pinned re-run existed and were not refreshed after it, while the verdict and inversion tables were computed later from the pinned file. The split pattern is consistent with that — a deliberate decision would more plausibly have been applied to all tables rather than to the three timing tables alone. An earlier revision of this amendment asserted "analyst error, not a design decision"; that asserted intent on no evidence and has been withdrawn.

**The correction is required under either explanation.** Even the most charitable reading — that the unpinned configuration is the one most practitioners actually run, so it is the more externally valid basis for a scaling table — leaves the manuscript internally inconsistent in three ways: §3.7 states that the unpinned duplicate "is only a robustness check, not a crossed-placement experiment"; §4.8's pinning comparison presupposes that the pinned run is the main dataset and the unpinned one the deviation; and the paper would then derive its 10.18× headline from the pinned `tau.csv` while deriving the 9.96× scaling figure beside it from unpinned timings, mixing the two sources inside a single argument. D10 mandated the re-run precisely because the unpinned session confounds placement policy with seed set. Recomputing from the pinned data is therefore the correction, whatever the original reason. **The error is self-contradicting:** §3.7 of the manuscript states that "an inadvertently unpinned duplicate of session 2 is only a robustness check, not a crossed-placement experiment," so the paper declares in its own methodology that these data are excluded from the main analysis and then uses them in three tables. **Affected:** Table 1 (wall-time scaling), Table 5 (completion within *T*), Table 6 (overrun ratios), and the realised-restart means in §4.6. **Not affected:** Table 2 (predicted *k*, including the 10.18× headline — identical under both files, because *k* derives from the pinned `tau.csv`), Tables 3, 4 and 7, the 27.0% pooled inversion rate, the Friedman result (χ²=16.75, *p*=0.00023), §4.1's flat CPSO scaling, the 312 excluded ties, and §4.8's pinning comparison (which is *intended* to use both files). **No conclusion changes under correction.** The 1→16-thread PSO speedup rises from 9.96× to 10.35×, marginally strengthening the parallelism argument; Rule B completion falls from 31.9% to 29.6%, still "just under three runs in ten"; the 142.45× maximum overrun is identical under both files. | **No — after submission** |

**D13 corrected values.** Mean total wall time (s), fixed-FE regime, pinned Session 2:

| | 1 | 2 | 4 | 6 | 8 | 12 | 16 | 1→16 |
|---|---|---|---|---|---|---|---|---|
| PSO | 1.2430 | 0.6477 | 0.3313 | 0.2263 | 0.1781 | 0.1389 | 0.1200 | **10.35×** |
| CPSO 3× | 4.2040 | 3.5457 | 3.2672 | 3.2542 | 3.1022 | 2.7313 | 2.6330 | 1.60× |
| CPSO 6× | 7.7354 | 7.2319 | 6.6838 | 6.5504 | 6.5711 | 6.1411 | 6.0779 | 1.27× |
| CPSO 12× | 13.9659 | 13.5692 | 13.1128 | 12.9992 | 12.6755 | 12.5199 | 13.0478 | 1.07× |

Completion within *T* (%): Rule A — PSO 16.3, CPSO 0.0 · Rule B — PSO 94.9, CPSO 29.6 · Rule C — PSO 94.4, CPSO 0.0.

Overrun ratio (median / p95 / mean / max): Rule A — PSO 1.14 / 1.46 / 1.15 / 1.58, CPSO 21.51 / 103.49 / 33.36 / 142.45 · Rule B — PSO 0.98 / 1.00 / 0.96 / 1.04, CPSO 1.04 / 1.79 / 1.15 / 2.39 · Rule C — PSO 0.96 / 1.00 / 0.92 / 1.07, CPSO 1.97 / 3.02 / 2.11 / 4.49.

Realised restarts at cost ratios 3 / 6 / 12: **10.39 / 25.31 / 54.69** (manuscript: 10.38 / 25.32 / 54.75) against predicted 10.80 / 26.05 / 56.00.

§4.8's figures are confirmed correct, and their estimator has been identified: *unpinned PSO 13.2% slower, unpinned CPSO 2.6% faster* is the **fixed-FE regime at sixteen threads, pooled totals** (n=60 matched cells). The manuscript does not state the regime, thread level or estimator; camera-ready should, because other defensible slices give PSO anywhere between −1.2% and +4.6%.

**Process change adopted:** the voided file is to be renamed `session2_VOID_unpinned.csv` in the release repository, with the reason recorded in the data README, so that the natural spelling cannot be loaded by accident. Analysis scripts are to assert an explicit expected file set rather than glob `session*.csv`.

Amendments D1–D8 were made before any main-matrix execution; D9 documents behaviour already exercised in Session 1 without altering it; D10 records an execution error and its resolution; D13 records an analyst error found after submission and before camera-ready.

**Process change adopted:** environment variables are set **inside** the run script rather than in the calling shell, and `mainmatrix.cpp` already echoes `OMP_PROC_BIND` / `OMP_PLACES` to the log — that echo is to be checked before leaving any session unattended. The first seven were made the same day as the original freeze, before any main-matrix execution, and are driven by calibration measurement rather than by outcomes. Original text is preserved above with strikethrough where superseded.

---

## 11b. Thermal record — Session 1 (18 Aug 2026)

Captured with HWiNFO at 10 s intervals; 1,079 samples. **Coverage: 179.7 min of the 258.8 min session (69%)** — warm-up plus roughly the first three hours. The final ~79 min is not captured; sessions 2–5 should log the full duration.

| Metric | Warm-up (0–10 min) | Steady state (>10 min) | sd | Max |
|---|---|---|---|---|
| CPU Tctl/Tdie | 67.5 °C | **68.1 °C** | 4.33 | 84.8 °C |
| CPU CCD1 (Tdie) | 64.5 °C | 66.6 °C | 4.51 | 83.0 °C |
| CPU CCD2 (Tdie) | 62.4 °C | 56.4 °C | 5.08 | 72.8 °C |
| Core clock (avg) | 3197 MHz | 3044 MHz | 749 | 4467 MHz |

Consistency check: CPU max ≥ CCD1 max ≥ CCD2 max (84.8 ≥ 83.0 ≥ 72.8), as expected given CCD1 is the hotter die.

**Throttling: none.** HTC events 0; PROCHOT (CPU) events 0 across all 1,079 samples. Thermal Limit mean 66.0%, max 79.8% of budget. Peak CPU temperature of 84.8 °C occurred at 91.7 min elapsed, with throttling flags clear at that sample. **No measurement in Session 1 was taken under thermal throttling.** This answers the principal internal-validity threat (§3.7, adversarial review B6) with measurement rather than argument.

**Warm-up protocol validated.** Equilibrium is reached inside the 600 s warm-up. Across the measurement phase, CPU mean by 20-minute bin is flat — 68.6, 69.2, 68.4, 67.6, 67.9, 68.5, 67.7, 67.1, 67.7 °C — with no monotonic trend over three hours. The steady-state sd of 4.33 °C reflects workload variation across interleaved conditions, not drift. The 600 s warm-up is sufficient and is retained unchanged.

**Finding: inter-CCD thermal asymmetry of 10.2 °C.** At steady state CCD1 averages 66.6 °C against CCD2's 56.4 °C. Two dies at materially different temperatures boost differently, so which CCD a thread is placed on affects its achievable clock. This is a concrete, measured mechanism for the pinning effect observed in pilot 3 — no benefit at 1–4 threads, 12–25% at 6–16 threads, precisely the range where the scheduler begins distributing work across both CCDs. Reported as supporting evidence for within-platform heterogeneity (RQ2) and as justification for the fixed pinning policy (D5).

**Ambient conditions.** All sessions are executed in an air-conditioned room maintained in the range **24–26 °C**. This is a controlled band rather than a per-session instrument reading; it was not logged for Session 1 and is reported as the stated room condition. The same setting applies to sessions 2–5, so ambient is held approximately constant across the experiment rather than varying between them.

*Limitation: ambient is reported as a controlled band, not a measured series, so short-term excursions within or beyond it (door opened, compressor cycling) are not captured. The die temperatures in the table above are the measured quantity on which the throttling analysis rests; ambient is reported for completeness and for comparability with the throttling literature, which conventionally states it.*

*Relevance of the 2 °C band: the edge-inference study cited in §3.7 reports throughput variation up to 27.7% across an ambient range of 22–36 °C — a 14 °C span. A 2 °C band is a small fraction of that, and Session 1 recorded zero throttling events with peak die temperature 84.8 °C against a ~95 °C limit, so ambient variation within this band is not expected to be a material source of timing variance. This is an argument from the measured die data, not an assumption about the room.*

*Correction note: two earlier drafts of this section were wrong. The first reported 70.5 min coverage; the second corrected coverage and peak CPU temperature but left CCD1 max, CCD2 max, and the sd values at their stale figures, producing an internal inconsistency in which CCD1's max sat below the overall CPU max. Root cause: the analysis filtered timestamps with a regex requiring zero-padded fields, silently discarding every sample whose second or month was a single digit (e.g. `23:40:0.740`, `18.8.2026`) — 775 of 1,079 rows. The table above is recomputed in full from the complete log. All conclusions are unchanged and now rest on 2.5× more data.*

---

## 11c. Pinning comparison — Session 2, pinned vs unpinned (19 Aug 2026)

Session 2 was first executed unpinned through operator error (D10) and subsequently re-run under the fixed policy. The two runs use **identical seeds (4–6) and an identical condition set**, differing only in `OMP_PROC_BIND`/`OMP_PLACES`. This is therefore a **controlled within-seed comparison**, not the confounded between-session comparison described in an earlier draft of this section.

**Method (stated explicitly so the figures are reproducible).** For each of the 3,120 conditions, mean `total_wall_s` is taken over the three seeds within a `(regime, rule, ratio, fid, algo, impl, threads, mb)` cell. Totals are then summed within each `(algo, threads)` group and the percentage change computed on the summed totals — a pooled ratio, not a mean of per-cell ratios, which would be biased by cells with small denominators.

**Wall-time change, unpinned relative to pinned (pooled totals):**

| Threads | PSO | SAPSO |
|---|---|---|
| 1 | −3.6% | −5.2% |
| 2 | −1.2% | −1.2% |
| 4 | −1.0% | −0.1% |
| 6 | −0.1% | −3.9% |
| 8 | −0.1% | −2.9% |
| 12 | +0.1% | −3.4% |
| 16 | **+0.3%** | **−3.1%** |

**Sensitivity to aggregation.** The direction is stable across estimators; the magnitude is not. At 16 threads PSO/SAPSO read +0.3%/−3.1% (pooled totals), +0.8%/−3.2% (median of cell ratios), and +4.4%/−6.0% (mean of cell ratios). The mean-of-ratios estimator is the most extreme and the least robust. **The pooled-totals figures in the table are the reported values**; the sensitivity is disclosed rather than suppressed.

**Direction of the effect.** From 6 threads upward, pinning and unpinning move the two algorithms in *opposite* directions: PSO is marginally slower unpinned while SAPSO is consistently ~3% faster. The proposed mechanism is that PSO parallelises across threads, so dispersed placement across two CCDs at differing temperatures (10.2 °C apart, §11b) costs it clock, whereas SAPSO's dominant cost is the serial O(m³) Cholesky refit, which benefits from the scheduler being free to migrate that single hot thread to a cooler, higher-boosting core.

**Evidential status.** The *direction* is consistent across three estimators and coheres with the pilot-3 affinity study and the measured CCD asymmetry. The *magnitudes* are small — a few percent — and materially smaller than the 12–25% reported in pilot 3, which used a different workload. The mechanism is offered as a plausible explanation supported by converging evidence, not as a demonstrated causal claim.

**Effect on verdicts: none.** Under the D11 tie rule the two runs give an **identical** inversion rate of 25.0% (n = 288 each), and FE-regime outcomes match seed-for-seed. Pinning shifts timing without changing which algorithm wins. For a study whose subject is wall-clock-budgeted comparison, that separation is worth stating: the execution policy perturbs the measured quantity while leaving the verdict intact at these magnitudes.

*Correction note: an earlier draft of this section reported a per-thread table with substantially larger magnitudes (−6.1% to +4.2%). Those figures compared Session 2 unpinned against Session 1, which differ in seed set as well as pinning policy, and were therefore confounded and not reproducible from any single stated method. The table above replaces them with the within-seed comparison that the Session 2 re-run made possible.*

---

## 11d. Interim results after four sessions — NOT FOR PUBLICATION

Project log only. Four of five sessions complete. All figures use the D11 tie rule and the D12 block-level unit of analysis.

**Headline.** Pooled across sessions 1–4: PSO wins **69.4%** of paired comparisons under FE-matched budgets and **96.0%** under time-fair Rule B. **Verdict inversion rate 26.6%.**

**Per-session stability** (condition-level rates, shown for comparison only — inference uses blocks):

| Session | Seeds | FE PSO wins | TF PSO wins | Inversions |
|---|---|---|---|---|
| S1 | 1–3 | 53.8% | 87.8% | 34.0% |
| S2 | 4–6 | 75.0% | 100.0% | 25.0% |
| S3 | 7–9 | 75.0% | 97.6% | 22.6% |
| S4 | 10–12 | 75.0% | 99.3% | 24.3% |

Sessions 2–4 agree closely; **Session 1 is the outlier** on both FE win rate and inversion rate. With three seeds per session this is within expected sampling variation, but it is worth watching in Session 5.

**Cost-ratio effect — significant under the correct test.** Block-level means (unit = problem × seed × ratio, 49 blocks per level):

| Cost ratio | Block-mean inversion rate |
|---|---|
| 3× | 23.7% |
| 6× | 27.0% |
| 12× | 29.1% |

**Friedman test across blocks: χ² = 13.04, p = 0.0015.** The trend is monotonic and significant once the pairing is respected. Note that the unpaired χ² on pooled conditions gave p = 0.23 and would have led to the wrong conclusion; that test is superseded by D12.

**Problem effect — not supported.** Earlier drafts reported a large, highly significant landscape effect. That was an artifact of pseudoreplication. Per-problem inversion rates across sessions:

| Problem | S1 | S2 | S3 | S4 |
|---|---|---|---|---|
| Rosenbrock | 11.1% | 33.3% | 30.6% | 0.0% |
| EllipsoidRotated | 26.4% | 66.7% | 26.4% | 0.0% |
| RastriginRotated | 33.3% | 0.0% | 0.0% | 66.7% |
| LunacekBiRastrigin | 76.4% | 0.0% | 33.3% | 30.6% |
| Sphere | 0.0% | — | — | — |

Ranges of 33–76 points across sessions on ~12 blocks per problem. **No claim about landscape is made.** The earlier statement that landscape dominates cost ratio is withdrawn.

**What this changes about the study's power.** Session 5 adds roughly 12 blocks, bringing the total to ≈61. That is the real sample size; the ~1,470 conditions are not independent observations. Reported intervals and tests must reflect this, and the paper should state the block structure explicitly rather than quoting condition counts.

---

## 12. Substantive observation recorded during calibration

**Gallagher101 exhibits ~20× platform-relative cost variation while four other functions do not.**

Measured single-threaded, dim 30, FE = 10,000×D, identical algorithm and parameters:

| Function | MSVC / Ryzen 9 8945HX | GCC / x86 container | Relative |
|---|---|---|---|
| Sphere | 0.292 s | ~0.30 s | ≈1× |
| EllipsoidRotated | 1.143 s | ~1.44 s | ≈1.3× |
| RastriginRotated | 1.602 s | ~1.83 s | ≈1.1× |
| Rosenbrock | 0.779 s | ~0.37 s | ≈0.5× |
| **Gallagher101** | **25.240 s** | **~1.25 s** | **≈20×** |

Supporting evidence for the mechanism: Gallagher21 (21 local optima) measured 5.518 s against Gallagher101 (101 optima) at 25.240 s — a ratio of 4.6× against a peak-count ratio of 4.8×. Cost scales with the number of stored local optima and their associated matrices, consistent with memory-bound rather than compute-bound behaviour.

**Why this matters to the study.** RQ2 asks whether platform differences perturb *rankings* rather than merely scaling *timings*. This is a direct instance: one function's cost relative to the rest changed by an order of magnitude across toolchain and platform, while four others held roughly constant. That is a non-uniform platform effect — precisely the phenomenon the reproducibility objection asserts and which has not previously been quantified for metaheuristic benchmarking.

The observation is exploratory (two platforms, one measurement each, differing compilers *and* hardware — the two are confounded here). It is reported as a motivating observation in the discussion, **not** as a controlled result. Disentangling compiler from hardware would require the second physical platform excluded in §9.
