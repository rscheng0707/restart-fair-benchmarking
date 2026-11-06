# Supplementary materials — *Internal Parallelism Breaks Wall-Clock Budget Invariance in Restart-Fair Metaheuristic Benchmarking*

Code, data, thermal logs, and the full pre-registration for the ICMLSC 2027 paper (ID M4338).

> **Citation.** *(fill in once the proceedings reference is assigned)*
> **Archived release.** *[![DOI](https://zenodo.org/badge/1448300647.svg)](https://doi.org/10.5281/zenodo.24015095)*

Everything in the paper regenerates from this repository with one command:

```bash
python analysis/reproduce_all.py
```

It prints all seven tables and every headline figure, and asserts the condition
counts (15,600 → 1,800 → 1,488 → 62 blocks) along the way. If an assertion
fails, the data and the paper have diverged.

---

## ⚠ Two things to know before loading any data

### 1. One session file is void and must not be analysed

`data/` contains **exactly five canonical session files**, `session1.csv` … `session5.csv`.

A sixth exists at `data/voided/session2_VOID_unpinned.csv`. **It is not part of the analysis.** Session 2 was first executed with thread pinning accidentally disabled (`OMP_PROC_BIND` and `OMP_PLACES` unset), which confounded placement policy with seed set. The session was re-run under the correct policy, and the re-run is what `data/session2.csv` holds. The voided run is published because it supports the pinning robustness check in §4.8 of the paper, and because discarding inconvenient data is not how this study was conducted.

It lives in a subdirectory precisely so `data/session*.csv` cannot pick it up. **Do not glob across `data/` recursively.** Deviation D10 in `PREREGISTRATION.md` records the original error; **D13 records that three tables in the submitted manuscript were computed from this voided file**, and gives the corrected values. `analysis/reproduce_all.py` asserts an explicit file list rather than globbing, for exactly this reason.

### 2. `tau.csv` has a non-obvious lookup rule

`data/tau.csv` is keyed on `(fid, ratio, threads, omp)`. The trap: **PSO's own τ is the `ratio = 0` row**, because PSO does not use the cost-load parameter, while **CPSO's τ is the row for its own ratio.** Under coupling rule B, `T_s` equals CPSO's τ for *both* algorithms in a pair, since the budget is set by the costed competitor.

A naive merge on `ratio` reconciles only half the rows and looks like corrupt data. Done correctly, all 3,600 time-fair rows reconcile exactly — `reproduce_all.py` asserts this at the end of its run.

---

## Contents

| Path | What it is |
|---|---|
| `data/session1.csv` … `session5.csv` | Main matrix, 3,120 conditions each, 15,600 total. One row per measured condition. |
| `data/voided/session2_VOID_unpinned.csv` | The voided unpinned session 2. Not analysed. |
| `data/tau.csv` | The 160-entry τ pilot table. Sets *T* and *k* for every configuration. See the lookup rule above. |
| `data/tau_recalibrated.csv` | Five τ entries re-measured after the Sphere evaluation-budget change (deviations D6/D8). |
| `data/calibration.csv` | Cost-load calibration: the `archive_cap` value chosen per problem and target ratio. |
| `data/costratio_verification.csv` | Achieved cost ratios, confirming all fifteen points land within 4.0% of target. |
| `thermal/session*_hwinfo.csv` | Raw HWiNFO64 logs, 10-second interval. `session1` covers 179.7 minutes / 1,079 samples and is the one analysed in the threats-to-validity section. |
| `code/` | The C++ harness and the Python helpers. See below. |
| `analysis/reproduce_all.py` | Regenerates every table and figure in the paper. |
| `PREREGISTRATION.md` | The dated methodology freeze with all 13 amendments. |

### `code/`

| File | Role |
|---|---|
| `mainmatrix.cpp` | The main experiment. Enumerates the crossed factor space and emits one `session*.csv`. |
| `taupilot.cpp` | The τ pilot. Produces `tau.csv`. Run before the main matrix; never reuses analysed runs. |
| `calibrate.cpp` | Chooses `archive_cap` per problem so the cost ratio hits 3×, 6×, 12×. |
| `costratio.cpp` | Verifies the achieved ratios against target. |
| `pilot3.cpp` | Timing pilot used to size the budget and check measurement precision. |
| `pso.hpp` | Canonical inertia-weight PSO. Per-particle RNG streams via a splitmix64 mixer, so serial and OpenMP results are bit-identical at one thread. |
| `sapso.hpp` | CPSO — the costed competitor. **See the naming note below.** |
| `bbob_adapter.hpp` | Wraps IOHexperimenter's BBOB problems. One problem instance per thread (they are stateful and not thread-safe). |
| `benchmarks.hpp`, `harness.hpp` | Shared problem definitions and the run harness. |
| `test_equivalence.cpp` | Asserts serial/OpenMP bit-identity at one thread. |
| `compute_T.py`, `pick_caps.py` | Derive *T* from `tau.csv`; select calibration caps. |
| `run_session.bat`, `Makefile`, `BUILD_IOH.md` | Build and run. **Thread placement is set inside `run_session.bat`, not in the calling shell** — that is the process change adopted after D10. |

---

## Data dictionary — `data/session*.csv`

| Column | Meaning |
|---|---|
| `session` | Session index 1–5. Sessions partition the 15 seeds into blocks of three. |
| `regime` | `fe` (evaluation-matched), `timefair`, `multife` (multi-budget at 0.1× and 1×). |
| `rule` | Coupling rule `A` (*T* = τ_PSO), `B` (*T* = τ_CPSO, primary), `C` (midpoint). Empty outside `timefair`. |
| `ratio` | Target per-iteration cost ratio: 3, 6, 12. `0` where the factor does not apply — **not** a ratio of zero. |
| `fid` | BBOB function id: 1 Sphere, 8 Rosenbrock, 10 EllipsoidRotated, 15 RastriginRotated, 24 LunacekBiRastrigin. |
| `problem` | Human-readable function name. |
| `algo` | `PSO`, or `SAPSO` — **the on-disk label for CPSO.** See below. |
| `impl` | `serial` or `omp`. |
| `threads` | 1, 2, 4, 6, 8, 12, 16. Always 1 for `serial`. |
| `mb` | Multi-budget multiplier (0.1 or 1) under `multife`; `1` elsewhere. |
| `seed` | 1–15. |
| `archive_cap` | Calibrated cost-load size *m* for this problem and ratio. Changes cost only — see the inertness note. |
| `fe_budget` | Evaluation budget: 10,000·D, or 30,000·D for Sphere (deviation D6). |
| `T_s` | Nominal wall-clock budget *T* in seconds. Empty outside `timefair`. |
| `tau_s` | Pilot-estimated mean single-run duration, from `tau.csv`. Never estimated from the analysed runs. |
| `k_predicted` | ⌊*T*/τ⌋ as the protocol specifies. |
| `runs_completed` | Restarts actually finished. Differs from `k_predicted` by 2.3–3.8%. |
| `best_fitness` | Best objective value returned. **The outcome variable.** |
| `total_wall_s` | Total wall-clock time for the condition. **The timing variable.** |
| `total_fes` | Function evaluations consumed. |
| `order_idx` | Randomised execution order within the session, recorded for thermal auditing. |

### `algo == "SAPSO"` means CPSO, and its cost load is behaviourally inert

The costed competitor was built as a surrogate-assisted PSO, and the label persists through the code and the CSVs. Post-collection inspection established that the surrogate never functions: the multiquadric kernel matrix is conditionally negative definite, Cholesky factorisation always fails, the model is never consulted, and candidate selection falls back to a fixed rule. `archive_cap` therefore manipulates **computational cost only**.

This is directly checkable in the released data. Across all 600 configurations matched by problem, implementation, thread count and seed, `best_fitness` is bit-identical at all three cost-ratio levels while `total_wall_s` roughly triples:

```python
import glob, pandas as pd
df = pd.concat(pd.read_csv(f) for f in sorted(glob.glob("data/session?.csv")))
g = df[(df.regime == "fe") & (df.algo == "SAPSO")]
g.groupby("ratio").total_wall_s.mean()   # 3.37 → 6.84 → 13.23
```

The paper therefore reports the algorithm as **CPSO (Costed PSO)** and restricts its claims accordingly: it tests how a wall-clock budget behaves when two otherwise identical algorithms differ in per-iteration cost, and makes no claim about genuine surrogate-assisted optimisers. The rename was not propagated into the data files, so that published figures reproduce from them without a translation step.

---

## Key results

| Quantity | Value |
|---|---|
| PSO wins, evaluation-matched | 69.4% |
| PSO wins, nominal-*T* rule B | 96.4% |
| Pooled inversion rate | 27.0% (n = 62 blocks) |
| By cost ratio | 24.2% / 27.4% / 29.4% — Friedman χ² = 16.75, *p* = 0.00023 |
| Blocks constant across the 4× cost-ratio span | 53 of 62; the 9 that move all rise (sign test *p* = 0.0039) |
| PSO 1→16 thread speedup | 10.35× |
| Mean predicted *k*, rule B, 1→16 threads | 6.5 → 66.2 (**10.18×**) |
| Maximum overrun of nominal *T* (rule A, CPSO) | 142.45× |

Inference uses the **problem × seed** block as the unit, not the condition — verdicts are constant within a block in 53 of 62 blocks, so condition-level tests are pseudoreplicated. Deviation D12 records this and the two conclusions it reversed.

---

## Build environment

Ryzen 9 8945HX (16 physical cores, SMT off for measurement), Windows 11, MSVC with `/openmp:llvm`. Thread placement fixed at `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, **set inside the run script rather than the calling shell** (the process change adopted after D10).

Two IOHexperimenter integration notes that cost real time:

- BBOB problem objects are **stateful and not thread-safe**. Each thread needs its own instance. `bbob_adapter.hpp` handles this.
- **Do not include the umbrella `ioh.hpp` header** — the CEC suites abort at static initialisation. Include the BBOB headers directly.

See `code/BUILD_IOH.md` for the full build procedure. IOHexperimenter is not vendored here; it is fetched at build time.

---

## Minting the DOI

1. Push to GitHub and tag a release, e.g. `v1.0-icmlsc2027`.
2. Enable the repository in Zenodo (Account → GitHub → toggle the repo on) **before** tagging. Zenodo only captures releases created after the toggle, so tagging first means re-tagging.
3. Collect the DOI badge Zenodo issues and paste it at the top of this file.
4. Put the **concept DOI** — the one that always resolves to the newest version — in the paper, not the version DOI. If the data is ever corrected, the paper's link stays valid.
