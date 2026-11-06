# Building pilot3 with IOHexperimenter

## 1. Get IOHexperimenter

The repo's submodule is configured for SSH, so `--recursive` fails on a plain
clone. Clone the dependency separately over HTTPS:

```
git clone --depth 1 https://github.com/IOHprofiler/IOHexperimenter.git ioh
cd ioh
git clone --depth 1 https://github.com/fmtlib/fmt.git external/fmt
cd ..
```

No CMake build or install is needed — pilot3 uses the headers directly.

## 2. Build (GCC / WSL2)

```
g++ -O2 -fopenmp -std=c++17 -I. \
    -Iioh/include -Iioh/external/fmt/include -DFMT_HEADER_ONLY \
    pilot3.cpp -o pilot3
```

## 3. Build (MSVC, matching the pilot 1 / pilot 2 stack)

```
cl /O2 /openmp:llvm /std:c++17 /EHsc /I. ^
   /Iioh\include /Iioh\external\fmt\include /DFMT_HEADER_ONLY ^
   pilot3.cpp /Fe:pilot3.exe
```

Keep the toolchain the same as pilots 1 and 2 for the diagnostic comparison to
hold. If you build with both, label the passes distinctly — compiler is an
experimental factor, not an implementation detail.

## Two integration notes

**Do not include `ioh.hpp`.** The umbrella header pulls in the CEC suites, whose
problem classes register at static-initialisation time and load external `.dat`
resource files; with `IOH_RESOURCE_DIR` unset the process aborts before `main()`
runs. `bbob_adapter.hpp` includes only `ioh/problem/bbob.hpp` and constructs
problem classes directly, avoiding the registry entirely.

**BBOB problem objects are stateful and not thread-safe.** Each instance holds
evaluation counters and transformation buffers mutated on every call. Sharing
one across OpenMP threads is a data race that will usually appear to work while
silently corrupting FE accounting. `BBOBProblem` allocates one instance per
thread; all are constructed with identical (fid, instance, dim), so they are
numerically identical and differ only in internal counters.

## 4. Three affinity passes

```
OMP_PROC_BIND=close OMP_PLACES=cores   ./pilot3 --label pinned_cores  > p3_pinned.csv   2> p3_pinned.log
OMP_PROC_BIND=false                    ./pilot3 --label unpinned      > p3_unpinned.csv 2> p3_unpinned.log
OMP_PROC_BIND=close OMP_PLACES=threads ./pilot3 --label pinned_smt --smt > p3_smt.csv   2> p3_smt.log

python analyze_pilot2.py p3_pinned.csv p3_unpinned.csv p3_smt.csv
```

The analysis script groups by `label`, so all three can be passed together.

## Expected run length

Measured single-thread, dim 30, 300k FE:

| fid | function          | time   |
|-----|-------------------|--------|
| 1   | Sphere            | 0.30 s |
| 8   | Rosenbrock        | 0.37 s |
| 10  | EllipsoidRotated  | 1.40 s |
| 15  | RastriginRotated  | 1.74 s |
| 21  | Gallagher101      | 2.27 s |

The rotated and multimodal functions reach the 1–10 s target that pilot 2 missed
by an order of magnitude. Sphere and Rosenbrock remain short; if their tau CI
stays above 5%, raise `--femult` for those two rather than for the whole suite.

The 7.5x per-evaluation cost spread across problems is itself useful — it is the
cost disparity the study is about, now present within the benchmark suite rather
than only between algorithms.
