#!/usr/bin/env python3
"""
compute_T.py - derive the cost-ratio table and the T value for each coupling
rule, from costratio.csv.

This produces two artefacts required by the methodology freeze:

  1. The PER-PROBLEM COST RATIO TABLE (SAPSO / PSO). Reported in the paper.
     The ratio is measured, not tuned: SAPSO parameters are fixed at freeze
     values. Because surrogate overhead is largely problem-independent while
     base evaluation cost is not, the ratio varies by problem, and that
     variation is itself a reported finding.

  2. The T VALUES for coupling rules A, B, C (freeze section 4.1):
       Rule B (primary)   T = SAPSO single-run time
       Rule A (robustness) T = PSO single-run time
       Rule C (control)    T = mean single-run time across the portfolio

Also reports the implied restart count k = floor(T / tau_PSO) per rule, which
is what determines whether an inversion is even possible: if k = 1 under some
rule, the baseline gets no restarts and the protocol's mechanism cannot operate.

Usage:
    python compute_T.py costratio.csv
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

df = pd.read_csv(sys.argv[1])
print(f"Loaded {len(df)} runs.")
print(f"Problems: {sorted(df.problem.unique())}")
print(f"Runs per (problem, algo): {df.groupby(['problem','algo']).size().min()}")

# ---- full precision throughout; round only for display ----
g = df.groupby(["fid", "problem", "algo"]).wall_s.agg(["mean", "std", "count"])
g = g.reset_index()

piv = g.pivot(index=["fid", "problem"], columns="algo", values="mean")
sd = g.pivot(index=["fid", "problem"], columns="algo", values="std")
n = g.pivot(index=["fid", "problem"], columns="algo", values="count")

print("\n" + "=" * 74)
print("1  COST RATIO TABLE  (SAPSO / PSO, single-threaded)")
print("=" * 74)
out = pd.DataFrame({
    "PSO_s":   piv["PSO"],
    "SAPSO_s": piv["SAPSO"],
    "ratio":   piv["SAPSO"] / piv["PSO"],
    "PSO_cv%":   sd["PSO"]   / piv["PSO"]   * 100,
    "SAPSO_cv%": sd["SAPSO"] / piv["SAPSO"] * 100,
})
print(out.to_string(float_format=lambda x: f"{x:.3f}"))
print(f"\nratio range: {out['ratio'].min():.2f}x - {out['ratio'].max():.2f}x "
      f"(spread {out['ratio'].max()/out['ratio'].min():.1f}x)")
print("Report this table. The ratio is measured at frozen SAPSO parameters,")
print("not tuned to a target.")

# ---- CI on tau, to check the freeze's 5% stability requirement ----
print("\n" + "=" * 74)
print("2  TAU PRECISION  (95% CI half-width as % of mean)")
print("=" * 74)
rows = []
for (fid, prob, algo), sub in df.groupby(["fid", "problem", "algo"]):
    t = sub.wall_s.values
    m, s, k = t.mean(), t.std(ddof=1), len(t)
    hw = stats.t.ppf(0.975, k - 1) * s / np.sqrt(k)
    rows.append(dict(fid=fid, problem=prob, algo=algo, tau_s=m,
                     ci_pct=hw / m * 100,
                     flag="UNSTABLE" if hw / m > 0.05 else "ok"))
tau = pd.DataFrame(rows)
print(tau.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
bad = (tau.flag == "UNSTABLE").sum()
if bad:
    print(f"\n{bad} configuration(s) exceed the 5% threshold set in the freeze.")
    print("Increase --runs for these before fixing T, or record an exclusion.")
else:
    print("\nAll within the 5% threshold. T may be fixed from these values.")

# ---- T per coupling rule ----
print("\n" + "=" * 74)
print("3  T VALUES PER COUPLING RULE")
print("=" * 74)
mean_all = (piv["PSO"] + piv["SAPSO"]) / 2.0
T = pd.DataFrame({
    "T_ruleB_s": piv["SAPSO"],          # primary: expensive algorithm
    "T_ruleA_s": piv["PSO"],            # robustness: cheap algorithm
    "T_ruleC_s": mean_all,              # control: portfolio mean
})
print(T.to_string(float_format=lambda x: f"{x:.4f}"))

print("\n" + "=" * 74)
print("4  IMPLIED RESTART COUNT  k = floor(T / tau_PSO)")
print("=" * 74)
K = pd.DataFrame({
    "k_ruleB": np.floor(piv["SAPSO"] / piv["PSO"]).astype(int),
    "k_ruleA": np.floor(piv["PSO"]   / piv["PSO"]).astype(int),
    "k_ruleC": np.floor(mean_all     / piv["PSO"]).astype(int),
})
print(K.to_string())
print("""
Reading this table:

  k = 1 means the baseline completes one run and gets NO restarts. Under that
  rule, on that problem, the protocol's restart mechanism cannot operate and no
  inversion is possible by construction. That is not a failed experiment - it is
  a boundary condition, and it is exactly why all three rules are reported.

  Rule A gives k = 1 everywhere by definition (T is the baseline's own run
  time). It is the control that shows what the protocol does when there is no
  restart budget.

  Under Rule B, k equals the floor of the cost ratio. Problems where the ratio
  is near 1 give k = 1 even under the primary rule - worth flagging in the
  results as the condition under which time-fairness reduces to FE-fairness.
""")

print("=" * 74)
print("Record these values in the freeze appendix. They are fixed from this")
print("point and are not recomputed from main-matrix runs.")
print("=" * 74)
