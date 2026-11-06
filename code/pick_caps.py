#!/usr/bin/env python3
"""
pick_caps.py - select archive_cap per (problem, target ratio) from calibrate.csv

Produces the parameter table required by the methodology freeze (section 4.4):
for each problem and each target cost ratio level (2x, 5x, 10x), the archive_cap
that comes closest, together with the ACHIEVED ratio.

Achieved ratios are reported alongside targets because discrete archive_cap
values will not land exactly on the target. The paper reports what was obtained.

Where a target is unreachable within the swept range, that is reported as
out of range rather than silently clamped to the nearest endpoint - a target
the design cannot reach is a fact about the design, not a rounding decision.

Usage:
    python pick_caps.py calibrate.csv [--targets 2,5,10]
"""
import sys
import numpy as np
import pandas as pd

TARGETS = [2.0, 5.0, 10.0]

args = [a for a in sys.argv[1:]]
if not args:
    print(__doc__)
    sys.exit(1)
path = args[0]
if "--targets" in args:
    TARGETS = [float(x) for x in args[args.index("--targets") + 1].split(",")]

df = pd.read_csv(path)

# Concatenating per-problem sweeps (e.g. Windows `copy /b a.csv + b.csv`) leaves
# a header row inside the data for every file after the first. Drop any row
# whose numeric columns did not parse, then coerce types.
before = len(df)
df = df[df.archive_cap.astype(str).str.strip().str.isdigit()].copy()
for c in ["fid", "dim", "swarm", "archive_cap", "n_cand", "refit_every", "n_runs"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
for c in ["pso_tau_s", "sapso_tau_s", "ratio", "pso_cv", "sapso_cv"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna(subset=["ratio", "archive_cap"])
if before != len(df):
    print(f"Dropped {before - len(df)} embedded header/malformed row(s) "
          f"from concatenated input.")

print(f"Loaded {len(df)} sweep points.")
print(f"Problems: {sorted(df.problem.unique())}")
print(f"archive_cap values swept: {sorted(df.archive_cap.unique())}")

print("\n" + "=" * 78)
print("1  RATIO vs archive_cap  (measured)")
print("=" * 78)
piv = df.pivot_table(index=["fid", "problem"], columns="archive_cap", values="ratio")
print(piv.to_string(float_format=lambda x: f"{x:.2f}"))

# CV sanity: high variance at a sweep point makes its ratio untrustworthy.
print("\nSAPSO CV% by sweep point (values above ~5% warrant more runs):")
cv = df.pivot_table(index=["fid", "problem"], columns="archive_cap", values="sapso_cv")
print(cv.to_string(float_format=lambda x: f"{x:.1f}"))

print("\n" + "=" * 78)
print("2  SELECTED archive_cap PER TARGET RATIO")
print("=" * 78)
rows = []
for (fid, prob), g in df.groupby(["fid", "problem"]):
    g = g.sort_values("archive_cap")
    lo, hi = g.ratio.min(), g.ratio.max()
    for t in TARGETS:
        if t < lo or t > hi:
            rows.append(dict(fid=fid, problem=prob, target=t, archive_cap=None,
                             achieved=None, err_pct=None,
                             note=f"OUT OF RANGE (swept {lo:.2f}x-{hi:.2f}x)"))
            continue
        i = (g.ratio - t).abs().idxmin()
        r = g.loc[i]
        err = (r.ratio - t) / t * 100
        note = ""
        if abs(err) > 20:
            note = "poor match - widen the cap ladder near this point"
        elif r.sapso_cv > 5:
            note = "high CV - increase --runs for this point"
        rows.append(dict(fid=fid, problem=prob, target=t,
                         archive_cap=int(r.archive_cap), achieved=r.ratio,
                         err_pct=err, note=note))
sel = pd.DataFrame(rows)
print(sel.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

print("\n" + "=" * 78)
print("3  IMPLIED RESTART COUNT  k = floor(achieved ratio)   [Rule B]")
print("=" * 78)
ok = sel.dropna(subset=["achieved"]).copy()
ok["k"] = np.floor(ok.achieved).astype(int)
k = ok.pivot(index=["fid", "problem"], columns="target", values="k")
print(k.to_string())
print("""
Under Rule B, T is SAPSO's single-run time, so k = floor(cost ratio): the
number of independent PSO runs that fit in the same budget. This is the whole
point of the factor. At the 2x level the baseline gets 2 restarts; at 10x it
gets 10. If any cell reads 1, the mechanism cannot operate there and the cap
ladder needs extending for that problem.
""")

miss = sel[sel.archive_cap.isna()]
if len(miss):
    print("=" * 78)
    print("UNREACHABLE TARGETS - extend the cap ladder and re-run these:")
    print(miss[["problem", "target", "note"]].to_string(index=False))
    caps = sorted(df.archive_cap.unique())
    print(f"\nSuggestion: --caps {','.join(str(c) for c in caps)},"
          f"{caps[-1]*2},{caps[-1]*3}")
    print("Note the O(m^3) fit cost: doubling archive_cap multiplies the")
    print("surrogate term by roughly 8. Use --maxsec to bound the sweep.")
else:
    print("All targets reachable within the swept range.")

print("\n" + "=" * 78)
print("Record the selected archive_cap values AND the achieved ratios in the")
print("freeze appendix (section 4.4). They are fixed from this point.")
print("=" * 78)
