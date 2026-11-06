#!/usr/bin/env python3
"""
Regenerate every table and headline figure in the paper from the released data.

Usage:  python analysis/reproduce_all.py            (run from the repository root)

Asserts an explicit file set rather than globbing. This is deliberate: the
repository contains a voided session (data/voided/) that must never enter the
analysis, and a recursive glob is exactly how it got into three tables of the
submitted manuscript. See deviation D13 in PREREGISTRATION.md.
"""
import sys
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SESSIONS = [os.path.join(ROOT, "data", f"session{i}.csv") for i in range(1, 6)]
TAU = os.path.join(ROOT, "data", "tau.csv")

TIE_RTOL, TIE_ATOL = 1e-9, 1e-12       # deviation D11
PAIR_KEY = ["fid", "impl", "threads", "seed", "ratio"]
PSO_KEY = ["fid", "impl", "threads", "seed"]


def load():
    missing = [p for p in SESSIONS + [TAU] if not os.path.exists(p)]
    if missing:
        sys.exit("missing required file(s):\n  " + "\n  ".join(missing))
    df = pd.concat([pd.read_csv(p) for p in SESSIONS], ignore_index=True)
    assert len(df) == 15600, f"expected 15,600 measured conditions, got {len(df):,}"
    return df, pd.read_csv(TAU)


def build_pairs(df):
    """Pair each CPSO condition with its PSO counterpart under both budget regimes."""
    fe = df[df.regime == "fe"]
    tf = df[(df.regime == "timefair") & (df.rule == "B")]

    pairs = (fe[fe.algo == "SAPSO"][PAIR_KEY + ["session", "best_fitness"]]
             .rename(columns={"best_fitness": "fe_cpso"})
             .merge(fe[fe.algo == "PSO"][PSO_KEY + ["best_fitness"]]
                    .rename(columns={"best_fitness": "fe_pso"}), on=PSO_KEY))
    idx = pd.MultiIndex.from_frame(pairs[PAIR_KEY])
    pairs["tf_pso"] = tf[tf.algo == "PSO"].set_index(PAIR_KEY).best_fitness.reindex(idx).values
    pairs["tf_cpso"] = tf[tf.algo == "SAPSO"].set_index(PAIR_KEY).best_fitness.reindex(idx).values
    pairs = pairs.dropna(subset=["tf_pso", "tf_cpso"]).copy()
    assert len(pairs) == 1800, f"expected 1,800 candidate pairs, got {len(pairs):,}"

    # A tie expresses no verdict, so it cannot be reversed (D11).
    tie = (np.isclose(pairs.fe_pso, pairs.fe_cpso, rtol=TIE_RTOL, atol=TIE_ATOL) |
           np.isclose(pairs.tf_pso, pairs.tf_cpso, rtol=TIE_RTOL, atol=TIE_ATOL))
    pairs["tie"] = tie
    a = pairs[~tie].copy()
    assert len(a) == 1488, f"expected 1,488 analysed comparisons, got {len(a):,}"

    a["fe_pso_wins"] = a.fe_pso < a.fe_cpso
    a["tf_pso_wins"] = a.tf_pso < a.tf_cpso
    a["inv"] = a.fe_pso_wins != a.tf_pso_wins
    return pairs, a


def h(title):
    print(f"\n{title}\n{'-' * len(title)}")


def main():
    df, tau = load()
    pairs, a = build_pairs(df)

    h("Condition accounting (Figure 1)")
    print(f"measured conditions      {len(df):>7,}")
    print(f"candidate paired         {len(pairs):>7,}")
    per_sess = pairs[pairs.tie].groupby("session").size()
    print(f"ties excluded (D11)      {pairs.tie.sum():>7,}   per session: "
          + " / ".join(str(int(n)) for n in per_sess.values))
    print(f"analysed comparisons     {len(a):>7,}")

    h("Table 1 — PSO wall-time scaling, evaluation-matched regime")
    p = df[(df.regime == "fe") & (df.algo == "PSO")].groupby("threads").total_wall_s.mean()
    base = p.loc[1]
    print(f"{'thr':>4} {'wall (s)':>10} {'speedup':>9} {'efficiency':>11}")
    for t in p.index:
        su = base / p.loc[t]
        print(f"{t:>4} {p.loc[t]:>10.4f} {su:>8.2f}x {su / t * 100:>10.1f}%")

    h("Table 2 — mean predicted k = floor(T/tau), rule B")
    tfb = df[(df.regime == "timefair") & (df.rule == "B")]
    kp = tfb[tfb.algo == "PSO"].groupby("threads").k_predicted.mean()
    kc = tfb[tfb.algo == "SAPSO"].groupby("threads").k_predicted.mean()
    print(f"{'thr':>4} {'PSO':>8} {'CPSO':>8}")
    for t in kp.index:
        print(f"{t:>4} {kp.loc[t]:>8.1f} {kc.loc[t]:>8.1f}")
    print(f"\nHEADLINE: allocated search changes {kp.loc[16] / kp.loc[1]:.2f}x "
          f"from one to sixteen threads while CPSO stays at {kc.loc[16]:.1f}")

    h("Table 3 — inversion rate by session")
    for s, g in a.groupby("session"):
        print(f"session {int(s)}   {g.inv.mean() * 100:>5.1f}%   n={len(g):>4}")

    h("Table 4 — inversion rate by cost ratio (block-level, D12)")
    blocks = a.groupby(["fid", "seed"])
    full = blocks.ratio.nunique()
    keep = set(full[full == 3].index)
    ab = a[a.set_index(["fid", "seed"]).index.isin(keep)]
    bl = ab.groupby(["fid", "seed", "ratio"]).inv.mean().unstack("ratio")
    print(f"{'ratio':>6} {'rate':>8}")
    for r in bl.columns:
        print(f"{r:>5}x {ab[ab.ratio == r].inv.mean() * 100:>7.1f}%")
    print(f"\nblocks contributing at every ratio: {len(bl)}")
    try:
        from scipy.stats import friedmanchisquare, binomtest
        st, pv = friedmanchisquare(*[bl[c].values for c in bl.columns])
        print(f"Friedman chi2 = {st:.2f}, p = {pv:.5f}")
        const = (bl.nunique(axis=1) == 1).sum()
        moving = bl[bl.nunique(axis=1) > 1]
        up = (moving[bl.columns[-1]] > moving[bl.columns[0]]).sum()
        print(f"constant across the 4x cost-ratio span: {const} of {len(bl)}")
        print(f"of the {len(moving)} that move, {up} rise with cost ratio "
              f"(sign test p = {binomtest(up, len(moving)).pvalue:.4f})")
    except ImportError:
        print("scipy not installed — skipping Friedman and sign test")

    h("Table 5 — runs completing within nominal T, by coupling rule")
    tfa = df[df.regime == "timefair"].copy()
    tfa["ok"] = tfa.total_wall_s <= tfa.T_s
    print(tfa.groupby(["rule", "algo"]).ok.mean().mul(100).round(1).unstack().to_string())

    h("Table 6 — realised wall time as a ratio to nominal T")
    tfa["ov"] = tfa.total_wall_s / tfa.T_s
    g = tfa.groupby(["rule", "algo"]).ov.agg(
        median="median", p95=lambda x: x.quantile(.95), mean="mean", max="max")
    print(g.round(2).to_string())

    h("Table 7 — inversion rate by coupling rule")
    for rule in ["A", "B", "C"]:
        tfr = df[(df.regime == "timefair") & (df.rule == rule)]
        idx = pd.MultiIndex.from_frame(a[PAIR_KEY])
        tp = tfr[tfr.algo == "PSO"].set_index(PAIR_KEY).best_fitness.reindex(idx).values
        tc = tfr[tfr.algo == "SAPSO"].set_index(PAIR_KEY).best_fitness.reindex(idx).values
        m = pd.DataFrame({"fe_w": a.fe_pso_wins.values, "tp": tp, "tc": tc}).dropna()
        print(f"rule {rule}   {((m.tp < m.tc) != m.fe_w).mean() * 100:>5.1f}%   n={len(m):>5}")

    h("Headline figures")
    print(f"PSO wins, evaluation-matched   {a.fe_pso_wins.mean() * 100:.1f}%")
    print(f"PSO wins, nominal-T rule B     {a.tf_pso_wins.mean() * 100:.1f}%")
    print(f"verdicts reversed              {a.inv.mean() * 100:.1f}%")
    print(f"max overrun (rule A, CPSO)     {g.loc[('A', 'SAPSO'), 'max']:.2f}x")

    h("tau.csv reconciliation")
    t = tfb.copy()
    t["omp"] = (t.impl == "omp").astype(int)
    # NOTE: PSO's own tau is the ratio=0 row; CPSO's is the row for its own ratio.
    tot = ok_n = 0
    for algo, use_own in [("PSO", False), ("SAPSO", True)]:
        s = t[t.algo == algo].copy()
        s["lr"] = s.ratio if use_own else 0
        m = s.merge(tau.rename(columns={"ratio": "lr", "tau_s": "ref"}),
                    on=["fid", "lr", "threads", "omp"], how="left")
        n = int(np.isclose(m.tau_s, m.ref, rtol=1e-9).sum())
        print(f"{algo:<6} tau_s matches tau.csv: {n}/{len(m)}")
        tot += len(m); ok_n += n
    assert ok_n == tot, "tau.csv does not reconcile with the session files"
    print(f"all {tot:,} rows reconcile")

    print("\nAll assertions passed.")


if __name__ == "__main__":
    main()
