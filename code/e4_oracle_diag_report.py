#!/usr/bin/env python3
"""EXPLORATORY report for the E4 oracle-arm diagnosis (spec E4_ORACLE_DIAGNOSIS_SPEC.md). Usage: e4_oracle_diag_report.py --in DIR --out DIR
Reads e4diag_{replay,fresh}[_perj]_cell{0..3}.csv produced by e4_oracle_diag.py and results/raw_amend3.csv (replay check)."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]; IN = args[args.index("--in") + 1]; OUT = args[args.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
NAMES = {0: "ring_128", 1: "ring_256", 2: "ring_512", 3: "ring_1024"}
LINES = []
def P(s=""): print(s); LINES.append(str(s))
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def f(k, n): lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)

P("EXPLORATORY E4 oracle-arm diagnosis (no gate). Denominator for coverage = non-abstained datasets; abstentions-as-misses shown beside it.")
# ---- replay check
r3 = pd.read_csv(os.path.join(HERE, "results", "raw_amend3.csv"), low_memory=False); r3 = r3[r3.est == "E4o"]
rows = []; summ = []
for c in range(4):
    d = pd.read_csv(os.path.join(IN, f"e4diag_replay_cell{c}.csv")); o = r3[r3.cell == c]
    m = d.merge(o, on=["cell", "ds"], suffixes=("", "_o"), how="outer", indicator=True)
    both = m[m._merge == "both"]
    ab_ok = int((both.abstain == both.abstain_o).sum())
    a = both[(both.abstain == 0) & (both.abstain_o == 0)]
    hw_ok = int((np.abs(a.hw - a.hw_o) < 1e-6).sum()); cov_ok = int((a.cov_simult == a.cov_simult_o).sum()); J_ok = int((a.J == a.J_o).sum())
    P(f"replay check {NAMES[c]}: rows {len(both)}/200; abstain equal {ab_ok}; answered {len(a)}: J equal {J_ok}, hw within 1e-6 {hw_ok}, cov_simult equal {cov_ok}")
    summ.append(dict(cell=NAMES[c], rows=len(both), abstain_equal=ab_ok, answered=len(a), J_equal=J_ok, hw_equal=hw_ok, cov_equal=cov_ok))
pd.DataFrame(summ).to_csv(os.path.join(OUT, "e4diag_replay_check.csv"), index=False)

def table(arm):
    out = []
    for c in range(4):
        d = pd.read_csv(os.path.join(IN, f"e4diag_{arm}_cell{c}.csv")); n = len(d); a = d[d.abstain == 0]; na = len(a)
        k = int(a.cov_simult.sum()); nm = int(a.miss.sum()); nb = int((a.miss_kind == "bound").sum()); nc_ = int((a.miss_kind == "calibration").sum())
        pj = pd.read_csv(os.path.join(IN, f"e4diag_{arm}_perj_cell{c}.csv"))
        r = dict(arm=arm, cell=NAMES[c], n=n, abstain=n - na, answered=na, covered=k, cov_cond=k / na, cov_lo=cp(k, na)[0], cov_hi=cp(k, na)[1], cov_strict=k / n,
                 misses=nm, miss_bound=nb, miss_calibration=nc_, datasets_with_any_bound_violation=int(a.bound_viol.sum()),
                 max_b_over_B=float(a.max_b_over_B.max()),
                 noise_cover=float(a.noise_ok.mean()), noise_cover_k=int(a.noise_ok.sum()),
                 mean_crit=float(a.crit.mean()), q95_tmax_noise=float(a.tmax_noise.quantile(.95)))
        # scale / centring from replicates (per-j, per J stratum with >=100 datasets): sd_emp over datasets vs bootstrap sd
        res = []
        for J, g in pj.groupby("J"):
            ids = g.ds.nunique()
            if ids < 100: continue
            for j, h in g.groupby("j"):
                res.append(dict(J=J, j=j, ids=ids, sd_ratio=h.sd.median() / h.n.std(ddof=1), mean_n_over_sd=h.n.mean() / h.n.std(ddof=1), b_over_B=(h.b.abs() / h.B).max()))
        res = pd.DataFrame(res)
        if len(res):
            res.assign(cell=NAMES[c], arm=arm).to_csv(os.path.join(OUT, f"e4diag_{arm}_{NAMES[c]}_perj_strata.csv"), index=False)
            r["median_sd_boot_over_sd_emp"] = float(res.sd_ratio.median()); r["min_sd_ratio"] = float(res.sd_ratio.min())
            r["mean_abs_centering_in_sd"] = float(res.mean_n_over_sd.abs().mean()); r["max_abs_centering_in_sd"] = float(res.mean_n_over_sd.abs().max())
        out.append(r)
    return pd.DataFrame(out)

T = pd.concat([table("replay"), table("fresh")])
T.to_csv(os.path.join(OUT, "e4diag_summary.csv"), index=False)
for r in T.itertuples():
    P(f"[{r.arm}] {r.cell}: coverage {f(r.covered, r.answered)}; abstained {r.abstain}/{r.n} (strict {100*r.cov_strict:.1f}%); misses {r.misses} = bound {r.miss_bound} + calibration {r.miss_calibration}; "
      f"datasets with any |b_j|>B_j: {r.datasets_with_any_bound_violation}; max |b|/B {r.max_b_over_B:.3f}; stochastic-part coverage P(max|n|/sd<=c) {100*r.noise_cover:.1f}%; "
      f"median sd_boot/sd_emp {getattr(r, 'median_sd_boot_over_sd_emp', float('nan')):.3f} (min {getattr(r, 'min_sd_ratio', float('nan')):.3f}); mean |centring| {getattr(r, 'mean_abs_centering_in_sd', float('nan')):.3f} sd")
# per-dataset attribution of the original 200 (replay)
parts = []
for c in range(4):
    d = pd.read_csv(os.path.join(IN, f"e4diag_replay_cell{c}.csv")); d["cell_name"] = NAMES[c]; parts.append(d[(d.abstain == 0) & (d.miss == 1)])
mm = pd.concat(parts)[["cell", "cell_name", "ds", "J", "M", "crit", "tmax_noise", "bound_viol", "max_b_over_B", "miss_kind", "miss_also_noise_excess"]]
mm.to_csv(os.path.join(OUT, "e4diag_replay_each_miss.csv"), index=False)
P(f"each original miss attributed: {len(mm)} misses -> " + str(mm.miss_kind.value_counts().to_dict()))
open(os.path.join(OUT, "e4diag_report.txt"), "w").write("\n".join(LINES))
