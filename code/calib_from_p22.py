"""Window-sufficiency calibration from p22 (PLAN Amendment 2, question c): for each estimator, the span log10(t_last/t_first) of the data window at the
minimum detectable L (power>=0.8) on the realised crossover families. Output results/causet_p22_calibration.csv"""
import numpy as np, pandas as pd
import p22lib as L
r = pd.read_csv("results/raw.csv"); a = pd.read_csv("results/raw_amend.csv"); m = pd.read_csv("results/cells_meta.csv")
r["src"] = "orig"; a["src"] = "amend"
d = pd.concat([r, a], ignore_index=True)
d = d.merge(m[["id", "fam", "L", "N", "grid"]], left_on="cell", right_on="id")
mins = {("E0", "aniso2_w16"): 192, ("E0", "aniso3_w8"): 96, ("E1", "aniso2_w16"): 256, ("E1", "aniso3_w8"): 96,
        ("E1x", "aniso2_w16"): 192, ("E1x", "aniso3_w8"): 96, ("E2", "aniso2_w16"): 192, ("E2", "aniso3_w8"): 128,
        ("E3", "aniso2_w16"): 128, ("E3", "aniso3_w8"): 48, ("E1q", "aniso2_w16"): 256, ("E1q", "aniso3_w8"): 96}
rows = []
for (est, fam), Lmin in mins.items():
    x = d[(d.est == est) & (d.fam == fam) & (d.L == Lmin) & (d.abstain == 0)]
    if not len(x): continue
    N = int(x.N.iloc[0])
    if est == "E3":
        span = np.nan; J = np.nan
    else:
        g = L.make_grid(24, 34, 8192)
        J = x.J.dropna()
        # E2's J counts secants; its last time is grid[J+1]; others grid[J-1]
        last = np.array([g[int(j) + 1] if est == "E2" else g[int(j) - 1] for j in J])
        span = float(np.median(np.log10(last / 24.0)))
    rows.append(dict(est=est, fam=fam, minL=Lmin, N=N, median_J=float(x.J.median()) if est != "E3" else np.nan, span_decades=span, n_datasets=len(x)))
out = pd.DataFrame(rows); out.to_csv("results/causet_p22_calibration.csv", index=False); print(out.to_string())
