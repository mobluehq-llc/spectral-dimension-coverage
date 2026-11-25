import pandas as pd, numpy as np
from scipy.stats import beta
def cp(k, n): return (beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(0.975, k + 1, n - k) if k < n else 1.0)
R = "results/"
ex = pd.read_csv(R + "mindrop_exact.csv"); print("== exact q (B1)")
print(ex[["fam", "L", "eps", "exact_drop", "D_lo", "D_hi", "status", "steps", "n_incon", "seconds"]].sort_values(["fam", "eps"]).to_string(index=False))
d = pd.read_csv(R + "mindrop_data.csv"); nd = pd.read_csv(R + "nodrop_data.csv")[["cell", "ds", "g_outcome"]]
d = d.merge(nd, on=["cell", "ds"], how="left")
for (f, Lv), g in d.groupby(["fam", "L"]):
    n = len(g); x = g.D_lo
    print(f"\n== data {f} {Lv}  n={n}  status {g.status.value_counts().to_dict()}  seconds mean {g.seconds.mean():.1f} max {g.seconds.max():.1f}")
    print("D_lo min/q1/med/q3/max/mean:", [round(v, 4) for v in (x.min(), x.quantile(.25), x.median(), x.quantile(.75), x.max(), x.mean())])
    for thr in (0.25, 0.5, 1.0):
        k = int((x > thr).sum()); print(f"  D_lo > {thr}: {k}/{n}  CP95 {tuple(round(v, 3) for v in cp(k, n))}")
    k = int((x >= 0.25).sum()); print("  D_lo >= 0.25 (bracket lower end):", k)
    print("  stored global DROP:", int((g.g_outcome == "DROP").sum()))
    a = int(((g.D_lo >= 0.25) & (g.g_outcome != "DROP")).sum()); b = int(((g.D_hi < 0.245) & (g.g_outcome == "DROP")).sum())
    st = int(((g.D_lo < 0.25) & (g.D_hi >= 0.245) & (g.status == "OK")).sum())
    print(f"  B2: D_lo>=0.25 but stored non-DROP: {a};  D_hi<0.245 but stored DROP: {b};  straddling 0.25: {st}")
    print("  stored DROP with D_lo==0 (status):", g[(g.g_outcome == 'DROP') & (g.D_lo == 0)].status.value_counts().to_dict())
    print("  datasets with an INCONCLUSIVE step:", int((g.n_incon > 0).sum()))
    ed = g.exact_drop.iloc[0]; k = int((g.D_lo <= ed).sum()); print(f"  B4: D_lo <= exact drop {ed:.4f}: {k}/{n} CP95 {tuple(round(v, 3) for v in cp(k, n))}")
    print("  D_lo histogram (0, (0,.25), [.25,.5), [.5,1), >=1):", [int((x == 0).sum()), int(((x > 0) & (x < .25)).sum()), int(((x >= .25) & (x < .5)).sum()), int(((x >= .5) & (x < 1)).sum()), int((x >= 1).sum())])
try:
    c = pd.read_csv(R + "mindrop_controls.csv"); print("\n== controls (B3)")
    for (f, Lv), g in c.groupby(["fam", "L"]):
        ed = g.exact_drop.iloc[0]; print(f, Lv, "n", len(g), "D_lo max", g.D_lo.max(), "exact drop", round(ed, 4), "violations of D_lo<=exact+0.005:", int((g.D_lo > ed + 0.005).sum()), g.status.value_counts().to_dict())
except FileNotFoundError: print("no controls file")
