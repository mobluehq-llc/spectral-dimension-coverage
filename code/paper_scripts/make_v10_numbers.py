"""v10 evidence lines: minimum compatible drop, continuum certificate, same-box comparison (stored per-dataset rows only; no simulation).
Run from p22_spectral_dimension_flow:  archive/.venv/bin/python paper/make_v10_numbers.py  -> results/v10_mindrop_cert_numbers.txt
Inputs: results/nodrop_controls.csv, results/mindrop_{exact,data,controls}.csv (rt_mindrop.py), results/certify_{2d,sentinels,sentinels_ctl}.csv (rt_certify.py),
results/nodrop_data.csv (rt_nodrop.py), results/rt_monolp_p2.csv (rt_monolp.py). PREREG_REDTEAM_FINAL.md, amendment 11-24 PM (8a5086d)."""
import pandas as pd
from scipy.stats import beta
def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1); hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return 100 * lo, 100 * hi
out = []
ex = pd.read_csv("results/mindrop_exact.csv")
for _, r in ex.sort_values(["fam", "L", "eps"]).iterrows():
    out.append(f"exact {r.fam} L={int(r.L)} eps {r.eps:g}: exact drop {r.exact_drop:.3f} | D_lo {r.D_lo:.3f} | D_hi {r.D_hi:.3f} | status {r.status}")
d = pd.read_csv("results/mindrop_data.csv"); nd = pd.read_csv("results/nodrop_data.csv")
p2 = pd.read_csv("results/rt_monolp_p2.csv"); p2 = p2[p2.variant == "B0"]
for (fam, L), g in d.groupby(["fam", "L"]):
    n = len(g); q = g.D_lo.quantile([0, .25, .5, .75, 1]).values
    k25 = int((g.D_lo > 0.25).sum()); k50 = int((g.D_lo > 0.5).sum()); lo, hi = cp(k25, n); lo5, hi5 = cp(k50, n)
    below = int((g.D_lo <= g.exact_drop).sum())
    out.append(f"data {fam} L={L} wpn32: n {n} | D_lo min {q[0]:.3f} q1 {q[1]:.3f} median {q[2]:.3f} q3 {q[3]:.3f} max {q[4]:.3f} | D_lo > 0.25: {k25} (CP {lo:.1f}-{hi:.1f}) | D_lo > 0.5: {k50} (CP {lo5:.1f}-{hi5:.1f}) | D_lo <= exact drop {g.exact_drop.iloc[0]:.3f}: {below}/{n}")
    m = g.merge(nd[(nd.fam == fam) & (nd.L == L)][["ds", "g_outcome", "j_drop", "e7_neg"]], on="ds", how="inner")
    assert len(m) == n, f"MISMATCH want {n} joined rows got {len(m)}"
    drop = m.g_outcome == "DROP"; ties = m[drop & (m.D_lo <= 0.25)]
    out.append(f"data {fam} L={L} wpn32: stored global DROP {int(drop.sum())} | D_lo >= 0.25 with stored non-DROP {int(((m.D_lo >= 0.25) & ~drop).sum())} | stored DROP with D_lo <= 0.25: {len(ties)}"
               + (f" (D_lo {ties.D_lo.min():.4f}-{ties.D_lo.max():.4f})" if len(ties) else "")
               + f" | datasets with an inconclusive bisection step {int((g.n_incon > 0).sum())} | LPs mean {g.nlp.mean():.0f} | seconds mean {g.seconds.mean():.1f}")
    s = p2[(p2.fam == fam) & (p2.L == L)].merge(m[["ds", "g_outcome", "j_drop", "e7_neg"]], on="ds", how="inner", suffixes=("_p2", ""))
    out.append(f"samebox {fam} L={L} wpn32: n {len(s)} | secant drop > 0.25: {int((s.sec_drop > 0.25).sum())} | endpoint {int((s.e7_neg > 0).sum())} | joint {int((s.j_drop > 0).sum())} | global DROP {int((s.g_outcome == 'DROP').sum())} | median pointwise width {s.width.median():.2f} | median secant drop {s.sec_drop.median():.2f}")
nc = pd.read_csv("results/nodrop_controls.csv"); mc = nc.merge(p2, on=["fam", "L", "ds"], how="inner", suffixes=("", "_p2"))
assert len(mc) == len(nc), f"MISMATCH want {len(nc)} joined control rows got {len(mc)}"
out.append(f"samebox controls wpn32: n {len(mc)} | box critical value, window and CP count identical in {int(((mc.c_box == mc.c_box_p2) & (mc.Jw == mc.Jw_p2) & (mc.n_cp == mc.n_cp_p2)).sum())} | secant drop > 0.25: {int((mc.sec_drop > 0.25).sum())} | global DROP {int((mc.g_outcome == 'DROP').sum())}")
tot = len(d); ok = int((d.D_lo <= d.exact_drop).sum())
out.append(f"data pooled: D_lo <= exact drop in {ok}/{tot}")
c = pd.read_csv("results/mindrop_controls.csv")
out.append(f"controls pooled: n {len(c)} | D_lo == 0 in {int((c.D_lo == 0).sum())}/{len(c)} | cells {c.groupby(['fam','L']).size().to_dict()}")
out.append(f"exact controls: D_lo == 0 in {int((ex[ex.fam.isin(['ring','torus2','torus3'])].D_lo == 0).sum())}/{int(ex.fam.isin(['ring','torus2','torus3']).sum())}")
z = pd.read_csv("results/certify_2d.csv")
out.append(f"certify 2D: declarations {len(z)} | CERTIFIED {int((z.verdict == 'CERTIFIED').sum())} | not certified {int((z.verdict != 'CERTIFIED').sum())} | solver failures {int(z.n_solver.sum())} | cells {int(z.ncells.sum())} certified {int(z.n_cert.sum())} (per dataset {int(z.ncells.min())}-{int(z.ncells.max())}, mean {z.ncells.mean():.1f}) | datasets with added atoms {int((z.added > 0).sum())} | atoms added {int(z.added.sum())} max {int(z.added.max())} | ball-verified {int((z.mp_verified.astype(str) == 'True').sum())} | float only {int((z.mp_verified.astype(str) != 'True').sum())} | seconds mean {z.seconds.mean():.1f}")
out.append("certify 2D ball-verified ds: " + " ".join(str(int(x)) for x in sorted(z[z.mp_verified.astype(str) == 'True'].ds)))
sn = pd.read_csv("results/certify_sentinels.csv")
e = sn[sn.kind == "exact"]; w = sn[sn.kind != "exact"]
out.append(f"sentinel exact: {len(e)} run | CERTIFIED {int((e.verdict == 'CERTIFIED').sum())} | cells {int(e.ncells.sum())} certified {int(e.n_cert.sum())}")
out.append(f"sentinel witness cells (2D WITNESS datasets): {len(w)} | witness {int((w.wit_status == 'witness').sum())} | certified {int((w.wit_status == 'certified').sum()) + int((w.verdict == 'CERTIFIED').sum())} | kinds {w.kind.value_counts().to_dict()}")
sc = pd.read_csv("results/certify_sentinels_ctl.csv")
out.append(f"sentinel control cells (ring 1024): {len(sc)} | witness {int((sc.wit_status == 'witness').sum())} | certified {int((sc.wit_status != 'witness').sum())}")
open("results/v10_mindrop_cert_numbers.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
