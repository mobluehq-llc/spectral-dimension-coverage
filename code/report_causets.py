"""Tables + figures for the causal-set section. Reads results/causet_estimates.csv, causet_curves.json, causet_graph_meta.csv, causet_p22_calibration.csv."""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
df = pd.read_csv("results/causet_estimates.csv"); cur = json.load(open("results/causet_curves.json")); meta = pd.read_csv("results/causet_graph_meta.csv")
# ---- graph diagnostics
g = meta.groupby(["dim", "N"]).agg(R=("rep", "size"), ordering_fraction=("rel_frac", "mean"), full_deg=("full_deg", "mean"), deg_k0=("deg_k0", "mean"), deg_k3=("deg_k3", "mean"),
                                    deg_k15=("deg_k15", "mean"), Neff_k0=("Neff_k0", "mean")).round(3)
g.to_csv("results/causet_graph_diag.csv"); print(g)
# ---- abstention and Delta table (primary: G4, thr=1.3 and G2/thr=6 side by side)
d = df.copy(); d["delta"] = np.where(d.est == "E3", d.d_late - d.d_early, np.where(d.decl_neg == 1, d.dmin, d.dmax))
d["flow"] = np.select([d.decl_neg == 1, d.decl_pos == 1], ["drop", "rise"], "none"); d.loc[d.abstain == 1, "flow"] = "abstain"
d.to_csv("results/causet_estimates_annotated.csv", index=False)
def cell(x):
    ok = x[x.abstain == 0]
    if not len(ok): return "abst"
    return f"{ok.delta.iloc[0]:+.2f} ({ok.flow.iloc[0][0]})"
rows = []
for (dim, con, est, grid, thr), x in d[d.est != "E3"].groupby(["dim", "con", "est", "grid", "thr"]):
    rows.append(dict(dim=dim, con=con, est=est, grid=grid, thr=thr, **{str(N): cell(x[x.N == N]) for N in sorted(d.N.unique())}))
for (dim, con, grid), x in d[d.est == "E3"].groupby(["dim", "con", "grid"]):
    rows.append(dict(dim=dim, con=con, est="E3", grid=grid, thr=np.nan, **{str(N): cell(x[x.N == N]) for N in sorted(d.N.unique())}))
tab = pd.DataFrame(rows); tab.to_csv("results/causet_delta_table.csv", index=False)
# abstention rate pooled per estimator x kind
d["kind"] = np.where(d.con.str.startswith("rgg"), "RGG control", "causet")
ab = d.groupby(["est", "kind", "thr"]).abstain.mean().unstack([1, 2]).round(2); ab.to_csv("results/causet_abstention.csv"); print(ab)
# ---- figures
cs = pd.DataFrame(cur)
if len(cs):
    for fig_i, (grid, thr) in enumerate([("G4", 1.3)]):
        fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharey="row")
        for r, dim in enumerate((1, 3)):
            for c, k in enumerate(("0", "3", "15")):
                ax = axs[r, c]
                for con, ls in ((f"k{k}", "-"), (f"rgg{k}", "--")):
                    sub = cs[(cs.dim == dim) & (cs.con == con) & (cs.grid == grid) & (cs.thr == thr) & (cs.est == "E1q")]
                    for _, row in sub.sort_values("N").iterrows():
                        t = np.array(row.t); dd = np.array(row.d)
                        ax.plot(t, dd, ls, marker="o", ms=2.5, label=f"{con} N={row.N}", color=plt.cm.viridis(np.log2(row.N / 500) / 4))
                ax.set_xscale("log"); ax.set_title(f"{'1+1' if dim == 1 else '3+1'}, k={k} (solid causet, dashed RGG)", fontsize=9)
                ax.set_xlabel("lazy-walk steps t"); ax.grid(alpha=.3)
            axs[r, 0].set_ylabel("E1q  d_s(t)")
        h = [plt.Line2D([], [], color=plt.cm.viridis(i / 4), label=f"N={500 * 2 ** i}") for i in range(5)]
        axs[0, 0].legend(handles=h, fontsize=7)
        fig.suptitle(f"E1q d_s(t) on causet graphs vs matched-degree RGG controls ({grid}, thr={thr}); missing curve = abstention")
        fig.tight_layout(); fig.savefig("results/fig5_causet_dcurves.png", dpi=130); plt.close(fig)
# Delta vs N
fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharey=True)
for r, dim in enumerate((1, 3)):
    for c, est in enumerate(("E0", "E1q", "E2")):
        ax = axs[r, c]
        for con in ("k0", "k3", "k15", "rgg0", "rgg3", "rgg15"):
            x = d[(d.dim == dim) & (d.con == con) & (d.est == est) & (d.grid == "G4") & (d.thr == 1.3) & (d.abstain == 0)].sort_values("N")
            if len(x):
                ax.plot(x.N, x.delta, "-" if con.startswith("k") else "--", marker="o", label=con)
        ax.axhspan(-0.25, 0.25, color="gray", alpha=.2); ax.set_xscale("log"); ax.set_title(f"{'1+1' if dim == 1 else '3+1'} {est}", fontsize=9); ax.set_xlabel("N"); ax.grid(alpha=.3)
    axs[r, 0].set_ylabel("Delta = d(late) - d(early)")
axs[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("Change statistic vs N (G4, thr=1.3; shaded = equivalence margin; gaps = abstention)"); fig.tight_layout(); fig.savefig("results/fig6_causet_delta_vs_N.png", dpi=130)
print(tab.to_string())

# ---- stability: d_hat(t=8) and its N-exponent (E1q, G4, thr 1.3)
rows = []
for (dim, con), x in cs[(cs.est == "E1q") & (cs.grid == "G4") & (cs.thr == 1.3)].groupby(["dim", "con"]):
    Ns_, ds_ = [], []
    for _, r in x.sort_values("N").iterrows():
        t = np.array(r.t); dd = np.array(r.d)
        if 8 in t: Ns_.append(r.N); ds_.append(dd[list(t).index(8)])
    if len(Ns_) >= 3:
        a = np.polyfit(np.log(Ns_), np.log(ds_), 1)[0]
        rows.append(dict(dim=dim, con=con, **{f"d(t=8),N={n}": round(v, 2) for n, v in zip(Ns_, ds_)}, exponent_dlogd_dlogN=round(a, 3)))
st = pd.DataFrame(rows); st.to_csv("results/causet_stability_t8.csv", index=False); print(st.to_string())
# ---- (c) window sufficiency: span log10(tlast/t_first) vs N, E1q rule, linear in log2 N, extrapolated to p22's required span
cal = pd.read_csv("results/causet_p22_calibration.csv")
need = {"E0 (E1x window)": ("E1x", ), "E1q": ("E1q",), "E2": ("E2",)}
rows = []
for dim in (1, 3):
    for con in ("k0", "k3", "k15", "rgg0", "rgg3", "rgg15"):
        for grid, t0 in (("G2", 2), ("G4", 4)):
            x = d[(d.est == "E1q") & (d.dim == dim) & (d.con == con) & (d.grid == grid) & (d.thr == 1.3) & (d.abstain == 0)].sort_values("N")
            if len(x) < 3: rows.append(dict(dim=dim, con=con, grid=grid, note="<3 non-abstaining N")); continue
            span = np.log10(x.tlast.values / t0); lg = np.log2(x.N.values)
            b, a0 = np.polyfit(lg, span, 1)
            rows.append(dict(dim=dim, con=con, grid=grid, span_N500=round(span[0], 2), span_N8000=round(span[-1], 2), slope_per_doubling=round(b, 3),
                             N_for_span_0p9=(float(2 ** ((0.903 - a0) / b)) if b > 0.003 else np.inf), N_for_span_1p66=(float(2 ** ((1.656 - a0) / b)) if b > 0.003 else np.inf)))
cc = pd.DataFrame(rows); cc.to_csv("results/causet_window_span_extrapolation.csv", index=False); print(cc.to_string())
