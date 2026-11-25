"""v8 Fig. 1 (presentation only: every plotted value is read from a stored result file; nothing is simulated).
Run from p22_spectral_dimension_flow:  .venv/bin/python paper/make_figs_v8.py
Writes paper/figs/fig1_budget_v8 (.pdf/.png) and paper/figs/fig1a_v8_numbers.csv. Panel (a) is the bias decomposition of LocQuad on
C_256 x C_16 from results/rt_rbc.csv (red-team round, 100 datasets per budget); panels (b)-(d) are those of make_figs_v7.py, unchanged."""
import re, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = "results"; A = "results_amend10_11"; O = "paper/figs"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 300})
NAME = {"E0": "TwoLine", "E1u": "LocLin-Under", "E1q": "LocQuad", "E2": "TraceSecant", "E4c": "BiasAware-CV", "E4r": "BiasAware-Rule",
        "E6": "StepParam", "Mb": "WML-Boot", "CL": "Collapse-Opt", "E7": "MonoLP"}
COL = {"E1q": "#0072B2", "E4c": "#D55E00", "E6": "#009E73", "E0": "#999999", "E1u": "#56B4E9", "E2": "#CC79A7", "E4r": "#E69F00",
       "Mb": "#7F7F7F", "CL": "#000000", "E7": "#6A3D9A"}
rows = []
def save(fig, n):
    fig.savefig(f"{O}/{n}.pdf", bbox_inches="tight"); fig.savefig(f"{O}/{n}.png", bbox_inches="tight", dpi=300); plt.close(fig)
def num(s): return float(str(s).split(" ")[0])

# ---- collect coverage against walker budget in the three confirmatory runs
W = [8, 32, 128]
from scipy.stats import beta
def cpint(k, n):
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1); hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return 100 * lo, 100 * hi
runs = {"replication": {}, "fixed": {}, "holdout": {}}   # run -> arm -> list over budgets of (covered, answered, generated)
T = pd.read_csv(f"{R}/amend9_b_table.csv")
for e in ("E1q", "E4c", "E7"):
    for w in W:
        m = re.search(r"(\d+)/(\d+) = ", T[(T.wpn == w) & (T.arm == e)].X_cov.iloc[0])
        runs["replication"].setdefault(e, []).append((int(m[1]), int(m[2]), 3600, f"{R}/amend9_b_table.csv"))
for line in open("paper/tab_r3_10c_pooled.tex"):
    c = [x.strip() for x in line.replace("\\\\", "").split("&")]
    if len(c) > 6 and c[1] == "fixed" and c[0] in ("E1q", "E4c"):
        runs["fixed"].setdefault(c[0], []).append((int(c[5]), int(c[4]), 3600, "results_round3/r3_10c_pooled.csv (paper/tab_r3_10c_pooled.tex)"))
a = pd.read_csv(f"{A}/agg_amend11_cells.csv").merge(pd.read_csv(f"{A}/cells_meta_amend11.csv")[["id", "kind"]], left_on="cell", right_on="id")
a = a[a.kind == "crossover"]
for e in ("E1q", "E4c", "E7"):
    for w in W:
        d = a[(a.est == e) & (a.wpn == w)]
        runs["holdout"].setdefault(e, []).append((int(d.k_cov.sum()), int(d.n_valid.sum()), int(d.n.sum()), f"{A}/agg_amend11_cells.csv"))
RUN = {"replication": "Frozen replication (own window)", "fixed": "Fixed windows", "holdout": "Holdout (own window)"}
tex = []
for r in ("replication", "fixed", "holdout"):
    for e in ("E1q", "E4c", "E7"):
        if e not in runs[r]: continue
        cells = []
        for w, (k, n, g, src) in zip(W, runs[r][e]):
            rows.append(("fig1b-d/tableIV", r, NAME[e], e, w, k, n, round(100 * k / n, 2), g, round(100 * k / g, 2), src))
            cells += [f"{n}", f"{100 * k / n:.1f}", f"{100 * k / g:.1f}"]
        tex.append(f"{RUN[r] if e == 'E1q' else ''} & {NAME[e]} & " + " & ".join(cells) + " \\\\")
    if r != "holdout": tex.append("\\hline")


# ---- Fig 1(a): bias against band half-width, LocQuad on C_256 x C_16, fixed window, stored per-dataset vectors
rt = pd.read_csv(f"{R}/rt_rbc.csv"); rt = rt[(rt.est == "LQ") & (rt.fam == "aniso2_w16") & (rt.L == 256) & (rt.abstain == 0)]
TG = []; t_ = [24, 34]
while len(TG) < 40: TG += t_; t_ = [2 * x for x in t_]
vec = lambda s: np.array([float(v) for v in s.split(";")])
fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.35), gridspec_kw={"width_ratios": [1.75, 1, 1, 1]})
ax = axs[0]; arows = []
for w, al in ((8, 0.13), (32, 0.26), (128, 0.5)):
    d = rt[rt.wpn == w]; J = int(d.J.iloc[0]); assert (d.J == J).all() and len(d) == 100
    B = np.array([vec(s) for s in d.b_vec]); HW = np.array([c * vec(s) for c, s in zip(d.crit, d.se_vec)]); t = np.array(TG[:J], float)
    b = B.mean(0); hw = np.median(HW, 0)
    ax.fill_between(t, -hw, hw, color=COL["E1q"], alpha=al, lw=0, label=f"half-width, {w}")
    ax.plot(t, b, color="k", lw=1.1, ls={8: ":", 32: "--", 128: "-"}[w], label=("smoothing bias (8, 32, 128)" if w == 128 else None))
    for tt, bb, hh in zip(t, b, hw): arows.append((w, int(tt), round(bb, 4), round(hh, 4), len(d), f"{R}/rt_rbc.csv"))
ax.axhline(0, color="0.4", lw=0.5); ax.set_xscale("log"); ax.set_xlabel("walk time $t$"); ax.set_ylabel("error in $d_s(t)$")
ax.set_ylim(-1.25, 0.72); ax.legend(frameon=False, loc="lower left", fontsize=6, handlelength=1.6, ncol=2, columnspacing=0.9)
ax.set_title("(a) $C_{256}\\times C_{16}$: bias stays, noise shrinks", fontsize=7.5)
pd.DataFrame(arows, columns=["wpn", "t", "mean_bias_dpop_minus_d", "median_band_halfwidth", "n_datasets", "source"]).to_csv(f"{O}/fig1a_v8_numbers.csv", index=False)

ttl = {"replication": "(b) frozen replication", "fixed": "(c) fixed windows", "holdout": "(d) holdout"}
MK = {"E1q": "o", "E4c": "s", "E7": "D"}
for ax, r in zip(axs[1:], ("replication", "fixed", "holdout")):
    for e in ("E7", "E4c", "E1q"):
        if e not in runs[r]: continue
        k = np.array([v[0] for v in runs[r][e]]); n = np.array([v[1] for v in runs[r][e]]); g = np.array([v[2] for v in runs[r][e]])
        p = 100 * k / n; ci = np.array([cpint(a_, b_) for a_, b_ in zip(k, n)])
        ax.errorbar(range(3), p, yerr=[p - ci[:, 0], ci[:, 1] - p], color=COL[e], marker=MK[e], ms=3.5, lw=1.2, capsize=2, label=NAME[e])
        if e != "E7": ax.plot(range(3), 100 * k / g, color=COL[e], marker=MK[e], ms=3.5, lw=0.9, ls="--", mfc="white")
        if e == "E1q":
            for i, q in enumerate(p): ax.annotate(f"{q:.1f}", (i, q), xytext=(4, 4), textcoords="offset points", fontsize=6, color=COL[e])
    ax.axhline(93, color="0.5", ls=":", lw=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(["8", "32", "128"]); ax.set_xlim(-0.3, 2.6); ax.set_title(ttl[r], fontsize=7.5); ax.set_ylim(0, 104)
    ax.set_xlabel("walkers per node")
    if r != "replication": ax.set_yticklabels([])
    if r == "fixed": ax.text(0.08, 0.60, "MonoLP not run\nin this comparison", transform=ax.transAxes, fontsize=6, color=COL["E7"])
axs[1].set_ylabel("crossover coverage (%)", labelpad=1)
from matplotlib.lines import Line2D
h = [Line2D([], [], color=COL[e], marker=MK[e], ms=3.5, lw=1.2, label=NAME[e]) for e in ("E1q", "E4c", "E7")]
h += [Line2D([], [], color="0.3", lw=1.2, label="% of answered datasets"), Line2D([], [], color="0.3", lw=0.9, ls="--", marker="o", mfc="white", ms=3.5, label="% of all generated datasets")]
fig.legend(handles=h, frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.62, -0.13), fontsize=6.5)
fig.subplots_adjust(wspace=0.12); p0 = axs[0].get_position(); axs[0].set_position([p0.x0 - 0.045, p0.y0, p0.width, p0.height]); save(fig, "fig1_budget_v8")
