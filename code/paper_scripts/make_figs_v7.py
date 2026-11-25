"""v7 figures (presentation only: every plotted value is read from a stored result file; nothing is simulated).
Run from p22_spectral_dimension_flow:  .venv/bin/python paper/make_figs_v7.py
Writes paper/figs/fig1_example_bands_v7, fig2_budget_v7, fig2_frontier_v7, fig3_step_validity_v7 (.pdf/.png)
and paper/figs/fig_v7_numbers.csv (every plotted number of Figs. 2 and 3 with its source file)."""
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
open("paper/tab_v7_joint.tex", "w").write("\n".join(tex) + "\n")

# ---- Fig 1: (a) one illustrative dataset per budget; (b-d) coverage against budget
z = np.load(f"{O}/fig1_budget_data.npz")
fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.35), gridspec_kw={"width_ratios": [1.75, 1, 1, 1]})
ax = axs[0]; m = z["tgrid"] < 2.5e3
ax.plot(z["tgrid"][m], z["dgrid"][m], "k-", lw=1.2, label="exact $d(t)$", zorder=5)
for w, al, ls in ((8, 0.16, "--"), (128, 0.45, "-")):
    ax.fill_between(z[f"t{w}"], z[f"q_lo{w}"], z[f"q_hi{w}"], color=COL["E1q"], alpha=al, lw=0)
    ax.plot(z[f"t{w}"], z[f"q_lo{w}"], color=COL["E1q"], lw=0.8, ls=ls); ax.plot(z[f"t{w}"], z[f"q_hi{w}"], color=COL["E1q"], lw=0.8, ls=ls, label=f"LocQuad, {w} walkers/node")
ax.plot(z["t128"], z["b_lo128"], color=COL["E4c"], lw=0.8, ls=":"); ax.plot(z["t128"], z["b_hi128"], color=COL["E4c"], lw=0.8, ls=":", label="BiasAware-CV, 128")
ax.set_xscale("log"); ax.set_xlabel("walk time $t$"); ax.set_ylabel("$d_s(t)$"); ax.set_ylim(-1.9, 3.6); ax.set_xlim(20, 2.5e3)
ax.legend(frameon=False, loc="lower left", fontsize=6, handlelength=1.6); ax.set_title("(a) $C_{256}\\times C_{16}$: one dataset per budget", fontsize=7.5)
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
fig.subplots_adjust(wspace=0.12); p0 = axs[0].get_position(); axs[0].set_position([p0.x0 - 0.045, p0.y0, p0.width, p0.height]); save(fig, "fig1_budget_v7")

# ---- Fig 3: coverage against detection size, per family
P = pd.read_csv(f"{R}/amend4_cell_summary.csv"); P = P[P.real == True]
P = P.assign(nv=(200 * (1 - P.abstain)).round(), k=(P["cov"] * 200 * (1 - P.abstain)).round())
fam = {"aniso2_w16": 0, "aniso3_w8": 1}
cov = {}
for e in ("E0", "E1u", "E1q", "E2", "E4c", "E4r", "E6"):
    for f, i in fam.items():
        d = P[(P.est == e) & (P.fam == f)]; cov[(e, i)] = (int(d.k.sum()), int(d.nv.sum()), f"{R}/amend4_cell_summary.csv")
r5 = pd.read_csv(f"{R}/raw_amend5.csv").merge(pd.read_csv(f"{R}/cells_meta_amend5.csv"), left_on="cell", right_on="id")
r5 = r5[(r5.kind == "crossover") & (r5.max_drop >= 0.6) & (r5.abstain == 0)]
for e in ("Mb", "CL"):
    for f, i in fam.items():
        d = r5[(r5.est == e) & (r5.fam == f)]; cov[(e, i)] = (int(d.cov_simult.sum()), len(d), f"{R}/raw_amend5.csv")
E = pd.read_csv(f"{R}/amend8_cell_table.csv"); cells = P[["fam", "L"]].drop_duplicates()   # the eight realised cells of the main benchmark
E = E.merge(cells, on=["fam", "L"])
for f, i in fam.items():
    d = E[E.fam == f]; cov[("E7", i)] = (int(round((d.E7_cov * 200).sum())), 200 * len(d), f"{R}/amend8_cell_table.csv")
F4 = pd.read_csv(f"{R}/amend4_final_table.csv").set_index("est"); F5 = pd.read_csv(f"{R}/amend5_passbar.csv").set_index("est")
ml = {e: [str(F4.loc[e, "minL_aniso2_w16"]), str(F4.loc[e, "minL_aniso3_w8"])] for e in ("E0", "E1u", "E1q", "E2", "E4c", "E4r", "E6")}
for e in ("Mb", "CL"): ml[e] = [str(F5.loc[e, "minL_aniso2_w16"]), str(F5.loc[e, "minL_aniso3_w8"])]
ml["E1u"][0] = ">256"; ml["E4c"][0] = ">768"; ml["E4r"][0] = "768"   # 2D thresholds from the extended sweep (Table III of the paper)
ml["E7"] = [">768", "256"]                                             # endpoint test (RESULTS.md, Amendment 8)
fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7), sharey=True)
MKS = {"E0": "o", "E1u": "v", "E1q": "o", "E2": "^", "E4c": "s", "E4r": "s", "E6": "P", "Mb": "X", "CL": "*", "E7": "D"}
JIT = {(0, "E0"): 0.95, (0, "E2"): 1.05, (0, "E6"): 0.95, (0, "Mb"): 1.05, (0, "E4c"): 0.93, (0, "E7"): 1.08, (1, "E0"): 0.94, (1, "E6"): 1.0, (1, "E1q"): 1.06, (1, "E1u"): 1.0,
       (1, "E2"): 0.95, (1, "Mb"): 1.05, (1, "E4c"): 0.95, (1, "E4r"): 1.05}
for k, (ax, tt) in enumerate(zip(axs, ("2D$\\to$1D ($C_L\\times C_{16}$)", "3D$\\to$1D ($C_L\\times C_8^2$)"))):
    for e in ml:
        x = ml[e][k]; arrow = x.startswith(">"); xv = float(x.lstrip(">")) * JIT.get((k, e), 1.0); kk, nn, src = cov[(e, k)]; y = 100 * kk / nn
        rows.append(("fig2", "2D" if k == 0 else "3D", NAME[e], e, 32, kk, nn, round(y, 2), "", "", src))
        ax.scatter(xv, y, s=34, color=COL[e], edgecolor="k", lw=0.4, zorder=3, marker=MKS[e], label=NAME[e] if k == 0 else None)
        if arrow: ax.annotate("", xy=(xv * 1.22, y), xytext=(xv * 1.04, y), arrowprops=dict(arrowstyle="->", lw=0.7, color=COL[e]))
    if k == 1:   # joint test on the same MonoLP constraint set (discretised measure; 3D only)
        kk, nn, _ = cov[("E7", 1)]; ax.scatter(160, 100 * kk / nn, s=34, facecolor="none", edgecolor=COL["E7"], lw=1.0, zorder=3, marker="D", label="MonoLP-Joint")
    ax.axhline(93, color="0.6", ls=":", lw=0.8); ax.set_xscale("log"); ax.set_title(tt, fontsize=8)
    ax.set_xlabel("smallest tested $L$ with power $\\geq 0.8$")
    ax.set_xticks([48, 64, 96, 128, 192, 256, 512, 768]); ax.set_xticklabels(["48", "64", "96", "128", "192", "256", "512", "768"]); ax.minorticks_off()
    ax.set_xlim(42, 1050)
axs[0].set_ylabel("coverage of exact $d(t)$ through realised\ncrossovers of that family (% of answered)"); axs[0].set_ylim(0, 106)
axs[0].text(0.02, 0.05, "arrow: power 0.8 not reached at the\nlargest $L$ tested for that estimator", transform=axs[0].transAxes, fontsize=6.5)
h0, l0 = axs[0].get_legend_handles_labels(); h1, l1 = axs[1].get_legend_handles_labels()
fig.legend(h0 + h1, l0 + l1, frameon=False, loc="center left", bbox_to_anchor=(0.905, 0.5), fontsize=6.5)
fig.subplots_adjust(wspace=0.06); save(fig, "fig2_frontier_v7")

# ---- Fig 4: validity boundary of the step assumption
S = pd.read_csv(f"{R}/amend6_separation_scan.csv")
fig, ax = plt.subplots(figsize=(3.4, 2.5))
for e in ("E6", "E4c"):
    d = S[S.est == e]; cv = d["cov"].map(num).values
    ax.errorbar(d.ratio, cv, yerr=[cv - 100 * d.cov_lo.values, 100 * d.cov_hi.values - cv], color=COL[e], marker="o", ms=4, lw=1, capsize=2, label=NAME[e])
ax.axhline(93, color="0.5", ls=":", lw=0.8); ax.axhline(95, color="0.8", lw=0.6)
ax.set_xlabel("step-separation ratio $b/8$"); ax.set_ylabel("simult. coverage of exact $d(t)$ (%)"); ax.set_ylim(0, 102); ax.legend(frameon=False, loc="lower left")
save(fig, "fig3_step_validity_v7")

pd.DataFrame(rows, columns=["figure", "panel", "name", "archive_code", "walkers_per_node", "covered", "answered", "pct_answered", "generated", "pct_all", "source"]).to_csv(f"{O}/fig_v7_numbers.csv", index=False)
print(pd.DataFrame(rows).to_string())
