"""Figures for the paper, from existing result CSVs (+ figs/fig1_data.npz from make_fig1_data.py). Run: .venv/bin/python paper/make_figs.py"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = "results"; O = "paper/figs"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 300})
COL = {"E1q": "#0072B2", "E4c": "#D55E00", "E6": "#009E73", "E0": "#999999", "E1u": "#56B4E9", "E2": "#CC79A7", "E4r": "#E69F00",
       "Mb": "#7F7F7F", "CL": "#000000", "E3": "#F0E442"}
def save(fig, n):
    fig.savefig(f"{O}/{n}.pdf", bbox_inches="tight"); fig.savefig(f"{O}/{n}.png", bbox_inches="tight", dpi=300); plt.close(fig)

# ---- Fig 1
z = np.load(f"{O}/fig1_data.npz")
fig, ax = plt.subplots(figsize=(3.4, 2.5))
m = z["tgrid"] < 2.5e3
ax.plot(z["tgrid"][m], z["dgrid"][m], "k-", lw=1.2, label="exact $d(t)$", zorder=5)
off = {"E1q": 0.985, "E4c": 1.0, "E6": 1.015}
for nm, lo, hi in (("E1q", z["e1q_lo"], z["e1q_hi"]), ("E4c", z["e4_lo"], z["e4_hi"]), ("E6", z["e6_lo"], z["e6_hi"])):
    t = z["t"] * off[nm]
    ax.fill_between(t, lo, hi, step=None, color=COL[nm], alpha=0.22, lw=0)
    ax.plot(t, lo, color=COL[nm], lw=0.8); ax.plot(t, hi, color=COL[nm], lw=0.8, label=nm.replace("E1q", "E1q (local-quadratic)").replace("E4c", "E4c (bias-aware)").replace("E6", "E6 (logistic step)"))
ax.set_xscale("log"); ax.set_xlabel("walk time $t$"); ax.set_ylabel("$d_s(t)$"); ax.set_ylim(-0.5, 3.4); ax.set_xlim(20, 2.5e3); ax.axhline(0, color="0.85", lw=0.5); ax.legend(frameon=False, loc="upper right")
ax.text(0.03, 0.05, r"$C_{256}\times C_{16}$ torus, one dataset", transform=ax.transAxes, fontsize=7)
save(fig, "fig1_example_bands")

# ---- Fig 2 frontier
P = pd.read_csv(f"{R}/amend4_final_table.csv").set_index("est")
Q = pd.read_csv(f"{R}/amend5_passbar.csv").set_index("est")
def num(s): return float(str(s).split(" ")[0])
def minl(s): return None if str(s).startswith("not") else float(s)
pts = {}
for e in ("E0", "E1u", "E1q", "E2", "E4c", "E4r", "E6"):
    pts[e] = (P.loc[e, "minL_aniso2_w16"], P.loc[e, "minL_aniso3_w8"], num(P.loc[e, "cov_cross_real"]))
for e in ("Mb", "CL"):
    pts[e] = (Q.loc[e, "minL_aniso2_w16"], Q.loc[e, "minL_aniso3_w8"], num(Q.loc[e, "cov_xover_real"]))
pts["E1u"] = (">256", pts["E1u"][1], pts["E1u"][2])
pts["E4c"] = (">768", pts["E4c"][1], pts["E4c"][2]); pts["E4r"] = ("768", pts["E4r"][1], pts["E4r"][2])  # 2D min L from Amendment 6 (extended sweep)
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
for k, (ax, ttl) in enumerate(zip(axs, ("2D$\\to$1D ($C_L\\times C_{16}$)", "3D$\\to$1D ($C_L\\times C_8^2$)"))):
    for e, v in pts.items():
        x = v[k]
        if x is None or str(x).startswith("not"): continue
        arrow = str(x).startswith(">"); xv = float(str(x).lstrip(">")); y = v[2]
        if k == 0 and e == "Mb": xv *= 1.07   # E6 and Mb coincide at (256, 99); Mb is displaced right for visibility (stated in the caption)
        ax.scatter(xv, y, s=28, color=COL[e], edgecolor="k", lw=0.4, zorder=3, marker=">" if arrow else "o")
        off = {"E6": (-14, 5), "Mb": (4, -9), "E4c": (-14, 5), "E4r": (4, -9), "E0": (4, -9), "E1u": (-16, 4)}.get(e, (4, 3))
        if k == 0 and e == "E1u": off = (6, -10)
        if k == 0 and e == "E2": off = (-14, 5)
        ax.annotate(e, (xv, y), xytext=off, textcoords="offset points", fontsize=7)
    ax.axhline(93, color="0.6", ls=":", lw=0.8); ax.set_xscale("log"); ax.set_title(ttl, fontsize=8)
    ax.set_xlabel("minimum $L$ with power $\\geq 0.8$")
    ax.set_xticks([48, 64, 96, 128, 192, 256, 512, 768]); ax.set_xticklabels(["48", "64", "96", "128", "192", "256", "512", "768"])
    ax.minorticks_off()
axs[0].set_ylabel("simultaneous coverage of exact $d(t)$\nthrough realised crossovers (%)")
axs[0].set_ylim(0, 104)
axs[0].text(0.02, 0.06, "triangle: not reached at the\nlargest $L$ tested for that estimator", transform=axs[0].transAxes, fontsize=6.5)
save(fig, "fig2_frontier")

# ---- Fig 3 separation ratio
S = pd.read_csv(f"{R}/amend6_separation_scan.csv")
fig, ax = plt.subplots(figsize=(3.4, 2.5))
for e, c in (("E6", COL["E6"]), ("E4c", COL["E4c"])):
    d = S[S.est == e]
    cv = d["cov"].map(num).values
    ax.errorbar(d.ratio, cv, yerr=[cv - 100 * d.cov_lo.values, 100 * d.cov_hi.values - cv], color=c, marker="o", ms=4, lw=1, capsize=2, label=e)
ax.axhline(93, color="0.5", ls=":", lw=0.8); ax.axhline(95, color="0.8", lw=0.6)
ax.set_xlabel("step-separation ratio $b/8$"); ax.set_ylabel("simult. coverage of exact $d(t)$ (%)"); ax.set_ylim(0, 102); ax.legend(frameon=False, loc="lower left")
save(fig, "fig3_e6_validity")

# ---- Fig 4 causal sets
D = pd.read_csv(f"{R}/causet_reference_peak.csv")
fig, axs = plt.subplots(1, 2, figsize=(3.4, 2.3), sharey=False)
for ax, dim, ttl in zip(axs, (1, 3), ("1+1", "3+1")):
    for con, c, lab in (("k0", "#0072B2", "causet (Hasse)"), ("rgg0", "#999999", "matched RGG")):
        d = D[(D.dim == dim) & (D.con == con)].sort_values("N")
        ax.plot(d.N, d.d_peak, "o-", ms=3, lw=1, color=c, label=lab)
    ax.set_xscale("log"); ax.set_title(ttl, fontsize=8); ax.set_xlabel("$N$"); ax.set_xticks([500, 2000, 8000]); ax.set_xticklabels(["500", "2000", "8000"]); ax.minorticks_off()
axs[0].set_ylabel("peak of $\\bar d_s(t)$ (ensemble mean)"); axs[0].legend(frameon=False, fontsize=6, loc="upper left")
fig.tight_layout(); save(fig, "fig4_causet_peak")
print("ok", pts)
