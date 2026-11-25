"""Audit figure for the paper (HadCET local DFA-2 slope), from audit/results. Run: .venv/bin/python paper/make_fig_audit.py"""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 300})
A = "audit/results"
real = pd.read_csv(f"{A}/cet_real.csv"); cur = json.load(open(f"{A}/cet_curves.json")); pop = json.load(open(f"{A}/cet_pop.json"))
arr = lambda s: np.array(str(s).split(";"), float)
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.axvspan(5, 20, color="0.9", zorder=0, lw=0)
for p, c, lab in (("P1", "#0072B2", "1851--1950"), ("P2", "#D55E00", "1951--2025")):
    s = np.array(cur[p]["scales"], float)
    r = real[(real.period == p) & (real.est == "E4c") & (real.delta == 0.25)].iloc[0]
    ax.fill_between(s, arr(r.lo), arr(r.hi), color=c, alpha=0.2, lw=0); ax.plot(s, arr(r.d), color=c, lw=1.3, label=f"{lab}, BiasAware-CV")
    r = real[(real.period == p) & (real.est == "E6") & (real.delta == 0.25)].iloc[0]
    ax.plot(s, arr(r.lo), color=c, ls=":", lw=0.9); ax.plot(s, arr(r.hi), color=c, ls=":", lw=0.9, label=f"{lab}, StepParam")
s = np.array(cur["P1"]["scales"], float)
ax.plot(s, pop["N0"]["d_coarse"], "k--", lw=1.0, label="fGn $H=0.65$ null (no crossover)")
ax.set_xscale("log"); ax.set_xlabel("DFA scale $s$ (days)"); ax.set_ylabel("local DFA-2 slope $H(s)$"); ax.set_ylim(-0.3, 1.8); ax.set_xlim(5, 1000)
ax.legend(frameon=False, loc="upper right")
fig.savefig("paper/figs/fig4_cet_slope_v7.pdf", bbox_inches="tight"); fig.savefig("paper/figs/fig4_cet_slope_v7.png", bbox_inches="tight", dpi=300)
