#!/usr/bin/env python3
"""P22 analysis: tables + plots from results/raw.csv and cells_meta.csv (pre-registered pass bar in PLAN.md)."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from p22lib import cp

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
raw = pd.read_csv(os.path.join(RES, "raw.csv"))
meta = pd.read_csv(os.path.join(RES, "cells_meta.csv"), keep_default_na=False)
meta["note"] = meta["note"].astype(str)
d = raw.merge(meta[["id", "fam", "L", "N", "kind", "note", "real"]], left_on="cell", right_on="id")
d["decl"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
d["decl_has"] = d.decl_neg.notna()
EST = ["E0", "E1", "E1x", "E2", "E2h", "E3"]
NAMES = {"E0": "E0 two-segment line", "E0g": "E0g global line", "E1": "E1 local-linear (split)", "E1x": "E1x no split", "E2": "E2 Hutchinson (bootstrap)", "E2h": "E2h Hutchinson (Hoeffding)", "E3": "E3 collapse (R&B-style)"}
COL = {"E0": "#7f7f7f", "E0g": "#bbbbbb", "E1": "#1f77b4", "E1x": "#9ecae1", "E2": "#d62728", "E2h": "#ff9896", "E3": "#2ca02c"}

# ---------------- per cell x estimator table
rows = []
for (cid, est), g in d.groupby(["cell", "est"]):
    m = meta[meta.id == cid].iloc[0]
    n = len(g); na = int(g.abstain.sum()); ok = g[g.abstain == 0]; k = len(ok)
    r = dict(cell=cid, fam=m.fam, L=m.L, N=m.N, kind=m.kind, note=m.note, real=m.real, est=est, n=n, abstain=na / n)
    if k:
        if "cov_simult" in ok and ok.cov_simult.notna().any():
            kk = int(ok.cov_simult.sum()); nn = int(ok.cov_simult.notna().sum()); r.update(cov_simult=kk / nn, cov_lo=cp(kk, nn)[0], cov_hi=cp(kk, nn)[1], n_cov=nn)
            r["cov_point"] = ok.cov_point.mean()
        if "cov_pop_simult" in ok and ok.cov_pop_simult.notna().any():
            r["cov_pop_simult"] = ok.cov_pop_simult.mean(); r["cov_pop_point"] = ok.cov_pop_point.mean()
        if ok.decl_neg.notna().any():
            dd = ok[ok.decl_neg.notna()]
            kf = int(dd.decl.sum()); r.update(decl_rate=kf / len(dd), decl_lo=cp(kf, len(dd))[0], decl_hi=cp(kf, len(dd))[1])
            r["decl_zero_rate"] = dd.decl_zero.mean()
            kp = int(dd.decl_neg.sum()); r["power_signed_nonabs"] = kp / len(dd)
            # power counts abstention as a miss
            kpa = int(dd.decl_neg.sum()); r["power"] = kpa / n; r["power_lo"], r["power_hi"] = cp(kpa, n)
        for c in ("bias", "mae", "hw", "se_mean", "h", "bias_pop", "J"):
            if c in ok and ok[c].notna().any(): r[c] = ok[c].mean()
    rows.append(r)
cells = pd.DataFrame(rows).sort_values(["kind", "fam", "L", "est"])
cells.to_csv(os.path.join(RES, "cell_summary.csv"), index=False)

# ---------------- estimator verdicts
V = []
ctrl = d[(d.kind == "control") & (d.note == "")]
cs = cells[(cells.kind == "control") & (cells.note == "")]
for est in EST:
    g = ctrl[(ctrl.est == est) & (ctrl.abstain == 0)]
    v = dict(est=est)
    gc = g[g.cov_simult.notna()]
    v["ctrl_n"] = len(gc); v["ctrl_abstain"] = float(ctrl[ctrl.est == est].abstain.mean())
    if len(gc):
        k = int(gc.cov_simult.sum()); v["cov_pooled"] = k / len(gc); v["cov_pooled_ci"] = "%.3f-%.3f" % cp(k, len(gc))
        v["cov_pointwise_pooled"] = gc.cov_point.mean()
        c_ = cs[(cs.est == est) & (cs.n_cov >= 20)]
        bad = c_[c_.cov_hi < 0.93]
        v["cov_bad_cells"] = ";".join(f"{a}{b}" for a, b in zip(bad.fam, bad.L)) or "none"
        v["C_pass"] = bool(v["cov_pooled"] >= 0.93 and len(bad) == 0)
    gf = g[g.decl_has]
    if len(gf):
        k = int(gf.decl.sum()); v["fpr_pooled"] = k / len(gf); v["fpr_ci"] = "%.3f-%.3f" % cp(k, len(gf)); v["fpr_zero_rule"] = gf.decl_zero.mean()
        f_ = cs[(cs.est == est) & cs.decl_rate.notna()]
        f_ = f_[(f_.n - (f_.abstain * f_.n)) >= 20]
        badf = f_[f_.decl_lo > 0.05]
        v["fpr_bad_cells"] = ";".join(f"{a}{b}" for a, b in zip(badf.fam, badf.L)) or "none"
        v["F_pass"] = bool(v["fpr_pooled"] <= 0.05 and len(badf) == 0)
    V.append(v)
V = pd.DataFrame(V)

# minimum detectable L
MIN = []
cr = cells[(cells.kind == "crossover") & (cells.real == True) & (cells.note == "")]
for est in EST:
    for fam in ("aniso2_w16", "aniso3_w8"):
        c = cr[(cr.est == est) & (cr.fam == fam)].sort_values("L").reset_index(drop=True)
        rec = dict(est=est, fam=fam, tested_real_L=",".join(map(str, c.L)))
        if c.power.isna().all() or len(c) == 0:
            rec.update(minL_point="n/a", minL_strict="n/a"); MIN.append(rec); continue
        pw = c.power.fillna(0).values; lo = c.power_lo.fillna(0).values
        def first(cond):
            for i in range(len(c)):
                if cond[i:].all(): return i
            return None
        ip = first(pw >= 0.8); isx = first(lo >= 0.8)
        rec["minL_point"] = int(c.L[ip]) if ip is not None else "not reached (>%d)" % c.L.max()
        rec["power_at_minL_point"] = float(pw[ip]) if ip is not None else np.nan
        rec["power_ci_at_minL_point"] = "%.2f-%.2f" % (c.power_lo[ip], c.power_hi[ip]) if ip is not None else ""
        rec["minL_strict_CPLB"] = int(c.L[isx]) if isx is not None else "not reached (>%d)" % c.L.max()
        rec["N_at_minL_point"] = int(c.N[ip]) if ip is not None else np.nan
        MIN.append(rec)
MIN = pd.DataFrame(MIN)
V = V.merge(MIN.pivot(index="est", columns="fam", values="minL_point").add_prefix("minL_point_").reset_index(), on="est", how="left")
V = V.merge(MIN.pivot(index="est", columns="fam", values="minL_strict_CPLB").add_prefix("minL_strict_").reset_index(), on="est", how="left")
V["PASS"] = V.C_pass.fillna(False) & V.F_pass.fillna(False)
V.to_csv(os.path.join(RES, "verdicts.csv"), index=False)
MIN.to_csv(os.path.join(RES, "min_detectable_L.csv"), index=False)

# power table incl. not-realised cells
pt = cells[(cells.kind == "crossover") & (cells.note == "")][["fam", "L", "N", "real", "est", "power", "power_lo", "power_hi", "abstain", "decl_zero_rate"]]
pt.to_csv(os.path.join(RES, "power_table.csv"), index=False)

# ablation table
ab = cells[cells.kind == "ablation"][["fam", "L", "note", "est", "abstain", "cov_simult", "cov_point", "bias", "mae", "decl_rate"]]
ab.to_csv(os.path.join(RES, "ablation_parity.csv"), index=False)
# descriptive: gasket, rrg
ds_ = cells[cells.kind == "descriptive"][["fam", "L", "N", "est", "abstain", "cov_simult", "cov_point", "cov_pop_simult", "bias", "mae"]]
ds_.to_csv(os.path.join(RES, "descriptive_gasket_rrg.csv"), index=False)

# convergence (d): first 100 vs last 100
cv = []
for est in EST:
    g = ctrl[(ctrl.est == est) & (ctrl.abstain == 0) & ctrl.cov_simult.notna()]
    a, b = g[g.ds < 100].cov_simult, g[g.ds >= 100].cov_simult
    if len(a) and len(b):
        pa, pb = a.mean(), b.mean(); se = np.sqrt(pa * (1 - pa) / len(a) + pb * (1 - pb) / len(b) + 1e-12)
        cv.append(dict(est=est, first100=pa, last100=pb, z=(pa - pb) / se))
pd.DataFrame(cv).to_csv(os.path.join(RES, "convergence_halves.csv"), index=False)

# ---------------- plots
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
# fig1 coverage vs N (controls)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
for a, fam in zip(ax, ("ring", "torus2", "torus3")):
    for est in EST:
        c = cells[(cells.fam == fam) & (cells.est == est) & (cells.note == "") & cells.cov_simult.notna() & (cells.abstain < 0.5)].sort_values("N")
        if len(c):
            a.errorbar(c.N, c.cov_simult, yerr=[c.cov_simult - c.cov_lo, c.cov_hi - c.cov_simult], marker="o", ms=3, lw=1, capsize=2, color=COL[est], label=NAMES[est])
    a.axhline(0.93, color="k", ls="--", lw=0.8); a.axhline(0.95, color="k", ls=":", lw=0.6)
    a.set_xscale("log"); a.set_title(f"{fam} (controls)"); a.set_xlabel("N nodes"); a.set_ylim(0, 1.03)
ax[0].set_ylabel("simultaneous coverage of exact d(t)\n(cells with abstention<50% only)"); ax[0].legend(fontsize=6.5, loc="lower left")
fig.suptitle("Fig 1. Simultaneous 95% band coverage of the exact d(t) on no-crossover controls (dashed = 93% bar)"); fig.tight_layout(); fig.savefig(os.path.join(RES, "fig1_coverage_controls.png"), dpi=140); plt.close(fig)

# fig2 power vs L
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
for a, fam in zip(ax, ("aniso2_w16", "aniso3_w8")):
    for est in EST:
        c = cells[(cells.fam == fam) & (cells.est == est) & (cells.note == "") & cells.power.notna()].sort_values("L")
        if not len(c): continue
        a.errorbar(c.L, c.power, yerr=[c.power - c.power_lo, c.power_hi - c.power], marker="o", ms=3, lw=1, capsize=2, color=COL[est], label=NAMES[est])
    nr = meta[(meta.fam == fam) & (meta.real == False) & (meta.note == "")]
    for L_ in nr.L: a.axvspan(L_ * 0.93, L_ * 1.07, color="0.9")
    a.axhline(0.8, color="k", ls="--", lw=0.8); a.set_xscale("log"); a.set_xlabel("long-axis length L (grey = crossover not realised in the exact window)")
    a.set_title({"aniso2_w16": "C_L x C_16  (2D -> 1D)", "aniso3_w8": "C_L x C_8 x C_8  (3D -> 1D)"}[fam]); a.set_ylim(-0.02, 1.03)
ax[0].set_ylabel("power (correct-sign change declared; abstention = miss)"); ax[0].legend(fontsize=6.5, loc="lower right")
fig.suptitle("Fig 2. Power to declare a dimensional change vs size (dashed = 0.8), Clopper-Pearson 95% bars"); fig.tight_layout(); fig.savefig(os.path.join(RES, "fig2_power_vs_L.png"), dpi=140); plt.close(fig)

# fig3 FPR
fig, ax = plt.subplots(figsize=(12, 3.8))
cc = cells[(cells.kind == "control") & (cells.note == "") & cells.decl_rate.notna() & cells.est.isin(["E0", "E1", "E1x", "E2", "E2h", "E3"])]
lab = cc[["fam", "L"]].drop_duplicates().sort_values(["fam", "L"]).reset_index(drop=True)
w = 0.13
for i, est in enumerate(["E0", "E1", "E1x", "E2", "E2h", "E3"]):
    x, y, lo, hi = [], [], [], []
    for j, r in lab.iterrows():
        c = cc[(cc.fam == r.fam) & (cc.L == r.L) & (cc.est == est)]
        if len(c) and c.abstain.iloc[0] < 0.9:
            x.append(j + (i - 2.5) * w); y.append(c.decl_rate.iloc[0]); lo.append(c.decl_rate.iloc[0] - c.decl_lo.iloc[0]); hi.append(c.decl_hi.iloc[0] - c.decl_rate.iloc[0])
    ax.errorbar(x, y, yerr=[lo, hi], fmt="o", ms=3, capsize=1.5, color=COL[est], label=NAMES[est])
ax.axhline(0.05, color="k", ls="--", lw=0.8); ax.set_xticks(range(len(lab))); ax.set_xticklabels([f"{a} {b}" for a, b in zip(lab.fam, lab.L)], rotation=60, ha="right")
ax.set_ylabel("false-positive rate (|Delta| band outside +-0.25)"); ax.legend(fontsize=6.5, ncol=3); ax.set_title("Fig 3. False-positive rate of declaring a dimensional change on controls (cells with abstention >= 90% omitted)")
fig.tight_layout(); fig.savefig(os.path.join(RES, "fig3_fpr_controls.png"), dpi=140); plt.close(fig)

# fig4 noise vs bias: coverage of exact d vs of own population functional, and leakage
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
for a, fam in zip(ax, ("torus2", "aniso2_w16")):
    for est, mk in (("E1", "o"), ("E1x", "s"), ("E2", "^")):
        c = cells[(cells.fam == fam) & (cells.est == est) & (cells.note == "") & cells.cov_simult.notna() & (cells.abstain < 0.5)].sort_values("N")
        if len(c):
            a.plot(c.N, c.cov_simult, marker=mk, color=COL[est], lw=1, label=f"{est} vs exact d(t)")
            a.plot(c.N, c.cov_pop_simult, marker=mk, color=COL[est], lw=1, ls=":", mfc="none", label=f"{est} vs own population functional")
    a.axhline(0.93, color="k", ls="--", lw=0.8); a.set_xscale("log"); a.set_ylim(0, 1.03); a.set_xlabel("N nodes"); a.set_title(fam)
ax[0].set_ylabel("simultaneous coverage"); ax[0].legend(fontsize=6.5)
fig.suptitle("Fig 4. Smoothing/selection bias vs noise: coverage of exact d(t) (solid) vs of the estimator's own limit (dotted); E1x = no sample split")
fig.tight_layout(); fig.savefig(os.path.join(RES, "fig4_bias_vs_noise.png"), dpi=140); plt.close(fig)

pd.set_option("display.width", 250, "display.max_columns", 50, "display.max_rows", 500)
print(V.to_string()); print(MIN.to_string()); print(pd.read_csv(os.path.join(RES, "convergence_halves.csv")).to_string())
