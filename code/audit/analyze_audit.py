"""Tables and figures for the real-data audit. usage: analyze_audit.py results_dir out_dir"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import audit_lib as AL, andi_lib as AN, dfa_lib as D
from p22lib import cp

RES, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
EST = ["E0", "E1q", "E4c", "E6"]; COL = dict(E0="#8c8c8c", E1q="#d95f02", E4c="#1b4f9c", E6="#1b9e77")
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130})
md = []


def pct(k, n):
    if n == 0: return "n/a"
    lo, hi = cp(int(k), int(n)); return f"{100*k/n:.1f} ({100*lo:.1f}-{100*hi:.1f}) n={int(n)}"


# =================================================================== Part A
a = pd.read_csv(os.path.join(RES, "andi_raw.csv"))
a = a[a.est.isin(EST)].copy()
a["decl_any"] = ((a.decl_neg == 1) | (a.decl_pos == 1)).astype(float)
a["decl_right"] = np.where(a.sign_lab > 0, a.decl_pos, a.decl_neg)
SS = a[a.cls == "SS"]; SW = a[a.cls == "SW"]
rows = []
for T in (100, 150, 200):
    for e in EST:
        s = SS[(SS["T"] == T) & (SS.est == e)]; n_all = len(s); s_na = s[s.abstain == 0]
        w = SW[(SW["T"] == T) & (SW.est == e)]; w_na = w[w.abstain == 0]
        r = dict(T=T, est=e, J=int(a[(a["T"] == T) & (a.est == "E4c")].J.median()), n_declare_any=0, SS_n=n_all, SS_nonabstain=len(s_na),
                 SS_cov_k=int(s_na.cov_simult.sum()) if len(s_na) else 0, SS_fc_k=int(s_na.decl_any.sum()) if len(s_na) else 0,
                 SS_hw_med=s_na.hw.median() if len(s_na) else np.nan,
                 SW_n=len(w), SW_any_k=int(w.decl_any.sum()) if len(w) else 0, SW_nonabstain=len(w_na), SW_S_med=w.S_lab.median() if len(w) else np.nan, SW_power_k=int(w.decl_right.fillna(0).sum()) if len(w) else 0,
                 SW_cov_k=int(w_na.cov_simult.sum()) if len(w_na) else 0, SW_hw_med=w_na.hw.median() if len(w_na) else np.nan)
        rows.append(r)
TA = pd.DataFrame(rows); TA.to_csv(os.path.join(OUT, "andi_summary.csv"), index=False)
md.append("## Part A tables (AnDi-2, simulated labelled tracks)\n")
md.append("| T | est | J | SS n (non-abst) | SS coverage of true alpha % (CP) | SS false-change % (CP) | median hw (SS) | C / F |\n|---|---|---|---|---|---|---|---|")
for _, r in TA.iterrows():
    n = r.SS_nonabstain
    if n:
        lo, hi = cp(int(r.SS_cov_k), int(n)); flo, fhi = cp(int(r.SS_fc_k), int(n))
        C = "pass" if (r.SS_cov_k / n >= .93 and hi >= .93) else "FAIL"; F = "pass" if (r.SS_fc_k / n <= .05 and flo <= .05) else "FAIL"
    else: C = F = "n/a"
    md.append(f"| {r['T']} | {r.est} | {r.J} | {r.SS_n} ({n}) | {pct(r.SS_cov_k, n)} | {pct(r.SS_fc_k, n)} | {r.SS_hw_med:.2f} | {C} / {F} |")
md.append("\n| T | est | SW n (non-abst) | median S_lab | any-sign declaration % (CP) | power (correct sign, abstain=miss) % (CP) | coverage of label-implied d_lab % | median hw (SW) |\n|---|---|---|---|---|---|---|---|")
for _, r in TA.iterrows():
    md.append(f"| {r['T']} | {r.est} | {r.SW_n} ({r.SW_nonabstain}) | {r.SW_S_med:.2f} | {pct(r.SW_any_k, r.SW_n)} | {pct(r.SW_power_k, r.SW_n)} | {pct(r.SW_cov_k, r.SW_nonabstain)} | {r.SW_hw_med:.2f} |")

# recipe comparison
md.append("\n## Minimum-size recipe vs AnDi (SW strata)\n")
md.append("| T | est | J | S (median S_lab) | recipe hw_law | measured median hw (SW) | hw ratio meas/law | recipe predicts detect | measured power | match (power>=0.8 vs predicted) | J needed for S=0.5 -> T needed |\n|---|---|---|---|---|---|---|---|---|---|---|")
rec = []
for T in (100, 150, 200):
    for e in EST:
        r = TA[(TA["T"] == T) & (TA.est == e)].iloc[0]
        if not r.SW_n: continue
        J = int(r.J); rv = AL.recipe_view(e, J, float(r.SW_S_med)); rv5 = AL.recipe_view(e, J, 0.5)
        pw = r.SW_power_k / r.SW_n
        Tn = None
        if rv5["J_needed"] is not None:
            for TT in range(24, 20000, 4):
                if len(AN.lag_grid(TT)) >= rv5["J_needed"]: Tn = TT; break
        hwr = r.SW_hw_med / rv["hw_law"] if rv["hw_law"] == rv["hw_law"] else np.nan
        match = "n/a" if rv["detect"] is None else ("yes" if (pw >= .8) == rv["detect"] else "NO")
        rec.append(dict(T=T, est=e, J=J, S=r.SW_S_med, hw_law=rv["hw_law"], hw_meas=r.SW_hw_med, ratio=hwr, pred=rv["detect"], power=pw, match=match, J_needed=rv5["J_needed"], T_needed=Tn))
        md.append(f"| {T} | {e} | {J} | {r.SW_S_med:.2f} | {rv['hw_law']:.2f} | {r.SW_hw_med:.2f} | {hwr:.2f} | {rv['detect']} | {pw:.2f} | {match} | {rv5['J_needed']} -> T~{Tn if Tn else '>20000'} |")
pd.DataFrame(rec).to_csv(os.path.join(OUT, "andi_recipe.csv"), index=False)

# figure A1: coverage and false-change
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
for k, (col, ttl, ref) in enumerate((("cov", "Coverage of the true alpha (simultaneous, 95% nominal)", .93), ("fc", "False-change rate (any-sign, margin 0.25)", .05))):
    for i, e in enumerate(EST):
        for j, T in enumerate((100, 150, 200)):
            r = TA[(TA["T"] == T) & (TA.est == e)].iloc[0]; n = r.SS_nonabstain
            if not n: continue
            kk = r.SS_cov_k if col == "cov" else r.SS_fc_k; lo, hi = cp(int(kk), int(n)); p = kk / n
            x = j + (i - 1.5) * 0.18
            ax[k].errorbar(x, p, yerr=[[p - lo], [hi - p]], fmt="o", color=COL[e], ms=4, capsize=2, label=e if j == 0 else None)
    ax[k].axhline(ref, color="k", ls="--", lw=0.8); ax[k].set_xticks(range(3)); ax[k].set_xticklabels(["T=100 (J=6)", "T=150 (J=7)", "T=200 (J=8)"]); ax[k].set_ylabel("proportion of single-state tracks"); ax[k].set_title(ttl, fontsize=9)
    ax[k].set_ylim(-0.02, 1.02)
ax[0].legend(frameon=False, loc="lower left", fontsize=8)
fig.suptitle("AnDi-2 simulated single-state tracks: the bands vs the label (CP 95% bars; dashed = p22 pass line)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figA1_andi_coverage_fpr.png")); plt.close(fig)

# figure A2: power on SW vs T with recipe markers
fig, ax = plt.subplots(figsize=(6.2, 4.0))
for i, e in enumerate(EST):
    for j, T in enumerate((100, 150, 200)):
        r = TA[(TA["T"] == T) & (TA.est == e)].iloc[0]
        if not r.SW_n: continue
        p = r.SW_power_k / r.SW_n; lo, hi = cp(int(r.SW_power_k), int(r.SW_n)); x = j + (i - 1.5) * 0.18
        ax.errorbar(x, p, yerr=[[p - lo], [hi - p]], fmt="o", color=COL[e], ms=5, capsize=2, label=e if j == 0 else None)
        rr = [q for q in rec if q["T"] == T and q["est"] == e]
        if rr and rr[0]["pred"] is not None: ax.text(x, 1.06, "Y" if rr[0]["pred"] else "n", ha="center", color=COL[e], fontsize=8)
ax.axhline(0.8, color="k", ls="--", lw=0.8); ax.set_ylim(-0.02, 1.15); ax.set_xticks(range(3)); ax.set_xticklabels(["T=100", "T=150", "T=200"]); ax.set_ylabel("power on switching tracks (S_lab >= 0.3)")
ax.legend(frameon=False, fontsize=8, loc="center left"); ax.set_title("Power vs track length; letters on top = recipe's predicted detect (Y) / not (n)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figA2_andi_power_vs_T.png")); plt.close(fig)

# figure A3: hw vs J with recipe law
fig, ax = plt.subplots(figsize=(6.2, 4.0)); Jg = np.linspace(6, 16, 50)
P = AL.recipe_params()
for e in EST:
    ax.plot(Jg, P[e].a * (Jg / 7) ** (-P[e].delta), color=COL[e], lw=1, ls="--")
    pts = TA[TA.est == e]
    ax.plot(pts.J, pts.SS_hw_med, "o", color=COL[e], label=e)
ax.set_xlabel("window points J (set by track length T)"); ax.set_ylabel("median band half-width (alpha units)"); ax.set_yscale("log")
ax.legend(frameon=False, fontsize=8); ax.set_title("Measured half-width on AnDi tracks (dots) vs the recipe's law a (J/7)^-delta (dashed)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figA3_andi_hw_vs_recipe.png")); plt.close(fig)

# =================================================================== Part B
real = pd.read_csv(os.path.join(RES, "cet_real.csv")); cur = json.load(open(os.path.join(RES, "cet_curves.json"))); meta = json.load(open(os.path.join(RES, "cet_meta.json")))
md.append("\n## Part B tables (HadCET)\n")
md.append(f"Record: {meta['n_days']} days {meta['first']} .. {meta['last']}, sha256 {meta['sha256'][:16]}...; training-residual sd {meta['train_resid_sd']:.2f} C.\n")
md.append("| period | est | margin | hw (H units) | declared neg / pos | contrast min..max (point) | E4c contrast bounds (min lower, max upper) | E6 location set (days) best | E6 step | E6 misfit |\n|---|---|---|---|---|---|---|---|---|---|")
for _, r in real.sort_values(["period", "est", "delta"]).iterrows():
    if r.est not in EST: continue
    loc = f"{np.exp(r.e6_m_lo):.1f}-{np.exp(r.e6_m_hi):.1f} ({np.exp(r.e6_m_best):.1f})" if r.est == "E6" and r.e6_m_lo == r.e6_m_lo else ""
    st = f"{r.e6_step_dB:+.2f}" if r.est == "E6" and r.e6_step_dB == r.e6_step_dB else ""
    bnd = f"{r.dlo_min:+.2f}, {r.dhi_max:+.2f}" if r.est == "E4c" else ""
    mf = ("ALARM" if r.lof_alarm == 1 else "ok") if r.est == "E6" else ""
    md.append(f"| {r.period} | {r.est} | {r.delta} | {r.hw:.2f} | {int(r.decl_neg)} / {int(r.decl_pos)} | {r.dmin:+.2f}..{r.dmax:+.2f} | {bnd} | {loc} | {st} | {mf} |")
g = lambda per, est, dl, k: real[(real.period == per) & (real.est == est) & (real.delta == dl)].iloc[0][k]


def verdict_for(dm):
    d4 = {p: (int(g(p, "E4c", dm, "decl_neg")), int(g(p, "E4c", dm, "decl_pos"))) for p in ("P1", "P2")}
    eq = {p: (g(p, "E4c", dm, "dlo_min") > -dm and g(p, "E4c", dm, "dhi_max") < dm) for p in ("P1", "P2")}
    loc = {p: (np.exp(g(p, "E6", dm, "e6_m_lo")), np.exp(g(p, "E6", dm, "e6_m_hi"))) for p in ("P1", "P2")}
    hit = {p: (loc[p][0] <= 20 and loc[p][1] >= 5) for p in loc}
    sup = d4["P1"][0] == 1 and d4["P2"][0] == 1 and hit["P1"]
    nsa = eq["P1"] and eq["P2"]
    nsb = d4["P1"][0] == 1 and d4["P2"][0] == 1 and (not hit["P1"]) and (not hit["P2"])
    V = "SUPPORTED" if sup else ("NOT SUPPORTED" if (nsa or nsb) else "UNDECIDABLE")
    txt = (f"margin {dm}: E4c declared-neg P1/P2 = {d4['P1'][0]}/{d4['P2'][0]} (pos {d4['P1'][1]}/{d4['P2'][1]}); equivalence inside +-{dm} P1/P2 = {eq['P1']}/{eq['P2']}; "
           f"E6 location set P1 {loc['P1'][0]:.1f}-{loc['P1'][1]:.1f} d, P2 {loc['P2'][0]:.1f}-{loc['P2'][1]:.1f} d (intersects 5-20 d: {hit['P1']}/{hit['P2']}) -> {V}")
    return V, txt


V25, t25 = verdict_for(0.25); V10, t10 = verdict_for(0.10)
md.append(f"\n**PRIMARY verdict (Amendment 1, margin 0.25): {V25}.** {t25}\n\nPre-registered original margin, confounded by DFA small-scale bias (null drift 0.17): {V10}. {t10}")
V = V25
# null-referenced contrast
pop0 = json.load(open(os.path.join(RES, "cet_pop.json")))["N0"]["d_coarse"]
for p in ("P1", "P2"):
    r = real[(real.period == p) & (real.est == "E4c") & (real.delta == 0.10)].iloc[0]; d = np.array(r.d.split(";"), float)
    md.append(f"\nE4c point slopes {p}: {np.round(d,2).tolist()}; minus N0 population slopes: {np.round(d-np.array(pop0),2).tolist()}")
# recipe view
md.append("\n## Recipe view (CET, J = 15 = extrapolation beyond fitted J 7-11)\n| est | hw_law | measured hw P1 / P2 | detect S=0.30 (margin 0.25)? | J needed for S=0.30 | smallest certifiable S at J=15 |\n|---|---|---|---|---|---|")
for e in ("E4c", "E6", "E1q", "E0"):
    rv = AL.recipe_view(e, 15, 0.30)
    md.append(f"| {e} | {rv['hw_law']:.2f} | {g('P1', e, 0.25, 'hw'):.2f} / {g('P2', e, 0.25, 'hw'):.2f} | {rv['detect']} (ratio {rv['ratio']:.2f} vs T {rv['T']:.2f}) | {rv['J_needed']} | {rv['S_min_at_J']:.2f} |")

md.append("\nRecipe with the MEASURED half-width (E4c/E6 at margin 0.25) and the observed contrast S_obs = max |point contrast| (E4c, E6):\n\n| period | est | measured hw | S_obs | ratio (S_obs-0.25)/hw | recipe T | recipe says detect | ratio for S_claim=0.30 | recipe says detect S_claim |\n|---|---|---|---|---|---|---|---|---|")
for p in ("P1", "P2"):
    for e in ("E4c", "E6"):
        r = real[(real.period == p) & (real.est == e) & (real.delta == 0.25)].iloc[0]; Tt = AL.recipe_params()[e]["T"]; So = max(abs(r.dmin), abs(r.dmax)); hw = r.hw
        md.append(f"| {p} | {e} | {hw:.2f} | {So:.2f} | {(So-0.25)/hw:.2f} | {Tt:.2f} | {(So-0.25)/hw >= Tt} | {(0.30-0.25)/hw:.2f} | {(0.30-0.25)/hw >= Tt} |")
# figure B1
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
for p, c in (("P1", "#1b4f9c"), ("P2", "#d95f02")):
    s = np.array(cur[p]["scales"]); F = np.array(cur[p]["F"]); sd = np.array(cur[p]["logF_sd"])
    ax[0].errorbar(s, F / s ** 0.5, yerr=[F / s ** 0.5 * (1 - np.exp(-1.96 * sd)), F / s ** 0.5 * (np.exp(1.96 * sd) - 1)], fmt="o-", color=c, ms=3, lw=1, label={"P1": "1851-1950", "P2": "1951-2025"}[p])
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("scale s (days)"); ax[0].set_ylabel("DFA-2 fluctuation F(s) / s^0.5"); ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title("Compensated DFA curve (flat = H 0.5)", fontsize=9)
for k, (p, c) in enumerate((("P1", "#1b4f9c"), ("P2", "#d95f02"))):
    for e, ls in (("E4c", "-"), ("E6", ":")):
        r = real[(real.period == p) & (real.est == e) & (real.delta == 0.10)].iloc[0]
        d = np.array(r.d.split(";"), float); lo = np.array(r.lo.split(";"), float); hi = np.array(r.hi.split(";"), float); s = np.array(cur[p]["scales"])
        if e == "E4c":
            ax[1].fill_between(s, lo, hi, color=c, alpha=0.18); ax[1].plot(s, d, color=c, lw=1.5, label=f"{p} E4c")
        else:
            ax[1].plot(s, lo, color=c, ls=":", lw=1.2); ax[1].plot(s, hi, color=c, ls=":", lw=1.2, label=f"{p} E6 band")
ax[1].axvspan(5, 20, color="0.85", zorder=0); ax[1].set_xscale("log"); ax[1].set_xlabel("scale s (days)  [grey = claimed crossover 5-20 d]"); ax[1].set_ylabel("local DFA slope H(s)")
ax[1].legend(frameon=False, fontsize=7, ncol=2); ax[1].set_title("Local slope: E4c honest band (shaded), E6 band (dotted)", fontsize=9)
# panel 3: contrast d_j - d_0 with E4c bounds approx via E4c band? show E6 location sets
for k, (p, c) in enumerate((("P1", "#1b4f9c"), ("P2", "#d95f02"))):
    r = real[(real.period == p) & (real.est == "E6") & (real.delta == 0.10)].iloc[0]
    ax[2].plot([np.exp(r.e6_m_lo), np.exp(r.e6_m_hi)], [k, k], color=c, lw=5, solid_capstyle="butt"); ax[2].plot(np.exp(r.e6_m_best), k, "k|", ms=14)
    ax[2].text(np.exp(r.e6_m_hi) * 1.15, k, f"step {r.e6_step_dB:+.2f}" + (", misfit alarm" if r.lof_alarm == 1 else ""), va="center", fontsize=8)
ax[2].axvspan(5, 20, color="0.85", zorder=0); ax[2].set_xscale("log"); ax[2].set_yticks([0, 1]); ax[2].set_yticklabels(["P1", "P2"]); ax[2].set_ylim(-0.6, 1.6); ax[2].set_xlabel("E6 logistic step location (days); lower end = edge of the search grid (2.6 d)"); ax[2].set_xticks([3, 5, 10, 20]); ax[2].set_xticklabels(["3", "5", "10", "20"]); ax[2].minorticks_off()
ax[2].set_title("Where E6 puts the crossover (99% set; tick = best fit)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figB1_cet_dfa.png")); plt.close(fig)

# synthetic calibration
sy = pd.read_csv(os.path.join(RES, "cet_synth.csv")); sy = sy[sy.est.isin(EST)].copy()
sy["decl_any"] = ((sy.decl_neg == 1) | (sy.decl_pos == 1)).astype(float)
pop = json.load(open(os.path.join(RES, "cet_pop.json")))
md.append("\n## Pipeline calibration on synthetic records with the CET layout (100 y, 100 replicates)")
for kind in ("N0", "N1"):
    d0 = np.array(pop[kind]["d_coarse"]); md.append(f"\n({kind} population local slope, 15 scales: {np.round(d0,2).tolist()}; max|d_j-d_0| = {np.abs(d0-d0[0]).max():.2f})\n\n| kind | est | margin | false-change / detection % (CP) | simultaneous coverage of population slope % (CP) | median hw |\n|---|---|---|---|---|---|")
    for e in EST:
        for dl in (0.10, 0.25):
            s = sy[(sy.kind == kind) & (sy.est == e) & (sy.delta == dl)]
            cv = s.cov_simult.dropna()
            what = "neg" if kind == "N1" else "any"
            k = s.decl_any.sum() if kind == "N0" else s.decl_neg.sum()
            md.append(f"| {kind} | {e} | {dl} | {pct(k, len(s))} ({'false-change' if kind=='N0' else 'neg declaration'}) | {pct(cv.sum(), len(cv))} | {s.hw.median():.2f} |")
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for j, (ttl, kind, col) in enumerate((("N0 (fGn, no crossover): false-change rate", "N0", "fc"), ("N0: coverage of population slope", "N0", "cov"), ("N1 (AR1+fGn, crossover): negative-declaration rate", "N1", "det"))):
    for i, e in enumerate(EST):
        for q, dl in enumerate((0.10, 0.25)):
            s = sy[(sy.kind == kind) & (sy.est == e) & (sy.delta == dl)]
            if col == "cov": s = s[s.cov_simult.notna()]; k = s.cov_simult.sum()
            elif col == "fc": k = s.decl_any.sum()
            else: k = s.decl_neg.sum()
            n = len(s)
            if not n: continue
            p = k / n; lo, hi = cp(int(k), int(n)); x = q + (i - 1.5) * 0.18
            ax[j].errorbar(x, p, yerr=[[p - lo], [hi - p]], fmt="o", color=COL[e], ms=4, capsize=2, label=e if q == 0 else None)
    ax[j].set_xticks([0, 1]); ax[j].set_xticklabels(["margin 0.10", "margin 0.25"]); ax[j].set_ylim(-0.02, 1.02); ax[j].set_title(ttl, fontsize=9)
    ax[j].axhline({"fc": .05, "cov": .93, "det": .8}[col], color="k", ls="--", lw=0.8)
ax[0].legend(frameon=False, fontsize=8); ax[0].set_ylabel("proportion of 100 synthetic records")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figB2_cet_synthetic_calibration.png")); plt.close(fig)
open(os.path.join(OUT, "audit_tables.md"), "w").write("\n".join(md)); print("\n".join(md))
