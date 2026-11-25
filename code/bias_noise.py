#!/usr/bin/env python3
"""P22 Amendment 14: bias-versus-noise decomposition of fixed-window E1q misses. Stored data only. See PLAN.md AMENDMENT 14."""
import os, sys, json, hashlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
ARGV = list(sys.argv)
sys.argv = [ARGV[0], "verify"]
import numpy as np, pandas as pd
import run_amend1011 as RA
import p22amend as A
sys.argv = ARGV
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "archive", "data", "raw_amend10c.csv")
OUT = os.path.join(HERE, "results")
H_QUAD = A.H_QUAD

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

raw = pd.read_csv(RAW)
print("raw sha256", sha(RAW))
e = raw[(raw.est == "E1q") & (raw["mode"] == "fixed")].copy()
meta = pd.read_csv(os.path.join(HERE, "results_amend10_11", "cells_meta_amend10c.csv"))
wins = json.load(open(os.path.join(HERE, "results_amend10_11", "amend10c_windows.json")))["cells"]
cc = RA.cells_10c()
spec = {}
for w in (8, 32, 128):
    for s in cc[w]:
        spec[s["id"]] = s
rows_h, rows_c, dropped = [], [], []
byid = {}
for cid, g in e.groupby("cell"):
    s = spec[cid]; byid[cid] = (s, g)
cells = {}
for cid, (s, g) in byid.items():
    cells.setdefault((s["fam"], int(s["L"])), {})[int(s["wpn"])] = (s, g)
lines = []
for (fam, L), d in sorted(cells.items()):
    ok = True; info = {}
    if set(d) != {8, 32, 128}: dropped.append((fam, L, "budgets", sorted(d))); continue
    for w, (s, g) in d.items():
        ga = g[g.abstain == 0]
        if len(ga) < 100: ok = False; dropped.append((fam, L, "answered<100", w, len(ga)))
    if not ok: continue
    J = int(wins[f"{fam}_{L}"]["exact_J_wpn32"])
    s0 = d[32][0]; ts = s0["ts_w"][:J].astype(float)
    x = np.log(ts); dt = s0["d_exact"](ts); q = s0["q_exact"](ts)
    # consistency with stored data
    ga32 = d[32][1][d[32][1].abstain == 0]
    if not (int(ga32.tlast.iloc[0]) == int(s0["ts_w"][J - 1]) and abs(ga32.true_dmin.iloc[0] - (dt[1:] - dt[0]).min()) < 1e-9):
        dropped.append((fam, L, "truth mismatch", int(ga32.tlast.iloc[0]), int(s0["ts_w"][J - 1]))); continue
    # deterministic bias per candidate h, population weights
    bh = {}
    for h in H_QUAD:
        _, sl = A.poly_rows(x, q * s0["N"], h, 2)
        bh[h] = -2 * sl @ np.log(q) - dt
        rows_h.append(dict(fam=fam, L=L, J=J, h=h, bias_mean=float(bh[h].mean()), bias_maxabs=float(np.abs(bh[h]).max())))
    for w in (8, 32, 128):
        s, g = d[w]; ga = g[g.abstain == 0]
        hh = ga.h.values
        bmean = np.array([bh[min(H_QUAD, key=lambda z: abs(z - v))].mean() for v in hh])
        bmax = np.array([np.abs(bh[min(H_QUAD, key=lambda z: abs(z - v))]).max() for v in hh])
        ident = (ga.bias - ga.bias_pop).values
        S = float(ga.bias_pop.std(ddof=1))
        hmod = pd.Series(hh).mode().iloc[0]
        bfrz = float(bh[hmod if False else 1.0].mean()) if False else None
        miss = int((ga.cov_simult == 0).sum())
        rows_c.append(dict(fam=fam, L=L, J=J, wpn=w, n_ans=len(ga), n_all=len(g), misses=miss, cov=float(ga.cov_simult.mean()),
                           S=S, se_mean=float(ga.se_mean.mean()), hw_med=float(ga.hw.median()),
                           bias_det=float(bmean.mean()), bias_det_absmean=float(np.abs(bmean).mean()), bias_max=float(bmax.mean()),
                           R1=abs(float(bmean.mean())) / S, R2=float(bmax.mean()) / float(ga.se_mean.mean()), R3=float(bmax.mean()) / float(ga.hw.median()),
                           ident_mean=float(ident.mean()), ident_vs_pop_medabsdiff=float(np.median(np.abs(ident - bmean))),
                           bias_pop_mean=float(ga.bias_pop.mean()), h_mean=float(hh.mean()), h_mode=float(hmod),
                           frac_h_le1=float((hh <= 1.0).mean())))
C = pd.DataFrame(rows_c); H = pd.DataFrame(rows_h)
os.makedirs(OUT, exist_ok=True)
C.to_csv(os.path.join(OUT, "bias_noise_cells.csv"), index=False); H.to_csv(os.path.join(OUT, "bias_noise_by_h.csv"), index=False)
# frozen-h partial arm: wpn-32 modal h held fixed
fz = []
for (fam, L), g in C.groupby(["fam", "L"]):
    g = g.set_index("wpn"); h32 = g.loc[32, "h_mode"]
    b = H[(H.fam == fam) & (H.L == L) & (np.isclose(H.h, h32))].bias_mean.iloc[0]
    for w in (8, 32, 128): fz.append(dict(fam=fam, L=L, wpn=w, h_frozen=h32, bias_frozen=b, R1_frozen_over_cvspread=abs(b) / g.loc[w, "S"]))
Z = pd.DataFrame(fz); Z.to_csv(os.path.join(OUT, "bias_noise_frozen_partial.csv"), index=False)
# rule
rise = []
for (fam, L), g in C.groupby(["fam", "L"]):
    g = g.set_index("wpn"); r = g.R1
    rise.append(dict(fam=fam, L=L, R8=r[8], R32=r[32], R128=r[128], rising=bool(r[8] < r[32] < r[128])))
Rr = pd.DataFrame(rise)
tot = C.misses.sum()
res = dict(n_cells=len(Rr), n_rising=int(Rr.rising.sum()), share_rising=float(Rr.rising.mean()), total_misses=int(tot),
           misses_R1_gt2=int(C[C.R1 > 2].misses.sum()), misses_R1_gt1=int(C[C.R1 > 1].misses.sum()),
           share_misses_R1_gt2=float(C[C.R1 > 2].misses.sum() / tot), share_misses_R1_gt1=float(C[C.R1 > 1].misses.sum() / tot),
           dropped=[str(z) for z in dropped], raw_sha256=sha(RAW), code_sha256=sha(__file__))
res["clause_i"] = res["share_rising"] >= 0.75; res["clause_ii"] = res["share_misses_R1_gt2"] >= 0.80
res["explained"] = bool(res["clause_i"] and res["clause_ii"])
json.dump(res, open(os.path.join(OUT, "bias_noise_rule.json"), "w"), indent=1, default=float)
Rr.to_csv(os.path.join(OUT, "bias_noise_rising.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(C[["fam", "L", "wpn", "n_ans", "misses", "cov", "S", "se_mean", "hw_med", "bias_det", "bias_max", "R1", "R2", "R3", "ident_mean", "ident_vs_pop_medabsdiff", "h_mean", "frac_h_le1"]].round(4).to_string())
print(Rr.round(3).to_string()); print(json.dumps(res, indent=1, default=float))
print(H.round(4).to_string())
print(Z.round(3).to_string())
