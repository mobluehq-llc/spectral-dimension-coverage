"""ILLUSTRATION ONLY (v7 Fig. 1a): one dataset of the cell C_256 x C_16 at 8 and at 128 walkers per node, chosen by a fixed rule
(dataset index 0, benchmark seed rule), with the LocQuad (archive code E1q) and BiasAware-CV (E4c) simultaneous bands.
Not part of any scored run; no number from it is reported. Same construction as make_fig1_data.py with the budget changed.
Run from the p22 folder: .venv/bin/python paper/make_fig1_budget_data.py"""
import sys, os
sys.path.insert(0, os.getcwd())
import numpy as np
import p22lib as L, p22amend as A, p22amend6 as A6, run_p22 as R
specs = [s for s in R.cell_specs() if s["fam"] == "aniso2_w16" and s["L"] == 256]
c = specs[0]; R.prep(c, np.random.default_rng(2200))
cap = {}; _orig = A.score
def score(dhat, Dstar, dtrue, tag, dpop=None):
    se = np.maximum(Dstar.std(0, ddof=1), 1e-9); cc = np.quantile((np.abs(Dstar - dhat) / se).max(1), 0.95)
    cap["dhat"] = dhat; cap["lo"] = dhat - cc * se; cap["hi"] = dhat + cc * se
    return _orig(dhat, Dstar, dtrue, tag, dpop)
A.score = score
pr = lambda s: np.array([float(v) for v in s.split(";")])
truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"]); ts, N = c["ts_w"], c["N"]
out = {}
for wpn in (8, 128):
    ds = 0
    ss = np.random.SeedSequence([22, c["id"], ds, wpn]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch[:6]]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = [C / m]
    r1 = A.est_E1b(U, m, ts, N, np.random.default_rng(ch[8]), truth, 2, 1.0, 200, c["thr"], "E1q")
    J = r1["J"]; out[f"t{wpn}"] = ts[:J].astype(float); out[f"q_lo{wpn}"] = cap["lo"].copy(); out[f"q_hi{wpn}"] = cap["hi"].copy()
    rows = A6.est_all(U, m, ts, N, rngs[2], rngs[3], rngs[4], truth, c["thr"], ("E4c",))
    e4 = [r for r in rows if r["est"] == "E4c"][0]
    out[f"b_lo{wpn}"] = pr(e4["e4lo"]); out[f"b_hi{wpn}"] = pr(e4["e4hi"])
    d = truth["d_exact"](ts[:J]); print(wpn, "J", J, "LocQuad covers:", bool(((d >= cap["lo"]) & (d <= cap["hi"])).all()), "len e4", len(out[f"b_lo{wpn}"]))
tg = np.logspace(np.log10(ts[0]), np.log10(ts[-1]), 400)
np.savez("paper/figs/fig1_budget_data.npz", tgrid=tg, dgrid=truth["d_exact"](tg), **out)
