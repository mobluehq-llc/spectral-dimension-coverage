"""Regenerate ONE illustrative dataset (cell aniso2_w16, L=256, ds=0, benchmark seeds) and save the exact d(t) and the
E1q / E4c / E6 simultaneous bands to paper/figs/fig1_data.npz. Run from the p22 folder with .venv/bin/python paper/make_fig1_data.py.
Bands for E4c/E6 come from p22amend6.est_all; E1q's band is captured from p22lib.score (simultaneous max-t band dhat +- c*se)."""
import sys, os
sys.path.insert(0, os.getcwd())
import numpy as np
import p22lib as L, p22amend as A, p22amend6 as A6, run_p22 as R

specs = [s for s in R.cell_specs() if s["fam"] == "aniso2_w16" and s["L"] == 256]
c = specs[0]; R.prep(c, np.random.default_rng(2200))
cap = {}
_orig = A.score
def score(dhat, Dstar, dtrue, tag, dpop=None):
    se = np.maximum(Dstar.std(0, ddof=1), 1e-9); cc = np.quantile((np.abs(Dstar - dhat) / se).max(1), 0.95)
    cap["dhat"] = dhat; cap["lo"] = dhat - cc * se; cap["hi"] = dhat + cc * se
    return _orig(dhat, Dstar, dtrue, tag, dpop)
A.score = score
ds = 0
ss = np.random.SeedSequence([22, c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch[:6]]
truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"]); ts, N = c["ts_w"], c["N"]
C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts, L.NBLOCK); U = [C / m]
r1 = A.est_E1b(U, m, ts, N, np.random.default_rng(ch[8]), truth, 2, 1.0, 200, c["thr"], "E1q")
rows = A6.est_all(U, m, ts, N, rngs[2], rngs[3], rngs[4], truth, c["thr"], ("E4c",))
pr = lambda s: np.array([float(v) for v in s.split(";")])
e6 = [r for r in rows if r["est"] == "E6"][0]; e4 = [r for r in rows if r["est"] == "E4c"][0]
J = r1["J"]; assert len(pr(e6["e6lo"])) == J, (J, len(pr(e6["e6lo"])))
tg = np.logspace(np.log10(ts[0]), np.log10(ts[-1]), 400)
np.savez("paper/figs/fig1_data.npz", t=ts[:J].astype(float), tgrid=tg, dgrid=truth["d_exact"](tg), dex=truth["d_exact"](ts[:J]),
         e1q_lo=cap["lo"], e1q_hi=cap["hi"], e1q=cap["dhat"], e4_lo=pr(e4["e4lo"]), e4_hi=pr(e4["e4hi"]),
         e6_lo=pr(e6["e6lo"]), e6_hi=pr(e6["e6hi"]), N=N, J=J)
for nm, lo, hi in (("E1q", cap["lo"], cap["hi"]), ("E4c", pr(e4["e4lo"]), pr(e4["e4hi"])), ("E6", pr(e6["e6lo"]), pr(e6["e6hi"]))):
    d = truth["d_exact"](ts[:J]); print(nm, "covers exact (simult.):", bool(((d >= lo) & (d <= hi)).all()), "mean hw %.2f" % ((hi - lo).mean() / 2))
print("J", J, "tlast", ts[J - 1])
