"""v8 evidence lines for the mechanism section (no simulation: stored per-dataset rows and exact spectra only).
Run from p22_spectral_dimension_flow:  .venv/bin/python paper/make_v8_numbers.py   -> results/v8_mechanism_numbers.txt
(1) LocQuad (E1q) at fixed windows, archive/data/raw_amend10c.csv: coverage of exact d and of its own smoothed target d_pop,
    bandwidth floor/ceiling shares, mean h.  (2) Noiseless LocQuad on C_256 x C_16 at h = 0.7 (frozen p22amend.poly_rows, weights
    proportional to the exact return probability, fixed window of 12 grid times) against the exact d(96).
(3) Log-convexity secant certificate on exact returns of C_768 x C_16."""
import sys, os
sys.path.insert(0, os.getcwd())
import numpy as np, pandas as pd
import p22amend as A
out = []
d = pd.read_csv("archive/data/raw_amend10c.csv", low_memory=False)
d = d[(d.est == "E1q") & (d["mode"] == "fixed") & (d.abstain == 0)]
for w in (8, 32, 128):
    g = d[d.wpn == w]; n = len(g)
    out.append(f"E1q fixed wpn {w}: answered {n} | cov exact {int(g.cov_simult.sum())}/{n} = {100 * g.cov_simult.mean():.1f}% | "
               f"cov own-target {int(g.cov_pop_simult.sum())}/{n} = {100 * g.cov_pop_simult.mean():.1f}% | "
               f"h at floor 0.7: {100 * (g.h == 0.7).mean():.1f}% | h at ceiling 4.0: {100 * (g.h == 4.0).mean():.1f}% | mean h {g.h.mean():.2f}")
def spec(Ls):
    lam = np.zeros(1)
    for L in Ls: lam = (lam[:, None] + np.cos(2 * np.pi * np.arange(L) / L)[None, :] / len(Ls)).ravel()
    return 0.5 + 0.5 * lam
def q_d(lam, t):
    t = np.atleast_1d(np.asarray(t, float))[:, None]; l = lam[lam > 0][None, :]; w = l ** t
    return w.mean(1), (-2 * t[:, 0] * (w * np.log(l)).sum(1) / w.sum(1))
lam = spec([256, 16]); ts = []; t_ = [24, 34]
while len(ts) < 12: ts += t_; t_ = [2 * x for x in t_]
ts = np.array(ts, float); q, dex = q_d(lam, ts)
_, sl = A.poly_rows(np.log(ts), q, 0.7, 2)
dpop = -2 * sl @ np.log(q); j = list(ts).index(96.0)
out.append(f"noiseless LocQuad C256xC16 h=0.7 t=96: d_pop {dpop[j]:.4f} | exact d {dex[j]:.4f} | bias {dpop[j] - dex[j]:.4f}")
lam = spec([768, 16]); tt = np.array([24, 34, 4352, 6144.0]); q, dex = q_d(lam, tt); lq = np.log(q)
lo24 = 2 * 24 * (lq[0] - lq[1]) / (34 - 24); hi6144 = 2 * 6144 * (lq[2] - lq[3]) / (6144 - 4352)
out.append(f"secant certificate C768xC16 exact returns: d(24) >= {lo24:.4f} | d(6144) <= {hi6144:.4f} | certified drop >= {lo24 - hi6144:.4f} | exact d(24) {dex[0]:.4f} d(6144) {dex[3]:.4f}")
for L in (96, 256, 768):
    _, dd = q_d(spec([L, 16]), [96, 163, 400]); out.append(f"exact d C{L}xC16 at t=96,163,400: {dd[0]:.3f} {dd[1]:.3f} {dd[2]:.3f}")
for L in (64, 192):
    _, dd = q_d(spec([L, 8, 8]), [34, 62, 192]); out.append(f"exact d C{L}xC8xC8 at t=34,62,192: {dd[0]:.3f} {dd[1]:.3f} {dd[2]:.3f}")
open("results/v8_mechanism_numbers.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
