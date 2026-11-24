#!/usr/bin/env python3
"""P22 convergence / sanity checks (PLAN.md section 'Convergence / sanity checks'). Writes results/checks.csv"""
import os
import numpy as np, pandas as pd
import p22lib as L

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
def rec(name, ok, detail):
    rows.append(dict(check=name, ok=bool(ok), detail=detail)); print(("PASS " if ok else "FAIL ") + name, detail, flush=True)

# (a) walker samplers vs exact q
def ring_step(N, lazy=0.5):
    nbr = np.stack([(np.arange(N) - 1) % N, (np.arange(N) + 1) % N], 1)
    return L.Graph("step", lazy, nbr=nbr, deg=np.full(N, 2))
ts = L.make_grid(24, 34, 200)
for name, g, lazy_ok in (("direct torus 24x24", L.Graph("torus", 0.5, dims=(24, 24)), True), ("direct ring 64", L.Graph("torus", 0.5, dims=(64,)), True),
                         ("direct 3D torus 8^3", L.Graph("torus", 0.5, dims=(8, 8, 8)), True), ("step ring 64", ring_step(64), True)):
    rng = np.random.default_rng(5)
    W = 400000
    C, m = g.sample_blocks(rng, W, ts, 200)
    tot = C.sum(0); Wn = C.shape[0] * m
    q, _ = g.exact(ts)
    z = (tot - Wn * q) / np.sqrt(Wn * q * (1 - q))
    rec(f"(a) sampler {name}", np.abs(z).max() < 4.0, f"max|z|={np.abs(z).max():.2f} over {len(ts)} times, W={Wn}")
# step sampler on gasket n=4 vs exact
g = L.gasket(4); tsg = L.make_grid(24, 34, 200); rng = np.random.default_rng(6)
C, m = g.sample_blocks(rng, 400000, tsg, 200); tot = C.sum(0); Wn = C.shape[0] * m; q, _ = g.exact(tsg)
z = (tot - Wn * q) / np.sqrt(Wn * q * (1 - q))
rec("(a) sampler step gasket n=4", np.abs(z).max() < 4.0, f"max|z|={np.abs(z).max():.2f}")

# (b) probes unbiased
g = L.Graph("torus", 0.5, dims=(16, 16)); rng = np.random.default_rng(7)
tsp = L.make_grid(24, 34, 100)
X = g.sample_probes(rng, 4000, tsp); q, _ = g.exact(tsp)
z = (X.mean(0) - q) / (X.std(0, ddof=1) / np.sqrt(X.shape[0]))
rec("(b) Hutchinson probes unbiased torus 16x16", np.abs(z).max() < 4.0, f"max|z|={np.abs(z).max():.2f}")
g = L.gasket(5); X = g.sample_probes(rng, 4000, tsp); q, _ = g.exact(tsp)
z = (X.mean(0) - q) / (X.std(0, ddof=1) / np.sqrt(X.shape[0]))
rec("(b) Hutchinson probes unbiased gasket n=5", np.abs(z).max() < 4.0, f"max|z|={np.abs(z).max():.2f}")

# (e) spectral truth vs dense matrix powers
for name, g in (("ring 40", L.Graph("torus", 0.5, dims=(40,))), ("torus 12x10", L.Graph("torus", 0.5, dims=(12, 10))), ("gasket n=5", L.gasket(5))):
    S = g.ssym().toarray(); q, _ = g.exact(np.array([1, 5, 24, 100]))
    Pw = [np.trace(np.linalg.matrix_power(S, t)) / g.N for t in (1, 5, 24, 100)]
    err = np.abs(q - np.array(Pw)).max()
    rec(f"(e) spectral q vs matrix power {name}", err < 1e-10, f"max err {err:.2e}")
    # derivative vs finite difference
    t0 = np.array([30.0]); h = 1e-3
    qa, da = g.exact(t0); qp, _ = g.exact(t0 + h); qm, _ = g.exact(t0 - h)
    fd = -2 * t0 * (np.log(qp) - np.log(qm)) / (2 * h)
    rec(f"(e) d(t) vs finite difference {name}", abs(fd - da).max() < 1e-6, f"|diff|={abs(fd-da).max():.2e}")

# (c) bootstrap B convergence on ring 256
c_truth = None
g = L.Graph("torus", 0.5, dims=(256,)); tsr = L.make_grid(24, 34, 1500)
truth = dict(q_exact=lambda t: g.exact(t)[0], d_exact=lambda t: g.exact(t)[1])
ratios = []
for ds in range(20):
    rng = np.random.default_rng([22, 999, ds])
    C, m = g.sample_blocks(rng, 32 * 256, tsr, 200)
    a = L.est_E1([C / m], m, tsr, 256, np.random.default_rng(1), truth, True, 200, L.THR)
    b = L.est_E1([C / m], m, tsr, 256, np.random.default_rng(2), truth, True, 1000, L.THR)
    if not a.get("abstain"): ratios.append(a["se_mean"] / b["se_mean"])
ratios = np.array(ratios)
rec("(c) bootstrap sd B=200 vs B=1000", np.abs(ratios - 1).max() < 0.10, f"ratio range [{ratios.min():.3f},{ratios.max():.3f}] over {len(ratios)} datasets")

# (f) gasket exact d vs asymptotic 2 ln3/ln5
ds_inf = 2 * np.log(3) / np.log(5)
for n in (4, 5, 6, 7):
    g = L.gasket(n); tt = L.make_grid(24, 34, 3072); q, d = g.exact(tt); win = (g.N * q >= 6) & (tt >= 24)
    if win.sum() >= 3:
        rec(f"(f) gasket n={n} exact window d", True, f"N={g.N} window t={tt[win][0]}..{tt[win][-1]} mean d={d[win].mean():.3f} (asymptotic {ds_inf:.3f}), min {d[win].min():.3f} max {d[win].max():.3f}")
    else:
        rec(f"(f) gasket n={n} exact window d", True, f"N={g.N} window too short (<3 grid points)")
pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "checks.csv"), index=False)
