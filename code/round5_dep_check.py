#!/usr/bin/env python3
"""ROUND 5 priority check (ChatGPT r5 s8): are the return counts at different times drawn from the SAME walkers (correlated) or independently?
Ring C_128, lazy walk, times 10, 24, 34. Exact per-walker cov(ret24, ret34) = q(24)[q(10)-q(34)] (vertex-transitive graph).
Run from the directory holding p22lib.py (host-D: <local path omitted>). Output: empirical vs exact, per walker and per block of m walkers."""
import sys, numpy as np
sys.path.insert(0, ".")
import p22lib as L
g = L.Graph("torus", lazy=0.5, dims=(128,)); ts = np.array([10, 24, 34]); q, _ = g.exact(ts); ex = q[1] * (q[0] - q[2])
W = 2_000_000; cs = []; m24 = []; m34 = []
for seed in range(5, 45):
    ret = g._direct(np.random.default_rng(seed), W, ts).astype(float); a, b = ret[:, 1], ret[:, 2]
    cs.append(np.mean(a * b) - a.mean() * b.mean()); m24.append(a.mean()); m34.append(b.mean())
cs = np.array(cs); se = cs.std(ddof=1) / np.sqrt(len(cs))
print(f"marginals: q(24) exact {q[1]:.6f} empirical {np.mean(m24):.6f} (se {np.std(m24,ddof=1)/np.sqrt(len(m24)):.6f}); q(34) exact {q[2]:.6f} empirical {np.mean(m34):.6f} (se {np.std(m34,ddof=1)/np.sqrt(len(m34)):.6f})")
print(f"q(10,24,34) = {q}\nexact per-walker cov = {ex:.6f}\nempirical (40 x {W} walkers) = {cs.mean():.6f}  se {se:.6f}  (z = {(cs.mean()-ex)/se:.2f}); independent draws would give 0")
C, m = g.sample_blocks(np.random.default_rng(11), 32 * 128 * 1000 * 4, ts, 200 * 1000 * 4)
bc = np.cov(C[:, 1], C[:, 2])[0, 1]; sen = np.sqrt(np.var(C[:, 1]) * np.var(C[:, 2]) / len(C))
print(f"blocks: m={m} walkers per block, {len(C)} blocks; empirical cov {bc:.4f} (se about {sen:.4f}); exact m*cov = {m*ex:.4f}; independent would be 0")
