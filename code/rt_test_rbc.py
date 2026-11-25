"""Unit tests for the RBC local-cubic slope (rt_rbc.py uses p22amend.poly_rows with deg=3)."""
import numpy as np, p22amend as A
rng = np.random.default_rng(1)
x = np.log(np.array([24, 34, 48, 68, 96, 136, 192, 272, 384, 544.0])); w = rng.uniform(5, 500, len(x))
# 1. exact on cubics under arbitrary weights and bandwidths
for h in (0.3, 0.7, 2.0, 4.0):
    c = rng.normal(size=4); y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3; dy = c[1] + 2*c[2]*x + 3*c[3]*x**2
    _, sl = A.poly_rows(x, w, h, 3)
    assert np.max(np.abs(sl @ y - dy)) < 1e-8, (h, np.max(np.abs(sl @ y - dy)))
# 2. quadratic fit is NOT exact on a cubic (bias exists), cubic fit is
_, sl2 = A.poly_rows(x, w, 0.7, 2); c = np.array([0, 1, 0, 1.0]); y = x**3 + x; dy = 3*x**2 + 1
assert np.max(np.abs(sl2 @ y - dy)) > 1e-3
# 3. convergence rates on a fine grid (interior points): quadratic-fit slope bias ~ h^2, cubic-fit slope bias ~ h^4
xs = np.linspace(2.0, 9.0, 1400); ws = np.ones_like(xs)
def err(h, deg):
    _, s = A.poly_rows(xs, ws, h, deg); m = slice(400, 1000)
    return np.abs(s @ np.sin(xs) - np.cos(xs))[m].max()
r2 = err(0.4, 2) / err(0.2, 2); r3 = err(0.4, 3) / err(0.2, 3)
assert 3.0 < r2 < 5.0, r2          # ~4
assert 12.0 < r3 < 20.0, r3        # ~16
# at the benchmark scale (h=0.7, spacing 0.35 grid) the cubic is better than the quadratic by a modest factor (recorded, not asserted at 5x)
xs = np.linspace(3.0, 7.0, 40); ws = np.ones_like(xs)
_, s2 = A.poly_rows(xs, ws, 0.7, 2); _, s3 = A.poly_rows(xs, ws, 0.7, 3)
e2 = np.abs(s2 @ np.sin(xs) - np.cos(xs))[5:-5].max(); e3 = np.abs(s3 @ np.sin(xs) - np.cos(xs))[5:-5].max()
assert e3 < e2
print("rt_test_rbc OK rates quad %.2f cubic %.2f; h=0.7 errors quad %.3e cubic %.3e (ratio %.1f)" % (r2, r3, e2, e3, e2 / e3))
