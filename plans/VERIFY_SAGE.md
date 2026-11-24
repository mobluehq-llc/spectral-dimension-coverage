# Independent verification in SageMath 10.10 (2026-11-21, run on host-T)
Scripts (this folder): `verify_sage_1_torus.sage`, `verify_sage_2_gasket.sage`, `verify_sage_3_parity.sage`, `verify_sage_4_p20.sage` (copy in `../p20_powerlaw_log_estimators/verify_sage_p20_closed_form.sage`). Run with `~/bin_sage.sh script.sage` with `p22lib.py` beside it (check 1 and 2 import p22lib only to compare; the truth values are recomputed independently). Jobs ran with <=6 threads.

## 1. Exact torus truth — PASS
Eigenvalues lambda = 1/2 + (1/2)(1/m) sum cos(2 pi j_a/n_a) in 120-bit arithmetic, P(t) = (1/N) sum lambda^t, d(t) = -2 t sum(lambda^t ln lambda)/sum(lambda^t), on the bench's own grid `make_grid(24,34,cap)` with cap = last t having N P >= 4.5 (as in `run_p22.prep`). Compared with `Graph.exact` (numpy):

| graph | N | grid pts | max abs diff d | max rel diff P | max abs diff of grid secant d |
|---|---|---|---|---|---|
| C64xC64 | 4096 | 10 | 4.6e-14 | 1.8e-14 | 3.5e-14 |
| C128xC128 | 16384 | 14 | 1.4e-13 | 6.9e-14 | 1.2e-13 |
| C64xC16 | 1024 | 6 | 1.1e-14 | 1.6e-15 | 9.1e-15 |
| C128xC16 | 2048 | 9 | 3.7e-14 | 1.6e-14 | 3.0e-14 |
| C64xC8xC8 | 4096 | 7 | 1.0e-14 | 5.5e-15 | 8.7e-15 |
| C128xC8xC8 | 8192 | 11 | 8.7e-14 | 3.3e-14 | 7.5e-14 |
| C256 | 256 | 11 | 3.0e-14 | 8.4e-15 | 2.5e-14 |

Worst: 1.35e-13 (d), 6.9e-14 (P relative), 1.17e-13 (secant). Numpy double precision is exact to ~1e-13; dropping lambda<=1e-14 eigenvalues has no effect. Sanity: d(t_min=24) = 2.0044 (2D), 3.0217 (3D); late-time d approaches 1 for the anisotropic cases.

## 2. Sierpinski gasket — PASS (with a refinement of the quoted numbers)
- Graph rebuilt independently (recursive midpoint subdivision, different code from `L.gasket`): N = (3^(n+1)+3)/2, edges 3^(n+1)... = 9,27,...,6561, identical degree multiset, and the lazy-walk spectrum matches `p22lib` to 0.0 for n=1..7 (N up to 3282). Exact rational check at n=2: tr(S^k) = 15/2, 153/32, 921/256, 21375/8192, 17102673/8388608 for k=1,2,3,5,8 equals the eigenvalue power sums to 1e-15. (Eigenvalues are double-precision LAPACK; the exact-rational trace check certifies them at n=2.)
- d* = 2 ln3/ln5 = 1.365212.
- Bench-grid (24,34), window N P >= 6, mean exact finite-graph d: n=5 1.303 (gap 0.062 below d*), n=6 1.319 (0.046), n=7 1.335 (0.030). **Confirms "0.03-0.06 below 1.365"**, and the gap shrinks with n. Individual grid points range 1.24-1.39 and some sit ABOVE d* (n=6: 1.373; n=7: 1.390, 1.372), so "all points below" would be false; it is the window mean that is below.
- Log-periodic wiggle: local extrema of d(t) at t = 18, 40, 90, 205, 443, 1079, 2142 (n=7): successive ratios ~2.2-2.5 ~ sqrt(5) = 2.236 (half period in ln t = ln5/2). Amplitude ~ +-0.04-0.05 about the trend (cubic-detrended residual max 0.043 n=7, 0.050 n=6 in the N P >= 6 window; 0.013 at n=5 whose window is only 0.8 decade). **Confirms "oscillates +-0.05, log-periodic".** Beyond the N P >= 6 window d falls away (to 0.88-1.19 near equilibration), so wider windows make the gap far larger.

## 3. Parity trap — PASS
- Symbolic: Sage simplifies cos(2 pi (j + n/2)/n) = -cos(2 pi j/n) to True, so the spectrum of an even cycle (and any Cartesian product with an even cycle factor, being bipartite) is symmetric under lambda -> -lambda; every odd power sum, hence tr(P^t), vanishes.
- Exact rational arithmetic: tr(P^t) = 0 for all odd t <= 21 on C_8, C_6xC_4, C_4^3 (<= 11 on C_64); even t strictly positive; C_7 (odd, non-bipartite) has tr(P^7)/7 = 1/64 != 0, so bipartiteness is what matters.
- Lazy walk: tr(((I+P)/2)^t) >= N 2^-t > 0 for every t on all of these (the identity term). Numerics C_256 at 100 bits: non-lazy P(35) = -1.3e-31 (zero to precision), P(34) = 0.1358; lazy P(35) = 0.0950, and the binomial-mixture identity P_lazy(t) = sum_k C(t,k) 2^-t P_k reproduces it to 1e-10.

## 4. p20 closed form — PASS (maximum error 0.94%, never above 1%)
See `../p20_powerlaw_log_estimators/VERIFY_SAGE_P20.md`.
