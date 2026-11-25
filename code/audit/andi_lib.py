"""AnDi-2 benchmark tracks -> TAMSD curves per coordinate (x = half A, y = half B), origin-block bootstrap, label-implied expected exponent."""
import numpy as np

LAG_BASE = 12      # Delta_max = floor(T / 12)
BLOCK_MIN = 4


def lag_grid(T):
    dmax = T // LAG_BASE
    g = sorted({int(round(2 ** (k / 2))) for k in range(0, 40) if round(2 ** (k / 2)) <= dmax})
    return np.array(g, dtype=int)


def block_tables(u, lags, b):
    """u: 1-D coordinate series. Returns S[nblock, K], n[nblock, K]: sum of squared displacements / counts for origins in blocks of length b."""
    T = len(u); nb = int(np.ceil((T - 1) / b))
    S = np.zeros((nb, len(lags))); n = np.zeros((nb, len(lags)))
    for k, d in enumerate(lags):
        t = np.arange(0, T - d); dx2 = (u[d:] - u[:T - d]) ** 2
        blk = t // b
        np.add.at(S[:, k], blk, dx2); np.add.at(n[:, k], blk, 1)
    return S, n


def curve_boot(S, n, B, rng):
    """y = -0.5 log TAMSD (so that d = slope of log TAMSD = alpha); bootstrap resamples origin blocks."""
    m = S.shape[0]
    y = -0.5 * np.log(S.sum(0) / n.sum(0))
    mult = rng.multinomial(m, np.full(m, 1.0 / m), size=B)
    Ss = mult @ S; ns = mult @ n
    return y, -0.5 * np.log(np.maximum(Ss, 1e-300) / np.maximum(ns, 1))


def track_inputs(xy, T, B, rng):
    """xy: (n,2) array of first T points. returns lags, (A, B, All) as (y, Ys) tuples; A = x coordinate, B = y coordinate."""
    lags = lag_grid(T)
    if len(lags) < 6: return lags, None
    b = int(max(lags[-1], BLOCK_MIN))
    Sx, nx = block_tables(xy[:T, 0], lags, b); Sy, ny = block_tables(xy[:T, 1], lags, b)
    yA, YsA = curve_boot(Sx, nx, B, rng); yB, YsB = curve_boot(Sy, ny, B, rng)
    mx = Sx.shape[0]
    # pooled: independent resampling of each coordinate's blocks, combined
    mult1 = rng.multinomial(mx, np.full(mx, 1.0 / mx), size=B); mult2 = rng.multinomial(mx, np.full(mx, 1.0 / mx), size=B)
    S = mult1 @ Sx + mult2 @ Sy; N = mult1 @ nx + mult2 @ ny
    yAll = -0.5 * np.log((Sx.sum(0) + Sy.sum(0)) / (nx.sum(0) + ny.sum(0)))
    YsAll = -0.5 * np.log(np.maximum(S, 1e-300) / np.maximum(N, 1))
    return lags, ((yA, YsA), (yB, YsB), (yAll, YsAll))


def label_expected_exponent(alpha, D, lags):
    """Label-implied local exponent d_lab(lag) = dlog E[TAMSD]/dlog lag with E[TAMSD](d) = sum_i n_i(d) D_i d^alpha_i / sum_i n_i(d),
    n_i(d) = max(len_i - d, 0): within-segment pairs only (pairs straddling a change are ignored: approximation, flagged)."""
    ch = np.nonzero((np.abs(np.diff(alpha)) > 1e-9) | (np.abs(np.diff(D)) > 1e-9))[0] + 1
    st = np.concatenate([[0], ch]); en = np.concatenate([ch, [len(alpha)]])
    segs = [(e - s, D[s], alpha[s]) for s, e in zip(st, en)]
    def logM(d):
        num = sum(max(L - d, 0.0) * Dd * d ** a for L, Dd, a in segs); den = sum(max(L - d, 0.0) for L, Dd, a in segs)
        return np.log(num / den)
    out = []
    for d in lags.astype(float):
        h = 0.02 * d
        out.append((logM(d + h) - logM(d - h)) / (np.log(d + h) - np.log(d - h)))
    return np.array(out), len(segs)
