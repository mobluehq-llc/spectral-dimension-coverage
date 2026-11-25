"""DFA fluctuation function with year-cluster bootstrap, for the CET audit. DFA-2 (quadratic detrending of the profile), boxes confined to 5-year blocks."""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

SCALES = np.array(sorted({int(round(7 * 2 ** (k / 2))) for k in range(0, 15)}))     # 7,10,14,20,28,40,57,79,113,160,226,320,452,640,905
ORDER = 2
BLOCK_YEARS = 5
_proj = {}


def proj(s):
    if s not in _proj:
        t = np.linspace(-1, 1, s); X = np.vander(t, ORDER + 1, increasing=True)
        _proj[s] = np.eye(s) - X @ np.linalg.pinv(X)
    return _proj[s]


def block_boxes(z, scales):
    """z: anomaly series of ONE block. returns per scale: array of box F^2 (mean squared DFA residual) and box start indices (50% overlap)."""
    Y = np.cumsum(z - z.mean())
    out = []
    for s in scales:
        if s > len(Y): out.append((np.zeros(0), np.zeros(0, int))); continue
        step = max(1, s // 2); starts = np.arange(0, len(Y) - s + 1, step)
        W = sliding_window_view(Y, s)[starts]
        R = W @ proj(s).T
        out.append(((R ** 2).mean(1), starts))
    return out


def cluster_tables(z, years, scales=SCALES):
    """z daily anomalies (len n), years int array (len n). Blocks = consecutive 5-year groups (by year index from first year).
    returns dict: block_id (per cluster year), arrays sumF2[year_idx, k], cnt[year_idx, k]; clusters are calendar years (box assigned to year of its first day)."""
    y0 = years.min(); yrs = np.unique(years); nyr = len(yrs)
    blk = (years - y0) // BLOCK_YEARS
    sumF = np.zeros((nyr, len(scales))); cnt = np.zeros((nyr, len(scales)))
    yidx = years - y0
    for b in np.unique(blk):
        idx = np.nonzero(blk == b)[0]
        if len(idx) < max(scales): continue
        res = block_boxes(z[idx], scales)
        for k, (f2, st) in enumerate(res):
            if len(f2) == 0: continue
            yy = yidx[idx[st]]
            np.add.at(sumF[:, k], yy, f2); np.add.at(cnt[:, k], yy, 1)
    return dict(sumF=sumF, cnt=cnt, blk=np.array([blk[np.nonzero(yidx == i)[0][0]] if (yidx == i).any() else -1 for i in range(nyr)]))


def split_halves(tab):
    """alternate 5-year blocks: A = even block ids, B = odd."""
    b = tab["blk"]; A = b % 2 == 0; B = (b % 2 == 1)
    sub = lambda m: dict(sumF=tab["sumF"][m], cnt=tab["cnt"][m])
    return sub(A), sub(B), dict(sumF=tab["sumF"], cnt=tab["cnt"])


def curve_and_boot(t, B, rng):
    """t: dict(sumF, cnt) rows = clusters (years). Returns y = -0.5*log F  (F = sqrt(F^2)), i.e. y = -0.25 log F^2, and B bootstrap copies (resample years)."""
    n = t["sumF"].shape[0]
    f2 = t["sumF"].sum(0) / np.maximum(t["cnt"].sum(0), 1)
    y = -0.25 * np.log(f2)
    mult = rng.multinomial(n, np.full(n, 1.0 / n), size=B)
    S = mult @ t["sumF"]; C = mult @ t["cnt"]
    Ys = -0.25 * np.log(np.maximum(S, 1e-300) / np.maximum(C, 1))
    return y, Ys


def fgn(n, H, rng):
    """Davies-Harte fractional Gaussian noise, unit variance."""
    k = np.arange(0, n + 1)
    g = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))
    r = np.concatenate([g, g[-2:0:-1]])
    lam = np.fft.fft(r).real
    lam = np.maximum(lam, 0)
    m = len(r)
    w = np.zeros(m, complex)
    w[0] = rng.standard_normal() * np.sqrt(lam[0] / m)
    w[m // 2] = rng.standard_normal() * np.sqrt(lam[m // 2] / m)
    a = rng.standard_normal(m // 2 - 1); b = rng.standard_normal(m // 2 - 1)
    w[1:m // 2] = (a + 1j * b) * np.sqrt(lam[1:m // 2] / (2 * m))
    w[m // 2 + 1:] = np.conj(w[1:m // 2][::-1])
    x = np.fft.fft(w).real
    return x[:n]
