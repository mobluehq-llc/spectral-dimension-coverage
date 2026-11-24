"""P22 library: graphs, exact spectral truth, walker + probe samplers, estimators E0/E1/E1x/E2/E2h/E3. See PLAN.md."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import scipy.sparse as sp
from scipy.stats import beta as _beta

DELTA = 0.25           # equivalence margin for the primary change rule
THR = 6.0              # N*P_hat window threshold
CMIN = 30              # min return count in window points
NBLOCK = 200           # walker blocks per dataset (all graphs together)
H_GRID = (0.5, 0.7, 1.0, 1.4, 2.0)


# ---------------------------------------------------------------- grid
def make_grid(a0, b0, tcap):
    ts, j = [], 0
    while a0 * 2 ** j <= tcap:
        for b in (a0, b0):
            if b * 2 ** j <= tcap:
                ts.append(b * 2 ** j)
        j += 1
    return np.array(sorted(ts), dtype=np.int64)


def cp(k, n, a=0.05):
    """Clopper-Pearson interval."""
    if n == 0:
        return (np.nan, np.nan)
    lo = 0.0 if k == 0 else _beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else _beta.ppf(1 - a / 2, k + 1, n - k)
    return (float(lo), float(hi))


# ---------------------------------------------------------------- graphs
class Graph:
    """kind 'torus' (dims, direct sampler) or 'step' (adjacency lists, step sampler). lazy = holding prob."""
    def __init__(self, kind, lazy=0.5, dims=None, nbr=None, deg=None, lam=None):
        self.kind, self.lazy, self.dims, self.nbr, self.deg = kind, lazy, dims, nbr, deg
        if kind == "torus":
            self.N = int(np.prod(dims))
        else:
            self.N = len(deg)
        self._lam = lam

    @property
    def lam(self):
        if self._lam is None:
            if self.kind == "torus":
                m = len(self.dims)
                g = np.meshgrid(*[np.cos(2 * np.pi * np.arange(n) / n) for n in self.dims], indexing="ij")
                self._lam = (self.lazy + (1 - self.lazy) * sum(g) / m).ravel()
            else:
                self._lam = np.linalg.eigvalsh(self.ssym().toarray())
        return self._lam

    def ssym(self):
        N = self.N
        if self.kind == "torus":
            m = len(self.dims)
            A = None
            for a, n in enumerate(self.dims):
                C = sp.diags([np.ones(n - 1), np.ones(n - 1), [1.0], [1.0]], [1, -1, n - 1, -(n - 1)], shape=(n, n)) if n > 2 else sp.csr_matrix(np.ones((2, 2)) - np.eye(2))
                f = [sp.identity(k, format="csr") for k in self.dims]
                f[a] = C
                T = f[0]
                for x in f[1:]:
                    T = sp.kron(T, x, format="csr")
                A = T if A is None else A + T
            return (self.lazy * sp.identity(N) + (1 - self.lazy) * A / (2 * m)).tocsr()
        rows = np.repeat(np.arange(N), self.nbr.shape[1])
        cols = self.nbr.ravel()
        ok = cols >= 0
        A = sp.csr_matrix((np.ones(ok.sum()), (rows[ok], cols[ok])), shape=(N, N))
        d = np.asarray(A.sum(1)).ravel()
        Dm = sp.diags(1 / np.sqrt(d))
        return (self.lazy * sp.identity(N) + (1 - self.lazy) * Dm @ A @ Dm).tocsr()

    def exact(self, ts, lam=None):
        """q(t)=mean(lam^t), d(t)=-2 t q'/q. lam<=0 eigenvalues dropped (lam^t -> 0 for the lazy walk; exact zeros contribute 0)."""
        lam = self.lam if lam is None else lam
        lp = lam[lam > 1e-14]
        ll = np.log(lp)
        ts = np.atleast_1d(np.asarray(ts, float))
        q = np.empty(len(ts)); dq = np.empty(len(ts))
        for i, t in enumerate(ts):
            a = np.exp(t * ll)
            q[i] = a.sum() / self.N
            dq[i] = (a * ll).sum() / self.N
        return q, -2 * ts * dq / q

    # ---- walkers: returns C (nblocks x K) return counts per block and walkers per block
    def sample_blocks(self, rng, W, ts, nb):
        W = int(W // nb * nb)
        K = len(ts)
        if self.kind == "torus":
            ret = self._direct(rng, W, ts)
        else:
            ret = self._steps(rng, W, ts)
        return ret.reshape(nb, W // nb, K).sum(1).astype(np.float64), W // nb

    def _direct(self, rng, W, ts):
        dims = self.dims; m = len(dims); l = self.lazy
        pos = np.zeros((W, m), dtype=np.int64)
        ret = np.zeros((W, len(ts)), dtype=np.int8)
        prev = 0
        for k, t in enumerate(ts):
            dt = int(t - prev); prev = int(t)
            rem = np.full(W, dt, dtype=np.int64)
            if l > 0:
                stay = rng.binomial(dt, l, size=W); rem = rem - stay
            for a in range(m):
                na = rem if a == m - 1 else rng.binomial(rem, 1.0 / (m - a))
                if a != m - 1:
                    rem = rem - na
                disp = 2 * rng.binomial(na, 0.5) - na
                pos[:, a] = (pos[:, a] + disp) % dims[a]
            ret[:, k] = (pos == 0).all(1)
        return ret

    def _steps(self, rng, W, ts):
        N = self.N; l = self.lazy
        start = rng.integers(0, N, size=W)
        pos = start.copy()
        tset = {int(t): k for k, t in enumerate(ts)}
        ret = np.zeros((W, len(ts)), dtype=np.int8)
        deg, nbr = self.deg, self.nbr
        for t in range(1, int(ts[-1]) + 1):
            r = rng.random(W)
            idx = np.minimum(((r - l) / (1 - l) * deg[pos]).astype(np.int64), deg[pos] - 1)
            idx = np.maximum(idx, 0)
            pos = np.where(r >= l, nbr[pos, idx], pos)
            k = tset.get(t)
            if k is not None:
                ret[:, k] = pos == start
        return ret

    # ---- probes
    def sample_probes(self, rng, P, ts):
        Z = rng.choice([-1.0, 1.0], size=(self.N, P))
        S = self.ssym()
        Y = Z.copy(); X = np.zeros((P, len(ts)))
        tset = {int(t): k for k, t in enumerate(ts)}
        for t in range(1, int(ts[-1]) + 1):
            Y = S @ Y
            k = tset.get(t)
            if k is not None:
                X[:, k] = (Z * Y).sum(0) / self.N
        return X


def gasket(n, lazy=0.5):
    tri = [(0, 0)]; s = 2 ** n
    for _ in range(n):
        s //= 2
        tri = [p for (a, b) in tri for p in ((a, b), (a + s, b), (a, b + s))]
    idx = {}; edges = set()
    def node(p):
        if p not in idx: idx[p] = len(idx)
        return idx[p]
    for (a, b) in tri:
        v = [node((a, b)), node((a + 1, b)), node((a, b + 1))]
        for i in range(3):
            for j in range(i + 1, 3):
                edges.add((min(v[i], v[j]), max(v[i], v[j])))
    N = len(idx)
    adj = [[] for _ in range(N)]
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    deg = np.array([len(a) for a in adj]); md = deg.max()
    nbr = -np.ones((N, md), dtype=np.int64)
    for i, a in enumerate(adj):
        nbr[i, :len(a)] = a
    return Graph("step", lazy, nbr=nbr, deg=deg)


def rrg(N, d, rng, lazy=0.5):
    from scipy.sparse.csgraph import connected_components
    while True:
        stubs = np.repeat(np.arange(N), d); rng.shuffle(stubs)
        a, b = stubs[0::2], stubs[1::2]
        if (a == b).any():
            continue
        e = np.sort(np.stack([a, b], 1), 1)
        if len(np.unique(e[:, 0] * N + e[:, 1])) < len(e):
            continue
        A = sp.csr_matrix((np.ones(len(e) * 2), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(N, N))
        if connected_components(A)[0] != 1:
            continue
        nbr = np.zeros((N, d), dtype=np.int64)
        Al = A.tolil().rows
        for i in range(N):
            nbr[i] = Al[i]
        return Graph("step", lazy, nbr=nbr, deg=np.full(N, d))


# ---------------------------------------------------------------- bootstrap of unit means (two-level: graphs x units)
def grand_mean(U):
    """U: list over graphs of (n_g x K) unit means -> equal-weight grand mean."""
    return np.mean([u.mean(0) for u in U], axis=0)


def boot_means(U, B, rng):
    G = len(U); K = U[0].shape[1]
    out = np.zeros((B, K))
    if G == 1:
        n = U[0].shape[0]
        mult = rng.multinomial(n, np.full(n, 1.0 / n), size=B)
        return (mult @ U[0]) / n
    for _ in range(G):
        gj = rng.integers(0, G, size=B)
        for g in range(G):
            rows = np.nonzero(gj == g)[0]
            if len(rows) == 0:
                continue
            n = U[g].shape[0]
            mult = rng.multinomial(n, np.full(n, 1.0 / n), size=len(rows))
            out[rows] += (mult @ U[g]) / n
    return out / G


# ---------------------------------------------------------------- band construction + scoring
def ll_rows(x, w, h):
    z = x[None, :] - x[:, None]
    u = np.exp(-0.5 * (z / h) ** 2) * w[None, :]
    S0 = u.sum(1); S1 = (u * z).sum(1); S2 = (u * z * z).sum(1)
    den = S0 * S2 - S1 ** 2
    slope = (S0[:, None] * u * z - S1[:, None] * u) / den[:, None]
    val = (S2[:, None] * u - S1[:, None] * u * z) / den[:, None]
    return val, slope


def score(dhat, Dstar, dtrue, tag, dpop=None):
    """bands and metrics for a vector of scale-wise estimates. dtrue may be None."""
    J = len(dhat)
    se = Dstar.std(0, ddof=1); se = np.maximum(se, 1e-9)
    c = np.quantile((np.abs(Dstar - dhat) / se).max(1), 0.95)
    r = dict(est=tag, abstain=0, J=J, hw=float(c * se.mean()), se_mean=float(se.mean()))
    if dtrue is not None:
        err = dhat - dtrue
        r.update(cov_simult=int((np.abs(err) <= c * se).all()), cov_point=float((np.abs(err) <= 1.96 * se).mean()),
                 bias=float(err.mean()), mae=float(np.abs(err).mean()))
    if dpop is not None:
        r["cov_pop_simult"] = int((np.abs(dhat - dpop) <= c * se).all())
        r["cov_pop_point"] = float((np.abs(dhat - dpop) <= 1.96 * se).mean())
        r["bias_pop"] = float((dhat - dpop).mean())
    if J < 2:
        return r
    # change test: Delta_j = d_j - d_0
    Dl = Dstar[:, 1:] - Dstar[:, [0]]; dl = dhat[1:] - dhat[0]
    sel = np.maximum(Dl.std(0, ddof=1), 1e-9)
    cl = np.quantile((np.abs(Dl - dl) / sel).max(1), 0.95)
    lo, hi = dl - cl * sel, dl + cl * sel
    r.update(decl_neg=int((hi < -DELTA).any()), decl_pos=int((lo > DELTA).any()), decl_zero=int(((hi < 0) | (lo > 0)).any()),
             dmin=float(dl.min()), dmax=float(dl.max()))
    if dtrue is not None:
        dt = dtrue[1:] - dtrue[0]
        r["delta_cov"] = int(((dt >= lo) & (dt <= hi)).all())
        r["true_dmin"] = float(dt.min())
    return r


def abst(tag, why=""):
    return dict(est=tag, abstain=1, why=why)


# ---------------------------------------------------------------- estimators
def window_len(cnt, p, N, thr, cmin=CMIN):
    ok = (cnt >= cmin) & (N * p >= thr)
    bad = np.nonzero(~ok)[0]
    return int(bad[0]) if len(bad) else len(ok)


def est_E1(Ublk, mwalk, ts, N, rng, truth, split=True, B=200, thr=THR, tag="E1"):
    """Ublk: list per graph (M x K) unit return probabilities; mwalk walkers per unit."""
    if split:
        UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    else:
        UA = UB = Ublk
    pA = grand_mean(UA)
    WA = sum(u.shape[0] for u in UA) * mwalk
    cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return abst(tag, "window")
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    best = None
    for h in H_GRID:
        val, _ = ll_rows(x, wA, h)
        hd = np.diag(val)
        cv = (wA * ((yA - val @ yA) / (1 - hd)) ** 2).sum()
        if best is None or cv < best[0]:
            best = (cv, h)
    h = best[1]
    _, sl = ll_rows(x, wA, h)
    WB = sum(u.shape[0] for u in UB) * mwalk
    floor = 0.5 / WB
    pB = np.maximum(grand_mean(UB), floor)[:J]
    dhat = -2 * sl @ np.log(pB)
    Ps = np.maximum(boot_means(UB, B, rng)[:, :J], floor)
    Ds = -2 * np.log(Ps) @ sl.T
    dtrue, dpop = truth["d_exact"](ts[:J]), None
    dpop = -2 * sl @ np.log(truth["q_exact"](ts[:J]))
    r = score(dhat, Ds, dtrue, tag, dpop)
    r.update(h=h, tlast=int(ts[J - 1]))
    return r


def est_E0(Ublk, mwalk, ts, N, rng, truth, B=200, thr=THR):
    p = grand_mean(Ublk)
    W = sum(u.shape[0] for u in Ublk) * mwalk
    cnt = p * W
    J = window_len(cnt, p, N, thr)
    if J < 6:
        return [abst("E0g", "window"), abst("E0", "window")]
    x = np.log(ts[:J].astype(float)); w = cnt[:J]
    def lsrow(idx):
        xi, wi = x[idx], w[idx]
        xb = (wi * xi).sum() / wi.sum()
        row = np.zeros(J); row[idx] = wi * (xi - xb) / (wi * (xi - xb) ** 2).sum()
        return row, xi.mean()
    half = (J + 1) // 2
    rg, xg = lsrow(np.arange(J)); r1, x1 = lsrow(np.arange(half)); r2, x2 = lsrow(np.arange(half, J))
    R = np.stack([r1, r2]); Rg = rg[None, :]
    floor = 0.5 / W
    y = np.log(np.maximum(p[:J], floor))
    Ps = np.maximum(boot_means(Ublk, B, rng)[:, :J], floor); Ys = np.log(Ps)
    out = []
    tref = np.exp([xg])
    d_g = -2 * (Rg @ y); D_g = -2 * Ys @ Rg.T
    qe = np.log(truth["q_exact"](ts[:J]))
    rg_ = score(d_g, D_g[:, :], truth["d_exact"](tref), "E0g", -2 * (Rg @ qe))
    out.append(rg_)
    tref2 = np.exp([x1, x2])
    out.append(score(-2 * (R @ y), -2 * Ys @ R.T, truth["d_exact"](tref2), "E0", -2 * (R @ qe)))
    return out


def est_E2(Xs, ts, N, rng, truth, B=400, thr=THR, tag="E2", step=2):
    """Xs: list per graph of (P x K) probe values. returns [E2, E2h]."""
    q = grand_mean(Xs)
    Jq = window_len(np.full(len(q), 1e9), np.maximum(q, 0) , N, thr)
    Jd = Jq - step
    if Jd < 6:
        return [abst("E2", "window"), abst("E2h", "window")]
    Ptot = sum(x.shape[0] for x in Xs)
    lr = np.log(ts[step:Jd + step] / ts[:Jd].astype(float))
    qe = np.maximum(q[:Jq], 1e-12)
    dhat = -2 * np.log(qe[step:Jd + step] / qe[:Jd]) / lr
    Qs = np.maximum(boot_means(Xs, B, rng)[:, :Jq], 1e-12)
    Ds = -2 * np.log(Qs[:, step:Jd + step] / Qs[:, :Jd]) / lr
    tref = np.sqrt(ts[:Jd] * ts[step:Jd + step].astype(float))
    dtrue = truth["d_exact"](tref)
    qx = truth["q_exact"](ts[:Jq])
    dpop = -2 * np.log(qx[step:Jd + step] / qx[:Jd]) / lr
    r = score(dhat, Ds, dtrue, "E2", dpop)
    # Hoeffding: per-probe terms in [0,1]; union bound over Jq grid points at total level 0.05
    eps = np.sqrt(np.log(2 * Jq / 0.05) / (2 * Ptot))
    ql, qh = np.maximum(qe - eps, 1e-12), qe + eps
    dlo = -2 * np.log(qh[step:Jd + step] / ql[:Jd]) / lr
    dhi = -2 * np.log(ql[step:Jd + step] / qh[:Jd]) / lr
    rh = dict(est="E2h", abstain=0, J=Jd, hw=float(np.mean(dhi - dlo) / 2), se_mean=float("nan"))
    rh["cov_simult"] = int(((dtrue >= dlo) & (dtrue <= dhi)).all())
    rh["cov_point"] = float(((dtrue >= dlo) & (dtrue <= dhi)).mean())
    rh["bias"] = float(((dlo + dhi) / 2 - dtrue).mean()); rh["mae"] = float(np.abs((dlo + dhi) / 2 - dtrue).mean())
    Dlo = dlo[1:] - dhi[0]; Dhi = dhi[1:] - dlo[0]
    rh.update(decl_neg=int((Dhi < -DELTA).any()), decl_pos=int((Dlo > DELTA).any()), decl_zero=int(((Dhi < 0) | (Dlo > 0)).any()))
    dt = dtrue[1:] - dtrue[0]
    rh["delta_cov"] = int(((dt >= Dlo) & (dt <= Dhi)).all())
    return [r, rh]


def _collapse_fit(P, cnt, Ns, ts, tsplit, region, dgrid):
    """P: list of arrays p_i(t_k) (window prefix), cnt: counts. returns d minimising pairwise collapse residual in region."""
    nS = len(P)
    ly = [np.log(Ns[i] * P[i]) for i in range(nS)]
    lt = [np.log(ts[i].astype(float)) for i in range(nS)]
    msk = [(ts[i] <= tsplit) if region == 0 else (ts[i] > tsplit) for i in range(nS)]
    def obj(d):
        tot, npt = 0.0, 0.0
        u = [lt[i] - (2.0 / d) * np.log(Ns[i]) for i in range(nS)]
        for i in range(nS):
            for j in range(nS):
                if i == j:
                    continue
                inr = msk[i] & (u[i] >= u[j][0]) & (u[i] <= u[j][-1])
                if inr.sum() == 0:
                    continue
                yy = np.interp(u[i][inr], u[j], ly[j])
                w = cnt[i][inr]
                tot += (w * (ly[i][inr] - yy) ** 2).sum(); npt += w.sum()
        return tot / npt if npt > 0 else np.inf, npt
    vals = np.array([obj(d)[0] for d in dgrid])
    if not np.isfinite(vals).any():
        return np.nan
    d0 = dgrid[np.argmin(vals)]
    fine = np.arange(d0 - 0.1, d0 + 0.1001, 0.02)
    fine = fine[fine > 0.2]
    vf = np.array([obj(d)[0] for d in fine])
    return float(fine[np.argmin(vf)])


def est_E3(fam, rng, truth_val, B=80, thr=1.05, t0=24.0, dmax=3.5):
    """fam: list of (U blocks (M x K) unit probs, mwalk, ts, N) for sizes largest..smallest. truth_val: plateau or None."""
    P, cnt, TS, Ns, Wt = [], [], [], [], []
    for (U, mw, ts, N) in fam:
        p = U.mean(0); W = U.shape[0] * mw; c = p * W
        J = window_len(c, p, N, thr)
        if J < 5:
            return abst("E3", "window")
        P.append(p[:J]); cnt.append(c[:J]); TS.append(ts[:J]); Ns.append(N); Wt.append(W)
    tsplit = np.sqrt(t0 * TS[0][-1])
    dgrid = np.arange(0.5, dmax + 0.0001, 0.1)
    def both(Ps):
        return [_collapse_fit(Ps, cnt, Ns, TS, tsplit, r, dgrid) for r in (0, 1)]
    est = both(P)
    if not np.all(np.isfinite(est)):
        return abst("E3", "nofit")
    boot = []
    for _ in range(B):
        Pb = []
        for (U, mw, ts, N), p0, Wi in zip(fam, P, Wt):
            n = U.shape[0]
            m = rng.multinomial(n, np.full(n, 1.0 / n))
            pb = (m @ U) / n
            Pb.append(np.maximum(pb[:len(p0)], 0.5 / Wi))
        b = both(Pb)
        if np.all(np.isfinite(b)):
            boot.append(b)
    boot = np.array(boot)
    if len(boot) < 0.8 * B:
        return abst("E3", "bootfail")
    dhat = np.array(est); se = np.maximum(boot.std(0, ddof=1), 1e-9)
    c = np.quantile((np.abs(boot - dhat) / se).max(1), 0.95)
    r = dict(est="E3", abstain=0, J=2, hw=float(c * se.mean()), se_mean=float(se.mean()), tsplit=float(tsplit))
    if truth_val is not None:
        err = dhat - truth_val
        r.update(cov_simult=int((np.abs(err) <= c * se).all()), cov_point=float((np.abs(err) <= 1.96 * se).mean()),
                 bias=float(err.mean()), mae=float(np.abs(err).mean()))
    dl = boot[:, 1] - boot[:, 0]; sD = max(dl.std(ddof=1), 1e-9); D = dhat[1] - dhat[0]
    lo, hi = D - 1.96 * sD, D + 1.96 * sD
    r.update(decl_neg=int(hi < -DELTA), decl_pos=int(lo > DELTA), decl_zero=int(hi < 0 or lo > 0), dmin=float(D), dmax=float(D))
    if truth_val is not None:
        r["delta_cov"] = int(lo <= 0 <= hi)
    r["d_early"], r["d_late"] = float(dhat[0]), float(dhat[1])
    return r
