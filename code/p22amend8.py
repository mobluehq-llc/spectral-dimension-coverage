"""P22 Amendment 8: E7 (moment-cone LP band), E8 (secant bracket), A1' alarm. See PLAN.md AMENDMENT 8.
Both arms use the full walker data of the dataset (no split), the E0-style window rule (count >= 30 and N*p >= thr, initial run), and the SAME
block-bootstrap draws (B = BOOT, rng supplied by the caller). Theorem used: for the lazy walk q(t)=sum w_k exp(-t s_k), w>=0, s>=0, w(s=0)>=1/N."""
import time
import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta
import p22lib as L

BOOT = 1000            # bootstrap draws (quantiles at 0.9875 / 0.95 need more than the B=200 of earlier arms; not specified in PLAN, stated)
A_BOX = 0.025          # total miss level of the q-box
GATE = 0.005
NATOMS = (400, 799, 1597)   # nested log grids: 400 -> 2*399+1 -> 2*798+1 (each contains the previous)
DELTA = L.DELTA


# ------------------------------------------------------------------ box
def make_box(U, m, rng, B=BOOT, a_box=A_BOX, qtrue=None):
    """U: (nb x K) block return probabilities, m walkers per block. Returns qhat, lo, hi, Ps (bootstrap draws), info.
    Points with total count >= 30: bootstrap max-t at level 1-a_box/2; others: exact Clopper-Pearson with Bonferroni over them at total a_box/2."""
    nb, K = U.shape; Wt = nb * m
    qhat = U.mean(0); cnt = qhat * Wt
    Ps = L.boot_means([U], B, rng)
    big = cnt >= L.CMIN
    lo = np.zeros(K); hi = np.zeros(K)
    if big.any():
        se = np.maximum(Ps[:, big].std(0, ddof=1), 1e-300)
        T = (np.abs(Ps[:, big] - qhat[big]) / se).max(1)
        c = float(np.quantile(T, 1 - a_box / 4))
        lo[big] = np.maximum(qhat[big] - c * se, 0.0); hi[big] = qhat[big] + c * se
    else:
        c = np.nan
    ncp = int((~big).sum())
    if ncp:
        al = (a_box / 2) / ncp
        k = np.rint(cnt[~big]).astype(float)
        lo[~big] = np.where(k > 0, beta.ppf(al / 2, np.maximum(k, 1e-9), Wt - k + 1), 0.0)
        hi[~big] = beta.ppf(1 - al / 2, k + 1, Wt - k)
    return qhat, lo, hi, Ps, dict(c_box=c, n_boot=int(big.sum()), n_cp=ncp)


# ------------------------------------------------------------------ LP
def atoms_grid(tmax, n):
    s = np.exp(np.linspace(np.log(0.02 / tmax), 0.0, n))
    return np.concatenate([[0.0], s, [np.inf]])        # lambda=1 atom, grid, lambda=0 atom


def lp_setup(ts, lo, hi, s_atoms, N):
    """returns scaled matrices reused for every target time. rows scaled by ref_i = (lo_i+hi_i)/2."""
    ts = np.asarray(ts, float)
    with np.errstate(invalid="ignore"):
        Aall = np.exp(-np.outer(ts, np.where(np.isinf(s_atoms), 0.0, s_atoms)))
    Aall[:, np.isinf(s_atoms)] = 0.0
    ref = np.maximum(0.5 * (lo + hi), 1e-300)
    As = Aall / ref[:, None]
    n = len(s_atoms); K = len(ts)
    # variables x = [y' (n), tau']
    rows = []
    rows.append(np.hstack([-As, (lo / ref)[:, None]]))      # lo*tau - A y <= 0
    rows.append(np.hstack([As, -(hi / ref)[:, None]]))      # A y - hi*tau <= 0
    e0 = np.zeros((1, n + 1)); e0[0, 0] = -1.0; e0[0, n] = 1.0 / N          # -y0 + tau/N <= 0
    rows.append(e0)
    Aub = np.vstack(rows); bub = np.zeros(Aub.shape[0])
    sfin = np.where(np.isinf(s_atoms), 0.0, s_atoms)
    return dict(As=As, Aub=Aub, bub=bub, sfin=sfin, n=n, ts=ts)


def lp_one(S, j):
    """min and max of d(t_j) under the box. returns (dmin, dmax, status) with status 0 ok, 2 infeasible, else other."""
    n = S["n"]
    Aeq = np.zeros((2, n + 1)); beq = np.array([1.0, 0.0])
    Aeq[0, :n] = S["As"][j]; Aeq[1, :n] = 1.0; Aeq[1, n] = -1.0
    cvec = np.zeros(n + 1); cvec[:n] = 2.0 * S["ts"][j] * S["sfin"] * S["As"][j]
    out = []
    st = 0
    for sgn in (1.0, -1.0):
        r = linprog(sgn * cvec, A_ub=S["Aub"], b_ub=S["bub"], A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
        if r.status != 0:
            return np.nan, np.nan, int(r.status)
        out.append(sgn * r.fun)
    return out[0], out[1], 0


def lp_band(ts, lo, hi, N, idx, tmax, n):
    S = lp_setup(ts, lo, hi, atoms_grid(tmax, n), N)
    res = [lp_one(S, j) for j in idx]
    return np.array([r[0] for r in res]), np.array([r[1] for r in res]), max(r[2] for r in res)


def e7_gated(ts, lo, hi, N, idx):
    """band with the atom-refinement validity gate. Returns dmin, dmax, status, info."""
    tmax = float(ts[-1]); info = dict(natoms=NATOMS[0], gate1=np.nan, gate2=np.nan, gate_fail=0)
    b0 = lp_band(ts, lo, hi, N, idx, tmax, NATOMS[0])
    b1 = lp_band(ts, lo, hi, N, idx, tmax, NATOMS[1])
    def diff(a, b):
        if a[2] != b[2]: return np.inf
        if a[2] != 0: return 0.0 if a[2] == b[2] else np.inf
        return float(max(np.abs(a[0] - b[0]).max(), np.abs(a[1] - b[1]).max()))
    d1 = diff(b0, b1); info["gate1"] = d1; fin = b1; info["natoms"] = NATOMS[1]
    if d1 > GATE or b1[2] != 0:                 # unstable or infeasible: refine once more
        b2 = lp_band(ts, lo, hi, N, idx, tmax, NATOMS[2]); d2 = diff(b1, b2)
        info["gate2"] = d2; fin = b2; info["natoms"] = NATOMS[2]
        info["gate_fail"] = int(d2 > GATE)
    return fin[0], fin[1], fin[2], info


# ------------------------------------------------------------------ E8
def e8_bands(ts, qhat, Ps, cnt, idx, floor, level=0.95):
    """secant bracket + max-t widening. idx: window indices to evaluate. Partners must have total count >= CMIN. Returns lo, hi (inf where open), info."""
    K = len(ts); lq = np.log(np.maximum(qhat, floor)); lP = np.log(np.maximum(Ps, floor)); t = ts.astype(float)
    stat = []   # (j, kind, partner)
    for j in idx:
        for p in (j + 1, j + 2):
            if p < K and cnt[p] >= L.CMIN: stat.append((j, 0, p))
        for p in (j - 1, j - 2):
            if p >= 0 and cnt[p] >= L.CMIN: stat.append((j, 1, p))
    def val(l, j, kind, p):
        return 2 * t[j] * (l[..., j] - l[..., p]) / (t[p] - t[j]) if kind == 0 else 2 * t[j] * (l[..., p] - l[..., j]) / (t[j] - t[p])
    S = np.array([val(lq, j, k, p) for j, k, p in stat])
    SB = np.stack([val(lP, j, k, p) for j, k, p in stat], 1)
    se = np.maximum(SB.std(0, ddof=1), 1e-12)
    c = float(np.quantile((np.abs(SB - S) / se).max(1), level))
    lo = {j: 0.0 for j in idx}; hi = {j: np.inf for j in idx}
    for (j, k, p), s_, e_ in zip(stat, S, se):
        if k == 0: lo[j] = max(lo[j], s_ - c * e_)
        else: hi[j] = min(hi[j], s_ + c * e_)
    return np.array([lo[j] for j in idx]), np.array([hi[j] for j in idx]), dict(c8=c, n_stat=len(stat))


# ------------------------------------------------------------------ scoring (common to E7/E8)
def score_band(tag, lo, hi, dtrue, extra=None):
    """lo, hi over window points (may contain +-inf / nan = infeasible). Simultaneous coverage = all points covered. Change test via interval arithmetic
    on Delta_j = d_j - d_0: [lo_j - hi_0, hi_j - lo_0]."""
    J = len(lo); r = dict(est=tag, abstain=0, J=J)
    if np.isnan(lo).any():
        r.update(infeasible=1, cov_simult=0, cov_point=0.0, hw=np.nan, decl_neg=0, decl_pos=0, decl_zero=0, delta_cov=0)
        if extra: r.update(extra)
        return r
    cov = (dtrue >= lo - 1e-12) & (dtrue <= hi + 1e-12)
    fin = np.isfinite(hi) & np.isfinite(lo)
    r.update(infeasible=0, cov_simult=int(cov.all()), cov_point=float(cov.mean()), hw=float(((hi - lo)[fin] / 2).mean()) if fin.any() else np.nan,
             n_open=int((~fin).sum()), width_mean=float((hi - lo)[fin].mean()) if fin.any() else np.nan)
    if J >= 2:
        dl_lo = lo[1:] - hi[0]; dl_hi = hi[1:] - lo[0]
        r.update(decl_neg=int((dl_hi < -DELTA).any()), decl_pos=int((dl_lo > DELTA).any()), decl_zero=int(((dl_hi < 0) | (dl_lo > 0)).any()))
        dt = dtrue[1:] - dtrue[0]; r["delta_cov"] = int(((dt >= dl_lo - 1e-12) & (dt <= dl_hi + 1e-12)).all())
    if extra: r.update(extra)
    return r


def run_e7_e8(U, m, ts, N, rng, truth, thr, J6=0, e6lo=None, e6hi=None):
    """U: (nb x K) block return probabilities. Returns list of rows [E7, E8] (+ alarm fields on E7 row when e6 band is given)."""
    nb, K = U.shape; Wt = nb * m
    p = U.mean(0); cnt = p * Wt
    Jw = L.window_len(cnt, p, N, thr)
    if Jw < 2:
        return [dict(est="E7", abstain=0, nowin=1, J=Jw), dict(est="E8", abstain=0, nowin=1, J=Jw)]
    Jev = max(Jw, J6)
    idx = list(range(Jev))
    qhat, lo, hi, Ps, binfo = make_box(U, m, rng)
    t0 = time.time()
    dmin, dmax, st, ginfo = e7_gated(ts, lo, hi, N, idx)
    tlp = time.time() - t0
    dtrue = truth["d_exact"](ts[:Jw])
    ex = dict(lp_seconds=tlp, nowin=0, tlast=int(ts[Jw - 1]), lp_status=st, **ginfo, **binfo)
    if st != 0:
        lo7 = np.full(Jw, np.nan); hi7 = lo7.copy()
    else:
        lo7, hi7 = dmin[:Jw], dmax[:Jw]
    r7 = score_band("E7", lo7, hi7, dtrue, ex)
    r7["e7lo"] = ";".join("%.5f" % v for v in (dmin if st == 0 else [])); r7["e7hi"] = ";".join("%.5f" % v for v in (dmax if st == 0 else []))
    if J6 and e6lo is not None:
        if st == 0:
            n = min(J6, len(dmin)); gap = np.maximum(e6lo[:n] - dmax[:n], dmin[:n] - e6hi[:n])
            r7["a1p_alarm"] = int((gap > 0).any()); r7["a1p_maxgap"] = float(gap.max()); r7["a1p_infeasible"] = 0
        else:
            r7["a1p_alarm"] = 0; r7["a1p_maxgap"] = np.nan; r7["a1p_infeasible"] = 1
    lo8, hi8, i8 = e8_bands(ts, qhat, Ps, cnt, list(range(Jw)), 0.5 / Wt)
    r8 = score_band("E8", lo8, hi8, dtrue, dict(nowin=0, tlast=int(ts[Jw - 1]), **i8))
    return [r7, r8]
