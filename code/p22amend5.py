"""P22 Amendment 5: E6 lack-of-fit diagnostic, Millan ML power-law fit (M, Mb), optimised finite-size collapse (CL), brush graph. See PLAN.md AMENDMENT 5."""
import numpy as np
from scipy.stats import chi2
from scipy.optimize import minimize_scalar
import p22lib as L
import p22shape as S
from p22lib import grand_mean, boot_means, abst, window_len, DELTA


# ------------------------------------------------------------------ graphs
def brush(Lb, Lf, lazy=0.5):
    """Dirac brush (Grimaldi et al. 2510.20520): periodic Lb x Lb base; a periodic ring fibre of Lf sites through every base node (fibre site 0)."""
    nb = Lb * Lb; N = nb * Lf
    adj = [[] for _ in range(N)]
    for b in range(nb):
        i, j = divmod(b, Lb)
        for s in range(Lf):
            adj[b * Lf + s] += [b * Lf + (s + 1) % Lf, b * Lf + (s - 1) % Lf]
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            adj[b * Lf].append((((i + di) % Lb) * Lb + (j + dj) % Lb) * Lf)
    deg = np.array([len(a) for a in adj]); nbr = -np.ones((N, deg.max()), dtype=np.int64)
    for i, a in enumerate(adj): nbr[i, :len(a)] = a
    return L.Graph("step", lazy, nbr=nbr, deg=deg)


# ------------------------------------------------------------------ E6 with lack-of-fit
def e6_lof(x, yA, SA, yB, Ys, truth_d, tag="E6"):
    J = len(x)
    SAs = (1 - S.SHRINK) * SA + S.SHRINK * np.diag(np.diag(SA)) + 1e-10 * np.trace(SA) / J * np.eye(J)
    SAi = np.linalg.inv(SAs)
    SB = np.cov(Ys, rowvar=False); SBs = SB + 1e-10 * np.trace(SB) / J * np.eye(J)
    SBi = np.linalg.inv(SBs)
    ms = np.linspace(x[0] - 1.0, x[-1] + 1.0, S.M_GRID_N)
    base = np.stack([np.ones(J), x], 1)
    fits = []
    for m in ms:
        for w in S.W_GRID:
            phi = 0.5 * w * S.softplus((x - m) / w)
            X = np.column_stack([base, phi])
            Kx = np.linalg.solve(X.T @ SAi @ X, X.T @ SAi)
            res = yA - X @ (Kx @ yA)
            fits.append((m, w, float(res @ SAi @ res)))
    rss = np.array([f[2] for f in fits]); keep = rss - rss.min() <= S.CHI2_2_99
    Lo = np.full(J - 1, np.inf); Hi = np.full(J - 1, -np.inf)
    cLo = np.full(J, np.inf); cHi = np.full(J, -np.inf)
    cr = np.sqrt(S.CHI2_2_96); Qmin = np.inf
    for (m, w, _), k in zip(fits, keep):
        if not k: continue
        phi = 0.5 * w * S.softplus((x - m) / w)
        X = np.column_stack([base, phi])
        Kf = np.linalg.solve(X.T @ SBi @ X, X.T @ SBi)
        rB = yB - X @ (Kf @ yB); Qmin = min(Qmin, float(rB @ SBi @ rB))
        K = Kf[1:]
        th = K @ yB; Sg = K @ SB @ K.T
        s = S.sigma((x - m) / w)
        V = np.stack([np.full(J, -2.0), -s], 1)
        cen = V @ th; sdv = np.sqrt(np.maximum(np.einsum("ji,ik,jk->j", V, Sg, V), 0))
        cLo = np.minimum(cLo, cen - cr * sdv); cHi = np.maximum(cHi, cen + cr * sdv)
        V2 = np.stack([np.zeros(J - 1), -(s[1:] - s[0])], 1)
        cen2 = V2 @ th; sd2 = np.sqrt(np.maximum(np.einsum("ji,ik,jk->j", V2, Sg, V2), 0))
        Lo = np.minimum(Lo, cen2 - cr * sd2); Hi = np.maximum(Hi, cen2 + cr * sd2)
    Bn = Ys.shape[0]
    thr = chi2.ppf(0.95, J - 3) * (Bn - 1) / (Bn - J)
    r = dict(est=tag, abstain=0, J=J, hw=float(np.mean(cHi - cLo) / 2), n_set=int(keep.sum()), n_grid=len(fits),
             decl_neg=int((Hi < -DELTA).any()), decl_pos=int((Lo > DELTA).any()), decl_zero=int(((Hi < 0) | (Lo > 0)).any()),
             dmin=float(Lo.min()), dmax=float(Hi.max()), lof_q=Qmin, lof_thr=float(thr), lof_alarm=int(Qmin > thr), lof_ratio=float(Qmin / thr))
    if truth_d is not None:
        r["cov_simult"] = int(((truth_d >= cLo) & (truth_d <= cHi)).all()); r["cov_point"] = float(((truth_d >= cLo) & (truth_d <= cHi)).mean())
        dt = truth_d[1:] - truth_d[0]
        r["delta_cov"] = int(((dt >= Lo) & (dt <= Hi)).all()); r["true_dmin"] = float(dt.min())
    return r


def est_E6(Ublk, mwalk, ts, N, rngB, rngA, truth, thr=L.THR, B=200):
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA); WA = sum(u.shape[0] for u in UA) * mwalk; cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return abst("E6", "window")
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J])
    WB = sum(u.shape[0] for u in UB) * mwalk; floor = 0.5 / WB
    yB = np.log(np.maximum(grand_mean(UB), floor)[:J])
    Ys = np.log(np.maximum(boot_means(UB, B, rngB)[:, :J], floor))
    PsA = np.maximum(boot_means(UA, B, rngA)[:, :J], 0.5 / WA)
    SA = np.cov(np.log(PsA), rowvar=False)
    r = e6_lof(x, yA, SA, yB, Ys, truth["d_exact"](ts[:J])); r["tlast"] = int(ts[J - 1])
    return r


# ------------------------------------------------------------------ M: Millan-style ML power-law fit (Poisson counts), CSN-style window choice by KS
MINPTS = 4


def _fit_pl(x, c):
    """Poisson ML of c_k ~ A t_k^-beta. returns beta, v (Var_w(x)), D (KS distance of normalised cumulatives), csum."""
    xb = (c * x).sum() / c.sum(); x0 = x - x.mean(); beta = 1.0
    for _ in range(40):
        w = np.exp(-beta * x0); w /= w.sum()
        m = (w * x0).sum(); v = (w * x0 * x0).sum() - m * m
        g = m - (xb - x.mean())
        step = g / max(v, 1e-12); beta = float(np.clip(beta + step, -2.0, 14.0))
        if abs(step) < 1e-10: break
    w = np.exp(-beta * x0); w /= w.sum(); m = (w * x0).sum(); v = (w * x0 * x0).sum() - m * m
    D = float(np.abs(np.cumsum(c) / c.sum() - np.cumsum(w)).max())
    return beta, v, D, float(c.sum())


def _best_window(x, c):
    best = None
    n = len(x)
    for i0 in range(n - MINPTS + 1):
        for i1 in range(i0 + MINPTS, n + 1):
            if c[i0:i1].sum() <= 0: continue
            beta, v, D, cs = _fit_pl(x[i0:i1], c[i0:i1])
            if best is None or D < best[0]:
                best = (D, 2 * beta, 2.0 / np.sqrt(max(cs * v, 1e-12)), float(x[i0:i1].mean()))
    return best  # D, d, sd_fisher, centre


def est_M(Ublk, mwalk, ts, N, rng, truth, thr=L.THR, B=80, plateau=None):
    p = grand_mean(Ublk); W = sum(u.shape[0] for u in Ublk) * mwalk; cnt = p * W
    J = window_len(cnt, p, N, thr)
    if J < 6: return [abst("M", "window"), abst("Mb", "window")]
    t = ts[:J].astype(float); x = np.log(t)
    tsplit = np.sqrt(24.0 * t[-1])
    ie = np.nonzero(t <= tsplit)[0]; il = np.nonzero(t > tsplit)[0]
    if len(ie) < MINPTS or len(il) < MINPTS: return [abst("M", "halves"), abst("Mb", "halves")]
    c0 = cnt[:J]
    res = [_best_window(x[ie], c0[ie]), _best_window(x[il], c0[il])]
    if any(r is None for r in res): return [abst("M", "nofit"), abst("Mb", "nofit")]
    dh = np.array([r[1] for r in res]); sd = np.array([r[2] for r in res]); cen = np.array([r[3] for r in res])
    dtr = truth["d_exact"](np.exp(cen))
    def base(tag):
        return dict(est=tag, abstain=0, J=J, tlast=int(t[-1]), d_early=float(dh[0]), d_late=float(dh[1]), win_e=float(res[0][3]), win_l=float(res[1][3]), ks_e=res[0][0], ks_l=res[1][0])
    # --- M: analytic (Fisher) error
    lo, hi = dh - 1.96 * sd, dh + 1.96 * sd
    r = base("M"); err = dh - dtr
    r.update(hw=float(1.96 * sd.mean()), se_mean=float(sd.mean()), cov_simult=int((np.abs(err) <= 1.96 * sd).all()), cov_point=float((np.abs(err) <= 1.96 * sd).mean()),
             bias=float(err.mean()), mae=float(np.abs(err).mean()))
    D = dh[1] - dh[0]; sD = float(np.sqrt((sd ** 2).sum()))
    r.update(decl_neg=int(D + 1.96 * sD < -DELTA), decl_pos=int(D - 1.96 * sD > DELTA), decl_zero=int((D + 1.96 * sD < 0) or (D - 1.96 * sD > 0)), dmin=float(D), dmax=float(D))
    dt = dtr[1] - dtr[0]; r["delta_cov"] = int(D - 1.96 * sD <= dt <= D + 1.96 * sD); r["true_dmin"] = float(dt)
    if plateau is not None: r["cov_simult_plateau"] = int((np.abs(dh - plateau) <= 1.96 * sd).all())
    out = [r]
    # --- Mb: block bootstrap, window re-optimised on every resample
    Ps = boot_means(Ublk, B, rng)[:, :J]; floor = 0.5 / W
    bs = []
    for b in range(B):
        cb = np.maximum(Ps[b], floor) * W
        rb = [_best_window(x[ie], cb[ie]), _best_window(x[il], cb[il])]
        if any(q is None for q in rb): continue
        bs.append([rb[0][1], rb[1][1]])
    bs = np.array(bs)
    if len(bs) < 0.8 * B: return [r, abst("Mb", "bootfail")]
    se = np.maximum(bs.std(0, ddof=1), 1e-9); c = np.quantile((np.abs(bs - dh) / se).max(1), 0.95)
    rb_ = base("Mb"); rb_.update(hw=float(c * se.mean()), se_mean=float(se.mean()), cov_simult=int((np.abs(err) <= c * se).all()),
                                  cov_point=float((np.abs(err) <= 1.96 * se).mean()), bias=float(err.mean()), mae=float(np.abs(err).mean()))
    sDb = max(float((bs[:, 1] - bs[:, 0]).std(ddof=1)), 1e-9)
    rb_.update(decl_neg=int(D + 1.96 * sDb < -DELTA), decl_pos=int(D - 1.96 * sDb > DELTA), decl_zero=int((D + 1.96 * sDb < 0) or (D - 1.96 * sDb > 0)), dmin=float(D), dmax=float(D))
    rb_["delta_cov"] = int(D - 1.96 * sDb <= dt <= D + 1.96 * sDb); rb_["true_dmin"] = float(dt)
    if plateau is not None: rb_["cov_simult_plateau"] = int((np.abs(dh - plateau) <= c * se).all())
    out.append(rb_)
    return out


# ------------------------------------------------------------------ CL: optimised collapse
DGRID = np.linspace(0.4, 8.0, 40)


def _hh_cost(d, lt, ly, cnt, lN, msk):
    nS = len(lt)
    u = [lt[i] - (2.0 / d) * lN[i] for i in range(nS)]
    tot, n = 0.0, 0
    for i in range(nS):
        ui, yi, ci = u[i][msk[i]], ly[i][msk[i]], cnt[i][msk[i]]
        if len(ui) == 0: continue
        for j in range(nS):
            if i == j: continue
            uj, yj, cj = u[j][msk[j]], ly[j][msk[j]], cnt[j][msk[j]]
            if len(uj) < 2: continue
            inr = (ui >= uj[0]) & (ui <= uj[-1])
            if not inr.any(): continue
            yy = np.interp(ui[inr], uj, yj); vj = np.interp(ui[inr], uj, 1.0 / cj)
            tot += (((yi[inr] - yy) ** 2) / (1.0 / ci[inr] + vj)).sum(); n += int(inr.sum())
    return tot / n if n >= 3 else np.inf


def _collapse_opt(P, cnt, Ns, TS, tsplit, region):
    nS = len(P)
    ly = [np.log(Ns[i] * P[i]) for i in range(nS)]; lt = [np.log(TS[i].astype(float)) for i in range(nS)]
    lN = [np.log(float(n)) for n in Ns]
    msk = [(TS[i] <= tsplit) if region == 0 else (TS[i] > tsplit) for i in range(nS)]
    f = lambda d: _hh_cost(d, lt, ly, cnt, lN, msk)
    vals = np.array([f(d) for d in DGRID])
    if not np.isfinite(vals).any(): return np.nan
    k = int(np.argmin(vals))
    lo, hi = DGRID[max(k - 1, 0)], DGRID[min(k + 1, len(DGRID) - 1)]
    g = lambda d: (lambda v: v if np.isfinite(v) else 1e12)(f(d))
    r = minimize_scalar(g, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-3))
    return float(r.x) if g(r.x) <= vals[k] else float(DGRID[k])


def est_CL(fam, rng, truth, plateau=None, B=80, thr=1.05, t0=24.0):
    P, cnt, TS, Ns, Wt = [], [], [], [], []
    for (U, mw, ts, N) in fam:
        p = U.mean(0); W = U.shape[0] * mw; c = p * W
        J = window_len(c, p, N, thr)
        if J < 5: return abst("CL", "window")
        P.append(p[:J]); cnt.append(c[:J]); TS.append(ts[:J]); Ns.append(N); Wt.append(W)
    tsplit = np.sqrt(t0 * TS[0][-1])
    def both(Ps): return [_collapse_opt(Ps, cnt, Ns, TS, tsplit, r) for r in (0, 1)]
    est = both(P)
    if not np.all(np.isfinite(est)): return abst("CL", "nofit")
    boot = []
    for _ in range(B):
        Pb = []
        for (U, mw, ts, N), p0, Wi in zip(fam, P, Wt):
            n = U.shape[0]; m = rng.multinomial(n, np.full(n, 1.0 / n)); pb = (m @ U) / n
            Pb.append(np.maximum(pb[:len(p0)], 0.5 / Wi))
        b = both(Pb)
        if np.all(np.isfinite(b)): boot.append(b)
    boot = np.array(boot)
    if len(boot) < 0.8 * B: return abst("CL", "bootfail")
    dh = np.array(est); se = np.maximum(boot.std(0, ddof=1), 1e-9)
    c = np.quantile((np.abs(boot - dh) / se).max(1), 0.95)
    t = TS[0].astype(float)
    cen = [np.exp(np.log(t[t <= tsplit]).mean()), np.exp(np.log(t[t > tsplit]).mean())] if (t <= tsplit).any() and (t > tsplit).any() else None
    if cen is None: return abst("CL", "halves")
    dtr = truth["d_exact"](np.array(cen)); err = dh - dtr
    r = dict(est="CL", abstain=0, J=2, hw=float(c * se.mean()), se_mean=float(se.mean()), tsplit=float(tsplit), d_early=float(dh[0]), d_late=float(dh[1]),
             cov_simult=int((np.abs(err) <= c * se).all()), cov_point=float((np.abs(err) <= 1.96 * se).mean()), bias=float(err.mean()), mae=float(np.abs(err).mean()))
    if plateau is not None: r["cov_simult_plateau"] = int((np.abs(dh - plateau) <= c * se).all())
    dl = boot[:, 1] - boot[:, 0]; sD = max(float(dl.std(ddof=1)), 1e-9); D = dh[1] - dh[0]
    r.update(decl_neg=int(D + 1.96 * sD < -DELTA), decl_pos=int(D - 1.96 * sD > DELTA), decl_zero=int((D + 1.96 * sD < 0) or (D - 1.96 * sD > 0)), dmin=float(D), dmax=float(D))
    dt = dtr[1] - dtr[0]; r["delta_cov"] = int(D - 1.96 * sD <= dt <= D + 1.96 * sD); r["true_dmin"] = float(dt)
    return r
