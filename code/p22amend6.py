"""P22 Amendment 6: E6 + E4 per-point bands, alarm A1 (E6 vs E4c band disagreement), fallback A2 (one-step vs two-step test on split A). See PLAN.md AMENDMENT 6."""
import numpy as np
from scipy.stats import chi2
import p22lib as L
import p22shape as S
import p22honest as H
import p22amend5 as A5
from p22lib import grand_mean, boot_means, abst, window_len, DELTA


def e6_bands(x, yA, SA, yB, Ys, truth_d):
    """Identical computation to p22amend5.e6_lof; additionally returns the per-point bands (cLo, cHi) and the shrunk half-A precision."""
    r = A5.e6_lof(x, yA, SA, yB, Ys, truth_d)
    J = len(x)
    SAs = (1 - S.SHRINK) * SA + S.SHRINK * np.diag(np.diag(SA)) + 1e-10 * np.trace(SA) / J * np.eye(J)
    SAi = np.linalg.inv(SAs)
    SB = np.cov(Ys, rowvar=False); SBs = SB + 1e-10 * np.trace(SB) / J * np.eye(J); SBi = np.linalg.inv(SBs)
    ms = np.linspace(x[0] - 1.0, x[-1] + 1.0, S.M_GRID_N); base = np.stack([np.ones(J), x], 1)
    rss = []; fits = []
    for m in ms:
        for w in S.W_GRID:
            phi = 0.5 * w * S.softplus((x - m) / w); X = np.column_stack([base, phi])
            Kx = np.linalg.solve(X.T @ SAi @ X, X.T @ SAi); res = yA - X @ (Kx @ yA)
            rss.append(float(res @ SAi @ res)); fits.append((m, w))
    rss = np.array(rss); keep = rss - rss.min() <= S.CHI2_2_99
    cLo = np.full(J, np.inf); cHi = np.full(J, -np.inf); cr = np.sqrt(S.CHI2_2_96)
    for (m, w), k in zip(fits, keep):
        if not k: continue
        phi = 0.5 * w * S.softplus((x - m) / w); X = np.column_stack([base, phi])
        K = np.linalg.solve(X.T @ SBi @ X, X.T @ SBi)[1:]
        th = K @ yB; Sg = K @ SB @ K.T; s = S.sigma((x - m) / w)
        V = np.stack([np.full(J, -2.0), -s], 1)
        cen = V @ th; sdv = np.sqrt(np.maximum(np.einsum("ji,ik,jk->j", V, Sg, V), 0))
        cLo = np.minimum(cLo, cen - cr * sdv); cHi = np.maximum(cHi, cen + cr * sdv)
    assert abs(float(np.mean(cHi - cLo) / 2) - r["hw"]) < 1e-9, "band recomputation differs from e6_lof"
    return r, cLo, cHi, SAi


def a2_test(x, yA, SAi, B=200):
    """One-step vs two-step on half A (GLS). Lambda = Q1 - Q2 >= 0; alarm iff Lambda > chi2_{3,.95}(B-1)/(B-J)."""
    J = len(x)
    Lc = np.linalg.cholesky(SAi)                       # SAi = Lc Lc^T ; whitening matrix is Lc^T
    base = np.stack([np.ones(J), x], 1)
    Bw = Lc.T @ base; P = np.eye(J) - Bw @ np.linalg.solve(Bw.T @ Bw, Bw.T)
    yt = P @ (Lc.T @ yA)
    ms = np.linspace(x[0] - 1.0, x[-1] + 1.0, S.M_GRID_N)
    phis = []; mi = []
    for i, m in enumerate(ms):
        for w in S.W_GRID:
            phis.append(0.5 * w * S.softplus((x - m) / w)); mi.append(i)
    Ph = P @ (Lc.T @ np.array(phis).T)                  # J x K, residualised
    mi = np.array(mi); yy = float(yt @ yt)
    G = Ph.T @ Ph; r = Ph.T @ yt; dg = np.diag(G)
    Q1 = yy - r ** 2 / np.maximum(dg, 1e-300)
    q1 = float(Q1.min())
    a = dg[:, None]; b = dg[None, :]; c = G
    det = a * b - c ** 2
    ok = (mi[:, None] < mi[None, :]) & (det > 1e-8 * a * b)
    with np.errstate(divide="ignore", invalid="ignore"):
        rss2 = yy - (b * r[:, None] ** 2 - 2 * c * r[:, None] * r[None, :] + a * r[None, :] ** 2) / det
    q2 = float(min(q1, np.where(ok, rss2, np.inf).min()))
    lam = max(q1 - q2, 0.0)
    thr = float(chi2.ppf(0.95, 3) * (B - 1) / (B - J))
    return dict(a2_lam=lam, a2_thr=thr, a2_alarm=int(lam > thr), a2_q1=q1, a2_q2=q2)


def e4_bands(UB, UA, mwalk, ts, N, rng, truth, J, pA_cnt, tags, B=200):
    """E4c/E4r only. Same computation as p22honest.est_E4 for those tags; also returns d-band arrays."""
    x = np.log(ts[:J].astype(float)); pA = grand_mean(UA)
    yA = np.log(pA[:J]); wA = pA_cnt[:J]
    Msel = dict(E4c=H.M_heldout(x, yA, wA), E4r=H.M_rule_of_thumb(x, yA, wA))
    WB = sum(u.shape[0] for u in UB) * mwalk; floor = 0.5 / WB
    yB = np.log(np.maximum(grand_mean(UB), floor)[:J])
    Ps = np.maximum(boot_means(UB, B, rng)[:, :J], floor); Ys = np.log(Ps); Sm = np.cov(Ys, rowvar=False)
    bt = -truth["d_exact"](ts[:J]) / 2.0
    out = {}
    for tag in tags:
        M = Msel[tag]
        if M is None: out[tag] = (abst(tag, "cvfail"), None, None); continue
        W = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J); bad = False
        for j in range(J):
            rr = H.honest_weights(x, x[j], Sm, M)
            if rr is None: bad = True; break
            W[j], sd[j], Bj[j] = rr
        if bad: out[tag] = (abst(tag, "solver"), None, None); continue
        sd = np.maximum(sd, 1e-12); bh = W @ yB; Bs = Ys @ W.T
        c = np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95); hwj = c * sd + Bj; err = bh - bt
        r = dict(est=tag, abstain=0, J=J, hw=float(2 * hwj.mean()), se_mean=float(2 * sd.mean()), bias_bound=float(2 * Bj.mean()), M=float(M), crit=float(c),
                 cov_simult=int((np.abs(err) <= hwj).all()), cov_point=float((np.abs(err) <= H.Z * sd + Bj).mean()), tlast=int(ts[J - 1]))
        Cc = W[1:] - W[[0]]; sdc = np.maximum(np.sqrt(np.einsum("ij,jk,ik->i", Cc, Sm, Cc)), 1e-12); Bc = Bj[1:] + Bj[0]
        dc = Cc @ yB; Dc = Ys @ Cc.T; cc = np.quantile((np.abs(Dc - dc) / sdc).max(1), 0.95); hc = cc * sdc + Bc
        dl = -2 * dc; hcd = 2 * hc; lo, hi = dl - hcd, dl + hcd
        r.update(decl_neg=int((hi < -DELTA).any()), decl_pos=int((lo > DELTA).any()), decl_zero=int(((hi < 0) | (lo > 0)).any()), dmin=float(dl.min()), dmax=float(dl.max()))
        dtt = -2 * (bt[1:] - bt[0]); r["delta_cov"] = int(((dtt >= lo) & (dtt <= hi)).all()); r["true_dmin"] = float(dtt.min())
        out[tag] = (r, -2 * bh - 2 * hwj, -2 * bh + 2 * hwj)
    return out


def est_all(Ublk, mwalk, ts, N, rngB, rngA, rngE4, truth, thr, e4tags=("E4c",), do_e6=True):
    """returns list of rows. E6 row carries a1/a2 alarms (when E4c present)."""
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA); WA = sum(u.shape[0] for u in UA) * mwalk; cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        rows = [abst("E6", "window")] if do_e6 else []
        return rows + [abst(t, "window") for t in e4tags]
    rows = []; e6lo = e6hi = None
    if do_e6:
        x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J])
        WB = sum(u.shape[0] for u in UB) * mwalk; floor = 0.5 / WB
        yB = np.log(np.maximum(grand_mean(UB), floor)[:J])
        Ys = np.log(np.maximum(boot_means(UB, 200, rngB)[:, :J], floor))
        PsA = np.maximum(boot_means(UA, 200, rngA)[:, :J], 0.5 / WA)
        SA = np.cov(np.log(PsA), rowvar=False)
        r6, e6lo, e6hi, SAi = e6_bands(x, yA, SA, yB, Ys, truth["d_exact"](ts[:J])); r6["tlast"] = int(ts[J - 1])
        r6.update(a2_test(x, yA, SAi))
        r6["e6lo"] = ";".join("%.5f" % v for v in e6lo); r6["e6hi"] = ";".join("%.5f" % v for v in e6hi)
        rows.append(r6)
    e4 = e4_bands(UB, UA, mwalk, ts, N, rngE4, truth, J, cntA, list(e4tags)) if e4tags else {}
    for tag, (r, lo, hi) in e4.items():
        if lo is not None:
            r["e4lo"] = ";".join("%.5f" % v for v in lo); r["e4hi"] = ";".join("%.5f" % v for v in hi)
            if tag == "E4c" and do_e6:
                gap = np.maximum(lo - e6hi, e6lo - hi)             # >0 means disjoint at that grid point
                rows[0]["a1_alarm"] = int((gap > 0).any()); rows[0]["a1_maxgap"] = float(gap.max())
        rows.append(r)
    if do_e6 and "a1_alarm" not in rows[0]:
        rows[0]["a1_alarm"] = 0; rows[0]["a1_maxgap"] = np.nan; rows[0]["a1_e4_abstain"] = 1
    return rows
