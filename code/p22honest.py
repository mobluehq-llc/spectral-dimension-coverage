"""P22 Amendment 3: E4 bias-aware honest fixed-length intervals (Armstrong & Kolesar style), variants E4o / E4c / E4r. See PLAN.md AMENDMENT 3."""
import numpy as np
from scipy.optimize import minimize
import p22lib as L
from p22lib import grand_mean, boot_means, abst, window_len, DELTA

Z = 1.96
FLOORM = 1e-3
H_Q4 = (1.0, 1.4, 2.0, 2.8, 4.0)
STATS = {"fail": 0, "solves": 0}


def exact_g3_cummax(graph, ts):
    """max |g'''(u)| (g = log q, u = log t) over [log ts[0], log ts[k]] for each k, analytic from the spectrum."""
    ll = np.log(graph.lam[graph.lam > 1e-14])
    fine = np.unique(np.concatenate([np.exp(np.linspace(np.log(ts[0]), np.log(ts[-1]), 600)), ts.astype(float)]))
    g3 = np.empty(len(fine))
    for i, t in enumerate(fine):
        w = np.exp(t * (ll - ll.max())); w /= w.sum()
        k1 = (w * ll).sum(); c = ll - k1
        k2 = (w * c ** 2).sum(); k3 = (w * c ** 3).sum()
        g3[i] = t * k1 + 3 * t ** 2 * k2 + t ** 3 * k3
    a = np.abs(g3); cm = np.maximum.accumulate(a)
    return np.array([cm[np.searchsorted(fine, t, side="right") - 1] for t in ts.astype(float)])


def M_rule_of_thumb(x, y, w):
    P = np.polynomial.Polynomial.fit(x, y, 4, w=np.sqrt(w))
    xs = np.linspace(x[0], x[-1], 400)
    return max(float(np.abs(P.deriv(3)(xs)).max()), FLOORM)


def M_heldout(x, y, w):
    J = len(x); best = None
    rows = {}
    for h in H_Q4:
        val = np.zeros((J, J)); c3 = np.zeros((J, J)); ok = True
        for j in range(J):
            z = x - x[j]; u = np.exp(-0.5 * (z / h) ** 2) * w
            X = np.vander(z, 5, increasing=True); A = X.T @ (u[:, None] * X)
            try:
                Mx = np.linalg.solve(A, X.T * u[None, :])
            except np.linalg.LinAlgError:
                ok = False; break
            val[j] = Mx[0]; c3[j] = Mx[3]
        if not ok: continue
        hd = np.diag(val)
        if np.any(np.abs(1 - hd) < 1e-6): continue
        cv = (w * ((y - val @ y) / (1 - hd)) ** 2).sum()
        if np.isfinite(cv) and (best is None or cv < best[0]): best = (cv, h, c3)
    if best is None:
        return None
    return max(float(np.abs(6 * best[2] @ y).max()), FLOORM)


def honest_weights(u, uj, S, M):
    """argmin 1.96 sqrt(a'Sa) + (M/6) sum |a||u-uj|^3  s.t. sum a=0, sum a d=1, sum a d^2=0. returns (a, sd, B) or None."""
    J = len(u); d = u - uj; w3 = np.abs(d) ** 3 * (M / 6.0)
    Aeq = np.vstack([np.ones(J), d, d ** 2])
    beq = np.array([0.0, 1.0, 0.0])
    # sd-minimal feasible start
    Sr = S + 1e-12 * np.trace(S) / J * np.eye(J)
    Si = np.linalg.solve(Sr, Aeq.T)
    a0 = Si @ np.linalg.solve(Aeq @ Si, beq)
    sc = np.sqrt(np.trace(S) / J)                      # scale so objective ~O(1)
    Ss = S / sc ** 2; w3s = w3 / sc
    def f(v):
        a = v[:J] - v[J:]; s = np.sqrt(max(a @ Ss @ a, 1e-300))
        return Z * s + w3s @ (v[:J] + v[J:])
    def g(v):
        a = v[:J] - v[J:]; s = np.sqrt(max(a @ Ss @ a, 1e-300)); gs = Z * (Ss @ a) / s
        return np.concatenate([gs + w3s, -gs + w3s])
    v0 = np.concatenate([np.maximum(a0, 0), np.maximum(-a0, 0)])
    cons = [dict(type="eq", fun=lambda v: Aeq @ (v[:J] - v[J:]) - beq, jac=lambda v: np.hstack([Aeq, -Aeq]))]
    STATS["solves"] += 1
    best = None
    for start in (v0, None):
        if start is None:
            start = np.concatenate([np.maximum(a0, 0), np.maximum(-a0, 0)]) * 0.5 + 1e-3
            start = np.concatenate([np.maximum(a0, 0), np.maximum(-a0, 0)])
        res = minimize(f, start, jac=g, bounds=[(0, None)] * (2 * J), constraints=cons, method="SLSQP", options=dict(maxiter=300, ftol=1e-11))
        a = res.x[:J] - res.x[J:]
        if np.abs(Aeq @ a - beq).max() < 1e-6:
            best = (res, a); break
        STATS["fail"] += 1
        v0 = np.concatenate([np.maximum(a0, 0), np.maximum(-a0, 0)]) + 1e-2
    if best is None:
        return None
    res, a = best
    # polish: tiny entries do not matter; recompute exact sd, B from a
    sd = float(np.sqrt(max(a @ S @ a, 0))); B = float(np.abs(a) @ w3)
    return a, sd, B


def est_E4(Ublk, mwalk, ts, N, rng, truth, Mo_cum, thr=L.THR, B=200):
    """returns list [E4o, E4c, E4r] result dicts (same keys as p22lib.score where meaningful)."""
    tags = ["E4o", "E4c", "E4r"]
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA); WA = sum(u.shape[0] for u in UA) * mwalk; cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return [abst(t, "window") for t in tags]
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    Ms = [max(float(Mo_cum[J - 1]), FLOORM), M_heldout(x, yA, wA), M_rule_of_thumb(x, yA, wA)]
    WB = sum(u.shape[0] for u in UB) * mwalk; floor = 0.5 / WB
    yB = np.log(np.maximum(grand_mean(UB), floor)[:J])
    Ps = np.maximum(boot_means(UB, B, rng)[:, :J], floor); Ys = np.log(Ps)
    S = np.cov(Ys, rowvar=False)
    bt = truth["d_exact"](ts[:J]); bt = -bt / 2.0
    out = []
    for tag, M in zip(tags, Ms):
        if M is None:
            out.append(abst(tag, "cvfail")); continue
        W = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J); bad = False
        for j in range(J):
            r = honest_weights(x, x[j], S, M)
            if r is None: bad = True; break
            W[j], sd[j], Bj[j] = r
        if bad:
            out.append(abst(tag, "solver")); continue
        sd = np.maximum(sd, 1e-12)
        bh = W @ yB; Bs = Ys @ W.T
        c = np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95)
        hwj = c * sd + Bj                                   # half-width in beta units
        err = bh - bt
        dh = -2 * bh; hwd = 2 * hwj
        r = dict(est=tag, abstain=0, J=J, hw=float(hwd.mean()), se_mean=float(2 * sd.mean()), bias_bound=float(2 * Bj.mean()), M=float(M), crit=float(c),
                 cov_simult=int((np.abs(err) <= hwj).all()), cov_point=float((np.abs(err) <= Z * sd + Bj).mean()),
                 bias=float((-2 * err).mean()), mae=float(np.abs(2 * err).mean()), tlast=int(ts[J - 1]))
        # change test: contrasts a^(j) - a^(0)
        Cc = W[1:] - W[[0]]
        sdc = np.maximum(np.sqrt(np.einsum("ij,jk,ik->i", Cc, S, Cc)), 1e-12)
        Bc = Bj[1:] + Bj[0]
        dc = Cc @ yB; Dc = Ys @ Cc.T
        cc = np.quantile((np.abs(Dc - dc) / sdc).max(1), 0.95)
        hc = cc * sdc + Bc
        dl = -2 * dc; hcd = 2 * hc
        lo, hi = dl - hcd, dl + hcd
        r.update(decl_neg=int((hi < -DELTA).any()), decl_pos=int((lo > DELTA).any()), decl_zero=int(((hi < 0) | (lo > 0)).any()),
                 dmin=float(dl.min()), dmax=float(dl.max()))
        dtt = -2 * (bt[1:] - bt[0])
        r["delta_cov"] = int(((dtt >= lo) & (dtt <= hi)).all()); r["true_dmin"] = float(dtt.min())
        out.append(r)
    return out
