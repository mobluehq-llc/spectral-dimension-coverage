"""P22 Amendment 10b: certified OUTER bounds for E7 / E7j (continuous-support relaxation). See PLAN.md AMENDMENT 10b.

Two outer relaxations of the set of ALL measures nu on s in [0, inf] (q(t)=int e^{-st} dnu, nu(s=0) >= 1/N, total mass 1) that put q inside the box:

  A  "literal" (the plan's wording): replace every mass by its nearest grid atom; this moves e^{-s t_j} by at most delta_j (exact max over the cells),
     so the box is widened by delta_j * (total mass = 1). The objective d(t_j) = 2 t_j N/D is enclosed by bounding numerator N (error eps_j) and denominator D
     (error delta_j) separately: d <= 2 t (N' + eps)/(D' - delta), d >= 2 t (N' - eps)/(D' + delta); both are linear-fractional, i.e. one Charnes-Cooper LP each.
  B  "cell-sandwich" (same idea, per cell instead of absolute x total mass): a measure puts mass m_k on cell k = [a_k, b_k]; for every t_j,
     q_j in [sum m_k e^{-b_k t_j}, sum m_k e^{-a_k t_j}] and N_j = sum m_k * int s e^{-s t_j} in [sum m_k Gm_kj, sum m_k Gp_kj] (Gm/Gp = min/max of s e^{-st} on the cell).
     The box constraint lo<=q<=hi is relaxed to  sum m Ap >= lo  and  sum m Am <= hi  (linear). d_j <= 2 t_j sum m Gp / max(lo_j, sum m Am) (two LPs: one per denominator),
     d_j >= 2 t_j sum m Gm / min(hi_j, sum m Ap). The relaxation is independent per constraint, hence an OUTER set.
All constants (e^{-s t}, s e^{-st}, their cell extrema) are computed with python-flint arb balls (precision 80 bits) and rounded outward to doubles. The LPs are solved by
HiGHS in double precision (tolerance ~1e-7 on scaled rows): this is not a verified LP solve; the decision gate MARGIN of Amendment 9 is the guard against that.
E7j (B only): certified upper bound U on sup (d_j - d_0) over the relaxed set by branch and bound on d_0 (as Amendment 9), with the lower bound L taken from the ORIGINAL
(discrete-atom) LP, whose optimal measure is a genuine measure in the true set. A drop is declared iff U < -0.25; certified NOT declarable iff L >= -0.25; otherwise ambiguous (not declared)."""
import time
import numpy as np
from scipy.optimize import linprog
try:
    import flint
    from flint import arb, ctx
    ctx.prec = 80
    HAVE_FLINT = True
except Exception:       # pragma: no cover
    HAVE_FLINT = False
import p22lib as L
import p22amend8 as A8
import p22amend9 as A9

DELTA = L.DELTA
NCERT = 1597
_CC = {}


def _dn(x): return float(np.nextafter(x, -np.inf))
def _up(x): return float(np.nextafter(x, np.inf))


def _arb_exp_neg(a, t):
    """enclosure [lo, hi] (doubles, outward) of exp(-a t); a double >=0 (inf allowed)."""
    if np.isinf(a): return 0.0, 0.0
    if not HAVE_FLINT:
        v = np.exp(-a * t); return _dn(_dn(v) * (1 - 1e-14)), _up(v * (1 + 1e-14))
    x = (-(arb(float(a)) * int(t))).exp()
    return _dn(float(x.lower())), _up(float(x.upper()))


def _arb_g(a, t):
    """enclosure of s e^{-s t} at s=a."""
    if np.isinf(a) or a == 0.0: return 0.0, 0.0
    if not HAVE_FLINT:
        v = a * np.exp(-a * t); return _dn(v * (1 - 1e-14)), _up(v * (1 + 1e-14))
    x = arb(float(a)) * (-(arb(float(a)) * int(t))).exp()
    return _dn(float(x.lower())), _up(float(x.upper()))


def cell_constants(ts, n=NCERT):
    """cells Z={0}, C_1..C_n (geometric partition of (0, 1]), INF=[1, inf]. Returns dict of arrays (K x nc) and representative atoms."""
    key = (tuple(int(t) for t in ts), n)
    if key in _CC: return _CC[key]
    s = A8.atoms_grid(float(ts[-1]), n)             # [0, s1..sn, inf]
    sg = s[1:-1]; e = np.sqrt(sg[:-1] * sg[1:])
    nc = n + 2
    lo_e = np.zeros(nc); hi_e = np.zeros(nc)
    for k in range(1, n + 1):
        lo_e[k] = 0.0 if k == 1 else e[k - 2]
        hi_e[k] = e[k - 1] if k < n else 1.0
    lo_e[n + 1] = 1.0; hi_e[n + 1] = np.inf
    K = len(ts)
    Am = np.zeros((K, nc)); Ap = np.zeros((K, nc)); Gm = np.zeros((K, nc)); Gp = np.zeros((K, nc))
    rA = np.zeros((K, nc)); rG = np.zeros((K, nc)); dA = np.zeros((K, nc)); dG = np.zeros((K, nc))
    for j, t in enumerate(ts):
        t = int(t); pk = 1.0 / t
        pkhi = _up((np.exp(-1.0) / t) * (1 + 1e-14))
        for k in range(nc):
            a, b = lo_e[k], hi_e[k]
            Ap[j, k] = _arb_exp_neg(a, t)[1]; Am[j, k] = _arb_exp_neg(b, t)[0]
            ga = _arb_g(a, t); gb = _arb_g(b, t)
            Gm[j, k] = min(ga[0], gb[0]); gp = max(ga[1], gb[1])
            if a <= pk * (1 + 1e-9) and b >= pk * (1 - 1e-9): gp = max(gp, pkhi)
            Gp[j, k] = gp
            sk = s[k]
            r = _arb_exp_neg(sk, t); rg = _arb_g(sk, t)
            rA[j, k] = 0.5 * (r[0] + r[1]); rG[j, k] = 0.5 * (rg[0] + rg[1])
            dA[j, k] = max(Ap[j, k] - r[0], r[1] - Am[j, k], 0.0) + 1e-15
            dG[j, k] = max(Gp[j, k] - rg[0], rg[1] - Gm[j, k], 0.0) + 1e-15
    C = dict(Am=Am, Ap=Ap, Gm=Gm, Gp=Gp, rA=rA, rG=rG, delta=dA.max(1), eps=dG.max(1), s=s, nc=nc, ts=np.asarray(ts, float))
    _CC[key] = C
    return C


# ------------------------------------------------------------------ LP helpers
def _solve(c, Aub, bub, Aeq, beq):
    last = None
    for meth in ("highs", "highs-ds", "highs-ipm"):
        r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method=meth)
        last = r
        if r.status in (0, 2, 3): break
    return last


def _relaxed_rows(C, lo, hi, N):
    """rows of the relaxed set in variables [y (nc), tau]. returns Aub, bub, Aeq, beq (Aeq: sum y - tau = 0)."""
    ts = C["ts"]; nc = C["nc"]; K = len(ts)
    ref = np.maximum(0.5 * (lo + hi), 1e-300)
    r1 = np.hstack([-C["Ap"] / ref[:, None], (lo / ref)[:, None]])        # lo tau - Ap y <= 0
    r2 = np.hstack([C["Am"] / ref[:, None], -(hi / ref)[:, None]])        # Am y - hi tau <= 0
    r3 = np.zeros((1, nc + 1)); r3[0, 0] = -1.0; r3[0, nc] = 1.0 / N
    Aub = np.vstack([r1, r2, r3]); bub = np.zeros(Aub.shape[0])
    Aeq = np.zeros((1, nc + 1)); Aeq[0, :nc] = 1.0; Aeq[0, nc] = -1.0
    return Aub, bub, Aeq, np.zeros(1), ref


def _empty(Aub, bub, Aeq0, nc):
    Aeq = np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]])
    return _solve(np.zeros(nc + 1), Aub, bub, Aeq, np.array([0.0, 1.0])).status == 2


def bound_B(C, lo, hi, N, j, side, extra_rows=None, rows=None):
    """certified bound on d(t_j) over the relaxed set. side=+1: upper bound, -1: lower bound. returns (value, status); status 0 ok, 2 empty relaxed set, 1 solver trouble.
    extra_rows: (Aextra, bextra) appended to the inequality rows (variables [y, tau])."""
    nc = C["nc"]; t = C["ts"][j]
    Aub, bub, Aeq0, beq0, ref = rows if rows is not None else _relaxed_rows(C, lo, hi, N)
    if extra_rows is not None:
        Aub = np.vstack([Aub, extra_rows[0]]); bub = np.concatenate([bub, extra_rows[1]])
    vals = []
    for variant in (1, 2):
        if variant == 1:                                               # denominator sum m Am (upper) or sum m Ap (lower) normalised to 1
            den = (C["Am"][j] if side > 0 else C["Ap"][j]) / ref[j]
            Aeq = np.vstack([Aeq0, np.append(den, 0.0)[None, :]]); beq = np.array([0.0, 1.0])
            cvec = np.zeros(nc + 1); cvec[:nc] = 2.0 * t * (C["Gp"][j] if side > 0 else C["Gm"][j]) / ref[j]
            scale = 1.0
        else:                                                          # denominator replaced by the box bound; total mass 1
            bd = lo[j] if side > 0 else hi[j]
            if not (bd > 0): continue
            Aeq = np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]]); beq = np.array([0.0, 1.0])
            cvec = np.zeros(nc + 1); cvec[:nc] = 2.0 * t * (C["Gp"][j] if side > 0 else C["Gm"][j]) / bd
            scale = 1.0
        r = _solve(-side * cvec, Aub, bub, Aeq, beq)
        if r.status == 2:
            if variant == 2 or _empty(Aub, bub, Aeq0, nc):
                return np.nan, 2
            vals.append(np.nan); continue                              # set nonempty but the normalisation excludes it (mass only in the INF cell): skip this variant
        if r.status == 3:
            vals.append(side * np.inf if variant == 1 else np.nan); continue
        if r.status != 0:
            return np.nan, 1
        vals.append(-side * r.fun * scale)
    v = [x for x in vals if not np.isnan(x)]
    if not v: return (np.inf if side > 0 else 0.0), 0
    out = min(v) if side > 0 else max(v)
    return (max(out, 0.0) if side < 0 else out), 0


def band_B(ts, lo, hi, N, idx, n=NCERT):
    """certified band for d(t_j), j in idx. returns dmin, dmax arrays, status (0 ok, 2 empty relaxed set = certified infeasible, 1 solver trouble)."""
    C = cell_constants(ts, n); rows = _relaxed_rows(C, lo, hi, N)
    dmin = np.zeros(len(idx)); dmax = np.zeros(len(idx))
    for i, j in enumerate(idx):
        u, su = bound_B(C, lo, hi, N, j, +1, rows=rows)
        if su != 0: return np.full(len(idx), np.nan), np.full(len(idx), np.nan), su
        l, sl = bound_B(C, lo, hi, N, j, -1, rows=rows)
        if sl != 0: return np.full(len(idx), np.nan), np.full(len(idx), np.nan), sl
        dmin[i], dmax[i] = l, u
    return dmin, dmax, 0


# ------------------------------------------------------------------ method A (literal)
def band_A(ts, lo, hi, N, idx, n=NCERT):
    """plan-literal outer bound: nearest-atom, box widened by delta_j (absolute, total mass 1), objective enclosure by separate numerator/denominator errors."""
    C = cell_constants(ts, n); K = len(ts); nc = C["nc"]; delta = C["delta"]; eps = C["eps"]
    lo2 = np.maximum(lo - delta, 0.0); hi2 = hi + delta
    ref = np.maximum(0.5 * (lo2 + hi2), 1e-300); R = C["rA"]
    r1 = np.hstack([-R / ref[:, None], (lo2 / ref)[:, None]]); r2 = np.hstack([R / ref[:, None], -(hi2 / ref)[:, None]])
    r3 = np.zeros((1, nc + 1)); r3[0, 0] = -1.0; r3[0, nc] = 1.0 / N
    Aub = np.vstack([r1, r2, r3]); bub = np.zeros(Aub.shape[0])
    Aeq0 = np.zeros((1, nc + 1)); Aeq0[0, :nc] = 1.0; Aeq0[0, nc] = -1.0
    dmin = np.zeros(len(idx)); dmax = np.zeros(len(idx)); status = 0
    for i, j in enumerate(idx):
        t = C["ts"][j]; g = C["rG"][j]
        # emptiness of the widened set
        # upper: (N'+eps tau)/(D' - delta tau), D'-delta tau = ref_j (normalised)
        Aeq = np.vstack([Aeq0, np.append(R[j] / ref[j], -delta[j] / ref[j])[None, :]]); beq = np.array([0.0, 1.0])
        cvec = np.append(2 * t * g / ref[j], 2 * t * eps[j] / ref[j])
        r = _solve(-cvec, Aub, bub, Aeq, beq)
        if r.status == 2:
            # either the set is empty, or D' <= delta is forced: test emptiness
            rr = _solve(np.zeros(nc + 1), Aub, bub, np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]]), np.array([0.0, 1.0]))
            if rr.status == 2: return np.full(len(idx), np.nan), np.full(len(idx), np.nan), 2
            dmax[i] = np.inf
        elif r.status == 3: dmax[i] = np.inf
        elif r.status == 0:
            dmax[i] = -r.fun
            # can the denominator reach <= delta (then the bound is +inf)?
            rm = _solve(np.append(R[j], 0.0), Aub, bub, np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]]), np.array([0.0, 1.0]))
            if rm.status == 0 and rm.fun <= delta[j]: dmax[i] = np.inf
        else: return np.full(len(idx), np.nan), np.full(len(idx), np.nan), 1
        Aeq = np.vstack([Aeq0, np.append(R[j] / ref[j], delta[j] / ref[j])[None, :]])
        cvec = np.append(2 * t * g / ref[j], -2 * t * eps[j] / ref[j])
        r = _solve(cvec, Aub, bub, Aeq, beq)
        if r.status == 0: dmin[i] = max(r.fun, 0.0)
        elif r.status == 3: dmin[i] = 0.0
        elif r.status == 2: dmin[i] = 0.0
        else: return np.full(len(idx), np.nan), np.full(len(idx), np.nan), 1
    return dmin, dmax, 0


# ------------------------------------------------------------------ E7j certified
def _u_rows(C, ua, ub):
    t0 = C["ts"][0]; nc = C["nc"]
    r1 = np.zeros(nc + 1); r1[:nc] = -(2 * t0 * C["Gp"][0] - ua * C["Am"][0])      # 2t0 sum Gp y - ua sum Am y >= 0
    r2 = np.zeros(nc + 1); r2[:nc] = 2 * t0 * C["Gm"][0] - ub * C["Ap"][0]         # 2t0 sum Gm y - ub sum Ap y <= 0
    return np.vstack([r1, r2]), np.zeros(2)


def sup_delta_cert(S, C, lo, hi, N, rows, j, ulo, uhi, sigma, thr=DELTA, tol=A9.TOL, maxlp=A9.MAXLP):
    """bounds on sup over the TRUE set of sigma*(d_j - d_0): L attained (original LP measure), U certified (relaxed set). Returns dict(L, U, n_lp, nbad, stop, nL)."""
    t0 = S["ts"][0]; T_ = -thr
    state = dict(nlp=0, nbad=0)
    def evL(ua, ub):
        st, h, Fj, F0 = A9._lp(S, j, ua, ub, t0, sigma)
        state["nlp"] += 1
        if st == 0 and np.isfinite(F0): return sigma * (Fj - F0)
        return -np.inf
    def evU(ua, ub):
        ex = _u_rows(C, ua, ub)
        h, st = bound_B(C, lo, hi, N, j, sigma, extra_rows=ex, rows=rows)
        state["nlp"] += 2
        if st == 2: return -np.inf, 2
        if st != 0: state["nbad"] += 1; return None, st
        return ((h - ua) if sigma > 0 else (ub - h)), 0
    edges = np.linspace(ulo - 1e-9, uhi + 1e-9, A9.N0 + 1)
    ints = [(edges[i], edges[i + 1]) for i in range(A9.N0)]
    mid = 0.5 * (ulo + uhi); ints.sort(key=lambda iv: abs(0.5 * (iv[0] + iv[1]) - mid))
    Lb = -np.inf; live = []
    def Umax_(): return max([x[0] for x in live], default=-np.inf)
    for iv in ints:
        Lb = max(Lb, evL(*iv))
        if Lb >= T_: return dict(L=Lb, U=np.inf, n_lp=state["nlp"], nbad=state["nbad"], stop="L")
        U_, st = evU(*iv)
        if U_ is None: live.append((np.inf, iv, True))
        elif U_ > -np.inf: live.append((U_, iv, False))
    while live and state["nlp"] < maxlp:
        live.sort(key=lambda x: -x[0]); Umax = live[0][0]
        if Umax < T_ or Umax - Lb <= tol: break
        U_, (a, b), fin = live.pop(0)
        if fin or b - a < 1e-7:
            live.insert(0, (U_, (a, b), True)); break
        mdl = 0.5 * (a + b)
        for iv in ((a, mdl), (mdl, b)):
            Lb = max(Lb, evL(*iv))
            if Lb >= T_: return dict(L=Lb, U=np.inf, n_lp=state["nlp"], nbad=state["nbad"], stop="L")
            U2, st = evU(*iv)
            if U2 is None: live.append((U_, iv, True))
            elif U2 > -np.inf: live.append((min(U2, U_), iv, False))
    Um = Umax_()
    stop = "U" if Um < T_ else ("tol" if Um - Lb <= tol else ("unres" if any(x[2] for x in live) else "maxlp"))
    return dict(L=Lb, U=Um, n_lp=state["nlp"], nbad=state["nbad"], stop=stop)


def run_cert(ts, lo, hi, N, Jw, dtrue, n=NCERT, do_j=True, do_A=True):
    """one dataset: certified bands (A and B) for the window, original gated E7 band, certified E7j decisions. Returns dict of flat fields."""
    out = {}
    idx = list(range(Jw)); C = cell_constants(ts, n)
    t1 = time.time(); o_lo, o_hi, o_st, o_info = A8.e7_gated(ts, lo, hi, N, idx); out["t_orig"] = time.time() - t1
    out.update(o_status=o_st, o_natoms=o_info["natoms"], o_gate_fail=o_info["gate_fail"])
    t1 = time.time(); b_lo, b_hi, b_st = band_B(ts, lo, hi, N, idx, n); out["t_B"] = time.time() - t1
    out["B_status"] = b_st
    out["o_lo"] = ";".join("%.6f" % v for v in o_lo) if o_st == 0 else ""; out["o_hi"] = ";".join("%.6f" % v for v in o_hi) if o_st == 0 else ""
    out["B_lo"] = ";".join("%.6f" % v for v in b_lo) if b_st == 0 else ""; out["B_hi"] = ";".join("%.6f" % v for v in b_hi) if b_st == 0 else ""
    if o_st == 0:
        out["o_cov"] = int(((dtrue >= o_lo - 1e-9) & (dtrue <= o_hi + 1e-9)).all()); out["o_width"] = float((o_hi - o_lo).mean())
    if b_st == 0:
        fin = np.isfinite(b_hi)
        out["B_cov"] = int(((dtrue >= b_lo - 1e-9) & (dtrue <= b_hi + 1e-9)).all()); out["B_width"] = float((b_hi - b_lo)[fin].mean()) if fin.any() else np.nan
        out["B_nopen"] = int((~fin).sum())
        if o_st == 0: out["B_contains_orig"] = int(((b_lo <= o_lo + 1e-6) & (b_hi >= o_hi - 1e-6)).all()); out["B_max_inward"] = float(max((b_lo - o_lo).max(), (o_hi - b_hi).max()))
        if Jw >= 2:
            dl_lo = b_lo[1:] - b_hi[0]; dl_hi = b_hi[1:] - b_lo[0]
            out["B_e7_neg"] = int((dl_hi < -DELTA).any()); out["B_e7_pos"] = int((dl_lo > DELTA).any())
    if o_st == 0 and Jw >= 2:
        out["o_e7_neg"] = int(((o_hi[1:] - o_lo[0]) < -DELTA).any()); out["o_e7_pos"] = int(((o_lo[1:] - o_hi[0]) > DELTA).any())
    if do_A:
        t1 = time.time(); a_lo, a_hi, a_st = band_A(ts, lo, hi, N, idx, n); out["t_A"] = time.time() - t1
        out["A_status"] = a_st
        if a_st == 0:
            fin = np.isfinite(a_hi)
            out["A_cov"] = int(((dtrue >= a_lo - 1e-9) & (dtrue <= a_hi + 1e-9)).all()); out["A_nopen"] = int((~fin).sum()); out["A_width"] = float((a_hi - a_lo)[fin].mean()) if fin.any() else np.nan
            if o_st == 0: out["A_contains_orig"] = int(((a_lo <= o_lo + 1e-6) & (a_hi >= o_hi - 1e-6)).all())
            if Jw >= 2:
                out["A_e7_neg"] = int(((a_hi[1:] - a_lo[0]) < -DELTA).any()); out["A_e7_pos"] = int(((a_lo[1:] - a_hi[0]) > DELTA).any())
    if do_j and b_st == 0 and Jw >= 2:
        t1 = time.time()
        S = A8.lp_setup(ts, lo, hi, A8.atoms_grid(float(ts[-1]), n), N); rows = _relaxed_rows(C, lo, hi, N)
        ulo, uhi = float(b_lo[0]), float(b_hi[0])
        if not np.isfinite(uhi): uhi = max(float(b_hi[0]), 10.0)
        dec = {+1: 0, -1: 0}; nlp = 0; nbad = 0; amb = 0; cnon = 0; Umin = {+1: np.inf, -1: np.inf}; Lmax = {+1: -np.inf, -1: -np.inf}
        for j in range(1, Jw):
            for sg in (+1, -1):
                r = sup_delta_cert(S, C, lo, hi, N, rows, j, ulo, uhi, sg)
                nlp += r["n_lp"]; nbad += r["nbad"]
                if r["U"] < -DELTA: dec[sg] = 1
                elif r["stop"] == "L": cnon += 1
                else: amb += 1
                Umin[sg] = min(Umin[sg], r["U"]); Lmax[sg] = max(Lmax[sg], r["L"])
        out.update(Bj_neg=dec[+1], Bj_pos=dec[-1], Bj_amb=amb, Bj_certnon=cnon, Bj_nlp=nlp, Bj_nbad=nbad, Bj_minU_drop=Umin[+1] if np.isfinite(Umin[+1]) else np.nan,
                   Bj_minU_rise=Umin[-1] if np.isfinite(Umin[-1]) else np.nan, t_Bj=time.time() - t1)
    return out
