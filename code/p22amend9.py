"""P22 Amendment 9: (9a) E7j joint change test; (9c) second-design patches. See PLAN.md AMENDMENT 9.
Imports the Amendment 8 code (frozen at git 1ba9ca6) unchanged.

9a method (exact up to the stated tolerance; the plan's "bisection with LP feasibility" cannot be used literally because Delta_j = d(t_j) - d(t_0)
is a difference of two linear-fractional terms with DIFFERENT denominators, so its feasibility set is not an LP for fixed value):
  Delta^max_j = max_u [ h_j(u) - u ],  h_j(u) = max d(t_j) over measures in the E7 set with d(t_0) = u.
  For fixed u both d(t_0) = u (2 t0 a0 - u q0 = 0, homogeneous) and the Charnes-Cooper objective are linear => one LP.
  The 1-D outer problem over u in [dmin_0, dmax_0] (E7's own band at t_0) is solved by RIGOROUS branch and bound: on a u-interval [ua, ub] the LP with
  ua <= d(t_0) <= ub gives h; then sup Delta <= h - ua (upper bound) and the LP's optimal measure gives an attained Delta (lower bound).
  A drop is declared at j iff the upper bound over all intervals is < -0.25 (certified); not declared as soon as an attained measure has Delta >= -0.25.
  Rise symmetric (sigma = -1: sup of -Delta, using the min-LP).  Tolerance of B&B 0.002 on U-L."""
import time
import numpy as np
from scipy.optimize import linprog
import p22lib as L
import p22amend8 as A8

DELTA = L.DELTA
TOL = 0.002
N0 = 6              # initial u-intervals
MAXLP = 400
MARGIN = 0.005      # gate: decision margin below which the finer atom grid is re-run


def _lp(S, j, ua, ub, t0, sigma):
    """LP on interval. returns (status, h, F_j(w), F_0(w)). sigma=+1: max F_j; sigma=-1: min F_j. Constrained: ua <= F_0 <= ub."""
    n = S["n"]; As = S["As"]; sf = S["sfin"]; ts = S["ts"]
    g0 = 2.0 * t0 * sf * As[0]; a0 = As[0]
    r1 = np.zeros(n + 1); r1[:n] = ua * a0 - g0
    r2 = np.zeros(n + 1); r2[:n] = g0 - ub * a0
    Aub = np.vstack([S["Aub"], r1, r2]); bub = np.concatenate([S["bub"], [0.0, 0.0]])
    Aeq = np.zeros((2, n + 1)); beq = np.array([1.0, 0.0])
    Aeq[0, :n] = As[j]; Aeq[1, :n] = 1.0; Aeq[1, n] = -1.0
    cv = np.zeros(n + 1); cv[:n] = 2.0 * ts[j] * sf * As[j]
    for meth in ("highs", "highs-ipm", "highs-ds"):
        r = linprog(-sigma * cv, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method=meth)
        if r.status in (0, 2): break
    if r.status == 2:
        return 2, np.nan, np.nan, np.nan
    if r.status != 0:
        return int(r.status), np.nan, np.nan, np.nan
    y = r.x[:n]
    qj = As[j] @ y; q0 = a0 @ y
    Fj = 2.0 * ts[j] * (sf * As[j]) @ y / qj
    F0 = 2.0 * t0 * (sf * a0) @ y / q0 if q0 > 0 else np.nan
    return 0, Fj, Fj, F0


def sup_delta(S, j, ulo, uhi, sigma, thr=DELTA, tol=TOL, maxlp=MAXLP):
    """Certified bounds on sup over the E7 set of sigma*(d_j - d_0). Returns dict(L, U, n_lp, nbad, stop).
    Intervals are (U, a, b, final): an interval whose LP stays numerically unresolved after three solvers inherits its parent's bound and is not split further."""
    t0 = S["ts"][0]; T_ = -thr if thr is not None else np.nan      # thr=None: no decision stop, run to tolerance
    def ev(ua, ub):
        st, h, Fj, F0 = _lp(S, j, ua, ub, t0, sigma)
        if st == 2: return -np.inf, -np.inf, 2
        if st != 0: return None, -np.inf, st
        ub_ = (h - ua) if sigma > 0 else (ub - h)
        lb_ = sigma * (Fj - F0) if np.isfinite(F0) else -np.inf
        return ub_, lb_, 0
    edges = np.linspace(ulo - 1e-9, uhi + 1e-9, N0 + 1)
    ints = [(edges[i], edges[i + 1]) for i in range(N0)]
    mid = 0.5 * (ulo + uhi)
    ints.sort(key=lambda iv: abs(0.5 * (iv[0] + iv[1]) - mid))
    Lb = -np.inf; nlp = 0; nbad = 0; live = []
    def Umax_():
        return max([x[0] for x in live], default=-np.inf)
    for iv in ints:
        U_, L_, st = ev(*iv); nlp += 1
        if st not in (0, 2): nbad += 1
        Lb = max(Lb, L_)
        if U_ is None: live.append((np.inf, iv, True))
        elif U_ > -np.inf: live.append((U_, iv, False))
        if Lb >= T_: return dict(L=Lb, U=Umax_(), n_lp=nlp, nbad=nbad, stop="L")
    while live and nlp < maxlp:
        live.sort(key=lambda x: -x[0]); Umax = live[0][0]
        if Umax < T_ or Umax - Lb <= tol: break
        U_, (a, b), fin = live.pop(0)
        if fin or b - a < 1e-7:
            live.insert(0, (U_, (a, b), True)); break
        mdl = 0.5 * (a + b)
        for iv in ((a, mdl), (mdl, b)):
            U2, L2, st = ev(*iv); nlp += 1
            if st not in (0, 2): nbad += 1
            Lb = max(Lb, L2)
            if U2 is None: live.append((U_, iv, True))
            elif U2 > -np.inf: live.append((min(U2, U_), iv, False))
        if Lb >= T_: return dict(L=Lb, U=Umax_(), n_lp=nlp, nbad=nbad, stop="L")
    Um = Umax_()
    stop = "U" if Um < T_ else ("tol" if Um - Lb <= tol else ("unres" if any(x[2] for x in live) else "maxlp"))
    return dict(L=Lb, U=Um, n_lp=nlp, nbad=nbad, stop=stop)


def e7j_decide(ts, lo, hi, N, Jw, natoms, thr=DELTA):
    """Joint change test on the box (lo, hi) with the E7 atom grid of size natoms. Returns dict of decisions and diagnostics.
    Uses E7's own d(t_0) band for the u-range; re-checks decisions with margin < MARGIN on a grid with twice the atoms (gate)."""
    tmax = float(ts[-1]); out = dict(j_infeasible=0, drop=0, rise=0, gate_flag=0, gate_rechecked=0, n_lp=0, nbad=0, n_amb=0, minU_drop=np.nan, minU_rise=np.nan,
                                     maxL_drop=np.nan, maxL_rise=np.nan)
    def build(n):
        return A8.lp_setup(ts, lo, hi, A8.atoms_grid(tmax, n), N)
    S = build(natoms)
    d0min, d0max, st = A8.lp_one(S, 0)
    if st != 0:
        out["j_infeasible"] = 1; return out
    S2 = [None]
    Us = {+1: [], -1: []}; Ls = {+1: [], -1: []}
    for j in range(1, Jw):
        for sg in (+1, -1):
            r = sup_delta(S, j, d0min, d0max, sg, thr)
            out["n_lp"] += r["n_lp"]; out["nbad"] += r["nbad"]
            dec = int(r["U"] < -thr)                               # certified declaration
            amb = int(r["stop"] in ("tol", "maxlp", "unres"))              # straddling thr at tolerance: not declared
            out["n_amb"] += amb
            marg = (-thr - r["U"]) if dec else (r["L"] + thr if r["stop"] == "L" else 0.0)
            if abs(marg) < MARGIN or amb:                           # gate: recheck at twice the atoms
                if S2[0] is None: S2[0] = build(2 * (natoms - 1) + 1)
                a0, a1, st2 = A8.lp_one(S2[0], 0)
                out["gate_rechecked"] += 1
                if st2 != 0: out["gate_flag"] += 1; dec = 0
                else:
                    r2 = sup_delta(S2[0], j, a0, a1, sg, thr); out["n_lp"] += r2["n_lp"]; out["nbad"] += r2["nbad"]
                    dec2 = int(r2["U"] < -thr)
                    if dec2 != dec: out["gate_flag"] += 1; dec = 0
            if dec:
                if sg > 0: out["drop"] = 1
                else: out["rise"] = 1
            Us[sg].append(r["U"]); Ls[sg].append(r["L"])
    if Us[1]: out["minU_drop"] = float(min(Us[1])); out["maxL_drop"] = float(max(Ls[1]))
    if Us[-1]: out["minU_rise"] = float(min(Us[-1])); out["maxL_rise"] = float(max(Ls[-1]))
    return out


def run_e7j(U, m, ts, N, rng, thr_w, J6=0):
    """Same box and E7 band as A8.run_e7_e8 (same rng usage), then the joint test. Returns one row dict (est E7j) plus the E7 interval-arithmetic decisions."""
    nb, K = U.shape; Wt = nb * m
    p = U.mean(0); cnt = p * Wt
    Jw = L.window_len(cnt, p, N, thr_w)
    if Jw < 2: return dict(est="E7j", abstain=0, nowin=1, J=Jw)
    qhat, lo, hi, Ps, binfo = A8.make_box(U, m, rng)
    t0 = time.time()
    dmin, dmax, st, ginfo = A8.e7_gated(ts, lo, hi, N, list(range(Jw)))
    tlp = time.time() - t0
    r = dict(est="E7j", abstain=0, nowin=0, J=Jw, e7_status=st, e7_natoms=ginfo["natoms"], e7_gate_fail=ginfo["gate_fail"], lp_seconds=tlp)
    if st != 0:
        r.update(infeasible=1, decl_neg=0, decl_pos=0, decl_any=0, e7_decl_neg=0, e7_decl_pos=0); return r
    dl_lo = dmin[1:Jw] - dmax[0]; dl_hi = dmax[1:Jw] - dmin[0]
    r["e7_decl_neg"] = int((dl_hi < -DELTA).any()); r["e7_decl_pos"] = int((dl_lo > DELTA).any())
    t1 = time.time()
    j = e7j_decide(ts, lo, hi, N, Jw, ginfo["natoms"])
    r.update(infeasible=int(j["j_infeasible"]), decl_neg=j["drop"], decl_pos=j["rise"], decl_any=int(j["drop"] or j["rise"]), j_seconds=time.time() - t1,
             **{k: j[k] for k in ("gate_flag", "gate_rechecked", "n_lp", "nbad", "n_amb", "minU_drop", "minU_rise", "maxL_drop", "maxL_rise")})
    r["viol_vs_e7"] = int(r["decl_neg"] < r["e7_decl_neg"] or r["decl_pos"] < r["e7_decl_pos"])   # E7j must declare whenever E7 does
    return r


# ---------------------------------------------------------------- 9c patches
def make_grid_quarter(a0, b0, tcap):
    """grid {a0 * 2^(j/4)} rounded to integer steps (b0 ignored), <= tcap."""
    ts = []; j = 0
    while a0 * 2 ** (j / 4) <= tcap + 0.5:
        v = int(round(a0 * 2 ** (j / 4)))
        if v <= tcap: ts.append(v)
        j += 1
    return np.array(sorted(set(ts)), dtype=np.int64)


def patch_design2(cmin=100, jmin=8):
    """Strict window: count >= cmin (via the default arg and L.CMIN) and J >= jmin else the window length is reported as 0 (=> abstain / nowin)."""
    import sys
    orig = L.window_len.__wrapped__ if hasattr(L.window_len, "__wrapped__") else L.window_len
    def wl(cnt, p, N, thr, cmin_=cmin):
        J = orig(cnt, p, N, thr, cmin_)
        return J if J >= jmin else 0
    wl.__wrapped__ = orig
    L.CMIN = cmin
    L.window_len = wl
    for name, mod in list(sys.modules.items()):
        if name.startswith("p22") and mod is not None and getattr(mod, "window_len", None) is orig:
            mod.window_len = wl
    L.make_grid = make_grid_quarter
    return wl
