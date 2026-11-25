#!/usr/bin/env python3
"""SA3 pick 1: continuum certificate of the stored 2D global declarations. Pre-registration: PREREG_REDTEAM_FINAL.md amendment 11-24 PM (P7).
Usage: rt_certify.py unit1 [ds] | run [--ncpu K] | sentinels [--ncpu K]"""
import os, sys, time, collections
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np, pandas as pd
import multiprocessing as mp
from scipy.optimize import linprog
import rt_nodrop as RT
import rt_mindrop as MD
import p22lib as L, p22amend8 as A8
import run_p22 as R
sys.argv = _argv
args = sys.argv[1:]
DELTA = 0.25; MARG = RT.MARG; PAD = 0.005
PHI_MIN = 1e-8; MAX_ROUNDS = 12; INFL = 1e-12; HMIN = 1e-13


# ---------------------------------------------------------------- cell LPs
class Cell:
    def __init__(self, ra, rb, thr, name):
        self.ra, self.rb, self.thr, self.name = ra, rb, thr, name   # ra/rb None = absent; thr None = no no-drop rows


class Ctx:
    def __init__(self, ts, lo, hi, N, Jw):
        self.ts, self.lo, self.hi, self.N, self.Jw = np.asarray(ts, float), lo, hi, N, Jw
        self.K = len(ts); self.ref = np.maximum(0.5 * (lo + hi), 1e-300)
        self.extra = []                                     # exchange atoms shared across the cells of one dataset
        self.tmax = float(ts[-1])

    def atoms(self):
        g = A8.atoms_grid(self.tmax, RT.NATOMS)             # [0, grid..., inf]
        mid = np.concatenate([g[1:-1], np.array(self.extra, float)]); mid = np.unique(mid)
        return np.concatenate([[0.0], mid, [np.inf]])


def cell_rows(ctx, S, cell):
    n = S["n"]; As = S["As"]; sf = S["sfin"]; ts = S["ts"]
    g = lambda j: 2.0 * ts[j] * sf * As[j]
    rows = []
    if cell.ra is not None: r = np.zeros(n + 1); r[:n] = cell.ra * As[0] - g(0); rows.append(r)
    if cell.rb is not None: r = np.zeros(n + 1); r[:n] = g(0) - cell.rb * As[0]; rows.append(r)
    if cell.thr is not None:
        for j in range(1, ctx.Jw):
            r = np.zeros(n + 1); r[:n] = cell.thr * As[j] - g(j); rows.append(r)
    return np.array(rows) if rows else np.zeros((0, n + 1))


def elastic(ctx, cell, atoms):
    """phase-1 LP. returns dict(status, phi, u, v, S, M, nbase)"""
    S = A8.lp_setup(ctx.ts, ctx.lo, ctx.hi, atoms, ctx.N); n = S["n"]
    M = np.vstack([S["Aub"], cell_rows(ctx, S, cell)]); m = M.shape[0]
    A_ub = np.hstack([M, -np.eye(m)])
    Aeq = np.zeros((2, n + 1 + m)); Aeq[0, :n] = 1.0; Aeq[0, n] = -1.0; Aeq[1, n] = 1.0; beq = np.array([0.0, 1.0])
    c = np.concatenate([np.zeros(n + 1), np.ones(m)])
    r = None
    for meth in ("highs-ds", "highs-ipm", "highs"):
        r = linprog(c, A_ub=A_ub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=(0, None), method=meth)
        if r.status == 0: break
    if r.status != 0: return dict(status="solver", S=S)
    u = np.maximum(-np.asarray(r.ineqlin.marginals), 0.0); v = -np.asarray(r.eqlin.marginals)
    return dict(status="ok", phi=float(r.fun), u=u, v=v, S=S, M=M, y=r.x[:n])


def coeffs(ctx, cell, u):
    """P_j, Q_j (j over all K times), u_N, tau-column sum  from the multiplier vector u (rows: lo(K), hi(K), N, extras)"""
    K = ctx.K; ts = ctx.ts; ref = ctx.ref
    P = np.zeros(K); Q = np.zeros(K)
    P += -u[:K] / ref + u[K:2 * K] / ref
    uN = u[2 * K]; k = 2 * K + 1
    taus = float((u[:K] * ctx.lo / ref).sum() - (u[K:2 * K] * ctx.hi / ref).sum() + uN / ctx.N)
    if cell.ra is not None:
        w = u[k]; k += 1; P[0] += w * cell.ra / ref[0]; Q[0] += -2 * ts[0] * w / ref[0]
    if cell.rb is not None:
        w = u[k]; k += 1; P[0] += -w * cell.rb / ref[0]; Q[0] += 2 * ts[0] * w / ref[0]
    if cell.thr is not None:
        for j in range(1, ctx.Jw):
            w = u[k]; k += 1; P[j] += w * cell.thr / ref[j]; Q[j] += -2 * ts[j] * w / ref[j]
    assert k == len(u)
    return P, Q, float(uN), taus


# ---------------------------------------------------------------- interval B&B over s in (0, 1] (float backend, 1e-12 relative inflation)
def gvals(P, Q, ts, s):
    e = np.exp(-np.outer(s, ts))
    lin = P[None, :] + np.outer(s, Q)
    f0 = lin * e
    f1 = (Q[None, :] - ts[None, :] * lin) * e
    return f0, f1, e


def bnb(P, Q, ts, c, thr, maxcells=200000):
    """lower-bound g(s) = c + sum_j (P_j + Q_j s) exp(-t_j s) over s in [0, 1]. returns (min_lb, violators[(s,g)], unresolved, ncell).
    pruned cells need lb >= thr. min_lb = min over pruned-cell lower bounds (with inflation)."""
    edges = np.concatenate([[0.0], np.geomspace(1e-7, 1.0, 90)])
    a = edges[:-1].copy(); b = edges[1:].copy()
    min_lb = np.inf; viol = []; unres = 0; ncell = 0
    while len(a):
        ncell += len(a)
        if ncell > maxcells: unres += len(a); break
        m = 0.5 * (a + b); h = 0.5 * (b - a)
        f0, f1, e = gvals(P, Q, ts, m)
        g = c + f0.sum(1); gp = f1.sum(1)
        mag = np.abs(c) + np.abs(f0).sum(1) + 1e-300
        err = INFL * (mag + np.abs(gp) * 0 )
        # |g''| bound on [a,b]: sum_j max over endpoints |(-2tQ + t^2 P + t^2 Q s)| * exp(-t a)
        t = ts[None, :]
        la = -2 * t * Q[None, :] + t * t * (P[None, :] + Q[None, :] * a[:, None]); lb_ = -2 * t * Q[None, :] + t * t * (P[None, :] + Q[None, :] * b[:, None])
        G2 = (np.maximum(np.abs(la), np.abs(lb_)) * np.exp(-t * a[:, None])).sum(1) * (1 + 1e-9)
        lb = g - np.abs(gp) * h - 0.5 * G2 * h * h - err - INFL * np.abs(gp) * h
        ok = lb >= thr
        if ok.any(): min_lb = min(min_lb, float(lb[ok].min()))
        bad = ~ok
        if not bad.any(): break
        # violators: midpoint itself below thr
        vv = bad & (g - err < thr)
        for i in np.where(vv)[0]: viol.append((float(m[i]), float(g[i])))
        small = bad & (h < HMIN + 1e-9 * m)
        unres += int(small.sum()); keep = bad & ~small
        am, bm, mm = a[keep], b[keep], m[keep]
        a = np.concatenate([am, mm]); b = np.concatenate([mm, bm])
    return min_lb, viol, unres, ncell


def certify_cell(ctx, cell, max_rounds=MAX_ROUNDS, log=None):
    """returns dict(status in certified|witness|margin|cap|solver, phi, eps, margin, rounds, added)"""
    added = 0
    for rnd in range(max_rounds + 1):
        E = elastic(ctx, cell, ctx.atoms())
        if E["status"] != "ok": return dict(status="solver", rounds=rnd, added=added)
        phi = E["phi"]
        if phi <= PHI_MIN: return dict(status="witness", phi=phi, rounds=rnd, added=added)
        u, v = E["u"], E["v"]; v1, v2 = float(v[0]), float(v[1])
        P, Q, uN, taus = coeffs(ctx, cell, u)
        tau_slack = taus - v1 + v2
        thr = -0.1 * phi
        min_lb, viol, unres, nc = bnb(P, Q, ctx.ts, v1, thr)
        tail = float(((np.abs(P) + np.abs(Q)) * np.exp(-ctx.ts)).sum()) * (1 + 1e-9)
        g0 = v1 + float(P.sum())
        scale = 1e-11 * (1.0 + float(np.abs(P).sum()))
        eps1 = max(0.0, -min_lb if np.isfinite(min_lb) else 0.0, tail - v1, -v1, uN - g0) + scale
        if unres == 0 and viol == []:
            eps2 = eps1 - tau_slack
            margin = phi - eps2
            if margin > 1e-9 * max(1.0, phi) and margin > 0:
                return dict(status="certified", phi=phi, eps1=eps1, margin=margin, rounds=rnd, added=added, P=P, Q=Q, v1=v1, v2=v2, uN=uN, taus=taus, thr=thr, ncell=nc)
            # no violator but eps too large relative to phi
            return dict(status="margin", phi=phi, eps1=eps1, rounds=rnd, added=added)
        if rnd == max_rounds: return dict(status="cap", phi=phi, rounds=rnd, added=added)
        if not viol: return dict(status="cap", phi=phi, rounds=rnd, added=added, note="unresolved")
        viol.sort(key=lambda x: x[1]); pts = []
        for s, gv in viol:
            if all(abs(s - p) > 1e-3 * p for p in pts): pts.append(s)
            if len(pts) >= 8: break
        ctx.extra += pts; added += len(pts)
    return dict(status="cap", rounds=max_rounds, added=added)


# ---------------------------------------------------------------- python-flint ARB ball re-verification of a certified cell (100-bit, rigorous enclosures)
def mp_verify(ts, cert, prec=100):
    from flint import arb, ctx
    ctx.prec = prec
    P, Q, v1, uN = cert["P"], cert["Q"], cert["v1"], cert["uN"]
    c = arb(float(v1)) + arb(float(cert["eps1"]))
    K = len(ts); Pi = [arb(float(x)) for x in P]; Qi = [arb(float(x)) for x in Q]; ti = [arb(float(x)) for x in ts]
    g0 = c
    for p in Pi: g0 = g0 + p
    if not (g0 - arb(float(uN)) >= 0) or not (c >= 0): return False, "atom"
    tail = arb(0)
    for j in range(K): tail = tail + (abs(Pi[j]) + abs(Qi[j])) * (-ti[j]).exp()
    if not (c - tail >= 0): return False, "tail"
    edges = [0.0] + list(np.geomspace(1e-7, 1.0, 90)); cells = list(zip(edges[:-1], edges[1:])); n = 0
    while cells:
        a, b = cells.pop(); n += 1
        if n > 400000: return False, "cap"
        ia, ib = arb(float(a)), arb(float(b)); m = (ia + ib) / 2; h = (ib - ia) / 2
        g = c; gp = arb(0); G2 = arb(0)
        for j in range(K):
            e = (-ti[j] * m).exp(); lin = Pi[j] + Qi[j] * m
            g = g + lin * e; gp = gp + (Qi[j] - ti[j] * lin) * e
            ea = (-ti[j] * ia).exp()
            l1 = -2 * ti[j] * Qi[j] + ti[j] * ti[j] * (Pi[j] + Qi[j] * ia); l2 = -2 * ti[j] * Qi[j] + ti[j] * ti[j] * (Pi[j] + Qi[j] * ib)
            G2 = G2 + (abs(l1) + abs(l2)) * ea
        lb = g - abs(gp) * h - G2 * h * h / 2
        if lb >= 0: continue
        if float(b) - float(a) < 1e-14 * max(1.0, float(b)): return False, "unres"
        mm = 0.5 * (a + b); cells += [(a, mm), (mm, b)]
    e2 = arb(float(cert["eps1"])) - arb(float(cert["taus"] - cert["v1"] + cert["v2"]))
    return bool(arb(float(cert["v2"])) + e2 < 0), "ok"


# ---------------------------------------------------------------- leaves of the r-branch-and-bound (grid, as find_nodrop)
def leaves(S, lo, hi, N, Jw, d0min, d0max):
    RT.DELTA = DELTA
    edges = np.linspace(d0min - 1e-9, d0max + 1e-9, RT.N0 + 1); ints = [(edges[i], edges[i + 1]) for i in range(RT.N0)]
    mid = 0.5 * (d0min + d0max); ints.sort(key=lambda iv_: abs(0.5 * (iv_[0] + iv_[1]) - mid)); dq = collections.deque(ints)
    out = []; nlp = 0
    while dq:
        if nlp >= RT.MAXLP: return None, "cap"
        ra, rb = dq.popleft()
        st, y = RT.feas(S, ra, rb, Jw, ra - DELTA - MARG); nlp += 1
        if st == 2: out.append((ra, rb)); continue
        if st != 0: return None, "badstat"
        ok, dj, md = RT.verify(S, y, lo, hi, N, Jw)
        if ok and md >= -DELTA: return None, "WITNESS"
        st, y = RT.feas(S, ra, rb, Jw, rb - DELTA); nlp += 1
        if st == 0:
            ok, dj, md = RT.verify(S, y, lo, hi, N, Jw)
            if ok and md >= -DELTA - 1e-9: return None, "WITNESS"
        if rb - ra < RT.WMIN: return None, "unres"
        m = 0.5 * (ra + rb); dq.append((ra, m)); dq.append((m, rb))
    return out, "DROP"


def certify_dataset(ts, lo, hi, N, Jw, verify_mp=False):
    t0 = time.time(); ctx = Ctx(ts, lo, hi, N, Jw)
    S = A8.lp_setup(ts, lo, hi, A8.atoms_grid(ctx.tmax, RT.NATOMS), N)
    d0min, d0max, st = A8.lp_one(S, 0)
    if st != 0: return dict(verdict="boxinfeasible", seconds=time.time() - t0)
    lv, why = leaves(S, lo, hi, N, Jw, d0min, d0max)
    if lv is None: return dict(verdict="regen_" + why, seconds=time.time() - t0)
    a0 = d0min - 1e-9; b0 = d0max + 1e-9
    cells = [Cell(ra, rb, ra - DELTA - MARG, "leaf") for ra, rb in lv]
    cells += [Cell(a0 - PAD, a0, a0 - PAD - DELTA - MARG, "pad_lo"), Cell(None, a0 - PAD, None, "end_lo"),
              Cell(b0, b0 + PAD, b0 - DELTA - MARG, "pad_hi"), Cell(b0 + PAD, None, b0 + PAD - DELTA - MARG, "end_hi")]
    # note: pad_lo / end_lo have ra / rb so that  d0 in [a0-PAD, a0]  and  d0 <= a0-PAD ; end_lo must have rb = a0-PAD
    cells[-3] = Cell(None, a0 - PAD, None, "end_lo")
    res = []; stat = collections.Counter(); minmargin = np.inf; nadd = 0; mpok = None; first_fail = None
    for cl in cells:
        r = certify_cell(ctx, cl); r["name"] = cl.name; res.append(r); stat[r["status"]] += 1; nadd += r.get("added", 0)
        if r["status"] == "certified": minmargin = min(minmargin, r["margin"] / r["phi"])
        elif first_fail is None: first_fail = cl.name + ":" + r["status"]
    certified = all(r["status"] == "certified" for r in res)
    if verify_mp and certified:
        mpok = all(mp_verify(ts, r)[0] for r in res)
    solver = stat.get("solver", 0) > 0
    verdict = "CERTIFIED" if certified else ("SOLVERFAIL" if solver else "NOTCERT")
    return dict(verdict=verdict, ncells=len(cells), n_cert=stat.get("certified", 0), n_witness=stat.get("witness", 0), n_margin=stat.get("margin", 0),
                n_cap=stat.get("cap", 0), n_solver=stat.get("solver", 0), first_fail=first_fail or "", added=nadd, atoms_total=len(ctx.extra),
                min_margin_frac=minmargin, mp_verified=mpok, seconds=time.time() - t0, d0min=d0min, d0max=d0max)


# ---------------------------------------------------------------- drivers
CELLS = None
def work(task):
    kind, cid, ds, flag = task; c = CELLS[cid]
    base = dict(kind=kind, cell=cid, ds=ds)
    if kind == "exact":
        ts, N = c["ts_w"], c["N"]; Jw = int(c["exact_J"]); q = c["graph"].exact(ts)[0]; eps = flag
        r = certify_dataset(ts, q * (1 - eps), q * (1 + eps), N, Jw); return dict(base, eps=eps, **r)
    b = MD.data_box(c, ds)
    ts, lo, hi, N, Jw = b
    if kind == "cert":
        return dict(base, **certify_dataset(ts, lo, hi, N, Jw, verify_mp=flag))
    if kind == "witcell":                       # sentinel: cell around the stored witness r must NOT certify
        rr = flag; ctx = Ctx(ts, lo, hi, N, Jw); cl = Cell(rr - 1e-3, rr + 1e-3, rr - 1e-3 - DELTA - MARG, "witcell")
        r = certify_cell(ctx, cl); return dict(base, wit_status=r["status"], phi=r.get("phi", np.nan))


def prep_cells():
    global CELLS
    sp_all = [s for s in RT.RI.core_specs() if s["kind"] in ("control", "crossover")]
    want = RT.TARGETS + RT.CONTROLS
    ids = {s["id"] for s in sp_all if (s["fam"], s["L"]) in want}
    sp = RT.specs(ids); CELLS = {s["id"]: s for s in sp}; return sp


def main():
    mode = args[0]; t0 = time.time(); sp = prep_cells(); c2 = [s for s in sp if s["fam"] == "aniso2_w16"][0]
    nd = pd.read_csv(os.path.join(R.RES, "nodrop_data.csv")); nd2 = nd[nd.cell == c2["id"]].sort_values("ds")
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 5; assert ncpu <= 5
    if mode == "unit1":
        ds = int(args[1]); r = work(("cert", c2["id"], ds, False)); print(r); return
    if mode in ("run", "run_mpall"):
        drops = nd2[nd2.g_outcome == "DROP"].ds.tolist(); print(len(drops), "stored DROP", flush=True)
        rng = np.random.default_rng(2211); mpset = set(rng.choice(drops, 10, replace=False).tolist())
        if mode == "run_mpall": mpset = set(drops)   # PREREG_REDTEAM_FINAL.md P8: Arb ball check on all 37
        tasks = [("cert", c2["id"], ds, ds in mpset) for ds in drops]; out = os.path.join(R.RES, "certify_2d.csv" if mode == "run" else "certify_2d_mpall.csv")
    elif mode == "sentinels":
        tasks = [("exact", c2["id"], -1, 1e-3)]
        wit = nd2[nd2.g_outcome == "WITNESS"]
        tasks += [("witcell", c2["id"], int(r.ds), float(r.g_r)) for r in wit.itertuples()]
        out = os.path.join(R.RES, "certify_sentinels.csv")
    elif mode == "ctlsent":
        nc = pd.read_csv(os.path.join(R.RES, "nodrop_controls.csv")); ring = [s for s in sp if s["fam"] == "ring"][0]["id"]
        w = nc[(nc.cell == ring) & (nc.g_outcome == "WITNESS")].sort_values("ds").head(10)
        tasks = [("witcell", ring, int(r.ds), float(r.g_r)) for r in w.itertuples()]; out = os.path.join(R.RES, "certify_sentinels_ctl.csv")
    rows = []
    with mp.get_context("fork").Pool(min(ncpu, len(tasks))) as pool:
        for k, r in enumerate(pool.imap_unordered(work, tasks)):
            rows.append(r); print(k + 1, len(tasks), {kk: r[kk] for kk in r if kk in ("kind", "ds", "verdict", "ncells", "n_cert", "first_fail", "min_margin_frac", "mp_verified", "seconds", "wit_status")}, flush=True)
            pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False); print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
