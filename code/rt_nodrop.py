#!/usr/bin/env python3
"""GLOBAL NO-DROP LP (PREREG_GLOBAL_NODROP.md). Frozen files are imported, never edited.
Usage: rt_nodrop.py exact | data [--ndata N] [--ncpu K] [--cells id,id] | unit
Question per box (lo, hi) on the grid ts, window length Jw: does some admissible measure (lazy walk: q(t)=sum w_k exp(-t s_k), w>=0, sum w = 1,
atom at s=0 of mass >= 1/N, lo <= q(ts) <= hi) have NO drop > 0.25, i.e.  d(t_j) >= d(t_0) - 0.25 for every j in 1..Jw-1 ?
  WITNESS      : yes, a verified measure is returned (the global null is NOT excluded)
  DROP         : no, for every r = d(t_0) in the E7 band of t_0 every relaxed LP is infeasible with margin MARG (certified drop), gated at 2x atoms
  INCONCLUSIVE : neither inside the budget (width WMIN, LP cap MAXLP, numerical status)
r-search: rigorous branch and bound on r-intervals [ra, rb] (d0 in [ra, rb]):
  OUTER (relaxation, superset): d_j >= ra - 0.25 - MARG for all j>=1.  Infeasible => no no-drop measure has d0 in [ra, rb] (certified exclusion).
  INNER (restriction, subset):  d_j >= rb - 0.25 for all j>=1.        Feasible and the returned measure verified => WITNESS.
Dual witness (drop witness): a verified measure in the box with d_j <= d_0 - 0.25 for some j (restriction LP over r-intervals)."""
import os, sys, time, json
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np, pandas as pd
import multiprocessing as mp
from scipy.optimize import linprog
import p22lib as L, p22amend8 as A8, p22amend9 as A9
import run_p22 as R, run_p22i as RI
sys.argv = _argv
args = sys.argv[1:]
DELTA = L.DELTA
MARG = 0.005          # certification margin: exclusion must hold even if drops up to DELTA+MARG were allowed
WMIN = 2e-4           # smallest r-interval
MAXLP = 400           # LP cap per dataset and atom grid
N0 = 6                # initial r-intervals
NATOMS = 799          # atom grid (gate: 1597)
VTOL = 1e-7           # verification tolerance (relative box, simplex)
EXACT_EPS = 1e-3      # relative half-width of the "exact q" box (ladder 1e-4, 1e-3, 1e-2 in exact mode; the atom grid cannot meet a narrower box)


def _solve(c, Aub, bub, Aeq, beq):
    r = None
    for meth in ("highs", "highs-ipm", "highs-ds"):
        r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method=meth)
        if r.status in (0, 2): break
    return r


def feas(S, ra, rb, Jw, thr_j, sense="nodrop"):
    """feasibility LP, tau fixed to 1 (so y = w). d0 in [ra, rb].  sense 'nodrop': d_j >= thr_j for all 1<=j<Jw. sense 'drop' : d_j <= thr_j for ONE j (S['jsel'])."""
    n = S["n"]; As = S["As"]; sf = S["sfin"]; ts = S["ts"]
    g = lambda j: 2.0 * ts[j] * sf * As[j]
    rows = []
    r1 = np.zeros(n + 1); r1[:n] = ra * As[0] - g(0); rows.append(r1)
    r2 = np.zeros(n + 1); r2[:n] = g(0) - rb * As[0]; rows.append(r2)
    if sense == "nodrop":
        for j in range(1, Jw):
            a = np.zeros(n + 1); a[:n] = thr_j * As[j] - g(j); rows.append(a)      # g_j y >= thr * a_j y
    else:
        j = S["jsel"]; a = np.zeros(n + 1); a[:n] = g(j) - thr_j * As[j]; rows.append(a)
    Aub = np.vstack([S["Aub"], np.array(rows)]); bub = np.concatenate([S["bub"], np.zeros(len(rows))])
    Aeq = np.zeros((2, n + 1)); Aeq[0, :n] = 1.0; Aeq[0, n] = -1.0; Aeq[1, n] = 1.0; beq = np.array([0.0, 1.0])
    r = _solve(np.zeros(n + 1), Aub, bub, Aeq, beq)
    if r.status == 2: return 2, None
    if r.status != 0: return int(r.status), None
    return 0, r.x[:n]


def curve(S, y):
    """d(t_j) of the measure y (= w) on the whole grid and the scaled q."""
    q = S["As"] @ y; g = 2.0 * S["ts"][:, None] * S["As"] * S["sfin"][None, :]
    return (g @ y) / q, q


def verify(S, y, lo, hi, N, Jw):
    """independent check of a returned measure: simplex, atom at 0, box (relative VTOL), and the actual drop. returns (ok, d, mindrop)"""
    if y.min() < -VTOL: return False, None, None
    yy = np.maximum(y, 0.0)
    if abs(yy.sum() - 1.0) > 1e-6: return False, None, None
    if yy[0] < 1.0 / N - 1e-9: return False, None, None
    d, qs = curve(S, yy)
    ref = np.maximum(0.5 * (lo + hi), 1e-300); q = qs * ref
    if (q < lo * (1 - VTOL) - 1e-300).any() or (q > hi * (1 + VTOL)).any(): return False, None, None
    dj = d[:Jw]
    return True, dj, float((dj[1:] - dj[0]).min()) if Jw > 1 else 0.0


def find_nodrop(S, lo, hi, N, Jw, d0min, d0max):
    """B&B over r. returns dict(outcome, nlp, unresolved, witness...)."""
    nlp = 0; unres = 0; badstat = 0; excluded = 0
    edges = np.linspace(d0min - 1e-9, d0max + 1e-9, N0 + 1)
    ints = [(edges[i], edges[i + 1]) for i in range(N0)]
    mid = 0.5 * (d0min + d0max); ints.sort(key=lambda iv: abs(0.5 * (iv[0] + iv[1]) - mid))
    from collections import deque
    dq = deque(ints)
    while dq:
        if nlp >= MAXLP:
            unres += len(dq); break
        ra, rb = dq.popleft()
        st, y = feas(S, ra, rb, Jw, ra - DELTA - MARG); nlp += 1          # OUTER
        if st == 2: excluded += 1; continue
        if st != 0: badstat += 1; unres += 1; continue
        ok, dj, md = verify(S, y, lo, hi, N, Jw)                           # the outer solution may itself be a witness
        if ok and md >= -DELTA:
            return dict(outcome="WITNESS", nlp=nlp, unres=unres, bad=badstat, r=float(dj[0]), mindelta=md, d=";".join("%.4f" % v for v in dj), y=y)
        st, y = feas(S, ra, rb, Jw, rb - DELTA); nlp += 1                  # INNER
        if st == 0:
            ok, dj, md = verify(S, y, lo, hi, N, Jw)
            if ok and md >= -DELTA - 1e-9:
                return dict(outcome="WITNESS", nlp=nlp, unres=unres, bad=badstat, r=float(dj[0]), mindelta=md, d=";".join("%.4f" % v for v in dj), y=y)
        if rb - ra < WMIN: unres += 1; continue
        m = 0.5 * (ra + rb); dq.append((ra, m)); dq.append((m, rb))
    out = "DROP" if (unres == 0 and not dq) else "INCONCLUSIVE"
    return dict(outcome=out, nlp=nlp, unres=unres, bad=badstat, r=np.nan, mindelta=np.nan, d="", y=None)


def find_dropwitness(S, lo, hi, N, Jw, d0min, d0max):
    """verified admissible measure with a drop > DELTA from t_0 to some t_j (restriction LP: d0 >= ra, d_j <= ra - DELTA). Late times first."""
    edges = np.linspace(d0min, d0max, N0 + 1); nlp = 0
    for j in range(Jw - 1, 0, -1):
        S["jsel"] = j
        for i in range(N0 - 1, -1, -1):
            ra, rb = edges[i], edges[i + 1]
            st, y = feas(S, ra, rb, Jw, ra - DELTA - 1e-4, sense="drop"); nlp += 1
            if st == 0:
                ok, dj, md = verify(S, y, lo, hi, N, Jw)
                if ok and md < -DELTA: return dict(found=1, nlp=nlp, maxdrop=-md, j=j)
    return dict(found=0, nlp=nlp, maxdrop=np.nan, j=-1)


def global_nodrop(ts, lo, hi, N, Jw, natoms=NATOMS):
    """full test with gate. returns row dict."""
    tmax = float(ts[-1]); t0 = time.time()
    def build(n): return A8.lp_setup(ts, lo, hi, A8.atoms_grid(tmax, n), N)
    S = build(natoms)
    d0min, d0max, st = A8.lp_one(S, 0)
    if st != 0:
        return dict(g_outcome="BOXINFEASIBLE", g_nlp=0, g_gate=0, g_seconds=time.time() - t0)
    res = find_nodrop(S, lo, hi, N, Jw, d0min, d0max); nlp = res["nlp"]; gate = 0
    if res["outcome"] == "DROP":                                           # gate: a certified drop must persist on twice the atoms
        S2 = build(2 * (natoms - 1) + 1)
        a0, a1, st2 = A8.lp_one(S2, 0)
        if st2 != 0: res = dict(res, outcome="INCONCLUSIVE"); gate = 1
        else:
            r2 = find_nodrop(S2, lo, hi, N, Jw, a0, a1); nlp += r2["nlp"]
            if r2["outcome"] != "DROP": res = r2; gate = 1
    row = dict(g_outcome=res["outcome"], g_nlp=nlp, g_unres=res["unres"], g_bad=res["bad"], g_gate=gate, g_r=res["r"], g_mindelta=res["mindelta"], g_d=res["d"],
               g_d0min=d0min, g_d0max=d0max)
    dw = find_dropwitness(S, lo, hi, N, Jw, d0min, d0max)
    row.update(dw_found=dw["found"], dw_maxdrop=dw["maxdrop"], dw_j=dw["j"], dw_nlp=dw["nlp"])
    row["g_seconds"] = time.time() - t0
    return row


# ------------------------------------------------------------------ data
def joint(ts, lo, hi, N, Jw):
    dmin, dmax, st, gi = A8.e7_gated(ts, lo, hi, N, list(range(Jw)))
    r = dict(e7_status=st)
    if st != 0: r.update(j_infeasible=1, j_drop=0, e7_neg=0); return r, gi
    r["e7_neg"] = int(((dmax[1:Jw] - dmin[0]) < -DELTA).any())
    j = A9.e7j_decide(ts, lo, hi, N, Jw, gi["natoms"])
    r.update(j_infeasible=j["j_infeasible"], j_drop=j["drop"], j_amb=j["n_amb"], j_gate=j["gate_flag"], j_maxL_drop=j["maxL_drop"], j_minU_drop=j["minU_drop"])
    return r, gi


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    p = U.mean(0); Wt = U.shape[0] * m; cnt = p * Wt
    Jw = L.window_len(cnt, p, N, thr)
    base = dict(cell=c["id"], fam=c["fam"], L=c["L"], ds=ds, Jw=Jw, mode="data")
    if Jw < 2: return dict(base, g_outcome="NOWINDOW")
    qhat, lo, hi, Ps, binfo = A8.make_box(U, m, rngs[7])
    r = dict(base, c_box=binfo["c_box"], n_cp=binfo["n_cp"])
    r.update(global_nodrop(ts, lo, hi, N, Jw))
    t1 = time.time(); jr, gi = joint(ts, lo, hi, N, Jw); r.update(jr); r["j_seconds"] = time.time() - t1
    return r


def run_exact(c, eps=EXACT_EPS):
    ts, N = c["ts_w"], c["N"]; Jw = int(c["exact_J"])
    q = c["graph"].exact(ts)[0]
    lo = q * (1 - eps); hi = q * (1 + eps)
    base = dict(cell=c["id"], fam=c["fam"], L=c["L"], ds=-1, Jw=Jw, mode="exact", eps=eps, d_start=float(c["d_start"]), d_end=float(c["d_end"]), max_drop_exact=float(c["max_drop"]))
    r = dict(base); r.update(global_nodrop(ts, lo, hi, N, Jw))
    jr, gi = joint(ts, lo, hi, N, Jw); r.update(jr)
    return r


CELLS = None
def worker(task):
    cid, a, b = task; c = CELLS[cid]
    return [run_dataset(c, ds) for ds in range(a, b)]


def specs(want_ids=None, fams=None):
    sp = [s for s in RI.core_specs() if s["kind"] in ("control", "crossover")]
    rr = np.random.default_rng(2200)
    sp = [s for s in sp if (want_ids is None or s["id"] in want_ids)]
    for s in sp: R.prep(s, rr)
    return sp


def pick(sp, fam, Lv):
    return [s for s in sp if s["fam"] == fam and s["L"] == Lv][0]


TARGETS = [("aniso2_w16", 768), ("aniso3_w8", 160)]
CONTROLS = [("ring", 1024), ("torus2", 96), ("torus3", 24)]


def main():
    global CELLS
    mode = args[0]; t0 = time.time()
    sp_all = [s for s in RI.core_specs() if s["kind"] in ("control", "crossover")]
    want = [(f, l) for f, l in (CONTROLS if "--controls" in args else TARGETS + (CONTROLS if mode == "exact" else []))]
    ids = {s["id"] for s in sp_all if (s["fam"], s["L"]) in want}
    sp = specs(ids); CELLS = {s["id"]: s for s in sp}
    print("prep %.0fs" % (time.time() - t0), [(s["id"], s["fam"], s["L"], s["N"], s["wpn"], s["exact_J"], round(float(s["max_drop"]), 4)) for s in sp], flush=True)
    outdir = R.RES
    if mode == "exact":
        rows = []
        for s in sp:
          for eps in (1e-4, 1e-3, 1e-2):
            r = run_exact(s, eps); rows.append(r)
            print(eps, {k: r.get(k) for k in ("fam", "L", "Jw", "g_outcome", "g_r", "g_mindelta", "g_nlp", "g_gate", "dw_found", "dw_maxdrop", "e7_neg", "j_drop", "max_drop_exact")}, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(outdir, "nodrop_exact.csv"), index=False); return
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else 100
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 5
    assert ncpu <= 5
    chunk = 2
    tasks = [(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)]
    rows = []; out = os.path.join(outdir, "nodrop_controls.csv" if "--controls" in args else "nodrop_data.csv")
    with mp.get_context("fork").Pool(ncpu) as pool:
        for k, r in enumerate(pool.imap_unordered(worker, tasks)):
            rows += r
            if k % 5 == 0:
                print(f"{k+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True); pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False); print("done", len(rows), "%.0fs" % (time.time() - t0), flush=True)


def unit():
    """synthetic unit test, no registered cell: a known measure with a known d(t); box = +-1e-6 around its q."""
    ts = np.array([24, 34, 48, 68, 96, 136, 192, 272, 384], float); N = 10000
    def box(ws, ss):
        q = np.array([(ws * np.exp(-t * ss)).sum() for t in ts]); return q * (1 - 2e-3), q * (1 + 2e-3), q
    # (a) two-scale measure: fast part decays, slow part gives a plateau: d falls (a drop) -> DROP expected
    ws = np.array([1e-4, 0.4999, 0.5]); ss = np.array([0.0, 0.0003, 0.05])
    ws = np.array([1e-4, 0.3, 0.7 - 1e-4]); ss = np.array([0.0, 0.00005, 0.08])
    lo, hi, q = box(ws, ss); d = np.array([2 * t * (ws * ss * np.exp(-t * ss)).sum() / (ws * np.exp(-t * ss)).sum() for t in ts])
    r = global_nodrop(ts, lo, hi, N, len(ts)); print("true d", np.round(d, 3), "drop", round(d[0] - d.min(), 3), "->", r["g_outcome"], r.get("dw_found"))
    # (b) single exponential: d(t) = 2 t s increases: no drop, WITNESS expected
    ws = np.array([1e-4, 1 - 1e-4]); ss = np.array([0.0, 0.01]); lo, hi, q = box(ws, ss)
    d = np.array([2 * t * (ws * ss * np.exp(-t * ss)).sum() / (ws * np.exp(-t * ss)).sum() for t in ts])
    r = global_nodrop(ts, lo, hi, N, len(ts)); print("true d", np.round(d, 3), "->", r["g_outcome"], "mindelta", r.get("g_mindelta"))


if __name__ == "__main__":
    unit() if args[0] == "unit" else main()
