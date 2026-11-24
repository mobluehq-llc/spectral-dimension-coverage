#!/usr/bin/env python3
"""Round-3 response, item 2: exact-rational verification of the Amendment 10b outer-relaxation (relaxation B) LPs.  Usage:
  round3_lp_verify.py select                         write results_round3/lp_verify_selection.json (rule + seed below), BEFORE any LP is solved
  round3_lp_verify.py run --out DIR [--ncpu K] [--n N]   verify the selected datasets (regenerates each box from its seed, exactly as run_p22j.box_for)
  round3_lp_verify.py tolerance --out DIR            decision margins of every stored E7j decision (needs only raw_amend10b.csv)

What is verified, and how.  For each dataset and each endpoint (window time j, side upper/lower) the program rebuilds the two LPs that
p22amend10.bound_B solves (variant 1: denominator normalised; variant 2: denominator replaced by the box bound), solves each with HiGHS
(dual simplex, presolve off, so a basis is returned) AND with HiGHS interior point (a second algorithm), and then checks the simplex
basis in EXACT rational arithmetic (python-flint fmpq, every double read as the exact rational it is):
  * primal feasibility of the basic solution (all basic structurals >= 0, all row activities inside their bounds),
  * dual feasibility (all reduced costs of the right sign, equality rows free),
  hence the vertex is an exact optimum of the stated LP data with a rational optimal value. No solver tolerance enters the verdict;
  the verdict is 'EXACT' or 'FAILED' (with the exact size of the violation).
Infeasibility (relaxed set empty) is certified by solving the always-feasible phase-1 LP  min sum(s) s.t. rows with artificial s>=0  and
verifying in the same exact way that its optimum is strictly positive.
The LP data are the double-precision constants of p22amend10.cell_constants (already rounded OUTWARD from python-flint arb balls), taken
as exact. So the claim certified is: the reported endpoint is the exact optimum of the relaxed problem whose data are those doubles."""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
args = sys.argv[1:]; MODE = args[0] if args else ""
OUTD = os.path.join(HERE, "results_round3")


def sel_rule():
    import pandas as pd
    b = pd.read_csv(os.path.join(HERE, "results_amend10_11", "raw_amend10b.csv"), low_memory=False)
    rs = np.random.default_rng(3101)
    sel = {}
    def add(rows, tag):
        for c, d in rows: sel.setdefault((int(c), int(d)), []).append(tag)
    prob = b[(b.B_status != 0) | (b.A_status != 0) | (b.o_status != 0)]
    add(prob[["cell", "ds"]].values, "problem(B/A/orig status != 0)")
    for cid in (702, 26, 703):
        x = b[b.cell == cid][["cell", "ds"]].values; add(x[np.sort(rs.choice(len(x), 8, replace=False))], f"minL3d cell {cid}")
    sens = b[(b.group.str.startswith("sens")) & (~b.cell.isin([702, 26, 703]))][["cell", "ds"]].values; add(sens[np.sort(rs.choice(len(sens), 8, replace=False))], "sens sample")
    st = b[(b.group.str.startswith("stable300")) & (~b.cell.isin([702, 26, 703]))][["cell", "ds"]].values; add(st[np.sort(rs.choice(len(st), 4, replace=False))], "stable sample")
    add(b[b.group.str.startswith("withheld")][["cell", "ds"]].values, "withheld E7j dataset")
    return sorted(sel.items())


if MODE == "select":
    s = sel_rule(); os.makedirs(OUTD, exist_ok=True)
    json.dump([dict(cell=c, ds=d, why=w) for (c, d), w in s], open(os.path.join(OUTD, "lp_verify_selection.json"), "w"), indent=1)
    print(len(s), "datasets selected"); sys.exit(0)

import pandas as pd
from fractions import Fraction


def fq(x):
    from flint import fmpq
    p, q = float(x).as_integer_ratio(); return fmpq(p, q)


def exact_check(c, Aub, bub, Aeq, beq, max_pivots=60, tight=False):
    """solve  min c.x  s.t. Aub x <= bub, Aeq x = beq, x >= 0  with HiGHS (dual simplex, presolve off) and certify the optimum EXACTLY:
    take the solver's basis, solve the basis system in rational arithmetic, test primal and dual feasibility exactly; if the solver's basis
    is exactly primal feasible but dual infeasible by a solver-tolerance amount, continue with exact rational primal-simplex pivots (Bland)
    until optimal.  Returns dict with exact_value (float of the exact rational optimum), exact flag, pivots, and the second-algorithm value."""
    import highspy
    from scipy.sparse import csc_matrix, vstack
    from flint import fmpq_mat, fmpq
    m1, m2 = Aub.shape[0], Aeq.shape[0]; n = len(c); m = m1 + m2
    A = vstack([csc_matrix(Aub), csc_matrix(Aeq)]).tocsc()
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("presolve", "off"); h.setOptionValue("solver", "simplex")
    if tight:
        h.setOptionValue("primal_feasibility_tolerance", 1e-10); h.setOptionValue("dual_feasibility_tolerance", 1e-10); h.setOptionValue("simplex_scale_strategy", 0)
    lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = m
    lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.zeros(n); lp.col_upper_ = np.full(n, highspy.kHighsInf)
    rl = np.concatenate([np.full(m1, -highspy.kHighsInf), np.asarray(beq, float)]); ru = np.concatenate([np.asarray(bub, float), np.asarray(beq, float)])
    lp.row_lower_ = rl; lp.row_upper_ = ru
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise; lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    h.passModel(lp); h.run()
    ms = h.getModelStatus(); out = dict(model_status=str(ms), tight=tight)
    if ms != highspy.HighsModelStatus.kOptimal:
        if not tight: return exact_check(c, Aub, bub, Aeq, beq, max_pivots=1000, tight=True)
        return out
    bas = h.getBasis(); out["solver_value"] = float(h.getInfo().objective_function_value)
    cs = list(bas.col_status); rs_ = list(bas.row_status); B = highspy.HighsBasisStatus
    F = lambda v: Fraction(*float(v).as_integer_ratio())
    Ad = A.toarray(); cols = [{r: F(Ad[r, i]) for r in np.nonzero(Ad[:, i])[0]} for i in range(n)]
    cF = [F(v) for v in c]
    lo_r = [None if not np.isfinite(rl[i]) else F(rl[i]) for i in range(m)]; hi_r = [None if not np.isfinite(ru[i]) else F(ru[i]) for i in range(m)]
    basic = [("x", i) for i in range(n) if cs[i] == B.kBasic] + [("r", i) for i in range(m) if rs_[i] == B.kBasic]
    if len(basic) != m: out["exact"] = "FAILED: basis size"; return out
    at = {}                                              # nonbasic row logicals: 'U' / 'L'
    for i in range(m):
        if rs_[i] != B.kBasic: at[i] = "U" if (rs_[i] == B.kUpper or (hi_r[i] is not None and lo_r[i] == hi_r[i])) else "L"
    def colvec(kind, idx):
        v = [Fraction(0)] * m
        if kind == "x":
            for r, a in cols[idx].items(): v[r] = a
        else: v[idx] = Fraction(-1)
        return v
    def to_mat(vs):
        M = fmpq_mat(m, len(vs))
        for cj, v in enumerate(vs):
            for r in range(m):
                if v[r] != 0: M[r, cj] = fmpq(v[r].numerator, v[r].denominator)
        return M
    def fr(q): return Fraction(int(q.p), int(q.q))
    pivots = 0
    while True:
        Bm = to_mat([colvec(k, i) for k, i in basic])
        rhs = fmpq_mat(m, 1)
        for i, st in at.items():
            val = hi_r[i] if st == "U" else lo_r[i]
            if val is not None: rhs[i, 0] = fmpq(val.numerator, val.denominator)
        try: z = Bm.solve(rhs); cB = fmpq_mat(m, 1)
        except Exception as e: out["exact"] = "FAILED: singular basis"; return out
        zb = [fr(z[k, 0]) for k in range(m)]
        for k, (kind, idx) in enumerate(basic): cB[k, 0] = fmpq(cF[idx].numerator, cF[idx].denominator) if kind == "x" else fmpq(0)
        y = [fr(v) for v in [Bm.transpose().solve(cB)[i, 0] for i in range(m)]]
        pv = Fraction(0)
        for k, (kind, idx) in enumerate(basic):
            lo, hi = (Fraction(0), None) if kind == "x" else (lo_r[idx], hi_r[idx])
            if lo is not None and zb[k] < lo: pv = max(pv, lo - zb[k])
            if hi is not None and zb[k] > hi: pv = max(pv, zb[k] - hi)
        if pv > 0:
            if not tight: return exact_check(c, Aub, bub, Aeq, beq, max_pivots=1000, tight=True)
            out.update(exact="FAILED: basis not exactly primal feasible", exact_primal_violation=float(pv), pivots=pivots); return out
        # reduced costs
        bset = {(k, i) for k, i in basic}; cand = None; dviol = Fraction(0)
        for i in range(n):
            if ("x", i) in bset: continue
            d = cF[i] - sum((a * y[r] for r, a in cols[i].items()), Fraction(0))
            if d < 0:
                dviol = max(dviol, -d)
                if cand is None: cand = ("x", i, +1)
        for i in range(m):
            if ("r", i) in bset or lo_r[i] == hi_r[i]: continue
            d = y[i]
            if at[i] == "U" and d > 0:
                dviol = max(dviol, d)
                if cand is None: cand = ("r", i, -1)
            if at[i] == "L" and d < 0:
                dviol = max(dviol, -d)
                if cand is None: cand = ("r", i, +1)
        if cand is None: break
        if pivots == 0: out["solver_basis_dual_violation"] = float(dviol)
        if pivots >= max_pivots:
            if not tight: return exact_check(c, Aub, bub, Aeq, beq, max_pivots=1000, tight=True)
            out.update(exact="FAILED: pivot limit", exact_dual_violation=float(dviol), pivots=pivots); return out
        kind, q, s = cand; aq = colvec(kind, q)
        wv = Bm.solve(to_mat([aq])); w = [fr(wv[k, 0]) for k in range(m)]
        best = None
        for k, (bk, bi) in enumerate(basic):
            dk = -s * w[k]
            if dk == 0: continue
            lo, hi = (Fraction(0), None) if bk == "x" else (lo_r[bi], hi_r[bi])
            if dk < 0 and lo is not None: th, to = (zb[k] - lo) / (-dk), "L"
            elif dk > 0 and hi is not None: th, to = (hi - zb[k]) / dk, "U"
            else: continue
            if best is None or th < best[0]: best = (th, k, to)
        if best is None: out["exact"] = "FAILED: unbounded"; return out
        th, k, to = best; lk = basic[k]
        basic[k] = (kind, q)
        if kind == "r": at.pop(q, None)
        if lk[0] == "r": at[lk[1]] = to
        pivots += 1
    val = sum((cF[i] * zb[k] for k, (kind, i) in enumerate(basic) if kind == "x"), Fraction(0))
    out.update(exact_value=float(val), exact_primal_violation=0.0, exact_dual_violation=0.0, exact="EXACT", pivots=pivots)
    h2 = highspy.Highs(); h2.setOptionValue("output_flag", False); h2.setOptionValue("solver", "ipm"); h2.passModel(lp); h2.run()
    out["ipm_value"] = float(h2.getInfo().objective_function_value) if h2.getModelStatus() == highspy.HighsModelStatus.kOptimal else None
    return out


def phase1_check(Aub, bub, Aeq, beq):
    """certify emptiness: min sum of artificials; verify exactly that the optimum is > 0."""
    m1, m2 = Aub.shape[0], Aeq.shape[0]; n = Aub.shape[1]
    # Aub x - s1 <= bub ; Aeq x + sp - sm = beq ; all >= 0 ; min sum s
    Ab = np.hstack([Aub, -np.eye(m1), np.zeros((m1, 2 * m2))]); Ae = np.hstack([Aeq, np.zeros((m2, m1)), np.eye(m2), -np.eye(m2)])
    c = np.concatenate([np.zeros(n), np.ones(m1 + 2 * m2)])
    r = exact_check(c, Ab, bub, Ae, beq)
    r["phase1_positive"] = bool(r.get("exact") == "EXACT" and r.get("exact_value", 0) > 0)
    return r


def lp_variant(C, lo, hi, N, j, side, variant, rows):
    import p22amend10 as A10
    nc = C["nc"]; t = C["ts"][j]; Aub, bub, Aeq0, beq0, ref = rows
    if variant == 1:
        den = (C["Am"][j] if side > 0 else C["Ap"][j]) / ref[j]
        Aeq = np.vstack([Aeq0, np.append(den, 0.0)[None, :]]); cvec = np.zeros(nc + 1); cvec[:nc] = 2.0 * t * (C["Gp"][j] if side > 0 else C["Gm"][j]) / ref[j]
    else:
        bd = lo[j] if side > 0 else hi[j]
        if not (bd > 0): return None
        Aeq = np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]]); cvec = np.zeros(nc + 1); cvec[:nc] = 2.0 * t * (C["Gp"][j] if side > 0 else C["Gm"][j]) / bd
    return -side * cvec, Aub, bub, Aeq, np.array([0.0, 1.0])


def verify_dataset(cid, ds):
    import run_p22j as RJ, p22amend10 as A10
    c = RJ.get_cell(cid); ts, N = c["ts_w"], c["N"]
    Jw, lo, hi = RJ.box_for(c, ds); idx = list(range(Jw))
    C = A10.cell_constants(ts, A10.NCERT); rows = A10._relaxed_rows(C, lo, hi, N)
    t0 = time.time(); bl, bh, st = A10.band_B(ts, lo, hi, N, idx)            # the reported certified band (recomputed here)
    rec = dict(cell=cid, ds=ds, J=Jw, B_status=st, n_lp=0, n_exact=0, n_failed=0, max_ds_ipm_gap=0.0, max_rel_gap_exact_vs_solver=0.0, max_primal_viol=0.0, max_dual_viol=0.0, endpoints=[])
    if st != 0:
        # certify emptiness with the exact phase-1 check on the variant-2 LP rows (rows alone define the set)
        Aub, bub, Aeq0, beq0, ref = rows; nc = C["nc"]
        Aeq = np.vstack([Aeq0, np.append(np.zeros(nc), 1.0)[None, :]])
        r = phase1_check(Aub, bub, Aeq, np.array([0.0, 1.0])); rec["empty_set_certified_exactly"] = bool(r["phase1_positive"]); rec["phase1"] = r
        return rec
    for jj in idx:
        for side in (+1, -1):
            vals = []; nvar = 0
            for variant in (1, 2):
                P = lp_variant(C, lo, hi, N, jj, side, variant, rows)
                if P is None: continue
                nvar += 1; r = exact_check(*P); rec["n_lp"] += 1
                if r.get("exact") == "EXACT":
                    rec["n_exact"] += 1; rec["n_pivoted"] = rec.get("n_pivoted", 0) + int(r.get("pivots", 0) > 0); rec["max_solver_basis_dual_viol"] = max(rec.get("max_solver_basis_dual_viol", 0.0), r.get("solver_basis_dual_violation", 0.0)); v = -side * r["exact_value"]
                    rec["max_primal_viol"] = max(rec["max_primal_viol"], r["exact_primal_violation"]); rec["max_dual_viol"] = max(rec["max_dual_viol"], r["exact_dual_violation"])
                    if r.get("ipm_value") is not None: rec["max_ds_ipm_gap"] = max(rec["max_ds_ipm_gap"], abs(r["ipm_value"] - r["exact_value"]))
                    vals.append(v)
                else:
                    rec["n_failed"] += 1; rec["endpoints"].append(dict(j=jj, side=side, variant=variant, status=r.get("exact", r.get("model_status")), pv=r.get("exact_primal_violation"), dv=r.get("exact_dual_violation")))
            if vals and len(vals) == nvar:
                ex = (min(vals) if side > 0 else max(vals)); ex = max(ex, 0.0) if side < 0 else ex
                rep = bh[jj] if side > 0 else bl[jj]
                rel = abs(ex - rep) / max(abs(rep), 1e-12); rec["max_rel_gap_exact_vs_solver"] = max(rec["max_rel_gap_exact_vs_solver"], rel)
    rec["seconds"] = time.time() - t0
    return rec


def _w(a): return verify_dataset(*a)


def run():
    import multiprocessing as mp
    out = args[args.index("--out") + 1]; os.makedirs(out, exist_ok=True)
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 4; assert ncpu <= 12
    sel = json.load(open(os.path.join(HERE, "results_round3", "lp_verify_selection.json")))
    if "--n" in args: sel = sel[: int(args[args.index("--n") + 1])]
    import run_p22j as RJ
    for cid in sorted({s["cell"] for s in sel}): RJ.get_cell(cid)
    tasks = [(s["cell"], s["ds"]) for s in sel]; why = {(s["cell"], s["ds"]): s["why"] for s in sel}
    res = []
    with mp.get_context("fork").Pool(ncpu) as pool:
        for k, r in enumerate(pool.imap_unordered(_w, tasks)):
            r["why"] = " | ".join(why[(r["cell"], r["ds"])]) if isinstance(why[(r["cell"], r["ds"])], list) else why[(r["cell"], r["ds"])]
            res.append(r); print(k + 1, "/", len(tasks), r["cell"], r["ds"], r.get("n_exact"), "/", r.get("n_lp"), flush=True)
            json.dump(res, open(os.path.join(out, "lp_verify_results.json"), "w"), indent=1, default=float)
    print("done")


def tolerance():
    out = args[args.index("--out") + 1]; os.makedirs(out, exist_ok=True)
    b = pd.read_csv(os.path.join(HERE, "results_amend10_11", "raw_amend10b.csv"), low_memory=False)
    ok = b[b.B_status == 0].copy()
    rows = []
    for col, nm in (("Bj_minU_drop", "drop"), ("Bj_minU_rise", "rise")):
        u = ok[col]; has = u.notna()
        # a decision 'declare' means U < -0.25. margin = distance of U from -0.25 (negative U = deeper = safer)
        marg = (u + 0.25).abs()
        for tol in (1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
            rows.append(dict(decision=nm, tol=tol, datasets_with_finite_U=int(has.sum()), within_tol_of_threshold=int((has & (marg <= tol)).sum()),
                             of_which_declared=int((has & (marg <= tol) & (u < -0.25)).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(out, "lp_tolerance_margins.csv"), index=False)
    # solver-trouble counters and ambiguity
    s = dict(datasets=len(b), B_ok=len(ok), Bj_nbad_gt0=int((ok.Bj_nbad.fillna(0) > 0).sum()), Bj_amb_gt0=int((ok.Bj_amb.fillna(0) > 0).sum()),
             Bj_decl_drop=int((ok.Bj_neg.fillna(0) == 1).sum()), Bj_decl_rise=int((ok.Bj_pos.fillna(0) == 1).sum()), total_lp=int(ok.Bj_nlp.fillna(0).sum()))
    json.dump(s, open(os.path.join(out, "lp_tolerance_counters.json"), "w"), indent=1)
    print(pd.DataFrame(rows).to_string()); print(s)


def retry():
    """re-verify only the datasets of a previous run that had a failed endpoint, with the final exact_check (tight retry inside)."""
    import multiprocessing as mp
    src = args[args.index("--in") + 1]; out = args[args.index("--out") + 1]; ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 4; assert ncpu <= 12
    prev = json.load(open(os.path.join(src, "lp_verify_results.json")))
    todo = [(r["cell"], r["ds"]) for r in prev if r.get("n_failed", 0) > 0 or r.get("n_exact", 0) + r.get("n_failed", 0) != r.get("n_lp", 0)]
    import run_p22j as RJ
    for cid in sorted({t[0] for t in todo}): RJ.get_cell(cid)
    new = {}
    with mp.get_context("fork").Pool(ncpu) as pool:
        for r in pool.imap_unordered(_w, todo): new[(r["cell"], r["ds"])] = r; print("retried", r["cell"], r["ds"], r["n_exact"], "/", r["n_lp"], flush=True)
    merged = []
    for r in prev:
        k = (r["cell"], r["ds"])
        if k in new: n = new[k]; n["why"] = r.get("why"); n["retried"] = True; merged.append(n)
        else: merged.append(r)
    json.dump(merged, open(os.path.join(out, "lp_verify_results_final.json"), "w"), indent=1, default=float); print("retry done", len(new))


if MODE == "run": run()
elif MODE == "retry": retry()
elif MODE == "tolerance": tolerance()
