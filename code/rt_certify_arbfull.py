#!/usr/bin/env python3
"""P9 (PREREG_REDTEAM_FINAL.md, amendment 11-24 NIGHT): full-chain Arb recheck of the 37 continuum certificates.
Differences from rt_certify.mp_verify (which is untouched): (1) the dual coefficients P_j, Q_j, the right-hand side taus and the
objective c - taus are RECOMPUTED in 100-bit ball arithmetic from the raw multipliers u, the box (lo, hi), N and the cell's r
thresholds, instead of being read as doubles computed in floating point; (2) the strictly negative dual objective c - taus is
tested in ball arithmetic; (3) the cell set is checked to tile the whole real line in r (exact double equality of shared endpoints).
Usage: rt_certify_arbfull.py run [--ncpu K]"""
import os, sys, time, collections
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np, pandas as pd
import multiprocessing as mp
import rt_certify as C
import rt_nodrop as RT
import rt_mindrop as MD
import p22lib as L, p22amend8 as A8
import run_p22 as R
sys.argv = _argv
args = sys.argv[1:]


def certify_cell2(ctx, cell, max_rounds=C.MAX_ROUNDS):
    """Copy of rt_certify.certify_cell that also returns the raw multiplier vector u of the certifying round."""
    added = 0
    for rnd in range(max_rounds + 1):
        E = C.elastic(ctx, cell, ctx.atoms())
        if E["status"] != "ok": return dict(status="solver", rounds=rnd, added=added)
        phi = E["phi"]
        if phi <= C.PHI_MIN: return dict(status="witness", phi=phi, rounds=rnd, added=added)
        u, v = E["u"], E["v"]; v1, v2 = float(v[0]), float(v[1])
        P, Q, uN, taus = C.coeffs(ctx, cell, u)
        tau_slack = taus - v1 + v2
        thr = -0.1 * phi
        min_lb, viol, unres, nc = C.bnb(P, Q, ctx.ts, v1, thr)
        tail = float(((np.abs(P) + np.abs(Q)) * np.exp(-ctx.ts)).sum()) * (1 + 1e-9)
        g0 = v1 + float(P.sum())
        scale = 1e-11 * (1.0 + float(np.abs(P).sum()))
        eps1 = max(0.0, -min_lb if np.isfinite(min_lb) else 0.0, tail - v1, -v1, uN - g0) + scale
        if unres == 0 and viol == []:
            eps2 = eps1 - tau_slack
            margin = phi - eps2
            if margin > 1e-9 * max(1.0, phi) and margin > 0:
                return dict(status="certified", phi=phi, eps1=eps1, margin=margin, rounds=rnd, added=added, u=np.array(u), v1=v1, v2=v2)
            return dict(status="margin", phi=phi, eps1=eps1, rounds=rnd, added=added)
        if rnd == max_rounds: return dict(status="cap", phi=phi, rounds=rnd, added=added)
        if not viol: return dict(status="cap", phi=phi, rounds=rnd, added=added, note="unresolved")
        viol.sort(key=lambda x: x[1]); pts = []
        for s, gv in viol:
            if all(abs(s - p) > 1e-3 * p for p in pts): pts.append(s)
            if len(pts) >= 8: break
        ctx.extra += pts; added += len(pts)
    return dict(status="cap", rounds=max_rounds, added=added)


def arb_full(ts, lo, hi, N, Jw, cell, cert, prec=100):
    """Rigorous (ball) check of one Farkas certificate. Returns (ok, reason, objective upper bound as float).
    Certificate: u >= 0 on the rows [lo(K), hi(K), N, cell rows]; G(s) = c + sum_j (P_j + Q_j s) exp(-t_j s) must be >= 0 for every
    atom s in [0, inf], and  c - taus < 0  (then 0 >= u^T M x = sum_k y_k (G_k - c) + tau*taus >= tau (taus - c) contradicts tau = 1)."""
    from flint import arb, ctx
    ctx.prec = prec
    A = lambda x: arb(float(x))
    K = len(ts); u = cert["u"]
    if not bool(np.all(u >= 0)): return False, "negative_u", np.nan
    ref = np.maximum(0.5 * (lo + hi), 1e-300)
    P = [(-A(u[j]) + A(u[K + j])) / A(ref[j]) for j in range(K)]; Q = [arb(0) for _ in range(K)]
    uN = A(u[2 * K])
    taus = arb(0)
    for j in range(K): taus = taus + A(u[j]) * A(lo[j]) / A(ref[j]) - A(u[K + j]) * A(hi[j]) / A(ref[j])
    taus = taus + uN / A(N)
    k = 2 * K + 1; t0 = A(ts[0])
    if cell.ra is not None:
        w = A(u[k]); k += 1; P[0] = P[0] + w * A(cell.ra) / A(ref[0]); Q[0] = Q[0] - 2 * t0 * w / A(ref[0])
    if cell.rb is not None:
        w = A(u[k]); k += 1; P[0] = P[0] - w * A(cell.rb) / A(ref[0]); Q[0] = Q[0] + 2 * t0 * w / A(ref[0])
    if cell.thr is not None:
        for j in range(1, Jw):
            w = A(u[k]); k += 1; P[j] = P[j] + w * A(cell.thr) / A(ref[j]); Q[j] = Q[j] - 2 * A(ts[j]) * w / A(ref[j])
    if k != len(u): return False, "row_count", np.nan
    c = A(cert["v1"]) + A(cert["eps1"])
    obj = c - taus
    objup = float(obj.mid()) + float(obj.rad())
    if not (obj < 0): return False, "objective", objup
    ti = [A(x) for x in ts]
    g0 = c
    for p in P: g0 = g0 + p
    if not (g0 - uN >= 0) or not (c >= 0): return False, "atom", objup
    tail = arb(0)
    for j in range(K): tail = tail + (abs(P[j]) + abs(Q[j])) * (-ti[j]).exp()
    if not (c - tail >= 0): return False, "tail", objup
    edges = [0.0] + list(np.geomspace(1e-7, 1.0, 90)); cells = list(zip(edges[:-1], edges[1:])); n = 0
    while cells:
        a, b = cells.pop(); n += 1
        if n > 400000: return False, "cap", objup
        ia, ib = arb(float(a)), arb(float(b)); m = (ia + ib) / 2; h = (ib - ia) / 2
        g = c; gp = arb(0); G2 = arb(0)
        for j in range(K):
            e = (-ti[j] * m).exp(); lin = P[j] + Q[j] * m
            g = g + lin * e; gp = gp + (Q[j] - ti[j] * lin) * e
            ea = (-ti[j] * ia).exp()
            l1 = -2 * ti[j] * Q[j] + ti[j] * ti[j] * (P[j] + Q[j] * ia); l2 = -2 * ti[j] * Q[j] + ti[j] * ti[j] * (P[j] + Q[j] * ib)
            G2 = G2 + (abs(l1) + abs(l2)) * ea
        lb = g - abs(gp) * h - G2 * h * h / 2
        if lb >= 0: continue
        if float(b) - float(a) < 1e-14 * max(1.0, float(b)): return False, "unres", objup
        mm = 0.5 * (a + b); cells += [(a, mm), (mm, b)]
    return True, "ok", objup


def cover_check(cells, a0, b0):
    """cells in certify_dataset order: leaves..., pad_lo, end_lo, pad_hi, end_hi. True iff their r-intervals tile (-inf, inf)."""
    leaves = sorted([c for c in cells if c.name == "leaf"], key=lambda c: c.ra)
    pad_lo, end_lo, pad_hi, end_hi = cells[-4:]
    assert [pad_lo.name, end_lo.name, pad_hi.name, end_hi.name] == ["pad_lo", "end_lo", "pad_hi", "end_hi"]
    ok = (end_lo.ra is None and end_lo.rb == pad_lo.ra and pad_lo.rb == a0 and leaves[0].ra == a0)
    ok = ok and all(leaves[i].rb == leaves[i + 1].ra for i in range(len(leaves) - 1))
    ok = ok and (leaves[-1].rb == b0 and pad_hi.ra == b0 and pad_hi.rb == end_hi.ra and end_hi.rb is None)
    return bool(ok)


def certify_dataset2(ts, lo, hi, N, Jw):
    t0 = time.time(); ctx = C.Ctx(ts, lo, hi, N, Jw)
    S = A8.lp_setup(ts, lo, hi, A8.atoms_grid(ctx.tmax, RT.NATOMS), N)
    d0min, d0max, st = A8.lp_one(S, 0)
    if st != 0: return dict(verdict="boxinfeasible")
    lv, why = C.leaves(S, lo, hi, N, Jw, d0min, d0max)
    if lv is None: return dict(verdict="regen_" + why)
    a0 = d0min - 1e-9; b0 = d0max + 1e-9
    cells = [C.Cell(ra, rb, ra - C.DELTA - C.MARG, "leaf") for ra, rb in lv]
    cells += [C.Cell(a0 - C.PAD, a0, a0 - C.PAD - C.DELTA - C.MARG, "pad_lo"), C.Cell(None, a0 - C.PAD, None, "end_lo"),
              C.Cell(b0, b0 + C.PAD, b0 - C.DELTA - C.MARG, "pad_hi"), C.Cell(b0 + C.PAD, None, b0 + C.PAD - C.DELTA - C.MARG, "end_hi")]
    cover = cover_check(cells, a0, b0)
    stat = collections.Counter(); nadd = 0; nfull = 0; reasons = collections.Counter(); maxobj = -np.inf; minrel = np.inf
    for cl in cells:
        r = certify_cell2(ctx, cl); stat[r["status"]] += 1; nadd += r.get("added", 0)
        if r["status"] == "certified":
            ok, why2, objup = arb_full(ts, lo, hi, N, Jw, cl, r)
            nfull += int(ok); reasons[why2] += 1
            if np.isfinite(objup): maxobj = max(maxobj, objup); minrel = min(minrel, -objup / r["phi"])
    certified = stat.get("certified", 0) == len(cells)
    return dict(verdict="CERTIFIED" if certified else "NOTCERT", ncells=len(cells), n_cert=stat.get("certified", 0), added=nadd,
                n_arbfull=nfull, arb_reasons=dict(reasons), cover_tiles_R=cover, max_obj_upper=maxobj, min_negobj_over_phi=minrel,
                seconds=time.time() - t0)


def work(task):
    cid, ds = task; c = C.CELLS[cid]; ts, lo, hi, N, Jw = MD.data_box(c, ds)
    return dict(cell=cid, ds=ds, **certify_dataset2(ts, lo, hi, N, Jw))


def main():
    sp = C.prep_cells(); c2 = [s for s in sp if s["fam"] == "aniso2_w16"][0]
    nd = pd.read_csv(os.path.join(R.RES, "nodrop_data.csv")); nd2 = nd[nd.cell == c2["id"]].sort_values("ds")
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 4; assert ncpu <= 4
    drops = nd2[nd2.g_outcome == "DROP"].ds.tolist(); print(len(drops), "stored DROP", flush=True)
    tasks = [(c2["id"], ds) for ds in drops]; out = os.path.join(R.RES, "certify_2d_arbfull.csv"); rows = []; t0 = time.time()
    with mp.get_context("fork").Pool(min(ncpu, len(tasks))) as pool:
        for k, r in enumerate(pool.imap_unordered(work, tasks)):
            rows.append(r); print(k + 1, len(tasks), r, flush=True)
            pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False); print("done %.0fs" % (time.time() - t0))


def falsify():
    """Falsifiers: a damaged certificate must FAIL arb_full; a gapped cell set must FAIL cover_check. Writes results/certify_arbfull_falsifier.log."""
    sp = C.prep_cells(); c2 = [s for s in sp if s["fam"] == "aniso2_w16"][0]
    nd = pd.read_csv(os.path.join(R.RES, "nodrop_data.csv")); ds = int(nd[(nd.cell == c2["id"]) & (nd.g_outcome == "DROP")].ds.iloc[0])
    ts, lo, hi, N, Jw = MD.data_box(C.CELLS[c2["id"]], ds); ctx = C.Ctx(ts, lo, hi, N, Jw)
    cl = C.Cell(None, 1.0e9, None, "huge_r")             # r <= 1e9: not a certifiable exclusion; use a real leaf instead
    S = A8.lp_setup(ts, lo, hi, A8.atoms_grid(ctx.tmax, RT.NATOMS), N); d0min, d0max, _ = A8.lp_one(S, 0)
    lv, why = C.leaves(S, lo, hi, N, Jw, d0min, d0max); ra, rb = lv[0]; cl = C.Cell(ra, rb, ra - C.DELTA - C.MARG, "leaf")
    r = certify_cell2(ctx, cl); assert r["status"] == "certified", r["status"]
    out = []; out.append(("intact", arb_full(ts, lo, hi, N, Jw, cl, r)))
    bad = dict(r); bad["eps1"] = r["eps1"] + 1e3; out.append(("eps1 + 1000 (objective must fail)", arb_full(ts, lo, hi, N, Jw, cl, bad)))
    bad = dict(r); bad["v1"] = r["v1"] - 1e3; out.append(("v1 - 1000 (g<0, must fail)", arb_full(ts, lo, hi, N, Jw, cl, bad)))
    bad = dict(r); u2 = r["u"].copy(); u2[0] += 1e3; bad["u"] = u2; out.append(("u[0] + 1000 (g < 0 somewhere, must fail)", arb_full(ts, lo, hi, N, Jw, cl, bad)))
    bad = dict(r); u2 = r["u"].copy(); u2[1] = -1e-3; bad["u"] = u2; out.append(("u[1] = -0.001 (negative multiplier, must fail)", arb_full(ts, lo, hi, N, Jw, cl, bad)))
    bad = dict(r); u2 = r["u"].copy(); u2[0] += 1e-3; bad["u"] = u2; out.append(("u[0] + 0.001 (any u >= 0 is self-consistent here; may pass, reported only)", arb_full(ts, lo, hi, N, Jw, cl, bad)))
    cells = [C.Cell(ra, rb, 0, "leaf") for ra, rb in lv] + [C.Cell(d0min - 1e-9 - C.PAD, d0min - 1e-9, 0, "pad_lo"), C.Cell(None, d0min - 1e-9 - C.PAD, None, "end_lo"),
             C.Cell(d0max + 1e-9, d0max + 1e-9 + C.PAD, 0, "pad_hi"), C.Cell(d0max + 1e-9 + C.PAD, None, 0, "end_hi")]
    a0 = d0min - 1e-9; b0 = d0max + 1e-9
    out.append(("tiling intact", cover_check(cells, a0, b0)))
    cells2 = list(cells); cells2[-1] = C.Cell(d0max + 1e-9 + 2 * C.PAD, None, 0, "end_hi")
    out.append(("tiling with a gap of one pad (must be False)", cover_check(cells2, a0, b0)))
    cells3 = [c for i, c in enumerate(cells) if i != 1]
    out.append(("tiling with a leaf removed (must be False)", cover_check(cells3, a0, b0)))
    with open(os.path.join(R.RES, "certify_arbfull_falsifier.log"), "w") as fh:
        for k, v in out: fh.write("%s: %s\n" % (k, v)); print(k, v)


if __name__ == "__main__":
    if args and args[0] == "run": main()
    if args and args[0] == "falsify": falsify()
