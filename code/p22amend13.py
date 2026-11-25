"""P22 Amendment 13: (a) independent-split bias-aware E4 variants (V2), (b) random-graph builders and ensemble truth.  See PLAN.md AMENDMENT 13.
Nothing here edits the frozen estimator files (p22lib.py, p22honest.py, p22amend*.py); it calls them.

E4 variants, all on the same walker data, same window rule, same half A / half B split as the frozen est_E4:
  V0 ("E4o", "E4c")        : frozen construction. S = bootstrap covariance of log p_hat on half B (draws Ys, rng B); weights from S; critical value from the SAME half-B draws.
  V2 ("E4o_V2", "E4c_V2")  : S = bootstrap covariance of log p_hat on half A (draws YsA, rng A, floor 0.5/W_A); weights and sd_j from that S; the critical value is the
                             0.95 quantile of max_j |a_j'(y* - y_B)| / sd_j over the half-B draws Ys (as V0), the point estimate and bounds are from half B.
  M: E4o = exact sup |g'''| over the window (oracle, needs the spectrum); E4c = held-out local-quartic rule on half A (p22honest.M_heldout). Identical to V0.
Change test for V2: contrast weights a_j - a_0, sd from S_A, max-t critical value from the half-B draws, bias bound B_j + B_0 (as V0)."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
import p22lib as L
import p22honest as H
from p22lib import grand_mean, boot_means, abst, window_len, DELTA


def _arm(tag, x, ts, J, M, S, Ys, yB, bt_list, Cc_S=None):
    """one (variant, M) arm. bt_list: list of (name, beta_true array) targets (beta = -d/2). Returns (row, dlo, dhi)."""
    W = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J)
    for j in range(J):
        r = H.honest_weights(x, x[j], S, M)
        if r is None:
            return abst(tag, "solver"), None, None
        W[j], sd[j], Bj[j] = r
    sd = np.maximum(sd, 1e-12)
    bh = W @ yB; Bs = Ys @ W.T
    c = float(np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95))
    hwj = c * sd + Bj
    row = dict(est=tag, abstain=0, J=J, hw=float(2 * hwj.mean()), se_mean=float(2 * sd.mean()), bias_bound=float(2 * Bj.mean()), M=float(M), crit=c, tlast=int(ts[J - 1]))
    names = [n for n, _ in bt_list]
    for n, bt in bt_list:
        err = bh - bt; sfx = "" if n == "S" else "_" + n
        row["cov_simult" + sfx] = int((np.abs(err) <= hwj).all()); row["cov_point" + sfx] = float((np.abs(err) <= H.Z * sd + Bj).mean())
        if n == "S":
            row["bias"] = float((-2 * err).mean()); row["mae"] = float(np.abs(2 * err).mean())
    Cc = W[1:] - W[[0]]
    sdc = np.maximum(np.sqrt(np.einsum("ij,jk,ik->i", Cc, S, Cc)), 1e-12); Bc = Bj[1:] + Bj[0]
    dc = Cc @ yB; Dc = Ys @ Cc.T
    cc = float(np.quantile((np.abs(Dc - dc) / sdc).max(1), 0.95)); hc = cc * sdc + Bc
    dl = -2 * dc; hcd = 2 * hc; lo, hi = dl - hcd, dl + hcd
    row.update(decl_neg=int((hi < -DELTA).any()), decl_pos=int((lo > DELTA).any()), decl_zero=int(((hi < 0) | (lo > 0)).any()), dmin=float(dl.min()), dmax=float(dl.max()))
    for n, bt in bt_list:
        sfx = "" if n == "S" else "_" + n
        dtt = -2 * (bt[1:] - bt[0]); row["delta_cov" + sfx] = int(((dtt >= lo) & (dtt <= hi)).all())
        if n == "S": row["true_dmin"] = float(dtt.min())
    return row, -2 * bh - 2 * hwj, -2 * bh + 2 * hwj


def est_E4v2(Ublk, mwalk, ts, N, rngB, rngA, truths, Mo_cum, thr=L.THR, B=200, tags=("E4o", "E4c"), variants=("V0", "V2"), keep_bands=False):
    """truths: list of (name, d_exact_function) ; name 'S' is the primary target (columns without suffix), others get a column suffix _<name>.
    Returns list of rows (one per tag x variant) in the order E4o, E4c, E4o_V2, E4c_V2 (those requested)."""
    out_tags = [t + ("" if v == "V0" else "_V2") for v in variants for t in tags]
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA); WA = sum(u.shape[0] for u in UA) * mwalk; cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return [abst(t, "window") for t in out_tags]
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    Ms = {}
    if "E4o" in tags: Ms["E4o"] = max(float(Mo_cum[J - 1]), H.FLOORM)
    if "E4c" in tags: Ms["E4c"] = H.M_heldout(x, yA, wA)
    WB = sum(u.shape[0] for u in UB) * mwalk; floorB = 0.5 / WB
    yB = np.log(np.maximum(grand_mean(UB), floorB)[:J])
    Ys = np.log(np.maximum(boot_means(UB, B, rngB)[:, :J], floorB))
    Ss = {"V0": np.cov(Ys, rowvar=False)}
    if "V2" in variants:
        YsA = np.log(np.maximum(boot_means(UA, B, rngA)[:, :J], 0.5 / WA)); Ss["V2"] = np.cov(YsA, rowvar=False)
    bts = [(n, -f(ts[:J]) / 2.0) for n, f in truths]
    rows = []
    for v in variants:
        for t in tags:
            tag = t + ("" if v == "V0" else "_V2"); M = Ms[t]
            if M is None:
                rows.append(abst(tag, "cvfail")); continue
            r, lo, hi = _arm(tag, x, ts, J, M, Ss[v], Ys, yB, bts)
            if keep_bands and lo is not None:
                r["e4lo"] = ";".join("%.5f" % q for q in lo); r["e4hi"] = ";".join("%.5f" % q for q in hi)
            rows.append(r)
    return rows


# ---------------------------------------------------------------------------- random graphs (arm b)
def _graph_from_adj(adj_sets):
    N = len(adj_sets); deg = np.array([len(a) for a in adj_sets]); md = int(deg.max())
    nbr = -np.ones((N, md), dtype=np.int64)
    for i, a in enumerate(adj_sets): nbr[i, :len(a)] = sorted(a)
    return L.Graph("step", 0.5, nbr=nbr, deg=deg)


def _largest_component(A):
    """A scipy sparse symmetric adjacency (csr). Returns relabelled adjacency sets of the largest connected component."""
    nc, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax(); keep = np.nonzero(lab == big)[0]
    sub = A[keep][:, keep].tocsr()
    return [set(sub.indices[sub.indptr[i]:sub.indptr[i + 1]].tolist()) for i in range(len(keep))]


def watts_strogatz(N, k, p, rng):
    """Newman-Watts-style standard WS: ring lattice, each node joined to k/2 neighbours each side; each lattice edge (i, i+s) is rewired with prob p (its far end moved to a
    uniform node, no self loops or duplicate edges); the largest connected component is returned (equal to the whole graph unless p is tiny and unlucky)."""
    adj = [set() for _ in range(N)]
    edges = []
    for s in range(1, k // 2 + 1):
        for i in range(N):
            j = (i + s) % N; adj[i].add(j); adj[j].add(i); edges.append((i, j))
    for (i, j) in edges:
        if rng.random() < p:
            for _ in range(50):
                h = int(rng.integers(0, N))
                if h != i and h not in adj[i]:
                    adj[i].discard(j); adj[j].discard(i); adj[i].add(h); adj[h].add(i); break
    rows = np.repeat(np.arange(N), [len(a) for a in adj]); cols = np.concatenate([np.fromiter(sorted(a), dtype=np.int64) for a in adj])
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N))
    return _graph_from_adj(_largest_component(A))


def site_percolation(Lside, pocc, rng):
    """supercritical 2D site percolation on an Lside x Lside square lattice, open boundaries, 4-neighbour clusters; largest cluster returned as a graph."""
    from scipy.ndimage import label
    occ = rng.random((Lside, Lside)) < pocc
    lab, n = label(occ)
    if n == 0: return None
    cnt = np.bincount(lab.ravel()); cnt[0] = 0; big = cnt.argmax()
    m = lab == big; idx = -np.ones((Lside, Lside), dtype=np.int64); idx[m] = np.arange(m.sum())
    adj = [set() for _ in range(int(m.sum()))]
    for dx, dy in ((1, 0), (0, 1)):
        a = m[:Lside - dx, :Lside - dy] & m[dx:, dy:]
        i = idx[:Lside - dx, :Lside - dy][a]; j = idx[dx:, dy:][a]
        for u, v in zip(i.tolist(), j.tolist()): adj[u].add(v); adj[v].add(u)
    return _graph_from_adj(adj)


def make_random_graph(fam, par, rng):
    if fam == "ws": return watts_strogatz(par["N"], par["k"], par["p"], rng)
    if fam == "perc": return site_percolation(par["Lside"], par["p"], rng)
    raise ValueError(fam)


def exact_curves(lam, N, ts):
    """q(t) = mean lam^t and dq(t) = mean lam^t ln lam on the grid (eigenvalues <= 1e-14 dropped, as Graph.exact)."""
    lp = lam[lam > 1e-14]; ll = np.log(lp); ts = np.asarray(ts, float)
    q = np.empty(len(ts)); dq = np.empty(len(ts))
    for i, t in enumerate(ts):
        a = np.exp(t * ll); q[i] = a.sum() / N; dq[i] = (a * ll).sum() / N
    return q, dq
