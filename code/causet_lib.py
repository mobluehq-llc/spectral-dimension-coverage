"""P22 causal-set extension: sprinkling into flat Minkowski causal diamonds (1+1, 3+1), graph constructions, exact lazy-walk q(t).
Sprinkling convention for 1+1 follows p14_lorentz_sprinkling.py (light-cone coordinates u,v; x<y iff u and v both increase); the 3+1 diamond and
all graph code are new. See PLAN.md AMENDMENT 2."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import p22lib as L

KS = (0, 3, 15)          # a pair is an edge iff it is causally related with <= k elements strictly between (k=0: Hasse diagram links)
T_CAP = 1024
GRID2 = L.make_grid(2, 3, T_CAP)                      # primary grid {2,3}*2^j  (2 t_k = t_{k+2})
GRID4 = np.array(sorted({b * 2 ** j for j in range(0, 9) for b in (4, 5, 6, 7) if b * 2 ** j <= T_CAP}), dtype=np.int64)  # dense grid (2 t_k = t_{k+4})
GRID = np.array(sorted(set(GRID2.tolist()) | set(GRID4.tolist())), dtype=np.int64)  # q is stored on the union


def sprinkle_relation(dim, N, rng):
    """returns C (N x N float32), C[i,j]=1 iff i precedes j, points sorted by time so C is strictly upper triangular."""
    if dim == 1:
        u = rng.random(N); v = rng.random(N)
        o = np.argsort(u); u, v = u[o], v[o]
        C = ((u[None, :] > u[:, None]) & (v[None, :] > v[:, None])).astype(np.float32)
        return C
    pts = []
    n = 0
    while n < N:                      # rejection sampling in the double cone |x| < 1-|t|, t in (-1,1)
        m = max(4 * (N - n), 1000)
        t = rng.uniform(-1, 1, m); x = rng.uniform(-1, 1, (m, 3))
        ok = np.linalg.norm(x, axis=1) < 1 - np.abs(t)
        pts.append(np.column_stack([t[ok], x[ok]])); n += ok.sum()
    P = np.concatenate(pts)[:N]
    P = P[np.argsort(P[:, 0])]
    C = np.zeros((N, N), dtype=np.float32)
    for a in range(0, N, 512):
        dt = P[None, :, 0] - P[a:a + 512, None, 0]
        dx = np.linalg.norm(P[None, :, 1:] - P[a:a + 512, None, 1:], axis=2)
        C[a:a + 512] = (dt > dx)
    return C


def lazy_q(A, lazy=0.5):
    """A symmetric scipy sparse or dense 0/1 adjacency (N x N). Largest component kept. returns (q on GRID, N_eff, mean degree)."""
    As = sp.csr_matrix(A)
    nc, lab = connected_components(As, directed=False)
    if nc > 1:
        big = np.argmax(np.bincount(lab)); keep = np.nonzero(lab == big)[0]
        As = As[keep][:, keep]
    N = As.shape[0]
    d = np.asarray(As.sum(1)).ravel()
    Dm = sp.diags(1.0 / np.sqrt(d))
    S = (lazy * np.eye(N) + (1 - lazy) * (Dm @ As @ Dm).toarray())
    lam = np.linalg.eigvalsh(S)
    lp = lam[lam > 1e-14]; ll = np.log(lp)
    q = np.array([np.exp(t * ll).sum() / N for t in GRID])
    return q, N, float(d.mean())


def rgg_adj(D, N, deg_target, rng):
    """random geometric graph in the unit D-cube (no wrap), radius bisected so the mean degree matches deg_target."""
    P = rng.random((N, D)); tr = cKDTree(P)
    lo, hi = 0.0, 1.0
    for _ in range(14):
        r = 0.5 * (lo + hi)
        deg = (tr.count_neighbors(tr, r) - N) / N
        lo, hi = (r, hi) if deg < deg_target else (lo, r)
    pr = tr.query_pairs(0.5 * (lo + hi), output_type="ndarray")
    A = sp.csr_matrix((np.ones(len(pr) * 2), (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(N, N))
    return A


def job(args):
    dim, N, rep = args
    rng = np.random.default_rng(np.random.SeedSequence([22, 900 + dim, N, rep]))
    C = sprinkle_relation(dim, N, rng)
    M = C @ C                                   # number of elements strictly between i and j (exact in float32 for N<2^24)
    out = dict(dim=dim, N=N, rep=rep)
    out["rel_frac"] = float(C.sum() / (N * (N - 1) / 2))
    out["full_deg"] = float(2 * C.sum() / N)
    for k in KS:
        A = (C > 0) & (M <= k)
        As = sp.csr_matrix(A.astype(np.float32)); As = ((As + As.T) > 0).astype(np.float64)
        q, Ne, md = lazy_q(As)
        out[f"q_k{k}"] = q; out[f"Neff_k{k}"] = Ne; out[f"deg_k{k}"] = md
        Ar = rgg_adj(2 if dim == 1 else 4, N, md, rng)
        qr, Ner, mdr = lazy_q(Ar)
        out[f"q_rgg{k}"] = qr; out[f"Neff_rgg{k}"] = Ner; out[f"deg_rgg{k}"] = mdr
    return out
