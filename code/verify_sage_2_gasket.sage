# Check 2: Sierpinski gasket level-n graph, independent construction; lazy walk S=(I+D^-1/2 A D^-1/2)/2.
import sys, os
sys.path.insert(0, os.getcwd())
import numpy as np
import p22lib as L
dstar = 2*log(3.0)/log(5.0)
print("2 ln3/ln5 = %.6f" % dstar)
def build(n):
    # independent: glue three copies of level n-1 at shared corners (vertex ids as frozen coordinate pairs in a 2^n integer triangular lattice via midpoint recursion)
    def tri(n, a, b, c):
        if n == 0:
            return [(a, b, c)]
        ab = tuple((x+y)//2 for x, y in zip(a, b)); bc = tuple((x+y)//2 for x, y in zip(b, c)); ca = tuple((x+y)//2 for x, y in zip(c, a))
        return tri(n-1, a, ab, ca) + tri(n-1, ab, b, bc) + tri(n-1, ca, bc, c)
    s = 2**n
    T = tri(n, (0, 0), (2*s, 0), (0, 2*s))   # coordinates doubled so midpoints stay integer
    E = set(); V = {}
    for t in T:
        for u in t:
            V.setdefault(u, len(V))
        for i in range(3):
            for j in range(i+1, 3):
                E.add(frozenset((t[i], t[j])))
    return V, E
def lam_of(n):
    V, E = build(n); N = len(V)
    import scipy.sparse as sp
    r = []; c = []
    for e in E:
        u, v = tuple(e); r += [V[u], V[v]]; c += [V[v], V[u]]
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(N, N))
    d = np.asarray(A.sum(1)).ravel()
    Dm = sp.diags(1/np.sqrt(d))
    S = 0.5*np.eye(N) + 0.5*(Dm @ A @ Dm).toarray()
    ev = np.linalg.eigvalsh(S)
    return N, len(E), ev, A, d
for n in range(1, 8):
    N, M, ev, A, d = lam_of(n)
    g = L.gasket(n)
    ok = (g.N == N) and (int(g.deg.sum())//2 == M) and (sorted(g.deg.tolist()) == sorted(d.astype(int).tolist()))
    ev2 = np.sort(g.lam)
    print("n=%d N=%d (formula %d) edges=%d degmult ok=%s  max|spec - p22lib spec|=%.2e  min eig=%.3e" % (n, N, (3**(n+1)+3)//2, M, ok, np.max(np.abs(np.sort(ev)-ev2)), ev.min()))
# independent exactness check of eigenvalues: trace identity tr(S^k)/N vs exact rational walk for n=2 (high precision)
n=2; N, M, ev, A, d = lam_of(n)
Aq = matrix(QQ, A.toarray().astype(int).tolist()); Dq = diagonal_matrix(QQ, [1/QQ(int(x)) for x in d]); Pq = Dq*Aq; Sq = (identity_matrix(QQ,N)+Pq)/2
Sk = Sq**1
for k in (1,2,3,5,8):
    tr = (Sq**k).trace()
    print("n=2 k=%d exact tr(S^k)=%s  float=%.15f  eig sum=%.15f" % (k, str(tr), float(tr), float(np.sum(ev**k))))
# d(t) on finite graph
print("\nlevel n: window = N P(t) >= 6 (bench thr), t from 24 (bench grid start); also thr 4.5 cap grid")
for n in (5, 6, 7):
    N, M, ev, A, d = lam_of(n)
    pos = ev[ev > 1e-14]; ll = np.log(pos)
    def dq(t):
        a = np.exp(t*ll); return a.sum()/N, -2*t*(a*ll).sum()/a.sum()
    ts = np.unique(np.round(np.geomspace(4, 20000, 4000)).astype(int))
    q = np.array([dq(t)[0] for t in ts]); dd = np.array([dq(t)[1] for t in ts])
    for thr in (6.0, 1.3):
        m = (N*q >= thr) & (ts >= 24)
        # contiguous from start
        idx = np.nonzero(m)[0]
        if len(idx) == 0: print("n=%d thr=%.1f empty" % (n, thr)); continue
        tt = ts[idx]; dw = dd[idx]; lt = np.log(tt)
        # log-periodic wiggle: remove smooth cubic in log t, amplitude of the remainder
        co = np.polyfit(lt, dw, 3); res = dw - np.polyval(co, lt)
        print("n=%d N=%d thr=%.1f window t=[%d,%d] (%.2f decades)  d range [%.4f,%.4f] mean %.4f  mean gap below d*=%.4f  gap range [%.4f,%.4f]  wiggle(resid cubic) +-%.4f" % (n, N, thr, tt[0], tt[-1], np.log10(tt[-1]/tt[0]), dw.min(), dw.max(), dw.mean(), dstar-dw.mean(), dstar-dw.max(), dstar-dw.min(), np.abs(res).max()))
    # peak-to-peak oscillation: local extrema of d(t) vs log t within thr 6 window and spacing (period ln5 expected in ln t)
    m = (N*q >= 1.3) & (ts >= 10); tt=ts[m]; dw=dd[m]
    ext = [i for i in range(1,len(dw)-1) if (dw[i]-dw[i-1])*(dw[i+1]-dw[i]) < 0 and abs(dw[i]-dw[i-1])>0]
    # collapse tiny jitter: keep extrema whose neighbours at +-5% differ
    print("   n=%d local extrema (t, d) in window thr1.3: %s" % (n, [(int(tt[i]), round(float(dw[i]),4)) for i in ext][:12]))
    # bench grid values
    g = L.gasket(n); big = L.make_grid(24, 34, 4096)
    qb, db = g.exact(big)
    mm = g.N*qb >= 6
    print("   n=%d bench-grid (24,34) d at Nq>=6: %s" % (n, np.round(db[mm], 3).tolist()))
