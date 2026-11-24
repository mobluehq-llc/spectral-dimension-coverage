# Check 1: exact (120-bit) torus truth vs numpy p22lib truth.  Run: sage verify_sage_1_torus.sage  (p22lib.py must sit beside it)
import sys, os
sys.path.insert(0, os.getcwd())
import numpy as np
import p22lib as L
R = RealField(120)
def eigs(dims, lazy=R(1)/2):
    m = len(dims)
    cs = [[R(cos(2*R(pi)*j/n)) for j in range(n)] for n in dims]
    out = [R(0)]
    out = None
    from itertools import product
    res = []
    for idx in product(*[range(n) for n in dims]):
        s = sum(cs[a][idx[a]] for a in range(m))
        res.append(lazy + (1-lazy)*s/m)
    return res
def exactP(lam, ts):
    N = len(lam)
    P=[]; D=[]
    lp=[x for x in lam if x > R(10)**(-40)]
    ll=[log(x) for x in lp]
    for t in ts:
        a=[x**int(t) for x in lp]
        S=sum(a); P.append(S/N)
        D.append(-2*t*sum(ai*li for ai,li in zip(a,ll))/S)
    return P, D
cases = [((64,64),'C64xC64'),((128,128),'C128xC128'),((64,16),'C64xC16'),((128,16),'C128xC16'),((64,8,8),'C64xC8xC8'),((128,8,8),'C128xC8xC8'),((256,),'C256')]
worst = 0; worstP=0; worstS=0
for dims, name in cases:
    lam = eigs(dims); N=len(lam)
    g = L.Graph("torus", 0.5, dims=tuple(dims))
    big = L.make_grid(24, 34, 3000000)
    # bench window: N*P >= thr-1.5 (=4.5) -> cap, as in run_p22.prep
    qn, dn = g.exact(big)
    idx = np.nonzero(N*qn >= 4.5)[0]
    cap = int(big[idx[-1]]) if len(idx) else int(big[0])
    ts = L.make_grid(24, 34, cap)
    P, D = exactP(lam, [int(t) for t in ts])
    qn, dn = g.exact(ts)
    dD = max(abs(float(D[i]) - dn[i]) for i in range(len(ts)))
    dP = max(abs(float(P[i]) - qn[i])/qn[i] for i in range(len(ts)))
    # secant d on the same grid (finite differences of log P vs log t)
    sec_s = [ -2*(log(P[i+1])-log(P[i]))/(log(R(ts[i+1]))-log(R(ts[i]))) for i in range(len(ts)-1)]
    sec_n = -2*np.diff(np.log(qn))/np.diff(np.log(ts.astype(float)))
    dS = max(abs(float(sec_s[i]) - sec_n[i]) for i in range(len(sec_s)))
    worst=max(worst,dD); worstP=max(worstP,dP); worstS=max(worstS,dS)
    print("%-12s N=%6d  ngrid=%3d tcap=%8d  max|d_exact-d_numpy|=%.3e  max relP=%.3e  max|secant diff|=%.3e  d[0]=%.6f d[-1]=%.6f" % (name,N,len(ts),cap,dD,dP,dS,float(D[0]),float(D[-1])))
print("WORST  d_analytic %.3e  relP %.3e  secant %.3e" % (worst,worstP,worstS))
