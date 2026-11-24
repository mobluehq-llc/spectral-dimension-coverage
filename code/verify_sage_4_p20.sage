# Check 4 (p20): Delta_p = -Cov_w(u, ln(1+q u))/Var_w(u)  vs exact weighted least squares of ln(lambda_i) on ln(x_c).
# lambda_i = int_{bin i} x^alpha (1+q ln x) dx (closed form in s=ln x), 20 log bins on [xmin,xmax], x0=1, u_i = ln x_c,i.
import itertools
R = RealField(100)
def lam_bins(al, q, xmin, xmax, NB=20):
    a=R(al); q=R(q); E=[R(xmin)*(R(xmax)/R(xmin))**(R(i)/NB) for i in range(NB+1)]
    b=a+1
    def F(s):   # antiderivative in s of exp(b s)(1+q s)
        if abs(b) < 1e-30: return s + q*s*s/2
        return exp(b*s)*((1+q*s)/b - q/b**2)
    S=[log(e) for e in E]
    lam=[F(S[i+1])-F(S[i]) for i in range(NB)]
    u=[(S[i]+S[i+1])/2 for i in range(NB)]
    return lam,u
def wls_slope(lam,u,w):
    W=sum(w); mu=sum(wi*ui for wi,ui in zip(w,u))/W; z=[log(l) for l in lam]; mz=sum(wi*zi for wi,zi in zip(w,z))/W
    return sum(wi*(ui-mu)*(zi-mz) for wi,ui,zi in zip(w,u,z))/sum(wi*(ui-mu)**2 for wi,ui in zip(w,u))
def formula(lam,u,q,w):
    W=sum(w); m=lambda f: sum(wi*fi for wi,fi in zip(w,f))/W
    g=[log(1+R(q)*ui) for ui in u]
    return -(m([a*b for a,b in zip(u,g)])-m(u)*m(g))/(m([a*a for a in u])-m(u)**2)
rows=[]; first_bad={}
for al in (-1.5,-2.5):
  for (xmin,xmax,nm) in ((1.5,150,'2dec'),(1.5,15,'1dec'),(1.5,1500,'3dec'),(1.1,110,'2dec_x1.1')):
    for wn in ('lam','uniform','lam2'):
      for q in (-0.15,-0.1,-0.05,0.02,0.05,0.2,0.5,1.0,2.0):
        lam,u=lam_bins(al,q,xmin,xmax)
        if min(1+R(q)*ui for ui in u) <= 0: continue
        w = lam if wn=='lam' else ([R(1)]*20 if wn=='uniform' else [l*l for l in lam])
        b=wls_slope(lam,u,w); exact=-(b-(al+1))   # Delta p = -(slope-1-alpha) (code convention: slope in ln lambda vs ln x is alpha+1 under unit-width log bins)
        f=formula(lam,u,q,w)
        # pointwise-model exact line (no bin-integration): lambda_i = xc^(al+1)(1+q u)
        lam_pt=[exp((R(al)+1)*ui)*(1+R(q)*ui) for ui in u]
        b2=wls_slope(lam_pt,u,w); exact_pt=-(b2-(al+1))
        rel=abs(f-exact)/abs(exact); relpt=abs(f-exact_pt)/abs(exact_pt)
        rows.append((al,nm,wn,q,float(f),float(exact),float(exact_pt),float(rel),float(relpt),float(-q)))
print("alpha range weights q | formula exact_binned exact_pointwise | rel err (binned) rel err (pointwise) | first-order -q")
for r in rows:
    if r[2]=='lam' or r[7]>0.01: print("%5.1f %-9s %-8s %6.2f | %9.5f %9.5f %9.5f | %8.2e %8.2e | %6.2f rel(-q) %8.2e" % (r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],abs(r[9]-r[5])/abs(r[5])))
import collections
print("\nSUMMARY: n=%d configs; formula vs exact binned WLS: max rel err %.3e; configs with err>1%%: %d" % (len(rows), max(r[7] for r in rows), sum(r[7]>0.01 for r in rows)))
print("formula vs pointwise-model exact WLS (algebraic identity expected to be exact when weights = lam_pt): max rel err %.3e" % max(r[8] for r in rows))
print("first-order Delta p = -q vs exact binned: configs with err>1%%: %d / %d, min err %.2e" % (sum(abs(r[9]-r[5])/abs(r[5])>0.01 for r in rows), len(rows), min(abs(r[9]-r[5])/abs(r[5]) for r in rows)))
print("by weights: " + ", ".join("%s max %.2e" % (w, max(r[7] for r in rows if r[2]==w)) for w in ('lam','uniform','lam2')))
print("by q (w=lam, 2dec): " + ", ".join("q=%.2f:%.1e" % (r[3],r[7]) for r in rows if r[2]=='lam' and r[1]=='2dec' and r[0]==-1.5))
