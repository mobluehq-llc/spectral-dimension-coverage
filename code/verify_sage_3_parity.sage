# Check 3: parity trap. Non-lazy walk on a bipartite graph has P(odd t)=0; the lazy walk does not.
# Symbolic part: eigenvalues of C_n (n even) come in +- pairs (cos(2 pi (j+n/2)/n) = -cos(2 pi j/n)), so every odd power sum vanishes.
var('n j'); 
print("symbolic: cos(2*pi*(j+n/2)/n) =", bool(cos(2*pi*(j+n/2)/n).simplify_full() == -cos(2*pi*j/n)))
# Exact (integer/rational) check on several bipartite graphs: tr(A^t)=0 for odd t; for the lazy walk tr((I+P)/2)^t) >= N 2^-t > 0.
def check(G, name, tmax=21):
    assert G.is_bipartite()
    A = G.adjacency_matrix(); N = G.order()
    d = G.degree(); regular = len(set(d)) == 1
    Dinv = diagonal_matrix(QQ, [1/QQ(x) for x in d]); P = Dinv*A
    I = identity_matrix(QQ, N); S = (I+P)/2
    bad_nl = [t for t in range(1, tmax+1, 2) if (P**t).trace() != 0]
    even_pos = all((P**t).trace() > 0 for t in range(2, tmax+1, 2))
    lazy_pos = all((S**t).trace() >= N/QQ(2)**t and (S**t).trace() > 0 for t in range(1, tmax+1))
    print("%-22s N=%3d  nonlazy tr(P^t)=0 for all odd t<=%d: %s | even t>0: %s | lazy P(t)>=2^-t>0 all t: %s" % (name, N, tmax, bad_nl == [], even_pos, lazy_pos))
check(graphs.CycleGraph(8), "C_8"); check(graphs.CycleGraph(64), "C_64", 11)
check(graphs.CycleGraph(6).cartesian_product(graphs.CycleGraph(4)), "C_6 x C_4")
check(graphs.GridGraph([4,4,4]) if False else graphs.CycleGraph(4).cartesian_product(graphs.CycleGraph(4)).cartesian_product(graphs.CycleGraph(4)), "C_4^3")
# counterexample: odd cycle (not bipartite) has odd-t returns
G = graphs.CycleGraph(7); P = (1/QQ(2))*G.adjacency_matrix()
print("C_7 (non-bipartite) tr(P^7)/7 =", (P**7).trace()/7, "(nonzero: parity trap needs bipartite)")
# Closed form on the cycle and torus, via eigenvalues, with high precision: odd/even and lazy removal
R = RealField(100)
nn = 256; lam = [R(cos(2*R(pi)*k/nn)) for k in range(nn)]
print("C_256 non-lazy P(t=35)=%s  P(t=34)=%.6f   lazy P(35)=%.6f" % (sum(x**35 for x in lam)/nn, float(sum(x**34 for x in lam)/nn), float(sum(((1+x)/2)**35 for x in lam)/nn)))
# P_lazy(t) = sum_k C(t,k) 2^-t P_k ; for bipartite nonlazy P_k=0 for odd k so lazy equals binomial mixture of even-time returns:
t=35; mix = sum(binomial(t,k)*R(2)**(-t)*sum(x**k for x in lam)/nn for k in range(0,t+1))
print("lazy binomial mixture identity at t=35: %.10f vs %.10f" % (float(mix), float(sum(((1+x)/2)**35 for x in lam)/nn)))
