"""Analysis of rt_monolp.csv and rt_rbc.csv (+ stored archive data). Prints markdown tables."""
import sys, numpy as np, pandas as pd
from scipy.stats import beta
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 50)
def cp(k, n, a=0.05):
    if n == 0: return (np.nan, np.nan)
    return (0.0 if k == 0 else beta.ppf(a/2, k, n-k+1), 1.0 if k == n else beta.ppf(1-a/2, k+1, n-k))
def fmt(k, n): lo, hi = cp(k, n); return f"{k}/{n} = {100*k/max(n,1):.1f}% ({100*lo:.1f}-{100*hi:.1f})"
which = sys.argv[1]
if which == "monolp":
    fs = sys.argv[2:] or ['results/rt_monolp_part1.csv', 'results/rt_monolp_p2.csv']
    d = pd.concat([pd.read_csv(f, low_memory=False) for f in fs]).drop_duplicates(['cell', 'ds', 'variant']); d = d[d.nowin == 0] if 'nowin' in d else d
    print("rows", len(d), "datasets", d.groupby(['cell','ds']).ngroups)
    d['dim'] = np.where(d.fam == 'aniso2_w16', '2D', np.where(d.fam == 'aniso3_w8', '3D', 'ctl'))
    ans = d[d.infeasible == 0]
    rows = []
    for (fam, Lx, v), g in d.groupby(['fam', 'L', 'variant']):
        a = g[g.infeasible == 0]
        rows.append(dict(fam=fam, L=Lx, var=v, n=len(g), infeas=int(g.infeasible.sum()), cov=a.cov_simult.mean() if len(a) else np.nan, width=a.width.mean(),
                         e7_neg=int(a.e7_neg.sum()), j_drop=int(a.j_drop.sum()), j_rise=int(a.j_rise.sum()), sec_drop_ge25=int((g.sec_drop > 0.25).sum()), gateflag=int(a.j_gate_flag.sum())))
    t = pd.DataFrame(rows); t.to_csv('results/rt_monolp_cells.csv', index=False)
    print(t.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    # pooled controls coverage / false joint decl
    for v in ('B0', 'C1', 'C2'):
        g = d[(d.variant == v)]; c = g[g.kind == 'control']; a = c[c.infeasible == 0]
        print(v, "controls: answered", len(a), "of", len(c), "cov", fmt(int(a.cov_simult.sum()), len(a)), "false joint any-sign", int(a.j_drop.sum() + a.j_rise.sum()), "e7 false neg", int(a.e7_neg.sum()))
        cr = g[(g.kind == 'crossover')]; a = cr[cr.infeasible == 0]
        print(v, "crossovers: answered", len(a), "of", len(cr), "cov", fmt(int(a.cov_simult.sum()), len(a)))
    # paired width change, 2D L>=192 and 3D
    for name, sel in (('2D L>=192', (d.fam == 'aniso2_w16') & (d.L >= 192)), ('2D all', d.fam == 'aniso2_w16'), ('3D L>=96', (d.fam == 'aniso3_w8') & (d.L >= 96)), ('controls', d.kind == 'control')):
        s = d[sel]; p = s.pivot_table(index=['cell', 'ds'], columns='variant', values='width')
        for v in ('C1', 'C2'):
            if v not in p.columns or 'B0' not in p.columns: continue
            r = (p[v] / p['B0']).dropna(); print(f"width ratio {v}/B0 {name}: median {r.median():.4f} mean {r.mean():.4f} n={len(r)}")
    # joint power change paired at 2D 768, 3D
    for (fam, Lx) in (('aniso2_w16', 768), ('aniso2_w16', 512), ('aniso3_w8', 128), ('aniso3_w8', 160)):
        s = d[(d.fam == fam) & (d.L == Lx)]
        if not len(s): continue
        p = s.pivot_table(index=['cell', 'ds'], columns='variant', values='j_drop')
        print(fam, Lx, {v: int(p[v].sum()) for v in p.columns}, "of", len(p))
if which == "rbc":
    d = pd.read_csv('results/rt_rbc.csv', low_memory=False)
    d = d[d.abstain == 0]
    print("rows", len(d), "cells", d.groupby(['fam', 'L']).ngroups, "datasets/cell/budget", d.groupby(['fam','L','wpn','est']).size().unique())
    lq = d[d.est == 'LQ']
    print("identity with frozen est_E1b (cov, h, hw):", int(lq.ident_frozen.sum()), "/", len(lq))
    for est in ('LQ', 'LQU', 'RBC'):
        for w in (8, 32, 128):
            g = d[(d.est == est) & (d.wpn == w)]; n = len(g)
            print(f"{est:4s} wpn {w:3d}: cov exact {fmt(int(g.cov_simult.sum()), n)} | cov own-target {fmt(int(g.cov_pop_simult.sum()), n)} | recentred {fmt(int(g.rec_cov.sum()), n)} | pointwise {g.cov_point.mean():.3f} | mean hw {g.hw.mean():.3f} | median worst|bias|/se {g.worst_ratio.median():.2f} | mean h {g.h.mean():.2f} | h at floor(0.7) {(g.h<=0.7+1e-9).mean():.2f} (LQ/RBC only)")
    # per cell for LQ vs RBC
    t = d.groupby(['fam', 'L', 'wpn', 'est']).cov_simult.mean().unstack('est').round(2); print(t.to_string())
    # worst-time bias / empirical sampling SD of (dhat - dpop) at that time, LQ
    out = []
    for (fam, Lx, w), g in d[d.est == 'LQ'].groupby(['fam', 'L', 'wpn']):
        B = np.array([[float(x) for x in s.split(';')] for s in g.b_vec]); 
        J = min(len(r) for r in B); 
        Dv = np.array([[float(x) for x in s.split(';')][:J] for s in g.dev_vec]); B = np.array([b[:J] for b in B])
        sd = Dv.std(0, ddof=1); bm = np.abs(B.mean(0)); k = (bm / sd).argmax(); out.append(dict(fam=fam, L=Lx, wpn=w, worst_bias=bm[k], sd_at_worst=sd[k], ratio=bm[k] / sd[k], t_index=k))
    o = pd.DataFrame(out); print(o.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(o.groupby('wpn').ratio.agg(['median', 'min', 'max']))
    # stored comparison on the same cells/datasets
    a = pd.read_csv('archive/data/raw_amend10c.csv', low_memory=False); a = a[(a.est == 'E1q') & (a['mode'] == 'fixed') & (a.abstain == 0)]
    key = lq[['cell', 'ds']]; s = a.merge(key, on=['cell', 'ds'])
    print("stored E1q on the same cell-datasets:", s.groupby('wpn').cov_simult.agg(['mean', 'size']).round(3).to_dict(), "regenerated LQ:", lq.groupby('wpn').cov_simult.mean().round(3).to_dict())
    m = lq.merge(s[['cell', 'ds', 'cov_simult', 'h', 'hw']], on=['cell', 'ds'], suffixes=('', '_st'))
    print("regeneration reproduces stored cov/h/hw:", int(((m.cov_simult == m.cov_simult_st) & (m.h == m.h_st) & ((m.hw - m.hw_st).abs() < 1e-9)).sum()), "/", len(m))
    # noise-free check numbers at t=96 on C256xC16 not in this file
