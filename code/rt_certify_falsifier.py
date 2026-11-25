import sys; sys.path.insert(0,'.'); sys.argv=['x','none']
import rt_certify as C, numpy as np
sp=C.prep_cells()
for c in [s for s in sp if s['fam'] in('aniso2_w16','aniso3_w8')]:
    ts,N=c['ts_w'],c['N']; Jw=int(c['exact_J']); q=c['graph'].exact(ts)[0]; eps=1e-3
    ctx=C.Ctx(ts,q*(1-eps),q*(1+eps),N,Jw); d0=float(c['d_start']); md=float(c['max_drop'])
    for dl in (md+0.1, md-0.2):
        cl=C.Cell(d0-1e-3,d0+1e-3,d0-1e-3-dl-C.MARG,'x'); r=C.certify_cell(ctx,cl); print(c['fam'],'true drop',round(md,3),'delta',round(dl,3),r['status'],r.get('phi'))
