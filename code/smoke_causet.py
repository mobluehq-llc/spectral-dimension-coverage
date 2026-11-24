import sys, time, numpy as np
import causet_lib as cl
for dim, N in [(1,500),(3,500),(1,2000),(3,2000),(1,8000)]:
    t=time.time(); o=cl.job((dim,N,0)); 
    print(dim,N,round(time.time()-t,1),'s', {k:round(o[k],1) for k in o if k.startswith('deg') or k=='full_deg'}, 'Nq(t=2,6,24):',[round(N*o['q_k0'][i],1) for i in (0,3,7)], flush=True)
