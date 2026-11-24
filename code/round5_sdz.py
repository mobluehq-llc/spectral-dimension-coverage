# ROUND 5 EXPLORATORY check (Gemini r5): SD of z_j=n_j/sd_boot_j from the stored round-4 fresh run (ring L=512, wpn 32, 1000 datasets, perj file on host-D).
import csv,numpy as np,collections
import sys
f=sys.argv[1] if len(sys.argv)>1 else 'markowitz_fresh_perj.csv'
D=collections.defaultdict(list)
for r in csv.DictReader(open(f)):
    D[(r['v'],int(r['J']),int(r['j']))].append((float(r['n']),float(r['sd'])))
for v in ('V0_baseline','V1_B1000','V2_split_A'):
    allz=[];sds=[];rat=[];per=[]
    for (vv,J,j),L in D.items():
        if vv!=v: continue
        a=np.array(L); z=a[:,0]/a[:,1]; allz+=list(z)
        if len(a)>=100:
            per.append((z.std(ddof=1), np.median(a[:,1])/a[:,0].std(ddof=1), z.mean(), len(a)))
    allz=np.array(allz); per=np.array(per)
    print(v,'N',len(allz),'pooled SD(z)=%.3f mean z=%.3f'%(allz.std(ddof=1),allz.mean()),'| strata',len(per),'median stratum SD(z)=%.3f (min %.3f max %.3f)'%(np.median(per[:,0]),per[:,0].min(),per[:,0].max()),'median sd-ratio %.3f'%np.median(per[:,1]), 'frac|z|>1.96: %.4f'%np.mean(np.abs(allz)>1.96))
