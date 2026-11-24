# ROUND 6 (ChatGPT r6 s3): define the 27 strata of the round-5 SD(z) check and report the range per variant.
# Stratum = (window length J, window point j) of the stored fresh-seed run (ring L=512, wpn 32, 1000 datasets), kept if >= 100 datasets.
# Only J=13 (359 datasets, 13 points) and J=14 (569 datasets, 14 points) qualify: 13+14 = 27 strata. Usage: round6_sdz_strata.py <markowitz_fresh_perj.csv>
import csv, sys, collections, numpy as np
D = collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[1])): D[(r['v'], int(r['J']), int(r['j']))].append((float(r['n']), float(r['sd'])))
for v in ('V0_baseline', 'V1_B1000', 'V2_split_A'):
    per = [(J, j, len(a), np.std(a[:, 0] / a[:, 1], ddof=1)) for (vv, J, j), L in sorted(D.items()) if vv == v for a in [np.array(L)] if len(a) >= 100]
    s = [p[3] for p in per]
    print(v, 'strata', len(per), 'J values', sorted(set(p[0] for p in per)), 'datasets per J', sorted(set((p[0], p[2]) for p in per)), 'SD(z) range %.3f-%.3f' % (min(s), max(s)), 'median %.3f' % np.median(s))
