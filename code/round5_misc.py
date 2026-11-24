#!/usr/bin/env python3
"""ROUND 5 small checks from stored results. (1) E2/E0 control abstention denominators (Gemini r5 item). (2) E7 endpoint test, relaxations A and B, and E7j (B only),
3D cells L=160/192/256 (ChatGPT r5). Usage: round5_misc.py OUTDIR"""
import os, sys, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); out = sys.argv[1]; os.makedirs(out, exist_ok=True); L = []
r = pd.read_csv(os.path.join(H, "results_round3", "r3_native_E0_E2_all_cells.csv"))
for e in ("E0", "E2"):
    c = r[(r.est == e) & (r.kind == "control")]; n = int(c.n.sum()); a = int(c.answered.sum()); k = int(c.k_cov_native.sum())
    L.append(f"{e} controls: {len(c)} cells, n={n}, answered={a}, covered={k}, abstained={n-a} = {100*(n-a)/n:.1f}%  (covered/answered {k}/{a}; the 1108 of the referee's reading is the COVERED count, not the answered count)")
d = pd.read_csv(os.path.join(H, "results_amend10_11", "raw_amend10b.csv"), low_memory=False)
for c in ("B_e7_neg", "B_e7_pos", "A_e7_neg", "A_e7_pos", "Bj_neg"): d[c] = d[c].fillna(0)
for Lq, g in d[(d.fam == "aniso3_w8") & d.L.isin([160, 192, 256])].groupby("L"):
    A = int(((g.A_status == 0) & ((g.A_e7_neg + g.A_e7_pos) > 0)).sum()); B = int(((g.B_status == 0) & ((g.B_e7_neg + g.B_e7_pos) > 0)).sum()); J = int(((g.B_status == 0) & (g.Bj_neg == 1)).sum())
    L.append(f"3D L={int(Lq)} n={len(g)}: endpoint test relaxation A {A}/{len(g)}, relaxation B {B}/{len(g)}; joint test E7j (B only) {J}/{len(g)}; no band A {int((g.A_status!=0).sum())}, B {int((g.B_status!=0).sum())}")
open(os.path.join(out, "round5_misc.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
