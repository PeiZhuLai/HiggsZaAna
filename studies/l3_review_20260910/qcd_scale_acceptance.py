#!/usr/bin/env python
"""L3 item 8: QCD-scale uncertainty on the SIGNAL ACCEPTANCE, per (m_a, era).

The AN currently quotes a flat 6 % taken from elsewhere. This measures it.

Method
    A_k / A_0 = (sum_sel w_k / sum_sel w_0) / (sum_gen w_k / sum_gen w_0)
The denominator matters: without it one measures the cross-section variation as
well, which is a separate uncertainty and would be double counted. It comes from
the Runs tree, where LHEScaleSumw is already normalized to genEventSumw (index 4
is exactly 1.0), so
    sum_gen w_k / sum_gen w_0 = sum_files(genEventSumw * LHEScaleSumw_k) / sum_files(genEventSumw)

Envelope: the 7-point convention drops the two combinations where muR and muF
move in opposite directions by a factor 4 (indices 2 and 6 of the 9).
Indices: 0..8 = (muR,muF) in {0.5,1,2} x {0.5,1,2}, muF fastest.

Output: JSON + a printed table.
"""
import json, os, sys, glob
import numpy as np
import pyarrow.parquet as pq

SIG = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/qcd_scale_acceptance.json"
NAMES = ["Zero","One","Two","Three","Four","Five","Six","Seven","Eight"]
SEVEN = [0,1,3,4,5,7,8]          # drop 2 and 6
GEN = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}

rows = {}
for d in sorted(glob.glob(os.path.join(SIG, "*"))):
    smp = os.path.basename(d)
    mp = os.path.join(d, "merged_nominal.parquet")
    if not os.path.exists(mp):
        continue
    t = pq.read_table(mp, columns=["weight_central"] + ["LHEScaleWeight_" + n for n in NAMES])
    w0 = np.asarray(t["weight_central"])
    if w0.sum() == 0:
        continue
    sel = np.array([float((w0 * np.asarray(t["LHEScaleWeight_" + n])).sum()) for n in NAMES])
    sel /= sel[4]
    gen = np.asarray(GEN.get(smp, [1.0] * 9), dtype=float)
    gen = gen / gen[4] if gen[4] else gen
    acc = sel / gen
    seven = acc[SEVEN]
    rows[smp] = {
        "sel_ratio": sel.tolist(), "gen_ratio": gen.tolist(), "acc_ratio": acc.tolist(),
        "acc_up": float(seven.max() - 1.0), "acc_dn": float(seven.min() - 1.0),
        "yield_up": float(sel[SEVEN].max() - 1.0), "yield_dn": float(sel[SEVEN].min() - 1.0),
        "has_gen": smp in GEN,
    }

print("%-28s %10s %10s   %10s %10s" % ("sample", "acc_up", "acc_dn", "yield_up", "yield_dn"))
for k in sorted(rows):
    r = rows[k]
    flag = "" if r["has_gen"] else "  (no gen sums -> yield only)"
    print("%-28s %+9.2f%% %+9.2f%%   %+9.2f%% %+9.2f%%%s"
          % (k, 100*r["acc_up"], 100*r["acc_dn"], 100*r["yield_up"], 100*r["yield_dn"], flag))
if rows:
    au = np.array([r["acc_up"] for r in rows.values()]); ad = np.array([r["acc_dn"] for r in rows.values()])
    print("\nacceptance envelope over %d (m_a, era) points: up mean %+.2f%% max %+.2f%% | dn mean %+.2f%% min %+.2f%%"
          % (len(rows), 100*au.mean(), 100*au.max(), 100*ad.mean(), 100*ad.min()))
json.dump(rows, open(OUT, "w"), indent=1)
print("\nwrote", OUT)
