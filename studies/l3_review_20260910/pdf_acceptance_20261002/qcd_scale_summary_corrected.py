"""QCD-scale summary over the 70 (m_a, era) points, as quoted in AN Sec. 9.1 / TWiki, before and
after replacing mA_M7_2024 by its NanoAOD-based recheck (results/qcd_scale_recheck_from_nano.json)."""
import json, numpy as np
Q = json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/qcd_scale_acceptance.json"))
F = json.load(open("results/qcd_scale_recheck_from_nano.json"))
def summ(rows, tag):
    au = np.array([r["acc_up"] for r in rows.values()]); ad = np.array([r["acc_dn"] for r in rows.values()])
    yu = np.array([r["yield_up"] for r in rows.values()]); yd = np.array([r["yield_dn"] for r in rows.values()])
    mx = np.maximum(np.abs(au), np.abs(ad)); k = list(rows)[int(mx.argmax())]
    print("%-10s n=%d  acc mean %+.2f%%/%+.2f%%  yield mean %+.2f%%/%+.2f%%  max|dA| median %.2f%% p90 %.2f%% worst %.2f%% (%s)"
          % (tag, len(rows), 100*au.mean(), 100*ad.mean(), 100*yu.mean(), 100*yd.mean(),
             100*np.median(mx), 100*np.percentile(mx, 90), 100*mx.max(), k))
summ(Q, "original")
C = dict(Q)
for s, r in F.items():
    C[s] = dict(Q[s], acc_up=r["acc_up"], acc_dn=r["acc_dn"], yield_up=r["yield_up"], yield_dn=r["yield_dn"])
summ(C, "corrected")
json.dump(C, open("results/qcd_scale_acceptance_corrected_mA7_2024.json", "w"), indent=1)
