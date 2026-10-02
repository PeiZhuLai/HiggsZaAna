"""Diagnose events whose parquet LHEPdfWeight_Up-Unit differs from the std recomputed from the NanoAOD members."""
import sys, glob, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdf_acceptance as pa
s = sys.argv[1]
p = pa.presel_ids(s)
k = pa.keyof(p["run"], p["luminosityBlock"], p["event"])
for u in pa.input_files(s):
    z = np.load(pa.cache_path(s, u))
    m = np.isin(k, z["key"])
    if not m.any():
        print(os.path.basename(u), "no presel events"); continue
    kk = k[m]; order = np.argsort(z["key"]); idx = order[np.searchsorted(z["key"], kk, sorter=order)]
    pdf = z["pdf"][idx].astype(float)
    d = np.abs((p["LHEPdfWeight_Up"][m] - p["LHEPdfWeight_Unit"][m]) - pdf[:, 1:].std(axis=1))
    print("%-45s nsel=%5d  maxdiff=%.2e  n(>1e-3)=%d  nevents=%d" % (os.path.basename(u), m.sum(), d.max(), (d > 1e-3).sum(), int(z["nevents"])))
    if (d > 1e-3).any():
        j = np.where(d > 1e-3)[0][:5]
        print("   parquet Up-Unit", (p["LHEPdfWeight_Up"][m] - p["LHEPdfWeight_Unit"][m])[j], " nano std", pdf[j, 1:].std(axis=1), "nano maxabs dev", np.abs(pdf[j,1:]-1).max(axis=1))
