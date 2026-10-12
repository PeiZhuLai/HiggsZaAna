"""Sanity checks of a Run 2 R1 output (see the Run 2 port plan)."""
import glob, os, sys
import numpy as np
import pyarrow.parquet as pq
out = sys.argv[1]
nom = pq.read_table(os.path.join(out, "output_job_1_nominal.parquet")).to_pandas()
print("nominal: %d events, columns %d" % (len(nom), nom.shape[1]))
w = [c for c in nom.columns if c.startswith("weight_") and c.endswith("_central")]
for c in sorted(w):
    v = nom[c].to_numpy()
    print("  %-62s mean %.4f  min %.4f  max %.4f" % (c, v.mean(), v.min(), v.max()))
for c in ("LHEScaleWeight_Four", "LHEScaleWeight_Zero", "LHEScaleWeight_Eight", "LHEPdfWeight_Up"):
    if c in nom:
        v = nom[c].to_numpy(); print("  %-28s mean %.4f  all==1: %s" % (c, v.mean(), bool(np.all(v == 1))))
dummy = [c for c in nom.columns if c.startswith("ALP_") and np.all(nom[c].to_numpy() == -999)]
print("  ALP columns all dummy (-999):", dummy)
for f in sorted(glob.glob(os.path.join(out, "*_up.parquet")) + glob.glob(os.path.join(out, "*_down.parquet"))):
    t = pq.read_table(f, columns=["Z_lead_lepton_pt", "Z_sublead_lepton_pt", "ALP_lead_photon_pt", "H_mass"]).to_pandas()
    print("  %-34s n=%5d  min lep pt %.2f  min pho pt %.2f  mean m_llgg %.3f" % (
        os.path.basename(f)[13:-8], len(t), min(t.Z_lead_lepton_pt.min(), t.Z_sublead_lepton_pt.min()),
        t.ALP_lead_photon_pt.min(), t.H_mass.mean()))
print("  nominal mean m_llgg %.3f" % nom.H_mass.mean())
