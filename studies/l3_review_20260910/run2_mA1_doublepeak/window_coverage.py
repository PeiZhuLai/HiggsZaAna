"""m_a = 1 GeV after the WP: fraction of the signal inside candidate signal-model windows."""
import glob
import numpy as np, uproot
w, m, c = [], [], []
for p in sorted(glob.glob("/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1/*.root")):
    a = uproot.open(p)["inclusive"].arrays(["factor", "H_m", "MVA_Score_mA_M1", "z_ee"], library="np")
    s = (a["MVA_Score_mA_M1"] > 0.963) & (a["H_m"] > 95) & (a["H_m"] < 180)
    w.append(a["factor"][s]); m.append(a["H_m"][s]); c.append(a["z_ee"][s])
w, m, c = map(np.concatenate, (w, m, c))
for lo, hi in ((115, 133), (115, 137), (115, 140), (115, 145), (115, 150), (110, 150), (105, 160), (95, 180)):
    print("%3d-%3d  all %.3f   ee %.3f   mumu %.3f" % (lo, hi, *[w[(m > lo) & (m < hi) & k].sum() / w[k].sum()
          for k in (np.ones_like(c, bool), c > 0, c == 0)]))
h, e = np.histogram(m, np.arange(110, 165, 2.5), weights=w)
print("\nm_llgg histogram (2.5 GeV, weighted, BDT WP):")
for lo, v in zip(e[:-1], h / h.sum()): print("  %6.1f %6.3f %s" % (lo, v, "#" * int(200 * v)))
