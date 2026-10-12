"""Run 3 m_a = 1 GeV: share of overlapping showers (gen dR(gg) < 0.05) at preselection and after the WP."""
import glob
import numpy as np
import uproot
BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1"
w, d, s, m = [], [], [], []
for f in sorted(glob.glob(BASE + "/*.root")):
    a = uproot.open(f)["inclusive"].arrays(["factor", "GenALP_dR_gg", "MVA_Score_mA_M1", "H_mass"], library="np")
    w.append(a["factor"]); d.append(a["GenALP_dR_gg"]); s.append(a["MVA_Score_mA_M1"]); m.append(a["H_mass"])
w, d, s, m = map(np.concatenate, (w, d, s, m))
inwin = (m > 95) & (m < 180)
for name, sel in (("preselection", inwin), ("BDT > 0.963", inwin & (s > 0.963))):
    print("%-14s raw %7d  weighted share gen dR<0.05 = %.3f   share in 115-133 = %.3f" % (
        name, sel.sum(), w[sel & (d < 0.05)].sum() / w[sel].sum(),
        w[sel & (m > 115) & (m < 133)].sum() / w[sel].sum()))
