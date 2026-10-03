"""Asimov Z (120-130 GeV) at the mA3 working-point-scan cuts of tab:closure_mA3_scan, with the
signal weights of scan_score_R_significance.py (HZA_SIGNAL_REWEIGHT=0: unweighted, reproduces the AN)."""
import os, scan_score_R_significance as S
bk = S.load_scored_bkg(3); sg = S.load_scored_sig(3)
tag = "rw" if os.environ.get("HZA_SIGNAL_REWEIGHT", "1") != "0" else "norw"
print(tag, " ".join("%.3f:%.1f" % (c, S._Z_direct(bk, sg, c)) for c in (0.980, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996)))
print(tag, "R", " ".join("%.3f:%.2f" % (c, S._R_direct(bk, c)) for c in (0.980, 0.984, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996)))
