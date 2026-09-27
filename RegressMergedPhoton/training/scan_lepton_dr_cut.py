#!/usr/bin/env python3
"""
Choose the photon-lepton dR radius for the merged-ML selection.

The merged-ML branch of za_tagger_merged.py has no photon-lepton cleaning, so
the electrons of the Z are reconstructed as merged photons: in DYto2E, 98.6% of
events passing the merged selection have dR(MLPhoton, nearest lepton) < 0.3, and
the fakes outnumber genuine candidates 69:1 on the signal region. The resolved
branch cleans at 0.3, but a merged photon is a single cluster rather than a
resolved pair, so the right radius here has to be measured rather than inherited.

This scans candidate radii and reports, per mass point, what fraction of SIGNAL
survives, against what fraction of the DY background survives -- the tradeoff
that fixes the radius.

Run with hza_ana (pyarrow >= 16; LCG_104's raises ArrowNotImplementedError on
these files, which looks like corruption but is not):

    export PATH=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin:$PATH
    python scan_lepton_dr_cut.py
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

if tuple(int(x) for x in pa.__version__.split(".")[:2]) < (16, 0):
    sys.exit("pyarrow %s too old; use hza_ana" % pa.__version__)

SIG = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_merged_DNA_v4/Sig_MC_MLNANO_all"
BKG = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML"
MASSES = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9"]
BKG_TAGS = ["Bkg_DYGto2LG_10to100_2024", "Bkg_DYJetsTo2E_2024",
            "Bkg_DYJetsTo2Mu_2024", "Bkg_DYJetsTo2Tau_2024"]
RADII = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

COLS = ["MLPhoton_lead_eta", "MLPhoton_lead_phi", "MLPhoton_lead_mass",
        "Z_lead_lepton_eta", "Z_lead_lepton_phi",
        "Z_sublead_lepton_eta", "Z_sublead_lepton_phi",
        "pass_allcuts_merged_ML"]


def delta_r(e1, p1, e2, p2):
    de = e1 - e2
    dp = (p1 - p2 + np.pi) % (2 * np.pi) - np.pi
    return np.sqrt(de * de + dp * dp)


def load(path):
    """-> (dR to nearest Z lepton, MLPhoton mass) for events passing merged-ML."""
    names = pq.read_schema(path).names
    cols = [c for c in COLS if c in names]
    if "pass_allcuts_merged_ML" not in cols:
        return None
    t = pq.read_table(path, columns=cols)
    g = lambda c: t[c].to_numpy().astype(float)
    keep = g("pass_allcuts_merged_ML") > 0.5
    me, mp, mm = g("MLPhoton_lead_eta"), g("MLPhoton_lead_phi"), g("MLPhoton_lead_mass")
    d = np.minimum(
        delta_r(me, mp, g("Z_lead_lepton_eta"), g("Z_lead_lepton_phi")),
        delta_r(me, mp, g("Z_sublead_lepton_eta"), g("Z_sublead_lepton_phi")),
    )
    ok = keep & np.isfinite(d) & (mm > -900)
    return d[ok], mm[ok]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radii", default=",".join(str(r) for r in RADII))
    a = ap.parse_args()
    radii = [float(x) for x in a.radii.split(",")]

    print("SIGNAL: fraction surviving dR(MLPhoton, nearest Z lepton) > cut")
    print("%-7s %9s " % ("mass", "n") + " ".join("%7.2f" % r for r in radii))
    sig_keep = {r: [0, 0] for r in radii}
    for m in MASSES:
        p = os.path.join(SIG, "mA_MLNANO_%s_2024" % m, "merged_nominal.parquet")
        if not os.path.exists(p):
            print("%-7s  (missing)" % m)
            continue
        got = load(p)
        if got is None:
            print("%-7s  (no merged-ML column)" % m)
            continue
        d, _ = got
        row = []
        for r in radii:
            f = float((d > r).mean()) if len(d) else 0.0
            row.append(f)
            sig_keep[r][0] += int((d > r).sum())
            sig_keep[r][1] += len(d)
        print("%-7s %9d " % (m, len(d)) + " ".join("%6.1f%%" % (100 * x) for x in row))

    print()
    print("BACKGROUND: fraction surviving the same cut")
    print("%-26s %9s " % ("sample", "n") + " ".join("%7.2f" % r for r in radii))
    bkg_keep = {r: [0, 0] for r in radii}
    for tag in BKG_TAGS:
        cand = glob.glob(os.path.join(BKG, tag, "*", "merged_nominal.parquet"))
        if not cand:
            print("%-26s  (missing)" % tag)
            continue
        got = load(cand[0])
        if got is None:
            continue
        d, _ = got
        row = []
        for r in radii:
            row.append(float((d > r).mean()) if len(d) else 0.0)
            bkg_keep[r][0] += int((d > r).sum())
            bkg_keep[r][1] += len(d)
        print("%-26s %9d " % (tag.replace("Bkg_", "").replace("_2024", ""), len(d))
              + " ".join("%6.1f%%" % (100 * x) for x in row))

    print()
    print("%-14s " % "combined" + " ".join("%7.2f" % r for r in radii))
    es = [sig_keep[r][0] / max(sig_keep[r][1], 1) for r in radii]
    eb = [bkg_keep[r][0] / max(bkg_keep[r][1], 1) for r in radii]
    print("%-14s " % "signal eff" + " ".join("%6.1f%%" % (100 * x) for x in es))
    print("%-14s " % "bkg eff" + " ".join("%6.1f%%" % (100 * x) for x in eb))
    print("%-14s " % "S/sqrt(B)" + " ".join("%7.2f" % (s / np.sqrt(b) if b > 0 else 0)
                                            for s, b in zip(es, eb)))
    print()
    print("S/sqrt(B) is relative to no cut (both effs are 1.0 there), so the")
    print("column with the largest value is the best radius on this metric.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
