#!/usr/bin/env python3
"""Prepare the background observable histogram for the merged low-mA fit.

Background modeling is data-driven: we fit the H_mass sideband+peak region of
DATA with a smooth function family (done later, in combine env, by
build_bkg_ws.py).  This step only fills the H_mass histogram = data_obs.

DATA STATUS (blocker): no merged-selection 2024 DATA parquet exists yet.  Until
it is produced, this falls back to the sum of DY MC (weighted, = Asimov
pseudo-data) and prints a LOUD warning.  The expected limits produced from MC
pseudo-data are a first-pass closure only; the real result needs DATA.

Run in the `higgs-alp-ana` env (pyarrow + ROOT).
"""
import argparse
import glob
import os

import numpy as np
import pyarrow.parquet as pq
import ROOT

ROOT.gROOT.SetBatch(True)

EOS_DEFAULT = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
DY_MC = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]


def find_data(eos):
    pats = ["*Data*", "*EGamma*", "*DoubleMuon*", "*DoubleEG*", "*Muon*"]
    paths = []
    for p in pats:
        paths += glob.glob(f"{eos}/parquet_friend/{p}/*/merged_nominal.parquet")
    return sorted(set(paths))


def load(path, cols):
    fs = glob.glob(path)
    if not fs:
        return None
    return pq.read_table(fs[0], columns=cols).to_pandas()


def fill_hist(h, arr_h, arr_w, selection_mask):
    for x, ww in zip(arr_h[selection_mask], arr_w[selection_mask]):
        h.Fill(float(x), float(ww))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eos", default=EOS_DEFAULT)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "output"))
    ap.add_argument("--selection", default="pass_allcuts_merged_AN2020")
    ap.add_argument("--fit-lo", type=float, default=105.0)
    ap.add_argument("--fit-hi", type=float, default=170.0)
    ap.add_argument("--bins", type=int, default=65)
    ap.add_argument("--force-mc", action="store_true",
                    help="ignore any DATA parquet and use DY MC pseudo-data")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    lo, hi = args.fit_lo, args.fit_hi
    hist = ROOT.TH1D("data_obs", "data_obs", args.bins, lo, hi)
    hist.Sumw2()

    data_paths = [] if args.force_mc else find_data(args.eos)
    is_pseudo = len(data_paths) == 0

    if not is_pseudo:
        print(f"[bkg] found {len(data_paths)} DATA parquet file(s):")
        for p in data_paths:
            print("   ", p)
        for p in data_paths:
            d = pq.read_table(p, columns=["H_mass", args.selection]).to_pandas()
            mask = (d[args.selection] == 1).to_numpy()
            mh = d["H_mass"].to_numpy()
            inwin = mask & (mh > lo) & (mh < hi) & np.isfinite(mh)
            for x in mh[inwin]:
                hist.Fill(float(x), 1.0)
    else:
        print("=" * 72)
        print("!!  WARNING: no merged-selection DATA parquet found.")
        print("!!  Falling back to DY MC (weighted) as Asimov PSEUDO-DATA.")
        print("!!  Expected limits from this are a CLOSURE ONLY -- not the result.")
        print("!!  Produce 2024 merged DATA friend parquet to unblock (stage 0).")
        print("=" * 72)
        for bk in DY_MC:
            d = load(f"{args.eos}/parquet_friend/Bkg_{bk}_2024/*/merged_nominal.parquet",
                     ["H_mass", "weight_central", args.selection])
            if d is None:
                print(f"   [mc missing] {bk}")
                continue
            mask = (d[args.selection] == 1).to_numpy()
            mh = d["H_mass"].to_numpy()
            ww = d["weight_central"].to_numpy()
            inwin = mask & (mh > lo) & (mh < hi) & np.isfinite(mh) & np.isfinite(ww)
            fill_hist(hist, mh, ww, inwin)
            print(f"   [mc] {bk}: +{ww[inwin].sum():.1f} weighted in window")

    nbkg = hist.Integral()
    print(f"[bkg] total background in [{lo},{hi}] = {nbkg:.1f}  "
          f"({'pseudo-data(MC)' if is_pseudo else 'DATA'})")

    out = os.path.join(args.outdir, "bkg_hist.root")
    f = ROOT.TFile(out, "RECREATE")
    hist.Write()
    meta = ROOT.TNamed("source", "pseudo_mc" if is_pseudo else "data")
    meta.Write()
    f.Close()
    print(f"[bkg] wrote {out}")
    # marker file so the driver / datacard can flag the result
    with open(os.path.join(args.outdir, "bkg_source.txt"), "w") as fh:
        fh.write("pseudo_mc\n" if is_pseudo else "data\n")


if __name__ == "__main__":
    main()
