#!/usr/bin/env python3
"""FSR recovery: dilepton and four-body mass before and after adding the FSR photon.

Answers the L3 question on Section 4.4 ("share before/after plots demonstrating
how adding up the FSR photon improves the mass resolution"). Also reports how
often the recovery actually fires, which is the quantity that turned out to be
the interesting one.

Input : /eos/home-p/pelai/HZa/parquet_DNA/Sig_MC/mA_M<ma>_<era>/merged_nominal.parquet
        (columns Z_mass / Z_noFSR_mass, H_mass / H_noFSR_mass, n_fsr, z_mumu)
Output: Plot/plots/fsrRecovery/<era>/fsr_<var>_mA_M<ma>_<era>.pdf
        Plot/plots/fsrRecovery/fsr_summary.txt

Usage : python plot_fsr_recovery.py [--ma 1,5,15,30] [--year 2024] [--fsr-only]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import ROOT  # noqa: E402
from _cms_style import LINE_WIDTH, cms_label, format_axes, make_legend  # noqa: E402

PARQUET_BASE = Path("/eos/home-p/pelai/HZa/parquet_DNA/Sig_MC")
OUT_BASE = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/fsrRecovery")
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]

COLUMNS = ["z_mumu", "n_fsr", "Z_mass", "Z_noFSR_mass",
           "H_mass", "H_noFSR_mass", "weight_central"]

VARIANTS = {
    "mll": dict(after="Z_mass", before="Z_noFSR_mass", nbins=60, lo=60.0, hi=120.0,
                title="m_{#font[12]{ll}} [GeV]"),
    "mllgg": dict(after="H_mass", before="H_noFSR_mass", nbins=42, lo=95.0, hi=180.0,
                  title="m_{#font[12]{ll}#gamma#gamma} [GeV]"),
}
_KEEP = []


def sigma_eff(values, weights=None):
    """Half-width of the narrowest interval containing 68.3% of the entries."""
    values = np.sort(np.asarray(values, dtype=float))
    n = len(values)
    if n < 10:
        return float("nan")
    k = int(round(0.683 * n))
    widths = values[k:] - values[:n - k]
    return float(np.min(widths) / 2.0)


def fill_hist(name, values, weights, cfg):
    hist = ROOT.TH1D(name, "", cfg["nbins"], cfg["lo"], cfg["hi"])
    hist.Sumw2()
    values = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    weights = np.ascontiguousarray(np.asarray(weights, dtype=np.float64))
    if len(values):
        hist.FillN(len(values), values, weights)
    hist.SetDirectory(0)
    return hist


def draw(h_before, h_after, cfg, header, out_path: Path, era: str):
    canvas = ROOT.TCanvas("c_fsr", "", 800, 700)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)
    canvas.SetRightMargin(0.05)

    for hist, color in ((h_before, ROOT.kBlack), (h_after, ROOT.kRed + 1)):
        hist.SetLineColor(color)
        hist.SetLineWidth(LINE_WIDTH)
    h_before.GetXaxis().SetTitle(cfg["title"])
    h_before.GetYaxis().SetTitle("A.U.")
    format_axes(h_before)
    h_before.GetYaxis().SetTitleOffset(1.20)
    h_before.SetMaximum(1.5 * max(h_before.GetMaximum(), h_after.GetMaximum()))
    h_before.Draw("HIST")
    h_after.Draw("HIST SAME")

    leg = make_legend(0.17, 0.72, 0.55, 0.86)
    leg.AddEntry(h_before, "Without FSR recovery", "l")
    leg.AddEntry(h_after, "With FSR recovery", "l")
    leg.Draw()

    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.040)
    txt.DrawLatex(0.60, 0.84, header)
    labels = cms_label(year=era, cms_text="CMS Simulation")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_path))
    canvas.Close()
    _KEEP.append((leg, labels, txt))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ma", default="all")
    parser.add_argument("--year", default="all")
    parser.add_argument("--fsr-only", action="store_true",
                        help="restrict to events where a photon was actually added")
    parser.add_argument("--out", default=str(OUT_BASE))
    args = parser.parse_args()

    masses = MASSES if args.ma == "all" else [int(x) for x in args.ma.split(",")]
    eras = ERAS if args.year == "all" else args.year.split(",")
    out_base = Path(args.out)

    lines = ["ma  era             N(mu)  n_fsr>0[%]  dressed[%]  "
             "sigma_eff(mll) off->on   sigma_eff(mllgg) off->on"]
    for ma in masses:
        for era in eras:
            path = PARQUET_BASE / f"mA_M{ma}_{era}" / "merged_nominal.parquet"
            if not path.exists():
                print(f"  [skip] {path}")
                continue
            table = pq.read_table(path, columns=COLUMNS)
            data = {c: table[c].combine_chunks().to_numpy(zero_copy_only=False)
                    for c in COLUMNS}
            mu = data["z_mumu"] == 1
            if mu.sum() == 0:
                continue
            # "dressed" = the FSR photon was actually added to a lepton.
            # The two Z candidates are built through separate code paths, so an
            # exact != comparison fires on ~1e-3 GeV float noise for every event;
            # require a physically meaningful shift instead.
            dressed = mu & (np.abs(data["Z_mass"] - data["Z_noFSR_mass"]) > 0.01)
            sel = dressed if args.fsr_only else mu

            row = [f"{ma:>2}  {era:<14} {int(mu.sum()):>6}",
                   f"{100.0 * (data['n_fsr'][mu] > 0).mean():10.2f}",
                   f"{100.0 * dressed.sum() / mu.sum():11.3f}"]
            for key, cfg in VARIANTS.items():
                before = data[cfg["before"]][sel]
                after = data[cfg["after"]][sel]
                weights = data["weight_central"][sel]
                s_before, s_after = sigma_eff(before), sigma_eff(after)
                row.append(f"  {s_before:6.3f}->{s_after:6.3f}")
                header = ("m_{a} = %d GeV, #mu#mu channel" % ma
                          + (", FSR events" if args.fsr_only else ""))
                tag = "fsrOnly_" if args.fsr_only else ""
                draw(fill_hist(f"hb_{key}", before, weights, cfg),
                     fill_hist(f"ha_{key}", after, weights, cfg),
                     cfg, header,
                     out_base / era / f"fsr_{tag}{key}_mA_M{ma}_{era}.pdf", era)
            lines.append(" ".join(row))
            print(lines[-1])

    out_base.mkdir(parents=True, exist_ok=True)
    name = "fsr_summary_fsrOnly.txt" if args.fsr_only else "fsr_summary.txt"
    (out_base / name).write_text("\n".join(lines) + "\n")
    print("wrote", out_base / name)


if __name__ == "__main__":
    main()
