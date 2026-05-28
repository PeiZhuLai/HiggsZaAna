"""
PyROOT version of the Z->ee tag-and-probe MLPhoton regressed-mass spectrum
(AN-20-142 Fig 7/8 analogue), split by probe |eta| bin.

Inputs:
    MLNanoAOD ROOT files produced by Prod_MLNanoAOD_run3_mc.py on DY -> ee MiniAOD.

Output:
    Plot/plots/tnp_mlphoton_mass.{png,pdf,root}
"""

from __future__ import annotations

import argparse
import array
import glob
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)


# ---- physics constants & cuts (AN-20-142 §5.2) ----
TIGHT_ELE_CB = 4
PROBE_MIN_PT = 7.0
PROBE_ABS_ETA_MAX = 1.44
ZMASS_LO, ZMASS_HI = 60.0, 120.0
ML_MATCH_DR = 0.05

ETA_BINS = [
    (0.0, 0.5, "central"),
    (0.5, 1.0, "middle"),
    (1.0, 1.44, "forward"),
]
COLOR_BY_NAME = {
    "central": ROOT.kAzure + 2,
    "middle": ROOT.kOrange + 1,
    "forward": ROOT.kGreen + 2,
}


def list_inputs(path: str) -> list[str]:
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, "*.root")))
    else:
        files = sorted(glob.glob(path))
    if not files:
        raise SystemExit(f"no .root inputs found at {path!r}")
    print(f"[plot_tnp] inputs: {len(files)} files")
    return files


def build_chain(files: list[str]) -> ROOT.TChain:
    chain = ROOT.TChain("Events")
    for f in files:
        chain.Add(f)
    return chain


def build_rdf(chain: ROOT.TChain):
    """Build RDataFrame and per-bin TH1F via a single event loop."""
    df = ROOT.RDataFrame(chain)

    # Per-event helper: build all probe rows (1 row per tag-probe pair).
    # Uses pure numpy through Define so the loop runs in compiled code.
    # We use Define-with-vector to flatten then Take/Histo through Snapshot of derived arrays.
    # Simpler: do a single-pass scan via PyROOT and fill TH1F objects directly. This is
    # easier to read and the dataset is small (<300k events).

    return df


def fill_histograms(chain: ROOT.TChain) -> dict[str, ROOT.TH1F]:
    """Single-pass per-event loop filling per-|eta|-bin TH1F."""
    histos: dict[str, ROOT.TH1F] = {}
    for lo, hi, name in ETA_BINS:
        h = ROOT.TH1F(
            f"h_m_{name}",
            f"{name} |#eta| in ({lo:.2f},{hi:.2f});m_{{#Gamma}}^{{regressed}} on Z#rightarrowee probe [GeV];probes / 10 MeV",
            60, 0.0, 0.6,
        )
        h.SetLineColor(COLOR_BY_NAME[name])
        h.SetLineWidth(2)
        histos[name] = h

    n_total = chain.GetEntries()
    print(f"[plot_tnp] entries: {n_total}")
    n_probes = 0

    for i, event in enumerate(chain):
        if i and i % 50_000 == 0:
            print(f"[plot_tnp] processed {i}/{n_total} events, accepted probes so far: {n_probes}")

        n_e = event.nElectron
        if n_e < 2:
            continue

        # snapshot lists once per event so we don't pay attribute cost in inner loop
        e_pt = list(event.Electron_pt)
        e_eta = list(event.Electron_eta)
        e_phi = list(event.Electron_phi)
        e_q = list(event.Electron_charge)
        e_cb = list(event.Electron_cutBased)

        # MLPhoton list
        n_ml = event.nMLPhoton
        if n_ml == 0:
            continue
        ml_eta = list(event.MLPhoton_eta)
        ml_phi = list(event.MLPhoton_phi)
        ml_mass = list(event.MLPhoton_mass)

        for ti in range(n_e):
            if e_cb[ti] < TIGHT_ELE_CB:
                continue
            for pi in range(n_e):
                if pi == ti:
                    continue
                if e_q[ti] * e_q[pi] != -1:
                    continue
                if e_pt[pi] < PROBE_MIN_PT:
                    continue
                if abs(e_eta[pi]) >= PROBE_ABS_ETA_MAX:
                    continue
                # m_ll
                dphi = e_phi[ti] - e_phi[pi]
                # cosh(deta) - cos(dphi)
                import math
                cd = math.cosh(e_eta[ti] - e_eta[pi]) - math.cos(dphi)
                if cd <= 0:
                    continue
                mll = math.sqrt(2.0 * e_pt[ti] * e_pt[pi] * cd)
                if not (ZMASS_LO < mll < ZMASS_HI):
                    continue
                # nearest MLPhoton to probe
                best_dr = 1e9
                best_k = -1
                for k in range(n_ml):
                    deta = ml_eta[k] - e_eta[pi]
                    dphi_k = ml_phi[k] - e_phi[pi]
                    while dphi_k > math.pi:
                        dphi_k -= 2 * math.pi
                    while dphi_k < -math.pi:
                        dphi_k += 2 * math.pi
                    dr = math.sqrt(deta * deta + dphi_k * dphi_k)
                    if dr < best_dr:
                        best_dr = dr
                        best_k = k
                if best_dr > ML_MATCH_DR or best_k < 0:
                    continue
                aeta = abs(e_eta[pi])
                for lo, hi, name in ETA_BINS:
                    if lo <= aeta < hi:
                        histos[name].Fill(ml_mass[best_k])
                        n_probes += 1
                        break

    print(f"[plot_tnp] selected probes total: {n_probes}")
    for name, h in histos.items():
        print(f"  {name}: N={int(h.GetEntries())}, mean={h.GetMean():.4f}, RMS={h.GetRMS():.4f}")
    return histos


def draw(histos: dict[str, ROOT.TH1F], out_stub: str) -> None:
    c = ROOT.TCanvas("c_tnp_mlphoton", "Z#rightarrowee tag-and-probe MLPhoton mass", 900, 700)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetBottomMargin(0.13)

    # Determine y-max
    ymax = max(h.GetMaximum() for h in histos.values()) * 1.15

    leg = ROOT.TLegend(0.55, 0.70, 0.92, 0.90)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    first = True
    for lo, hi, name in ETA_BINS:
        h = histos[name]
        if first:
            h.SetMaximum(ymax)
            h.GetXaxis().SetTitleOffset(1.05)
            h.GetYaxis().SetTitleOffset(1.40)
            h.Draw("HIST")
            first = False
        else:
            h.Draw("HIST SAME")
        leg.AddEntry(h,
                     f"|#eta| {lo:.1f}-{hi:.2f} ({name}, N={int(h.GetEntries())})",
                     "l")

    leg.Draw()

    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.DrawLatex(0.15, 0.93, "Z #rightarrow ee tag-and-probe (Run3 2024 MC)")
    txt.SetTextSize(0.032)
    txt.DrawLatex(0.15, 0.88,
                  "tight tag (cutBased#geq4), opp.-charge probe p_{T}>7, 60<m_{ll}<120, |#eta_{probe}|<1.44")

    for ext in ("png", "pdf", "root"):
        path = f"{out_stub}.{ext}"
        c.SaveAs(path)
        print(f"[plot_tnp] saved {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "input", nargs="?",
        default="/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/DYto2E_MLL50",
        help="dir or glob of MLNanoAOD ROOT files",
    )
    ap.add_argument(
        "--output-stub", "-o",
        default="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/tnp_mlphoton_mass",
        help="output filename without extension",
    )
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.output_stub), exist_ok=True)

    files = list_inputs(args.input)
    chain = build_chain(files)
    histos = fill_histograms(chain)
    draw(histos, args.output_stub)


if __name__ == "__main__":
    main()
