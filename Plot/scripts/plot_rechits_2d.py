"""
2D EB RecHit (i#eta, i#phi) energy heatmaps from a pre-dumped ROOT TTree
(`EBRecHitDumper/rechits`), all in PyROOT.

Two modes:

  --mode event   per-event 30x30 crystal window centered on the highest-energy
                 EB RecHit. Saves the first N events to PNG/PDF (one file each).

  --mode aggregate
                 sum across N events, each centered on its own seed crystal.
                 Equivalent to the AN-20-142 Fig 2 average shower template.

Input: a ROOT file produced by dump_rechits_cfg.py, stored under
       /eos/home-p/pelai/HZa/root_rechit/
       (TTree path: 'dumper/rechits' under TFileService default folder).

Output: PNG/PDF only (no ROOT) under
        /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/rechits_2d/.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from array import array

import ROOT

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

N_CONTOURS = 255


def set_jet_palette(n_contours: int = N_CONTOURS) -> None:
    """Install the classic 'jet' rainbow palette used in the EGM-20-001 /
    end-to-end mass-regression figures (HIGApp_2021MAR04_v2.pdf, slide 6):
    dark-navy (zero) -> blue -> cyan -> green -> yellow -> red (-> dark red).

    A custom gradient table is used instead of a ROOT built-in so the look is
    reproducible across ROOT versions. Empty (content==0) bins stay white in
    COLZ, matching the reference figure.
    """
    stops = array("d", [0.00, 0.125, 0.375, 0.625, 0.875, 1.00])
    red = array("d", [0.00, 0.00, 0.00, 1.00, 1.00, 0.50])
    green = array("d", [0.00, 0.00, 1.00, 1.00, 0.00, 0.00])
    blue = array("d", [0.50, 1.00, 1.00, 0.00, 0.00, 0.00])
    ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, n_contours)
    ROOT.gStyle.SetNumberContours(n_contours)


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
set_jet_palette()

DEFAULT_INPUT = "/eos/home-p/pelai/HZa/root_rechit/mA_M1_rechits_100ev_numEvent100.root"
DEFAULT_OUT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/rechits_2d"
# 2024 integrated luminosity (HiggsDNA constants.py LUMI["2024"]).
DEFAULT_LUMI = "109.82 fb^{-1}"

ISIZE = 30  # window size (Cluster::isize in the MLPhoton producer)


def _make_th2(name: str, title: str) -> ROOT.TH2F:
    h = ROOT.TH2F(
        name,
        title + ";i#eta - i#eta_{seed};i#phi - i#phi_{seed};Energy [GeV]",
        ISIZE, -ISIZE / 2.0, ISIZE / 2.0,
        ISIZE, -ISIZE / 2.0, ISIZE / 2.0,
    )
    for ax in (h.GetXaxis(), h.GetYaxis(), h.GetZaxis()):
        ax.SetTitleSize(0.055)
        ax.SetLabelSize(0.05)
    h.GetXaxis().SetTitleOffset(1.00)
    h.GetYaxis().SetTitleOffset(1.10)
    h.GetZaxis().SetTitleOffset(1.20)
    h.SetContour(N_CONTOURS)
    return h


def _find_seed(ieta, iphi, energy):
    seed_e = -1.0
    seed_ieta = 0
    seed_iphi = 0
    n = len(energy)
    for k in range(n):
        e = energy[k]
        if e <= 0:
            continue
        if e > seed_e:
            seed_e = e
            seed_ieta = int(ieta[k])
            seed_iphi = int(iphi[k])
    return seed_ieta, seed_iphi, seed_e


def _fill_into(h: ROOT.TH2F, ieta, iphi, energy, seed_ieta: int, seed_iphi: int):
    half = ISIZE / 2.0
    n = len(energy)
    for k in range(n):
        e = energy[k]
        if e <= 0:
            continue
        dx = int(ieta[k]) - seed_ieta
        dy = int(iphi[k]) - seed_iphi
        if dy > 180:
            dy -= 360
        elif dy < -180:
            dy += 360
        if abs(dx) >= half or abs(dy) >= half:
            continue
        h.Fill(dx, dy, e)


def _mass_label_from_path(path: str) -> str | None:
    """Parse the nominal a-boson mass (GeV) from an 'mA_M<...>' tag in a path.

    Examples: mA_M0p1 -> "0.1", mA_M0p5 -> "0.5", mA_M1 -> "1.0".
    Returns None if no mass tag is found.
    """
    import re

    m = re.search(r"mA_M([0-9]+(?:p[0-9]+)?)", path)
    if not m:
        return None
    val = m.group(1).replace("p", ".")
    if "." not in val:
        val += ".0"
    return val


def _draw_save(h: ROOT.TH2F, out_stub: str, *, logz: bool = True,
               mass_label: str | None = None, lumi_label: str = DEFAULT_LUMI):
    c = ROOT.TCanvas(f"c_{h.GetName()}", h.GetTitle(), 800, 660)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.18)
    c.SetTopMargin(0.10)
    c.SetBottomMargin(0.14)
    if logz and h.GetMaximum() > 0:
        c.SetLogz(True)
        h.SetMinimum(1e-3)
    h.Draw("COLZ")

    # Above the frame: "CMS Simulation" (top-left) + lumi/energy (top-right),
    # following the standard CMS header convention.
    # CMS-standard: bold "CMS" + italic suffix (font 42 + #bf/#it inline).
    txt_left = ROOT.TLatex()
    txt_left.SetNDC(True)
    txt_left.SetTextFont(42)
    txt_left.SetTextSize(0.045)
    txt_left.SetTextAlign(11)
    txt_left.DrawLatex(0.13, 0.92, "#bf{CMS} #it{Simulation}")

    lumi_text = "%s (13.6 TeV)" % lumi_label if lumi_label else "(13.6 TeV)"
    txt_lumi = ROOT.TLatex()
    txt_lumi.SetNDC(True)
    txt_lumi.SetTextFont(42)
    txt_lumi.SetTextSize(0.045)
    txt_lumi.SetTextAlign(31)
    txt_lumi.DrawLatex(0.82, 0.92, lumi_text)

    # Inside the frame, top-right corner: the a-boson mass label.
    txt_right = None
    if mass_label is not None:
        txt_right = ROOT.TLatex()
        txt_right.SetNDC(True)
        txt_right.SetTextFont(42)
        txt_right.SetTextSize(0.045)
        txt_right.SetTextAlign(33)
        txt_right.DrawLatex(0.80, 0.87, "m_{a} = %s GeV" % mass_label)

    os.makedirs(os.path.dirname(out_stub), exist_ok=True)
    for ext in ("png", "pdf"):
        path = f"{out_stub}.{ext}"
        c.SaveAs(path)
        print(f"[plot_rechits_2d] saved {path}")


def _open_tree(inputs: list[str]) -> ROOT.TChain:
    chain = ROOT.TChain("dumper/rechits")
    n = 0
    for f in inputs:
        n_added = chain.Add(f)
        n += n_added
    if n == 0 or chain.GetEntries() == 0:
        raise SystemExit(f"no entries in 'dumper/rechits' from {inputs!r}")
    print(f"[plot_rechits_2d] entries: {chain.GetEntries()} from {n} file(s)")
    return chain


def _entry_arrays(tree, i):
    tree.GetEntry(i)
    return (
        list(tree.rh_ieta),
        list(tree.rh_iphi),
        list(tree.rh_energy),
        int(tree.event),
        int(tree.run),
        int(tree.lumi),
        float(tree.genALP_mass) if tree.genALP_found else float("nan"),
        float(tree.genALP_pt) if tree.genALP_found else float("nan"),
        float(tree.genALP_eta) if tree.genALP_found else float("nan"),
    )


def mode_event(tree: ROOT.TChain, out_dir: str, max_events: int,
               mass_label: str | None = None, lumi_label: str = DEFAULT_LUMI):
    n_total = tree.GetEntries()
    n_draw = min(max_events, n_total)
    print(f"[plot_rechits_2d] event mode, drawing {n_draw} events")
    for i in range(n_draw):
        ieta, iphi, energy, evt, run, lumi, gen_mass, gen_pt, gen_eta = _entry_arrays(tree, i)
        sieta, siphi, se = _find_seed(ieta, iphi, energy)
        if se <= 0:
            continue
        title = (
            f"ev #{i} run/lumi/evt {run}/{lumi}/{evt} | "
            f"seed E={se:.2f} GeV @ i#eta={sieta}, i#phi={siphi} | "
            f"gen m_a={gen_mass:.3f}, p_T={gen_pt:.1f}, #eta={gen_eta:.2f}"
        )
        h = _make_th2(f"h_ev_{i:03d}", title)
        _fill_into(h, ieta, iphi, energy, sieta, siphi)
        stub = os.path.join(out_dir, "event", f"event_{i:04d}")
        _draw_save(h, stub, logz=True, mass_label=mass_label, lumi_label=lumi_label)


def mode_aggregate(tree: ROOT.TChain, out_dir: str, max_events: int,
                   mass_label: str | None = None, lumi_label: str = DEFAULT_LUMI):
    n_total = tree.GetEntries()
    n_use = min(max_events, n_total)
    h = _make_th2("h_aggregate", "Average EB RecHit energy (seed-centered)")
    n_used = 0
    for i in range(n_use):
        ieta, iphi, energy, *_ = _entry_arrays(tree, i)
        sieta, siphi, se = _find_seed(ieta, iphi, energy)
        if se <= 0:
            continue
        _fill_into(h, ieta, iphi, energy, sieta, siphi)
        n_used += 1
    h.SetTitle(f"Sum of EB RecHit energy on {n_used} events, seed-centered")
    stub = os.path.join(out_dir, "aggregate", f"aggregate_{n_used}ev")
    _draw_save(h, stub, logz=True, mass_label=mass_label, lumi_label=lumi_label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "input", nargs="*", default=[DEFAULT_INPUT],
        help="ROOT file(s) from dump_rechits_cfg.py",
    )
    ap.add_argument("--mode", choices=("event", "aggregate"), default="event")
    ap.add_argument("--max-events", "-n", type=int, default=10)
    ap.add_argument("--out-dir", "-o", default=DEFAULT_OUT_DIR)
    ap.add_argument(
        "--mass", default=None,
        help="a-boson mass (GeV) for the in-frame top-right label; "
             "auto-derived from the 'mA_M<...>' tag in the input/out-dir if omitted",
    )
    ap.add_argument(
        "--lumi", default=DEFAULT_LUMI,
        help="luminosity string for the top-right CMS header (before '(13.6 TeV)'); "
             f"default '{DEFAULT_LUMI}' (2024). Pass '' to show only the energy.",
    )
    args = ap.parse_args()

    # Top-right m_a label: explicit --mass wins, else parse the mA_M<...> tag
    # from the input file name, else from the output directory.
    mass_label = args.mass
    if mass_label is None:
        for src in list(args.input) + [args.out_dir]:
            mass_label = _mass_label_from_path(src)
            if mass_label is not None:
                break
    print(f"[plot_rechits_2d] m_a label: {mass_label}")

    tree = _open_tree(args.input)
    if args.mode == "event":
        mode_event(tree, args.out_dir, args.max_events, mass_label, args.lumi)
    else:
        mode_aggregate(tree, args.out_dir, args.max_events, mass_label, args.lumi)


if __name__ == "__main__":
    main()
