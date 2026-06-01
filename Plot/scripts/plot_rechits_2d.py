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

import ROOT

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from _cms_style import cms_label

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(64)

DEFAULT_INPUT = "/eos/home-p/pelai/HZa/root_rechit/mA_M1_rechits_100ev_numEvent100.root"
DEFAULT_OUT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/rechits_2d"

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
    h.SetContour(64)
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


def _draw_save(h: ROOT.TH2F, out_stub: str, *, logz: bool = True):
    c = ROOT.TCanvas(f"c_{h.GetName()}", h.GetTitle(), 800, 660)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.18)
    c.SetTopMargin(0.10)
    c.SetBottomMargin(0.14)
    if logz and h.GetMaximum() > 0:
        c.SetLogz(True)
        h.SetMinimum(1e-3)
    h.Draw("COLZ")

    cms_lab = cms_label(year="2024", x_left=0.13, x_right=0.82)

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


def mode_event(tree: ROOT.TChain, out_dir: str, max_events: int):
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
        _draw_save(h, stub, logz=True)


def mode_aggregate(tree: ROOT.TChain, out_dir: str, max_events: int):
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
    _draw_save(h, stub, logz=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "input", nargs="*", default=[DEFAULT_INPUT],
        help="ROOT file(s) from dump_rechits_cfg.py",
    )
    ap.add_argument("--mode", choices=("event", "aggregate"), default="event")
    ap.add_argument("--max-events", "-n", type=int, default=10)
    ap.add_argument("--out-dir", "-o", default=DEFAULT_OUT_DIR)
    args = ap.parse_args()

    tree = _open_tree(args.input)
    if args.mode == "event":
        mode_event(tree, args.out_dir, args.max_events)
    else:
        mode_aggregate(tree, args.out_dir, args.max_events)


if __name__ == "__main__":
    main()
