"""
2D ECAL RecHit (i#eta, i#phi) energy maps from MiniAOD events, in PyROOT.

Reproduces the style of AN-20-142 Fig 2 (right column): the ECAL shower
pattern around a merged a -> gamma gamma candidate, drawn as a 2D heatmap of
energy on a 30x30 crystal window centered on the highest-energy crystal.

Two modes:

  --mode event   single-event display: draws N event-by-event maps in a
                 multi-page PDF + single PNGs. Useful for sanity checks
                 and for picking interesting events.

  --mode aggregate
                 sum of energy maps across many events, centered on the
                 leading-energy crystal of each event. Shows the average
                 shower shape (akin to a "template").

Inputs: any MiniAOD file (PAT or RECO), reading
        reducedEgamma:reducedEBRecHits (default tag).

CMSSW env required (we read EcalRecHit + EBDetId):
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=el9_amd64_gcc12
    cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src
    eval `scramv1 runtime -sh`

Output:
    Plot/plots/rechits_2d/{mode}/{stub}.{png,pdf,root}
"""

from __future__ import annotations

import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(64)

# Default MiniAODv6 inputs (Run3 2024) — small subset is enough.
DEFAULT_INPUT = (
    "root://cms-xrd-global.cern.ch//store/mc/RunIII2024Summer24MiniAODv6/"
    "HZa-Zto2L-ato2G_Par-M-0p4_TuneCP5_13p6TeV_madgraph-pythia8/MINIAODSIM/"
    "150X_mcRun3_2024_realistic_v2-v2/2560000/"
    "d8063b96-2b41-4a53-9236-c6417c2de7e9.root"
)
DEFAULT_OUT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/rechits_2d"

EB_INPUT_TAG = ("reducedEgamma", "reducedEBRecHits")
ISIZE = 30  # window half-size used by the MLPhoton classifier (Cluster::isize)


def _load_cmssw_libraries():
    """Make sure the EDM data formats we need are loaded under PyROOT."""
    needed = [
        "FWCoreFWLite",
        "DataFormatsFWLite",
        "DataFormatsEcalRecHit",
        "DataFormatsEcalDetId",
    ]
    for lib in needed:
        if ROOT.gSystem.Load(f"lib{lib}") < 0:
            raise SystemExit(
                f"failed to load CMSSW lib {lib}; did you cmsenv inside CMSSW_15_0_14?"
            )
    ROOT.FWLiteEnabler.enable()


def _open_minilike(files: list[str]):
    """Open MiniAOD inputs as an Events handle (FWLite-style)."""
    from DataFormats.FWLite import Events, Handle  # type: ignore
    events = Events(files)
    handle = Handle("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >")
    return events, handle


def _make_th2(name: str, title: str) -> ROOT.TH2F:
    h = ROOT.TH2F(
        name, title,
        ISIZE, -ISIZE / 2, ISIZE / 2,   # x: i#eta offset
        ISIZE, -ISIZE / 2, ISIZE / 2,   # y: i#phi offset
    )
    h.GetXaxis().SetTitle("i#eta - i#eta_{seed}")
    h.GetYaxis().SetTitle("i#phi - i#phi_{seed}")
    h.GetZaxis().SetTitle("Energy [GeV]")
    h.SetContour(64)
    return h


def _find_seed(rechits) -> tuple[int, int, float]:
    """Return (ieta_seed, iphi_seed, max_energy)."""
    seed_e = -1.0
    seed_ieta = 0
    seed_iphi = 0
    for rh in rechits:
        e = rh.energy()
        if e <= 0.0:
            continue
        # EBDetId
        det = ROOT.EBDetId(rh.detid())
        if e > seed_e:
            seed_e = e
            seed_ieta = det.ieta()
            seed_iphi = det.iphi()
    return seed_ieta, seed_iphi, seed_e


def _fill_into(h: ROOT.TH2F, rechits, seed_ieta: int, seed_iphi: int, *, add: bool = False):
    """Fill h with (i#eta - seed, i#phi - seed, energy)."""
    half = ISIZE / 2
    for rh in rechits:
        e = rh.energy()
        if e <= 0.0:
            continue
        det = ROOT.EBDetId(rh.detid())
        dx = det.ieta() - seed_ieta
        dy = det.iphi() - seed_iphi
        # wrap φ
        if dy > 180:
            dy -= 360
        elif dy < -180:
            dy += 360
        if abs(dx) >= half or abs(dy) >= half:
            continue
        if add:
            h.Fill(dx, dy, e)
        else:
            h.Fill(dx, dy, e)


def _draw_save(h: ROOT.TH2F, out_stub: str, *, logz: bool = True):
    c = ROOT.TCanvas(f"c_{h.GetName()}", h.GetTitle(), 720, 640)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.15)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    if logz and h.GetMaximum() > 0:
        c.SetLogz(True)
        # Avoid zero-bins killing the log scale
        h.SetMinimum(1e-3)
    h.Draw("COLZ")

    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)
    txt.DrawLatex(0.12, 0.94, h.GetTitle())

    os.makedirs(os.path.dirname(out_stub), exist_ok=True)
    for ext in ("png", "pdf", "root"):
        c.SaveAs(f"{out_stub}.{ext}")
        print(f"[plot_rechits_2d] saved {out_stub}.{ext}")


def mode_event(events, handle, out_dir: str, max_events: int):
    written = 0
    for i, ev in enumerate(events):
        if written >= max_events:
            break
        ev.getByLabel(EB_INPUT_TAG[0], EB_INPUT_TAG[1], handle)
        if not handle.isValid():
            continue
        rh = handle.product()
        seed_ieta, seed_iphi, seed_e = _find_seed(rh)
        if seed_e <= 0.0:
            continue
        run = ev.eventAuxiliary().run()
        lumi = ev.eventAuxiliary().luminosityBlock()
        evt = ev.eventAuxiliary().event()
        title = (f"ev #{i} run/lumi/evt {run}/{lumi}/{evt}: "
                 f"seed E={seed_e:.2f} GeV @ i#eta={seed_ieta}, i#phi={seed_iphi}")
        h = _make_th2(f"h_ev_{i}", title)
        _fill_into(h, rh, seed_ieta, seed_iphi)
        stub = os.path.join(out_dir, "event", f"event_{i:04d}")
        _draw_save(h, stub, logz=True)
        written += 1
    if written == 0:
        print("[plot_rechits_2d] no qualifying events found")


def mode_aggregate(events, handle, out_dir: str, max_events: int):
    h = _make_th2("h_aggregate", f"Sum of EB RecHit energy maps (N#leq{max_events})")
    n_used = 0
    for i, ev in enumerate(events):
        if n_used >= max_events:
            break
        ev.getByLabel(EB_INPUT_TAG[0], EB_INPUT_TAG[1], handle)
        if not handle.isValid():
            continue
        rh = handle.product()
        seed_ieta, seed_iphi, seed_e = _find_seed(rh)
        if seed_e <= 0.0:
            continue
        _fill_into(h, rh, seed_ieta, seed_iphi, add=True)
        n_used += 1
    h.SetTitle(f"Sum of EB RecHit energy on {n_used} events centered on seed crystal")
    stub = os.path.join(out_dir, "aggregate", f"aggregate_{n_used}ev")
    _draw_save(h, stub, logz=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "input", nargs="*", default=[DEFAULT_INPUT],
        help="MiniAOD input file(s) — root:// or local path. Default: one mA=0.4 file.",
    )
    ap.add_argument(
        "--mode", choices=("event", "aggregate"), default="event",
        help="event: per-event 2D map, default. aggregate: sum over events.",
    )
    ap.add_argument(
        "--max-events", "-n", type=int, default=10,
        help="cap on events processed (default 10).",
    )
    ap.add_argument(
        "--out-dir", "-o", default=DEFAULT_OUT_DIR,
        help="output directory under Plot/plots/.",
    )
    args = ap.parse_args()

    _load_cmssw_libraries()
    events, handle = _open_minilike(args.input)

    print(f"[plot_rechits_2d] inputs: {args.input}")
    print(f"[plot_rechits_2d] mode={args.mode}, max_events={args.max_events}")
    print(f"[plot_rechits_2d] out_dir={args.out_dir}")

    if args.mode == "event":
        mode_event(events, handle, args.out_dir, args.max_events)
    else:
        mode_aggregate(events, handle, args.out_dir, args.max_events)


if __name__ == "__main__":
    main()
