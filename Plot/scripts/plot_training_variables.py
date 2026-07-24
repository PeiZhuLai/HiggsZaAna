#!/usr/bin/env python3
"""
Plot the BDT training-variable signal vs. background distributions with PyROOT.

This is the PyROOT + project-convention rewrite of the matplotlib block in
    HZaMVA/scripts/run3_Za_BDT.py
(the `for hlf, xlabel_hlf, fn in zip(variables+['param']+mass_variables, ...)` loop).

The distributions are reproduced *exactly as the BDT sees them*:
  - same input ROOT files / `train` trees / H_m in (95, 180) selection,
  - same per-event weight `factor`,
  - the SAME iterative sideband reweight applied to BOTH signal and background
    (reuses HZaMVA/scripts/sideband_reweight.py, so it stays in sync with training),
  - density-normalised curves (matplotlib density=True -> ROOT Integral("width")=1).

Two figures are produced per variable:
  1. <var>_SigBkg.pdf   : all-mass signal (red) vs. background (blue)  [= original plot]
  2. <var>_perMass.pdf  : background + one signal curve per mA hypothesis (colour gradient)

Convention (via Plot/scripts/_cms_style.py):
  - top-left  "#bf{CMS} #it{Preliminary}", top-right "172.13 fb^{-1} (13.6 TeV)"
  - axis title size 0.055 / label size 0.05.

Env:  higgs-zg-plot_python  (PyROOT).  Run from anywhere.

Usage:
    python plot_training_variables.py                 # all variables, default OUTDIR
    OUTDIR=/some/dir python plot_training_variables.py
    python plot_training_variables.py pho1R9 var_dR_Za   # only these variables
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import uproot          # read trees WITHOUT ROOT/cling JIT (robust on conda-ROOT)
import ROOT            # plotting only (TH1/TCanvas/TLatex need no JIT)

# --- project style helpers (Plot/scripts/_cms_style.py) ----------------------
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from _cms_style import (  # noqa: E402
    AXIS_LABEL_SIZE,
    AXIS_TITLE_SIZE,
    LINE_WIDTH,
    LUMI_BY_YEAR,
    cms_label,
    format_axes,
    make_legend,
)

# --- reuse the training's sideband reweighter --------------------------------
_HZAMVA_SCRIPTS = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts"
sys.path.append(_HZAMVA_SCRIPTS)
from sideband_reweight import load_sideband_reweighter  # noqa: E402

ROOT.gROOT.SetBatch(True)

# ---------------------------------------------------------------------------
# Configuration (mirrors run3_Za_BDT.py)
# ---------------------------------------------------------------------------
FILE_PATH = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
BKG_DIR = "All_Bkg"
TRAIN_TREE = "train"          # BDT train sample, same as run3_Za_BDT.py
TRAIN_MASS_LOW = 95
TRAIN_MASS_HIGH = 180
SELECTION = "H_m>{} && H_m<{}".format(TRAIN_MASS_LOW, TRAIN_MASS_HIGH)
YEAR = "run3"
LUMI = "172.13"               # exact label used across the analysis

MASS_LIST = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 15., 20., 25., 30.]

# 14 base BDT variables + the one kept derived feature (pho_pt_asym).
BASE_VARIABLES = [
    "pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr",
    "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55", "pho2PIso_noCorr",
    "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2", "var_dR_g1Z",
    "var_PtaOverMh", "H_pt_oHm",
]
DERIVED_VARIABLES = ["pho_pt_asym"]          # computed, NOT a ROOT branch
VARIABLES = BASE_VARIABLES + DERIVED_VARIABLES
MASS_VARIABLES = ["ALP_m", "H_m"]

# Variables to plot, in the numbered order used by the dataVmc plots
# (Plot/plots/variables_dataVmc/plot_UL_mva/nominal): the filename number of a
# given variable matches that convention so plot sets line up across studies.
PLOT_VARIABLES = [
    "pho1Pt_oHm",              # 1
    "pho1R9",                  # 2
    "pho1IetaIeta55",          # 3
    "pho2Pt_oHm",              # 4
    "pho2R9",                  # 5
    "pho2IetaIeta55",          # 6
    "pho1PIso_noCorr",         # 7  (pho1ECALIso)
    "pho2PIso_noCorr",         # 8  (pho2ECALIso)
    "ALP_calculatedPhotonIso", # 9
    "var_dR_Za",               # 10
    "var_dR_g1g2",             # 11
    "var_dR_g1Z",              # 12
    "var_PtaOverMh",           # 13
    "H_pt_oHm",                # 14 (H_pt)
    "pho_pt_asym",             # 15
    "param",                   # 16
    "ALP_m",                   # 17
    "H_m",                     # 18
]
# 1-based number prefix per variable (no zero padding, matching the reference).
VAR_NUM = {v: i + 1 for i, v in enumerate(PLOT_VARIABLES)}

# ROOT/TLatex axis titles (PyROOT translation of the matplotlib latex labels).
XLABELS = {
    "pho1Pt_oHm": "#gamma_{Leading} p_{T} / m_{H}",
    "pho1R9": "#gamma_{Leading} R9",
    "pho1IetaIeta55": "#gamma_{Leading} #sigma_{i#etai#eta}^{5x5}",
    "pho1PIso_noCorr": "#gamma_{Leading} PF_{#gamma} Iso",
    "pho2Pt_oHm": "#gamma_{Subleading} p_{T} / m_{H}",
    "pho2R9": "#gamma_{Subleading} R9",
    "pho2IetaIeta55": "#gamma_{Subleading} #sigma_{i#etai#eta}^{5x5}",
    "pho2PIso_noCorr": "#gamma_{Subleading} PF_{#gamma} Iso",
    "ALP_calculatedPhotonIso": "#gamma#gamma Iso",
    "var_dR_Za": "#DeltaR(Z,a)",
    "var_dR_g1g2": "#DeltaR(#gamma,#gamma)",
    "var_dR_g1Z": "#DeltaR(#gamma_{Leading}, Z)",
    "var_PtaOverMh": "p_{T,a} / m_{H}",
    "H_pt_oHm": "p_{T,H} / m_{H}",
    "pho_pt_asym": "(p_{T,1}-p_{T,2})/(p_{T,1}+p_{T,2})",
    "param": "(m_{a} - m_{a,hyp}) / m_{H}",
    "ALP_m": "m_{a} [GeV]",
    "H_m": "m_{H} [GeV]",
}

X_LIMITS = {
    "pho1Pt_oHm": (0.0, 0.6),
    "pho1R9": (0.15, 1.3),
    "pho1IetaIeta55": (0.003, 0.04),
    "pho1PIso_noCorr": (-2, 25),
    "pho2Pt_oHm": (0.0, 0.5),
    "pho2R9": (0.15, 1.3),
    "pho2IetaIeta55": (0.003, 0.04),
    "pho2PIso_noCorr": (-2, 25),
    "ALP_calculatedPhotonIso": (1, 130),
    "var_dR_Za": (-0.6, 6),
    "var_dR_g1g2": (-0.2, 4.2),
    "var_dR_g1Z": (-0.6, 6),
    "var_PtaOverMh": (-0.1, 0.9),
    "H_pt_oHm": (0.0, 2.5),
    "pho_pt_asym": (0.0, 1.0),
    "param": (-1, 1),
    "ALP_m": (0, 55),
    "H_m": (110, 180),
}

BIN_SIZES = {
    "pho1Pt_oHm": 0.01,
    "pho1R9": 0.025,
    "pho1IetaIeta55": 0.001,
    "pho1PIso_noCorr": 1.,
    "pho2Pt_oHm": 0.01,
    "pho2R9": 0.025,
    "pho2IetaIeta55": 0.001,
    "pho2PIso_noCorr": 1.,
    "ALP_calculatedPhotonIso": 5,
    "var_dR_Za": 0.2,
    "var_dR_g1g2": 0.1,
    "var_dR_g1Z": 0.2,
    "var_PtaOverMh": 0.02,
    "H_pt_oHm": 0.04,
    "pho_pt_asym": 0.05,
    "param": 0.05,
    "ALP_m": 1,
    "H_m": 0.5,
}

OUTDIR = os.environ.get(
    "OUTDIR",
    "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/train_vars_sigVbkg",
)

# Pad margins.  The CMS/lumi labels start from the pad margins (see cms_label
# calls below): CMS label begins at the LEFT margin so it never runs off the
# frame, lumi ends at (1 - right margin).
LEFT_MARGIN = 0.15
RIGHT_MARGIN = 0.05
TOP_MARGIN = 0.08
BOTTOM_MARGIN = 0.16   # extra room so the x-axis title isn't clipped

# Per-mass colours (user-specified palette) and a 4-style line cycle so the
# 14 signal curves stay distinguishable in the per-mass overlay.
MASS_COLORS_HEX = {
    1: "#000000", 2: "#DC2626", 3: "#059669", 4: "#7C3AED", 5: "#84CC16",
    6: "#0F766E", 7: "#BE185D", 8: "#4B5563", 9: "#FACC15", 10: "#EF476F",
    15: "#991B1B", 20: "#6FFFE9", 25: "#78350F", 30: "#ff7f00",
}
# m_a=1 (black) and m_a=25 (brown) moved off blue so they don't clash with the
# blue background curve.  Two-style cycle keeps neighbouring curves separable.
LINE_STYLE_CYCLE = [1, 2]

# Branches to read from ROOT (pho_pt_asym & param are computed, not read).
_READ_BRANCHES = BASE_VARIABLES + ["ALP_m", "H_m", "factor", "Z_m", "event"]

SIDEBAND_REWEIGHTER = load_sideband_reweighter()


# ---------------------------------------------------------------------------
# Data loading (faithful to run3_Za_BDT.py)
# ---------------------------------------------------------------------------
def _read_frame(path):
    """Read the `train` tree of `path` with the training selection into a DataFrame.

    Uses uproot (pure python) so we never touch ROOT's cling JIT, which is
    unreliable in the conda-ROOT env; ROOT is used only for drawing.
    """
    with uproot.open(path) as f:
        tree = f[TRAIN_TREE]
        available = set(tree.keys())
        cols = [b for b in _READ_BRANCHES if b in available]
        arrays = tree.arrays(cols, library="np")
    frame = pd.DataFrame({c: np.asarray(arrays[c], dtype=float) for c in cols})
    # Training selection: H_m in (95, 180).
    frame = frame[(frame["H_m"] > TRAIN_MASS_LOW) & (frame["H_m"] < TRAIN_MASS_HIGH)]
    frame = frame.reset_index(drop=True)
    # pho_pt_asym: derived exactly as add_derived_features() in the training script.
    pt1 = frame["pho1Pt_oHm"].astype(float)
    pt2 = frame["pho2Pt_oHm"].astype(float)
    frame["pho_pt_asym"] = (pt1 - pt2) / (pt1 + pt2 + 1e-6)
    return frame


def load_background():
    path = "{}/{}/run3.root".format(FILE_PATH, BKG_DIR)
    print("[bkg] reading", path)
    frame = _read_frame(path)
    if SIDEBAND_REWEIGHTER is not None:
        # Background gets its param from the event-hash mass hypothesis, then the
        # same reweight as training multiplies `factor`.
        SIDEBAND_REWEIGHTER.ensure_param(frame, output_col="param", mass_output_col="mass")
        frame["factor_nominal"] = frame["factor"]
        SIDEBAND_REWEIGHTER.apply_to_dataframe(
            frame, base_weight_col="factor_nominal", output_weight_col="factor")
    else:
        SIDEBAND_REWEIGHTER_none_warn()
        frame["param"] = np.nan
    print("[bkg] entries={}  sum(factor)={:.6g}".format(len(frame), frame["factor"].sum()))
    return frame


def load_signal_by_mass():
    """Return {mass: DataFrame} plus the concatenated all-mass frame."""
    per_mass = {}
    for mass in MASS_LIST:
        path = "{}/mA_M{}/run3.root".format(FILE_PATH, int(mass))
        print("[sig] reading", path)
        frame = _read_frame(path)
        frame["mass"] = mass
        frame["param"] = (frame["ALP_m"] - frame["mass"]) / frame["H_m"]
        if SIDEBAND_REWEIGHTER is not None:
            # Same reweight as background, using the TRUE-mass param set above.
            frame["factor_nominal"] = frame["factor"]
            SIDEBAND_REWEIGHTER.apply_to_dataframe(
                frame, base_weight_col="factor_nominal", output_weight_col="factor")
        per_mass[mass] = frame
        print("    entries={}  sum(factor)={:.6g}".format(len(frame), frame["factor"].sum()))
    all_sig = pd.concat([per_mass[m] for m in MASS_LIST], ignore_index=True)
    return per_mass, all_sig


def SIDEBAND_REWEIGHTER_none_warn():
    print("WARNING: sideband reweight JSON not found; using nominal `factor` "
          "(distributions will differ from training).")


# ---------------------------------------------------------------------------
# Histogram helpers
# ---------------------------------------------------------------------------
def _binning(var):
    x_min, x_max = X_LIMITS[var]
    bin_size = BIN_SIZES.get(var, (x_max - x_min) / 40.)
    n_bins = int((x_max - x_min) / float(bin_size))
    return n_bins, float(x_min), float(x_max), bin_size


def _make_density_hist(name, values, weights, var):
    """Weighted, density-normalised TH1 (matplotlib density=True equivalent)."""
    n_bins, x_min, x_max, _ = _binning(var)
    h = ROOT.TH1D(name, "", n_bins, x_min, x_max)
    h.Sumw2()
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    good = np.isfinite(v) & np.isfinite(w)
    for xi, wi in zip(v[good], w[good]):
        h.Fill(xi, wi)
    integral = h.Integral("width")
    if integral != 0.0:
        h.Scale(1.0 / integral)
    return h


def _mass_color(mass):
    """ROOT colour index for a given mA hypothesis (user-specified palette)."""
    return ROOT.TColor.GetColor(MASS_COLORS_HEX[int(mass)])


def _style_frame(h, var, bin_size, y_max):
    h.GetXaxis().SetTitle(XLABELS[var])
    h.GetYaxis().SetTitle("A.U. / {:g}".format(bin_size))
    format_axes(h)
    h.SetMinimum(0.0)
    h.SetMaximum(1.35 * y_max)
    h.GetXaxis().SetNdivisions(505)


def _canvas():
    c = ROOT.TCanvas("c", "c", 800, 600)
    c.SetLeftMargin(LEFT_MARGIN)
    c.SetRightMargin(RIGHT_MARGIN)
    c.SetTopMargin(TOP_MARGIN)
    c.SetBottomMargin(BOTTOM_MARGIN)
    return c


def _draw_cms_label():
    """CMS Preliminary starting at the pad LEFT margin (so it never overruns the
    frame), lumi ending at the pad RIGHT edge."""
    cms_label(cms_text="CMS Preliminary", lumi=LUMI,
              x_left=LEFT_MARGIN, x_right=1.0 - RIGHT_MARGIN, y=0.94)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_sig_vs_bkg(var, sig_frame, bkg_frame):
    """All-mass signal (red) vs. background (blue) — the original combined plot."""
    _, _, _, bin_size = _binning(var)
    h_sig = _make_density_hist("h_sig_" + var, sig_frame[var], sig_frame["factor"], var)
    h_bkg = _make_density_hist("h_bkg_" + var, bkg_frame[var], bkg_frame["factor"], var)

    h_sig.SetLineColor(ROOT.kRed + 1)
    h_sig.SetLineStyle(1)
    h_sig.SetLineWidth(3)
    h_bkg.SetLineColor(ROOT.kBlue + 1)
    h_bkg.SetLineStyle(9)          # dash-dot-like, matching the matplotlib '-.'
    h_bkg.SetLineWidth(LINE_WIDTH)

    y_max = max(h_sig.GetMaximum(), h_bkg.GetMaximum())

    c = _canvas()
    _style_frame(h_sig, var, bin_size, y_max)
    h_sig.Draw("HIST")
    h_bkg.Draw("HIST SAME")

    leg = make_legend(0.62, 0.74, 0.93, 0.88)
    leg.AddEntry(h_sig, "Signal", "l")
    leg.AddEntry(h_bkg, "Background", "l")
    leg.Draw()

    _draw_cms_label()
    c.RedrawAxis()
    out = os.path.join(OUTDIR, "{}_{}_SigBkg.pdf".format(VAR_NUM[var], var))
    c.SaveAs(out)
    print("[plot]", out)


def plot_per_mass(var, per_mass, bkg_frame):
    """Background + one signal density curve per mA hypothesis (colour gradient)."""
    _, _, _, bin_size = _binning(var)

    h_bkg = _make_density_hist("hpm_bkg_" + var, bkg_frame[var], bkg_frame["factor"], var)
    h_bkg.SetLineColor(ROOT.kBlue + 1)   # match the background colour in the SigBkg plot
    h_bkg.SetLineStyle(1)
    h_bkg.SetLineWidth(LINE_WIDTH + 1)

    masses = MASS_LIST
    sig_hists = []
    for i, mass in enumerate(masses):
        frame = per_mass[mass]
        h = _make_density_hist("hpm_sig_{}_{}".format(var, int(mass)),
                               frame[var], frame["factor"], var)
        h.SetLineColor(_mass_color(mass))
        h.SetLineWidth(3)
        h.SetLineStyle(LINE_STYLE_CYCLE[i % len(LINE_STYLE_CYCLE)])
        sig_hists.append((mass, h))

    y_max = max([h_bkg.GetMaximum()] + [h.GetMaximum() for _, h in sig_hists])

    c = _canvas()
    _style_frame(h_bkg, var, bin_size, y_max)
    h_bkg.SetMaximum(1.55 * y_max)   # extra headroom for the right-side legend
    h_bkg.Draw("HIST")
    for _, h in sig_hists:
        h.Draw("HIST SAME")
    h_bkg.Draw("HIST SAME")   # keep background on top for reference

    # 15 entries -> compact 3-column legend, shifted to the right so it clears
    # left-peaking variables (e.g. ALP_calculatedPhotonIso).
    leg = make_legend(0.42, 0.72, 0.95, 0.90, text_size=0.028, ncol=3)
    leg.AddEntry(h_bkg, "Bkg", "l")
    for mass, h in sig_hists:
        leg.AddEntry(h, "m_{{a}}={:g} GeV".format(mass), "l")
    leg.Draw()

    _draw_cms_label()
    c.RedrawAxis()
    out = os.path.join(OUTDIR, "{}_{}_perMass.pdf".format(VAR_NUM[var], var))
    c.SaveAs(out)
    print("[plot]", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    requested = [v for v in sys.argv[1:] if not v.startswith("-")]
    variables = requested if requested else PLOT_VARIABLES
    for v in variables:
        if v not in X_LIMITS:
            print("[skip] {}: no x_limits/label defined".format(v))
    variables = [v for v in variables if v in X_LIMITS]

    os.makedirs(OUTDIR, exist_ok=True)
    print("OUTDIR =", OUTDIR)
    if SIDEBAND_REWEIGHTER is not None:
        print("Sideband reweight:", SIDEBAND_REWEIGHTER.source_path)

    bkg_frame = load_background()
    per_mass, all_sig = load_signal_by_mass()

    for v in variables:
        plot_sig_vs_bkg(v, all_sig, bkg_frame)
        plot_per_mass(v, per_mass, bkg_frame)

    print("Done. {} variables -> {} PDFs in {}".format(
        len(variables), 2 * len(variables), OUTDIR))


if __name__ == "__main__":
    main()
