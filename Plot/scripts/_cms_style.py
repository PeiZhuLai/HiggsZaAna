"""
CMS-style helpers for PyROOT plots in this project.

- CMS Simulation label top-left
- "<lumi> fb^{-1} (13.6 TeV)" top-right
- Line width 3 (override hist line widths if you call apply_line_width)
- Legend default text size 0.045
- Axis title size 0.055 / label size 0.05 (matches the rest of the project)

Use:
    from _cms_style import cms_label, make_legend, COLORS, format_axes, LINE_WIDTH
"""

from __future__ import annotations

import os
import sys

# Make this importable from anywhere
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import ROOT

LINE_WIDTH = 3
LEGEND_TEXT_SIZE = 0.045
AXIS_TITLE_SIZE = 0.055
AXIS_LABEL_SIZE = 0.05

# Default lumi for the standard data-taking years used in this analysis.
LUMI_BY_YEAR = {
    "2022preEE": 7.99,
    "2022postEE": 26.68,
    "2023preBPix": 17.96,
    "2023postBPix": 9.68,
    "2024": 109.82,
    "run3": 172.13,
}


def _setup_globals():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)


_setup_globals()


def format_axes(h, *, title_size: float = AXIS_TITLE_SIZE,
                label_size: float = AXIS_LABEL_SIZE,
                x_title_offset: float = 1.05,
                y_title_offset: float = 1.30) -> None:
    """Common axis formatting; works for both TH1 and TH2."""
    for ax_name, off in [("GetXaxis", x_title_offset),
                          ("GetYaxis", y_title_offset),
                          ("GetZaxis", None)]:
        ax = getattr(h, ax_name)()
        ax.SetTitleSize(title_size)
        ax.SetLabelSize(label_size)
        if off is not None:
            ax.SetTitleOffset(off)


def apply_line_width(h, width: int = LINE_WIDTH) -> None:
    h.SetLineWidth(width)


def cms_label(*, lumi: float | str | None = None,
              year: str | None = None,
              cms_text: str = "CMS Simulation",
              energy: str = "13.6 TeV",
              x_left: float = 0.13,
              x_right: float = 0.95,
              y: float = 0.92,
              text_size: float = 0.045) -> tuple:
    """Draw CMS top-left + lumi top-right labels on the current pad.

    Returns the TLatex objects so the caller can keep refs alive.
    """
    if lumi is None and year is not None and year in LUMI_BY_YEAR:
        lumi = LUMI_BY_YEAR[year]

    # CMS-standard label: bold "CMS" + italic suffix (e.g. "CMS Simulation" ->
    # "#bf{CMS} #it{Simulation}"). Use font 42 and inline #bf/#it so the suffix
    # is italic, not bold. (matches plot_MVASigEffVmA.py)
    parts = cms_text.split(None, 1)
    if parts and parts[0] == "CMS":
        label_latex = "#bf{CMS}" + (" #it{%s}" % parts[1] if len(parts) > 1 else "")
    else:
        label_latex = "#bf{%s}" % cms_text
    txt_left = ROOT.TLatex()
    txt_left.SetNDC(True)
    txt_left.SetTextFont(42)
    txt_left.SetTextSize(text_size)
    txt_left.SetTextAlign(11)
    txt_left.DrawLatex(x_left, y, label_latex)

    txt_right = ROOT.TLatex()
    txt_right.SetNDC(True)
    txt_right.SetTextFont(42)
    txt_right.SetTextSize(text_size)
    txt_right.SetTextAlign(31)
    if lumi is not None:
        try:
            lumi_str = f"{float(lumi):.1f} fb^{{-1}} ({energy})"
        except (TypeError, ValueError):
            lumi_str = f"{lumi} fb^{{-1}} ({energy})"
    else:
        lumi_str = f"({energy})"
    txt_right.DrawLatex(x_right, y, lumi_str)
    return txt_left, txt_right


def make_legend(x1: float, y1: float, x2: float, y2: float,
                *, text_size: float = LEGEND_TEXT_SIZE,
                ncol: int = 1) -> "ROOT.TLegend":
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(text_size)
    leg.SetTextFont(42)
    if ncol > 1:
        leg.SetNColumns(ncol)
    return leg
