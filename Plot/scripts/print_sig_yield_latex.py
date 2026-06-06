#!/usr/bin/env python3
"""
Read the signal-fit workspaces and print the expected signal yield per
(mA, year, channel), applying the 1/TEST_TREE_FRAC = 1/0.3 = 3.3333... rescaling
that brings the test-tree sumEntries back to the full inclusive expectation.

Per ws,  yield = TEST_TREE_INV * sumEntries(workspace_dataset)
which is algebraically equal to  xs * br * ea * lumi(pb^-1)  because the
ea spline already absorbs the same 1/TEST_TREE_FRAC factor.

Outputs:
  1) A LaTeX table with columns:
       mA | year (one per row) | ele | mu | ele+mu
     plus a TOTAL row per mA (sum over years) and an overall total.
  2) A second, compact LaTeX table with one row per mA: ele, mu, ele+mu.
"""
from __future__ import print_function

import argparse
import os
import sys
from collections import OrderedDict

import ROOT  # noqa: E402

ROOT.gErrorIgnoreLevel = ROOT.kError

TEST_TREE_FRAC = 0.3
TEST_TREE_INV  = 1.0 / TEST_TREE_FRAC

DEFAULT_BASE = (
    "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/"
    "flashggFinalFit/Signal"
)

DEFAULT_MA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
DEFAULT_YEARS    = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
DEFAULT_CHANNELS = ["ele", "mu"]


def _open_workspace(path):
    if not os.path.exists(path):
        return None, None
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return None, None
    ws = None
    for k in f.GetListOfKeys():
        obj = k.ReadObj()
        if obj.InheritsFrom("RooWorkspace"):
            ws = obj
            break
    if ws is None:
        f.Close()
        return None, None
    return f, ws


def _sum_entries(ws):
    """Sum entries of the nominal signal mass dataset stored in the workspace."""
    for ds in ws.allData():
        name = ds.GetName()
        if "sig_mass" in name and name.endswith("_cat0"):
            return float(ds.sumEntries())
    # fall-back: first dataset
    for ds in ws.allData():
        return float(ds.sumEntries())
    return 0.0


def gather_yields(base, mAs, years, channels):
    yields = OrderedDict()
    for ch in channels:
        for m in mAs:
            for yr in years:
                path = (
                    "{base}/outdir_{ch}/signalFit/output/"
                    "{m}_CMS-HGG_sigfit_{yr}_{ch}_Hm125.root"
                ).format(base=base, ch=ch, m=m, yr=yr)
                f, ws = _open_workspace(path)
                if ws is None:
                    yields[(m, ch, yr)] = None
                    continue
                yields[(m, ch, yr)] = TEST_TREE_INV * _sum_entries(ws)
                f.Close()
    return yields


def _fmt(v, width=10, decimals=3):
    if v is None:
        return "{:>{w}s}".format("--", w=width)
    return "{:>{w}.{d}f}".format(v, w=width, d=decimals)


def latex_table_per_year(yields, mAs, years, channels, out_path=None):
    col_spec = "c" + "c" * (len(years) * len(channels) + len(years) + 2)
    head_year = "".join([" & \\multicolumn{{{}}}{{c}}{{{}}}".format(len(channels) + 1, yr) for yr in years])
    head_ch = "".join([" & " + " & ".join(channels + ["e+\\mu"]) for _ in years])
    lines = [
        "% Expected signal yield = (1/0.3) * sumEntries(test tree workspace)",
        "\\begin{tabular}{" + col_spec + "}",
        "\\hline",
        "$m_a$ [GeV]" + head_year + " & \\multicolumn{3}{c}{Total} \\\\",
        " " + head_ch + " & ele & $\\mu$ & e+$\\mu$ \\\\",
        "\\hline",
    ]
    total_year = {(yr, ch): 0.0 for yr in years for ch in channels}
    total_year_lep = {yr: 0.0 for yr in years}
    grand_ele = grand_mu = 0.0
    for m in mAs:
        row = ["{}".format(m)]
        m_ele = m_mu = 0.0
        for yr in years:
            ele_v = yields.get((m, "ele", yr))
            mu_v  = yields.get((m, "mu",  yr))
            lep_v = (ele_v or 0.0) + (mu_v or 0.0) if (ele_v is not None or mu_v is not None) else None
            row.append("{:.3f}".format(ele_v) if ele_v is not None else "--")
            row.append("{:.3f}".format(mu_v) if mu_v is not None else "--")
            row.append("{:.3f}".format(lep_v) if lep_v is not None else "--")
            if ele_v is not None: total_year[(yr, "ele")] += ele_v; m_ele += ele_v
            if mu_v  is not None: total_year[(yr, "mu")]  += mu_v;  m_mu  += mu_v
            if lep_v is not None: total_year_lep[yr] += lep_v
        row += ["{:.3f}".format(m_ele), "{:.3f}".format(m_mu), "{:.3f}".format(m_ele + m_mu)]
        grand_ele += m_ele; grand_mu += m_mu
        lines.append(" & ".join(row) + " \\\\")
    lines.append("\\hline")
    tot = ["Total"]
    for yr in years:
        tot.append("{:.3f}".format(total_year[(yr, "ele")]))
        tot.append("{:.3f}".format(total_year[(yr, "mu")]))
        tot.append("{:.3f}".format(total_year_lep[yr]))
    tot += ["{:.3f}".format(grand_ele), "{:.3f}".format(grand_mu), "{:.3f}".format(grand_ele + grand_mu)]
    lines.append(" & ".join(tot) + " \\\\")
    lines += ["\\hline", "\\end{tabular}"]
    text = "\n".join(lines)
    if out_path:
        with open(out_path, "w") as fp:
            fp.write(text + "\n")
        print("[write] {}".format(out_path))
    return text


def latex_table_per_mA(yields, mAs, years, out_path=None):
    lines = [
        "% Expected signal yield per mA (years summed). yield = (1/0.3) * sumEntries.",
        "\\begin{table}[h]",
        "\\begin{center}",
        "\\topcaption{Expected signal yield ($\\mathrm{ggH}\\!\\rightarrow\\!Z\\alpha\\!\\rightarrow\\!2\\ell 2\\gamma$) "
        "in the electron and muon channels after the BDT category selection, "
        "with the assumed cross section 0.1pb and branching ratio 100\\% of $\\mathrm{ggH}\\!\\rightarrow\\!Z\\alpha\\!\\rightarrow\\!2\\ell 2\\gamma$, summed over Run3 eras (172.13~fb$^{-1}$).}",
        "\\label{tab:expected_sig_yield_perMA}",
        "\\small",
        "\\begin{tabular}{|c|c|c|c|}",
        "\\hline",
        "$m_a$ [GeV] & $N_{\\mathrm{sig}}^{ele}$ & $N_{\\mathrm{sig}}^{\\mu}$ & $N_{\\mathrm{sig}}^{e+\\mu}$ \\\\",
        "\\hline",
    ]
    grand_e = grand_m = 0.0
    for m in mAs:
        e = sum((yields.get((m, "ele", yr)) or 0.0) for yr in years)
        u = sum((yields.get((m, "mu",  yr)) or 0.0) for yr in years)
        grand_e += e; grand_m += u
        lines.append("{} & {:.3f} & {:.3f} & {:.3f} \\\\".format(m, e, u, e + u))
    lines.append("\\hline")
    lines.append("Total & {:.3f} & {:.3f} & {:.3f} \\\\".format(grand_e, grand_m, grand_e + grand_m))
    lines += ["\\hline", "\\end{tabular}", "\\end{center}", "\\end{table}"]
    text = "\n".join(lines)
    if out_path:
        with open(out_path, "w") as fp:
            fp.write(text + "\n")
        print("[write] {}".format(out_path))
    return text


def text_table(yields, mAs, years, channels):
    lines = []
    header = ["mA", "channel"] + years + ["sum_yr"]
    lines.append(" ".join("{:>11s}".format(h) for h in header))
    lines.append("-" * (12 * len(header)))
    per_total = {}
    for m in mAs:
        for ch in channels:
            row_vals = [yields.get((m, ch, yr)) for yr in years]
            s = sum(v for v in row_vals if v is not None)
            per_total[(m, ch)] = s
            cells = ["{:>11d}".format(m), "{:>11s}".format(ch)] + [_fmt(v, 11) for v in row_vals] + ["{:>11.3f}".format(s)]
            lines.append(" ".join(cells))
    lines.append("")
    lines.append("Per-mA totals (years summed)")
    lines.append("{:>5s} {:>12s} {:>12s} {:>12s}".format("mA", "ele", "mu", "ele+mu"))
    for m in mAs:
        e = per_total.get((m, "ele"), 0.0)
        u = per_total.get((m, "mu"),  0.0)
        lines.append("{:>5d} {:>12.3f} {:>12.3f} {:>12.3f}".format(m, e, u, e + u))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default=DEFAULT_BASE, help="flashggFinalFit/Signal base path")
    ap.add_argument("--mas", default=",".join(str(m) for m in DEFAULT_MA), help="comma-separated mA list")
    ap.add_argument("--years", default=",".join(DEFAULT_YEARS), help="comma-separated years")
    ap.add_argument("--channels", default=",".join(DEFAULT_CHANNELS), help="comma-separated channels")
    ap.add_argument("--latex-per-year", default=None, help="path to write per-year LaTeX table")
    ap.add_argument("--latex-per-mA",   default=None, help="path to write per-mA LaTeX table (years summed)")
    args = ap.parse_args()

    mAs       = [int(x) for x in args.mas.split(",") if x.strip()]
    years     = [x.strip() for x in args.years.split(",") if x.strip()]
    channels  = [x.strip() for x in args.channels.split(",") if x.strip()]

    yields = gather_yields(args.base, mAs, years, channels)
    print("\n=== Plain-text yields (expected events = sumEntries / TEST_TREE_FRAC) ===")
    print(text_table(yields, mAs, years, channels))

    print("\n=== LaTeX table (per mA, years summed) ===")
    print(latex_table_per_mA(yields, mAs, years, args.latex_per_mA))

    print("\n=== LaTeX table (per mA per year) ===")
    print(latex_table_per_year(yields, mAs, years, channels, args.latex_per_year))


if __name__ == "__main__":
    main()
