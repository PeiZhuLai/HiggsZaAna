#!/usr/bin/env python3
"""L3 review S1 plots (PyROOT, house style via HZgamma plot/scripts/hzg_style.py, read-only import).
Run with LCG_104:  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
"""
import sys
import numpy as np, pandas as pd
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDFs; use LCG_104")
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/plot/scripts")
from hzg_style import (canvas_single, style_axes, cms_header, panel_label, house_legend,  # noqa
                       color, PETROFF, MARKERS, save, _keep, below_frame)
from lumi_constants import load_lumi  # noqa

ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
S = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927"
OUT = f"{S}/plots"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
LUMI = load_lumi()
LUMI_RUN3 = sum(LUMI[e] for e in ERAS)
ALTS = [("zmhsb", "m_{ll} 50-80, 100-120 & m_{ll#gamma#gamma} SB"),
        ("zmnarrowhsb", "m_{ll} 80-86, 96-102 & m_{ll#gamma#gamma} SB"),
        ("zm", "m_{ll} 50-80, 100-120 (incl. SR)"),
        ("zmnarrow", "m_{ll} 80-86, 96-102 (incl. SR)")]

summ = pd.read_csv(f"{S}/results/s1_altsideband_run3_summary.csv")
full = pd.read_csv(f"{S}/results/s1_altsideband_yields.csv")


def frame(name, ylo, yhi, ytitle):
    c = canvas_single(name)
    h = ROOT.TH1F(name + "_f", "", 100, 0, 31.5)
    h.SetMinimum(ylo); h.SetMaximum(yhi)
    style_axes(h, c, "m_{a} [GeV]", ytitle)
    h.Draw("axis"); _keep.append(h)
    ln = ROOT.TLine(0, 0, 31.5, 0); ln.SetLineStyle(2); ln.SetLineColor(ROOT.kGray + 2); ln.Draw(); _keep.append(ln)
    return c


def graph(x, y, ey, i, style_idx=None):
    g = ROOT.TGraphErrors(len(x), np.array(x, float), np.array(y, float), np.zeros(len(x)), np.array(ey, float))
    k = i if style_idx is None else style_idx
    g.SetMarkerStyle(MARKERS[k % len(MARKERS)]); g.SetMarkerSize(1.4)
    g.SetMarkerColor(color(PETROFF[k])); g.SetLineColor(color(PETROFF[k])); g.SetLineWidth(3)
    _keep.append(g)
    return g


# 1) Run-3 relative yield difference vs mA, with the current nuisance as a band
a = summ[summ.ch == "all"].sort_values("mA")
m = a.mA.to_numpy(float)
c = frame("c_run3", -14, 34, "#DeltaN_{sig}/N_{sig} [%]")
band = ROOT.TGraphAsymmErrors(len(m), m, np.zeros(len(m)), np.full(len(m), 0.4), np.full(len(m), 0.4),
                              -100 * a.cur_dn.to_numpy(), 100 * a.cur_up.to_numpy())
band.SetFillColor(ROOT.kGray); band.SetFillStyle(1001); band.SetLineColor(ROOT.kGray)
band.Draw("2 same"); _keep.append(band)
entries, objs = [], []
for i, (k, lab) in enumerate(ALTS):
    g = graph(m + (i - 1.5) * 0.18, 100 * a[f"dy_{k}"], 100 * a[f"dyerr_{k}"], i)
    if k in ("zm", "zmnarrow"):
        g.SetMarkerStyle([24, 25][i - 2])
    g.Draw("P same"); entries.append(lab); objs.append(g)
entries.append("current CMS_hza_mva_reweight"); objs.append(band)
cms_header(c, LUMI_RUN3, "13.6", extra="Preliminary")
panel_label(c, "Signal, Run 3, ee+#mu#mu", notes=["alt. sideband vs nominal, after BDT WP"])
leg = house_legend(c, entries, size=0.034, top=below_frame(0.06 + 0.055 + 0.035, c))
for o, t in zip(objs, entries):
    leg.AddEntry(o, t, "f" if o is band else "pe")
leg.Draw()
c.SetGrid(0, 0)
save(c, f"{OUT}/s1_dyield_vs_mA_run3")

# 2) per era for the two sideband-excluding alternatives
for k, lab in ALTS[:2]:
    lo, hi = (-18, 16)
    c = frame("c_era_" + k, lo, hi, "#DeltaN_{sig}/N_{sig} [%]")
    ents, objs = [], []
    for i, era in enumerate(ERAS):
        e = full[(full.era == era) & (full.ch == "all")].sort_values("mA")
        g = graph(e.mA.to_numpy(float) + (i - 2) * 0.15, 100 * e[f"dy_{k}"], 100 * e[f"dyerr_{k}"], i)
        g.Draw("P same"); ents.append(era); objs.append(g)
    cms_header(c, LUMI_RUN3, "13.6", extra="Preliminary")
    panel_label(c, "Signal, ee+#mu#mu", notes=[lab + " vs nominal"])
    leg = house_legend(c, ents, ncols=2, size=0.036, top=below_frame(0.06 + 0.055 + 0.035, c))
    for o, t in zip(objs, ents):
        leg.AddEntry(o, t, "pe")
    leg.Draw()
    save(c, f"{OUT}/s1_dyield_vs_mA_{k}_per_era")

# 3) shape: sigma_eff relative change and mean shift vs nominal (Run 3)
for q, ytit, scale, rng in (("rseff", "#Delta#sigma_{eff}/#sigma_{eff} [%]", 100, (-8, 17)),
                            ("dmean", "#Delta mean [GeV]", 1, (-0.35, 0.75))):
    c = frame("c_" + q, rng[0], rng[1], ytit)
    ents, objs = [], []
    for i, (k, lab) in enumerate(ALTS):
        g = graph(m + (i - 1.5) * 0.18, scale * a[f"{q}_{k}"], np.zeros(len(m)), i)
        if k in ("zm", "zmnarrow"):
            g.SetMarkerStyle([24, 25][i - 2])
        g.Draw("P same"); ents.append(lab); objs.append(g)
    cms_header(c, LUMI_RUN3, "13.6", extra="Preliminary")
    panel_label(c, "Signal, Run 3, ee+#mu#mu", notes=["m_{ll#gamma#gamma} in fit window, after BDT WP"])
    leg = house_legend(c, ents, size=0.034, top=below_frame(0.06 + 0.055 + 0.035, c))
    for o, t in zip(objs, ents):
        leg.AddEntry(o, t, "pe")
    leg.Draw()
    save(c, f"{OUT}/s1_shape_{q}_vs_mA_run3")
