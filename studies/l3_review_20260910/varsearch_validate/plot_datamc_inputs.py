#!/usr/bin/env python3
"""Data/MC comparison of the input variables of the two BDTs selected by the variable search.

Region "sideband": the preselected signal-region events with m_llgg in 95-115 or 135-180 GeV (the population the BDT
  sees as background). Data vs DY simulation (DYGto2LG + DYJets, 'factor' weight), MC normalized to data. Signal
  (test trees, normalized) overlaid for m_a = 1, 3, 10, 30.
Region "zmmg": probe photons of Z -> mu mu gamma (parquet_tnp_zmmg: tight cut-based ID, 80 < m_mumugamma < 100),
  data vs DY (+ ttbar) simulation, barrel and endcap separately. The simulation is reweighted to the data in
  (pT, |eta|) so the comparison tests the shower shape and isolation modeling, not the kinematics. Photon-level
  variables only (event-level ones are not defined there).
HZa applies no shower-shape correction to the simulation (no flow / quantile correction in HiggsDNA), so the
nominal simulation is the only version.

Per variable: chi2/ndf of data vs normalized MC (data + MC statistical errors), and the largest |data/MC - 1| over
bins with >= 1% of the events. Env: higgs-alp-ana (ROOT 6.24).
Output: <out>/{sideband,zmmg_EB,zmmg_EE}/<var>.{png,pdf}; eval/datamc_summary.txt
"""
import argparse
import ast
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bdt_varsearch"))
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts")
import search_job as sj  # noqa: E402
import ROOT  # noqa: E402
from _cms_style import format_axes  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
# the authoritative lumi table; only the LUMI literal is taken (the module itself does not import: ERA uses sets as keys)
_CONST = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
_LUMI = next(ast.literal_eval(n.value) for n in ast.parse(open(_CONST).read()).body
             if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "LUMI")
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
LUMI = sum(_LUMI[e] for e in ERAS)
LUMI_LABEL = "%.1f fb^{-1} (13.6 TeV)" % LUMI

SEL = json.load(open(os.path.join(HERE, "selected.json")))
TNP = "/eos/home-p/pelai/HZa/parquet_tnp_zmmg"
TNP_MC = ["DYJetsToLL_2022preEE", "DYJetsToLL_2022postEE", "DYJetsToLL_2023preBPix", "DYJetsToLL_2023postBPix",
          "DYJetsTo2Mu_2024", "TT_2022preEE", "TT_2022postEE", "TT_2023preBPix", "TT_2023postBPix", "TT_2024"]
# analysis input -> Z->mumugamma probe column (lead/sublead both map to the probe photon)
ZMMG = {"sieip": "probe_photon_sieip", "etaWidth": "probe_photon_etaWidth", "phiWidth": "probe_photon_phiWidth",
        "sieie": "probe_photon_sieie", "trkSumPtHollowConeDR03": "probe_photon_trkSumPtHollowConeDR03",
        "pfChargedIsoWorstVtx": "probe_photon_pfChargedIsoWorstVtx", "chiso": "probe_photon_chiso",
        "ecalPFClusterIso": "probe_photon_ecalPFClusterIso", "relEerr": None, "pt": "probe_photon_pt",
        "abseta": None}
TITLES = {"min_lg_deltaR": "min #DeltaR(#gamma,l)", "pho2_relEerr": "#sigma_{E}/E (sublead #gamma)",
          "pho1Pt": "p_{T}^{#gamma1} [GeV]", "H_pt": "p_{T}^{ll#gamma#gamma} [GeV]", "var_PtaOverMa": "p_{T}^{a}/m_{a}",
          "var_dEta_g1Z": "#Delta#eta(#gamma1,Z)", "pho1IetaIeta55": "#sigma_{i#etai#eta}^{5x5} (lead #gamma)",
          "pho1PIso_noCorr": "ECAL PF cluster iso (lead #gamma) [GeV]", "ALP_calculatedPhotonIso": "a-candidate photon iso",
          "relEerr": "#sigma_{E}/E", "pt": "p_{T}^{#gamma} [GeV]", "abseta": "|#eta^{#gamma}|"}


def title(v):
    if v in TITLES:
        return TITLES[v]
    for w in ("lead", "sublead"):
        p = "ALP_%s_photon_" % w
        if v.startswith(p):
            return "%s (%s #gamma)" % (v[len(p):], w)
    return v


def th1(name, vals, w, edges):
    h = ROOT.TH1D(name, "", len(edges) - 1, np.asarray(edges, float))
    h.Sumw2()
    c, _ = np.histogram(np.clip(vals, edges[0], edges[-1] - 1e-9 * (edges[-1] - edges[0])), edges, weights=w)
    e, _ = np.histogram(np.clip(vals, edges[0], edges[-1] - 1e-9 * (edges[-1] - edges[0])), edges, weights=w * w)
    for i in range(len(c)):
        h.SetBinContent(i + 1, c[i]); h.SetBinError(i + 1, np.sqrt(e[i]))
    return h


def edges_for(x):
    x = x[np.isfinite(x)]
    lo, hi = np.percentile(x, [0.5, 99.5])
    if np.isclose(lo, hi):
        hi = lo + 1
    return np.linspace(lo, hi, 41)


def draw(out, var, data, wdata, mc, wmc, sigs, label, extra=""):
    edges = edges_for(data)
    hd = th1("d", data, wdata, edges); hm = th1("m", mc, wmc, edges)
    hm.Scale(hd.Integral() / hm.Integral())
    dcont = np.array([hd.GetBinContent(i + 1) for i in range(len(edges) - 1)])
    derr = np.array([hd.GetBinError(i + 1) for i in range(len(edges) - 1)])
    mcont = np.array([hm.GetBinContent(i + 1) for i in range(len(edges) - 1)])
    merr = np.array([hm.GetBinError(i + 1) for i in range(len(edges) - 1)])
    ok = (derr ** 2 + merr ** 2) > 0
    chi2 = float(np.sum((dcont[ok] - mcont[ok]) ** 2 / (derr[ok] ** 2 + merr[ok] ** 2)))
    ndf = int(ok.sum() - 1)
    big = (dcont >= 0.01 * dcont.sum()) & (mcont > 0)
    maxdev = float(np.max(np.abs(dcont[big] / mcont[big] - 1))) if big.any() else float("nan")
    meanshift = (np.average(data, weights=wdata) - np.average(mc, weights=wmc)) / max(np.sqrt(np.cov(data, aweights=wdata)), 1e-12)

    c = ROOT.TCanvas("c", "", 800, 800)
    p1 = ROOT.TPad("p1", "", 0, 0.3, 1, 1); p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.15); p1.SetRightMargin(0.05); p1.SetTopMargin(0.08)
    p2 = ROOT.TPad("p2", "", 0, 0, 1, 0.3); p2.SetTopMargin(0.03); p2.SetBottomMargin(0.38); p2.SetLeftMargin(0.15); p2.SetRightMargin(0.05)
    p1.Draw(); p2.Draw(); p1.cd()
    hm.SetFillColor(ROOT.kOrange - 3); hm.SetLineColor(ROOT.kOrange + 2)
    hd.SetMarkerStyle(20); hd.SetMarkerSize(0.9); hd.SetLineColor(ROOT.kBlack)
    hs = []
    cols = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]
    for k, (lab, x, w) in enumerate(sigs):
        h = th1("s%d" % k, x, w, edges); h.Scale(hd.Integral() / h.Integral() if h.Integral() > 0 else 1)
        h.SetLineColor(cols[k % 4]); h.SetLineWidth(3); h.SetLineStyle(2); hs.append((lab, h))
    ymax = max([hd.GetMaximum(), hm.GetMaximum()] + [h.GetMaximum() for _, h in hs]) * 1.45
    hm.SetMaximum(ymax); hm.SetMinimum(0)
    format_axes(hm); hm.GetXaxis().SetLabelSize(0); hm.GetYaxis().SetTitle("Events / bin (MC norm. to data)")
    hm.GetYaxis().SetTitleSize(0.055); hm.GetYaxis().SetLabelSize(0.050)
    hm.Draw("hist"); hme = hm.Clone("hme"); hme.SetFillColor(ROOT.kGray + 2); hme.SetFillStyle(3354); hme.SetMarkerSize(0); hme.Draw("e2 same")
    for _, h in hs:
        h.Draw("hist same")
    hd.Draw("e1 same")
    leg = ROOT.TLegend(0.55, 0.60, 0.93, 0.89); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.040)
    leg.AddEntry(hd, "Data", "pe"); leg.AddEntry(hm, label, "f")
    for lab, h in hs:
        leg.AddEntry(h, lab, "l")
    leg.Draw()
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.050); lat.SetTextFont(42)
    lat.DrawLatex(0.15, 0.935, "#bf{CMS} #it{Preliminary}")
    lat.SetTextAlign(31); lat.DrawLatex(0.95, 0.935, LUMI_LABEL)
    lat.SetTextAlign(11); lat.SetTextSize(0.038)
    lat.DrawLatex(0.19, 0.84, extra)
    lat.DrawLatex(0.19, 0.78, "#chi^{2}/ndf = %.1f/%d" % (chi2, ndf))
    p2.cd()
    r = hd.Clone("r"); r.Divide(hm)
    for i in range(1, r.GetNbinsX() + 1):
        r.SetBinError(i, hd.GetBinError(i) / hm.GetBinContent(i) if hm.GetBinContent(i) > 0 else 0)
    band = hm.Clone("band")
    for i in range(1, band.GetNbinsX() + 1):
        b = band.GetBinContent(i); band.SetBinError(i, band.GetBinError(i) / b if b > 0 else 0); band.SetBinContent(i, 1)
    band.SetFillColor(ROOT.kGray + 1); band.SetFillStyle(3354); band.SetMarkerSize(0); band.SetLineColor(ROOT.kGray + 1)
    format_axes(band, x_title_offset=1.0, y_title_offset=0.5)
    for ax in (band.GetXaxis(), band.GetYaxis()):
        ax.SetTitleSize(0.055 * 0.7 / 0.3); ax.SetLabelSize(0.050 * 0.7 / 0.3)
    band.GetYaxis().SetNdivisions(505); band.GetYaxis().SetTitle("Data/MC"); band.GetYaxis().SetTitleOffset(0.5)
    band.GetXaxis().SetTitle(title(var)); band.SetMinimum(0.5); band.SetMaximum(1.5)
    band.Draw("e2"); r.Draw("e1 same")
    one = ROOT.TLine(edges[0], 1, edges[-1], 1); one.SetLineStyle(2); one.Draw()
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "pdf"):
        c.SaveAs(os.path.join(out, "%s.%s" % (var, ext)))
    c.Close()
    return chi2, ndf, maxdev, float(meanshift)


def sideband(out, cache, rows):
    S = {n: dict(np.load(os.path.join(cache, n + ".npz"))) for n in ("data_sb", "dy", "test_sig_m1", "test_sig_m3", "test_sig_m10", "test_sig_m30")}
    feats = []
    for model in ("low", "high"):
        for f in SEL[model]["feats"]:
            if f != "param_rel" and f not in feats:
                feats.append(f)
    d, m = S["data_sb"], S["dy"]
    for f in feats:
        used = "/".join(k for k in ("low", "high") if f in SEL[k]["feats"])
        x_d = sj.column(d, f, 1.0); x_m = sj.column(m, f, 1.0)
        okm = np.isfinite(x_m) & np.isfinite(m["factor"])
        sigs = [("m_{a} = %d GeV (shape)" % k, sj.column(S["test_sig_m%d" % k], f, float(k)), S["test_sig_m%d" % k]["factor"])
                for k in (1, 3, 10, 30)]
        sigs = [(l, x[np.isfinite(x)], w[np.isfinite(x)]) for l, x, w in sigs]
        r = draw(out, f, x_d[np.isfinite(x_d)], np.ones(np.isfinite(x_d).sum()), x_m[okm], m["factor"][okm], sigs,
                 "DY simulation", "m_{ll#gamma#gamma} sidebands, model: %s" % used)
        rows.append(("sideband", f, used) + r)
        print("  sideband %-45s chi2/ndf %7.1f/%d  max|r-1| %.2f  mean shift %+.2f sd" % ((f,) + r), flush=True)


def load_tnp(names, kind):
    cols = sorted({c for c in ZMMG.values() if c} | {"probe_photon_eta", "probe_photon_energyErr", "probe_photon_isScEtaEB",
                                                    "probe_photon_isScEtaEE", "weight_central"})
    acc = {}
    for n in names:
        f = "%s/%s/%s/merged_nominal.parquet" % (TNP, kind, n)
        t = pq.read_table(f, columns=cols).to_pandas()
        for c in cols:
            acc.setdefault(c, []).append(t[c].to_numpy().astype(np.float64))
    a = {c: np.concatenate(v) for c, v in acc.items()}
    a["relEerr"] = a["probe_photon_energyErr"] / (a["probe_photon_pt"] * np.cosh(a["probe_photon_eta"]))
    a["abseta"] = np.abs(a["probe_photon_eta"])
    return a


def zmmg(out, rows):
    d = load_tnp(["Data_tnp_zmmg_%s" % e for e in ERAS], "data")
    m = load_tnp(TNP_MC, "mc")
    for reg, key in (("EB", "probe_photon_isScEtaEB"), ("EE", "probe_photon_isScEtaEE")):
        sd, sm = d[key] > 0.5, m[key] > 0.5
        wd = np.ones(sd.sum()); wm = m["weight_central"][sm]
        # (pT, |eta|) reweighting of the simulation to the data
        pe = np.array([10, 12, 15, 20, 25, 30, 40, 60, 1e9]); ee = np.quantile(d["abseta"][sd], np.linspace(0, 1, 6))
        ee[0], ee[-1] = 0, 3
        hd, _, _ = np.histogram2d(d["probe_photon_pt"][sd], d["abseta"][sd], [pe, ee])
        hm, _, _ = np.histogram2d(m["probe_photon_pt"][sm], m["abseta"][sm], [pe, ee], weights=wm)
        ratio = np.where(hm > 0, (hd / hd.sum()) / np.where(hm > 0, hm / hm.sum(), 1), 1.0)
        ip = np.clip(np.digitize(m["probe_photon_pt"][sm], pe) - 1, 0, len(pe) - 2)
        ie = np.clip(np.digitize(m["abseta"][sm], ee) - 1, 0, len(ee) - 2)
        wrw = wm * ratio[ip, ie]
        for v, col in ZMMG.items():
            xd = d[col][sd] if col else d[v][sd]
            xm = m[col][sm] if col else m[v][sm]
            kin = v in ("pt", "abseta")
            w_use = wm if kin else wrw
            ok = np.isfinite(xm) & np.isfinite(w_use)
            r = draw(os.path.join(out, "zmmg_%s" % reg), v, xd[np.isfinite(xd)], wd[np.isfinite(xd)], xm[ok], w_use[ok], [],
                     "Z#rightarrow#mu#mu#gamma sim." + ("" if kin else " (p_{T},|#eta| rw)"),
                     "Z#rightarrow#mu#mu#gamma probe #gamma, %s" % reg)
            rows.append(("zmmg_%s" % reg, v, "probe") + r)
            print("  zmmg %s %-25s chi2/ndf %7.1f/%d  max|r-1| %.2f  mean shift %+.2f sd" % ((reg, v) + r), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="/eos/home-p/pelai/HZa/sandbox_varsearch_v1/eval_samples64")
    ap.add_argument("--out", default="/eos/home-p/pelai/HZa/sandbox_varsearch_v1/datamc")
    ap.add_argument("--regions", default="sideband,zmmg")
    a = ap.parse_args()
    rows = []
    if "sideband" in a.regions:
        sideband(os.path.join(a.out, "sideband"), a.cache, rows)
    if "zmmg" in a.regions:
        zmmg(a.out, rows)
    with open(os.path.join(HERE, "eval", "datamc_summary_%s.txt" % a.regions.replace(",", "_")), "w") as f:
        f.write("# region variable model chi2 ndf max|data/MC-1| (bins>=1%%) mean-shift[sd]; plots in %s\n" % a.out)
        for r in rows:
            f.write("%-10s %-45s %-9s %9.1f %3d %6.3f %+6.3f\n" % r)


if __name__ == "__main__":
    main()
