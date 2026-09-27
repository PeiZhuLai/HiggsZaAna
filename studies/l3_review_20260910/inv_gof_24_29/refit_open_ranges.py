"""In-memory refit of the mA24 envelope members with the turn-on/width/sigma ranges opened
(nothing is written back). Also draws old-vs-new normalized spectra for mA 20, 24, 25, 28."""
import ROOT, math
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/inv_gof_24_29"
N = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3"
O = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/snapshot_preFSRfix_20260924/fit_results_run3"
OPEN = {"turnon": (95., 125.), "step": (95., 125.), "width": (0.1, 50.), "stepWidth": (0.1, 50.), "sigma": (0.05, 20.), "gsigma": (0.05, 20.)}
f = ROOT.TFile.Open(N + "/24/CMS-HGG_mva_13p6TeV_multipdf.root"); w = f.Get("multipdf")
x = w.var("CMS_hza_mass"); d = w.data("roohist_data_mass_cat0")
mp = [p for p in w.allPdfs() if p.InheritsFrom("RooMultiPdf")][0]
for i in range(mp.getNumPdfs()):
    pdf = mp.getPdf(i)
    nll0 = pdf.createNLL(d).getVal()
    for p in pdf.getParameters(ROOT.RooArgSet(x)):
        key = [k for k in OPEN if p.GetName().split("_")[-2 if p.GetName().split("_")[-1].startswith("p") else -1] == k or p.GetName().endswith("_" + k)]
        if key: p.setRange(*OPEN[key[0]])
    r = pdf.fitTo(d, ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Minimizer("Minuit2", "migrad"), ROOT.RooFit.Strategy(1))
    nll1 = pdf.createNLL(d).getVal()
    vals = {p.GetName(): round(p.getVal(), 2) for p in pdf.getParameters(ROOT.RooArgSet(x)) if any(k in p.GetName() for k in ("turnon", "step", "width", "sigma"))}
    print("mA24 %-6s NLL bounded %.1f -> open %.1f (dNLL %.1f) status %d  %s" % (pdf.GetName(), nll0, nll1, nll0 - nll1, r.status(), vals))

c = ROOT.TCanvas("c", "", 1200, 900); c.Divide(2, 2); keep = []
for k, m in enumerate((20, 24, 25, 28)):
    c.cd(k + 1)
    hs = []
    for tag, base, col in (("pre-FSR-fix", O, ROOT.kBlue + 1), ("FSR fix", N, ROOT.kRed + 1)):
        ff = ROOT.TFile.Open("%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (base, m)); ww = ff.Get("multipdf")
        h = ww.data("roohist_data_mass_cat0").createHistogram("h%s%d" % (tag[:3], m), ww.var("CMS_hza_mass"), ROOT.RooFit.Binning(17, 95, 180))
        h.SetDirectory(0); n = h.Integral(); h.Scale(1. / n); h.SetLineColor(col); h.SetMarkerColor(col); h.SetLineWidth(2)
        h.SetTitle(";m_{#ell#ell#gamma#gamma} [GeV];fraction / 5 GeV"); h.GetXaxis().SetTitleSize(0.055); h.GetYaxis().SetTitleSize(0.055)
        h.GetXaxis().SetLabelSize(0.05); h.GetYaxis().SetLabelSize(0.05); hs.append((h, "%s (N=%d)" % (tag, n))); keep += [h, ff]
    hs[0][0].SetMaximum(1.3 * max(h.GetMaximum() for h, _ in hs)); hs[0][0].Draw("E"); hs[1][0].Draw("E SAME")
    leg = ROOT.TLegend(0.45, 0.7, 0.88, 0.88); leg.SetBorderSize(0); leg.SetHeader("m_{a} = %d GeV" % m)
    for h, lab in hs: leg.AddEntry(h, lab, "l")
    leg.Draw(); keep.append(leg)
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.05); lat.DrawLatex(0.12, 0.92, "#bf{CMS} #it{Preliminary}"); keep.append(lat)
c.SaveAs(OUT + "/spectra_old_vs_new_mA20_24_25_28.pdf")
