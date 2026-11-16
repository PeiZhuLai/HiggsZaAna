#!/usr/bin/env python3
"""Post-fit distributions on the background-only Asimov dataset (Combine review: pre/post-fit with bands).

Input : Combine/root_t2w_xsec/<m>_Datacard_leptons.root, the b-only Asimov dataset
        (Checks/review_round2/asimov/higgsCombine_asimov_mA<m>...root, toys/toy_asimov) and the s+b fit of it
        (Checks/review_round2/fitdiag/xsec_mA<m>/multidimfit_sb.root, RooFitResult fit_mdf, Hesse covariance).
Output: /eos/home-p/pelai/HZa/output_plots/combine_review_round2/postfit_asimov_mA<m>.{pdf,png}
Top: Asimov dataset in 1 GeV bins, post-fit background model (envelope member chosen by the fit) with a band
from 500 parameter sets drawn from the fit covariance, and the signal model at the median expected limit (r = q50).
The band is the RMS of the sampled curves around the post-fit curve in each bin, the convention of the Combine
FitDiagnostics shapes. The 16-84% quantiles of the sampled curves (used before) are not centered on the post-fit curve
when a parameter sits at its range limit with a large Hesse error (mA20: width 0.13 +- 41, range 0.1-50, covariance
forced positive definite): the draws are clipped at the limit and the quantile band moves away from the curve. Bottom: Asimov minus post-fit background, with the same band and the signal.
Nothing here uses the observed data.
"""
import ast, math, os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)  # NaN warnings of the sampled sets fill GBs of log
ROOT.RooAbsReal.setEvalErrorLoggingMode(ROOT.RooAbsReal.Ignore)
FF = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
CHK = FF + "/Combine/Checks/review_round2"
OUT = "/eos/home-p/pelai/HZa/output_plots/combine_review_round2"
W = os.path.dirname(os.path.abspath(__file__))
NSAMP = 500


def read_lumi():
    src = open("/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py").read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "LUMI":
            return ast.literal_eval(node.value)


LUMI = read_lumi()
LUMI_LABEL = "%.0f fb^{-1} (13.6 TeV)" % (LUMI["2022"] + LUMI["2023"] + LUMI["2024"])
Q50 = {int(l.split()[0]): float(l.split()[1]) for l in open(W + "/q50.txt")}


def main():
    os.makedirs(OUT, exist_ok=True)
    for m in (1, 10, 20, 30):
        fw = ROOT.TFile.Open("%s/Combine/root_t2w_xsec/%d_Datacard_leptons.root" % (FF, m)); w = fw.Get("w")
        fr = ROOT.TFile.Open("%s/fitdiag/xsec_mA%d/multidimfit_sb.root" % (CHK, m)); res = fr.Get("fit_mdf")
        fa = ROOT.TFile.Open("%s/asimov/higgsCombine_asimov_mA%d.GenerateOnly.mH125.38.123456.root" % (CHK, m))
        asimov = fa.Get("toys/toy_asimov")
        x = w.var("CMS_hza_mass"); lo, hi = x.getMin(), x.getMax(); nb = int(round(hi - lo))
        # post-fit values, including the envelope index chosen by the fit
        w.allVars().assign(res.floatParsFinal())
        idx = res.constPars().find("pdfindex_cat0_13p6TeV")
        w.cat("pdfindex_cat0_13p6TeV").setIndex(idx.getIndex() if idx else 0)
        bpdf = w.pdf("shapeBkg_bkg_mass_cat0")
        bnorm = w.function("n_exp_final_bincat0_proc_bkg_mass")
        procs = sorted(p.GetName()[len("shapeSig_"):-len("_cat0")] for p in w.allPdfs()
                       if p.GetName().startswith("shapeSig_") and p.GetName().endswith("_cat0"))
        obs = ROOT.RooArgSet(x)
        centers = [lo + i + 0.5 for i in range(nb)]

        def bkg_curve():
            out = []
            n = bnorm.getVal()
            for c in centers:  # 1 GeV bins: density x bin width x yield
                x.setVal(c); out.append(n * bpdf.getVal(obs))
            return out

        nom = bkg_curve()
        # Sample only the parameters of the selected background function and its normalization, with the
        # covariance of the fit reduced to them (the full matrix contains zero-width parameters of the
        # unselected envelope members, for which RooFitResult::randomizePars returns NaN).
        import numpy as np
        sel = bpdf.getCurrentPdf().getParameters(obs)
        names = [q.GetName() for q in res.floatParsFinal() if q.getError() > 0 and
                 (sel.find(q.GetName()) or q.GetName() == "shapeBkg_bkg_mass_cat0__norm")]
        lst = ROOT.RooArgList()
        for n in names:
            lst.add(res.floatParsFinal().find(n))
        cov = res.reducedCovarianceMatrix(lst)
        C = np.array([[cov[i][j] for j in range(len(names))] for i in range(len(names))])
        mu = np.array([res.floatParsFinal().find(n).getVal() for n in names])
        rng = np.random.default_rng(12345)
        samples = []
        tries = 0
        while len(samples) < NSAMP and tries < 5 * NSAMP:
            tries += 1
            vals = rng.multivariate_normal(mu, C)
            for n, v in zip(names, vals):
                w.var(n).setVal(float(v))  # setVal clips to the parameter range
            cur = bkg_curve()
            if all(math.isfinite(v) for v in cur):
                samples.append(cur)
        print("mA%d: sampled %s" % (m, ",".join(names)))
        w.allVars().assign(res.floatParsFinal())
        print("mA%d: %d finite parameter sets out of %d drawn" % (m, len(samples), tries))
        band = []
        for i in range(nb):
            rms = math.sqrt(sum((s[i] - nom[i]) ** 2 for s in samples) / len(samples))
            band.append((nom[i] - rms, nom[i] + rms))
        # signal at r = q50, nuisances at the post-fit values (all zero on this Asimov)
        w.var("r").setVal(Q50[m])
        sig = []
        norms = [w.function("n_exp_final_bincat0_proc_" + p) for p in procs]
        pdfs = [w.pdf("shapeSig_%s_cat0" % p) for p in procs]
        for c in centers:
            x.setVal(c); sig.append(sum(nn.getVal() * pp.getVal(obs) for nn, pp in zip(norms, pdfs)))
        # Asimov dataset (weighted) in 1 GeV bins
        hA = ROOT.TH1D("hA%d" % m, "", nb, lo, hi); hA.Sumw2()
        for i in range(asimov.numEntries()):
            row = asimov.get(i); hA.Fill(row.getRealValue("CMS_hza_mass"), asimov.weight())
        for i in range(1, nb + 1):
            hA.SetBinError(i, hA.GetBinContent(i) ** 0.5)
        gB = ROOT.TGraph(nb); gBand = ROOT.TGraphAsymmErrors(nb); gS = ROOT.TGraph(nb)
        gR = ROOT.TGraphAsymmErrors(nb); gRB = ROOT.TGraphAsymmErrors(nb); gRS = ROOT.TGraph(nb)
        for i, c in enumerate(centers):
            gB.SetPoint(i, c, nom[i]); gS.SetPoint(i, c, nom[i] + sig[i])
            gBand.SetPoint(i, c, nom[i]); gBand.SetPointError(i, 0.5, 0.5, nom[i] - band[i][0], band[i][1] - nom[i])
            d = hA.GetBinContent(i + 1)
            gR.SetPoint(i, c, d - nom[i]); gR.SetPointError(i, 0, 0, d ** 0.5, d ** 0.5)
            gRB.SetPoint(i, c, 0); gRB.SetPointError(i, 0.5, 0.5, nom[i] - band[i][0], band[i][1] - nom[i])
            gRS.SetPoint(i, c, sig[i])
        c1 = ROOT.TCanvas("c%d" % m, "", 900, 950)
        top = ROOT.TPad("t", "", 0, 0.32, 1, 1); bot = ROOT.TPad("b", "", 0, 0, 1, 0.32)
        top.SetLeftMargin(0.17); top.SetRightMargin(0.04); top.SetTopMargin(0.09 / 0.68); top.SetBottomMargin(0.02)
        bot.SetLeftMargin(0.17); bot.SetRightMargin(0.04); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.42)
        top.Draw(); bot.Draw(); top.cd()
        ymax = 2.0 * max(max(hA.GetBinContent(i) for i in range(1, nb + 1)), max(b[1] for b in band), max(nom[i] + sig[i] for i in range(nb)))
        f1 = top.DrawFrame(lo, -0.02 * ymax, hi, ymax)
        f1.GetYaxis().SetTitle("Events / GeV"); f1.GetXaxis().SetLabelSize(0)
        f1.GetYaxis().SetTitleSize(0.055 / 0.68); f1.GetYaxis().SetLabelSize(0.050 / 0.68); f1.GetYaxis().SetTitleOffset(1.05)
        gBand.SetFillColor(ROOT.kGreen - 8); gBand.SetLineColor(0); gBand.Draw("2")
        gB.SetLineColor(ROOT.kBlue + 1); gB.SetLineWidth(3); gB.Draw("L")
        gS.SetLineColor(ROOT.kRed + 1); gS.SetLineWidth(3); gS.SetLineStyle(2); gS.Draw("L")
        hA.SetMarkerStyle(20); hA.SetMarkerSize(0.9); hA.SetLineColor(ROOT.kBlack); hA.Draw("E0 P same")
        ytop = 1 - top.GetTopMargin()
        lab = ROOT.TLatex(); lab.SetNDC(); lab.SetTextFont(42); lab.SetTextSize(0.040 / 0.68)
        lab.DrawLatex(0.21, ytop - 0.09 / 0.68, "m_{a} = %d GeV" % m)
        lab.DrawLatex(0.21, ytop - 0.15 / 0.68, "B-only Asimov")
        leg = ROOT.TLegend(0.48, ytop - 0.02 / 0.68 - 4 * 0.065 / 0.68, 0.93, ytop - 0.02 / 0.68)
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.040 / 0.68); leg.SetMargin(0.22)
        leg.AddEntry(hA, "Asimov (A)", "pe")
        leg.AddEntry(gB, "Post-fit bkg. (B)", "l")
        leg.AddEntry(gBand, "#pm1#sigma bkg.", "f")
        leg.AddEntry(gS, "B + signal (r = %.3g)" % Q50[m], "l")
        leg.Draw()
        t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.050 / 0.68)
        t.DrawLatex(0.17, ytop + 0.02 / 0.68, "#bf{CMS} #it{Preliminary}")
        t.SetTextAlign(31); t.DrawLatex(0.96, ytop + 0.02 / 0.68, LUMI_LABEL)
        bot.cd()
        rr = 1.3 * max(max(hA.GetBinContent(i + 1) ** 0.5 for i in range(nb)), max(sig))
        f2 = bot.DrawFrame(lo, -rr, hi, rr)
        f2.GetXaxis().SetTitle("m_{ll#gamma#gamma} [GeV]"); f2.GetYaxis().SetTitle("A #minus B")
        for ax in (f2.GetXaxis(), f2.GetYaxis()):
            ax.SetTitleSize(0.055 / 0.32); ax.SetLabelSize(0.050 / 0.32)
        f2.GetYaxis().SetTitleOffset(0.46); f2.GetYaxis().CenterTitle(); f2.GetYaxis().SetNdivisions(505); f2.GetXaxis().SetTitleOffset(1.05)
        f2.GetXaxis().SetLabelOffset(0.015)
        gRB.SetFillColor(ROOT.kGreen - 8); gRB.SetLineColor(0); gRB.Draw("2")
        gRS.SetLineColor(ROOT.kRed + 1); gRS.SetLineWidth(3); gRS.SetLineStyle(2); gRS.Draw("L")
        gR.SetMarkerStyle(20); gR.SetMarkerSize(0.9); gR.Draw("P")
        ln = ROOT.TLine(lo, 0, hi, 0); ln.SetLineStyle(3); ln.Draw()
        for ext in ("pdf", "png"):
            c1.SaveAs("%s/postfit_asimov_mA%d.%s" % (OUT, m, ext))
        for f in (fw, fr, fa):
            f.Close()


if __name__ == "__main__":
    main()
