#!/usr/bin/env python3
"""Summary of the Asimov (r = 0) fits of the Combine-review checks.

Input : <checks>/fitdiag/<card>_mA<m>/multidimfit_{sb,b}.root   (RooFitResult fit_mdf)
        the packaged datacards, to know which parameters carry a constraint (lnN / param)
Output: <W>/fit_summary.txt
        <plots>/corr_<card>_mA<m>.pdf  correlation matrix of the floating parameters (s+b fit)
Per card and mass point: fit status and covariance quality of both fits, r_hat +- sigma_r, the
largest |pull| and the smallest post-fit/pre-fit width among the constrained nuisances (both are
in units of the pre-fit width, which is 1 for lnN and param nuisances in Combine), and the
parameter most correlated with r.
"""
import os, sys
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDF output")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
W = os.path.dirname(os.path.abspath(__file__))
CHK = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/review_round2"
DC = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons"
PLOTS = "/eos/home-p/pelai/HZa/output_plots/combine_review_round2"
CORR_MASSES = (1, 10, 20, 30)


def constrained(card_path):
    names = set()
    for line in open(card_path):
        t = line.split()
        if len(t) > 1 and t[1] in ("lnN", "param"):
            names.add(t[0])
    return names


def read_lumi_label():
    import ast
    src = open("/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py").read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "LUMI":
            lumi = ast.literal_eval(node.value)
            return "%.0f fb^{-1} (13.6 TeV)" % (lumi["2022"] + lumi["2023"] + lumi["2024"])
    raise RuntimeError("LUMI not found in constants.py")


# Axis labels: the parameter names of the workspace, except the background normalization, whose name is too long
# for the axis (the AN caption gives the mapping).
LABEL = {"shapeBkg_bkg_mass_cat0__norm": "N_{bkg}"}


def corr_plot(res, out):
    allp = [p.GetName() for p in res.floatParsFinal() if p.getError() > 0]
    # keep only parameters correlated with at least one other one; the rest is the identity
    pars = [a for a in allp if a == "r" or any(abs(res.correlation(a, b)) > 0.05 for b in allp if b != a)]
    n = len(pars)
    h = ROOT.TH2D("h", "", n, 0, n, n, 0, n)
    for i, a in enumerate(pars):
        h.GetXaxis().SetBinLabel(i + 1, LABEL.get(a, a)); h.GetYaxis().SetBinLabel(i + 1, LABEL.get(a, a))
        for j, b in enumerate(pars):
            h.SetBinContent(i + 1, j + 1, res.correlation(a, b))
    c = ROOT.TCanvas("c", "", 1000, 900)
    c.SetLeftMargin(0.25); c.SetBottomMargin(0.25); c.SetRightMargin(0.20); c.SetTopMargin(0.08)
    h.SetMinimum(-1); h.SetMaximum(1)
    h.GetXaxis().SetLabelSize(0.050); h.GetYaxis().SetLabelSize(0.050); h.GetXaxis().LabelsOption("v")
    h.GetZaxis().SetLabelSize(0.045); h.GetZaxis().SetTitle("Correlation"); h.GetZaxis().SetTitleSize(0.050)
    h.GetZaxis().SetTitleOffset(1.30)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetNumberContours(100)
    h.Draw("colz")
    # cell values drawn by hand: white on the dark (negative) end of kBird, black elsewhere
    t = ROOT.TLatex(); t.SetTextFont(42); t.SetTextAlign(22); t.SetTextSize(min(0.040, 0.22 / n))
    for i in range(n):
        for j in range(n):
            v = h.GetBinContent(i + 1, j + 1)
            t.SetTextColor(ROOT.kWhite if v < -0.35 else ROOT.kBlack)
            t.DrawLatex(i + 0.5, j + 0.5, "%.2f" % v)
    hdr = ROOT.TLatex(); hdr.SetNDC(); hdr.SetTextFont(42); hdr.SetTextSize(0.040)
    hdr.DrawLatex(c.GetLeftMargin(), 1 - c.GetTopMargin() + 0.02, "#bf{CMS} #it{Preliminary}")
    hdr.SetTextAlign(31); hdr.DrawLatex(1 - c.GetRightMargin(), 1 - c.GetTopMargin() + 0.02, read_lumi_label())
    c.SaveAs(out)


def main():
    os.makedirs(PLOTS, exist_ok=True)
    L = ["card   mA  | s+b: status covQual  r_hat +- sigma_r | b-only: status covQual | "
         "max|pull| (param) | min constraint (param) | max|corr(r,x)| (param)"]
    for card in ("xsec", "wilson"):
        cons = None
        for m in CORR_MASSES:
            d = "%s/fitdiag/%s_mA%d" % (CHK, card, m)
            if cons is None:
                cons = constrained("%s/%d_pruned_datacard_leptons%s.txt" % (DC, m, "_xsec" if card == "xsec" else ""))
            try:
                fsb = ROOT.TFile.Open(d + "/multidimfit_sb.root"); rsb = fsb.Get("fit_mdf")
                fb = ROOT.TFile.Open(d + "/multidimfit_b.root"); rb = fb.Get("fit_mdf")
            except Exception:
                rsb = rb = None
            if not rsb or not rb:
                L.append("%-6s %2d  | missing fit result" % (card, m)); continue
            fl = {p.GetName(): p for p in rsb.floatParsFinal()}
            r = fl["r"]
            cp = [(abs(p.getVal()), n) for n, p in fl.items() if n in cons]
            cw = [(p.getError(), n) for n, p in fl.items() if n in cons and p.getError() > 0]
            cr = [(abs(rsb.correlation("r", n)), n) for n, p in fl.items() if n != "r" and p.getError() > 0]
            mp, mc, mr = max(cp), min(cw), max(cr)
            L.append("%-6s %2d  | %d %d  %+.5f +- %.5f | %d %d | %.3f (%s) | %.3f (%s) | %.3f (%s)" % (
                card, m, rsb.status(), rsb.covQual(), r.getVal(), r.getError(), rb.status(), rb.covQual(),
                mp[0], mp[1], mc[0], mc[1], mr[0], mr[1]))
            if True:
                corr_plot(rsb, "%s/corr_%s_mA%d.pdf" % (PLOTS, card, m)); corr_plot(rsb, "%s/corr_%s_mA%d.png" % (PLOTS, card, m))
    open(W + "/fit_summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
