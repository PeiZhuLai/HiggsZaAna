#!/usr/bin/env python3
"""Signal shape variations of the datacard model (Combine review: "plots of the shape variations").

Input : Combine/root_t2w_xsec/<m>_Datacard_leptons.root (production workspace)
Output: /eos/home-p/pelai/HZa/output_plots/combine_review_round2/sigshape_syst_mA<m>.{pdf,png}
For each mass point the total signal model (sum over the ten year x lepton-flavor processes, each weighted
by its expected yield) is drawn at the nominal values of the nuisances and with all energy-scale nuisances
(electron, muon, photon, FNUF, material) shifted together by +1 and -1 sigma (left), and with all resolution
nuisances shifted together by +1 and -1 sigma (right). A coherent shift of all nuisances of a group is the
largest variation the model allows for that group; the individual nuisances are smaller.
"""
import os, sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def read_lumi():
    """LUMI table of the authoritative constants.py (the module itself does not import in the CMSSW python)."""
    import ast
    src = open("/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py").read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "LUMI":
            return ast.literal_eval(node.value)
    raise RuntimeError("LUMI not found in constants.py")


LUMI = read_lumi()

FF = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
OUT = "/eos/home-p/pelai/HZa/output_plots/combine_review_round2"
LUMI_LABEL = "%.0f fb^{-1} (13.6 TeV)" % (LUMI["2022"] + LUMI["2023"] + LUMI["2024"])
XMIN, XMAX, NPT = 115.0, 135.0, 400


def curve(w, procs, x, shifts):
    nus = w.genobj("ModelConfig").GetNuisanceParameters()
    for v in nus:
        v.setVal(0.0)
    for n, val in shifts.items():
        w.var(n).setVal(val)
    g = ROOT.TGraph(NPT)
    norms = [w.function("n_exp_final_bincat0_proc_" + p).getVal() for p in procs]
    pdfs = [w.pdf("shapeSig_%s_cat0" % p) for p in procs]
    obs = ROOT.RooArgSet(x)
    for i in range(NPT):
        xv = XMIN + (XMAX - XMIN) * i / (NPT - 1)
        x.setVal(xv)
        g.SetPoint(i, xv, sum(nn * pp.getVal(obs) for nn, pp in zip(norms, pdfs)))
    return g


def main():
    os.makedirs(OUT, exist_ok=True)
    for m in (1, 10, 20, 30):
        f = ROOT.TFile.Open("%s/Combine/root_t2w_xsec/%d_Datacard_leptons.root" % (FF, m))
        w = f.Get("w")
        w.var("MH").setVal(125.38)
        w.var("r").setVal(1.0)
        x = w.var("CMS_hza_mass")
        procs = sorted(p.GetName()[len("shapeSig_"):-len("_cat0")] for p in w.allPdfs()
                       if p.GetName().startswith("shapeSig_") and p.GetName().endswith("_cat0"))
        nus = [v.GetName() for v in w.genobj("ModelConfig").GetNuisanceParameters()]
        groups = {"Energy scale": [n for n in nus if "_scale_" in n],
                  "Resolution": [n for n in nus if "_res_" in n]}
        c = ROOT.TCanvas("c%d" % m, "", 1600, 860)
        c.Divide(2, 1)
        keep = []
        for ip, (label, names) in enumerate(groups.items()):
            c.cd(ip + 1)
            top = ROOT.TPad("t%d%d" % (m, ip), "", 0, 0.32, 1, 1); bot = ROOT.TPad("b%d%d" % (m, ip), "", 0, 0, 1, 0.32)
            top.SetLeftMargin(0.17); top.SetRightMargin(0.04); top.SetTopMargin(0.09); top.SetBottomMargin(0.02)
            bot.SetLeftMargin(0.17); bot.SetRightMargin(0.04); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.45)
            top.Draw(); bot.Draw()
            nom = curve(w, procs, x, {})
            up = curve(w, procs, x, {n: 1.0 for n in names})
            dn = curve(w, procs, x, {n: -1.0 for n in names})
            ymax = max(nom.GetY()[i] for i in range(NPT))
            ratios = []
            for g in (up, dn):
                rg = ROOT.TGraph()  # only where the nominal model is above 2% of its peak
                for i in range(NPT):
                    if nom.GetY()[i] > 2e-2 * ymax:
                        rg.SetPoint(rg.GetN(), g.GetX()[i], g.GetY()[i] / nom.GetY()[i])
                ratios.append(rg)
            for g in (nom, up, dn):  # unit peak of the nominal curve
                for i in range(NPT):
                    g.SetPoint(i, g.GetX()[i], g.GetY()[i] / ymax)
            for g, col, sty in ((nom, ROOT.kBlack, 1), (up, ROOT.kRed + 1, 2), (dn, ROOT.kAzure + 2, 2),
                                (ratios[0], ROOT.kRed + 1, 1), (ratios[1], ROOT.kAzure + 2, 1)):
                g.SetLineColor(col); g.SetLineWidth(3); g.SetLineStyle(sty)
            top.cd()
            fr = top.DrawFrame(XMIN, 0, XMAX, 1.4)
            fr.GetYaxis().SetTitle("A.U. (nominal peak = 1)")
            fr.GetXaxis().SetLabelSize(0)
            fr.GetYaxis().SetTitleSize(0.055 / 0.68); fr.GetYaxis().SetLabelSize(0.050 / 0.68); fr.GetYaxis().SetTitleOffset(1.05)
            dn.Draw("L"); up.Draw("L"); nom.Draw("L")
            leg = ROOT.TLegend(0.20, 0.55, 0.62, 0.87)
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.055)
            leg.SetHeader("m_{a} = %d GeV, %s (%d nuisances)" % (m, label.lower(), len(names)))
            leg.AddEntry(nom, "Nominal", "l")
            leg.AddEntry(up, "All +1#sigma", "l")
            leg.AddEntry(dn, "All #minus1#sigma", "l")
            leg.Draw()
            t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.065)
            t.DrawLatex(0.17, 0.93, "#bf{CMS} #it{Preliminary}")
            t.SetTextAlign(31); t.DrawLatex(0.96, 0.93, LUMI_LABEL)
            bot.cd()
            dev = max(abs(rg.GetY()[i] - 1) for rg in ratios for i in range(rg.GetN()) if 119 < rg.GetX()[i] < 131)
            lim = next(v for v in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5) if v >= min(0.5, 1.3 * dev + 0.005))  # round range: labels 1-lim, 1, 1+lim
            fb = bot.DrawFrame(XMIN, 1 - 1.3 * lim, XMAX, 1 + 1.3 * lim)  # labels at 1 +- lim stay off the frame edges
            fb.GetXaxis().SetTitle("m_{ll#gamma#gamma} [GeV]"); fb.GetYaxis().SetTitle("Ratio")
            for ax in (fb.GetXaxis(), fb.GetYaxis()):
                ax.SetTitleSize(0.055 / 0.32); ax.SetLabelSize(0.050 / 0.32)
            fb.GetXaxis().SetNdivisions(505); fb.GetXaxis().SetTitleOffset(1.05); fb.GetXaxis().SetLabelOffset(0.02)
            fb.GetYaxis().SetTitleOffset(0.48); fb.GetYaxis().SetNdivisions(502); fb.GetYaxis().CenterTitle(); fb.GetYaxis().SetLabelSize(0.042 / 0.32)
            # y labels only at 1 - lim, 1, 1 + lim, drawn by an axis inside the frame so the lowest one stays off the
            # x labels (ROOT puts its labels on the frame edges)
            fb.GetYaxis().SetLabelSize(0)
            ax = ROOT.TGaxis(XMIN, 1 - lim, XMIN, 1 + lim, 1 - lim, 1 + lim, 2, "S")
            ax.SetLabelFont(42); ax.SetLabelSize(0.042 / 0.32); ax.SetTickSize(0); ax.SetLabelOffset(0.01); ax.Draw()
            keep.append(ax)
            line = ROOT.TLine(XMIN, 1, XMAX, 1); line.SetLineStyle(3); line.Draw()
            ratios[0].Draw("L"); ratios[1].Draw("L")
            keep += [nom, up, dn, leg, fr, fb, t, top, bot, line] + ratios
        for ext in ("pdf", "png"):
            c.SaveAs("%s/sigshape_syst_mA%d.%s" % (OUT, m, ext))
        f.Close()


if __name__ == "__main__":
    main()
