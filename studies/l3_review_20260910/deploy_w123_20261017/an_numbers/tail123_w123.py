#!/usr/bin/env python3
"""m_a = 1 signal model vs simulation below 123 GeV inside the 115-133 GeV fit window, per era and channel
(Signal/outdir_<ch>/signalFit/output/1_CMS-HGG_sigfit_<era>_<ch>_Hm125.root: dataset = in-window simulation,
pdf hggpdfsmrel_* at its stored MH). Output: an_numbers/tail123_w123.txt"""
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
F = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
out = []; worst = (0, None)
for ch in ("ele", "mu"):
    for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        f = ROOT.TFile.Open("%s/Signal/outdir_%s/signalFit/output/1_CMS-HGG_sigfit_%s_%s_Hm125.root" % (F, ch, e, ch))
        w = f.Get("wsig_13p6TeV"); x = w.var("CMS_hza_mass")
        d = w.data("sig_mass_m125_GG2H_%s_cat0_13p6TeV" % e); pdf = w.pdf("hggpdfsmrel_GG2H_%s_cat0_13p6TeV" % e)
        tot = d.sumEntries("CMS_hza_mass>=115 && CMS_hza_mass<=133"); low = d.sumEntries("CMS_hza_mass>=115 && CMS_hza_mass<123")
        x.setRange("win", 115, 133); x.setRange("low", 115, 123)
        ip = pdf.createIntegral(ROOT.RooArgSet(x), ROOT.RooFit.NormSet(ROOT.RooArgSet(x)), ROOT.RooFit.Range("low")).getVal() / \
             pdf.createIntegral(ROOT.RooArgSet(x), ROOT.RooFit.NormSet(ROOT.RooArgSet(x)), ROOT.RooFit.Range("win")).getVal()
        mc = low / tot
        out.append("%-4s %-13s MCfrac<123=%.3f PDFfrac<123=%.3f diff=%+.1f points" % (ch, e, mc, ip, 100 * (ip - mc)))
        if abs(ip - mc) > abs(worst[0]): worst = (ip - mc, (ch, e, mc, ip))
        f.Close()
ch, e, mc, ip = worst[1]
out.append("largest |diff| %.1f points: %s %s model %.1f%% vs simulation %.1f%%" % (100 * abs(worst[0]), ch, e, 100 * ip, 100 * mc))
open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017/an_numbers/tail123_w123.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
