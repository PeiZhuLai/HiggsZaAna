"""Recompute GOF for the envelope members with a Poisson likelihood-ratio (Baker-Cousins) chi2
and occupancy stats; flag parameters at their range limits. Reads the multipdf workspaces only."""
import ROOT, math, sys
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
N = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3"
O = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/snapshot_preFSRfix_20260924/fit_results_run3"
for tag, base in (("new", N), ("old", O)):
    for m in [int(x) for x in sys.argv[1:]]:
        f = ROOT.TFile.Open("%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (base, m)); w = f.Get("multipdf")
        x = w.var("CMS_hza_mass"); d = w.data("roohist_data_mass_cat0")
        mp = [p for p in w.allPdfs() if p.InheritsFrom("RooMultiPdf")][0]
        h = d.createHistogram("h%s%d" % (tag, m), x, ROOT.RooFit.Binning(85, 95, 180))
        ntot = h.Integral(); low = sum(1 for i in range(1, 86) if h.GetBinContent(i) < 5)
        print("%s mA%d N=%d bins<5: %d/85  mean/bin %.2f" % (tag, m, ntot, low, ntot / 85.))
        for i in range(mp.getNumPdfs()):
            pdf = mp.getPdf(i)
            pars = pdf.getParameters(ROOT.RooArgSet(x)); npar = pars.getSize() + 1
            rail = [p.GetName() for p in pars if p.InheritsFrom("RooRealVar") and not p.isConstant()
                    and (abs(p.getVal() - p.getMin()) < 1e-3 * max(1, abs(p.getMax() - p.getMin())) or abs(p.getVal() - p.getMax()) < 1e-3 * max(1, abs(p.getMax() - p.getMin())))]
            bc = 0.; ney = 0.
            for b in range(1, 86):
                lo, hi = h.GetBinLowEdge(b), h.GetBinLowEdge(b) + h.GetBinWidth(b)
                x.setRange("b", lo, hi)
                mu = ntot * pdf.createIntegral(ROOT.RooArgSet(x), ROOT.RooFit.NormSet(ROOT.RooArgSet(x)), ROOT.RooFit.Range("b")).getVal()
                n = h.GetBinContent(b)
                bc += 2 * (mu - n + (n * math.log(n / mu) if n > 0 else 0))
                if n > 0: ney += (n - mu) ** 2 / n
            ndf = 85 - npar
            print("   %-6s npar=%d  BakerCousins chi2=%.1f p=%.3g | Neyman(sqrt n) chi2=%.1f p=%.3g %s" % (
                pdf.GetName(), npar, bc, ROOT.TMath.Prob(bc, ndf), ney, ROOT.TMath.Prob(ney, ndf), ("RAIL:" + ",".join(rail)) if rail else ""))
