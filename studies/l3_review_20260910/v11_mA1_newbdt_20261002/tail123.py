import ROOT, sys
ROOT.gROOT.SetBatch(True); ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
B="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_%s/signalFit/output/1_CMS-HGG_sigfit_%s_%s_Hm125.root"
rows=[]
for ch in ("ele","mu"):
  for era in ("2022preEE","2022postEE","2023preBPix","2023postBPix","2024"):
    f=ROOT.TFile.Open(B%(ch,era,ch)); w=f.Get("wsig_13p6TeV")
    x=w.var("CMS_hza_mass"); w.var("MH").setVal(125.)
    pdf=w.pdf("hggpdfsmrel_GG2H_%s_cat0_13p6TeV"%era)
    d=w.data("sig_mass_m125_GG2H_%s_cat0_13p6TeV"%era)
    lo,hi=115.,133.
    tot=d.sumEntries(); 
    below=d.sumEntries("CMS_hza_mass<123"); inwin=d.sumEntries("CMS_hza_mass>=115&&CMS_hza_mass<=133")
    x.setRange("win",lo,hi); x.setRange("low",lo,123.)
    pw=pdf.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.NormSet(ROOT.RooArgSet(x)),ROOT.RooFit.Range("win")).getVal()
    pl=pdf.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.NormSet(ROOT.RooArgSet(x)),ROOT.RooFit.Range("low")).getVal()
    norm=w.function("hggpdfsmrel_GG2H_%s_cat0_13p6TeV_normThisLumi"%era).getVal()
    lumi=w.var("IntLumi").getVal()
    rows.append((ch,era,tot,inwin,below/tot,pl/pw,norm,x.getMin(),x.getMax(),lumi))
T=sum(r[6] for r in rows)
for r in rows: print("%-4s %-13s dsum=%.4f inwin=%.4f MCfrac<123=%.3f PDFfrac<123=%.3f norm=%.3e share=%.3f range=[%g,%g] lumi=%g"%(r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[6]/T,r[7],r[8],r[9]))
