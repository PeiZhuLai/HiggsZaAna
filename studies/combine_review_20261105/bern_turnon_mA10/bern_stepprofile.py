"""Profile the NLL of an envelope member in its turn-on position (others refit from several starts); also print every member's NLL."""
import ROOT, sys, random
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
mA, pdfname, stepvar = int(sys.argv[1]), sys.argv[2], sys.argv[3]
f = ROOT.TFile.Open("flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % mA)
w = f.Get("multipdf"); data = w.data("Data_13p6TeV")
mp = [p for p in w.allPdfs() if p.ClassName() == "RooMultiPdf"][0]
for i in range(mp.getNumPdfs()):
    p = mp.getPdf(i); n = p.createNLL(data)
    k = len([v for v in p.getParameters(data) if not v.isConstant()])
    print("  %-6s NLL %.3f  npar %d" % (p.GetName(), n.getVal(), k))
pdf = w.pdf(pdfname); pars = [v for v in pdf.getParameters(data) if not v.isConstant()]
prod = {v.GetName(): v.getVal() for v in pars}
nll = pdf.createNLL(data); nll0 = nll.getVal()
step = w.var(stepvar)
def fit():
    m = ROOT.RooMinimizer(nll); m.setPrintLevel(-1); m.setStrategy(1); m.setEps(1e-6)
    m.setMaxFunctionCalls(20000); m.setMaxIterations(20000)
    return m.minimize("Minuit2", "migrad"), nll.getVal()
random.seed(2)
print("profile in %s (%s at mA%d), 2*dNLL vs production:" % (stepvar, pdfname, mA))
for s in [95.3, 96, 97, 98, 99, 100, 101, 102]:
    best = None
    for k in range(6):
        for v in pars: v.setVal(prod[v.GetName()])
        if k:
            for v in pars:
                if v.GetName() == stepvar: continue
                if "Width" in v.GetName() or "_width" in v.GetName(): v.setVal(random.choice([0.5, 1, 2, 4, 7]))
                elif "sigma" in v.GetName(): v.setVal(random.choice([0.5, 1, 2, 4]))
                elif v.GetName().startswith("Bern_b0"): v.setVal(random.uniform(0.05, 1.0))
        step.setVal(s); step.setConstant(True)
        st, val = fit()
        if st == 0 and (best is None or val < best[0]):
            best = (val, {v.GetName(): v.getVal() for v in pars})
    step.setConstant(False)
    o = best[1]
    print("  step %6.1f  2dNLL %+7.3f   " % (s, 2 * (best[0] - nll0)) + " ".join("%s=%.3g" % (n.replace("Bern_", ""), o[n]) for n in sorted(o) if n != stepvar))
