"""Refit one envelope member of an mA multipdf workspace from many starting points (same data, range, NLL as fTest)."""
import ROOT, sys, itertools, random
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
mA, pdfname = int(sys.argv[1]), sys.argv[2]
f = ROOT.TFile.Open("flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % mA)
w = f.Get("multipdf"); data = w.data("Data_13p6TeV"); pdf = w.pdf(pdfname); mass = w.var("CMS_hza_mass")
pars = [v for v in pdf.getParameters(data) if not v.isConstant()]
prod = {v.GetName(): v.getVal() for v in pars}
nll = pdf.createNLL(data)
def fit():
    m = ROOT.RooMinimizer(nll); m.setPrintLevel(-1); m.setStrategy(1); m.setEps(1e-6)
    m.setMaxFunctionCalls(20000); m.setMaxIterations(20000)
    st = m.minimize("Minuit2", "migrad")
    return st, nll.getVal()
def setp(d):
    for v in pars:
        if v.GetName() in d: v.setVal(d[v.GetName()])
setp(prod); nll0 = nll.getVal()
st, nllp = fit(); refit_prod = {v.GetName(): v.getVal() for v in pars}
print("production point NLL %.4f ; refit from it: status %d NLL %.4f" % (nll0, st, nllp))
fam = pdfname.rstrip("0123456789")
stepn = {"Bern": ("Bern_step", "Bern_stepWidth", "Bern_gsigma"), "Pow": ("Pow_turnon_p1", "Pow_width_p1", "Pow_sigma_p1")}[fam]
random.seed(1)
best = (nllp, refit_prod)
rows = []
for s, wd, sg in itertools.product([97, 99, 101, 103, 105], [0.5, 2, 4, 7], [0.5, 2, 4]):
    for k in range(2):
        start = dict(prod); start[stepn[0]] = s; start[stepn[1]] = wd; start[stepn[2]] = sg
        if k == 1:
            for n in start:
                if n.startswith("Bern_b0"): start[n] = random.uniform(0.05, 1.0)
        setp(start)
        st, v = fit()
        vals = {p.GetName(): p.getVal() for p in pars}
        rows.append((v, st, s, wd, sg, k, vals))
        if st == 0 and v < best[0]: best = (v, vals)
rows.sort(key=lambda r: r[0])
print("lowest 8 minima (status, start step/width/sigma/randcoef -> dNLL vs production, fitted step/width/sigma):")
for v, st, s, wd, sg, k, vals in rows[:8]:
    print("  st%d  start %5.1f %4.1f %4.1f %d -> dNLL %+8.4f   %7.2f %6.2f %6.2f" % (st, s, wd, sg, k, v - nll0, vals[stepn[0]], vals[stepn[1]], vals[stepn[2]]))
import collections
c = collections.Counter(round(r[0] - nll0, 2) for r in rows if r[1] == 0)
print("distinct converged minima (dNLL: count):", sorted(c.items())[:12])
print("best:", " ".join("%s=%.4g" % kv for kv in sorted(best[1].items())), " dNLL %+.4f" % (best[0] - nll0))
