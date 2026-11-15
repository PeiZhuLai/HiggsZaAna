#!/usr/bin/env python3
"""Lists, for every member of the production envelope at the 30 mass points, the floating shape parameters whose
stored best-fit value is at (within 0.1% of the range from) an edge of its range, plus the range of the turn-on,
width and smearing parameters. Read-only.
Input : flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/<m>/CMS-HGG_mva_13p6TeV_multipdf.root
Output: stdout (redirect to params_at_limit.txt)"""
import ROOT
ROOT.gROOT.SetBatch(True)
FR = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3'
tot = nb = 0; per_kind = {}
for m in range(1, 31):
    f = ROOT.TFile.Open('%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root' % (FR, m))
    ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
    mass = ws.var('CMS_hza_mass')
    mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i); name = pdf.GetName().split('_')[-1] if False else pdf.GetName()
        at = []; rng = []
        for v in pdf.getParameters(ROOT.RooArgSet(mass)):
            if v.isConstant() or not v.InheritsFrom('RooRealVar'):
                continue
            n = v.GetName(); lo, hi, x = v.getMin(), v.getMax(), v.getVal(); tol = 1e-3 * (hi - lo)
            tot += 1
            kind = 'turnon' if ('turnon' in n or n == 'Bern_step') else 'width' if ('width' in n or 'stepWidth' in n) \
                else 'sigma' if ('sigma' in n) else 'other'
            per_kind.setdefault(kind, [0, 0]); per_kind[kind][0] += 1
            if kind != 'other':
                rng.append('%s[%g,%g]' % (n, lo, hi))
            if x - lo < tol or hi - x < tol:
                at.append('%s=%.4g(%s)' % (n, x, 'lo' if x - lo < tol else 'hi')); nb += 1; per_kind[kind][1] += 1
        print('mA%-2d %-28s at_limit: %-60s ranges: %s' % (m, pdf.GetName(), ' '.join(at) or '-', ' '.join(rng)))
    f.Close()
print('TOTAL floating shape parameters %d, at a limit %d' % (tot, nb))
for k, (n, b) in sorted(per_kind.items()):
    print('  %-7s %3d of %3d at a limit' % (k, b, n))
