#!/usr/bin/env python3
"""Staging (opened ranges) vs production envelope, per mass: members, default index, GOF, every shape parameter
(value and range) and how many sit at a range limit (within 1e-3 of the range width).
mA19 is the regression point (ranges unchanged): its parameters must agree with production.
Input : fit_results/<m>/ (staging), flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/<m>/ (production)
Output: compare.txt"""
import os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
H = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get('RES', 'fit_results')  # staging result directory, e.g. fit_results_ms
P = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3'
masses = [int(a) for a in sys.argv[1:]] or [19, 10, 15, 21, 22, 23]

def env(d):
    out = {}
    p = os.path.join(d, 'HZAmassInde_fTest', 'EnvelopeResults.txt')
    if os.path.exists(p):
        for l in open(p):
            if 'pdf :' in l:
                f = [x.strip() for x in l.split(',')]
                out[f[1].split(':')[1].strip()] = (float(f[2].split(':')[1]), f[3].split(':')[1].strip())
    return out

def load(d):
    f = ROOT.TFile.Open(os.path.join(d, 'CMS-HGG_mva_13p6TeV_multipdf.root'))
    if not f or f.IsZombie(): return None
    w = f.Get('multipdf'); m = [p for p in w.allPdfs() if p.InheritsFrom('RooMultiPdf')][0]
    obs = w.var('CMS_hza_mass') or w.var('CMS_hgg_mass')
    idx = [c for c in w.allCats()][0].getIndex(); pdfs = {}
    for i in range(m.getNumPdfs()):
        pdf = m.getPdf(i); pars = {}
        for v in pdf.getParameters(ROOT.RooArgSet(obs)):
            if v.InheritsFrom('RooRealVar') and not v.isConstant():
                pars[v.GetName()] = (v.getVal(), v.getMin(), v.getMax())
        pdfs[pdf.GetName()] = pars
    names = [m.getPdf(i).GetName() for i in range(m.getNumPdfs())]
    f.Close()
    return names, idx, pdfs

def atlim(v, lo, hi):
    e = 1e-3*(hi-lo); return abs(v-lo) < e or abs(hi-v) < e

out = []
for mA in masses:
    s, p = load(os.path.join(H, RES, str(mA))), load(os.path.join(P, str(mA)))
    es, ep = env(os.path.join(H, RES, str(mA))), env(os.path.join(P, str(mA)))
    out.append('=== mA%d' % mA)
    if not s: out.append('  staging MISSING'); continue
    out.append('  production: %s default %s | staging: %s default %s' % (','.join(p[0]), p[0][p[1]], ','.join(s[0]), s[0][s[1]]))
    out.append('  GOF (truth flag) production: %s' % ', '.join('%s %.2f(%s)' % (k, v[0], v[1]) for k, v in ep.items()))
    out.append('  GOF (truth flag) staging   : %s' % ', '.join('%s %.2f(%s)' % (k, v[0], v[1]) for k, v in es.items()))
    nlp = sum(atlim(*t) for d in p[2].values() for t in d.values()); nls = sum(atlim(*t) for d in s[2].values() for t in d.values())
    out.append('  params at a range limit: production %d, staging %d' % (nlp, nls))
    maxdiff = 0.0
    for name in s[0]:
        for k, (v, lo, hi) in sorted(s[2][name].items()):
            pv = p[2].get(name, {}).get(k)
            ptxt = '%.4g [%.4g, %.4g]%s' % (pv[0], pv[1], pv[2], ' *' if atlim(*pv) else '') if pv else '-'
            if pv: maxdiff = max(maxdiff, abs(v - pv[0])/max(1e-9, pv[2]-pv[1]))
            out.append('    %-6s %-16s prod %-32s stg %.4g [%.4g, %.4g]%s' % (name, k, ptxt, v, lo, hi, ' *' if atlim(v, lo, hi) else ''))
    out.append('  max |staging - production| / production range width: %.2e' % maxdiff)
txt = '\n'.join(out) + '\n'
open(os.path.join(H, 'compare.txt' if RES == 'fit_results' else 'compare_%s.txt' % RES), 'w').write(txt); print(txt)
