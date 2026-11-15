#!/usr/bin/env python3
"""Truth models for the bias screen of one staging envelope: every envelope member at its stored (fTest) values,
plus the wide-truth model of bias_wide_truth/truth.json (best of 16 starts with the opened ranges, picked by
min(2NLL + npar)) when that function is in the envelope. Injection r = 0 and q50 of q50_at_bias_submission.txt.
Usage : make_truths.py <m> <tag>
Input : cards/<tag>/<m>_Datacard_leptons.root, ../../bias_wide_truth/truth.json, ../../q50_at_bias_submission.txt
Output: cards/<tag>/<m>_truths.json  {key: {pdf, index, set, q50}}"""
import json, os, sys, ROOT
ROOT.gROOT.SetBatch(True)
m, tag = int(sys.argv[1]), sys.argv[2]
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, '..', '..')
q50 = {int(l.split()[0]): float(l.split()[1]) for l in open(os.path.join(R, 'q50_at_bias_submission.txt')) if l.strip()}
f = ROOT.TFile.Open(os.path.join(H, 'cards', tag, '%d_Datacard_leptons.root' % m)); w = f.Get('w')
mp = [p for p in w.allPdfs() if p.InheritsFrom('RooMultiPdf')][0]; obs = w.var('CMS_hza_mass')
out = {}
for i in range(mp.getNumPdfs()):
    pdf = mp.getPdf(i)
    pars = [v for v in pdf.getParameters(ROOT.RooArgSet(obs)) if v.InheritsFrom('RooRealVar') and not v.isConstant()]
    out[pdf.GetName()] = dict(pdf=pdf.GetName(), index=i, q50=q50[m], set=','.join('%s=%.10g' % (v.GetName(), v.getVal()) for v in pars))
wt = json.load(open(os.path.join(R, 'bias_wide_truth', 'truth.json'))).get(str(m))
names = [mp.getPdf(i).GetName() for i in range(mp.getNumPdfs())]
if wt and wt['pdf'] in names:
    out['wide' + wt['pdf']] = dict(pdf=wt['pdf'], index=names.index(wt['pdf']), q50=q50[m], set=wt['set'])
json.dump(out, open(os.path.join(H, 'cards', tag, '%d_truths.json' % m), 'w'), indent=1)
for k, v in out.items(): print(k, v['index'], v['set'])
