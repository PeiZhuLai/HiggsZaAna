#!/usr/bin/env python3
"""Is the staging fTest minimum the global one? NLL of each staging envelope member at its stored values, on
roohist_data_mass_cat0 of the staging multipdf (createNLL on the whole dataset, as fTest runFit and
range_full_refit/refit_limit.py), against nll_wide of range_full_refit (best of 16 starts, same opened ranges).
Input : fit_results/<m>/CMS-HGG_mva_13p6TeV_multipdf.root, ../range_full_refit/results/mA<m>.json
Output: nll_check.txt"""
import json, os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
H = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get('RES', 'fit_results')  # staging result directory, e.g. fit_results_ms
out = ['mA  pdf    NLL staging  NLL 16-start(wide)  diff(staging-16start)']
for m in [int(a) for a in sys.argv[1:]] or [10, 15, 21, 22, 23]:
    f = ROOT.TFile.Open(os.path.join(H, RES, str(m), 'CMS-HGG_mva_13p6TeV_multipdf.root'))
    w = f.Get('multipdf'); data = w.data('roohist_data_mass_cat0')
    mp = [p for p in w.allPdfs() if p.InheritsFrom('RooMultiPdf')][0]
    ref = {p['pdf']: p['nll_wide'] for p in json.load(open(os.path.join(H, '..', 'range_full_refit', 'results', 'mA%d.json' % m)))['pdfs']}
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i); v = pdf.createNLL(data).getVal(); r = ref.get(pdf.GetName())
        out.append('%-3d %-6s %.3f   %s   %s' % (m, pdf.GetName(), v, '%.3f' % r if r else '-', '%+.3f' % (v - r) if r else '-'))
    f.Close()
txt = '\n'.join(out) + '\n'; open(os.path.join(H, 'nll_check.txt' if RES == 'fit_results' else 'nll_check_%s.txt' % RES), 'w').write(txt); print(txt)
