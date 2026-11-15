#!/usr/bin/env python3
"""Stored (fTest) parameter values of every envelope member, with the range and a flag when within 1% of the range
of an edge. Input: <multipdf root> ... Output: stdout"""
import sys, ROOT
ROOT.gROOT.SetBatch(True)
for fn in sys.argv[1:]:
    f = ROOT.TFile.Open(fn); ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
    mass = ws.var('CMS_hza_mass'); mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
    print('==', fn)
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i); out = []
        for v in pdf.getParameters(ROOT.RooArgSet(mass)):
            if v.isConstant() or not v.InheritsFrom('RooRealVar'): continue
            lo, hi, x = v.getMin(), v.getMax(), v.getVal(); e = 0.01*(hi - lo)
            b = ' AT_LO' if x - lo < e else (' AT_HI' if hi - x < e else '')
            if any(t in v.GetName() for t in ('turnon', 'width', 'sigma', 'step', 'Width')):
                out.append('%s=%.3g [%.3g,%.3g]%s' % (v.GetName().split('_', 1)[1], x, lo, hi, b))
        print('  %-6s %s' % (pdf.GetName(), '  '.join(out)))
