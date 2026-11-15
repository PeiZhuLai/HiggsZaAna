#!/usr/bin/env python3
"""Every param nuisance of a card must act on the model: in the text2workspace output it must be a variable of the
workspace, and moving it from 0 to +1 must change the signal+background model pdf (model_s) at some m_llgg.
A param whose name matches no workspace variable is created free-standing by text2workspace (only its Gaussian
constraint), which combine accepts silently.
Usage : check_params.py <card.txt> <t2w.root>
Output: one line per param, then 'ALL OK' or 'PROBLEMS n'"""
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
card, ws = sys.argv[1], sys.argv[2]
params = [l.split()[0] for l in open(card) if len(l.split()) > 1 and l.split()[1] == 'param']
f = ROOT.TFile.Open(ws); w = f.Get('w'); mod = w.pdf('model_s'); m = w.var('CMS_hza_mass')
obs = ROOT.RooArgSet(m); bad = 0
def vals():
    out = []
    for x in (110., 118., 122., 125., 128., 132.):
        m.setVal(x); out.append(mod.getVal(obs))
    return out
r = w.var('r'); r.setVal(1.)  # signal on, so the signal shape enters model_s
base = vals()
for p in params:
    v = w.var(p)
    if not v: print('[ERROR] %s: not in the workspace' % p); bad += 1; continue
    v0 = v.getVal(); v.setVal(v0 + 1.); d = max(abs(a - b) / max(abs(b), 1e-300) for a, b in zip(vals(), base)); v.setVal(v0)
    ok = d > 1e-9
    print('%s %-50s max rel change of model_s at +1: %.2e' % ('[OK]   ' if ok else '[ERROR]', p, d)); bad += 0 if ok else 1
print('ALL OK' if not bad else 'PROBLEMS %d' % bad)
