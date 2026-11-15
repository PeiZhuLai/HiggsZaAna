#!/usr/bin/env python3
"""Fixed-index impacts (IMPACT_FIX_INDEX=1, production ranges) against the free-index ones of 11-06, r = q50.
Input : Combine/output_impacts{_xsec,}_q50{,_fixidx_prodrange}/<m>_impacts.json, logs/impact_fixidx/*.out (chosen index)
Output: compare_fixidx.txt
Columns: sigma_r (half width of the 68% interval) free and fixed; chosen function; top three impacts in units of the
fixed sigma_r; background shape parameters of the chosen function with an impact, and how many have a one-sided
interval (best fit at a range limit, i.e. the fit value equals one end of its interval)."""
import glob, json, os, re
C = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine'
H = os.path.dirname(os.path.abspath(__file__)); L = os.path.join(H, '..', 'logs', 'impact_fixidx')
BKG = re.compile(r'^(Exp_|Pow_|Lau_|Bern_)')
def load(p):
    return json.load(open(p)) if os.path.exists(p) else None
def sig(d):
    r = [p for p in d['POIs'] if p['name'] == 'r'][0]['fit']; return 0.5*(r[2]-r[0])
def imp(p): return max(abs(p['r'][2]-p['r'][1]), abs(p['r'][1]-p['r'][0]))
out = ['card   mA  sig_free  sig_fix  ratio | chosen | top 3 (impact/sigma_fix) | bkg shape: with impact / one-sided | nparams']
for card, tag in (('xsec', '_xsec'), ('wilson', '')):
    for m in range(1, 31):
        fr = load('%s/output_impacts%s_q50/%d_impacts.json' % (C, tag, m)); fx = load('%s/output_impacts%s_q50_fixidx_prodrange/%d_impacts.json' % (C, tag, m))
        if not fx: out.append('%-6s %-3d MISSING' % (card, m)); continue
        chosen = '?'
        for lg in glob.glob(os.path.join(L, '%s_mA%d.*.out' % (card, m))):
            for l in open(lg):
                if 'chosen pdfindex' in l: chosen = l.split('=')[1].split(';')[0].strip()
        s = sig(fx); P = sorted(fx['params'], key=lambda p: -imp(p))
        top = ', '.join('%s %.2f' % (p['name'].replace('CMS_hza_', '').replace('shapeBkg_bkg_mass_cat0__norm', 'bkg_norm'), imp(p)/s) for p in P[:3])
        bk = [p for p in fx['params'] if BKG.match(p['name'])]
        withimp = [p for p in bk if imp(p) > 0]; one = [p for p in bk if abs(p['fit'][1]-p['fit'][0]) < 1e-9 or abs(p['fit'][2]-p['fit'][1]) < 1e-9]
        sf = sig(fr) if fr else float('nan')
        out.append('%-6s %-3d %.5f  %.5f  %.3f | %-2s | %s | %d/%d / %d | %d' % (card, m, sf, s, s/sf, chosen, top, len(withimp), len(bk), len(one), len(fx['params'])))
open(os.path.join(H, 'compare_fixidx.txt'), 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
