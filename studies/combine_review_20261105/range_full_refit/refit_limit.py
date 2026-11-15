#!/usr/bin/env python3
"""What do the expected limits become if the background functions are refitted with the wider turn-on ranges?

Fit: every envelope function is refitted to the full 95-180 GeV binned data, as the production fit does
(fTest_ALP_turnOn.cpp runFit: createNLL on the whole dataset, Minuit2 migrad, strategy 1, then 2 if it fails), once
with the production ranges and once with the wider ranges of limit_wide_ranges.py (turn-on 95-125 GeV, width 0.1-50,
smearing 0.05-20 GeV). The fit uses the events in 115-135 GeV, exactly as production does; no number restricted to
that window is computed or printed. The workspaces are opened read-only.

Limits (AsymptoticLimits --run blind, i.e. Asimov from the parameter values set here):
  prod      production workspace as it is
  refitP    parameters refitted with the production ranges (closure: should reproduce prod)
  wideP     wider ranges, refitted parameters, production pdfindex
  wideB     wider ranges, refitted parameters, pdfindex re-chosen by the production rule min(2NLL + npar)

Input : Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root,
        Combine/root_t2w_xsec/<mA>_Datacard_leptons.root (both read only)
Output: results/mA<mA>.json, limits/mA<mA>.txt, combine files in work/mA<mA>/
"""
import json, os, subprocess, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'range_fit_quality'))
import fit_quality as q   # wide_range, at_bound, sideband_bins, stat (sideband-only statistic)
FF = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit'
T2W = FF + '/Combine/root_t2w_xsec/%d_Datacard_leptons.root'
LIM = '--run blind --rMax 2 --cminDefaultMinimizerStrategy 0 --rAbsAcc 1e-5 --rRelAcc 1e-5 --freezeParameters MH'.split()
IDX = 'pdfindex_cat0_13p6TeV'

def nll_val(pdf, data):
    return pdf.createNLL(data).getVal()

def one_fit(pdf, data):
    nll = pdf.createNLL(data, ROOT.RooFit.Offset(True), ROOT.RooFit.Optimize(2))
    m = ROOT.RooMinimizer(nll); m.setPrintLevel(-1); m.setStrategy(1); m.setOffsetting(True); m.optimizeConst(1)
    m.setEps(1e-6); m.setMaxFunctionCalls(20000); m.setMaxIterations(20000)
    st = m.minimize('Minuit2', 'migrad')
    if st != 0:
        m.setStrategy(2); st = m.minimize('Minuit2', 'migrad')
    return st

def best_fit(pdf, data, pars, starts):
    """Fit from every start; keep the lowest NLL among converged fits (any fit if none converged)."""
    res = []
    for s in starts:
        for v in pars: v.setVal(min(max(s[v.GetName()], v.getMin()), v.getMax()))
        st = one_fit(pdf, data); res.append((st != 0, nll_val(pdf, data), st, {v.GetName(): v.getVal() for v in pars}))
    res.sort(key=lambda r: (r[0], r[1])); b = res[0]
    for v in pars: v.setVal(b[3][v.GetName()])
    return b[1], b[2], b[3], [r[2] for r in res]

def refit(m):
    f = ROOT.TFile.Open(q.BASE % m); ws = f.Get('multipdf')
    mass = ws.var('CMS_hza_mass'); data = ws.data('roohist_data_mass_cat0')
    centers, counts = q.sideband_bins(data, mass)
    mp = ws.pdf('CMS_hgg_cat0_13p6TeV_bkgshape'); idx = ws.cat(IDX); prod_index = idx.getIndex()
    out = []
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i); pars = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
        stored = {v.GetName(): v.getVal() for v in pars}; narrow = {v.GetName(): (v.getMin(), v.getMax()) for v in pars}
        nll_s = nll_val(pdf, data); S_s = q.stat(pdf, mass, centers, counts)[0]; bound_s = [v.GetName() for v in pars if q.at_bound(v)]
        nll_p, st_p, val_p, sts_p = best_fit(pdf, data, pars, [stored])
        for v in pars:
            w = q.wide_range(v.GetName())
            if w: v.setRange(*w)
        step = [k for k in stored if 'turnon' in k or k == 'Bern_step']; soft = [k for k in stored if 'width' in k or k == 'Bern_stepWidth']
        starts = [val_p] + [dict(val_p, **{k: t for k in step}, **{k: w for k in soft}) for t in (97.0, 101.0, 105.0, 110.0, 115.0) for w in (0.5, 5.0, 20.0)]
        nll_w, st_w, val_w, sts_w = best_fit(pdf, data, pars, starts)
        S_w = q.stat(pdf, mass, centers, counts)[0]; bound_w = [v.GetName() for v in pars if q.at_bound(v)]
        wide = {v.GetName(): (v.getMin(), v.getMax()) for v in pars}
        for v in pars: v.setRange(*narrow[v.GetName()]); v.setVal(stored[v.GetName()])
        r = dict(i=i, pdf=pdf.GetName(), npar=len(pars), nll_stored=nll_s, nll_prod=nll_p, nll_wide=nll_w,
                 status_prod=st_p, status_wide=st_w, statuses_wide=sts_w, values_stored=stored, values_prod=val_p,
                 values_wide=val_w, ranges_wide=wide, bound_stored=bound_s, bound_wide=bound_w, S_sb_stored=S_s, S_sb_wide=S_w)
        out.append(r)
        print('mA%-2d %-6s st %d/%d bound %d->%d S_sb %.2f->%.2f closure 2dNLL(stored-refitP) %.3f' % (
            m, r['pdf'], st_p, st_w, len(bound_s), len(bound_w), S_s, S_w, 2*(nll_s-nll_p)), flush=True)
    f.Close()
    score = lambda r, k: 2*r[k] + r['npar']
    choice = dict(prod_index=prod_index,
                  best_stored=min(out, key=lambda r: score(r, 'nll_stored'))['i'],
                  best_refitP=min(out, key=lambda r: score(r, 'nll_prod'))['i'],
                  best_wide=min(out, key=lambda r: score(r, 'nll_wide'))['i'])
    return out, choice

def limit(m, tag, extra):
    wd = os.path.join(HERE, 'work', 'mA%d' % m); os.makedirs(wd, exist_ok=True)
    cmd = ['combine', '-M', 'AsymptoticLimits', T2W % m, '-m', '125.38', '-n', '.%s%d' % (tag, m)] + LIM + extra
    with open(os.path.join(wd, '%s.log' % tag), 'w') as lg: subprocess.run(cmd, cwd=wd, stdout=lg, stderr=subprocess.STDOUT)
    fn = os.path.join(wd, 'higgsCombine.%s%d.AsymptoticLimits.mH125.38.root' % (tag, m))
    f = ROOT.TFile.Open(fn); x = [e.limit for e in f.Get('limit')][:5]; f.Close(); return x

def setp(vals, index):
    return ['--setParameters', ','.join(['MH=125.38', '%s=%d' % (IDX, index)] + ['%s=%.10g' % kv for kv in vals.items()])]

def main():
    m = int(sys.argv[1]); os.makedirs(os.path.join(HERE, 'results'), exist_ok=True); os.makedirs(os.path.join(HERE, 'limits'), exist_ok=True)
    out, ch = refit(m)
    json.dump(dict(mA=m, choice=ch, pdfs=out), open(os.path.join(HERE, 'results', 'mA%d.json' % m), 'w'), indent=1)
    allP = {}; allW = {}; rng = []
    for r in out:
        allP.update(r['values_prod']); allW.update(r['values_wide'])
        rng += ['%s=%g,%g' % (k, a, b) for k, (a, b) in r['ranges_wide'].items() if q.wide_range(k)]
    L = dict(prod=limit(m, 'prod', ['--setParameters', 'MH=125.38']),
             refitP=limit(m, 'refitP', setp(allP, ch['prod_index'])),
             wideP=limit(m, 'wideP', setp(allW, ch['prod_index']) + ['--setParameterRanges', ':'.join(rng)]),
             wideB=limit(m, 'wideB', setp(allW, ch['best_wide']) + ['--setParameterRanges', ':'.join(rng)]))
    name = {r['i']: r['pdf'] for r in out}
    with open(os.path.join(HERE, 'limits', 'mA%d.txt' % m), 'w') as o:
        o.write('mA%d prod_index %s best_stored %s best_refitP %s best_wide %s\n' % (m, name[ch['prod_index']], name[ch['best_stored']], name[ch['best_refitP']], name[ch['best_wide']]))
        for k, v in L.items(): o.write('%-7s %s  q50/prod %.4f\n' % (k, ' '.join('%.5f' % x for x in v), v[2]/L['prod'][2]))
    print(open(os.path.join(HERE, 'limits', 'mA%d.txt' % m)).read(), flush=True)

if __name__ == '__main__': main()
