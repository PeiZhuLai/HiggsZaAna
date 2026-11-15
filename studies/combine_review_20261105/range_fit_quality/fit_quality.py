#!/usr/bin/env python3
"""Do the sideband data want the wider turn-on ranges? Refit every envelope function with the production ranges and
with the ranges of limit_wide_ranges.py (turn-on 95-125, width 0.1-50, smearing 0.05-20 GeV), sidebands only.

The bins in 115-135 GeV are skipped when the data are read: no fit and no number here uses them. Only the shape
parameters of the function are floated, as in the production fit; the normalization follows the sideband yield.
Workspaces are opened read-only; nothing is written back.

Input : Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root
Output: <outdir>/mA<mA>.json, one record per function
  S        Baker-Cousins statistic on the sideband bins, 2*sum[mu - n + n*ln(n/mu)], minimized directly; mu = pdf shape
           scaled to the sideband yield, so a difference of S is 2*delta(NLL)
  dS       S(production ranges) - S(wider ranges); >= 0 when the wider ranges contain the production ones
  sr_ratio background predicted in 115-135 GeV per sideband event, wider over production ranges (uses the pdf only)
"""
import argparse, json, math, os, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
BASE = ('/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/'
        'ALP_BkgModel_ReReco/fit_results_run3/%d/CMS-HGG_mva_13p6TeV_multipdf.root')
BLIND = (115.0, 135.0)

def wide_range(name):
    if 'turnon' in name or name == 'Bern_step': return (95.0, 125.0)
    if 'width' in name or name == 'Bern_stepWidth': return (0.1, 50.0)
    if 'sigma' in name: return (0.05, 20.0)
    return None

def at_bound(v):
    span = v.getMax() - v.getMin()
    return v.getVal() - v.getMin() < 1e-3*span or v.getMax() - v.getVal() < 1e-3*span

def sideband_bins(hist, mass):
    """Bin centers and contents of the binned data outside the blinded window."""
    centers, counts = [], []
    for i in range(hist.numEntries()):
        c = hist.get(i).getRealValue(mass.GetName()); n = hist.weight()
        if BLIND[0] < c < BLIND[1]: continue
        centers.append(c); counts.append(n)
    return centers, counts

def shape_at(pdf, mass, centers):
    ns = ROOT.RooArgSet(mass); vals = []
    for c in centers:
        mass.setVal(c); vals.append(pdf.getVal(ns))
    return vals

def stat(pdf, mass, centers, counts):
    f = shape_at(pdf, mass, centers); tot = sum(counts); fs = sum(f)
    if fs <= 0: return float('inf'), float('inf')
    S = P = 0.0
    for fi, n in zip(f, counts):
        mu = tot*fi/fs
        if mu <= 0: return float('inf'), float('inf')
        S += 2*(mu - n + (n*math.log(n/mu) if n > 0 else 0.0)); P += (n - mu)**2/mu
    return S, P

def sr_per_sideband(pdf, mass, lo, hi):
    ns = ROOT.RooArgSet(mass); sr = sb = 0.0; c = lo + 0.5
    while c < hi:
        mass.setVal(c); v = pdf.getVal(ns)
        if BLIND[0] < c < BLIND[1]: sr += v
        else: sb += v
        c += 1.0
    return sr/sb if sb > 0 else float('nan')

def minimize(pdf, mass, centers, counts, pars):
    """Minuit2 on S itself. RooFit's Range("a,b") is not usable here: for RooFFTConvPdf its NLL stays normalized over
    the whole mass window, which rewards a pdf that empties the blinded region (checked at mA10, Pow1)."""
    def fcn(x):
        for i, v in enumerate(pars): v.setVal(x[i])
        S = stat(pdf, mass, centers, counts)[0]
        return S if S < 1e9 else 1e9
    mini = ROOT.Math.Factory.CreateMinimizer('Minuit2', 'Migrad'); func = ROOT.Math.Functor(fcn, len(pars))
    mini.SetFunction(func); mini.SetPrintLevel(-1); mini.SetStrategy(1); mini.SetTolerance(0.01); mini.SetErrorDef(1.0); mini.SetMaxFunctionCalls(5000)
    for i, v in enumerate(pars):
        span = v.getMax() - v.getMin(); mini.SetLimitedVariable(i, v.GetName(), v.getVal(), 0.02*span, v.getMin(), v.getMax())
    mini.Minimize(); xs = mini.X()
    for i, v in enumerate(pars): v.setVal(xs[i])
    return mini.Status()

def best_fit(pdf, mass, centers, counts, pars, starts):
    """Fit from each start, keep the lowest S; the parameters are left at that fit."""
    best = None; statuses = []
    for s in starts:
        for v in pars: v.setVal(min(max(s[v.GetName()], v.getMin()), v.getMax()))
        st = minimize(pdf, mass, centers, counts, pars); S, P = stat(pdf, mass, centers, counts); statuses.append(st)
        if best is None or S < best[0]:
            best = (S, P, st, {v.GetName(): v.getVal() for v in pars})
    for v in pars: v.setVal(best[3][v.GetName()])
    return best, statuses

def fit_function(pdf, mass, centers, counts, lo, hi):
    pars = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
    stored = {v.GetName(): v.getVal() for v in pars}; narrow = {v.GetName(): (v.getMin(), v.getMax()) for v in pars}
    mid = {k: 0.5*(a+b) for k, (a, b) in narrow.items()}
    (S0, P0, st0, val0), sts0 = best_fit(pdf, mass, centers, counts, pars, [stored, dict(stored, **{k: mid[k] for k in mid if wide_range(k)})])
    bound0 = [v.GetName() for v in pars if at_bound(v)]; sr0 = sr_per_sideband(pdf, mass, lo, hi)
    changed = False
    for v in pars:
        w = wide_range(v.GetName())
        if w and (v.getMin(), v.getMax()) != w: v.setRange(*w); changed = True
    if changed:
        step = [k for k in stored if 'turnon' in k or k == 'Bern_step']; soft = [k for k in stored if 'width' in k or k == 'Bern_stepWidth']
        starts = [val0] + [dict(val0, **{k: t for k in step}, **{k: w for k in soft}) for t in (97.0, 101.0, 105.0, 110.0, 115.0) for w in (0.5, 5.0, 20.0)]
        (S1, P1, st1, val1), sts1 = best_fit(pdf, mass, centers, counts, pars, starts)
    else:
        S1, P1, st1, val1, sts1 = S0, P0, st0, val0, [st0]
    bound1 = [v.GetName() for v in pars if at_bound(v)]; sr1 = sr_per_sideband(pdf, mass, lo, hi)
    nested = all(wide_range(k) is None or (wide_range(k)[0] <= a and b <= wide_range(k)[1]) for k, (a, b) in narrow.items())
    for v in pars:
        v.setRange(*narrow[v.GetName()]); v.setVal(stored[v.GetName()])
    return dict(pdf=pdf.GetName(), npar=len(pars), opened=changed, nested=nested,
                S_prod=S0, P_prod=P0, status_prod=st0, statuses_prod=sts0, bound_prod=bound0, values_prod=val0,
                S_wide=S1, P_wide=P1, status_wide=st1, statuses_wide=sts1, bound_wide=bound1, values_wide=val1,
                dS=S0-S1, sr_ratio=sr1/sr0 if sr0 > 0 else float('nan'))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('mass', type=int); ap.add_argument('outdir'); ap.add_argument('--pdf', default='')
    a = ap.parse_args()
    f = ROOT.TFile.Open(BASE % a.mass)
    ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
    mass = ws.var('CMS_hza_mass'); lo, hi = mass.getMin(), mass.getMax()
    centers, counts = sideband_bins(ws.data('roohist_data_mass_cat0'), mass)
    assert not any(BLIND[0] < c < BLIND[1] for c in centers)
    mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]; out = []
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i)
        if a.pdf and pdf.GetName() != a.pdf: continue
        r = fit_function(pdf, mass, centers, counts, lo, hi); r.update(mA=a.mass, nbins=len(centers)); out.append(r)
        print('mA%-2d %-6s S %.2f -> %.2f dS %.3f st %d/%d bound %d->%d sr_ratio %.3f' % (a.mass, r['pdf'], r['S_prod'], r['S_wide'], r['dS'], r['status_prod'], r['status_wide'], len(r['bound_prod']), len(r['bound_wide']), r['sr_ratio']), flush=True)
    os.makedirs(a.outdir, exist_ok=True); json.dump(out, open(os.path.join(a.outdir, 'mA%d.json' % a.mass), 'w'), indent=1)

if __name__ == '__main__': main()
