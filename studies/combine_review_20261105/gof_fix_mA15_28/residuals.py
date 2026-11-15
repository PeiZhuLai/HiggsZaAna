#!/usr/bin/env python3
"""Where does the sideband fit of mA28 / mA15 fail? Fit every envelope member of the deployed multipdf on the sideband
bins (production ranges, multi-start as in the deployment) and print the per-bin Pearson residuals (n - mu)/sqrt(mu),
plus the residuals summed in 5 GeV groups.

The bins in 115-135 GeV are dropped when the data are read (fit_quality.sideband_bins); nothing here uses them.
Input : Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root (read-only)
Output: stdout only
"""
import sys, math, ROOT
sys.path.insert(0, '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_fit_quality')
import fit_quality as fq

def main(m):
    f = ROOT.TFile.Open(fq.BASE % m)
    ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
    mass = ws.var('CMS_hza_mass')
    centers, counts = fq.sideband_bins(ws.data('roohist_data_mass_cat0'), mass)
    assert not any(fq.BLIND[0] < c < fq.BLIND[1] for c in centers)
    mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
    print('mA%d  sideband bins %d  events %d' % (m, len(centers), sum(counts)))
    for i in range(mp.getNumPdfs()):
        pdf = mp.getPdf(i)
        pars = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
        stored = {v.GetName(): v.getVal() for v in pars}
        step = [k for k in stored if 'turnon' in k or k == 'Bern_step']; soft = [k for k in stored if 'width' in k or k == 'Bern_stepWidth']
        starts = [stored] + [dict(stored, **{k: t for k in step}, **{k: w for k in soft}) for t in (97.0, 101.0, 105.0, 110.0) for w in (0.5, 5.0, 20.0)]
        (S, P, st, val), _ = fq.best_fit(pdf, mass, centers, counts, pars, starts)
        f_ = fq.shape_at(pdf, mass, centers); tot = sum(counts); fs = sum(f_)
        mus = [tot*x/fs for x in f_]
        r = [(n - mu)/math.sqrt(mu) for n, mu in zip(counts, mus)]
        print('\n%-6s npar %d  S %.2f  Pearson %.2f  status %d  pars %s' % (pdf.GetName(), len(pars), S, P, st,
              ' '.join('%s=%.3g' % (k.split('_')[-2] + '_' + k.split('_')[-1] if k.count('_') >= 2 else k, v) for k, v in val.items())))
        # 5 GeV groups: sum of n, sum of mu, (N-M)/sqrt(M)
        groups = {}
        for c, n, mu in zip(centers, counts, mus):
            g = int((c - 95) // 5); groups.setdefault(g, [0, 0.0]); groups[g][0] += n; groups[g][1] += mu
        print('  5GeV groups (lo: N/M pull):', '  '.join('%d:%d/%.1f %+.1f' % (95 + 5*g, N, M, (N - M)/math.sqrt(M)) for g, (N, M) in sorted(groups.items())))
        big = [(c, n, mu, x) for c, n, mu, x in zip(centers, counts, mus, r) if abs(x) > 2.0]
        print('  |pull|>2 bins:', '  '.join('%.1f:%d/%.1f(%+.1f)' % b for b in big))
        # runs test on the sign of the residuals
        signs = [1 if x > 0 else -1 for x in r]; runs = 1 + sum(1 for a, b in zip(signs, signs[1:]) if a != b)
        npos = signs.count(1); nneg = signs.count(-1); n = npos + nneg
        er = 1 + 2.0*npos*nneg/n; vr = 2.0*npos*nneg*(2*npos*nneg - n)/(n*n*(n - 1))
        print('  runs %d (expected %.1f +- %.1f)  sum pull^2 low(<115) %.1f high(>135) %.1f' % (runs, er, math.sqrt(vr),
              sum(x*x for c, x in zip(centers, r) if c < 115), sum(x*x for c, x in zip(centers, r) if c > 135)))

if __name__ == '__main__':
    for m in map(int, sys.argv[1:]): main(m)
