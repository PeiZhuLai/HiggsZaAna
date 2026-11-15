#!/usr/bin/env python3
"""Sideband goodness-of-fit table for the AN from the toys of run_gof_toys.sh.
For every mass point and envelope member: Baker-Cousins statistic S of the sideband fit to the data, mean and
standard deviation of S over the toys, the p-value (fraction of toys with S_toy >= S_obs) and the number of toys
whose fit converged. Only sideband bins enter (115-135 GeV dropped when the data are read).
Input : toys/mA<m>_<pdf>_lo<lo>_an.json, logs/members.txt
Output: /eos/home-p/pelai/HZa/output_tables/bkg_envelope/sideband_gof.{tex,csv}; summary on stdout"""
import csv, json, os
import numpy as np
W = os.path.dirname(os.path.abspath(__file__)); OUT = '/eos/home-p/pelai/HZa/output_tables/bkg_envelope'
ORDER = {'Exp': 0, 'Pow': 1, 'Lau': 2, 'Bern': 3}
rows = []
for line in open(os.path.join(W, 'logs', 'members.txt')):
    m, pdf, lo = line.split(); lo = float(lo)
    d = json.load(open(os.path.join(W, 'toys', 'mA%s_%s_lo%g_an.json' % (m, pdf, lo))))
    s = np.array([t['S'] if isinstance(t, dict) else t for t in d['toys']], dtype=float)
    rows.append(dict(mA=int(m), pdf=pdf, lo=lo, nbins=d['nbins'], npar=d['npar'], S_obs=d['S_obs'],
                     S_mean=float(s.mean()), S_sd=float(s.std(ddof=1)), p=d['p'], ntoys=len(s),
                     bound=' '.join(d.get('bound', []))))
rows.sort(key=lambda r: (r['mA'], ORDER[''.join(c for c in r['pdf'] if c.isalpha())]))
with open(os.path.join(OUT, 'sideband_gof.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def cell(r):
    return '%s & $%.1f$ & $%.1f\\pm%.1f$ & %.3f' % (r['pdf'], r['S_obs'], r['S_mean'], r['S_sd'], r['p'])


tex = []
for lo_m, hi_m in ((1, 10), (11, 20), (21, 30)):
    sub = [r for r in rows if lo_m <= r['mA'] <= hi_m]
    tex.append(r'''\begin{table}[htbp]
  \begin{center}
    \topcaption{Goodness of fit of the envelope members in the sidebands, $%d \le m_{a} \le %d~\GeV$. Each member is fitted to the data in the bins outside $115<m_{\ell\ell\gamma\gamma}<135~\GeV$ (from 95~\GeV, 100~\GeV at $m_{a}=1~\GeV$, with the parameter ranges of the fit to the data and several starting points of the turn-on position and width), and the Baker--Cousins likelihood-ratio statistic $S$ is compared with its distribution in 300 pseudo-datasets generated from the fitted function in the same bins and fitted in the same way. The $p$-value is the fraction of pseudo-datasets with a larger $S$.}
    \label{tab:sideband_gof_mA%dto%d}
    \begin{tabular}{|c|c|c|c|c|} \hline
      $m_{a}$ (\GeV) & Function & $S$ (data) & $S$ (toys, mean $\pm$ s.d.) & $p$-value \\ \hline
''' % (lo_m, hi_m, lo_m, hi_m))
    for m in range(lo_m, hi_m + 1):
        mr = [r for r in sub if r['mA'] == m]
        for i, r in enumerate(mr):
            tex.append('      %s & %s \\\\\n' % (str(m) if i == 0 else '', cell(r)))
        tex.append('      \\hline\n')
    tex.append('    \\end{tabular}\n  \\end{center}\n\\end{table}\n\n')
open(os.path.join(OUT, 'sideband_gof.tex'), 'w').write(''.join(tex))
ps = [r['p'] for r in rows]
print('%d members at %d mass points; p min %.3f max %.3f; p<0.01: %s; p<0.05: %s; min ntoys %d' % (
    len(rows), len({r['mA'] for r in rows}), min(ps), max(ps),
    [(r['mA'], r['pdf'], round(r['p'], 3)) for r in rows if r['p'] < 0.01],
    [(r['mA'], r['pdf'], round(r['p'], 3)) for r in rows if 0.01 <= r['p'] < 0.05], min(r['ntoys'] for r in rows)))
