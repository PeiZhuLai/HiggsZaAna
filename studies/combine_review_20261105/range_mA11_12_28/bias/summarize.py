#!/usr/bin/env python3
"""Pulls of the bias screen, same definition as bias_wide_truth/summarize.py and RunBiasStudy.py:
(r_hat - r_inj) / sigma with sigma on the side facing the injected value; toys without a valid interval dropped.
Gaussian-fit mean on a TH1 with 80 bins in [-4, 4] (AN convention), sample mean and median.
Usage : summarize.py <tag> [<tag> ...]
Input : out/<tag>/mA<m>_<truth>_<inj>/chunk*.root, cards/<tag>/<m>_truths.json
Output: summary_<tag>.txt"""
import glob, json, math, os, sys, ROOT
ROOT.gROOT.SetBatch(True)
H = os.path.dirname(os.path.abspath(__file__))
for tag in sys.argv[1:]:
    out = ['mA  truth      inj  r_inj   | ntoys  mean_pull +- err   gaus_mean +- err   median | chunks']
    for tj in sorted(glob.glob(os.path.join(H, 'cards', tag, '*_truths.json'))):
        m = os.path.basename(tj).split('_')[0]; truths = json.load(open(tj))
        for tr in truths:
            for inj in ('0', 'q50'):
                r_inj = truths[tr]['q50'] if inj == 'q50' else 0.0
                files = sorted(glob.glob(os.path.join(H, 'out', tag, 'mA%s_%s_%s' % (m, tr, inj), 'chunk*.root')))
                pulls = []
                for fn in files:
                    f = ROOT.TFile.Open(fn); t = f.Get('limit')
                    for i in range(t.GetEntries() // 3):
                        t.GetEntry(3*i); bf, q0 = t.r, t.quantileExpected
                        t.GetEntry(3*i+1); lo, q1 = t.r, t.quantileExpected
                        t.GetEntry(3*i+2); hi, q2 = t.r, t.quantileExpected
                        if q0 != -1 or abs(q1+0.32) > 1e-3 or abs(q2-0.32) > 1e-3: continue
                        u = (hi - bf) if bf < r_inj else (bf - lo)
                        if u > 0: pulls.append((bf - r_inj)/u)
                    f.Close()
                n = len(pulls)
                if n == 0: out.append('%-3s %-10s %-4s %.4f | 0 | %d' % (m, tr, inj, r_inj, len(files))); continue
                mu = sum(pulls)/n; sd = math.sqrt(sum((p-mu)**2 for p in pulls)/max(n-1, 1))
                h = ROOT.TH1F('h%s%s%s%s' % (tag, m, tr, inj), '', 80, -4, 4)
                for p in pulls: h.Fill(p)
                h.Fit('gaus', 'Q0'); g = h.GetFunction('gaus')
                gm, ge = (g.GetParameter(1), g.GetParError(1)) if g else (float('nan'), float('nan'))
                out.append('%-3s %-10s %-4s %.4f | %5d  %+.3f +- %.3f   %+.3f +- %.3f   %+.3f | %d' % (
                    m, tr, inj, r_inj, n, mu, sd/math.sqrt(n), gm, ge, sorted(pulls)[n//2], len(files)))
    txt = '\n'.join(out) + '\n'; open(os.path.join(H, 'summary_%s.txt' % tag), 'w').write(txt); print(tag); print(txt)
