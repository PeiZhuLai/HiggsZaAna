#!/usr/bin/env python3
"""Pulls of the wide-truth bias test, same definition as RunBiasStudy.py: (r_hat - r_inj) / sigma, with sigma the
uncertainty on the side facing the injected value; toys without a valid interval are dropped. Mean and its error
from the sample (and the Gaussian-fit mean as in the AN, on a TH1 with 80 bins in [-4, 4]).
Input : out/mA<m>_<inj>/chunk*_{narrow,wide}.root, truth.json. Output: summary.txt"""
import glob, json, math, os, ROOT
ROOT.gROOT.SetBatch(True)
H = os.path.dirname(os.path.abspath(__file__)); truth = json.load(open(os.path.join(H, 'truth.json')))
out = ['mA  truth  inj    r_inj    | fit     ntoys  mean_pull +- err   gaus_mean +- err   median | chunks']
for m in sorted(truth, key=int):
    for inj in ('0', 'q50'):
        r_inj = truth[m]['q50'] if inj == 'q50' else 0.0
        for cfg in ('narrow', 'wide'):
            files = sorted(glob.glob(os.path.join(H, 'out', 'mA%s_%s' % (m, inj), 'chunk*_%s.root' % cfg)))
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
            if n == 0: out.append('%-3s %-6s %-4s %.4f | %-6s  0' % (m, truth[m]['pdf'], inj, r_inj, cfg)); continue
            mu = sum(pulls)/n; sd = math.sqrt(sum((p-mu)**2 for p in pulls)/max(n-1, 1))
            h = ROOT.TH1F('h%s%s%s' % (m, inj, cfg), '', 80, -4, 4)
            for p in pulls: h.Fill(p)
            h.Fit('gaus', 'Q0'); g = h.GetFunction('gaus')
            gm, ge = (g.GetParameter(1), g.GetParError(1)) if g else (float('nan'), float('nan'))
            med = sorted(pulls)[n//2]
            out.append('%-3s %-6s %-4s %.4f | %-6s  %5d  %+.3f +- %.3f   %+.3f +- %.3f   %+.3f | %d' % (
                m, truth[m]['pdf'], inj, r_inj, cfg, n, mu, sd/math.sqrt(n), gm, ge, med, len(files)))
open(os.path.join(H, 'summary.txt'), 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
