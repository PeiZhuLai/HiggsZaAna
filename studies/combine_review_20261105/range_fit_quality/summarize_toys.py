#!/usr/bin/env python3
"""Tail probabilities from toys_dS.py. Input: toys/*.json. Output: toys_summary.txt"""
import glob, json
out = ['mA pdf ntoys | dS: obs, toy median, q84, q95, N(toy>=obs) | S_prod: obs, N(toy>=obs) | toys with a fit status != 0']
for f in sorted(glob.glob('toys/*.json'), key=lambda x: int(x.split('mA')[1].split('_')[0])):
    d = json.load(open(f)); t = d['toys']; n = len(t); s = sorted(x['dS'] for x in t)
    out.append('%-3d %-5s %d | %.2f %.2f %.2f %.2f %d | %.2f %d | %d' % (d['mA'], d['pdf'], n, d['obs_dS'], s[n//2], s[int(.84*n)], s[int(.95*n)], sum(x >= d['obs_dS'] for x in s), d['obs_S_prod'], sum(x['S_prod'] >= d['obs_S_prod'] for x in t), sum(any(x['st']) for x in t)))
open('toys_summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
