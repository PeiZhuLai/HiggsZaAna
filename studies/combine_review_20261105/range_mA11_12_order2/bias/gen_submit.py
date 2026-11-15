#!/usr/bin/env python3
"""Condor submit file for the bias screen of the order-2 envelopes: every truth of cards/<tag>/<m>_truths.json x
injection q50 only (the bias requirement is judged at q50) x 10 chunks of 1000 toys (10000 toys per truth).
Usage : gen_submit.py <tag> <m> [<m> ...]
Input : cards/<tag>/<m>_truths.json
Output: bias_screen_<tag>.submit"""
import json, os, sys
B = os.path.dirname(os.path.abspath(__file__)); tag, masses = sys.argv[1], [int(x) for x in sys.argv[2:]]
rows = []
for m in masses:
    truths = json.load(open(os.path.join(B, 'cards', tag, '%d_truths.json' % m)))
    rows += ['%s %d %s q50 %d' % (tag, m, t, c) for t in truths for c in range(10)]
sub = '''universe       = vanilla
executable     = {B}/runjob.sh
request_memory = 3000
request_cpus   = 1
+JobFlavour    = "workday"
batch_name     = bias_screen_{tag}
on_exit_hold   = (ExitBySignal == True) || (ExitCode != 0)
periodic_release = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
output         = /dev/null
error          = {B}/logs/$(tag)_mA$(m)_$(tr)_$(inj)_c$(c).$(ClusterId).err
log            = {B}/logs/bias_screen.$(ClusterId).log
arguments      = $(tag) $(m) $(tr) $(inj) $(c) 1000
queue tag, m, tr, inj, c from (
{rows}
)
'''.format(B=B, tag=tag, rows='\n'.join(rows))
open(os.path.join(B, 'bias_screen_%s.submit' % tag), 'w').write(sub); print(len(rows), 'jobs')
