#!/usr/bin/env python3
"""Condor submit file for the bias screen of one staging envelope: every truth of cards/<tag>/<m>_truths.json x
r = 0, q50 x 10 chunks of 1000 toys (10000 toys per truth and injection).
Usage : [EXE=runjob_v2.sh] gen_submit.py <tag> <m>    (executable runjob.sh unless EXE is set)
Input : cards/<tag>/<m>_truths.json
Output: bias_screen_<tag>_mA<m>.submit"""
import json, os, sys
B = os.path.dirname(os.path.abspath(__file__)); tag, m = sys.argv[1], int(sys.argv[2])
truths = json.load(open(os.path.join(B, 'cards', tag, '%d_truths.json' % m)))
rows = ['%s %d %s %s %d' % (tag, m, t, inj, c) for t in truths for inj in ('0', 'q50') for c in range(10)]
sub = '''universe       = vanilla
executable     = {B}/{exe}
request_memory = 3000
request_cpus   = 1
+JobFlavour    = "workday"
batch_name     = bias_screen_{tag}_mA{m}
on_exit_hold   = (ExitBySignal == True) || (ExitCode != 0)
periodic_release = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
output         = /dev/null
error          = {B}/logs/$(tag)_mA$(m)_$(tr)_$(inj)_c$(c).$(ClusterId).err
log            = {B}/logs/bias_screen.$(ClusterId).log
arguments      = $(tag) $(m) $(tr) $(inj) $(c) 1000
queue tag, m, tr, inj, c from (
{rows}
)
'''.format(exe=os.environ.get('EXE', 'runjob.sh'), B=B, tag=tag, m=m, rows='\n'.join(rows))
open(os.path.join(B, 'bias_screen_%s_mA%d.submit' % (tag, m)), 'w').write(sub); print(len(rows), 'jobs')
