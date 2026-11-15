#!/usr/bin/env python3
"""Writes chain_deploy.sh for the mA1, 2 deployment from deploy_20261113/chain_deploy.sh; every edit must match once.
Input : ../deploy_20261113/chain_deploy.sh
Output: chain_deploy.sh
Usage : make_chain.py"""
import os
D = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(D, '..', 'deploy_20261113', 'chain_deploy.sh')).read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, (old, s.count(old))
    s = s.replace(old, new)


rep('''# Production deployment of the per-mass range decision of the Combine review: the step/smearing ranges are opened
# (turn-on 95-125, width 0.1-50, smearing 0.05-20 GeV) and the data are fitted from several starting points at
# mA9, 11, 12, 13, 16, 17, 20, as already done at mA10, 15, 21-29.
# Reason: with the bias fit profiling the production ranges, the narrow ranges fail the bias requirement at these
# points (range_per_mass/logs/bias_compare_consistent.txt) and, except mA9, the wide ranges pass the wide-truth test
# (range_per_mass/summary.txt); turn-ons of the wide fits are 97-103 GeV, outside the signal region.''',
    '''# Production deployment of the mA1 and mA2 range decision of the Combine review: the step/smearing ranges are opened
# with the turn-on capped below the signal region (turn-on 95-113, width 0.1-50, smearing 0.05-20 GeV) and the data
# are fitted from several starting points.
# Reason: with the per-mass windows the bias fit (production ranges) fails the bias requirement (mA1 Pow1 +0.212,
# mA2 Exp1 +0.218, Pow3 +0.234); the fully opened ranges put the turn-on inside the signal region. The capped ranges
# pass the bias requirement and the sideband goodness of fit (range_mA1_2/bias/summary_cap113.txt,
# range_mA1_2/gof/gof_summary.txt). mA1 is fitted from 100 GeV, as in production.''')
rep('''#  1  build: PdfModelBuilder.cc openStepGausRanges() and the fTest multiStartMass() get the 7 masses; library and
#     fTest rebuilt (fTest source = deploy_20261111/src/fTest_ALP_turnOn.deploy.cpp, the source of the production
#     binary). Regression: mA19 and mA10 (not redeployed) must reproduce production exactly.
#  2  staging fTest of the 7 masses -> work/fit_results/<m> (production untouched); report envelope and parameters.''',
    '''#  1  build: PdfModelBuilder.cc openStepGausRanges() gets mA1, 2 and their turn-on upper edge is 113 GeV; the fTest
#     multiStartMass() gets mA1, 2; library and fTest rebuilt (fTest source = deploy_20261113/src/
#     fTest_ALP_turnOn.deploy.cpp, the source of the production binary). Regression: mA19, mA10 and mA9 (not
#     redeployed) must reproduce production exactly.
#  2  staging fTest of mA1, 2 -> work/fit_results/<m> (production untouched); must be identical to the side study
#     range_mA1_2/fit_results/<m>.''')
rep('#  4  text2workspace + blind expected limits + check_params for the 7 masses, both cards; old files backed up.',
    '#  4  text2workspace + blind expected limits + check_params for mA1, 2, both cards; old files backed up.')
rep('D=$ST/deploy_20261113;', 'D=$ST/deploy_mA1_2;')
rep('DIED deploy_20261113', 'DIED deploy_mA1_2')
rep('DONE deploy_20261113 stages', 'DONE deploy_mA1_2 stages')
rep('bak_deploy_20261113', 'bak_deploy_mA1_2', 4)
rep('SRC0=$ST/deploy_20261111/src/fTest_ALP_turnOn.deploy.cpp',
    'SRC0=$ST/deploy_20261113/src/fTest_ALP_turnOn.deploy.cpp\nSIDE=$ST/range_mA1_2/fit_results')
rep('MASSES="9 11 12 13 16 17 20"', 'MASSES="1 2"')
rep('''ftest_one() {  # <m> <out dir>; production command line of shellScripts/bkg/fit_bkg.sh (mhLow 95 for m > 1, 100 GOF toys)
  local m=$1 P=$2''', '''ftest_one() {  # <m> <out dir>; production command line of shellScripts/bkg/fit_bkg.sh (mhLow 100 for mA1, 95 otherwise, 100 GOF toys)
  local m=$1 P=$2 lo=95
  [ "$m" = 1 ] && lo=100''')
rep('--mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 )',
    '--mhLow $lo --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 )')
rep('''old = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 10 || massALP == 15 || (massALP >= 21 && massALP <= 29); }'
new = ('inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 9 || massALP == 10 || massALP == 11 || massALP == 12 || massALP == 13 || massALP == 15 || massALP == 16 || massALP == 17 || massALP == 20 || (massALP >= 21 && massALP <= 29); }')
if new in s: print('already edited')
else:
    assert s.count(old) == 1, 'openStepGausRanges line not found'
    s = s.replace(old, new)
    s = s.replace('// mA10, 15, 21-23 opened as well', '// mA9, 11, 12, 13, 16, 17, 20 opened as well (Combine review: with the bias fit profiling the production ranges\\n// the per-mass windows fail the bias requirement there; the wide ranges pass the wide-truth test).\\n// mA10, 15, 21-23 opened as well', 1)
    open(p, 'w').write(s); print('edited')''',
    '''old = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 9 ||'
new = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 1 || massALP == 2 || massALP == 9 ||'
old2 = '  openRange(turnon, turnonLo, turnonHi, 95.0, 125.0);'
new2 = '  openRange(turnon, turnonLo, turnonHi, 95.0, (massALP == 1 || massALP == 2) ? 113.0 : 125.0);'
c = '// mA9, 11, 12, 13, 16, 17, 20 opened as well'
if new in s and new2 in s: print('already edited')
else:
    assert s.count(old) == 1, 'openStepGausRanges line not found'
    assert s.count(old2) == 1, 'turn-on openRange line not found'
    assert s.count(c) == 1, 'comment anchor not found'
    s = s.replace(old, new).replace(old2, new2)
    s = s.replace(c, '// mA1, 2 opened as well, with the turn-on capped at 113 GeV below the signal region (Combine review: the\\n// per-mass windows fail the bias requirement there; the fully opened turn-on moves into the signal region).\\n' + c, 1)
    open(p, 'w').write(s); print('edited')''')
rep('''old = 'static bool multiStartMass(int m) { return m == 10 || m == 15 || (m >= 21 && m <= 29); }'
new = 'static bool multiStartMass(int m) { return m == 9 || m == 10 || m == 11 || m == 12 || m == 13 || m == 15 || m == 16 || m == 17 || m == 20 || (m >= 21 && m <= 29); }'
''', '''old = 'static bool multiStartMass(int m) { return m == 9 ||'
new = 'static bool multiStartMass(int m) { return m == 1 || m == 2 || m == 9 ||'
''')
rep('''hb "    library and fTest rebuilt; regression mA19 and mA10 against production"
rm -rf $D/regress/19 $D/regress/10
printf '19 %s\\n10 %s\\n' $D/regress/19 $D/regress/10 | xargs -P2 -L1 bash -c 'ftest_one $0 $1'
grep -q 'rc=0' $D/regress/19/.ftest_rc && grep -q 'rc=0' $D/regress/10/.ftest_rc || die "regression fTest failed"
for m in 19 10; do''', '''hb "    library and fTest rebuilt; regression mA19, mA10 and mA9 against production"
rm -rf $D/regress/19 $D/regress/10 $D/regress/9
printf '19 %s\\n10 %s\\n9 %s\\n' $D/regress/19 $D/regress/10 $D/regress/9 | xargs -P3 -L1 bash -c 'ftest_one $0 $1'
for m in 19 10 9; do grep -q 'rc=0' $D/regress/$m/.ftest_rc || die "regression fTest mA$m failed"; done
for m in 19 10 9; do''')
rep('hb "    regression: mA19 and mA10 IDENTICAL to production"', 'hb "    regression: mA19, mA10 and mA9 IDENTICAL to production"')
rep('''  hb "    mA$m $(envelope $D/work/fit_results/$m/HZAmassInde_fTest/EnvelopeResults.txt)"''',
    '''  hb "    mA$m $(envelope $D/work/fit_results/$m/HZAmassInde_fTest/EnvelopeResults.txt)"
  python3 $STG/compare_ws.py $SIDE/$m/$M $D/work/fit_results/$m/$M > $L/compare_side_$m.txt 2>&1
  tail -1 $L/compare_side_$m.txt | grep -qx IDENTICAL || die "staging mA$m differs from the side study ($L/compare_side_$m.txt)"''')
rep('''  local m=$1 P=$FR/$1
''', '''  local m=$1 P=$FR/$1 lo=95
  [ "$m" = 1 ] && lo=100
''')
rep('--mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi 172.13',
    '--mhLow $lo --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi 172.13')
rep('hb "    14 workspaces, 14 blind expected limits, check_params ALL OK"',
    'hb "    4 workspaces, 4 blind expected limits, check_params ALL OK"')
open(os.path.join(D, 'chain_deploy.sh'), 'w').write(s)
print('written')
