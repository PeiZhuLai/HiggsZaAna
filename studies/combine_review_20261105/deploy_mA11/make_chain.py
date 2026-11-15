#!/usr/bin/env python3
"""Writes chain_deploy.sh for the mA11 deployment (1 GeV floor of the step width and smearing) from
deploy_mA1_2/chain_deploy.sh; every edit must match the expected number of times.
Input : ../deploy_mA1_2/chain_deploy.sh
Output: chain_deploy.sh
Usage : make_chain.py"""
import os
D = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(D, '..', 'deploy_mA1_2', 'chain_deploy.sh')).read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, (old, s.count(old))
    s = s.replace(old, new)


rep('''# Production deployment of the mA1 and mA2 range decision of the Combine review: the step/smearing ranges are opened
# with the turn-on capped below the signal region (turn-on 95-113, width 0.1-50, smearing 0.05-20 GeV) and the data
# are fitted from several starting points.
# Reason: with the per-mass windows the bias fit (production ranges) fails the bias requirement (mA1 Pow1 +0.212,
# mA2 Exp1 +0.218, Pow3 +0.234); the fully opened ranges put the turn-on inside the signal region. The capped ranges
# pass the bias requirement and the sideband goodness of fit (range_mA1_2/bias/summary_cap113.txt,
# range_mA1_2/gof/gof_summary.txt). mA1 is fitted from 100 GeV, as in production.''',
    '''# Production deployment of the mA11 decision of the Combine review: the lower edge of the step width and of the
# Gaussian smearing is raised from 0.1 / 0.05 GeV to 1 GeV, as already done for massALP = 0.
# Reason: with the opened ranges the Pow1 truth fails the bias requirement (-0.223); the data fit puts the smearing on
# its 0.05 GeV floor, narrower than the FFT binning and far below the ~2 GeV m_llgg resolution. With the 1 GeV floor
# every truth passes (Pow1 -0.193, range_mA11_12_28/bias/summary_flr.txt) and the sideband goodness of fit passes
# (p 0.67-0.78, range_mA11_12_28/gof/gof_summary.txt).''')
rep('''#  1  build: PdfModelBuilder.cc openStepGausRanges() gets mA1, 2 and their turn-on upper edge is 113 GeV; the fTest
#     multiStartMass() gets mA1, 2; library and fTest rebuilt (fTest source = deploy_20261113/src/
#     fTest_ALP_turnOn.deploy.cpp, the source of the production binary). Regression: mA19, mA10 and mA9 (not
#     redeployed) must reproduce production exactly.
#  2  staging fTest of mA1, 2 -> work/fit_results/<m> (production untouched); must be identical to the side study
#     range_mA1_2/fit_results/<m>.''',
    '''#  1  build: PdfModelBuilder.cc wsFloor gets mA11; library and fTest rebuilt (fTest source unchanged =
#     deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp, the source of the production binary). Regression: mA19, mA12 and
#     mA1 (not redeployed) must reproduce production exactly.
#  2  staging fTest of mA11 -> work/fit_results/11 (production untouched); must be identical to the side study
#     range_mA11_12_28/fit_results/11.''')
rep('check_params for mA1, 2, both cards', 'check_params for mA11, both cards')
rep('D=$ST/deploy_mA1_2;', 'D=$ST/deploy_mA11;')
rep('DIED deploy_mA1_2', 'DIED deploy_mA11')
rep('DONE deploy_mA1_2 stages', 'DONE deploy_mA11 stages')
rep('bak_deploy_mA1_2', 'bak_deploy_mA11', 4)
rep('SRC0=$ST/deploy_20261113/src/fTest_ALP_turnOn.deploy.cpp\nSIDE=$ST/range_mA1_2/fit_results',
    'SRC0=$ST/deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp\nSIDE=$ST/range_mA11_12_28/fit_results')
rep('MASSES="1 2"', 'MASSES="11"')
rep('''old = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 9 ||'
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
    open(p, 'w').write(s); print('edited')''',
    '''old = '  const double wsFloor = (massALP == 0) ? 1.0 : 0.0;'
new = ('  // mA11 (Combine review): with the 0.05 GeV smearing floor the Pow1 truth fails the bias requirement; the 1 GeV\\n'
       '  // floor is narrower than the m_llgg resolution and every truth passes.\\n'
       '  const double wsFloor = (massALP == 0 || massALP == 11) ? 1.0 : 0.0;')
if new in s: print('already edited')
else:
    assert s.count(old) == 1, 'wsFloor line not found'
    s = s.replace(old, new)
    open(p, 'w').write(s); print('edited')''')
rep('''old = 'static bool multiStartMass(int m) { return m == 9 ||'
new = 'static bool multiStartMass(int m) { return m == 1 || m == 2 || m == 9 ||'
assert s.count(old) == 1, 'multiStartMass line not found'
open(sys.argv[2], 'w').write(s.replace(old, new))''',
    '''assert 'static bool multiStartMass(int m) { return m == 1 || m == 2 || m == 9 ||' in s, 'unexpected fTest source'
open(sys.argv[2], 'w').write(s)''')
rep('''hb "    library and fTest rebuilt; regression mA19, mA10 and mA9 against production"
rm -rf $D/regress/19 $D/regress/10 $D/regress/9
printf '19 %s\\n10 %s\\n9 %s\\n' $D/regress/19 $D/regress/10 $D/regress/9 | xargs -P3 -L1 bash -c 'ftest_one $0 $1'
for m in 19 10 9; do grep -q 'rc=0' $D/regress/$m/.ftest_rc || die "regression fTest mA$m failed"; done
for m in 19 10 9; do''', '''hb "    library and fTest rebuilt; regression mA19, mA12 and mA1 against production"
rm -rf $D/regress/19 $D/regress/12 $D/regress/1
printf '19 %s\\n12 %s\\n1 %s\\n' $D/regress/19 $D/regress/12 $D/regress/1 | xargs -P3 -L1 bash -c 'ftest_one $0 $1'
for m in 19 12 1; do grep -q 'rc=0' $D/regress/$m/.ftest_rc || die "regression fTest mA$m failed"; done
for m in 19 12 1; do''')
rep('hb "    regression: mA19, mA10 and mA9 IDENTICAL to production"', 'hb "    regression: mA19, mA12 and mA1 IDENTICAL to production"')
rep('hb "    4 workspaces, 4 blind expected limits, check_params ALL OK"',
    'hb "    2 workspaces, 2 blind expected limits, check_params ALL OK"')
assert 'mA1, 2' not in s and 'deploy_mA1_2/logs' not in s, 'leftover mA1, 2 text'
open(os.path.join(D, 'chain_deploy.sh'), 'w').write(s)
print('written')
