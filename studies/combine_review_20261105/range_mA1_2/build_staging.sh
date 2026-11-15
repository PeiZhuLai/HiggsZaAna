#!/usr/bin/env bash
# Side study for mA1 and mA2: with the per-mass windows the bias fit (production ranges) fails the bias requirement
# (mA1 Pow1 +0.212, mA2 Exp1 +0.218 and Pow3 +0.234), and the fully opened ranges put the turn-on inside the signal
# region (116.6, 119.5 GeV). Here the step/smearing ranges are opened with the turn-on capped below the signal region:
# turn-on 95-113 GeV, width 0.1-50 GeV, smearing 0.05-20 GeV, and the data are fitted from several starting points.
# Production Background/ is not touched.
#  - library: copy of production Background/{src,interface,makefile}; openStepGausRanges() gets mA1, 2 and the
#    turn-on upper edge is 113 GeV for them (125 GeV elsewhere, unchanged).
#  - fTest: source of the production binary (deploy_20261113/src/fTest_ALP_turnOn.deploy.cpp) with mA1, 2 added to
#    multiStartMass(); rpath -> the staging library.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_20261113/src/fTest_ALP_turnOn.deploy.cpp
# Output: Background/{src,interface,lib,obj,bin} under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC=$D/../deploy_20261113/src/fTest_ALP_turnOn.deploy.cpp
S=$D/Background
[ -f $D/../deploy_20261113/logs/stage1.DONE ] || { echo "production build of deploy_20261113 not done"; exit 1; }
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
python3 - "$S/src/PdfModelBuilder.cc" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
old = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 9 ||'
new = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 1 || massALP == 2 || massALP == 9 ||'
assert s.count(old) == 1, 'openStepGausRanges line not found exactly once'
s = s.replace(old, new)
old = '  openRange(turnon, turnonLo, turnonHi, 95.0, 125.0);'
new = '  openRange(turnon, turnonLo, turnonHi, 95.0, (massALP == 1 || massALP == 2) ? 113.0 : 125.0);'
assert s.count(old) == 1, 'turn-on openRange line not found exactly once'
open(p, 'w').write(s.replace(old, new))
EOF
grep -n 'inline bool openStepGausRanges\|openRange(turnon' $S/src/PdfModelBuilder.cc
python3 - $SRC $S/fTest_ALP_turnOn.build.cpp <<'EOF'
import sys
s = open(sys.argv[1]).read()
old = 'static bool multiStartMass(int m) { return m == 9 ||'
new = 'static bool multiStartMass(int m) { return m == 1 || m == 2 || m == 9 ||'
assert s.count(old) == 1, 'multiStartMass line not found exactly once'
open(sys.argv[2], 'w').write(s.replace(old, new))
EOF
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $S/fTest_ALP_turnOn.build.cpp -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
ls -l $S/lib $S/bin
