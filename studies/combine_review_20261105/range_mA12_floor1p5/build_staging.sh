#!/usr/bin/env bash
# Side study (Combine review, option a for mA12): the 1 GeV floor of the step width and smearing improves the Pow1
# truth bias of mA12 (-0.240 -> -0.216) but does not pass (range_mA11_12_28/bias/summary_flr.txt), and every member
# of the envelope puts the width and smearing on the 1 GeV floor. Here the floor of mA12 is 1.5 GeV (scan 1.5, 2, 3 GeV around the
# ~2 GeV m_llgg resolution). Staging build only; production Background/ is not touched.
#  - library: copy of production Background/{src,interface,makefile} (mA11 1 GeV floor already deployed); mA12 floor 1.5 GeV.
#  - fTest: source of the production binary (deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp, unchanged); rpath -> staging library.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
# Output: Background/{src,interface,lib,obj,bin} under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC=$D/../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
S=$D/Background
[ -f $D/../deploy_mA11/logs/stage1.DONE ] || { echo "production build of mA11 not done"; exit 1; }
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
python3 - "$S/src/PdfModelBuilder.cc" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
old = '  const double wsFloor = (massALP == 0 || massALP == 11) ? 1.0 : 0.0;'
new = '  const double wsFloor = (massALP == 0 || massALP == 11) ? 1.0 : (massALP == 12 ? 1.5 : 0.0);'
assert s.count(old) == 1, 'wsFloor line not found exactly once'
open(p, 'w').write(s.replace(old, new))
PY
grep -n 'wsFloor = ' $S/src/PdfModelBuilder.cc
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $SRC -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
ls -l $S/lib $S/bin
