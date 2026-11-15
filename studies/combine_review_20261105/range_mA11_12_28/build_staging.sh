#!/usr/bin/env bash
# Side studies (option B of the Combine review discussion), staging build only; production Background/ is not touched.
#  mA11, 12: after the wide-range deployment the Pow1 truth still fails the bias requirement (-0.223, -0.240); the
#    Pow1 and Lau1 data fits put the smearing on its 0.05 GeV floor with a 3.6-5.1 GeV step width, while the
#    envelope picks Lau1 or Exp1 on 20-30% of the Pow1 toys (deploy_20261113/diag_pow1). Here width and smearing get
#    the 1 GeV floor already used for massALP = 0 (a 0.05 GeV smearing is narrower than the FFT binning and far below
#    the ~2 GeV m_llgg resolution).
#  mA28: the sideband GOF of the envelope (Exp1, Pow1, Lau1) fails; Bernstein is re-enabled up to order 6 as in
#    gof_fix_mA15_28 (F-test picks Bern5).
#  - library: copy of production Background/{src,interface,makefile}; the width/smearing floor is 1 GeV for mA0, 11, 12.
#  - fTest: source of the production binary (deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp) with case 28 -> order 6
#    in the Bernstein switch; rpath -> the staging library.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
# Output: Background/{src,interface,lib,obj,bin} and fTest_ALP_turnOn.build.cpp under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC=$D/../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
S=$D/Background
[ -f $D/../deploy_mA1_2/logs/chain.DONE ] || { echo "production deployment of mA1, 2 not done"; exit 1; }
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
python3 - "$S/src/PdfModelBuilder.cc" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
old = '  const double wsFloor = (massALP == 0) ? 1.0 : 0.0;'
new = '  const double wsFloor = (massALP == 0 || massALP == 11 || massALP == 12) ? 1.0 : 0.0;'
assert s.count(old) == 1, 'wsFloor line not found exactly once'
open(p, 'w').write(s.replace(old, new))
EOF
grep -n 'wsFloor' $S/src/PdfModelBuilder.cc
python3 - $SRC $D/fTest_ALP_turnOn.build.cpp <<'EOF'
import sys
s = open(sys.argv[1]).read()
old = '    case 29: return 6; case 30: return 6;'
new = '    case 28: return 6; case 29: return 6; case 30: return 6;'
assert s.count(old) == 1, 'Bernstein switch line not found exactly once'
open(sys.argv[2], 'w').write(s.replace(old, new))
EOF
grep -n 'case 28' $D/fTest_ALP_turnOn.build.cpp
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $D/fTest_ALP_turnOn.build.cpp -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
ls -l $S/lib $S/bin
