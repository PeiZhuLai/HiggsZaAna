#!/usr/bin/env bash
# Side study (Combine review, mA11): is the 1 GeV lower limit of the step width and smearing still needed once the
# multi-start also scans the exponent? Staging build with the production sources, except that the mA11 floor is
# removed from PdfModelBuilder.cc (wsFloor only for massALP = 0, i.e. mA11 back to 0.1 / 0.05 GeV like the other
# opened masses). fTest = production source (deploy_mA11_12_ps/src/fTest_ALP_turnOn.deploy.cpp, exponent scan at
# mA11 and mA12 written in the code), linked to the staging library. Production Background/ is not touched.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_mA11_12_ps/src/fTest_ALP_turnOn.deploy.cpp
# Output: Background/{src,interface,lib,obj,bin} under this directory; logs/build.log, logs/PdfModelBuilder.diff
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC=$D/../deploy_mA11_12_ps/src/fTest_ALP_turnOn.deploy.cpp; S=$D/Background
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
python3 - $S/src/PdfModelBuilder.cc <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
old = '  const double wsFloor = (massALP == 0 || massALP == 11) ? 1.0 : 0.0;'
new = '  const double wsFloor = (massALP == 0) ? 1.0 : 0.0;  // side study: mA11 floor removed'
assert s.count(old) == 1, 'wsFloor line not found'
open(p, 'w').write(s.replace(old, new)); print('patched')
PY
diff $F/Background/src/PdfModelBuilder.cc $S/src/PdfModelBuilder.cc > $D/logs/PdfModelBuilder.diff || true
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt); LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $SRC -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
date > $D/logs/build.DONE
