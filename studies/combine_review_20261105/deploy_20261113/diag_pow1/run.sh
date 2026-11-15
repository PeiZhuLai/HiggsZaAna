#!/usr/bin/env bash
# Diagnostic for the Pow1 truth bias after the wide-range deployment (mA11 -0.223, mA12 -0.240, mA9 -0.209; mA17
# -0.192 as a reference). For each mass: 300 toys generated from Pow1 (production card, best-fit snapshot, injection
# as in the production highstat run), then fitted twice with the production bias-fit options:
#   env  : discrete profiling over the envelope, the selected index saved per toy
#   pow1 : index frozen to Pow1 (truth function only)
# Unbiased pow1 with biased env means the envelope picks other functions on these toys.
# Input : Combine/root_t2w/<m>_Datacard_leptons.root (cross-section card, read only)
# Output: out/mA<m>/{toys,env,pow1}.root and logs; logs/run.DONE
# Usage : run.sh   (4 processes in the background on lxplus)
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd)
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
mkdir -p $D/logs; rm -f $D/logs/run.DONE
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
source $FF/Combine/Checks/Bias_nominal/bias_fit_config.sh
OPTS="$(resolve_bias_fit_combine_options)"
IDX=pdfindex_cat0_13p6TeV
one() {
  local m=$1 ip=$2 r=$3 o=$D/out/mA$1 card=$FF/Combine/root_t2w/$1_Datacard_leptons.root
  rm -rf $o; mkdir -p $o; cd $o
  combine -M GenerateOnly -d $card -m 125.38 --expectSignal $r --saveToys -t 300 -s $((910000 + m)) -n _gen \
    --setParameters $IDX=$ip --freezeParameters $IDX,MH > gen.log 2>&1 || { echo "rc=gen_fail" > .rc; return; }
  mv higgsCombine_gen.GenerateOnly.mH125.38.*.root toys.root
  eval combine -M MultiDimFit -d $card -m 125.38 -P r --algo singles $OPTS -t 300 --toysFile toys.root \
    --saveSpecifiedIndex $IDX -n _env > env.log 2>&1
  mv higgsCombine_env.MultiDimFit.mH125.38.*.root env.root
  eval combine -M MultiDimFit -d $card -m 125.38 -P r --algo singles ${OPTS/--freezeParameters MH/--freezeParameters MH,$IDX} \
    --setParameters $IDX=$ip -t 300 --toysFile toys.root -n _pow1 > pow1.log 2>&1
  mv higgsCombine_pow1.MultiDimFit.mH125.38.*.root pow1.root
  echo "rc=0" > .rc
}
one 9 0 0.0168 & one 11 2 0.0223 & one 12 3 0.0238 & one 17 1 0.0285 & wait
for m in 9 11 12 17; do echo "mA$m $(cat $D/out/mA$m/.rc 2>/dev/null || echo rc=missing)"; done > $D/logs/run.DONE
