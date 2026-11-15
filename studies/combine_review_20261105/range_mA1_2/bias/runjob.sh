#!/bin/bash
# One chunk of the bias screen on a staging envelope: generate toys from one truth (cards/<tag>/<m>_truths.json),
# fit them with the envelope at the stored ranges and the minimizer options of bias_fit_config.sh (as the "narrow"
# fit of bias_wide_truth/runjob.sh). At the opened masses the stored ranges are the opened ones.
# Usage : runjob.sh TAG MASS TRUTH INJ CHUNK NTOYS     INJ = 0 or q50
# Input : cards/<tag>/<m>_Datacard_leptons.root, cards/<tag>/<m>_truths.json
# Output: out/<tag>/mA<m>_<truth>_<inj>/chunk<c>.root (MultiDimFit tree), chunk<c>.done
set -eo pipefail
tag="$1"; m="$2"; tr="$3"; inj="$4"; c="$5"; n="$6"
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_mA1_2/bias
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
card="${H}/cards/${tag}/${m}_Datacard_leptons.root"
out="${H}/out/${tag}/mA${m}_${tr}_${inj}"; mkdir -p "$out"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
ulimit -s unlimited
wd="${_CONDOR_SCRATCH_DIR:-${H}/work}/${tag}_mA${m}_${tr}_${inj}_c${c}_$$"; mkdir -p "$wd"; cd "$wd"
read -r idx setp q50 < <(python3 -c "
import json; t=json.load(open('${H}/cards/${tag}/${m}_truths.json'))['${tr}']; print(t['index'], t['set'], t['q50'])")
[[ "$inj" == "q50" ]] && r_inj="$q50" || r_inj=0
seed=$(( 800000 + m*1000 + c )); [[ "$inj" == "q50" ]] && seed=$(( seed + 500 ))
echo "tag=${tag} mA=${m} truth=${tr} (index ${idx}) inj=${inj} r=${r_inj} chunk=${c} ntoys=${n} seed=${seed}"
combine -M GenerateOnly -d "$card" -m 125.38 --expectSignal "$r_inj" --saveToys -t "$n" -s "$seed" -n _gen \
  --setParameters "MH=125.38,pdfindex_cat0_13p6TeV=${idx},${setp}" --freezeParameters MH,pdfindex_cat0_13p6TeV > gen.log 2>&1
toys=$(ls higgsCombine_gen.GenerateOnly.mH125.38.${seed}.root)
source "${FF}/Combine/Checks/Bias_nominal/bias_fit_config.sh"
wide_opts="$(default_bias_fit_combine_options)"
opts="${wide_opts%%--setParameterRanges*}--setParameterRanges r=-5,5"
eval combine -M MultiDimFit -d "$card" -m 125.38 -P r --algo singles $opts -t "$n" -n _fit --toysFile="$toys" > fit.log 2>&1
cp -p higgsCombine_fit.MultiDimFit.mH125.38.*.root "${out}/chunk${c}.root"
date > "${out}/chunk${c}.done"
echo "[OK] ${tag} mA=${m} truth=${tr} inj=${inj} chunk=${c}"
