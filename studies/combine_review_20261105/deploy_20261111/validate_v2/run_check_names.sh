#!/usr/bin/env bash
# Official check_names.py (cms-analysis/general/systematics, copy in ../../naming_check/systematics) on the deployed
# cards of representative mass points: mA1 (dropped unwired param, muon trigger present), mA10 (boundary bias point),
# mA15 and mA22 (envelope changed by the opened ranges), Wilson and cross-section cards.
# Uses the card copies with relative shapes paths made by run_all.sh (cards/<m>_<card>.txt) and systematics_draft.yml.
# Input : cards/{1,10,15,22}_{wilson,xsec}.txt, ~/doc/HZa/review_comments/Combine_review/systematics_draft.yml
# Output: check_names/<m>_<card>.log, check_names_summary.txt
set -o pipefail
V=$(cd "$(dirname "$0")" && pwd); S=$V/../../naming_check; C=$HOME/doc/HZa/review_comments/Combine_review
IMG=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/combine-container:CMSSW_14_1_0_pre4-combine_v10.6.1-harvester_v3.1.0
mkdir -p $V/check_names
apptainer exec -B /afs -B /eos -B /cvmfs $IMG bash -c "source /cvmfs/cms.cern.ch/cmsset_default.sh 2>/dev/null; cd /home/cmsusr/CMSSW_14_1_0_pre4/src && eval \$(scramv1 runtime -sh) 2>/dev/null; cd $S
for m in 1 10 15 22 29; do for c in wilson xsec; do
  python3 systematics/check_names.py --input $V/cards/\${m}_\${c}.txt --systematics-dict $C/systematics_draft.yml --global-dict systematics/systematics_master.yml > $V/check_names/\${m}_\${c}.log 2>&1
  echo \"\${m}_\${c} rc=\$? \$(grep 'Summary:' $V/check_names/\${m}_\${c}.log)\"
done; done" > $V/check_names_summary.txt 2>&1
cat $V/check_names_summary.txt
