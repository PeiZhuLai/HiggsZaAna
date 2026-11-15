#!/usr/bin/env bash
# End-to-end test of the renaming (Wilson card, production otherwise unchanged):
#  A  production workspace (old names)                              -> blind limits
#  B  apply_naming.py --no-split + renamed signal workspaces        -> must equal A (pure renaming)
#  C  apply_naming.py (trigger split) + renamed signal workspaces   -> effect of the electron/muon split
#  plus check_params.py on B: every param nuisance of the card is a variable of the workspace that changes the
#  signal model when moved to +1 (a name mismatch would leave it free-standing and silently ineffective).
# Input : Datacard/output_Datacard_leptons/<m>_pruned_datacard_leptons.txt, Combine/root_t2w_wilson/<m>_Datacard_leptons.root (read only)
# Output: e2e/<m>_{B,C}.{txt,root}, e2e/logs/, e2e/limits.txt, e2e/params_<m>.txt, e2e/.done
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); E=$S/e2e; C=~/doc/HZa/review_comments/Combine_review
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit; DC=$FF/Datacard/output_Datacard_leptons
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh); mkdir -p $E/logs; rm -f $E/.done $E/limits.txt
LIM="-M AsymptoticLimits -m 125.38 --run blind --rMax 2 --cminDefaultMinimizerStrategy 0 --rAbsAcc 1e-6 --rRelAcc 1e-6 --freezeParameters MH"
for m in ${MASSES:-1 10 11}; do
  W=$DC/${m}_pruned_datacard_leptons.txt
  python3 $S/rename_ws.py $S/ws_named $(awk '$1=="shapes" && $2 ~ /^ggH_/ {print $4}' $W) > $E/logs/${m}_rename.log 2>&1 || { echo "mA$m rename failed"; continue; }
  python3 $C/apply_naming.py $W $E/${m}_B.txt --signal-ws-dir $S/ws_named --no-split
  python3 $C/apply_naming.py $W $E/${m}_C.txt --signal-ws-dir $S/ws_named
  for t in B C; do text2workspace.py $E/${m}_$t.txt -o $E/${m}_$t.root -m 125.38 higgsMassRange=115,135 > $E/logs/${m}_${t}_t2w.log 2>&1; done
  python3 $S/check_params.py $E/${m}_B.txt $E/${m}_B.root > $E/params_${m}.txt 2>&1
  cd $E/logs
  combine $LIM -d $FF/Combine/root_t2w_wilson/${m}_Datacard_leptons.root -n _A_mA$m > ${m}_A.log 2>&1
  for t in B C; do combine $LIM -d $E/${m}_$t.root -n _${t}_mA$m > ${m}_$t.log 2>&1; done
  for t in A B C; do echo "mA$m $t $(grep -E 'Expected (2.5|16.0|50.0|84.0|97.5)%' ${m}_$t.log | awk '{print $NF}' | tr '\n' ' ')"; done >> $E/limits.txt
  cd $FF/..
done
date > $E/.done
