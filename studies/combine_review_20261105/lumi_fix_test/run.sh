#!/usr/bin/env bash
# Effect of the LUM POG lumi fix on the blind expected limits (cross-section card): replace the 8 lumi lnN lines
# (5 Uncorrelated with the yearly totals + 3 "Correlated" lines holding the combined-dataset totals 1.01/1.41/1.20%)
# by the 3 nuisances of LumiRecommendationsRun3 (2022+2023+2024 reduction), then text2workspace and AsymptoticLimits
# --run blind with the options of range_full_refit/refit_limit.py; reference: production root_t2w_xsec, same options.
# Input : Datacard/output_Datacard_leptons/<m>_pruned_datacard_leptons_xsec.txt, Combine/root_t2w_xsec/<m>_Datacard_leptons.root (read only)
# Output: cards/<m>_lumifix.{txt,root}, logs/<m>_{old,new}.log, limits.txt, .done
set -o pipefail
H=$(cd "$(dirname "$0")" && pwd); FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh); mkdir -p $H/cards $H/logs; rm -f $H/.done
LIM="-M AsymptoticLimits -m 125.38 --run blind --rMax 2 --cminDefaultMinimizerStrategy 0 --rAbsAcc 1e-5 --rRelAcc 1e-5 --freezeParameters MH"
for m in ${MASSES:-1 10 20 30}; do
  python3 - $FF/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt $H/cards/${m}_lumifix.txt <<'PY'
import re, sys
NEW = {'lumi_13p6TeV_222324': {'2022': '1.0138', '2023': '1.0017', '2024': '1.0020'},
       'lumi_13p6TeV_2324': {'2023': '1.0127', '2024': '1.0068'},
       'lumi_13p6TeV_2024': {'2024': '1.0144'}}
out, procs, done = [], None, False
for l in open(sys.argv[1]):
    f = l.split()
    if f and f[0] == 'process' and procs is None: procs = f[1:]
    if len(f) > 1 and f[1] == 'lnN' and f[0].startswith('lumi_13p6TeV_'):
        if not done:
            for n, v in NEW.items():
                cols = []
                for p in procs:
                    mm = re.match(r'ggH_(\d{4})', p)
                    cols.append(v.get(mm.group(1), '-') if mm else '-')
                out.append('%-51s %-13s ' % (n, 'lnN') + ' '.join('%-15s' % c for c in cols) + '\n')
            done = True
        continue
    out.append(l)
assert done
open(sys.argv[2], 'w').writelines(out)
PY
  text2workspace.py $H/cards/${m}_lumifix.txt -o $H/cards/${m}_lumifix.root -m 125.38 higgsMassRange=115,135 > $H/logs/${m}_t2w.log 2>&1
  cd $H/logs
  combine $LIM -d $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -n _old_mA$m > ${m}_old.log 2>&1
  combine $LIM -d $H/cards/${m}_lumifix.root -n _new_mA$m > ${m}_new.log 2>&1
  for t in old new; do echo "mA$m $t $(grep 'Expected 50.0%' ${m}_${t}.log | awk '{print $NF}') $(grep 'Expected 16.0%' ${m}_${t}.log | awk '{print $NF}') $(grep 'Expected 84.0%' ${m}_${t}.log | awk '{print $NF}')"; done >> $H/limits.txt
  cd $FF/..
done
date > $H/.done
