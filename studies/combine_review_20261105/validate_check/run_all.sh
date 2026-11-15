#!/usr/bin/env bash
# ValidateDatacards.py (CI container) on all 30 mass points x both cards (Wilson, cross-section) at the default
# thresholds (-c 0.1, -s 0.001). The tool joins the card directory to the shapes paths, so each card is copied here
# with relative shapes paths; nuisance names are left as in production.
# For a parametric model only three checks apply (large normalization effect, empty process, small signal);
# 10_card.txt with -c 0.01 -s 0.2 (lowthr.json) showed that they do run on these cards.
# Input : flashggFinalFit/Datacard/output_Datacard_leptons/<m>_pruned_datacard_leptons{,_xsec}.txt (read only)
# Output: cards/<m>_<card>.txt, json/<m>_<card>.json, logs/<m>_<card>.log, summary.txt, .done
set -o pipefail
V=$(cd "$(dirname "$0")" && pwd)
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons
IMG=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/combine-container:CMSSW_14_1_0_pre4-combine_v10.6.1-harvester_v3.1.0
mkdir -p $V/cards $V/json $V/logs; rm -f $V/.done
for m in $(seq 1 30); do
  for c in wilson xsec; do
    [ $c = xsec ] && src=$FF/${m}_pruned_datacard_leptons_xsec.txt || src=$FF/${m}_pruned_datacard_leptons.txt
    python3 - $src $V/cards/${m}_${c}.txt <<'EOF'
import os, sys
d = os.path.dirname(sys.argv[2]); out = []
for l in open(sys.argv[1]):
    f = l.split()
    if f and f[0] == 'shapes' and f[3].startswith('/'): l = l.replace(f[3], os.path.relpath(f[3], d), 1)
    out.append(l)
open(sys.argv[2], 'w').writelines(out)
EOF
  done
done
apptainer exec -B /afs -B /eos -B /cvmfs $IMG bash -c "source /cvmfs/cms.cern.ch/cmsset_default.sh 2>/dev/null; cd /home/cmsusr/CMSSW_14_1_0_pre4/src && eval \$(scramv1 runtime -sh) 2>/dev/null; cd $V/cards
for f in *.txt; do b=\${f%.txt}; ValidateDatacards.py \$f --jsonFile $V/json/\$b.json -p 2 > $V/logs/\$b.log 2>&1; echo \"\$b rc=\$? \$(grep -c . $V/json/\$b.json)\"; done" > $V/summary_rc.txt 2>&1
python3 - $V <<'EOF'
import glob, json, os, sys
V = sys.argv[1]; out = []
for j in sorted(glob.glob(os.path.join(V, 'json', '*.json')), key=lambda p: (int(os.path.basename(p).split('_')[0]), p)):
    d = json.load(open(j)); b = os.path.basename(j)[:-5]
    if not d: out.append('%-10s null (no warnings)' % b); continue
    for k, v in d.items():
        out.append('%-10s %s: %s' % (b, k, ', '.join(sorted(v)) if isinstance(v, dict) else v))
open(os.path.join(V, 'summary.txt'), 'w').write('\n'.join(out) + '\n')
EOF
date > $V/.done
