#!/usr/bin/env bash
# Reproduction test of the production low-mass trainer (HZaMVA/scripts/train_lowmass_final.py, w123 recipe):
# retrain into retrain_check/ (production model files untouched) and compare with the sandbox w123 model (varsearch_validate/models/model_low_relE2_w123.pkl)
# : booster bytes, pkl md5, scores on the training matrix.
# Output: retrain_check/model_Za_BDT_lowmass_run3.{pkl,json,meta.json}, retrain_check/compare.txt
# Log: logs/retrain_check.log; completion marker: retrain_check/DONE (or FAIL)
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
SCR=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
U=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/models/model_low_relE2_w123.pkl
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
O=$D/retrain_check; rm -f $O/DONE $O/FAIL
cd $SCR
env -i HOME=$HOME PATH=/usr/bin:/bin HZA_LOWMODEL_OUT=$O/model_Za_BDT_lowmass_run3 $PY train_lowmass_final.py || { touch $O/FAIL; exit 1; }
$PY - $O/model_Za_BDT_lowmass_run3.pkl $U > $O/compare.txt 2>&1 <<'PY'
import hashlib, pickle, sys
new, old = (pickle.load(open(p, "rb")) for p in sys.argv[1:3])
md5 = [hashlib.md5(open(p, "rb").read()).hexdigest() for p in sys.argv[1:3]]
rb = [m.get_booster().save_raw("json") for m in (new, old)]
print("pkl md5      new %s  sandbox %s  equal %s" % (md5[0], md5[1], md5[0] == md5[1]))
print("booster json equal %s (%d / %d bytes)" % (rb[0] == rb[1], len(rb[0]), len(rb[1])))
print("RESULT", "BIT-EXACT" if rb[0] == rb[1] else "DIFFERENT")
PY
cat $O/compare.txt
grep -q "RESULT BIT-EXACT" $O/compare.txt && touch $O/DONE || touch $O/FAIL
