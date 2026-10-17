#!/usr/bin/env bash
# Step A: build the production low-mass model set (pkl / xgboost json / meta.json) from the variable-search model
# low_relE2_w123 (relE2 inputs, m_a 1-3 all trained inside 115-133 GeV) and evaluate it with the production diagnostics (AUC: plot_mva_diagnostics.py; KS:
# plot_overtrain_lowmass.py), before anything is installed.
# Input : varsearch_validate/models/model_low_relE2_w123.{pkl,meta.json}
# Output: model_set/model_Za_BDT_lowmass_run3.{pkl,json,meta.json}, model_set/diag/ (plots), logs/model_set_*.log
set -euo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
SRC=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/models/model_low_relE2_w123
SCR=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
O=$D/model_set; B=$O/model_Za_BDT_lowmass_run3
cp $SRC.pkl $B.pkl
$PY - $B.pkl $SRC.meta.json $B.json $B.meta.json <<'PY'
import json, pickle, sys
clf = pickle.load(open(sys.argv[1], "rb")); src = json.load(open(sys.argv[2]))
clf.save_model(sys.argv[3])
meta = {"features": src["features"], "masses": src["masses"], "window": src["mA1_window"], "window_masses": src["window_masses"], "mA1_window": src["mA1_window"],
        "model": "XGB depth4/300/lr0.1/subsample0.8/mcw5", "name": src["model"], "inputs": src["inputs"],
        "dtype": src["dtype"], "xgboost": src["xgboost"],
        "training": "HZaMVA/scripts/train_lowmass_final.py (reproduces varsearch_validate/models/model_low_relE2_w123 bit-exact, retrain_check/compare.txt)"}
json.dump(meta, open(sys.argv[4], "w"), indent=2)
assert clf.n_features_in_ == len(meta["features"]), (clf.n_features_in_, meta["features"])
PY
cd $SCR
env -i HOME=$HOME PATH=/usr/bin:/bin $PY plot_mva_diagnostics.py --model $B.pkl --meta $B.meta.json --outdir $O/diag > $D/logs/model_set_diag.log 2>&1
env -i HOME=$HOME PATH=/usr/bin:/bin HZA_LOWMODEL_PKL=$B.pkl HZA_LOWMODEL_META=$B.meta.json HZA_LOWMODEL_PLOT=$O/diag/overtrain_lowmass_run3.pdf \
  $PY plot_overtrain_lowmass.py > $D/logs/model_set_overtrain.log 2>&1
$PY - $B.meta.json $D/logs/model_set_diag.log $D/logs/model_set_overtrain.log <<'PY'
import json, re, sys
meta = json.load(open(sys.argv[1]))
m = re.search(r"Test AUC=([0-9.]+), Validation AUC=([0-9.]+)", open(sys.argv[2]).read())
k = re.search(r"capped-KS p: sig=([0-9.e+-]+)\s+bkg=([0-9.e+-]+)", open(sys.argv[3]).read())
meta.update(auc_test=float(m.group(1)), auc_validation=float(m.group(2)), ks_sig=float(k.group(1)), ks_bkg=float(k.group(2)))
json.dump(meta, open(sys.argv[1], "w"), indent=2)
print(json.dumps(meta, indent=2))
PY
touch $B.pkl $B.json   # the three files of the set carry the same write time (install_lowhigh_models.sh check)
echo "MODEL SET DONE"
