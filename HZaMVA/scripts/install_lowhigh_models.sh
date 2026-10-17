#!/usr/bin/env bash
# Install the low-/high-mass BDTs from HZaMVA/scripts/ into HZaMVA/using/ as a SET:
#   model_Za_BDT_<m>_run3.pkl, .json, .meta.json
# 2026-10-01: the retrain chain copied only the .pkl, so using/ held a new model next to the 09-24
# meta.json; plot_mva_diagnostics.py reads the meta and nothing complained. Here all three files are
# copied together and checked: same bytes as the source, written by the same training (mtimes within
# MAXDT seconds) and the pickle's n_features_in_ equals len(meta["features"]).
# Before anything is copied, meta.json "features" of each model must equal, name by name and in order, the feature
# list the scorer feeds it (LOW_/HIGH_FEATURE_COLUMNS in Parquet2Rootfile/Parque2Root_BDT.py). A model of the same
# size but different inputs (e.g. the old pho1Pt_oHm low-mass set) would otherwise be scored without any error.
# Usage: install_lowhigh_models.sh [backup_tag]     (python with xgboost/sklearn must be on PATH)
set -uo pipefail
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
U=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using
P2R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py
TAG=${1:-bak_$(date +%Y%m%d_%H%M)}
MAXDT=600
python - $P2R $S/model_Za_BDT_lowmass_run3.meta.json $S/model_Za_BDT_highmass_run3.meta.json <<'PY' || exit 1
import ast, json, sys
cols = {}
for node in ast.parse(open(sys.argv[1]).read()).body:
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and getattr(node.targets[0], "id", "") in ("LOW_FEATURE_COLUMNS", "HIGH_FEATURE_COLUMNS"):
        cols[node.targets[0].id] = ast.literal_eval(node.value)
ok = True
for key, meta_path in (("LOW_FEATURE_COLUMNS", sys.argv[2]), ("HIGH_FEATURE_COLUMNS", sys.argv[3])):
    feats = json.load(open(meta_path)).get("features")
    if key not in cols or feats != cols[key]:
        ok = False
        print("FEATURE MISMATCH %s\n  meta   : %s\n  scorer : %s (%s)" % (meta_path, feats, cols.get(key), key))
    else:
        print("  features of %s match %s" % (meta_path.split("/")[-1], key))
sys.exit(0 if ok else 1)
PY
for m in lowmass highmass; do
  b=model_Za_BDT_${m}_run3
  for e in pkl json meta.json; do
    [ -s $S/$b.$e ] || { echo "MISSING $S/$b.$e"; exit 1; }
  done
  t0=$(stat -c %Y $S/$b.pkl)
  for e in json meta.json; do
    dt=$(( $(stat -c %Y $S/$b.$e) - t0 )); dt=${dt#-}
    [ $dt -le $MAXDT ] || { echo "$b.$e is ${dt}s away from $b.pkl -- not from the same training"; exit 1; }
  done
  for e in pkl json meta.json; do
    [ -f $U/$b.$e ] && cp -n $U/$b.$e $U/$b.$e.$TAG
    cp -f $S/$b.$e $U/$b.$e
    cmp -s $S/$b.$e $U/$b.$e || { echo "copy of $b.$e differs"; exit 1; }
  done
  python - $U/$b.pkl $U/$b.meta.json <<'PY' || exit 1
import json, pickle, sys
m = pickle.load(open(sys.argv[1], "rb")); meta = json.load(open(sys.argv[2]))
nf, nm = getattr(m, "n_features_in_", None), len(meta["features"])
print("  %s: n_features=%s meta features=%d auc_test=%.4f" % (sys.argv[1].split("/")[-1], nf, nm, meta.get("auc_test", float("nan"))))
sys.exit(0 if nf == nm else 1)
PY
done
echo "installed low/high-mass models as complete sets (backups *.$TAG)"
