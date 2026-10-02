#!/usr/bin/env bash
# Install the low-/high-mass BDTs from HZaMVA/scripts/ into HZaMVA/using/ as a SET:
#   model_Za_BDT_<m>_run3.pkl, .json, .meta.json
# 2026-10-01: the retrain chain copied only the .pkl, so using/ held a new model next to the 09-24
# meta.json; plot_mva_diagnostics.py reads the meta and nothing complained. Here all three files are
# copied together and checked: same bytes as the source, written by the same training (mtimes within
# MAXDT seconds) and the pickle's n_features_in_ equals len(meta["features"]).
# Usage: install_lowhigh_models.sh [backup_tag]     (python with xgboost/sklearn must be on PATH)
set -uo pipefail
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
U=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using
TAG=${1:-bak_$(date +%Y%m%d_%H%M)}
MAXDT=600
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
print("  %s: n_features=%s meta features=%d auc_test=%.4f" % (sys.argv[1].split("/")[-1], nf, nm, meta["auc_test"]))
sys.exit(0 if nf == nm else 1)
PY
done
echo "installed low/high-mass models as complete sets (backups *.$TAG)"
