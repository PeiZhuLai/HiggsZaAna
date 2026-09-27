#!/usr/bin/env bash
# Repoint HZaMVA from run3_bdt_inputs_nominal to the FSR-fix ROOT files.
#
# Nine files hardcode the nominal directory. Rather than swap the literal (which makes
# it impossible to go back and impossible to tell which production a result came from),
# each becomes an env-overridable default: HZA_P2ROOT_BASE wins, otherwise fsrfix.
# Set HZA_P2ROOT_BASE=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal to
# reproduce anything from the pre-FSR-fix era.
set -uo pipefail
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA
OLD=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal
NEW=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix
PYFILES="scripts/1_make_sideband_reweight.py scripts/sculpt_diag.py scripts/run3_Za_BDT.py
         scripts/sculpt_retrain_test.py scripts/hza_features.py scripts/sculpt_param_test.py
         scripts/disco_prep.py"
for f in $PYFILES; do
  [ -f "$f" ] || { echo "skip (absent): $f"; continue; }
  cp -n "$f" "${f}.bak_nominal_20260921"
  # ensure `os` is importable where we inject os.environ
  grep -q '^import os' "$f" || sed -i '0,/^import /s//import os\nimport /' "$f"
  sed -i "s|\"$OLD|os.environ.get(\"HZA_P2ROOT_BASE\", \"$NEW\") + \"|g; s|'$OLD|os.environ.get('HZA_P2ROOT_BASE', '$NEW') + '|g" "$f"
  python3 -c "import ast,sys;ast.parse(open('$f').read())" \
    && echo "patched  $f" || { echo "SYNTAX BROKEN: $f -- restoring"; cp "${f}.bak_nominal_20260921" "$f"; }
done
for f in scripts/run_apply_nn_inputs.sh; do
  [ -f "$f" ] || continue
  cp -n "$f" "${f}.bak_nominal_20260921"
  sed -i "s|INBASE=\"$OLD\"|INBASE=\"\${HZA_P2ROOT_BASE:-$NEW}\"|" "$f"
  bash -n "$f" && echo "patched  $f"
done
echo "=== remaining hardcoded nominal references ==="
grep -rn "run3_bdt_inputs_nominal" --include=*.py --include=*.sh . 2>/dev/null | grep -v '\.bak' || echo "  (none)"
