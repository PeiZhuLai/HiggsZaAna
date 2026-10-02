#!/usr/bin/env bash
# Stage 4b: retrain the two models the SCORING stage actually uses.
#
# Parque2Root_BDT.add_mva_scores() routes ma 1-3 to model_Za_BDT_lowmass_run3.pkl and
# ma 4-30 to model_Za_BDT_highmass_run3.pkl. It never calls get_model(), so the single
# model_Za_BDT_run3.pkl that 2_train.sh produces is NOT what scoring uses -- that one is
# consumed by Plot/lib/Analyzer_Configs.py for dataVmc. All three must be retrained on
# the FSR-fix samples, and it is easy to retrain only the first and think you are done.
#
# These scripts write their pkl into scripts/ and do NOT copy to using/. The copy is
# explicit here, after the models have been verified.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
U=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using
START=$(date +%s)

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
export PYTHONPATH="${PYTHONPATH:-}:$S"
echo "[env] ROOT_DIR seen by hza_features: $(cd $S && python -c 'from hza_features import ROOT_DIR; print(ROOT_DIR)')"

for m in lowmass highmass; do
  for ext in pkl json meta.json; do
    cp -n "$S/model_Za_BDT_${m}_run3.${ext}" "$S/model_Za_BDT_${m}_run3.${ext}.bak_preFSRfix_20260922" 2>/dev/null
  done
  cp -n "$U/model_Za_BDT_${m}_run3.pkl" "$U/model_Za_BDT_${m}_run3.pkl.bak_preFSRfix_20260922" 2>/dev/null
done
echo "[backup] July low/high models saved with .bak_preFSRfix_20260922"

cd "$S"
echo; echo "=============== [1/4] train low-mass (ma 1-3) ==============="
python train_lowmass_final.py
echo "raw exit: $? (not a verdict)"

echo; echo "=============== [2/4] train high-mass (ma 4-30) ==============="
python train_highmass_final.py
echo "raw exit: $? (not a verdict)"

echo; echo "=============== [3/4] install into using/ ==============="
for m in lowmass highmass; do
  src="$S/model_Za_BDT_${m}_run3.pkl"
  if [ ! -s "$src" ]; then echo "MISSING $src -- refusing to install"; exit 1; fi
  if [ "$(stat -c %Y "$src")" -lt "$START" ]; then echo "STALE $src -- refusing to install"; exit 1; fi
done
# 2026-10-01: install pkl + json + meta.json together (only the pkl used to be copied)
bash "$S/install_lowhigh_models.sh" "bak_install_$(date +%Y%m%d_%H%M)" || { echo "model install failed"; exit 1; }

echo; echo "=============== [4/4] verify all THREE models ==============="
python - <<PY
import os, pickle, sys, time, json
U="$U"; S="$S"; START=$START
expect = {"model_Za_BDT_run3.pkl": 16,
          "model_Za_BDT_lowmass_run3.pkl": 5,
          "model_Za_BDT_highmass_run3.pkl": 16}
ok = True
for name, nfeat_exp in expect.items():
    p = os.path.join(U, name)
    if not os.path.exists(p):
        print("MISSING", p); ok = False; continue
    mt = os.path.getmtime(p)
    fresh = mt >= (START - 6*3600)   # the single model was trained earlier this morning
    try:
        with open(p, "rb") as f: m = pickle.load(f)
        nf = getattr(m, "n_features_in_", None)
    except Exception as e:
        print("UNREADABLE %s -> %s" % (p, str(e)[:80])); ok = False; continue
    good = fresh and (nf == nfeat_exp)
    print("%-34s %-6s mtime=%s n_features=%s (expect %d)"
          % (name, "OK" if good else "BAD", time.strftime('%F %T', time.localtime(mt)), nf, nfeat_exp))
    if not good: ok = False
print()
for m in ("lowmass", "highmass"):
    f = os.path.join(S, "model_Za_BDT_%s_run3.meta.json" % m)
    if os.path.exists(f):
        d = json.load(open(f))
        print("%-9s AUC(test)=%.4f  ks_sig=%.3g ks_bkg=%.3g" % (m, d["auc_test"], d["ks_sig"], d["ks_bkg"]))
        print("          R: " + "  ".join("mA%s=%.2f" % (k, v) for k, v in d["R"].items()))
sys.exit(0 if ok else 1)
PY
rc=$?
echo "VERDICT: $([ $rc -eq 0 ] && echo PASSED || echo FAILED)"
exit $rc
