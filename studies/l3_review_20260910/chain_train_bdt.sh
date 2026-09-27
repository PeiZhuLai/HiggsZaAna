#!/usr/bin/env bash
# Stage 4: retrain the BDT on the FSR-fix samples, then save the model into HZaMVA/using.
#
# Retraining is the right call here, not an optional extra: the signal m_llgg
# distribution moved (sigma_eff improved 2.69% in the muon channel, 0.00% in the
# electron control), the 2022 DYGto2LG background is a different MC sample with 24-31%
# more statistics, and the sideband reweight the trainer applies now comes from the
# regenerated JSON. Scoring the new samples with the July model would train and apply on
# different distributions.
#
# The old model is kept, not overwritten in place -- 3_save_model.sh copies into
# ../using/, so the previous one is backed up first.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
U=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using
START=$(date +%s)

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"

cp -n "$S/model_Za_BDT_run3.pkl" "$S/model_Za_BDT_run3.pkl.bak_preFSRfix_20260922" 2>/dev/null
cp -n "$U/model_Za_BDT_run3.pkl" "$U/model_Za_BDT_run3.pkl.bak_preFSRfix_20260922" 2>/dev/null
echo "[backup] previous models saved with .bak_preFSRfix_20260922"

echo "=============== [1/3] train ==============="
cd "$S"
bash 2_train.sh
echo "raw exit: $? (not a verdict)"

echo "=============== [2/3] save model ==============="
bash 3_save_model.sh
echo "raw exit: $? (not a verdict)"

echo "=============== [3/3] verify the product, not the exit code ==============="
python - <<PY
import os, pickle, sys, time
S="$S"; U="$U"; START=$START
ok=True
for p in (S+"/model_Za_BDT_run3.pkl", U+"/model_Za_BDT_run3.pkl"):
    if not os.path.exists(p):
        print("MISSING", p); ok=False; continue
    mt=os.path.getmtime(p)
    fresh = mt >= START
    try:
        with open(p,"rb") as f: m=pickle.load(f)
        kind=type(m).__name__
        nfeat=getattr(m,"n_features_in_",None)
    except Exception as e:
        print("UNREADABLE %s -> %s"%(p,str(e)[:80])); ok=False; continue
    print("%-58s %s  mtime=%s fresh=%s  type=%s n_features=%s"
          % (p, "OK" if fresh else "STALE", time.strftime('%F %T',time.localtime(mt)), fresh, kind, nfeat))
    if not fresh: ok=False
sys.exit(0 if ok else 1)
PY
rc=$?
echo "VERDICT: $([ $rc -eq 0 ] && echo PASSED || echo FAILED)"
exit $rc
