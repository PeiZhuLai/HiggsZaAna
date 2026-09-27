#!/bin/bash
# Full v4 validation chain: closure -> signal MLNanoAOD -> HiggsDNA -> bias/ROI.
#
# Everything the earlier rounds taught, applied:
#   * wait on PRODUCTS, never on `pgrep -f <script>` (that also matches the
#     interactive shells used to check progress, so the wait never ends)
#   * validate the pilot's COLUMNS before launching the rest -- a missing friend
#     join writes a perfectly valid parquet with no MLPhoton_* in it
#   * one mass point failing must not skip the others (DAS is flaky)
#   * outputs go to _v4 paths; nothing existing is overwritten
#
# Runs detached:
#     setsid nohup bash deploy_chain_v4.sh > /dev/null 2>&1 &
#
# Log     : <EOS>/deploy_chain_v4.log
# Products: <EOS>/bias_scan_v4.txt , <EOS>/roi_windows_v4.txt
set -o pipefail

EOSB=/eos/cms/store/group/phys_susy/pelai/HZa_merged
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
TRAIN="${REPO}/RegressMergedPhoton/training"
CONDOR="${REPO}/RegressMergedPhoton/condor"
HDNA="${REPO}/HiggsDNA"
ML_V4="${EOSB}/MLNanoAOD_v4"
PQ_V4="${EOSB}/parquet_merged_DNA_v4/Sig_MC_MLNANO_all"
OLD_PQ="${EOSB}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all"
V3_PQ="${EOSB}/parquet_merged_DNA_v3/Sig_MC_MLNANO_all"
LOG="${EOSB}/deploy_chain_v4.log"
PY_ANA=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
MASSES="M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9"

say() { echo "[$(date '+%F %T')] $*" >> "${LOG}"; }
say "=========== v4 chain started ==========="

# ---- 0. closure: does CMSSW actually load the new regressor? ---------------
say "step 0: preprocessing closure with the installed v4 models"
bash "${TRAIN}/run_closure.sh" M1 200 > "${EOSB}/closure_v4.log" 2>&1
if grep -q "RESULT: PASS" "${EOSB}/closure_v4.log"; then
    say "closure PASS"
else
    say "FATAL: closure did not pass -- see ${EOSB}/closure_v4.log"
    exit 1
fi

# ---- 1. signal MLNanoAOD ---------------------------------------------------
say "step 1: re-producing signal MLNanoAOD -> ${ML_V4}"
OUT_EOS_BASE="${ML_V4}" bash "${CONDOR}/submit_all_signal_v3.sh" \
    >> "${EOSB}/mlnano_v4_submit.log" 2>&1
say "submitted; waiting for the queue to drain"

released=0
while true; do
    st=$(condor_q -constraint 'regexp("MLPhoton_Regressing_", JobBatchName)' -af JobStatus 2>/dev/null)
    run=$(echo "${st}" | grep -c '^2$' || true)
    idle=$(echo "${st}" | grep -c '^1$' || true)
    held=$(echo "${st}" | grep -c '^5$' || true)
    n=$(find "${ML_V4}" -name '*.root' 2>/dev/null | wc -l)
    say "  mlnano: running=${run} idle=${idle} held=${held} files=${n}"
    if [ "${held}" -gt 0 ] && [ "${released}" -eq 0 ]; then
        condor_release -constraint 'regexp("MLPhoton_Regressing_", JobBatchName)' >> "${LOG}" 2>&1
        released=1
        say "  released ${held} held job(s) (one-shot)"
    fi
    [ "${run}" -eq 0 ] && [ "${idle}" -eq 0 ] && break
    sleep 600
done
n_ml=$(find "${ML_V4}" -name '*.root' 2>/dev/null | wc -l)
say "step 1 done: ${n_ml} MLNanoAOD files (the v3 round produced 529)"
if [ "${n_ml}" -lt 500 ]; then
    say "FATAL: MLNanoAOD is short (${n_ml} < 500); not continuing"
    exit 1
fi

# ---- 2. HiggsDNA, pilot first ---------------------------------------------
say "step 2a: HiggsDNA pilot (M0p5)"
ML_BASE="${ML_V4}" OUTDIR="${PQ_V4}" \
    bash "${HDNA}/scripts/run_merged_signal_v3.sh" M0p5 \
    > "${EOSB}/higgsdna_v4_M0p5.log" 2>&1

pilot="${PQ_V4}/mA_MLNANO_M0p5_2024/merged_nominal.parquet"
prev=-1
while true; do
    if [ -f "${pilot}" ]; then
        cur=$(stat -c%s "${pilot}" 2>/dev/null || echo 0)
        [ "${cur}" -gt 0 ] && [ "${cur}" = "${prev}" ] && break
        prev="${cur}"
    fi
    sleep 60
done
say "pilot parquet settled ($(stat -c%s "${pilot}") bytes)"

"${PY_ANA}" - "${pilot}" <<'PYEOF' >> "${LOG}" 2>&1
import sys, numpy as np, pyarrow.parquet as pq
p = sys.argv[1]
names = pq.read_schema(p).names
need = ["MLPhoton_lead_mass", "pass_allcuts_merged_ML", "MergedML_mass"]
missing = [c for c in need if c not in names]
print(f"[validate] columns={len(names)} missing={missing}")
if missing:
    raise SystemExit(2)
t = pq.read_table(p, columns=["MLPhoton_lead_mass", "pass_allcuts_merged_ML"])
v = t["MLPhoton_lead_mass"].to_numpy(); s = t["pass_allcuts_merged_ML"].to_numpy()
m = s & np.isfinite(v) & (v > -100)
print(f"[validate] rows={len(v)} selected={int(m.sum())} median={np.median(v[m]):.4f}")
if m.sum() < 20:
    raise SystemExit(3)
PYEOF
if [ $? -ne 0 ]; then
    say "FATAL: pilot validation failed -- not launching the rest"
    exit 1
fi
say "pilot validated"

say "step 2b: HiggsDNA for the remaining mass points"
ML_BASE="${ML_V4}" OUTDIR="${PQ_V4}" \
    bash "${HDNA}/scripts/run_merged_signal_v3.sh" M0p1 M0p2 M0p3 M0p4 M0p6 M0p7 M0p8 M0p9 \
    > "${EOSB}/higgsdna_v4_rest.log" 2>&1
say "step 2b exit: $?"
n_dirs=$(ls -d "${PQ_V4}"/mA_MLNANO_M0p*_2024 2>/dev/null | wc -l)
say "mass-point dirs: ${n_dirs}/9"
grep -E "FAILED mass points|all requested" "${EOSB}/higgsdna_v4_rest.log" >> "${LOG}" 2>&1

# ---- 3. the deliverables ---------------------------------------------------
say "step 3: bias_scan + ROI"
"${PY_ANA}" "${TRAIN}/bias_scan.py" --base "${PQ_V4}" > "${EOSB}/bias_scan_v4.txt" 2>&1
"${PY_ANA}" "${TRAIN}/derive_roi_windows.py" --base "${PQ_V4}" --compare "${OLD_PQ}" \
    > "${EOSB}/roi_windows_v4.txt" 2>&1
# v4 vs v3 directly -- the more useful comparison now that both exist
"${PY_ANA}" "${TRAIN}/derive_roi_windows.py" --base "${PQ_V4}" --compare "${V3_PQ}" \
    > "${EOSB}/roi_windows_v4_vs_v3.txt" 2>&1
say "chain finished"
say "  bias : ${EOSB}/bias_scan_v4.txt"
say "  ROI  : ${EOSB}/roi_windows_v4.txt (vs old), roi_windows_v4_vs_v3.txt (vs v3)"
