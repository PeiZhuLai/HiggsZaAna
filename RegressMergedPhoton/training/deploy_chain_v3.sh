#!/bin/bash
# Finish the v3 deployment: validate the M0p5 pilot, run the remaining sub-GeV
# mass points through HiggsDNA, then produce the two numbers that decide whether
# the retraining was worth it -- the bias_scan table and the re-derived ROI
# windows.
#
# Runs detached so it survives the session:
#     setsid nohup bash deploy_chain_v3.sh > /dev/null 2>&1 &
#
# Log     : <EOS>/deploy_chain_v3.log
# Products: parquet_merged_DNA_v3/Sig_MC_MLNANO_all/
#           <EOS>/bias_scan_v3.txt , <EOS>/roi_windows_v3.txt
set -o pipefail

EOSB=/eos/cms/store/group/phys_susy/pelai/HZa_merged
HDNA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
TRAIN=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/training
NEW="${EOSB}/parquet_merged_DNA_v3/Sig_MC_MLNANO_all"
OLD="${EOSB}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all"
LOG="${EOSB}/deploy_chain_v3.log"
PY_ANA=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
REST="${REST:-M0p1 M0p2 M0p3 M0p4 M0p6 M0p7 M0p8 M0p9}"

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${LOG}"; }
say "chain started"

# ---- 1. wait for the M0p5 pilot -------------------------------------------
# Wait on the PRODUCT, not on a process name. `pgrep -f run_merged_signal_v3.sh`
# also matches any interactive shell whose command line happens to contain that
# string -- including the ones used to inspect progress -- so the loop never
# exited even though the wrapper had long finished. Waiting for the merged
# parquet to appear and stop growing has no such ambiguity.
pilot="${NEW}/mA_MLNANO_M0p5_2024/merged_nominal.parquet"
prev=-1
while true; do
    if [ -f "${pilot}" ]; then
        cur=$(stat -c%s "${pilot}" 2>/dev/null || echo 0)
        # same size on two consecutive polls -> the merge has finished writing
        [ "${cur}" -gt 0 ] && [ "${cur}" = "${prev}" ] && break
        prev="${cur}"
    fi
    sleep 60
done
say "pilot parquet settled ($(stat -c%s "${pilot}") bytes)"
if [ ! -f "${pilot}" ]; then
    say "FATAL: pilot produced no merged_nominal.parquet at ${pilot}"
    say "       check ${EOSB}/higgsdna_v3_M0p5.log"
    exit 1
fi

# The parquet must actually carry the MLPhoton columns the ROI depends on --
# a run that silently drops the friend would still write a file.
"${PY_ANA}" - "${pilot}" <<'PYEOF' >> "${LOG}" 2>&1
import sys, numpy as np, pyarrow.parquet as pq
p = sys.argv[1]
names = pq.read_schema(p).names
need = ["MLPhoton_lead_mass", "pass_allcuts_merged_ML", "MergedML_mass"]
missing = [c for c in need if c not in names]
print(f"[validate] {p}")
print(f"[validate] columns={len(names)} missing={missing}")
if missing:
    raise SystemExit(2)
t = pq.read_table(p, columns=["MLPhoton_lead_mass", "pass_allcuts_merged_ML"])
v = t["MLPhoton_lead_mass"].to_numpy(); s = t["pass_allcuts_merged_ML"].to_numpy()
m = s & np.isfinite(v) & (v > -100)
print(f"[validate] rows={len(v)} selected={int(m.sum())} "
      f"MLPhoton_lead_mass median={np.median(v[m]):.4f}")
if m.sum() < 20:
    raise SystemExit(3)
PYEOF
rc=$?
if [ ${rc} -ne 0 ]; then
    say "FATAL: pilot validation failed (rc=${rc}); not launching the rest"
    exit 1
fi
say "pilot validated"

# ---- 2. the remaining sub-GeV points --------------------------------------
say "running HiggsDNA for: ${REST}"
bash "${HDNA}/scripts/run_merged_signal_v3.sh" ${REST} \
    >> "${EOSB}/higgsdna_v3_rest.log" 2>&1
say "HiggsDNA rest exit code: $?"

produced=$(ls -d "${NEW}"/mA_MLNANO_M0p*_2024 2>/dev/null | wc -l)
say "mass-point dirs produced: ${produced}/9"

# ---- 3. the two deliverables ----------------------------------------------
say "bias_scan on the new production"
"${PY_ANA}" "${TRAIN}/bias_scan.py" --base "${NEW}" \
    > "${EOSB}/bias_scan_v3.txt" 2>&1
say "bias_scan (old production, same script) for reference"
"${PY_ANA}" "${TRAIN}/bias_scan.py" --base "${OLD}" \
    > "${EOSB}/bias_scan_old.txt" 2>&1

say "re-deriving ROI windows"
"${PY_ANA}" "${TRAIN}/derive_roi_windows.py" --base "${NEW}" --compare "${OLD}" \
    > "${EOSB}/roi_windows_v3.txt" 2>&1

say "chain finished"
say "  bias (new): ${EOSB}/bias_scan_v3.txt"
say "  bias (old): ${EOSB}/bias_scan_old.txt"
say "  ROI       : ${EOSB}/roi_windows_v3.txt"
