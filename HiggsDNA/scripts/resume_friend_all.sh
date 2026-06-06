#!/usr/bin/env bash
# Resume the friend-tree HiggsDNA batch after the original driver
# (run_merged_friend_all.sh, PID 1715011) died mid-mA_M0p4 on 2026-06-02.
#
# mA_M0p1/M0p2/M0p3 are already done+merged. mA_M0p4 has 22 jobs still idle in
# the condor queue plus an intact analysis_manager.pkl, so we RECONNECT to those
# jobs by re-invoking run_analysis.py WITHOUT deleting the pkl (do NOT use the
# run_merged_friend_all.sh wrapper for M0p4 — it rm's the pkl and would orphan
# the 22 jobs and resubmit duplicates).
#
# After M0p4 merges, hand off to the normal driver for the remaining tags.

set -uo pipefail

repo_dir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
cd "${repo_dir}"

unset PYTHONPATH
export PYTHONPATH="${repo_dir}"

OUTDIR_BASE="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend"

# --- Step 1: reconnect to the in-flight mA_M0p4 jobs (pkl preserved) ---
m0p4_out="${OUTDIR_BASE}/mA_M0p4"
echo "=== RESUME mA_M0p4 (reconnect to existing pkl + idle condor jobs) ==="
python scripts/run_analysis.py \
    --config metadata/friend_tmp_configs/mA_M0p4.json \
    --log-level INFO \
    --n_cores 4 \
    --output_dir "${m0p4_out}" \
    --batch_system condor \
    --unretire_jobs \
    --merge_outputs \
    --fpo 1
m0p4_rc=$?
echo "=== mA_M0p4 run_analysis.py exit code: ${m0p4_rc} ==="

# Sanity: confirm M0p4 produced merged parquet before moving on.
n_merged=$(ls "${m0p4_out}"/*/merged_*.parquet 2>/dev/null | wc -l)
echo "=== mA_M0p4 merged parquet count: ${n_merged} ==="
if [ "${n_merged}" -eq 0 ]; then
    echo "[WARN] mA_M0p4 produced no merged parquet; continuing to remaining tags anyway."
fi

# --- Step 2: remaining tags via the normal driver (fresh, pkl-rm is fine) ---
REMAINING="mA_M0p5,mA_M0p6,mA_M0p7,mA_M0p8,mA_M0p9,mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,Bkg_DYGto2LG_10to100_2024,Bkg_DYJetsTo2E_2024,Bkg_DYJetsTo2Mu_2024,Bkg_DYJetsTo2Tau_2024"
echo "=== HANDOFF to driver for remaining tags: ${REMAINING} ==="
TAGS="${REMAINING}" bash scripts/run_merged_friend_all.sh

echo "=== resume_friend_all.sh finished ==="
