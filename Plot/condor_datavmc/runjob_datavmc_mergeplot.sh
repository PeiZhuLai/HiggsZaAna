#!/bin/bash
# HTCondor worker job: merge (hadd) the 6 sample partials for one (region, tag) and draw
# the data/MC plots. hadd is done on the worker's LOCAL scratch (AFS has pathological
# per-object latency from remote workers: hadd of ~2900 histograms over AFS took >1h;
# local hadd is seconds). Only bulk file copies touch AFS.
# args: REGION_KEY(SR|CR|mva) FINAL_TAG(nominal|sideband_rwgt)
set -eo pipefail

region_key="$1"; final_tag="$2"
if [[ -z "$region_key" || -z "$final_tag" ]]; then
    echo "[ERROR] usage: $0 REGION TAG"; exit 2
fi

PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
scriptsDir="$PLOT/scripts"
V="$PLOT/plots/variables_dataVmc"

case "$region_key" in
    SR)  suffix="UL_SR";  region_arg=(--region 1) ;;
    CR)  suffix="UL_CR";  region_arg=(--region 2) ;;
    mva) suffix="UL_mva"; region_arg=(-b) ;;
    *)   echo "[ERROR] unknown region $region_key"; exit 2 ;;
esac

source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
export PYTHONPATH="${PYTHONPATH:-}:$PLOT/lib"

scratch="${_CONDOR_SCRATCH_DIR:-/tmp}/dvmc_${suffix}_${final_tag}_$$"
mkdir -p "$scratch"
trap 'rm -rf "$scratch"' EXIT

# --- bulk-copy the 6 partials to local scratch, then hadd locally (low latency) ---
local_inputs=()
for st in data dyll dyg sig_m1_m5 sig_m6_m10 sig_m15_m30; do
    src="$V/ALP_plot_run3_${suffix}_${final_tag}_part_${st}.root"
    if [[ ! -s "$src" ]]; then echo "[ERROR] missing partial $src"; exit 3; fi
    cp "$src" "$scratch/part_${st}.root"
    local_inputs+=("$scratch/part_${st}.root")
done

local_merged="$scratch/ALP_plot_run3_${suffix}_${final_tag}.root"
echo "[mergeplot] local hadd -> $local_merged"
hadd -f "$local_merged" "${local_inputs[@]}"

# --- publish merged back to AFS (single bulk copy) ---
target="$V/ALP_plot_run3_${suffix}_${final_tag}.root"
cp "$local_merged" "$target"
echo "[mergeplot] published merged -> $target"

# --- draw (reads the small merged file; writes ~19 pdfs to AFS) ---
cd "$PLOT"
echo "[mergeplot] draw region=$region_key tag=$final_tag"
python3 "$scriptsDir/2_plot_dataVmc.py" -y run3 -m --ln \
    --inputTag "$final_tag" --outputTag "$final_tag" "${region_arg[@]}"
echo "[mergeplot][DONE] $(date '+%F %T')"
