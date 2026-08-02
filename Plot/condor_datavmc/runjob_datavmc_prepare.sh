#!/bin/bash
# HTCondor worker job: one dataVmc prepare (1_prepare_dataVmc.py) for a
# (region, tag, sample-group). Mirrors run_plot_task in A_run_dataVmc.sh but for batch,
# so the heavy RDataFrame event loop runs on a worker node (not the lxplus login arbiter).
# args: REGION_KEY(SR|CR|mva) FINAL_TAG(nominal|sideband_rwgt) SAMPLE_TAG SAMPLES
set -eo pipefail

region_key="$1"; final_tag="$2"; sample_tag="$3"; samples="$4"
if [[ -z "$region_key" || -z "$final_tag" || -z "$sample_tag" || -z "$samples" ]]; then
    echo "[ERROR] usage: $0 REGION TAG SAMPLE_TAG SAMPLES"; exit 2
fi

PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
scriptsDir="$PLOT/scripts"
sidebandReweightJson=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json

# --- environment: higgs-alp-ana conda (uproot + ROOT + RDataFrame) ---
source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
export PYTHONPATH="${PYTHONPATH:-}:$PLOT/lib"
# keep the RDataFrame event loop single-threaded (request_cpus=1)
export ALP_RDF_THREADS=0

partial_tag="${final_tag}_part_${sample_tag}"
cmd=(python3 "$scriptsDir/1_prepare_dataVmc.py" -y run3 -m --ln --histOnly \
        --samples "$samples" --outputTag "$partial_tag")
case "$region_key" in
    SR)  cmd+=(--region 1) ;;
    CR)  cmd+=(--region 2) ;;
    mva) cmd+=(-b) ;;
    *)   echo "[ERROR] unknown region $region_key"; exit 2 ;;
esac
if [[ "$final_tag" == "sideband_rwgt" ]]; then
    cmd+=(--useSidebandReweight --sidebandReweightJson "$sidebandReweightJson")
fi

cd "$PLOT"
echo "[runjob] host=$(hostname) region=$region_key tag=$final_tag sample_tag=$sample_tag"
printf '[runjob][CMD]'; printf ' %q' "${cmd[@]}"; echo
"${cmd[@]}"
echo "[runjob][DONE] $(date '+%F %T')"
