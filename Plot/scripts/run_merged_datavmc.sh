#!/usr/bin/env bash
# dataVmc after the MERGED-photon selection, for all three regions.
#
# LCG_104 (ROOT 6.28), not hza_ana and never system ROOT:
#   * hza_ana is ROOT 6.34, which rejects `from ROOT import *` -- Analyzer_Configs
#     does exactly that, so the plotter dies on import;
#   * system ROOT 6.40 breaks RDataFrame.AsNumpy, which the fast-fill path uses,
#     and the failure is silent (plots simply do not update);
#   * ROOT 6.34 also drops histograms when writing PDFs.
#
# Usage: bash run_merged_datavmc.sh [region ...]     (default 0 1 2)

set -o pipefail

REGIONS="${*:-0 1 2}"
IN=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_mergedflag
OUT_BASE=/eos/home-p/pelai/HZa/output_plots/merged_dataVmc
LOG_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged

# Sub-GeV signal points overlaid on the plots (flashgg spelling: 0p5 == 0.5 GeV).
SIG_OVERLAY="${SIG_OVERLAY:-M0p1,M0p2,M0p5,M0p9}"

# linear | log | both.  The sub-GeV signal is 3-56 weighted events against a
# background peaking in the thousands, so on a linear axis the overlay is drawn
# but unreadable -- the resolved dataVmc plots have always been log for the same
# reason.  Log output is suffixed `_log` and does not overwrite the linear PDFs.
YSCALE="${YSCALE:-both}"
case "${YSCALE}" in
    linear) MODES=("") ;;
    log)    MODES=("--ln") ;;
    both)   MODES=("" "--ln") ;;
    *) echo "YSCALE must be linear|log|both, got '${YSCALE}'" >&2; exit 2 ;;
esac

declare -A NAME=([0]=all [1]=SR [2]=CR)

set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
set -u

cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna

for r in ${REGIONS}; do
    out="${OUT_BASE}/region${r}_${NAME[$r]}"
    log="${LOG_BASE}/datavmc_region${r}.log"
    mkdir -p "${out}"
    echo "=== region ${r} (${NAME[$r]}) -> ${out}"
    for mode in "${MODES[@]}"; do
        tag="linear"; [ -n "${mode}" ] && tag="log"
        mlog="${log%.log}_${tag}.log"
        python Plot/scripts/plot_fast_variable_dataVmc.py \
            --mergedOnly --region "${r}" \
            --sig-overlay "${SIG_OVERLAY}" ${mode} \
            --input-dir "${IN}" \
            --out-dir "${out}" > "${mlog}" 2>&1
        rc=$?
        n=$(ls "${out}" 2>/dev/null | wc -l)
        echo "    ${tag}: rc=${rc}  files_in_dir=${n}  log=${mlog}"
    done
done
