#!/bin/bash
# =============================================================================
#  1_grand_merged.sh  --  end-to-end pipeline for the LOW-mA MERGED H->Z a
#                         analysis  (m_a = 0.1 .. 0.9 GeV, a -> single merged
#                         photon),  built on the flashggFinalFit framework
#                         (same machinery as the resolved 1_grand.sh).
#
#  Physics (per-mA ROI, doc sec.12): fit observable = CMS_hza_mass =
#  H_merged_ML_mass = m(ll + MLPhoton_lead), recomputed in merged_p2root.py
#  (peaks ~125). Selection = pass_allcuts_merged_ML. Each m_a opens a ROI window
#  on MLPhoton_lead_mass (reco ML merged-photon mass; tracks m_a) -> a DIFFERENT
#  data sub-sample -> a PER-mA background (like resolved binning by m_gg~m_a).
#  Signal = double-CB with lepton/photon shape systematics; background =
#  data-driven discrete-profiling envelope (RooMultiPdf), fit PER m_a.
#
#  Lumi: only 2024 data exists today; the datacard projects to the full Run-3
#  lumi (2022+2023+2024 = 172.13 fb^-1) via makeYields --lumi (single x1.567).
#
#  Chain (mirrors doc/HZa_merged/low_mA_merged_progress.md sec.11):
#    [0] (opt) produce 2024 merged-selection DATA friend parquet
#    [1] merged_p2root.py   : friend parquet -> flashgg DiphotonTree ROOT trees
#    [2] Tree2WS            : signal (per lep) + data -> RooWorkspace
#    [3] signal model       : fTest -> calcPhotonSyst -> signalFit (DCB + syst)
#    [4] background          : fTest_ALP_turnOn -> PER-mA multipdf (per-mA ROI)
#    [5] datacard           : makeYields(--mergedLowMA --lumi) -> makeDatacard
#                             -> text2workspace
#    [6] combine             : AsymptoticLimits per m_a
#    [7] plot                : makeLimitsPlot_lowmA.py (log-x Brazil band)
#
#  Env split: p2root + limit plot use conda `higgs-alp-ana` (uproot/PyROOT);
#  Tree2WS/signal/bkg/datacard/combine use the flashgg cmsenv.
#
#  The earlier standalone-RooFit prototype (fit_signal.py / build_bkg_ws.py /
#  ...) still lives under MergedAna/ for fast iteration; it is NOT used here.
#
#  Usage:  bash 1_grand_merged.sh
#  Toggle stages with RUN_* env vars (defaults below). Subset masses with MASSES.
# =============================================================================
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
MERGED_DIR="${MERGED_DIR:-${PROJECT_DIR}/MergedAna}"
EOS="${EOS:-/eos/project/h/htozg-dy-privatemc/pelai/HZa}"
CMSSW_SRC="${CMSSW_SRC:-/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src}"
FLASHGG="${FLASHGG:-${CMSSW_SRC}/flashggFinalFit}"
ROOT_MVACUT="${ROOT_MVACUT:-/eos/home-p/pelai/HZa/root_MVAcut}"
# Signal ML parquet (all 9 sub-GeV mass points, nominal+16 syst; MLPhoton_lead_*).
SIG_ML_DIR="${SIG_ML_DIR:-/eos/home-p/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all}"

# Sub-GeV mass points (flashgg label form, e.g. 0p5 == m_a 0.5 GeV) and channels.
MASSES="${MASSES:-0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9}"
LEPS="${LEPS:-ele mu}"
YEAR="${YEAR:-2024}"

# Run-3 lumi projection (172.13 fb^-1) applied at the datacard; 0 => use 2024 only.
RUN3_LUMI="${RUN3_LUMI:-172.13}"
ASSUME_XS="${ASSUME_XS:-100.0}"          # reference xs [fb] for the limit-plot y-axis
GOF_TOYS="${GOF_TOYS:-50}"               # bkg fTest goodness-of-fit toys
MAX_PARALLEL="${MAX_PARALLEL:-4}"        # cap for the per-mass signal/Tree2WS loops
# Full-range overlay: stitch the merged sub-GeV limits with the RESOLVED (m_a>=1)
# combine roots into one 0.1-30 GeV plot. Set RESOLVED_LIMITS_DIR='' to skip it.
RESOLVED_LIMITS_DIR="${RESOLVED_LIMITS_DIR:-${FLASHGG}/Combine/output_combine_results}"

# flashgg outputs
SIG_DIR="${FLASHGG}/Signal"
BKG_WSDIR="${BKG_WSDIR:-${FLASHGG}/Background/ALP_BkgModel_merged}"
# per-mA ROI: bkg multipdf + data live per-mA under ${BKG_WSDIR}/<m>/ and
# ${ROOT_MVACUT}/data/Data_merged_M<m>/ (Stage 4 / Stage 2). BKG_MERGED/DATA_WS
# (the old single shared paths) are retained only for backward-compat references.
BKG_MERGED="${BKG_WSDIR}/merged"
DATA_WS="${ROOT_MVACUT}/data/Data_merged/ws/run3.root"
DC_DIR="${FLASHGG}/Datacard"
DC_OUT="${DC_DIR}/output_Datacard_leptons"
LIMDIR="${LIMDIR:-${FLASHGG}/Combine/merged_limits}"
PLOTDIR="${PLOTDIR:-${FLASHGG}/Plots/plot_limits}"

ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
ANACONDA_CONDA_SH="${ANACONDA_CONDA_SH:-/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh}"
PYROOT_ENV="${PYROOT_ENV:-higgs-alp-ana}"

# ---- stage switches ---------------------------------------------------------
RUN_MERGED_DATA_PRODUCTION="${RUN_MERGED_DATA_PRODUCTION:-0}"   # blocker; off by default
RUN_P2ROOT="${RUN_P2ROOT:-1}"
RUN_TREE2WS="${RUN_TREE2WS:-1}"
RUN_SIG_MODEL="${RUN_SIG_MODEL:-1}"
RUN_BKG="${RUN_BKG:-1}"
RUN_DATACARD="${RUN_DATACARD:-1}"
RUN_COMBINE="${RUN_COMBINE:-1}"
RUN_PLOT_LIMITS="${RUN_PLOT_LIMITS:-1}"

echo_step() { echo; echo "==== $* ===="; }
echo_skip() { echo_step "Skip: $*"; }

# ---- env helpers ------------------------------------------------------------
init_conda_shell() {
    if declare -F conda >/dev/null 2>&1; then return 0; fi
    local c
    for c in "$ANACONDA_SETUP" "$ANACONDA_CONDA_SH"; do
        [[ -r "$c" ]] || continue
        set +u; source "$c" >/dev/null 2>&1 || true; set -u
        if declare -F conda >/dev/null 2>&1; then return 0; fi
        if command -v conda >/dev/null 2>&1; then
            set +u; eval "$(conda shell.bash hook 2>/dev/null)"; set -u
            declare -F conda >/dev/null 2>&1 && return 0
        fi
    done
    return 1
}

activate_pyroot_env() {
    if [[ "${CONDA_DEFAULT_ENV:-}" == "$PYROOT_ENV" ]]; then return 0; fi
    init_conda_shell || { echo "[ERROR] cannot init conda" >&2; exit 1; }
    set +u; conda activate "$PYROOT_ENV"; set -u
    [[ "${CONDA_DEFAULT_ENV:-}" == "$PYROOT_ENV" ]] || {
        echo "[ERROR] failed to activate ${PYROOT_ENV}" >&2; exit 1; }
    echo "[ENV] conda env = ${CONDA_DEFAULT_ENV}"
}

deactivate_conda() {
    if declare -F conda >/dev/null 2>&1 && [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
        set +u; conda deactivate || true; set -u
    fi
}

activate_cmsenv() {
    [[ -n "${CMSSW_BASE:-}" && "${CMSSW_BASE}" == "${CMSSW_SRC%/src}" ]] && return 0
    deactivate_conda
    set +u
    source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
    pushd "$CMSSW_SRC" >/dev/null
    eval "$(scramv1 runtime -sh)"
    popd >/dev/null
    set -u
    export PYTHONPATH="${FLASHGG}/tools:${FLASHGG}/Signal/tools:${PYTHONPATH:-}"
    echo "[ENV] CMSSW_BASE=${CMSSW_BASE:-<unset>}"
}

# Run a list of "cmd args" jobs with a parallelism cap. Each arg is a single
# string passed to bash -c.
run_parallel() {
    local -a pids=()
    local job
    for job in "$@"; do
        bash -c "$job" &
        pids+=($!)
        if (( ${#pids[@]} >= MAX_PARALLEL )); then wait "${pids[@]}"; pids=(); fi
    done
    [[ ${#pids[@]} -gt 0 ]] && wait "${pids[@]}"
    return 0
}

# =============================================================================
#  Stage 0 -- produce 2024 merged-selection DATA friend parquet  (BLOCKER)
# =============================================================================
if [[ "$RUN_MERGED_DATA_PRODUCTION" == "1" ]]; then
    echo_step "merged: produce ${YEAR} DATA friend parquet"
    DATA_SCRIPT="${PROJECT_DIR}/HiggsDNA/scripts/run_merged_data_2024.sh"
    if [[ -x "$DATA_SCRIPT" ]]; then
        ( cd "${PROJECT_DIR}/HiggsDNA" && bash "$DATA_SCRIPT" )
    else
        echo "[WARN] ${DATA_SCRIPT} not found." >&2
    fi
else
    echo_skip "merged: produce DATA friend parquet (using existing parquet)"
fi

# =============================================================================
#  Stage 1 -- merged_p2root: friend parquet -> flashgg DiphotonTree ROOT
#             (conda higgs-alp-ana; uproot)
# =============================================================================
if [[ "$RUN_P2ROOT" == "1" ]]; then
    echo_step "merged: parquet -> flashgg ROOT trees (signal + data)"
    activate_pyroot_env
    sig_csv=$(echo "$MASSES" | sed -E 's/([0-9p]+)/M\1/g; s/ +/,/g')
    # per-mA ROI ML path: signal from Sig_MC_MLNANO_all, data per-mA (merged_p2root
    # uses its DATA_ML_GLOB default = parquet_friend/Data_2024_*/). Writes
    # sig/mA_M<m>/output_<year>.root and data/Data_merged_M<m>/run3.root (PER-mA).
    python3 "${MERGED_DIR}/merged_p2root.py" \
        --indir "${SIG_ML_DIR}" --outdir "${ROOT_MVACUT}" \
        --signals "$sig_csv" --year "$YEAR" --do-signal --do-data
else
    echo_skip "merged: parquet -> ROOT"
fi

# ----- everything below needs the flashgg cmsenv -----------------------------
if [[ "$RUN_TREE2WS$RUN_SIG_MODEL$RUN_BKG$RUN_DATACARD$RUN_COMBINE" != "00000" ]]; then
    echo_step "Switch to flashgg combine environment (cmsenv)"
    activate_cmsenv
fi

# =============================================================================
#  Stage 2 -- Tree2WS: signal (per lepton) + data
# =============================================================================
if [[ "$RUN_TREE2WS" == "1" ]]; then
    echo_step "merged: Tree2WS (signal per-lepton + data)"
    cd "${FLASHGG}/Trees2WS"
    mkdir -p logs
    jobs=()
    for m in $MASSES; do
        for lep in $LEPS; do
            jobs+=("cd '${FLASHGG}/Trees2WS' && python3 trees2ws.py --inputConfig config.py \
--inputTreeFile '${ROOT_MVACUT}/sig/mA_M${m}/output_${YEAR}.root' --inputMass 125 \
--productionMode ggh --year ${YEAR} --lepton ${lep} --doSystematics \
> logs/sig_mA_M${m}_${YEAR}_${lep}.log 2>&1")
        done
    done
    run_parallel "${jobs[@]}"
    # data: PER-mA now (each m_a has its own ROI sub-sample) ->
    #   data/Data_merged_M<m>/run3.root  ->  data/Data_merged_M<m>/ws/run3.root
    for m in $MASSES; do
        dtree="${ROOT_MVACUT}/data/mA_M${m}/run3.root"
        [[ -f "$dtree" ]] || { echo "  [WARN] missing per-mA data tree M${m}: $dtree"; continue; }
        python3 trees2ws_data.py --inputConfig config.py \
            --inputTreeFile "$dtree" > "logs/data_mA_M${m}.log" 2>&1
    done
    echo "  Tree2WS done (signal + per-mA data)"
else
    echo_skip "merged: Tree2WS"
fi

# =============================================================================
#  Stage 3 -- signal model: fTest -> calcPhotonSyst -> signalFit (DCB + 8 syst)
# =============================================================================
if [[ "$RUN_SIG_MODEL" == "1" ]]; then
    echo_step "merged: signal model (fTest -> calcPhotonSyst -> signalFit)"
    cd "${SIG_DIR}/scripts"
    mkdir -p "${SIG_DIR}/logs_merged"
    jobs=()
    for m in $MASSES; do
        for lep in $LEPS; do
            ws="${ROOT_MVACUT}/sig/mA_M${m}/ws_Tree2WS"
            lg="${SIG_DIR}/logs_merged/${m}_${lep}.log"
            jobs+=("cd '${SIG_DIR}/scripts' && \
python3 fTest.py          --mass_ALP ${m} --year ${YEAR} --mass 125 --channel ${lep} --inputWSDir '${ws}'  > '${lg}' 2>&1 && \
python3 calcPhotonSyst.py --mass_ALP ${m} --year ${YEAR}            --channel ${lep} --inputWSDir '${ws}' >> '${lg}' 2>&1 && \
python3 signalFit.py      --mass_ALP ${m} --year ${YEAR}            --channel ${lep} --inputWSDir '${ws}' --doSystematics >> '${lg}' 2>&1")
        done
    done
    run_parallel "${jobs[@]}"
    n=0; for m in $MASSES; do for lep in $LEPS; do
        [[ -f "${SIG_DIR}/outdir_${lep}/signalFit/output/${m}_CMS-HGG_sigfit_${YEAR}_${lep}_Hm125.root" ]] && n=$((n+1)) || echo "  [WARN] missing signal model M${m} ${lep}"
    done; done
    echo "  signal models built: ${n}"
else
    echo_skip "merged: signal model"
fi

# =============================================================================
#  Stage 4 -- background: PER-mA discrete-profiling fit on each m_a data ROI WS
# =============================================================================
if [[ "$RUN_BKG" == "1" ]]; then
    echo_step "merged: background discrete-profiling envelope (PER-mA ROI)"
    cd "${FLASHGG}/Background"
    if [[ ! -x ./bin/fTest_ALP_turnOn ]]; then
        echo "[ERROR] ./bin/fTest_ALP_turnOn missing (compile Background first)" >&2; exit 1
    fi
    # PER-mA: each m_a selects a different MLPhoton_lead_mass ROI sub-sample -> its
    # own data spectrum -> its own multipdf. Output dir = ${BKG_WSDIR}/<m> so it
    # matches makeYields bkg_mass_for_io=<mass_ALP> (Stage 5). Unlike the old shared
    # fit, the per-mA bkgplots are genuinely different (different data sub-samples).
    n_bkg=0
    for m in $MASSES; do
        data_ws="${ROOT_MVACUT}/data/mA_M${m}/ws/run3.root"
        if [[ ! -f "$data_ws" ]]; then
            echo "  [WARN] M${m} data WS not found: ${data_ws}; skip"; continue
        fi
        bkg_m="${BKG_WSDIR}/${m}"
        mkdir -p "${bkg_m}/HZAmassInde_fTest"
        ./bin/fTest_ALP_turnOn -i "${data_ws}" \
            --saveMultiPdf "${bkg_m}/CMS-HGG_mva_13p6TeV_multipdf.root" \
            -D "${bkg_m}/HZAmassInde_fTest" \
            --mass_ALP 1 -c 1 --isFlashgg 0 --isData 0 -f data, \
            --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
            --gtoys "$GOF_TOYS" > "${bkg_m}/ftest.log" 2>&1
        if [[ ! -f "${bkg_m}/CMS-HGG_mva_13p6TeV_multipdf.root" ]]; then
            echo "[ERROR] M${m} bkg multipdf not produced; see ${bkg_m}/ftest.log" >&2; continue
        fi
        echo "  M${m} bkg multipdf: ${bkg_m}/CMS-HGG_mva_13p6TeV_multipdf.root"
        n_bkg=$((n_bkg+1))
        # per-mA banded Data+envelope plot (real per-mA, NOT a shared symlink)
        if [[ -x ./bin/makeBkgPlots_ALP ]]; then
            AFR="${bkg_m}/AllFitResults"; mkdir -p "$AFR"
            ma_gev="${m//p/.}"   # flashgg label 0p5 -> physical mA 0.5 GeV for the plot label
            ./bin/makeBkgPlots_ALP -b "${bkg_m}/CMS-HGG_mva_13p6TeV_multipdf.root" \
                -d "${bkg_m}/BkgPlots" --total_OutDir "$AFR" -o "${bkg_m}/BkgPlots.root" \
                --sqrts 13p6TeV --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal "${ma_gev}" \
                --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
                --intLumi "$RUN3_LUMI" -c 0 --isFlashgg 0 --doBands > "${bkg_m}/bkgplot.log" 2>&1 || true
            # makeBkgPlots names outputs bkgplot_<%.0f>.pdf (0 or 1 for sub-GeV); glob picks the single file per per-mA dir
            for base in bkgplot allPdfs; do
                f=$(ls "$AFR"/${base}_[0-9]*.pdf 2>/dev/null | head -1)
                [[ -n "$f" ]] && mv -f "$f" "$AFR/${base}_M${m}.pdf"
            done
        fi
    done
    echo "  per-mA bkg fits done: ${n_bkg} under ${BKG_WSDIR}/<m>/"
else
    echo_skip "merged: background fit"
fi

# =============================================================================
#  Stage 5 -- datacard: makeYields(--mergedLowMA --lumi) -> makeDatacard
#                       -> text2workspace
# =============================================================================
if [[ "$RUN_DATACARD" == "1" ]]; then
    echo_step "merged: datacard (Run-3 lumi=${RUN3_LUMI}) per m_a"
    export PYTHONPATH="${FLASHGG}/tools:${PYTHONPATH:-}"
    lumi_args=(); [[ "$RUN3_LUMI" != "0" ]] && lumi_args=(--lumi "$RUN3_LUMI")
    mkdir -p "$LIMDIR"
    for m in $MASSES; do
        echo "  -- M${m}"
        cd "$DC_DIR"
        python3 makeYields.py --inputWSDirMap "${YEAR}=${ROOT_MVACUT}" \
            --mass_ALP ${m} --channel leptons --mergedLowMA "${lumi_args[@]}" \
            --bkgModelWSDir "${BKG_WSDIR}" \
            --doSystematics --ignore-warnings > "${LIMDIR}/makeYields_${m}.log" 2>&1
        python3 makeDatacard.py --mass_ALP ${m} --channel leptons --years "${YEAR}" \
            --doSystematics --output "mergedDatacard_${m}" > "${LIMDIR}/makeDatacard_${m}.log" 2>&1
        dc="${DC_OUT}/${m}_pruned_datacard_leptons.txt"
        [[ -f "$dc" ]] || { echo "  [WARN] M${m} no datacard; skip"; continue; }
        cd "$DC_OUT"
        text2workspace.py "${m}_pruned_datacard_leptons.txt" -o "ws_${m}_merged.root" -m 125 \
            > "${LIMDIR}/t2w_${m}.log" 2>&1
    done
else
    echo_skip "merged: datacard"
fi

# =============================================================================
#  Stage 6 -- combine: AsymptoticLimits per m_a (MH frozen, rMax 1000)
# =============================================================================
if [[ "$RUN_COMBINE" == "1" ]]; then
    echo_step "merged: combine AsymptoticLimits per m_a"
    mkdir -p "$LIMDIR"
    cd "$DC_OUT"
    for m in $MASSES; do
        ws="ws_${m}_merged.root"
        [[ -f "$ws" ]] || { echo "  [WARN] M${m} no workspace; skip"; continue; }
        combine -M AsymptoticLimits "$ws" -m 125 -n "_merged_flashgg_M${m}" \
            --freezeParameters MH --rMin 0 --rMax 1000 --rRelAcc 0.01 \
            > "${LIMDIR}/combine_${m}.log" 2>&1
        mv "higgsCombine_merged_flashgg_M${m}.AsymptoticLimits.mH125.root" "${LIMDIR}/" 2>/dev/null || true
        med=$(grep "Expected 50" "${LIMDIR}/combine_${m}.log" | sed 's/.*< //')
        echo "  M${m}  median expected r < ${med:-?}"
    done
else
    echo_skip "merged: combine limits"
fi

# =============================================================================
#  Stage 7 -- plot limits (conda higgs-alp-ana; log-x Brazil band)
# =============================================================================
if [[ "$RUN_PLOT_LIMITS" == "1" ]]; then
    echo_step "merged: plot expected limit vs m_a (log-x)"
    if command -v scram >/dev/null 2>&1 && [[ -n "${CMSSW_BASE:-}" ]]; then
        eval "$(scram unsetenv -sh)" || true
    fi
    activate_pyroot_env
    masses_csv=$(echo "$MASSES" | sed -E 's/([0-9])p([0-9])/0.\2/g; s/ +/,/g')
    plot_lumi="$RUN3_LUMI"; [[ "$plot_lumi" == "0" ]] && plot_lumi=109.82
    cd "${FLASHGG}/Plots"
    python3 makeLimitsPlot_lowmA.py \
        --combine-dir "${LIMDIR}" \
        --pattern "higgsCombine_merged_flashgg_M0p{i}.AsymptoticLimits.mH125.root" \
        --masses "$masses_csv" --assume-xs "$ASSUME_XS" \
        --lumi "$plot_lumi" --tag lowmA_flashgg_run3 --outdir "${PLOTDIR}"
    # Full-range overlay (merged sub-GeV + resolved >=1 GeV) on one log-x canvas.
    if [[ -n "$RESOLVED_LIMITS_DIR" && -d "$RESOLVED_LIMITS_DIR" ]]; then
        python3 makeLimitsPlot_full.py \
            --merged-dir "${LIMDIR}" --resolved-dir "${RESOLVED_LIMITS_DIR}" \
            --merged-masses "$masses_csv" --assume-xs "$ASSUME_XS" \
            --lumi "$plot_lumi" --tag full_0p1_30 --outdir "${PLOTDIR}"
    else
        echo "  [INFO] RESOLVED_LIMITS_DIR unset/missing -> skip full-range overlay"
    fi
else
    echo_skip "merged: plot limits"
fi

echo_step "merged pipeline done"
echo "  signal models : ${SIG_DIR}/outdir_<lep>/signalFit/output/"
echo "  bkg multipdf  : ${BKG_WSDIR}/<m>/CMS-HGG_mva_13p6TeV_multipdf.root (per-mA ROI)"
echo "  datacards     : ${DC_OUT}/<m>_pruned_datacard_leptons.txt"
echo "  limits        : ${LIMDIR}/higgsCombine_merged_flashgg_M0p*.root"
echo "  limit plot    : ${PLOTDIR}/Limits_XS_lowmA_flashgg_run3.{pdf,png}"
