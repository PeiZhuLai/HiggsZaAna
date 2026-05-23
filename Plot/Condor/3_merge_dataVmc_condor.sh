#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PLOT_DIR="${PROJECT_DIR}/Plot"
SCRIPTS_DIR="${PLOT_DIR}/scripts"
OUTPUT_DIR="${PLOT_DIR}/plots"
VARIABLES_DIR="${OUTPUT_DIR}/variables_dataVmc"
LOG_DIR="${OUTPUT_DIR}/logs_split"
JOBS_FILE="${PROJECT_DIR}/Plot/Condor/dataVmc_jobs.txt"
SUBMIT_FILE="${PROJECT_DIR}/Plot/Condor/dataVmc.submit"
JOBS_GENERATOR="${PROJECT_DIR}/Plot/Condor/make_dataVmc_condor_jobs.py"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RUN_MERGE_PLOTS="${RUN_MERGE_PLOTS:-1}"
RUN_DATAVMC_PLOTS="${RUN_DATAVMC_PLOTS:-1}"
RUN_OPTIMIZATION="${RUN_OPTIMIZATION:-1}"
CHECK_ROOT_KEYS="${CHECK_ROOT_KEYS:-1}"
CHECK_SYSTEMATICS="${CHECK_SYSTEMATICS:-1}"
REMAKE_JOBS="${REMAKE_JOBS:-0}"

mkdir -p "$LOG_DIR" "$VARIABLES_DIR"
export PYTHONPATH="${PYTHONPATH:-}:${PLOT_DIR}/lib:${PROJECT_DIR}/HZaMVA/scripts"

cd "$PLOT_DIR"

echo "[Config] RUN_MERGE_PLOTS=$RUN_MERGE_PLOTS"
echo "[Config] RUN_DATAVMC_PLOTS=$RUN_DATAVMC_PLOTS"
echo "[Config] RUN_OPTIMIZATION=$RUN_OPTIMIZATION"
echo "[Config] CHECK_SYSTEMATICS=$CHECK_SYSTEMATICS"
echo "[Config] VARIABLES_DIR=$VARIABLES_DIR"
echo "[Config] LOG_DIR=$LOG_DIR"
echo "[Config] JOBS_FILE=$JOBS_FILE"

if [[ "$RUN_DATAVMC_PLOTS" == "1" && "$RUN_MERGE_PLOTS" != "1" ]]; then
    echo "[Info] RUN_DATAVMC_PLOTS=1 with RUN_MERGE_PLOTS=$RUN_MERGE_PLOTS; redraw uses existing merged ROOT files in $VARIABLES_DIR"
fi

jobs_stale=0
if [[ -s "$JOBS_FILE" && "$JOBS_FILE" -ot "$JOBS_GENERATOR" ]]; then
    jobs_stale=1
fi
if [[ -s "$SUBMIT_FILE" && "$SUBMIT_FILE" -ot "$JOBS_GENERATOR" ]]; then
    jobs_stale=1
fi

if [[ "$REMAKE_JOBS" == "1" || "$jobs_stale" == "1" || ! -s "$JOBS_FILE" || ! -s "$SUBMIT_FILE" ]]; then
    if [[ "$jobs_stale" == "1" && "$REMAKE_JOBS" != "1" ]]; then
        echo "[Info] Existing submit/jobs are older than $JOBS_GENERATOR; regenerating."
    fi
    make_jobs_args=()
    if [[ -n "${DATA_VMC_MAKE_JOBS_ARGS:-}" ]]; then
        # shellcheck disable=SC2206
        make_jobs_args=(${DATA_VMC_MAKE_JOBS_ARGS})
    fi
    "$PYTHON_BIN" "$JOBS_GENERATOR" --condor-dir "${PROJECT_DIR}/Plot/Condor" "${make_jobs_args[@]}"
fi

region_suffix() {
    local region_key="$1"
    case "$region_key" in
        SR)  echo "UL_SR" ;;
        CR)  echo "UL_CR" ;;
        mva) echo "UL_mva" ;;
        *)
            echo "[ERROR] Unknown region key: $region_key" >&2
            return 1
            ;;
    esac
}

root_has_keys() {
    local path="$1"

    if [[ "$CHECK_ROOT_KEYS" != "1" ]]; then
        return 0
    fi

    "$PYTHON_BIN" - "$path" <<'PY' >/dev/null 2>&1
import sys

try:
    import ROOT
except Exception:
    sys.exit(2)

path = sys.argv[1]
ROOT.gROOT.SetBatch(True)
root_file = ROOT.TFile.Open(path, "READ")
if not root_file or root_file.IsZombie():
    sys.exit(1)

keys = root_file.GetListOfKeys()
sys.exit(0 if keys and keys.GetEntries() > 0 else 1)
PY
}

root_has_nontrivial_systematics() {
    local path="$1"

    if [[ "$CHECK_SYSTEMATICS" != "1" ]]; then
        return 0
    fi

    "$PYTHON_BIN" - "$path" <<'PY'
import sys

try:
    import ROOT
except Exception:
    sys.exit(2)

path = sys.argv[1]
samples = ("DYJetsToLL", "DYGto2LG")

ROOT.gROOT.SetBatch(True)
root_file = ROOT.TFile.Open(path, "READ")
if not root_file or root_file.IsZombie():
    sys.exit(1)

raw_dir = root_file.Get("raw_plots")
sys_dir = root_file.Get("sys_dir")
if not raw_dir or not sys_dir:
    print("[ERROR] Missing raw_plots or sys_dir in %s" % path)
    sys.exit(1)

checks = 0
different = 0
examples = []

for key in sys_dir.GetListOfKeys():
    sys_hist_name = key.GetName()
    sample = None
    marker = None
    for candidate in samples:
        candidate_marker = "_" + candidate + "_"
        if candidate_marker in sys_hist_name:
            sample = candidate
            marker = candidate_marker
            break
    if sample is None:
        continue

    var_name, syst_name = sys_hist_name.rsplit(marker, 1)
    if not (syst_name.endswith("_up") or syst_name.endswith("_down")):
        continue

    nominal = raw_dir.Get(var_name + "_" + sample)
    shifted = sys_dir.Get(sys_hist_name)
    if not nominal or not shifted:
        continue

    checks += 1
    diff = 0.0
    for i_bin in range(0, nominal.GetNbinsX() + 2):
        diff += abs(float(shifted.GetBinContent(i_bin)) - float(nominal.GetBinContent(i_bin)))

    scale = max(abs(float(nominal.Integral())), 1.0)
    if diff > max(1e-9, 1e-8 * scale):
        different += 1
        if len(examples) < 3:
            examples.append("%s diff=%.6g" % (sys_hist_name, diff))

root_file.Close()

if checks == 0:
    print("[ERROR] No background systematic histograms found in %s" % path)
    sys.exit(1)

if different == 0:
    print("[ERROR] All %d checked background systematic histograms match nominal in %s" % (checks, path))
    print("[ERROR] This usually means the partial ROOT files were produced with --skipSystematics or from stale outputs.")
    sys.exit(1)

print("[OK] Systematics check: %d/%d background systematic histograms differ from nominal" % (different, checks))
for example in examples:
    print("[OK]   " + example)
PY
}

collect_sample_tags() {
    local region_key="$1"
    local final_tag="$2"

    if [[ ! -s "$JOBS_FILE" ]]; then
        echo "[ERROR] Missing jobs file: $JOBS_FILE" >&2
        return 1
    fi

    awk -v region_key="$region_key" -v final_tag="$final_tag" '
        $1 == region_key && $2 == final_tag { print $3 }
    ' "$JOBS_FILE"
}

merge_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local suffix target input sample_tag
    local inputs=()
    local sample_tags=()

    suffix="$(region_suffix "$region_key")"
    target="${VARIABLES_DIR}/ALP_plot_run3_${suffix}_${final_tag}.root"

    while IFS= read -r sample_tag; do
        [[ -n "$sample_tag" ]] || continue
        sample_tags+=("$sample_tag")
    done < <(collect_sample_tags "$region_key" "$final_tag")

    if [[ "${#sample_tags[@]}" -eq 0 ]]; then
        echo "[ERROR] No sample tags found in $JOBS_FILE for tag=${final_tag} region=${region_key}" >&2
        return 1
    fi

    for sample_tag in "${sample_tags[@]}"; do
        input="${VARIABLES_DIR}/ALP_plot_run3_${suffix}_${final_tag}_part_${sample_tag}.root"
        if [[ ! -s "$input" ]]; then
            echo "[ERROR] Missing partial ROOT file: $input" >&2
            return 1
        fi
        if ! root_has_keys "$input"; then
            echo "[ERROR] Partial ROOT file has no keys or is unreadable: $input" >&2
            return 1
        fi
        inputs+=("$input")
    done

    echo "[hadd] $target"
    hadd -f "$target" "${inputs[@]}"
    if ! root_has_keys "$target"; then
        echo "[ERROR] Merged ROOT file has no keys or is unreadable: $target" >&2
        return 1
    fi
    if ! root_has_nontrivial_systematics "$target"; then
        echo "[ERROR] Merged ROOT file has no non-trivial background systematics: $target" >&2
        echo "[ERROR] Re-run Plot/Condor/1_submit_dataVmc_condor.sh without SKIP_SYSTEMATICS=1, then merge again." >&2
        return 1
    fi
}

draw_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local log_file="${LOG_DIR}/draw_${final_tag}_${region_key}.log"
    local cmd=(
        "$PYTHON_BIN" "$SCRIPTS_DIR/2_plot_dataVmc.py"
        -y run3
        -m
        --ln
        --inputTag "$final_tag"
        --outputTag "$final_tag"
    )

    case "$region_key" in
        SR)  cmd+=(--region 1) ;;
        CR)  cmd+=(--region 2) ;;
        mva) cmd+=(-b) ;;
        *)
            echo "[ERROR] Unknown region key: $region_key" >&2
            return 1
            ;;
    esac

    {
        echo "[START] $(date '+%F %T') draw tag=${final_tag} region=${region_key}"
        printf "[CMD]"
        printf " %q" "${cmd[@]}"
        echo
    } > "$log_file"

    if "${cmd[@]}" >> "$log_file" 2>&1; then
        echo "[DONE] $(date '+%F %T')" >> "$log_file"
    else
        local status=$?
        echo "[FAILED] $(date '+%F %T') exit=${status}" >> "$log_file"
        return "$status"
    fi
}

if [[ "$RUN_MERGE_PLOTS" == "1" ]]; then
    merge_pids=()
    merge_labels=()

    for final_tag in nominal sideband_rwgt; do
        for region_key in SR CR mva; do
            echo "[submit hadd] tag=${final_tag} region=${region_key}"
            merge_plot_output "$region_key" "$final_tag" &
            merge_pids+=("$!")
            merge_labels+=("${final_tag}/${region_key}")
        done
    done

    merge_status=0
    for i in "${!merge_pids[@]}"; do
        if wait "${merge_pids[$i]}"; then
            echo "[done hadd] ${merge_labels[$i]}"
        else
            echo "[ERROR] hadd failed: ${merge_labels[$i]}" >&2
            merge_status=1
        fi
    done

    if [[ "$merge_status" != "0" ]]; then
        exit "$merge_status"
    fi
else
    echo "[Info] RUN_MERGE_PLOTS=$RUN_MERGE_PLOTS; skip hadd merge"
fi

if [[ "$RUN_DATAVMC_PLOTS" == "1" ]]; then
    for final_tag in nominal sideband_rwgt; do
        for region_key in SR CR mva; do
            draw_plot_output "$region_key" "$final_tag"
        done
    done
else
    echo "[Info] RUN_DATAVMC_PLOTS=$RUN_DATAVMC_PLOTS; skip 2_plot_dataVmc.py"
fi

if [[ "$RUN_OPTIMIZATION" == "1" ]]; then
    for final_tag in nominal sideband_rwgt; do
        "$PYTHON_BIN" "$SCRIPTS_DIR/ALP_Optimization.py" \
            -y run3 \
            -o "${OUTPUT_DIR}/optimize_run3UL_${final_tag}" \
            --region 1 \
            -p --sigVSscore -s --doOpt -c 1 \
            --inputTag "$final_tag"
    done
else
    echo "[Info] RUN_OPTIMIZATION=$RUN_OPTIMIZATION; skip ALP_Optimization.py"
fi
