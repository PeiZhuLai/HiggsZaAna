#!/usr/bin/env bash
#
# Submit the bkg and data re-production in waves, with the concurrency bounded.
#
# WHY
#   The signal stage was submitted as one block of 1104 jobs. The farm started
#   ~1089 of them at once, they all read the hza_ana conda env off /eos at the
#   same moment, and EOS gave up:
#       749x  OSError: [Errno 5] Input/output error
#        98x  ImportError: /eos/.../hza_ana/.../*.so: cannot read file data
#   66 % of the 2304 job attempts failed; only HiggsDNA's retry loop got the
#   stage to finish. Measured later in the same run, once the farm had drained
#   to ~211 concurrent jobs: 68 attempts, ZERO I/O errors. The failure rate is a
#   function of concurrency, nothing else.
#
# STRATEGY (two knobs, no new machinery)
#   1. Raise fpo. The conda-env read is a fixed cost per JOB, so packing 4x more
#      input files into each job cuts the number of EOS env reads by 4x for the
#      same CPU work. bkg 3287 -> ~822 jobs, data 3475 -> ~434 jobs.
#   2. Submit one (stage, year) wave at a time and wait for it. run_ana_*.sh
#      blocks until its jobs are done, and each wave is sized to land near the
#      ~250-300 concurrent jobs that were measured to be clean.
#   Waves also mean a bad wave costs one wave, not the whole stage.
#
#   If a wave still shows many "Input/output error" failures, the fallback is
#   B-mode (HZA_BMODE_PACK, conda-pack tarball to node-local scratch) -- the
#   wrapper already supports it, it just needs a tarball built.
#
# USAGE
#   bash submit_bkg_data_waves.sh             # dry run, prints the plan
#   DRY_RUN=0 bash submit_bkg_data_waves.sh   # submit
#   DRY_RUN=0 STAGES=bkg bash submit_bkg_data_waves.sh
#
set -uo pipefail

D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
REPO="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
STAGING="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix"
LOGDIR="${D}/logs_fsrfix"
DRY_RUN="${DRY_RUN:-1}"
STAGES="${STAGES:-bkg data}"
# Space-separated "stage:year" entries to skip, for waves already done outside
# this script (e.g. the manual verification wave).
SKIP="${SKIP:-}"

ENVDIR="/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana"
export CONDA_PREFIX="${ENVDIR}"
export PATH="${ENVDIR}/bin:${PATH}"
export PYTHONPATH="${REPO}"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"

# wave := "<stage> <year> <fpo> <expected jobs> <sample list>"
#
# The sample list is PER YEAR and must contain only samples that actually exist
# that year. Passing the full 7-sample list to a single year makes HiggsDNA
# build a Task with zero jobs for the absent ones, and Task.summarize() then
# divides by counts["all"] == 0:
#     progress_bar.py:18  ZeroDivisionError: division by zero
# That killed every bkg wave on the first attempt.
#
# expected jobs = (job dirs of the previous fpo=4 production) * 4 / fpo
# 2022 used to be covered by two PTG-binned DY+gamma samples, but CMS has
# invalidated every 2022 version of both (PTG-10to50 and PTG-50to100, v11/v12/v13,
# preEE and postEE alike) -- DAS hides INVALID datasets, so the sample manager
# just reported "No files with instance=prod/global" and built a Task with zero
# jobs, which crashed the driver on ZeroDivisionError. The single-range
# PTG-10to100 sample is VALID for 2022 (341 + 872 files) and its cross section
# is already the sum of the two (124 + 2.088 = 126.1 ~ 126.6 pb), so 2022 now
# uses the same binning as 2023. The 2022 paths were added to
# metadata/samples/zgamma_tutorial.json.
B22="DYGto2LG_10to100,DYJetsToLL"
B23="DYGto2LG_10to100,DYJetsToLL"
B24="DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau"
WAVES=(
  "bkg  2022preEE      4  166 ${B22}"
  "bkg  2022postEE     4  492 ${B22}"
  "bkg  2023preBPix    4  187 ${B23}"
  "bkg  2023postBPix   4  108 ${B23}"
  "bkg  2024           4  382 DYGto2LG_10to100"
  "bkg  2024           4  748 DYJetsTo2E"
  "bkg  2024           4  744 DYJetsTo2Mu"
  "bkg  2024           4  460 DYJetsTo2Tau"
  "data 2022preEE      4  169 Data"
  "data 2022postEE     4  340 Data"
  "data 2023preBPix    4  294 Data"
  "data 2023postBPix   4  145 Data"
  "data 2024           4 2527 Data"
)

# fpo stays at 4 -- the value the original production ran at.
#
# Raising it was a mistake. The idea was to cut the number of per-job conda-env
# reads off /eos, but what actually fixes the EOS storm is bounding CONCURRENCY,
# which the wave structure already does. Meanwhile HiggsDNA accumulates all of a
# job's input files in memory before writing its single parquet, so memory scales
# with fpo and the tail explodes:
#     fpo=1  (signal, n=397): median 35 MB, p90 74 MB, max 4151 MB
#     fpo=16 (measured here): 23060-23699 MB  -> held against a 12000 MB request
# fpo=32 would have needed ~46 GB/job, which essentially never matches a slot.
# So: keep memory at the proven 3000 MB and control concurrency with waves only.
#
# 2024 is split by SAMPLE instead of by fpo, which keeps every wave under ~750
# jobs without touching memory at all. data 2024 is the one sample that cannot be
# split that way (2527 jobs); it is throttled at the condor level instead.
WAVE_MEMORY="${WAVE_MEMORY:-3000}"

# Each wave rebuilds analysis_manager.pkl from scratch (CLEAN_ANALYSIS_STATE=1).
# Keeping the pickle across waves looks tidier but is wrong: AnalysisManager
# takes the PICKLED sample_list/fpo over the ones passed on the command line --
#   "kwarg 'sample_list' ... will be ignored in favor of the pickled value"
# -- so a stale pickle silently reimposes the full 7-sample list and its stale
# (empty) file lists, and every Task is built with 0 input files:
#   Task 'DYGto2LG_10to50_2022preEE' : splitting 0 input files into 0 jobs
#   ZeroDivisionError: division by zero
# Rebuilding is safe because each wave is restricted to its own SAMPLE_LIST and
# YEARS, the waves run serially, and jobs whose output already exists are skipped.

script_for() { [ "$1" = "bkg" ] && echo run_ana_bkgmc.sh || echo run_ana_data.sh; }
outdir_for() { [ "$1" = "bkg" ] && echo "${STAGING}/Bkg_MC" || echo "${STAGING}/Data"; }

echo "=== bkg + data re-production, wave schedule ==="
printf '%-6s %-14s %5s %8s\n' stage year fpo jobs
total=0
for w in "${WAVES[@]}"; do
    read -r st yr fpo nj sl <<<"$w"
    case " ${STAGES} " in *" ${st} "*) ;; *) continue ;; esac
    printf '%-6s %-14s %5s %8s\n' "$st" "$yr" "$fpo" "$nj"
    total=$((total + nj))
done
echo "total jobs across the selected waves: ${total}  (was 6762 at fpo=4)"
echo "dry run: ${DRY_RUN}"
echo

mkdir -p "${LOGDIR}"
for w in "${WAVES[@]}"; do
    read -r st yr fpo nj sl <<<"$w"
    case " ${STAGES} " in *" ${st} "*) ;; *) continue ;; esac

    case " ${SKIP} " in *" ${st}:${yr} "*) echo "[wave] skip ${st} ${yr} (already done)"; continue ;; esac

    # Never start a wave while another driver still owns this output directory:
    # two run_analysis processes on the same outdir race on analysis_manager.pkl.
    if [ "${DRY_RUN}" != "1" ]; then
        waited=0
        while ps -u "$(id -un)" -o cmd= 2>/dev/null | grep -q "[r]un_analysis.py --config metadata/za_$( [ "$st" = bkg ] && echo bkgmc || echo data )_run3.json"; do
            [ "$waited" -eq 0 ] && echo "[wave] waiting for the previous ${st} driver to finish"
            waited=$((waited+1)); sleep 120
            [ "$waited" -gt 180 ] && { echo "[wave] ABORT: previous ${st} driver still running after 6 h"; exit 1; }
        done
    fi

    left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
    if [ "${left:-0}" -lt 7200 ] && [ "${DRY_RUN}" != "1" ]; then
        echo "[wave] ABORT before ${st}/${yr}: proxy under 2 h (${left}s)."
        echo "       voms-proxy-init --rfc --voms cms -valid 192:00, then rerun."
        exit 1
    fi

    log="${LOGDIR}/wave_${st}_${yr}.log"
    echo "[wave] ${st} ${yr} (fpo=${fpo}, ~${nj} jobs) -> ${log}"
    if [ "${DRY_RUN}" = "1" ]; then
        echo "       would run: OUTDIR=$(outdir_for "$st") YEARS=${yr} FPO=${fpo} SAMPLE_LIST=${sl} CONDOR_REQ_MEMORY=${WAVE_MEMORY}"
        continue
    fi

    ( cd "${REPO}" && env \
        "CONDA_PREFIX=${ENVDIR}" "PATH=${ENVDIR}/bin:${PATH}" "PYTHONPATH=${REPO}" \
        "OUTDIR=$(outdir_for "$st")" "YEARS=${yr}" "FPO=${fpo}" \
        "SAMPLE_LIST=${sl}" "CONDOR_REQ_MEMORY=${WAVE_MEMORY}" \
        "CLEAN_ANALYSIS_STATE=1" \
        bash "scripts/$(script_for "$st")" ) > "${log}" 2>&1
    rc=$?
    echo "[wave] ${st} ${yr} driver returned ${rc} at $(date '+%F %T')"

    # The driver's exit status means nothing -- it counts retired jobs as
    # complete and returned 0 for a stage where every job had died. Judge the
    # wave by what landed on disk, and window the I/O error count to this wave
    # (a cumulative grep over eos_logs just reports the signal-stage history
    # forever and never changes).
    got=$(find "$(outdir_for "$st")" -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' \
          -newermt "-${WAVE_WINDOW_MIN:-600} minutes" 2>/dev/null | wc -l | tr -d ' ')
    io=$(find "${REPO}/eos_logs" -name '*.err' -newermt "-${WAVE_WINDOW_MIN:-600} minutes" 2>/dev/null \
         | xargs -r grep -l "Input/output error" 2>/dev/null | wc -l | tr -d ' ')
    echo "[wave] ${st} ${yr}: chunks written in this window=${got} (expected ~${nj}), EOS I/O failures=${io}"
    if [ "${got:-0}" -lt $(( nj / 2 )) ]; then
        echo "[wave] ALERT ${st} ${yr} produced ${got} of ~${nj} expected chunks -- stopping the wave train"
        echo "[wave] inspect ${log} and ${REPO}/eos_logs before continuing"
        exit 1
    fi
done
echo "=== waves finished (dry run=${DRY_RUN}) ==="
