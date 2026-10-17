#!/usr/bin/env bash
# Detached watchdog for the HZa_merged lepton-cleaning limit chain and the v5
# MLPhoton regressor trainings. Lands results in the registry; writes numbers
# only, no interpretation.
#   chain   : bash 1_grand_merged.sh (pid CHAIN_PID), done marker "merged pipeline done"
#   training: condor clusters, done marker "[train] best model" in <run>/train.log
set -uo pipefail
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs
REG=$LOGD/PENDING_RESULTS.txt
STATUS=$LOGD/watch_lepclean_v5.status
CHAIN_PID=871964
CHAIN_LOG=$LOGD/grand_merged_lepclean_20261013_013135.log
RUNS=/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_runs
declare -A CL=([v5a1]=12907724 [v5a2]=12907725 [v5a3]=12907726)
declare -A FINAL=()
chain_final=""

for _ in $(seq 1 576); do                       # 48 h ceiling at 5-min polls
    line="poll $(date '+%F %H:%M:%S')"
    # ---- chain
    if [ -z "$chain_final" ]; then
        if grep -q 'merged pipeline done' "$CHAIN_LOG" 2>/dev/null; then
            chain_final=DONE
            echo "$(date '+%F %H:%M') DONE  merged_limit_chain_lepclean  -> $CHAIN_LOG" >> "$REG"
        elif ! kill -0 "$CHAIN_PID" 2>/dev/null; then
            chain_final=DIED
            echo "$(date '+%F %H:%M') DIED  merged_limit_chain_lepclean  last: $(grep '====' "$CHAIN_LOG" | tail -1)" >> "$REG"
        fi
    fi
    line+=" chain=${chain_final:-running($(grep -c 'bkg multipdf:' "$CHAIN_LOG" 2>/dev/null)/9 bkg)}"
    # ---- trainings
    for t in v5a1 v5a2 v5a3; do
        if [ -z "${FINAL[$t]:-}" ]; then
            st=$(condor_q "${CL[$t]}" -af JobStatus 2>/dev/null | head -1)
            if [ -z "$st" ]; then                 # left the queue: judge by the artifact
                if grep -q '^\[train\] best model' "$RUNS/regressor_$t/train.log" 2>/dev/null; then
                    FINAL[$t]=DONE
                    echo "$(date '+%F %H:%M') DONE  regressor_$t  -> $RUNS/regressor_$t/train.log" >> "$REG"
                else
                    FINAL[$t]=DIED
                    echo "$(date '+%F %H:%M') DIED  regressor_$t  cluster ${CL[$t]} left queue without best-model line" >> "$REG"
                fi
            elif [ "$st" = 5 ]; then
                line+=" $t=HELD"
                continue
            fi
        fi
        line+=" $t=${FINAL[$t]:-status$st}"
    done
    echo "$line" > "$STATUS"
    n=0; for t in v5a1 v5a2 v5a3; do [ -n "${FINAL[$t]:-}" ] && n=$((n+1)); done
    [ -n "$chain_final" ] && [ "$n" = 3 ] && { echo "$line ALL_FINAL" > "$STATUS"; exit 0; }
    sleep 300
done
echo "$(date '+%F %H:%M') TIMEOUT watch_lepclean_v5 after 48h" >> "$REG"
