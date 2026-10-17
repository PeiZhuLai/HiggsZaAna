#!/usr/bin/env bash
# CERN-side loop for the IHEP part of the HZa_merged lepclean re-production.
# Runs as systemd --user unit hza-lc-ihep (IHEP login nodes kill long local loops).
#
# Every POLL s: ssh one cycle of _ihep_lc_cycle.sh (held memory ladder, missing-idx
# resubmission, status) and append its output to ihep_cycle.log. STOP lines (held at
# 20 GB, idx failed 10 times) go to the registry once.
# When all EXPECT_TAGS tags are complete and nothing is queued:
#   1. submit _ihep_lc_merge_job.sh once (hep_sub, 8000 MB) and wait for MERGE_DONE;
#   2. rc=0 -> rsync every <tag>/merged/merged_nominal.parquet back to CERN
#      parquet_friend_lepclean/<tag>/merged/, re-check each file's row count at CERN
#      against the IHEP report, registry DONE; otherwise registry FAIL + report.
set -uo pipefail
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
LOG=$LOGD/ihep_cycle.log
IR=/publicfs/cms/user/laipeizhu/HZa/HiggsZaAna/HiggsDNA
IOUT=/publicfs/cms/user/laipeizhu/Eos/HZa_merged/parquet_friend_lepclean
COUT=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_lepclean
POLL="${POLL:-600}"   # throttled submission needs a tighter poll than the 30-min report
EXPECT_TAGS="${EXPECT_TAGS:-36}"
export PATH="/usr/bin:/bin"
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

issh() {  # try the login nodes in turn
    for h in ihep1 ihep2 ihep3 ihep4; do
        out=$(timeout 1800 ssh -o BatchMode=yes -o ConnectTimeout=30 "$h" "$@" 2>/dev/null) && { echo "$out"; return 0; }
    done
    return 1
}

say "loop start (poll ${POLL}s, expect ${EXPECT_TAGS} tags)"
merge_submitted=0
while true; do
    if ! c=$(issh "THROTTLE_SUBMIT=1 MAXQ=${MAXQ:-600} IHEP_MEM=${IHEP_MEM:-8000} bash $IR/scripts/_ihep_lc_cycle.sh"); then
        say "SSH_FAILED all login nodes"; sleep 300; continue
    fi
    echo "$c" >> "$LOG"
    echo "$c" | grep ' STOP ' | while read -r l; do
        grep -qF "$l" "$REG" || echo "$(date '+%F %H:%M') ALERT prod_lepclean IHEP $l" >> "$REG"
    done
    st=$(echo "$c" | grep ' STATUS ' | tail -1)
    tot=$(echo "$st" | sed -n 's/.*tags_complete=[0-9]*\/\([0-9]*\).*/\1/p')
    if echo "$c" | grep -q ALL_FILES_DONE && [ "${tot:-0}" -ge "$EXPECT_TAGS" ]; then
        if [ "$merge_submitted" -eq 0 ]; then
            cmd="hep_sub $IR/scripts/_ihep_lc_merge_job.sh -g cms -m 8000 -o $IR/logs/ihep_friend_lepclean/merge_job.out -e $IR/logs/ihep_friend_lepclean/merge_job.err"
            r=$(issh "export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:\$PATH; rm -f $IR/logs/ihep_friend_lepclean/merge_report.txt; $cmd")
            say "MERGE submitted :: $cmd :: $r"; merge_submitted=1
        fi
        rep=$(issh "cat $IR/logs/ihep_friend_lepclean/merge_report.txt 2>/dev/null")
        if echo "$rep" | grep -q '^MERGE_DONE'; then
            echo "$rep" > "$LOGD/ihep_merge_report.txt"
            if echo "$rep" | grep -q '^MERGE_DONE rc=0'; then
                bad=0
                for tag in $(echo "$rep" | awk '$1=="PASS"{print $2}'); do
                    mkdir -p "$COUT/$tag/merged"
                    rsync -a "ihep1:$IOUT/$tag/merged/merged_nominal.parquet" "$COUT/$tag/merged/" 2>/dev/null \
                      || rsync -a "ihep2:$IOUT/$tag/merged/merged_nominal.parquet" "$COUT/$tag/merged/" 2>/dev/null
                    want=$(echo "$rep" | awk -v t="$tag" '$2==t{for(i=1;i<=NF;i++) if($i ~ /^merged_rows=/){sub("merged_rows=","",$i); print $i}}')
                    got=$($PY -c "import pyarrow.parquet as pq;print(pq.ParquetFile('$COUT/$tag/merged/merged_nominal.parquet').metadata.num_rows)" 2>/dev/null || echo -1)
                    [ "$got" = "$want" ] || { bad=$((bad+1)); say "SYNC_MISMATCH $tag ihep=$want cern=$got"; }
                done
                n=$(echo "$rep" | grep -c '^PASS')
                if [ "$bad" -eq 0 ] && [ "$n" -ge "$EXPECT_TAGS" ]; then
                    say "SYNCED $n tags to $COUT, rows verified"
                    echo "$(date '+%F %H:%M') DONE  prod_lepclean IHEP 2022/23 merged+synced ($n tags) -> $LOGD/ihep_merge_report.txt" >> "$REG"
                else
                    echo "$(date '+%F %H:%M') FAIL  prod_lepclean IHEP sync: pass=$n mismatch=$bad -> $LOG" >> "$REG"
                fi
            else
                echo "$(date '+%F %H:%M') FAIL  prod_lepclean IHEP merge -> $LOGD/ihep_merge_report.txt" >> "$REG"
            fi
            say "loop end"; exit 0
        fi
        sleep 300; continue
    fi
    sleep "$POLL"
done
