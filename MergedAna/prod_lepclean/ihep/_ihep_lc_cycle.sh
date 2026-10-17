#!/bin/bash
# ONE cycle of the IHEP lepclean friend production watch. Runs on an IHEP login node,
# called over ssh by the CERN-side loop (ihep_loop.sh): IHEP login nodes kill long
# local loops, so no loop lives here; state is kept in files.
#
# Per cycle:
#  1. held jobs: memory -> hep_edit -m 8000 -> 16000 -> 20000 + hep_release (held-job
#     rule); still held at 20000 -> left held, reported (STOP line).
#     Non-memory holds are not touched, only reported.
#  2. when no job of a tag is in the queue, every file idx must have a chunk or an
#     EMPTY marker; a missing idx (job failed: cceos / transient I/O / strict join) is
#     resubmitted as a single job, at most 10 times per idx, same memory as its last
#     known request. Over 10 -> reported, not resubmitted.
#  3. one status line; "ALL_FILES_DONE" when every tag is complete and nothing queued.
# Every action is printed with job id, reason, old/new memory and the exact command.
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
R=/publicfs/cms/user/laipeizhu/HZa/HiggsZaAna/HiggsDNA
W=$R/scripts/_ihep_lc_file_worker.sh
CAT=$R/metadata/samples/za_merged_data_2223_perds.json
OUT=/publicfs/cms/user/laipeizhu/Eos/HZa_merged/parquet_friend_lepclean
LOGDIR=$R/logs/ihep_friend_lepclean
ST=$LOGDIR/state; mkdir -p "$ST"
MEMF=$ST/mem.txt; RESUB=$ST/resub.txt; touch "$MEMF" "$RESUB"
BASEMEM=${IHEP_MEM:-$(awk '/SUBMIT/{m=$NF} END{sub("mem=","",m); print m+0}' "$LOGDIR/submissions.txt")}
ts=$(date '+%m-%d %H:%M')

# One cycle at a time: the CERN loop retries on another login node when ssh times out,
# while the first cycle may still be resubmitting there -> duplicate jobs. mkdir lock,
# considered stale after 40 min.
LOCK=$ST/cycle.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  age=$(( $(date +%s) - $(stat -c %Y "$LOCK" 2>/dev/null || date +%s) ))
  if [ "$age" -lt 2400 ]; then echo "$ts LOCKED (another cycle running for ${age}s)"; exit 0; fi
  echo "$ts STALE_LOCK removed (${age}s)"; rmdir "$LOCK" 2>/dev/null; mkdir "$LOCK" || exit 0
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

q=$(hep_q -u laipeizhu 2>/dev/null | grep '_ihep_lc_file_worker')
[ $? -gt 1 ] && { echo "$ts HEP_Q_FAILED"; exit 0; }

# ---- 1. held ----
hold=$(hep_q -u laipeizhu -hold 2>/dev/null)
for jid in $(echo "$q" | awk '$6=="H"{print $1}'); do
  reason=$(echo "$hold" | awk -v j="$jid" '$1""==j""' | head -1 | cut -c1-300)
  key=$(echo "$q" | awk -v j="$jid" '$1""==j""{print $(NF-1)"/"$NF}')
  if echo "$reason" | grep -qiE 'memory|mem|oom'; then
    cur=$(awk -v k="$key" '$1==k{m=$2} END{print m+0}' "$MEMF"); [ "$cur" -eq 0 ] && cur=$BASEMEM
    if   [ "$cur" -lt 8000 ];  then new=8000
    elif [ "$cur" -lt 16000 ]; then new=16000
    elif [ "$cur" -lt 20000 ]; then new=20000
    else echo "$ts STOP held at 20000 MB: $jid $key | $reason"; continue; fi
    echo "$ts HELD $jid $key reason=[$reason] mem $cur -> $new :: hep_edit -m $new $jid ; hep_release $jid"
    hep_edit -m "$new" "$jid" >/dev/null 2>&1 && hep_release "$jid" >/dev/null 2>&1
    echo "$key $new" >> "$MEMF"
  else
    echo "$ts HELD_OTHER $jid $key (not touched) reason=[$reason]"
  fi
done

# ---- 2. per-tag completeness, resubmit missing idx of tags with nothing queued ----
python3 - "$CAT" "$OUT" "$q" "$LOGDIR/submissions.txt" <<'PY' > "$ST/missing.txt"
import json, os, sys
cat = json.load(open(sys.argv[1])); out = sys.argv[2]
queued = {}
for line in sys.argv[3].splitlines():
    p = line.split()
    if len(p) >= 2 and p[-1].isdigit(): queued.setdefault(p[-2], set()).add(int(p[-1]))
submitted = {l.split()[3] for l in open(sys.argv[4]) if " SUBMIT " in l}
for tag in sorted(cat):
    if tag not in submitted:
        continue
    n = len(list(cat[tag]["files"].values())[0])
    d = os.path.join(out, tag, "chunks")
    done = set()
    # A worker wipes its f<idx> dir at start, so a dir can vanish between listdir and
    # the inner listdir; one such race used to kill the whole scan (10-14 21:52 the
    # cycle reported 4 of 12 tags). Treat any per-dir error as "not done yet".
    try:
        entries = os.listdir(d) if os.path.isdir(d) else []
    except OSError:
        entries = []
    for f in entries:
        if not f.startswith("f"): continue
        fd = os.path.join(d, f)
        try:
            if os.path.exists(os.path.join(fd, "EMPTY")) or any(
                    os.path.exists(os.path.join(fd, s, "merged_nominal.parquet"))
                    for s in os.listdir(fd) if os.path.isdir(os.path.join(fd, s))):
                done.add(int(f[1:]))
        except (OSError, ValueError):
            pass
    q = queued.get(tag, set())
    miss = sorted(set(range(n)) - done - q)
    print(tag, n, len(done), len(q), ",".join(map(str, miss)))
PY

ntag=0; ncomplete=0; nq=$(echo "$q" | grep -c .); nresub=0; nfiles=0; ndone=0
while read -r tag n done queued miss; do
  ntag=$((ntag+1)); nfiles=$((nfiles+n)); ndone=$((ndone+done))
  [ "$done" -eq "$n" ] && { ncomplete=$((ncomplete+1)); continue; }
  # missing already excludes the idx still in the queue, so no need to wait for the
  # whole tag to drain (that left ENV_FETCH_FAIL holes waiting for hours)
  for i in ${miss//,/ }; do
    [ "$nresub" -ge "${MAXRESUB:-40}" ] && break   # cap reached: stop resubmitting this cycle but keep
                                                     # counting the remaining tags (break 2 truncated STATUS)
    k="$tag/$i"
    c=$(grep -c "^$k$" "$RESUB")
    if [ "$c" -ge 10 ]; then
      grep -q "^GAVEUP $k$" "$RESUB" || { echo "GAVEUP $k" >> "$RESUB"; echo "$ts STOP $k failed 10 resubmissions"; }
      continue
    fi
    mem=$(awk -v k="$k" '$1==k{m=$2} END{print m+0}' "$MEMF"); [ "$mem" -eq 0 ] && mem=$BASEMEM
    err=$(tail -3 "$LOGDIR/${tag}_f${i}.err" 2>/dev/null | tr '\n' ' ' | cut -c1-200)
    cmd="hep_sub $W -g cms -m $mem -argu $tag $i -o $LOGDIR/${tag}_f${i}.out -e $LOGDIR/${tag}_f${i}.err"
    echo "$ts RESUB#$((c+1)) $k mem=$mem :: $cmd | last err: $err"
    $cmd >/dev/null 2>&1 && echo "$k" >> "$RESUB"
    nresub=$((nresub+1))
  done
done < "$ST/missing.txt"

# ---- 3. throttled submission of the remaining tags ----
# cceos reads fail en masse at high concurrency (~3000 concurrent -> ~99% failures;
# feedback_ihep_higgsdna_must_throttle), so whole tags (one array each) are only
# released while the queue holds fewer than MAXQ jobs. Largest tags first.
MAXQ=${MAXQ:-600}
nsubmitted_new=0
if [ "${THROTTLE_SUBMIT:-0}" = "1" ]; then
  python3 - "$CAT" "$LOGDIR/submissions.txt" <<'PY' > "$ST/unsubmitted.txt"
import json, sys
cat = json.load(open(sys.argv[1]))
done = {l.split()[3] for l in open(sys.argv[2]) if " SUBMIT " in l}
rest = [t for t in cat if t not in done]
for t in sorted(rest, key=lambda t: -len(list(cat[t]["files"].values())[0])):
    print(t, len(list(cat[t]["files"].values())[0]))
PY
  inq=$nq
  while read -r tag n; do
    [ $((inq + n)) -gt "$MAXQ" ] && [ "$inq" -gt 0 ] && break
    echo "$ts SUBMIT_TAG $tag N=$n (queue $inq -> $((inq+n)), MAXQ=$MAXQ)"
    APPLY=1 MEM=$BASEMEM bash $R/scripts/_ihep_lc_array_submit.sh "$tag" | sed "s/^/$ts   /"
    inq=$((inq + n)); nsubmitted_new=$((nsubmitted_new+1))
  done < "$ST/unsubmitted.txt"
  nleft=$(( $(wc -l < "$ST/unsubmitted.txt") - nsubmitted_new ))
else
  nleft=0
fi

nr=$(echo "$q" | awk '$6=="R"' | grep -c .); ni=$(echo "$q" | awk '$6=="I"' | grep -c .); nh=$(echo "$q" | awk '$6=="H"' | grep -c .)
echo "$ts STATUS tags_complete=$ncomplete/$ntag files_done=$ndone/$nfiles queue=$nq run=$nr idle=$ni held=$nh resub_this_cycle=$nresub new_tags=$nsubmitted_new unsubmitted_tags=$nleft"
if [ "$ntag" -gt 0 ] && [ "$ncomplete" -eq "$ntag" ] && [ "$nq" -eq 0 ] && [ "$nleft" -eq 0 ] && [ "$nsubmitted_new" -eq 0 ]; then echo "$ts ALL_FILES_DONE"; fi
