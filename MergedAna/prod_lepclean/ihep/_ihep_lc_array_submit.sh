#!/bin/bash
# Submit the lepclean friend re-production as per-FILE hep_sub arrays, one per tag.
#   _ihep_lc_array_submit.sh <tag> [<tag> ...]     (ALL = every tag in the catalog)
# Each tag with N files -> `hep_sub -n N`, job idx = %{ProcId}. N_CORES=1, no -wt.
# Dry-run unless APPLY=1: prints the full hep_sub command lines first.
set -uo pipefail
export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
R=/publicfs/cms/user/laipeizhu/HZa/HiggsZaAna/HiggsDNA
WORKER=$R/scripts/_ihep_lc_file_worker.sh
CATALOG=$R/metadata/samples/za_merged_data_2223_perds.json
MEM=${MEM:-3500}
LOGDIR=$R/logs/ihep_friend_lepclean
mkdir -p "$LOGDIR"; cd "$LOGDIR"

if [ "${1:-}" = "ALL" ]; then
  tags=$(python3 -c "import json;print(' '.join(sorted(json.load(open('$CATALOG')).keys())))")
else
  tags="$*"
fi
[ -z "$tags" ] && { echo "no tags"; exit 1; }

total=0
for tag in $tags; do
  N=$(python3 -c "
import json
c=json.load(open('$CATALOG')); t=c['$tag']; f=list(t['files'].values())[0]
assert isinstance(f,list) and f and f[0].endswith('.root')
print(len(f))") || { echo "SKIP $tag (catalog not a PFN list)"; continue; }
  cmd=(hep_sub "$WORKER" -g cms -m "$MEM" -n "$N" -argu "$tag" "%{ProcId}"
       -o "$LOGDIR/${tag}_f%{ProcId}.out" -e "$LOGDIR/${tag}_f%{ProcId}.err")
  echo "${cmd[*]}"
  if [ "${APPLY:-0}" = "1" ]; then
    "${cmd[@]}" 2>&1 | tail -1
    echo "$(date '+%F %T') SUBMIT $tag N=$N mem=$MEM" >> "$LOGDIR/submissions.txt"
  fi
  total=$((total+N))
done
echo "=== $total file-jobs (mem=${MEM} MB)$([ "${APPLY:-0}" = "1" ] || echo ' -- dry-run, APPLY=1 to submit') ==="
