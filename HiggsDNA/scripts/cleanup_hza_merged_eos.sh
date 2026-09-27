#!/bin/bash
# =============================================================================
# cleanup_hza_merged_eos.sh -- tiered cleanup of HZa_merged on eoscms
#
# WHY: the zh group quota on /eos/cms/store/group/phys_susy is exceeded
#      (357.43 / 360 TB raw, 99.29%). HZa_merged holds 1.58 TB logical
#      (= 3.16 TB raw, replica x2). The full re-production with photon-lepton
#      cleaning needs headroom.
#
# DEFAULT IS DRY-RUN. Nothing is deleted unless you pass DRY_RUN=0.
#
#   bash cleanup_hza_merged_eos.sh                 # dry-run, tier A
#   TIER=B bash cleanup_hza_merged_eos.sh          # dry-run, tier A+B
#   TIER=A DRY_RUN=0 bash cleanup_hza_merged_eos.sh   # actually delete tier A
#
# Tiers
#   A  stray logs, empty dirs, smoke/probe leftovers   ~36 MB logi
#   B  superseded v3 products and training scratch     ~23 GB logi
#   C  large friend/merged products (REVIEW FIRST)     ~224 GB logi
#
# NEVER TOUCHED (hard-coded protection, see PROTECT below):
#   MLNanoAOD/              1232 GB  upstream ML photon reconstruction;
#                                    the tagger re-production does NOT regenerate it
#   parquet_merged_DNA_tmp/    4 GB  despite the name this is the LIVE default
#                                    SIG_ML_DIR in 1_grand_merged.sh:48
#   root_MVAcut/               8 GB  live ROOT_MVACUT in 1_grand_merged.sh:46
#
# input : none (reads EOS listings)
# output: none in dry-run; deletes under $BASE when DRY_RUN=0
# =============================================================================
set -uo pipefail

BASE="${BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged}"
MGM="${MGM:-root://eoscms.cern.ch}"
DRY_RUN="${DRY_RUN:-1}"
TIER="${TIER:-A}"

PROTECT="MLNanoAOD parquet_merged_DNA_tmp root_MVAcut"

freed=0

is_protected() {
    local name="$1"
    for p in $PROTECT; do [ "$name" = "$p" ] && return 0; done
    return 1
}

dirsize() {
    eos "$MGM" ls -l "$BASE/" 2>/dev/null | awk -v n="$1" '$NF==n{print $(NF-4); exit}'
}

drop_dir() {
    local name="$1" why="$2"
    if is_protected "$name"; then
        echo "  [PROTECTED] $name -- refusing ($why)"
        return
    fi
    local sz; sz=$(dirsize "$name")
    if [ -z "$sz" ]; then
        echo "  [absent]    $name"
        return
    fi
    printf "  [%s] %-40s %9.2f GB   %s\n" \
        "$([ "$DRY_RUN" = 0 ] && echo DELETE || echo "dry-run")" \
        "$name/" "$(echo "$sz/1073741824" | bc -l)" "$why"
    freed=$((freed + sz))
    if [ "$DRY_RUN" = 0 ]; then
        eos "$MGM" rm -r "$BASE/$name" || echo "  !! rm failed: $name"
    fi
}

drop_file() {
    local name="$1"
    local sz; sz=$(eos "$MGM" ls -l "$BASE/" 2>/dev/null | awk -v n="$name" '$NF==n{print $(NF-4); exit}')
    [ -z "$sz" ] && return
    freed=$((freed + sz))
    [ "$DRY_RUN" = 0 ] && { eos "$MGM" rm "$BASE/$name" || echo "  !! rm failed: $name"; }
}

echo "=================================================================="
echo " HZa_merged EOS cleanup   BASE=$BASE"
echo " DRY_RUN=$DRY_RUN   TIER=$TIER"
echo "=================================================================="

# ---------------------------------------------------------------- TIER A ----
echo
echo "--- Tier A: stray logs / empty dirs / smoke leftovers ---"

echo "  stray top-level log+txt files:"
nlog=0; slog=0
while read -r sz nm; do
    [ -z "$nm" ] && continue
    case "$nm" in
        *.log|*.txt)
            nlog=$((nlog+1)); slog=$((slog+sz))
            [ "$DRY_RUN" = 0 ] && eos "$MGM" rm "$BASE/$nm" >/dev/null 2>&1
            ;;
    esac
done < <(eos "$MGM" ls -l "$BASE/" 2>/dev/null | awk '$1!~/^d/{print $(NF-4), $NF}')
printf "    %s %d files, %.2f MB\n" \
    "$([ "$DRY_RUN" = 0 ] && echo deleted || echo "would delete")" "$nlog" "$(echo "$slog/1048576"|bc -l)"
freed=$((freed + slog))

echo "  zero-byte sample dirs under MLNanoAOD/ (never produced):"
nempty=0
while read -r nm; do
    [ -z "$nm" ] && continue
    nempty=$((nempty+1))
    echo "    $nm"
    [ "$DRY_RUN" = 0 ] && eos "$MGM" rm -r "$BASE/MLNanoAOD/$nm" >/dev/null 2>&1
done < <(eos "$MGM" ls -l "$BASE/MLNanoAOD/" 2>/dev/null | awk '$1~/^d/ && $(NF-4)==0 {print $NF}')
echo "    ($nempty dirs)"

drop_dir gpuprobe                          "one-off GPU probe, no references"
drop_dir parquet_friend_ML_smoke           "smoke test superseded by parquet_friend_ML"
drop_dir parquet_friend_ML_lepclean_smoke  "lepton-cleaning smoke test, conclusion recorded in doc"
drop_dir parquet_merged_DNA_lepclean_smoke "lepton-cleaning smoke test, conclusion recorded in doc"
drop_dir parquet_merged_DNA_lepclean_full  "aborted lepclean run (0.3 MB)"

# ---------------------------------------------------------------- TIER B ----
if [ "$TIER" = "B" ] || [ "$TIER" = "C" ]; then
echo
echo "--- Tier B: superseded v3 products / training scratch ---"
drop_dir MLNanoAOD_v3          "superseded by MLNanoAOD_v4 (same 9 mass points, 2026-08-18)"
drop_dir parquet_merged_DNA_v3 "superseded by parquet_merged_DNA_v4"
drop_dir train_packed          "superseded by train_packed_v4"
drop_dir train_packed_merged   "superseded by train_packed_v4"
drop_dir train_npz_r3          "training scratch, no references"
drop_dir train_dump            "training scratch, regenerable from MLNanoAOD"
fi

# ---------------------------------------------------------------- TIER C ----
if [ "$TIER" = "C" ]; then
echo
echo "--- Tier C: LARGE. Only after the re-production lands. ---"
echo "    Read the notes before enabling; these are NOT obviously dead."
drop_dir parquet_friend    "184 GB. resolved-branch friend; still referenced by DATA_ML_GLOB"
drop_dir parquet_friend_ML "21 GB. background ML friend; obsoleted by the cleaned re-production"
fi

echo
printf "=== %s: %.2f GB logical (%.2f GB raw, replica x2) ===\n" \
    "$([ "$DRY_RUN" = 0 ] && echo freed || echo "would free")" \
    "$(echo "$freed/1073741824"|bc -l)" "$(echo "$freed*2/1073741824"|bc -l)"
[ "$DRY_RUN" = 1 ] && echo "(dry-run -- nothing was deleted; pass DRY_RUN=0 to act)"
