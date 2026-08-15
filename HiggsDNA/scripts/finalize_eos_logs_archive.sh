#!/bin/bash
# Finalize the eos_logs .out archival:
#   wait for tar -> verify -> only then delete the archived originals.
#
# The deletion is gated on THREE checks, all of which must pass:
#   1. tar -tzf exits 0                       (archive is readable end to end)
#   2. .out count inside == keep_out.txt count (nothing was dropped while taring)
#   3. every path in keep_out.txt is present inside the archive (exact set match)
# Any failure leaves the originals untouched and records the reason.

set -uo pipefail

H="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
TARBALL="$H/eos_logs_out_archive_20260811.tar.gz"
KEEP="/tmp/pelai/claude-175325/-afs-cern-ch-user-p-pelai/9a28ea12-d78b-4d45-9ea9-99e7dff90fa4/scratchpad/keep_out.txt"
STATUS="$H/scripts/.eos_logs_archive.status"
LIST="$H/scripts/.eos_logs_archive_contents.txt"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$STATUS"; }

: > "$STATUS"
say "waiting for tar to finish"
while pgrep -x tar > /dev/null; do sleep 20; done
say "tar process gone; archive size: $(stat -c %s "$TARBALL" 2>/dev/null) bytes"

# --- gate 1: archive readable ---
if ! tar -tzf "$TARBALL" > "$LIST" 2>/dev/null; then
  say "FAIL gate1: tar -tzf returned non-zero (archive truncated or corrupt). Originals kept."
  exit 1
fi
say "gate1 OK: archive readable"

# --- gate 2: counts match ---
n_tar=$(grep -c '\.out$' "$LIST")
n_keep=$(wc -l < "$KEEP")
say "gate2: .out in archive = $n_tar ; keep list = $n_keep"
if [ "$n_tar" -ne "$n_keep" ]; then
  say "FAIL gate2: count mismatch. Originals kept."
  exit 1
fi
say "gate2 OK"

# --- gate 3: exact set match (paths in archive are relative to $H) ---
missing=$(comm -23 \
  <(sed "s|^$H/||" "$KEEP" | sort) \
  <(grep '\.out$' "$LIST" | sort) | wc -l)
say "gate3: keep-list entries absent from archive = $missing"
if [ "$missing" -ne 0 ]; then
  say "FAIL gate3: archive does not contain every kept .out. Originals kept."
  exit 1
fi
say "gate3 OK: all three gates passed"

# --- delete originals ---
before=$(du -sb "$H/eos_logs" | cut -f1)
xargs -a "$KEEP" -d '\n' -r rm -f
after=$(du -sb "$H/eos_logs" | cut -f1)
remaining=$(find "$H/eos_logs" -name '*.out' | wc -l)
say "deleted archived .out; eos_logs $((before/1024/1024)) MB -> $((after/1024/1024)) MB ; .out remaining = $remaining"
say "DONE"
