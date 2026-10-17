#!/bin/bash
# hep_sub job: merge + reconcile every lepclean friend tag on IHEP (ihep_lc_merge.py).
# Writes the per-tag verdict lines to logs/ihep_friend_lepclean/merge_report.txt and a
# final MERGE_DONE rc=<rc> line (the CERN-side loop waits for it).
set -o pipefail
R=/publicfs/cms/user/laipeizhu/HZa/HiggsZaAna/HiggsDNA
PACK=/publicfs/cms/user/laipeizhu/App/Conda/pack/hza_ana_pack.tar.gz
REP=$R/logs/ihep_friend_lepclean/merge_report.txt
ENVROOT="${TMPDIR:-/tmp}/hza_lc_merge_$$"
mkdir -p "$ENVROOT"; cp "$PACK" "$ENVROOT/p.tgz" && tar xzf "$ENVROOT/p.tgz" -C "$ENVROOT" && rm -f "$ENVROOT/p.tgz"
source "$ENVROOT/bin/activate"; "$ENVROOT/bin/conda-unpack" 2>/dev/null || true
"$ENVROOT/bin/python" $R/scripts/ihep_lc_merge.py $R/metadata/samples/za_merged_data_2223_perds.json \
    /publicfs/cms/user/laipeizhu/Eos/HZa_merged/parquet_friend_lepclean "$@" > "$REP.tmp" 2>&1
rc=$?
echo "MERGE_DONE rc=$rc $(date '+%F %T')" >> "$REP.tmp"
mv "$REP.tmp" "$REP"
rm -rf "$ENVROOT"
exit $rc
