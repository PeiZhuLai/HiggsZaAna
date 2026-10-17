#!/bin/bash
# IHEP per-FILE worker for the HZa_merged lepton-cleaned (dR>0.3) friend re-production
# of 2022/2023 data. Args: $1 = tag, $2 = file index (0-based in the PFN catalog).
#
# Same flow as _ihep_friend_file_worker.sh (August production), with:
#   * output to parquet_friend_lepclean/ (the old parquet_friend/ is left untouched);
#   * HZA_MLPHOTON_STRICT=1: a failed MLPhoton join fails the job instead of writing
#     events with no MLPhoton (the 2024 Muon1 RunD failure mode);
#   * the per-file output dir is wiped first: a job OOM-killed while writing
#     analysis_manager.pkl leaves a truncated pickle, and every retry then dies with
#     EOFError on it;
#   * a file with zero selected events gets an EMPTY marker and exits 0, so it is
#     not retried forever (run_analysis writes no parquet for it).
# Mode-B env (node-local conda-pack), nano via cceos, MLNanoAOD from scratchfs2.
set -o pipefail
TAG="${1:?usage: _ihep_lc_file_worker.sh <tag> <idx>}"
IDX="${2:?usage: _ihep_lc_file_worker.sh <tag> <idx>}"

REPO=/publicfs/cms/user/laipeizhu/HZa/HiggsZaAna/HiggsDNA
MLNANO=/scratchfs2/cms/laipeizhu/Eos/HZa_merged/MLNanoAOD
OUTBASE=/publicfs/cms/user/laipeizhu/Eos/HZa_merged/parquet_friend_lepclean
CATALOG=metadata/samples/za_merged_data_2223_perds.json
PACK=/publicfs/cms/user/laipeizhu/App/Conda/pack/hza_ana_pack.tar.gz
PROXY=/publicfs/cms/user/laipeizhu/x509up_laipeizhu

echo "==== IHEP lepclean FILE worker: $TAG idx=$IDX  @ $(hostname) $(date) ===="

ENVROOT="${TMPDIR:-/tmp}/hza_lc_env_${TAG}_${IDX}_$$"
# Some nodes have a full local scratch (10-14: hxmt127/bws0976/hxmt116, 476 jobs died in
# seconds with "No space left on device" copying the 1.5 GB pack). Failing fast turns
# such a node into a black hole that eats the queue, so hold the slot for a while
# before failing; the cycle script resubmits the idx and it lands elsewhere.
free_kb=$(df -Pk "${TMPDIR:-/tmp}" | awk 'NR==2{print $4}')
if [ "${free_kb:-0}" -lt 8000000 ]; then
  echo "SCRATCH_FULL $(hostname) free=$((free_kb/1024)) MB -- sleeping 1200 s before failing"
  sleep 1200; exit 75
fi
mkdir -p "$ENVROOT"
cp "$PACK" "$ENVROOT/pack.tar.gz" || { echo "ENV_FETCH_FAIL"; exit 110; }
tar xzf "$ENVROOT/pack.tar.gz" -C "$ENVROOT" || { echo "ENV_EXTRACT_FAIL"; exit 111; }
rm -f "$ENVROOT/pack.tar.gz"
source "$ENVROOT/bin/activate"
"$ENVROOT/bin/conda-unpack" 2>/dev/null || true

export PYTHONPATH="$REPO"
export X509_USER_PROXY="$PROXY"
export HDNA_REDIRECTOR="root://cceos.ihep.ac.cn:1094//eos/ihep/cms"
export PYTHONNOUSERSITE=1
export HOME=/publicfs/cms/user/laipeizhu/dasmaps_home
export HZA_MLPHOTON_STRICT=1
# IHEP workers have no xrdcp on PATH (staging fails with rc=32512 = exit 127); stream
# from cceos as in the August production.
export HZA_STAGE_INPUTS="${HZA_STAGE_INPUTS:-0}"
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1 || true
source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh >/dev/null 2>&1 || true
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
export X509_VOMS_DIR=/cvmfs/grid.cern.ch/etc/grid-security/vomsdir
cd "$REPO"

CFGDIR=metadata/friend_tmp_configs_ihep_lc
mkdir -p "$CFGDIR"
subcat="$CFGDIR/${TAG}_f${IDX}.catalog.json"
cfg="$CFGDIR/${TAG}_f${IDX}.json"
"$ENVROOT/bin/python" - "$TAG" "$IDX" "$cfg" "$subcat" "$CATALOG" "$MLNANO/$TAG" "metadata/parent_maps/${TAG}.json" <<'PY'
import json, os, sys
tag, idx, out, subcat, catalog, ml_dir, pm = sys.argv[1:]
idx = int(idx)
cat = json.load(open(catalog))
era = list(cat[tag]["files"].keys())[0]
pfns = cat[tag]["files"][era]
if not (isinstance(pfns, list) and pfns and pfns[0].endswith(".root")):
    raise SystemExit("ERR: catalog[%s][%s] is not an explicit PFN list" % (tag, era))
if idx >= len(pfns):
    raise SystemExit("ERR: idx %d out of range (n=%d)" % (idx, len(pfns)))
sub = {tag: {"fpo": cat[tag].get("fpo", 1), "files": {era: [pfns[idx]]}}}
json.dump(sub, open(subcat, "w"), indent=2)
d = json.load(open("metadata/za_merged_data_ML_2024.json"))
d["samples"]["catalog"] = subcat
d["samples"]["sample_list"] = [tag]
d["samples"]["years"] = [era]
d["mlphoton_friend"] = True
d["mlphoton_parent_map"] = os.path.abspath(pm)
d["mlphoton_dir"] = ml_dir
d["mlphoton_tag"] = tag
json.dump(d, open(out, "w"), indent=2)
print("idx=%d era=%s file=%s" % (idx, era, os.path.basename(pfns[idx])))
PY
rc=$?
[ $rc -ne 0 ] && { echo "CONFIG_BUILD_FAIL rc=$rc"; rm -rf "$ENVROOT"; exit $rc; }

outdir="$OUTBASE/$TAG/chunks/f${IDX}"
rm -rf "$outdir"
mkdir -p "$outdir"
"$ENVROOT/bin/python" scripts/run_analysis.py \
  --config "$cfg" --log-level INFO --n_cores 1 \
  --output_dir "$outdir" --batch_system local --merge_outputs
rc=$?
echo "==== run_analysis rc=$rc  @ $(date) ===="
if [ $rc -eq 0 ]; then
  chunk=$(find "$outdir" -name merged_nominal.parquet 2>/dev/null | head -1)
  if [ -n "$chunk" ]; then
    echo "CHUNK_OK $chunk"
  elif find "$outdir" -name '*summary_job*.json' 2>/dev/null | grep -q .; then
    touch "$outdir/EMPTY"; echo "CHUNK_EMPTY (job completed, 0 selected events)"
  else
    echo "CHUNK_MISSING"; rc=42
  fi
fi
rm -rf "$ENVROOT" 2>/dev/null || true
exit $rc
