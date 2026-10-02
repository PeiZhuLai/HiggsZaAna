#!/usr/bin/env bash
#
# V2 reviewer study, extension to the Appendix B mass points m_a >= 2 GeV.
# Identical to ../v2_trigdenom_20260927/run_trigdenom_driver.sh (m_a = 1 GeV) except for the
# sample catalog and the output/log directory names:
#   config  metadata/za_signal_run3_trigdenom_mAge2.json  (same tagger option
#           zgammas.study_lep_trigger_eff_otherleg = true; catalog za_trigdenom_mAge2.json)
#   output  /eos/home-p/pelai/HZa/trigdenom_V2_mAge2_20261002/Sig_MC_trigdenomV2mAge2[TEST]
#   logs    HiggsDNA/eos_logs/Sig_MC_trigdenomV2mAge2[TEST]  (CondorManager uses the last
#           component of the output dir, so production eos_logs/Sig_MC is never touched)
#
# MODE=test -> --short, SAMPLES (default mA_M2), YEARS (default 2024): one job per sample/era
# MODE=full -> every file of the 13 mass points x 5 eras (1044 NanoAOD files = 1044 jobs at fpo=1)
# DRY_RUN=1 (default) prints the HiggsDNA command without submitting.
#
# Launched as a systemd --user unit (lxplus el9 kills setsid/nohup on logout).
set -uo pipefail
MODE="${MODE:-test}"
DRY_RUN="${DRY_RUN:-1}"
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_mAge2_20261002
BASE=/eos/home-p/pelai/HZa/trigdenom_V2_mAge2_20261002
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
CFG=metadata/za_signal_run3_trigdenom_mAge2.json
ALL=mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30
case "$MODE" in
  test) OUT=$BASE/Sig_MC_trigdenomV2mAge2TEST; SMP="${SAMPLES:-mA_M2}"; YRS="${YEARS:-2024}"; SHORTF=1 ;;
  full) OUT=$BASE/Sig_MC_trigdenomV2mAge2; SMP="${SAMPLES:-$ALL}"; YRS="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"; SHORTF=0 ;;
  *) echo "bad MODE=$MODE"; exit 2 ;;
esac

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1   # dasgoclient, voms-proxy-info
set -u
export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$REPO"
export X509_USER_PROXY=/tmp/x509up_u175325
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
export HZA_STAGE_INPUTS=1
export HIGGSDNA_CONDOR_REQ_MEMORY=8000

echo "[$(date '+%F %T')] MODE=$MODE DRY_RUN=$DRY_RUN OUT=$OUT CFG=$CFG SAMPLES=$SMP YEARS=$YRS host=$(hostname)"
left=$(voms-proxy-info -file "$X509_USER_PROXY" -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
[ "${left:-0}" -lt 86400 ] && { echo "[ERROR] proxy under 24 h"; exit 1; }
if [ /tmp/x509up_u175325 -nt "$REPO/x509up_u175325" ]; then cp -p "$X509_USER_PROXY" "$REPO/x509up_u175325"; fi
openssl x509 -noout -enddate -in "$REPO/x509up_u175325"

mkdir -p "$OUT"
cd "$REPO"
OUTDIR="$OUT" CONFIG="$CFG" FPO=1 CLEAN_ANALYSIS_STATE=1 UNRETIRE_JOBS=1 MERGE_OUTPUTS=0 \
  N_CORES=1 SHORT=$SHORTF CONDOR_REQ_MEMORY=8000 SAMPLE_LIST="$SMP" YEARS="$YRS" DRY_RUN="$DRY_RUN" \
  bash scripts/run_ana_signal.sh
rc=$?
echo "[$(date '+%F %T')] run_ana_signal.sh rc=$rc (NOT a verdict; check summaries + .out payloads)"
[ "$DRY_RUN" = "0" ] && echo "DRIVER_DONE rc=$rc $(date '+%F %T')" > "$D/logs/driver_${MODE}.done"
exit $rc
