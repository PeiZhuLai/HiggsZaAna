#!/usr/bin/env bash
#
# V2 reviewer study (Vinay; also Eiko/Sam on Fig. 5): dilepton-trigger efficiency vs
# leading (subleading) lepton pT with the OTHER leg required above its own double-lepton
# threshold. Signal m_a = 1 GeV only, all five Run 3 eras.
#
# WHAT IS DIFFERENT FROM THE PRODUCTION SIGNAL RUN
#   Only the tagger option zgammas.study_lep_trigger_eff_otherleg = true
#   (metadata/za_signal_run3_trigdenom.json). It ADDS trigeff curves with ord labels
#   leadN2/subleadN2 (>= 2 selected same-flavor leptons) and leadOL/subleadOL (>= 2 and the
#   other leg above 25/15 GeV for e, 20/10 GeV for mu). The old lead/sublead curves are still
#   produced in the same job, so old and new come from identical events.
#   Environment mirrors the FSR-fix production (submit_fsr_reproduction.sh):
#   hza_ana env, B-mode pack, copy-then-read staging, fpo=1, RequestMemory 3000.
#
# WHY a separate catalog (metadata/samples/za_trigdenom_mA1.json)
#   SampleManager rewrites <catalog>_sample_manager_full.json; running with zgamma_tutorial.json
#   would rewrite the validated production's file.
# WHY the last output-dir component is Sig_MC_trigdenomV2[TEST]
#   CondorManager writes logs to HiggsDNA/eos_logs/<last component of output_dir>;
#   ".../Sig_MC" would write into the production eos_logs/Sig_MC.
#
# MODE=test -> --short, 2024 + 2023preBPix (1 job each), output .../Sig_MC_trigdenomV2TEST
# MODE=full -> all jobs, five eras,                        output .../Sig_MC_trigdenomV2
#
# Launched as a systemd --user unit (lxplus el9 kills setsid/nohup on logout).
set -uo pipefail
MODE="${MODE:-test}"
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_20260927
BASE=/eos/home-p/pelai/HZa/trigdenom_V2_20260927
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
CFG=metadata/za_signal_run3_trigdenom.json
case "$MODE" in
  test) OUT=$BASE/Sig_MC_trigdenomV2TEST; YRS=2024,2023preBPix; SHORTF=1 ;;
  full) OUT=$BASE/Sig_MC_trigdenomV2; YRS=2022preEE,2022postEE,2023preBPix,2023postBPix,2024; SHORTF=0 ;;
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
export HIGGSDNA_CONDOR_REQ_MEMORY=3000

echo "[$(date '+%F %T')] MODE=$MODE OUT=$OUT CFG=$CFG YEARS=$YRS host=$(hostname)"
left=$(voms-proxy-info -file "$X509_USER_PROXY" -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
[ "${left:-0}" -lt 86400 ] && { echo "[ERROR] proxy under 24 h"; exit 1; }
# The repo copy is what ships with each job (jobs.py GRID_PROXY); /tmp being fresh is not enough.
# Only refresh it if /tmp is newer (another session copies the same file).
if [ /tmp/x509up_u175325 -nt "$REPO/x509up_u175325" ]; then cp -p "$X509_USER_PROXY" "$REPO/x509up_u175325"; fi
openssl x509 -noout -enddate -in "$REPO/x509up_u175325"

mkdir -p "$OUT"
cd "$REPO"
OUTDIR="$OUT" CONFIG="$CFG" FPO=1 CLEAN_ANALYSIS_STATE=1 UNRETIRE_JOBS=1 MERGE_OUTPUTS=0 \
  N_CORES=1 SHORT=$SHORTF CONDOR_REQ_MEMORY=3000 SAMPLE_LIST=mA_M1 YEARS="$YRS" DRY_RUN="${DRY_RUN:-0}" \
  bash scripts/run_ana_signal.sh
rc=$?
echo "[$(date '+%F %T')] run_ana_signal.sh rc=$rc (NOT a verdict; check summaries + .out payloads)"
echo "DRIVER_DONE rc=$rc $(date '+%F %T')" > "$D/logs/driver_${MODE}.done"
exit $rc
