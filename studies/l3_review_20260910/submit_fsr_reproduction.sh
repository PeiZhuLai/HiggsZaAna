#!/usr/bin/env bash
#
# HZa Run 3 full re-production after the FSR-recovery fixes (+ LHE weight wiring).
#
# WHY
#   Two bugs were found while producing the FSR before/after plot the new L3
#   convener asked for, and both are in the code the current samples were made
#   with:
#     A) higgs_dna/selections/photon_selections.py: the three delta_R calls in
#        select_resolved_fsr_photons were not wrapped in ak.fill_none(...,True).
#        The selected-electron collection arrives as an option type, so the mask
#        came back as ?bool; FsrPhoton[mask] became an option array, ak.num
#        counted the None placeholders as photons (n_fsr > 0 in 12-19 % of muon
#        events) and the assignment found nothing to add -- the FSR recovery
#        never ran (dressing rate 0.008-0.054 %).
#     B) Once A was fixed, the recovery started eating the SUB-LEADING ALP
#        photon: the signal-photon veto only covered photons_sorted[:, :1]. At
#        m_a = 30 GeV, 72 % of the dressed candidates were the sub-leading ALP
#        photon itself and sigma_eff(m_llgg) got 16 % WORSE. The veto now covers
#        every photon passing the signal criteria, which is the two-photon
#        generalization of nano2pico's FsrSeparationReq = 0.2.
#   With A+B the recovery behaves: dressing 1.9 % (m_a=1) / 3.0 % (m_a=30), no
#   signal photon absorbed, residual FSR photon pT median 6.8 GeV (soft), and
#   sigma_eff(m_llgg) improves 2.8 % at m_a=30. The fit observable is the plain
#   H_mass (no refit branch reaches the fit inputs), so this propagates to the
#   signal model and the limits -- hence a full re-production, data included:
#   the muon-channel m(llgg) moves in data too.
#
#   The same pass also wires up the LHE weights (analysis.py attach_lhe_weights
#   + 12 fields in metadata/za_signal_run3.json). LHEScaleWeight_* / LHEPdfWeight_*
#   never existed in the output before -- load_events silently drops branches
#   that are not in the NanoAOD, and those flat names only exist as arrays there.
#   Doing it now avoids a second signal-only pass for the QCD-scale acceptance
#   study (L3 item 8).
#
# SAFETY
#   Writes to a NEW staging directory, so the current production -- the one the
#   AN is based on -- is untouched until the new one is reconciled.
#
# USAGE
#   bash submit_fsr_reproduction.sh            # dry run, prints the commands
#   DRY_RUN=0 bash submit_fsr_reproduction.sh  # actually submit
#   STAGE=signal DRY_RUN=0 bash submit_fsr_reproduction.sh   # one stage only
#
set -euo pipefail

DRY_RUN="${DRY_RUN:-1}"
STAGE="${STAGE:-all}"          # all | signal | bkg | data

REPO="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"

# HiggsDNA bakes the worker-side interpreter into the condor wrapper as
# "$CONDA_PREFIX/bin/python" (jobs.py:414 via misc_utils.get_HiggsDNA_conda,
# which is literally `echo $CONDA_PREFIX`). Exporting PATH alone is NOT enough:
# with CONDA_PREFIX unset the wrapper ends up with <repo>/bin/python, every job
# dies in ~2 s with exit 127 ("No such file or directory"), and the queue looks
# empty while nothing is produced. run_analysis.py still returns 0 in that case
# -- it counts retired jobs as complete -- so the only reliable check is the
# output itself. All three variables are required.
ENVDIR="${ENVDIR:-/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana}"
export CONDA_PREFIX="${ENVDIR}"
export PATH="${ENVDIR}/bin:${PATH}"
export PYTHONPATH="${REPO}"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"

# B-mode: every worker fetches the conda-pack tarball to node-local scratch
# instead of reading the env off /eos. Without it, ~1089 concurrent jobs made
# EOS return I/O errors for 66 % of attempts (749x OSError Errno 5, 98x
# ImportError on .so files). The wrapper (condor/lxplus/exe_template.sh) already
# has the branch; it only needs this variable, and submit_template.txt is
# getenv=True so it reaches the job. higgs_dna itself still comes from AFS via
# PYTHONPATH, so local code changes take effect -- the tarball only carries
# python/ROOT/numpy. See ref_hza_bmode_pack.
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"

# fpo and memory stay at the values the original production used. Raising fpo to
# cut EOS env reads was the wrong tool -- B-mode fixes that directly, while fpo
# drives memory (fpo=16 measured 23 GB/job against a 12 GB request).
export HIGGSDNA_CONDOR_REQ_MEMORY="${HIGGSDNA_CONDOR_REQ_MEMORY:-3000}"

if [[ ! -x "${ENVDIR}/bin/python" ]]; then
    echo "[ERROR] ${ENVDIR}/bin/python is missing -- the condor wrapper would bake a bad interpreter." >&2
    exit 1
fi
TAG="${TAG:-fsrfix}"
STAGING="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_${TAG}"
LOGDIR="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_${TAG}"

# Expected job counts, taken from the job directories of the current production.
# The sample manager recomputes them at submit time; these are for reporting.
declare -A EXPECTED=( [signal]=1104 [bkg]=3287 [data]=3475 )

echo "=== HZa Run 3 re-production (${TAG}) ==="
echo "repo        : ${REPO}"
echo "staging out : ${STAGING}/{Sig_MC,Bkg_MC,Data}"
echo "job logs    : ${LOGDIR}"
echo "driver logs : ${LOGDIR}/<stage>.log"
echo "stages      : ${STAGE}"
echo "expected    : signal ${EXPECTED[signal]} + bkg ${EXPECTED[bkg]} + data ${EXPECTED[data]} jobs"
echo "dry run     : ${DRY_RUN}"
echo

# Create the output tree from the login node. Condor workers cannot create
# parent directories under /eos/project ([3010] Operation not permitted); they
# can write files into directories that already exist.
make_dirs() {
    for d in "${STAGING}/Sig_MC" "${STAGING}/Bkg_MC" "${STAGING}/Data" "${LOGDIR}"; do
        if [[ "${DRY_RUN}" == "1" ]]; then
            echo "mkdir -p ${d}"
        else
            mkdir -p "${d}"
        fi
    done
}

run_stage() {
    local name="$1" script="$2" outdir="$3"
    if [[ "${STAGE}" != "all" && "${STAGE}" != "${name}" ]]; then
        echo "--- skip ${name}"
        return 0
    fi
    echo "--- ${name}: expecting ~${EXPECTED[$name]} condor jobs -> ${outdir}"
    # higgs_dna is not pip-installed into the conda env, and run_analysis.py
    # lives in scripts/ so Python puts scripts/ on sys.path, not the repo root.
    # Without this the driver dies with ModuleNotFoundError: No module named
    # 'higgs_dna' before a single job is submitted.
    local cmd=(env "CONDA_PREFIX=${ENVDIR}" "PATH=${ENVDIR}/bin:${PATH}"
               "PYTHONPATH=${REPO}" "HZA_BMODE_PACK=${HZA_BMODE_PACK}"
               "CONDOR_REQ_MEMORY=${HIGGSDNA_CONDOR_REQ_MEMORY}"
               "OUTDIR=${outdir}" "DRY_RUN=${DRY_RUN}"
               bash "${REPO}/scripts/${script}")
    if [[ "${DRY_RUN}" == "1" ]]; then
        printf '    '; printf '%q ' "${cmd[@]}"; printf '\n'
        ( cd "${REPO}" && "${cmd[@]}" ) | sed 's/^/    /'
    else
        ( cd "${REPO}" && "${cmd[@]}" ) 2>&1 | tee "${LOGDIR}/${name}.log"
    fi
    echo
}

make_dirs
echo
# Signal first: smallest, and the only stage that carries the new LHE fields.
run_stage signal run_ana_signal.sh "${STAGING}/Sig_MC"
run_stage bkg    run_ana_bkgmc.sh  "${STAGING}/Bkg_MC"
run_stage data   run_ana_data.sh   "${STAGING}/Data"

echo "=== done (dry run=${DRY_RUN}) ==="
if [[ "${DRY_RUN}" == "1" ]]; then
    echo "Nothing was submitted. Re-run with DRY_RUN=0 to submit."
fi
