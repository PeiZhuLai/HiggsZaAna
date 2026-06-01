#!/bin/bash
#
# For each low-mA signal mass point (M0p1..M0p9, M1), dump 10 events of
# EB RecHits from one MiniAOD file via cmsRun, then run the PyROOT 2D plotter
# in event mode to write 10 per-event heatmaps under
# Plot/plots/rechits_2d/{mass}/event/.
#
# Usage:
#     bash dump_and_plot_rechits_all_masses.sh
#

set -e

PROJECT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
CMSSW=${PROJECT}/CMSSW_15_0_14
CFG=${CMSSW}/src/RegressMergedPhoton/RecHitDumper/test/dump_rechits_cfg.py
DUMP_DIR=/eos/home-p/pelai/HZa/root_rechit   # ROOT dumps live on EOS
PLOT_BASE=${PROJECT}/Plot/plots/rechits_2d   # PNG/PDF plots on AFS
PLOTTER=${PROJECT}/Plot/scripts/plot_rechits_2d.py
N_EVENTS=10

# CMSSW + conda
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd ${CMSSW}/src
eval `scramv1 runtime -sh`
cd ${PROJECT}

declare -A DATASET
DATASET[M0p1]=/HZa-Zto2L-ato2G_Par-M-0p1_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p2]=/HZa-Zto2L-ato2G_Par-M-0p2_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p3]=/HZa-Zto2L-ato2G_Par-M-0p3_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p4]=/HZa-Zto2L-ato2G_Par-M-0p4_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p5]=/HZa-Zto2L-ato2G_Par-M-0p5_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p6]=/HZa-Zto2L-ato2G_Par-M-0p6_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p7]=/HZa-Zto2L-ato2G_Par-M-0p7_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p8]=/HZa-Zto2L-ato2G_Par-M-0p8_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M0p9]=/HZa-Zto2L-ato2G_Par-M-0p9_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
DATASET[M1]=/HZa-Zto2L-ato2G_Par-M-1_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM

mkdir -p ${DUMP_DIR}

for tag in M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9 M1; do
    ds=${DATASET[$tag]}
    echo "===== ${tag}: ${ds} ====="
    f=$(dasgoclient -query "file dataset=${ds}" | head -1)
    if [ -z "${f}" ]; then
        echo "[skip] no DAS hits for ${ds}"
        continue
    fi
    out=${DUMP_DIR}/mA_${tag}_rechits_${N_EVENTS}ev.root
    if [ ! -f "${out%.root}_numEvent${N_EVENTS}.root" ]; then
        echo "[dump] cmsRun → ${out}"
        cmsRun ${CFG} \
            inputFiles=root://cms-xrd-global.cern.ch/${f} \
            outputFile=${out} \
            maxEvents=${N_EVENTS} 2>&1 | tail -3
    else
        echo "[skip dump] ${out%.root}_numEvent${N_EVENTS}.root exists"
    fi
done

# Plot stage runs in a fresh subshell so the CMSSW PyROOT / conda PyROOT
# environments do not collide.
for tag in M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9 M1; do
    dump=${DUMP_DIR}/mA_${tag}_rechits_${N_EVENTS}ev_numEvent${N_EVENTS}.root
    out_dir=${PLOT_BASE}/mA_${tag}
    if [ ! -f "${dump}" ]; then
        echo "[skip plot] missing dump ${dump}"
        continue
    fi
    echo "===== plot ${tag}: ${dump} → ${out_dir} ====="
    env -i HOME=${HOME} PATH=/usr/bin:/bin bash -lc "
        source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
        conda activate hza_ana
        python ${PLOTTER} ${dump} --mode event -n ${N_EVENTS} --out-dir ${out_dir}
    "
done

echo "===== all done ====="
