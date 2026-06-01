#!/bin/bash
#
# Render ROI residual + m_Gamma plots for all available mA points in
#   /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all/
# (falls back to Sig_MC_MLNANO/ if the _all directory does not exist yet)
# Output goes to /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/roi_residual/
# as two separate PDF/PNG files per mass:
#   mA_<X>_mGamma.{png,pdf,root}
#   mA_<X>_residual.{png,pdf,root}
#
set -e

PLOTTER=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts/plot_roi_residual.py
OUT_DIR=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/roi_residual

PARQUET_BASE_A=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all
PARQUET_BASE_B=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO

declare -A GENMA
GENMA[M0p1]=0.1
GENMA[M0p2]=0.2
GENMA[M0p3]=0.3
GENMA[M0p4]=0.4
GENMA[M0p5]=0.5
GENMA[M0p6]=0.6
GENMA[M0p7]=0.7
GENMA[M0p8]=0.8
GENMA[M0p9]=0.9
GENMA[M1]=1.0

mkdir -p "${OUT_DIR}"

for tag in M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9 M1; do
    parquet=""
    for base in "${PARQUET_BASE_A}" "${PARQUET_BASE_B}"; do
        cand="${base}/mA_MLNANO_${tag}_2024/merged_nominal.parquet"
        if [ -f "${cand}" ]; then
            parquet="${cand}"
            break
        fi
    done
    if [ -z "${parquet}" ]; then
        echo "[skip ${tag}] parquet not found"
        continue
    fi
    gen_ma=${GENMA[$tag]}
    tag_label=$(echo "mA_${gen_ma}" | sed 's/\./p/')
    out_stub="${OUT_DIR}/${tag_label}"
    echo "===== ${tag} (gen ma=${gen_ma}) → ${out_stub}_{mGamma,residual} ====="
    env -i HOME=${HOME} PATH=/usr/bin:/bin bash -lc "
        source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
        conda activate hza_ana
        python ${PLOTTER} ${parquet} --gen-ma ${gen_ma} --output-stub ${out_stub} > /dev/null 2>&1
    "
done

echo "===== done ====="
