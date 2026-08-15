#!/bin/bash
# 2022/2023 friend parquet: submit=eos-home hza_ana(ROOT正常);worker=B-mode tarball(零eos home儲存);
# FPO=50 大幅降 job/cluster 數→condor_q 撐得住、監控正常(避免 FPO=2 的 schedd 死亡螺旋)。
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
FPO=50 CATALOG=metadata/samples/za_merged_data_2223_perds.json \
bash scripts/orchestrate_data_ml_production.sh
