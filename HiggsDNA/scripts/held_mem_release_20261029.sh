#!/usr/bin/env bash
# HZa phidzee 生產：記憶體超限被 held 的 job，RequestMemory 依 3000 -> 6000 -> 9000 調高後 release（qedit 不觸發湊核，仍是 1 核）。
# 到 9000 仍 held 就不再動，印出來給人看。只處理輸出在 parquet_DNA_tmp_phidzee / parquet_DNA_zee_zmmg_ctrl 的 job。
C='Owner=="pelai" && JobStatus==5 && regexp("parquet_DNA_tmp_phidzee|parquet_DNA_zee_zmmg_ctrl",JobSubmitFile) && regexp("memory",HoldReason)'
condor_q -const "$C" -af ClusterId ProcId RequestMemory 2>/dev/null | while read c p m; do
  if [ "$m" -lt 6000 ]; then n=6000; elif [ "$m" -lt 9000 ]; then n=9000; else echo "STUCK $c.$p at ${m}MB"; continue; fi
  condor_qedit $c.$p RequestMemory $n >/dev/null && condor_release $c.$p >/dev/null && echo "released $c.$p ${m}->${n}MB"
done
