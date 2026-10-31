#!/usr/bin/env bash
# =============================================================================
# prestage_stuck_inputs_phidzee_20261030.sh — 把反覆卡在 xrdcp 的輸入檔先拉到 eoscms，再讓 job 重排
# =============================================================================
# 為什麼：job 內 analysis.py 的 stage_and_open 用沒有 timeout 的 xrdcp 從 xrootd-cms.infn.it 拉輸入，
#   部分 worker 上會卡住不動（10-30 有 28 個卡 2-6 h，vacate 重排後 10 個又卡在同一個檔）；
#   同一個檔從 lxplus 拉只要數十秒。所以從 lxplus 拉到 eoscms，用 DAS 的 adler32 核對，
#   改寫該 job 的 config json（輸入 URL 換成 eoscms 副本，原檔備份 .bak_global），再 condor_vacate_job。
#   job 重啟時 stage_and_open 改從 eoscms 拉。
#
# 用法：bash prestage_stuck_inputs_phidzee_20261030.sh <cluster.proc> ... [--go]   （預設 dry-run）
# 輸入：condor job 的 JobSubmitFile 所在目錄 + Cmd 的 <task>/job_N → *config*.json
#       （parquet_DNA_tmp_phidzee 與 parquet_DNA_zee_zmmg_ctrl 都適用）
# 輸出：/eos/cms/store/group/phys_susy/pelai/HZa/prestage_phidzee/<LFN>；改寫的 config json（備份 .bak_global）
# log ：HiggsDNA/signal_rerun_logs/phidzee_20261029/prestage_20261030.log
set -uo pipefail
GO=0; IDS=()
for a in "$@"; do [ "$a" = "--go" ] && GO=1 || IDS+=("$a"); done
export X509_USER_PROXY=/tmp/x509up_u$(id -u)
DAS=/cvmfs/cms.cern.ch/common/dasgoclient
AFS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs
OUTB=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_phidzee
PS=/eos/cms/store/group/phys_susy/pelai/HZa/prestage_phidzee
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/signal_rerun_logs/phidzee_20261029/prestage_20261030.log
for id in "${IDS[@]}"; do
  if [ -d "$id" ]; then
    # 直接給 job 目錄（job 已不在 queue，例如重試用盡被 driver 放棄）：只 prestage + 改寫 config，不 vacate
    jd=${id%/}
  else
    read -r cmd subf <<< "$(condor_q "$id" -af Cmd JobSubmitFile 2>/dev/null)"
    [ -n "$cmd" ] || { echo "!! $id not in queue"; continue; }
    # job 目錄 = submit 檔所在的 unit 目錄 + Cmd 的 <task>/job_N（主分析 MC 與 control 都適用）
    jd=$(dirname "$subf")/$(basename "$(dirname "$(dirname "$cmd")")")/$(basename "$(dirname "$cmd")")
  fi
  cfg=$(ls "$jd"/*_config_job*.json 2>/dev/null | head -1)
  [ -f "$cfg" ] || { echo "!! $id no config json in $jd"; continue; }
  for u in $(grep -o 'root://[^"]*\.root' "$cfg" | sort -u); do
    lfn=/store/${u#*//store/}
    dest=$PS$lfn
    want=$($DAS -query="file file=$lfn" -json 2>/dev/null | grep -oE '"adler32":"[0-9a-f]+"' | head -1 | cut -d'"' -f4)
    [ -n "$want" ] || { echo "!! $id no DAS adler32 for $lfn"; continue; }
    echo "$id  $lfn  adler32=$want"
    [ "$GO" -eq 1 ] || continue
    if [ ! -s "$dest" ] || [ "$(xrdadler32 "$dest" 2>/dev/null | awk '{print $1}')" != "$want" ]; then
      mkdir -p "$(dirname "$dest")"
      okc=0
      for r in xrootd-cms.infn.it cms-xrd-global.cern.ch cmsxrootd.fnal.gov; do
        timeout 900 xrdcp -f --nopbar --cksum adler32:$want "root://$r/$lfn" "root://eoscms.cern.ch/$dest" >/dev/null 2>&1 && { okc=1; break; }
      done
      [ "$okc" -eq 1 ] || { echo "!! $id copy failed from all redirectors: $lfn" | tee -a "$LOG"; continue 2; }
    fi
    [ -f "$cfg.bak_global" ] || cp -p "$cfg" "$cfg.bak_global"
    python3 - "$cfg" "$u" "root://eoscms.cern.ch/$dest" <<'EOF'
import sys
fn, old, new = sys.argv[1:]
s = open(fn).read()
assert s.count(old) >= 1, "url not found"
open(fn, "w").write(s.replace(old, new))
EOF
    echo "$(date +%F_%T) $id staged $lfn -> $dest" >> "$LOG"
  done
  if [ "$GO" -eq 1 ]; then
    grep -q 'root://xrootd-cms.infn.it' "$cfg" && { echo "!! $id config still has infn URL"; continue; }
    [ -d "$id" ] && { echo "$(date +%F_%T) $id config rewritten (not in queue, no vacate)" | tee -a "$LOG"; continue; }
    condor_vacate_job "$id" >/dev/null 2>&1 && echo "$(date +%F_%T) $id vacated" | tee -a "$LOG"
  fi
done
[ "$GO" -eq 1 ] || echo "DRY-RUN (add --go)"
