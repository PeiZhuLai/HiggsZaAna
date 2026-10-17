#!/usr/bin/env bash
# Step B: put the low_swapPt_relE2 scores (MVA_Score_mA_M1..3) into every production scored file, in place.
# Per file (restartable; state markers in state/):
#   1. first time: download the production file, adler32 = server checksum, upload an untouched copy to the backup
#      tree (xrdcp --cksum checks the destination), write state/<s>__<e>.orig (the original checksum);
#   2. patch_scores.py on the backup copy (verified against .orig), with --ref = the sandbox chain copy when it exists
#      (m_a 1-3 signal + data/DY: the new scores must equal the ones the validated chain used);
#   3. upload to <file>.relE2partial (checksum checked), xrdfs mv over the production file (atomic overwrite), check
#      the production checksum = patched checksum, write state/<s>__<e>.ok.
# Order: the 273 files the m_a 1-3 chain reads first (data, DY, mA_M1-3 + systematic variations), then the rest.
# Input : /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix (1208 *.root)
# Backup: /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_before_relE2_20261016
# Output: the production files in place; logs/patch/<s>__<e>.log; logs/patch_progress.txt; logs/patch.DONE|FAILED
# Usage: run_patch.sh [smoke|first|rest|all] (default all; smoke = Data 2022preEE + mA_M1 2022preEE)  -- PATCH_P files at a time (default 6; lxplus front-end limit 6 processes in total)
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_relE2_20261016
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
XR=root://eosuser.cern.ch/
P=/eos/user/p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix
B=/eos/user/p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_before_relE2_20261016
R=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_low_swapPt_relE2
OLD=$D/model_set/previous_lowmass_run3.pkl
NEW=$D/model_set/model_Za_BDT_lowmass_run3.pkl
WHICH=${1:-all}
mkdir -p $D/state $D/logs/patch $D/work
[ -s $NEW ] && [ -s $OLD ] || { echo "model set missing"; exit 1; }
cmp -s $NEW /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/models/model_low_swapPt_relE2.pkl || { echo "NEW is not the relE2 model"; exit 1; }

ck() { xrdfs eosuser.cern.ch query checksum "$1" 2>/dev/null | awk '{print $2}'; }
one() {
  local s=$1 e=$2 k=${1}__${2} W L o r pck
  W=$D/work/$k; L=$D/logs/patch/$k.log
  [ -e $D/state/$k.ok ] && return 0
  mkdir -p $W
  {
  echo "[$(date '+%F %T')] $s $e"
  if [ ! -s $D/state/$k.orig ]; then
    o=$(ck $P/$s/$e.root); [ -n "$o" ] || { echo "no server checksum"; return 1; }
    xrdcp -f -s $XR$P/$s/$e.root $W/in.root || return 1
    [ "$(xrdadler32 $W/in.root | awk '{print $1}')" = "$o" ] || { echo "download checksum mismatch"; return 1; }
    xrdfs eosuser.cern.ch mkdir -p $B/$s
    xrdcp -f -s --cksum adler32:$o $W/in.root $XR$B/$s/$e.root || { echo "backup upload failed"; return 1; }
    echo $o > $D/state/$k.orig
  else
    o=$(cat $D/state/$k.orig)
    [ "$(ck $B/$s/$e.root)" = "$o" ] || { echo "backup checksum differs from the recorded original"; return 1; }
    xrdcp -f -s $XR$B/$s/$e.root $W/in.root || return 1
    [ "$(xrdadler32 $W/in.root | awk '{print $1}')" = "$o" ] || { echo "backup download checksum mismatch"; return 1; }
  fi
  ref=""
  if xrdfs eosuser.cern.ch stat $R/$s/$e.root >/dev/null 2>&1; then
    xrdcp -f -s $XR$R/$s/$e.root $W/ref.root || return 1; ref="--ref $W/ref.root"
  fi
  env -i HOME=$HOME PATH=/usr/bin:/bin $PY $D/patch_scores.py --src $W/in.root --out $W/out.root --old-low $OLD --new-low $NEW $ref || return 1
  pck=$(xrdadler32 $W/out.root | awk '{print $1}')
  xrdcp -f -s --cksum adler32:$pck $W/out.root $XR$P/$s/$e.root.relE2partial || { echo "upload failed"; return 1; }
  xrdfs eosuser.cern.ch mv $P/$s/$e.root.relE2partial $P/$s/$e.root || { echo "mv failed"; return 1; }
  [ "$(ck $P/$s/$e.root)" = "$pck" ] || { echo "production checksum after mv differs"; return 1; }
  echo "$o $pck ${ref:+ref}" > $D/state/$k.ok
  echo "[$(date '+%F %T')] OK orig=$o patched=$pck ${ref:+(equal to the sandbox chain scores)}"
  } > $L 2>&1
  local rc=$?
  rm -rf $W
  return $rc
}
try() {   # 3 attempts per file
  local a
  for a in 1 2 3; do one $1 $2 && { echo "$(date '+%F %T') OK   $1 $2" >> $D/logs/patch_progress.txt; return 0; }
    grep -q "GATE\|READBACK" $D/logs/patch/${1}__${2}.log && break; sleep 30; done
  echo "$(date '+%F %T') FAIL $1 $2 (logs/patch/${1}__${2}.log)" >> $D/logs/patch_progress.txt; return 1
}
export -f one try ck; export D PY XR P B R OLD NEW

xrdfs eosuser.cern.ch ls -R $P 2>/dev/null | grep '\.root$' | sed "s#^$P/##" | awk -F/ 'NF==2{sub(/\.root$/,"",$2); print $1, $2}' | sort > $D/logs/patch_all.txt
n=$(wc -l < $D/logs/patch_all.txt); [ "$n" -eq 1208 ] || { echo "$(date '+%F %T') expected 1208 files, listed $n" > $D/logs/patch.FAILED; exit 1; }
grep -E '^(Data|DY[A-Za-z0-9_]*|mA_M[123]|mA_M[123]_[A-Za-z_]+) ' $D/logs/patch_all.txt > $D/logs/patch_first.txt
# rest: nominal signal dirs first (dataVmc reads them), then the systematic-variation dirs
grep -vxFf $D/logs/patch_first.txt $D/logs/patch_all.txt | grep -E '^mA_M[0-9]+ ' > $D/logs/patch_rest.txt
grep -vxFf $D/logs/patch_first.txt $D/logs/patch_all.txt | grep -vE '^mA_M[0-9]+ ' >> $D/logs/patch_rest.txt
[ $(wc -l < $D/logs/patch_first.txt) -eq 273 ] || { echo "$(date '+%F %T') first batch is not 273 files" > $D/logs/patch.FAILED; exit 1; }
grep -E '^(Data 2022preEE|mA_M1 2022preEE)$' $D/logs/patch_first.txt > $D/logs/patch_smoke.txt
case $WHICH in smoke) L=$D/logs/patch_smoke.txt;; first) L=$D/logs/patch_first.txt;; rest) L=$D/logs/patch_rest.txt;; *) L=$D/logs/patch_first.txt;; esac
rm -f $D/logs/patch.DONE $D/logs/patch.FAILED
xargs -P${PATCH_P:-6} -L1 bash -c 'try $0 $1' < $L
[ $WHICH = all ] && xargs -P${PATCH_P:-6} -L1 bash -c 'try $0 $1' < $D/logs/patch_rest.txt
nok=$(ls $D/state/*.ok 2>/dev/null | wc -l)
want=1208; [ $WHICH = first ] && want=273; [ $WHICH = smoke ] && want=2
if [ "$nok" -ge "$want" ]; then echo "$(date '+%F %T') $nok files patched" > $D/logs/patch.DONE
else echo "$(date '+%F %T') only $nok/$want files patched (logs/patch_progress.txt)" > $D/logs/patch.FAILED; fi
