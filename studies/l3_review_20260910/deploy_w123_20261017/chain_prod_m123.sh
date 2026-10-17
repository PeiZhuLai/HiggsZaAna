#!/usr/bin/env bash
# Step C+D: deploy the low-mass BDT low_relE2_w123 (relE2 inputs, m_a 1-3 all trained inside 115-133 GeV) in production
# for m_a = 1, 2, 3 (high-mass model and m_a >= 4 untouched). Same recipe as the validated sandbox chain
# varsearch_validate/chainlow_relE2w123 (chainlow.sh), now on the
# production code and paths; each stage is compared with the sandbox result it must reproduce.
#   0  wait for the patched scored files (run_patch.sh: the 273 files of data, DY, mA_M1-3), fingerprint production
#   1  working points mA1-3 in Plot/output/MVAcut_points_run3.json: 0.98115 0.98837 0.99226 (x0.5 rule)
#   2  S1 lnN for mA1-3 (eval_alt_reweight S1_RW_TAG=fsrfix_dyveto + make_s1_lnn_json --reference nominal), merged
#      into the production default JSON; must equal the sandbox values, mA4-30 unchanged
#   3  MVA cut data + signal mA1-3 (EOS root_MVAcut backed up first); entries equal to the sandbox chain
#   4  Tree2WS mA1-3; data workspace counts equal to the sandbox
#   5  background envelopes mA1-3 (production binaries, fixed window 95-180); envelope + GOF equal to the sandbox
#   6  signal model mA1-3 (sig_sys 1-4)
#   7  datacards mA1-3 (wilson + xsec); identical to the sandbox cards up to the code/EOS path prefixes
#   8  text2workspace + expected limits (both cards; equal to the sandbox), observed limits (both), significance
#   9  interpolation-efficiency JSON (signal_eff_sumw.py): only the mA1-3 entries may change
#  10  bias: the sandbox high-stat results (same workspaces, 10 x 1000 toys, xsec card) copied into Bias_nominal
#  11  impacts mA1-3, both cards (condor; stuck jobs are re-run locally)
#  12  regression: every changed production file belongs to mA1-3; limits of mA4-30 unchanged
# Backups: AFS products -> backup_prod/<tag>.tar (whole directories), EOS root_MVAcut mA1-3 ->
#          /eos/home-p/pelai/HZa/root_MVAcut_before_w123_20261017, bias mA1-3 -> bias_outputs_highstat_variants/
# Markers: logs/chain.DONE or logs/chain.FAILED; heartbeat logs/chain_heartbeat.txt
# At most 4 local processes (the rest of the score patch runs beside it with 2; lxplus limit 6).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
SBX=$CH/flashggFinalFit
SMV=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
DC=$F/Datacard/output_Datacard_leptons
MV=/eos/home-p/pelai/HZa/root_MVAcut
MVB=/eos/home-p/pelai/HZa/root_MVAcut_before_w123_20261017
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
WPJ=$PLOT/output/MVAcut_points_run3.json
S1D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927
LNN=$S1D/results/s1_lnn_scheme2_sigrw_dyveto.json
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
TAG=before_w123_20261017
L=$D/logs; DR=$D/drivers; BK=$D/backup_prod
START=${START_STAGE:-0}
mkdir -p $L $DR $BK
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/chain_heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain.FAILED; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
rm -f $L/chain.DONE $L/chain.FAILED
AFSDIRS="Signal/outdir_ele Signal/outdir_mu Datacard/yields_leptons Datacard/output_Datacard_leptons Combine/root_t2w_xsec
 Combine/root_t2w_wilson Combine/output_combine_results_xsec Combine/output_combine_results_wilson
 Combine/output_combine_results_observed Combine/output_combine_results_observed_wilson Combine/output_significance
 Combine/output_impacts Combine/output_impacts_xsec Background/ALP_BkgModel_ReReco/fit_results_run3"
manifest() { ( cd $F && find $AFSDIRS -type f ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 md5sum ); }

if run 0; then
hb "[0] wait for the 273 patched scored files"
while true; do
  n=0; while read -r s e; do [ -e $D/state/${s}__${e}.ok ] && n=$((n + 1)); done < $L/patch_first.txt
  [ $n -eq 273 ] && break
  [ -e $L/patch.FAILED ] && die "score patch failed: $(cat $L/patch.FAILED)"
  grep -q " FAIL " $L/patch_progress.txt && die "score patch has failed files (patch_progress.txt)"
  hb "    patched $n/273"; sleep 600
done
hb "    273/273 patched (data, DY, mA_M1-3)"
[ -s $L/manifest.before ] || manifest > $L/manifest.before
hb "    production manifest: $(wc -l < $L/manifest.before) files"
fi

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
export PYTHONPATH="$F/tools:$F/Signal/tools"
set -u

if run 1; then
hb "[1] working points mA1-3"
cp -n $WPJ $WPJ.$TAG
python3 - $WPJ <<'PY' | tee -a $L/chain_heartbeat.txt || die "working-point JSON update"
import json, sys
p = sys.argv[1]; d = json.load(open(p)); new = {1: 0.98115, 2: 0.98837, 3: 0.99226}
for e in d["results"]:
    m = int(e["mA"])
    if m in new:
        print("    mA%d %.5f -> %.5f" % (m, float(e["MVAcut"]), new[m]))
        e["MVAcut"] = new[m]; e["Significance"] = None
        e["note"] = "low_relE2_w123: x0.5 rule of the closure scan (closure + data GOF pass), not the optimization"
json.dump(d, open(p, "w"), indent=1)
r = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(p))["results"]}
assert [r[m] for m in (1, 2, 3)] == [0.98115, 0.98837, 0.99226] and len(r) == 14   # anchors only
PY
python3 - $WPJ $WPJ.$TAG <<'PY' || die "working points of mA4-30 changed"
import json, sys
a, b = ({int(e["mA"]): e for e in json.load(open(p))["results"]} for p in sys.argv[1:3])
assert sorted(a) == sorted(b) and all(a[m] == b[m] for m in b if m > 3)
PY
fi

if run 2; then
hb "[2] S1 lnN mA1-3"
CSV=$S1D/results/s1_altsideband_yields_dyveto_w123_m123.csv
( cd $S1D && env -i HOME=$HOME PATH=/usr/bin:/bin S1_RW_TAG=fsrfix_dyveto $ALP_PY eval_alt_reweight.py --masses 1,2,3 --out $CSV ) > $L/s1_eval.log 2>&1 \
  || die "eval_alt_reweight failed (s1_eval.log)"
grep -q "variants: \['nominal', 'zmhsb', 'zmnarrowhsb'\]" $L/s1_eval.log || die "S1 eval did not load the three dyveto reweights"
SRC=$($ALP_PY -c "import json; print(' '.join(json.load(open('$LNN'))['source']))")
NEWJ=$S1D/results/s1_lnn_scheme2_sigrw_dyveto_w123.json
( cd $S1D/results && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY ../make_s1_lnn_json.py --reference nominal $NEWJ $SRC $CSV ) > $L/s1_json.log 2>&1 || die "make_s1_lnn_json failed"
$ALP_PY - $NEWJ $LNN $CH/s1/s1_lnn_scheme2_chain4.json <<'PY' | tee -a $L/chain_heartbeat.txt || die "S1 lnN differs from the sandbox or mA4-30 changed"
import json, sys
new, old, sbx = (json.load(open(p))["values"] for p in sys.argv[1:4])
bad = [m for m in old if int(m) > 3 and new.get(m) != old[m]] + ([] if sorted(new) == sorted(old) else ["keys"])
bad += [m for m in "123" if new[m] != sbx[m]]
for m in "123": print("    mA%s %s (was %s)" % (m, new[m], old[m]))
sys.exit(1 if bad else 0)
PY
cp -n $LNN $LNN.$TAG
cp $NEWJ $LNN
fi
export HZA_MVA_REWEIGHT_LNN_JSON=$LNN

if run 3; then
hb "[3] MVA cut mA1-3 (backup to $MVB)"
# file by file (xrdcp -r on EOS fails on the directory creation); destination checksum checked
( cd $MV && find data/mA_M1 data/mA_M2 data/mA_M3 sig/mA_M1 sig/mA_M2 sig/mA_M3 -type f ) > $L/mvacut_backup_list.txt
while read -r f; do
  x=/eos/user/p/pelai/HZa/root_MVAcut/$f; y=/eos/user/p/pelai/HZa/root_MVAcut_before_w123_20261017/$f
  ck=$(xrdfs eosuser.cern.ch query checksum $x | awk '{print $2}')
  [ "$(xrdfs eosuser.cern.ch query checksum $y 2>/dev/null | awk '{print $2}')" = "$ck" ] && continue
  xrdfs eosuser.cern.ch mkdir -p $(dirname $y)
  xrdcp -f -s --cksum adler32:$ck root://eosuser.cern.ch/$x root://eosuser.cern.ch/$y || die "backup of $f"
done < $L/mvacut_backup_list.txt
n0=$(find $MV/data/mA_M{1,2,3} $MV/sig/mA_M{1,2,3} -type f | wc -l); n1=$(find $MVB/data/mA_M{1,2,3} $MVB/sig/mA_M{1,2,3} -type f | wc -l)
[ $n0 -eq $n1 ] || die "root_MVAcut backup has $n1 files, production $n0"
hb "    backup: $n1 files"
sed 's/^mAs = \[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30\]$/mAs = [1,2,3]/' \
  $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py > $DR/apply_bdt_data_m123.py
grep -q '^mAs = \[1,2,3\]$' $DR/apply_bdt_data_m123.py || die "apply_bdt_data mass list not reduced"
T=$(date +%s)
( ok=0; for a in 1 2 3; do
    ( timeout 5400 python3 $DR/apply_bdt_data_m123.py 2>&1; echo "CHAIN_RC=$?" ) | head -c 50000000 > $L/mva_data_try$a.log
    grep -q '^CHAIN_RC=0$' $L/mva_data_try$a.log && { ok=1; break; }; sleep 60; done; [ $ok = 1 ] ) &
pd=$!
sig_one() { local s=$1 y=$2 a; for a in 1 2 3; do
    python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples $s --years $y > $L/mva_sig_${s}_$y.log 2>&1 && return 0; sleep 60; done; return 1; }
export -f sig_one; export F L
for m in 1 2 3; do for y in $YEARS; do echo mA_M$m $y; done; done | xargs -P3 -L1 bash -c 'sig_one $0 $1 || echo "FAILED $0 $1"' > $L/mva_sig_status.txt
wait $pd || die "apply_bdt_data failed 3 times (mva_data_try*.log)"
grep -q FAILED $L/mva_sig_status.txt && die "apply_bdt_sig failed: $(grep FAILED $L/mva_sig_status.txt | tr '\n' ' ')"
for m in 1 2 3; do fresh $MV/data/mA_M$m/run3.root $T || die "data mA$m stale"; for y in $YEARS; do fresh $MV/sig/mA_M$m/output_$y.root $T || die "sig mA$m $y stale"; done; done
python3 - $MV $SMV/root_MVAcut $LNN <<'PY' | tee $L/mva_compare.txt || die "MVA cut outputs differ from the sandbox chain"
import json, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
P, S, J = sys.argv[1], sys.argv[2], json.load(open(sys.argv[3]))["values"]; bad = 0
def trees(path):
    out = {}
    f = ROOT.TFile.Open(path)
    def walk(d, pre):
        for k in d.GetListOfKeys():
            o = k.ReadObj()
            if o.InheritsFrom("TDirectory"): walk(o, pre + o.GetName() + "/")
            elif o.InheritsFrom("TTree"):
                n = o.GetEntries(); sw = 0.0
                if o.GetBranch("weight"):
                    for e in o: sw += e.weight
                out[pre + o.GetName()] = (n, round(sw, 6))
    walk(f, ""); f.Close(); return out
for m in (1, 2, 3):
    for rel in ["data/mA_M%d/run3.root" % m] + ["sig/mA_M%d/output_%s.root" % (m, y) for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")]:
        a, b = trees("%s/%s" % (P, rel)), trees("%s/%s" % (S, rel))
        ok = a == b
        if rel.startswith("sig"):
            f = ROOT.TFile.Open("%s/%s" % (P, rel)); t = f.Get("DiphotonTree/ggh_125_Za_ele_13p6TeV_cat0")
            if t and t.GetEntries():
                t.GetEntry(0); ok = ok and abs(t.weight_mva_reweight_Up - J[str(m)]["up"]) < 1e-6 and abs(t.weight_mva_reweight_Down - J[str(m)]["down"]) < 1e-6
            ok = ok and len(a) == 34
            f.Close()
        print("  %-36s %d trees %s" % (rel, len(a), "= sandbox" if ok else "DIFFERS"))
        bad += not ok
sys.exit(1 if bad else 0)
PY
fi

if run 4; then
hb "[4] Tree2WS mA1-3"
sed -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(1 2 3)/' \
    -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(1 2 3)/' \
    -e 's/^max_parallel=6$/max_parallel=4/' \
    $F/Trees2WS/run_tree2ws.sh > $DR/run_tree2ws_m123.sh
grep -q '^max_parallel=4$' $DR/run_tree2ws_m123.sh || die "tree2ws parallelism not reduced"
grep -q '^mAs_sig=(1 2 3)$' $DR/run_tree2ws_m123.sh && grep -q '^mAs_data=(1 2 3)$' $DR/run_tree2ws_m123.sh || die "tree2ws mass lists not reduced"
T=$(date +%s)
( cd $F/Trees2WS && bash $DR/run_tree2ws_m123.sh ) > $L/tree2ws.log 2>&1
for m in 1 2 3; do
  fresh $MV/data/mA_M$m/ws/run3.root $T || die "data ws mA$m stale"
  n=$(find $MV/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 10 ] || die "signal ws mA$m: $n fresh files"
done
python3 - $MV $SMV/root_MVAcut <<'PY' | tee -a $L/chain_heartbeat.txt || die "data workspaces differ from the sandbox"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
bad = 0
for m in (1, 2, 3):
    n = []
    for b in sys.argv[1:3]:
        f = ROOT.TFile.Open("%s/data/mA_M%d/ws/run3.root" % (b, m)); n.append(f.Get("CMS_hza_workspace").data("Data_13p6TeV").sumEntries()); f.Close()
    print("    mA%d data ws %d (sandbox %d)" % (m, n[0], n[1])); bad += n[0] != n[1]
sys.exit(bad)
PY
fi

if run 5; then
hb "[5] background envelopes mA1-3"
[ -s $BK/fit_results_run3_m123.tar ] || tar -C $FR -cf $BK/fit_results_run3_m123.tar 1 2 3 AllFitResults_blind AllFitResults_unblind || die "backup fit_results_run3"
T=$(date +%s)
bkg_one() {
  local m=$1 P=$FR/$1
  rm -rf $P; mkdir -p $P/HZAmassInde_fTest
  ( cd $F/Background && ./bin/fTest_ALP_turnOn -i $MV/data/mA_M$m/ws/run3.root \
      --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root -D $P/HZAmassInde_fTest \
      --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  for b in blind unblind; do
    extra=""; [ $b = unblind ] && extra="--unblind"
    ( cd $F/Background && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_$b \
        --total_OutDir $FR/AllFitResults_$b -o $P/BkgPlots_$b.root --sqrts 13p6TeV --isMultiPdf --useBinnedData \
        --massStep 2.5 --mhVal 125.0 --maVal $m --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
        --intLumi 172.13 -c 0 --isFlashgg 0 --doBands $extra ) >> $P/makeBkgPlots.stdout.log 2>&1
  done
  echo "    [$(date '+%H:%M:%S')] bkg mA$m done"
}
export -f bkg_one; export F FR MV
printf '1\n2\n3\n' | xargs -P3 -I{} bash -c 'bkg_one {}' | tee -a $L/chain_heartbeat.txt
python3 - $FR $SBX/Background/ALP_BkgModel_ReReco/fit_results_run3 $T <<'PY' | tee $L/bkg_compare.txt || die "background envelopes differ from the sandbox"
import os, re, sys
P, S, T = sys.argv[1], sys.argv[2], int(sys.argv[3]); bad = 0
for m in (1, 2, 3):
    env = []
    for b in (P, S):
        txt = open("%s/%d/HZAmassInde_fTest/EnvelopeResults.txt" % (b, m)).read()
        env.append(re.findall(r"pdf : ([A-Za-z0-9]+) , gof : ([0-9.e-]+)", txt))
    ok = env[0] == env[1] and os.path.getmtime("%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (P, m)) >= T
    print("    mA%d %s %s" % (m, " ".join("%s:%s" % e for e in env[0]), "= sandbox" if ok else "DIFFERS from %s" % env[1]))
    bad += not ok
sys.exit(bad)
PY
fi

if run 6; then
hb "[6] signal model mA1-3"
[ -s $BK/signal_outdirs.tar ] || tar -C $F/Signal -cf $BK/signal_outdirs.tar outdir_ele outdir_mu || die "backup signal outdirs"
T=$(date +%s)
mkdir -p $DR/sig_single
sigm() {
  local m=$1 s
  for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
    # the RunPlotter script also replots the 0pX (HZa_merged, 2024) models in a second loop: cut it (mAs_lowMA=())
    sed -e "s/^mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )$/mAs=( $m )/" -e "s/^mAs_lowMA=( 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 )$/mAs_lowMA=()/" \
        $F/shellScripts/sig_sys/$s.sh > $DR/sig_single/${s}_m$m.sh
    grep -q "^mAs=( $m )$" $DR/sig_single/${s}_m$m.sh || { echo "  copy of $s for mA$m failed"; return 1; }
    bash $DR/sig_single/${s}_m$m.sh > $L/sig_${s}_m$m.log 2>&1
  done
  echo "    [$(date '+%H:%M:%S')] signal mA$m done"
}
export -f sigm; export F DR L
printf '1\n2\n3\n' | xargs -P3 -I{} bash -c 'sigm {}' | tee -a $L/chain_heartbeat.txt
for m in 1 2 3; do for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || die "sigfit mA$m $y $ch stale (sig_*_m$m.log)"
done; done; done
hb "    signal fits 30/30 fresh"
fi

if run 7; then
hb "[7] datacards mA1-3"
[ -s $BK/datacards.tar ] || tar -C $F/Datacard -cf $BK/datacards.tar yields_leptons output_Datacard_leptons || die "backup datacards"
for f in 1_runjob_gen_datacard_makeYields.sh 2_runjob_gen_datacard_makeDatacard.sh; do
  sed 's/^mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )$/mAs=( 1 2 3 )/' $F/Datacard/$f > $DR/${f%.sh}_m123.sh
  grep -q '^mAs=( 1 2 3 )$' $DR/${f%.sh}_m123.sh || die "$f mass list not reduced"
done
T=$(date +%s)
( cd $F/Datacard && sh $DR/1_runjob_gen_datacard_makeYields_m123.sh ) > $L/dc_1.log 2>&1
( cd $F/Datacard && sh $DR/2_runjob_gen_datacard_makeDatacard_m123.sh ) > $L/dc_2.log 2>&1
for m in 1 2 3; do for s in "" _xsec; do
  fresh $DC/${m}_pruned_datacard_leptons$s.txt $T || die "datacard mA$m$s stale (dc_*.log)"
  sed -e "s#$F#FFF#g" -e "s#$MV#MVMV#g" $DC/${m}_pruned_datacard_leptons$s.txt > $L/dc_prod_${m}$s.txt
  sed -e "s#$SBX#FFF#g" -e "s#$SMV/root_MVAcut#MVMV#g" $SBX/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$s.txt > $L/dc_sbx_${m}$s.txt
  diff -q $L/dc_prod_${m}$s.txt $L/dc_sbx_${m}$s.txt > /dev/null || die "datacard mA$m$s differs from the sandbox (diff $L/dc_prod_${m}$s.txt $L/dc_sbx_${m}$s.txt)"
done; done
hb "    6 datacards identical to the sandbox up to the path prefixes"
fi

lim_one() {   # same definition as inside stage 8; here so that a restart at stage 9 has it
  local m=$1 k=$2 sfx="" W
  [ $k = xsec ] && sfx=_xsec
  W=$D/work_combine/${k}_$m; rm -rf $W; mkdir -p $W
  ( cd $C && python3 RunText2Workspace.py --batch local --queue workday --mA $m --card $k ) > $L/t2w_${k}_m$m.log 2>&1
  ( cd $W && combine -M AsymptoticLimits $C/root_t2w_$k/${m}_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_${k}_m$m.log 2>&1
    mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_$k/
    combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$sfx.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/observed_${k}_m$m.log 2>&1
    if [ $k = xsec ]; then mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed/
    else mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed_wilson/; fi
    if [ $k = xsec ]; then
      combine -M Significance $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 -n $m > $L/significance_m$m.log 2>&1
      mv higgsCombine$m.Significance.mH125.38.root $C/output_significance/
    fi )
}
export -f lim_one; export C D L F

if run 8; then
hb "[8] text2workspace + limits mA1-3"
[ -s $BK/combine.tar ] || tar -C $C -cf $BK/combine.tar root_t2w_xsec root_t2w_wilson output_combine_results_xsec output_combine_results_wilson \
  output_combine_results_observed output_combine_results_observed_wilson output_significance || die "backup combine outputs"
T=$(date +%s)
lim_one() {
  local m=$1 k=$2 sfx="" W
  [ $k = xsec ] && sfx=_xsec
  W=$D/work_combine/${k}_$m; rm -rf $W; mkdir -p $W
  ( cd $C && python3 RunText2Workspace.py --batch local --queue workday --mA $m --card $k ) > $L/t2w_${k}_m$m.log 2>&1
  ( cd $W && combine -M AsymptoticLimits $C/root_t2w_$k/${m}_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_${k}_m$m.log 2>&1
    mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_$k/
    combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$sfx.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/observed_${k}_m$m.log 2>&1
    if [ $k = xsec ]; then mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed/
    else mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed_wilson/; fi
    if [ $k = xsec ]; then
      combine -M Significance $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 -n $m > $L/significance_m$m.log 2>&1
      mv higgsCombine$m.Significance.mH125.38.root $C/output_significance/
    fi )
}
export -f lim_one; export C D L F
for m in 1 2 3; do for k in xsec wilson; do echo "$m $k"; done; done | xargs -P4 -L1 bash -c 'lim_one $0 $1'
for m in 1 2 3; do
  for k in xsec wilson; do fresh $C/root_t2w_$k/${m}_Datacard_leptons.root $T && fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "expected limit $k mA$m"; done
  fresh $C/output_combine_results_observed/higgsCombine$m.AsymptoticLimits.mH125.38.root $T && fresh $C/output_combine_results_observed_wilson/higgsCombine$m.AsymptoticLimits.mH125.38.root $T \
    && fresh $C/output_significance/higgsCombine$m.Significance.mH125.38.root $T || die "observed/significance mA$m"
done
python3 - $C $SBX/Combine $BK <<'PY' | tee $L/limits_compare.txt || die "expected limits differ from the sandbox"
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def q(path):
    f = ROOT.TFile.Open(path); out = {}
    for e in f.Get("limit"): out[round(e.quantileExpected, 3)] = e.limit
    f.Close(); return out
bad = 0
print("    expected limit on r (median [-1s, +1s]); sandbox = chainlow_relE2w123")
for k in ("xsec", "wilson"):
    for m in (1, 2, 3):
        p = q("%s/output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (sys.argv[1], k, m))
        s = q("%s/output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (sys.argv[2], k, m))
        same = all(abs(p[x] / s[x] - 1) < 1e-6 for x in (0.025, 0.16, 0.5, 0.84, 0.975))
        print("    %-6s mA%d %.4f [%.4f, %.4f]  sandbox %.4f  %s" % (k, m, p[0.5], p[0.16], p[0.84], s[0.5], "same" if same else "DIFFERS"))
        bad += not same
sys.exit(bad)
PY
fi

if run 9; then
hb "[9] interpolation-efficiency JSON"
for f in $PLOT/output/sigEfficiencyVmA_*_quadratic_interp_ma_points.json; do cp -n $f $f.$TAG; done
T=$(date +%s)
( cd $PLOT/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin PYTHONPATH=$PLOT/lib $ALP_PY signal_eff_sumw.py ) > $L/signal_eff_sumw.log 2>&1 || die "signal_eff_sumw.py failed"
python3 - $PLOT/output $TAG <<'PY' | tee -a $L/chain_heartbeat.txt || die "interpolation JSON: entries other than mA1-3 changed"
import glob, json, sys
bad = 0
for p in sorted(glob.glob("%s/sigEfficiencyVmA_*_quadratic_interp_ma_points.json" % sys.argv[1])):
    n, o = json.load(open(p))["values"], json.load(open(p + "." + sys.argv[2]))["values"]
    moved = sorted({int(m) for y in n for m in n[y] if abs(float(n[y][m]) - float(o[y][m])) > 1e-12 * max(1.0, abs(float(o[y][m])))})
    rel = max([abs(float(n[y][m]) / float(o[y][m]) - 1) for y in n for m in n[y] if int(m) > 3] or [0.0])
    print("    %s: changed mA %s; largest relative change at mA>=4 %.2e (global quadratic spline through the 14 anchors)"
          % (p.split("/")[-1], moved, rel))
    bad += rel > 1e-3
sys.exit(bad)
PY

INTERP="11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29"
hb "[9b] datacards + limits of the 16 interpolated masses (their efficiency rateParam comes from the JSON)"
for f in 1_runjob_gen_datacard_makeYields.sh 2_runjob_gen_datacard_makeDatacard.sh; do
  sed "s/^mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )$/mAs=( $INTERP )/" $F/Datacard/$f > $DR/${f%.sh}_interp.sh
  grep -q "^mAs=( $INTERP )$" $DR/${f%.sh}_interp.sh || die "$f interpolated mass list"
done
T=$(date +%s)
( cd $F/Datacard && sh $DR/1_runjob_gen_datacard_makeYields_interp.sh ) > $L/dc_interp_1.log 2>&1
( cd $F/Datacard && sh $DR/2_runjob_gen_datacard_makeDatacard_interp.sh ) > $L/dc_interp_2.log 2>&1
for m in $INTERP; do
  for s in "" _xsec; do fresh $DC/${m}_pruned_datacard_leptons$s.txt $T || die "interpolated datacard mA$m$s stale"; done
  [ $(grep -c "^interpolate_yield_.*rateParam" $DC/${m}_pruned_datacard_leptons.txt) -eq 10 ] || die "mA$m: interpolation rateParam missing"
done
for m in $INTERP; do for k in xsec wilson; do echo "$m $k"; done; done | xargs -P4 -L1 bash -c 'lim_one $0 $1'
for m in $INTERP; do for k in xsec wilson; do fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit $k mA$m"; done; done
mkdir -p $D/work_combine/old && tar -C $D/work_combine/old -xf $BK/combine.tar output_combine_results_xsec output_combine_results_wilson
python3 - $C $D/work_combine/old "$INTERP" <<'PY' | tee -a $L/limits_compare.txt || die "interpolated limits moved by more than 1e-3"
import sys, ROOT
ROOT.gROOT.SetBatch(True)
C, O, ms = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split()]
def q(f):
    f = ROOT.TFile.Open(f); out = {round(e.quantileExpected, 3): e.limit for e in f.Get("limit")}; f.Close(); return out
worst = 0
for k in ("xsec", "wilson"):
    for m in ms:
        name = "output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (k, m)
        o, n = q("%s/%s" % (O, name)), q("%s/%s" % (C, name))
        r = max(abs(n[x] / o[x] - 1) for x in (0.025, 0.16, 0.5, 0.84, 0.975)); worst = max(worst, r)
        print("    %-6s mA%-2d median %.5f -> %.5f (largest quantile change %.1e)" % (k, m, o[0.5], n[0.5], r))
print("    interpolated masses: largest relative limit change %.1e" % worst)
sys.exit(1 if worst > 1e-3 else 0)
PY
fi

if run 10; then
hb "[10] bias mA1-3: sandbox high-stat results (same workspaces) into Bias_nominal"
mkdir -p $BN/bias_outputs_highstat_variants
for m in 1 2 3; do
  [ -d $BN/bias_outputs_highstat_variants/mA_${m}_$TAG ] || mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG || die "move old bias mA$m"
  rm -rf $BN/bias_outputs_highstat/mA_$m
  cp -a $SMV/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat/mA_$m || die "copy bias mA$m"
  diff -rq $SMV/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat/mA_$m > /dev/null || die "bias copy mA$m incomplete"
done
cat > $BN/bias_outputs_highstat/README_w123_m123.txt <<EOF
mA_1, mA_2, mA_3: copied from the sandbox chain varsearch_validate/chainlow_relE2w123 (EOS
$SMV/bias_outputs_highstat), 10 chunks x 1000 toys, cross-section card. The production workspaces of
$C/root_t2w_xsec for mA1-3 come from datacards identical to the sandbox ones (up to the path prefixes) and give the
same expected limits (deploy_w123_20261017/logs/limits_compare.txt). mA_1 = envelope Exp1/Lau1/Pow1.
Previous low-mass model results: bias_outputs_highstat_variants/mA_<m>_$TAG.
EOF
( BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh ) > $L/bias_collect.log 2>&1 || hb "    collect_bias_results.sh rc!=0 (bias_collect.log)"
python3 - $BN/bias_outputs_highstat <<'PY' | tee -a $L/chain_heartbeat.txt
import json, sys
for m in (1, 2, 3):
    fr = json.load(open("%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (sys.argv[1], m, m)))["fit_results"]
    print("    mA%d bias " % m + " ".join("%s:%+.3f" % (k, v["mean"]) for k, v in sorted(fr.items(), key=lambda x: -abs(x[1]["mean"]))))
PY
fi

if run 11; then
hb "[11] impacts mA1-3 (wilson + xsec)"
[ -s $BK/impacts.tar ] || tar -C $C -cf $BK/impacts.tar output_impacts output_impacts_xsec || die "backup impacts"
T0=$(date +%s)
out=$( (IMPACT_CARD=wilson bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh 1 2 3; IMPACT_CARD=xsec bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh 1 2 3) 2>&1)
echo "$out" > $L/impact_submit.log
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$out" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "impact submission failed (impact_submit.log)"
hb "    impact clusters $CL"
echo "$CL" > $L/impact_clusters.txt
start=$(date +%s)
while true; do
  q=$(condor_q $CL -af ProcId JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(awk '$2==1' <<<"$q" | wc -l); runn=$(awk '$2==2' <<<"$q" | wc -l); held=$(awk '$2==5' <<<"$q" | wc -l)
  hb "    impacts idle=$idle running=$runn held=$held"
  [ $((idle + runn + held)) -eq 0 ] && break
  if [ $runn -eq 0 ] && [ $(( $(date +%s) - start )) -gt 7200 ]; then
    hb "    impact jobs idle/held for 2 h: remove them and run locally (one at a time)"
    condor_rm $CL > /dev/null 2>&1; sleep 30
    for k in wilson xsec; do for m in 1 2 3; do
      j=$C/$( [ $k = xsec ] && echo output_impacts_xsec || echo output_impacts )/${m}_impacts.json
      fresh $j $T0 && continue
      IMPACT_CARD=$k bash $F/shellScripts/impact/Condor/runjob_expectedImpact.sh $m $F ${CMSSW_SRC%/src} > $L/impact_local_${k}_m$m.log 2>&1 || hb "    local impact $k mA$m rc!=0"
    done; done
    break
  fi
  sleep 1800
done
python3 - $C $T0 <<'PY' | tee $L/impact_check.txt || die "impacts missing or incomplete"
import json, os, sys
C, T0 = sys.argv[1], float(sys.argv[2]); bad = 0
for k, d in (("wilson", "output_impacts"), ("xsec", "output_impacts_xsec")):
    for m in (1, 2, 3):
        p = "%s/%s/%d_impacts.json" % (C, d, m)
        if not (os.path.exists(p) and os.path.getmtime(p) >= T0): print("    %s mA%d MISSING/stale" % (k, m)); bad += 1; continue
        j = json.load(open(p)); r = [x for x in j["POIs"] if x["name"] == "r"][0]["fit"]
        print("    %-6s mA%d r = %.3f -%.3f/+%.3f, %d parameters" % (k, m, r[1], r[1] - r[0], r[2] - r[1], len(j["params"])))
sys.exit(bad)
PY
fi

if run 12; then
hb "[12] regression check"
manifest > $L/manifest.after
python3 - $L/manifest.before $L/manifest.after <<'PY' | tee $L/regression.txt
import re, sys
def rd(p): return {l[34:].strip(): l[:32] for l in open(p)}
a, b = rd(sys.argv[1]), rd(sys.argv[2])
changed = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
pat = re.compile(r"(^|/)(11|12|13|14|16|17|18|19|21|22|23|24|26|27|28|29)_(pruned|datacard|Datacard)|higgsCombine(11|12|13|14|16|17|18|19|21|22|23|24|26|27|28|29)\.|(^|/)(1|2|3)(_|/)|_(1|2|3)_(20|run3|ERA)|_(1|2|3)\.(pdf|png|root|txt|json|C)$|higgsCombine(1|2|3)\.|mA_?(1|2|3)(_|\.|/)|_m(1|2|3)(_|\.)")
other = [k for k in changed if not pat.search(k)]
print("    changed files: %d, not recognizably mA1-3 or an interpolated-mass datacard/limit: %d" % (len(changed), len(other)))
for k in other[:50]: print("      ", k)
sys.exit(1 if other else 0)
PY
[ ${PIPESTATUS[0]} -eq 0 ] || die "files outside mA1-3 changed (regression.txt)"
fi

hb "CHAIN DONE"
date '+%F %T' > $L/chain.DONE
echo "$(date '+%F %H:%M') DONE  deploy_w123 production chain mA1-3 -> $L/chain_heartbeat.txt" >> /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/logs/PENDING_RESULTS.txt
