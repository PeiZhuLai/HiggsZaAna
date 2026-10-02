#!/usr/bin/env bash
# 2024 DY+jets overlap-veto retrain, user decisions 2026-10-01:
#  * mA3 working point 0.992 -> 0.988: after the retrain the closure at 0.992 is Z = 1.59 (bootstrap
#    40 replicas: Z 1.62 +- 1.17, 42% above 2); 0.988 is the most stable (bootstrap Z 1.03 +- 1.07,
#    15% above 2; data-envelope GOF >= 0.92; R 1.16, lowest).
#  * envelope per AN Sec 7.4 after the high-stat bias (> 0.2): mA25 Bernstein capped 6 -> 5
#    (Bern6 +0.212), Laurent removed at mA24 (Lau1 -0.232) and mA25 (Lau1 -0.215) -- Lau1 is the
#    lowest order, as for mA8 Pow1 on 2026-09-25.
# Adapted from chain_wp_mA3_0992.sh: mA3 MVA cut -> ws -> bkg fit; bkg refit of mA24/25 with the
# rebuilt fTest; signal model mA3; all datacards + limits; bias + impact for every mass whose limit
# moved (must include 3, 24, 25); then wait, merge and summarize those biases.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_dyveto_20260930
BK=$F/Background
SRCF=$BK/test/fTest_ALP_turnOn.cpp
HB=$L/dyveto_wp3_env2425_heartbeat.txt
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
NEWCUT=0.988
TAG=before_dyveto_wp0988_env2425_$(date +%Y%m%d_%H%M)
S1=$D/wp3_single; mkdir -p $S1
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

echo "=============== [1] backups ($TAG) ==============="
V=$C/output_combine_results_variants; mkdir -p $V
# 1_makeLimits.sh does rm -fr output_combine_results: move every non-standard file out first
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
cp -p $PLOT/output/MVAcut_points_run3.json $PLOT/output/MVAcut_points_run3.json.$TAG
cp -p $PLOT/scripts/scan_score_R_significance.py $PLOT/scripts/scan_score_R_significance.py.$TAG
for j in $PLOT/output/sigEfficiencyVmA_{ele,muon}_byYear_5years_quadratic_interp_ma_points.json; do cp -p $j $j.$TAG; done
for d in data sig; do [ -d $EOSMVA/$d/mA_M3 ] && cp -rp $EOSMVA/$d/mA_M3 $EOSMVA/${d}_variants_mA_M3_$TAG; done
# variants live OUTSIDE ALP_BkgModel_ReReco: sync_figures.sh copies that whole tree into the AN
mkdir -p $F/Background/fit_results_variants && for m in 3 24 25; do cp -rp $FR/$m $F/Background/fit_results_variants/${m}_$TAG; done
cp -p $SRCF $SRCF.$TAG
mkdir -p $F/Signal/sigfit_variants_$TAG && cp -p $F/Signal/outdir_{ele,mu}/signalFit/output/3_* $F/Signal/sigfit_variants_$TAG/ 2>/dev/null
cp -rp $F/Datacard/output_Datacard_leptons $F/Datacard/output_Datacard_leptons_$TAG
cp -p $FZ/3_Datacard_leptons.root $FZ/3_Datacard_leptons.root.$TAG
mkdir -p $BN/bias_outputs_highstat_variants && [ -d $BN/bias_outputs_highstat/mA_3 ] && mv $BN/bias_outputs_highstat/mA_3 $BN/bias_outputs_highstat_variants/mA_3_$TAG
cp -p $C/output_impacts/3_impacts.json $C/output_impacts/3_impacts.json.$TAG 2>/dev/null
echo "  done"

echo "=============== [2] working point: FIXED_WP[3] -> $NEWCUT, rewrite JSON ==============="
python3 - $PLOT/scripts/scan_score_R_significance.py $NEWCUT <<'PY' || exit 1
import sys
p, cut = sys.argv[1], sys.argv[2]; s = open(p).read()
old = "FIXED_WP = {2: 0.975, 3: 0.992}"
new = ("FIXED_WP = {2: 0.975, 3: %s}   # mA3 0.992 -> %s (2026-10-01): after the 2024 DY+jets overlap-veto\n"
       "                                  # retrain the closure at 0.992 was Z = 1.59 (bootstrap 1.62 +- 1.17); %s is\n"
       "                                  # the most stable cut (bootstrap Z 1.03 +- 1.07, data GOF >= 0.92). Earlier:") % (cut, cut, cut)
assert s.count(old) == 1, "FIXED_WP line not found"
open(p, "w").write(s.replace(old, new)); print("  FIXED_WP updated")
PY
cd $PLOT && PYTHONPATH=$PLOT/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts $ALP_PY scripts/scan_score_R_significance.py --write-json > $L/dyv_wp3_scan.log 2>&1
grep -E '^ +[0-9]+ +0\.' $L/dyv_wp3_scan.log
python3 - $PLOT/output/MVAcut_points_run3.json $NEWCUT <<'PY' || die "JSON not updated"
import json, sys
r = {e["mA"]: e["MVAcut"] for e in json.load(open(sys.argv[1]))["results"]}
print("  JSON mA1-4:", {m: r[m] for m in (1, 2, 3, 4)})
sys.exit(0 if abs(r[3] - float(sys.argv[2])) < 1e-6 and r[2] == 0.975 else 1)
PY

echo "=============== [2b] efficiency JSON (signal_eff_sumw.py in a clean env: under cmsenv cling clashes) ==============="
T=$(date +%s)
( cd $PLOT && env -i HOME=$HOME PATH=/usr/bin:/bin PYTHONPATH=$PLOT/lib $ALP_PY scripts/signal_eff_sumw.py > $L/dyv_wp3_signal_eff_sumw.log 2>&1 ) || die "signal_eff_sumw failed (wp3_signal_eff_sumw.log)"
$ALP_PY - $T $PLOT/output/sigEfficiencyVmA_{ele,muon}_byYear_5years_quadratic_interp_ma_points.json <<'PY' || die "efficiency JSON stale"
import json, os, sys
T = float(sys.argv[1]); ok = True
for p in sys.argv[2:]:
    d = json.load(open(p)); v = d.get("values", {})
    nz = sum(1 for y in v for m in range(1, 31) if float(v[y].get(str(m), 0)) > 0)
    good = os.path.getmtime(p) >= T and d.get("meta", {}).get("input_base", "").endswith("run3_bdt_scored_fsrfix") and nz == 150
    print("  %s: %s (%d/150 non-zero)" % (os.path.basename(p), "OK" if good else "STALE", nz)); ok &= good
sys.exit(0 if ok else 1)
PY

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"

echo "=============== [3] MVA cut mA3: data (full-read verified) + signal 5 eras ==============="
bash $D/regen_mvacut_data.sh 3 > $L/dyv_wp3_mva_data.log 2>&1 || die "data MVA cut failed (wp3_mva_data.log)"
tail -2 $L/dyv_wp3_mva_data.log
T=$(date +%s)
for y in $YEARS; do
  python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M3 --years $y > $L/dyv_wp3_mva_sig_$y.log 2>&1 &
done
wait
for y in $YEARS; do fresh $EOSMVA/sig/mA_M3/output_$y.root $T || die "signal MVA cut mA3 $y stale/missing"; done
python3 - $EOSMVA $YEARS <<'PY' || die "signal MVA-cut trees fail a full read"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B, years = sys.argv[1], sys.argv[2:]
for y in years:
    f = ROOT.TFile.Open("%s/sig/mA_M3/output_%s.root" % (B, y)); n = 0
    for k in f.GetListOfKeys():
        o = k.ReadObj()
        if o.InheritsFrom("TDirectory"):
            for kk in o.GetListOfKeys():
                t = kk.ReadObj()
                if t.InheritsFrom("TTree"):
                    for i in range(t.GetEntries()):
                        if t.GetEntry(i) <= 0: print("FAIL", y, t.GetName(), i); sys.exit(1)
                    n += t.GetEntries()
    print("  sig mA3 %-13s full read OK (%d entries over all trees)" % (y, n))
PY

echo "=============== [4] Tree2WS mA3 ==============="
sed -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(3)/' \
    -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(3)/' \
    $F/Trees2WS/run_tree2ws.sh > $S1/run_tree2ws_m3.sh
grep -q '^mAs_sig=(3)$' $S1/run_tree2ws_m3.sh && grep -q '^mAs_data=(3)$' $S1/run_tree2ws_m3.sh || die "tree2ws single-mass copy failed"
T=$(date +%s)
cd $F/Trees2WS && bash $S1/run_tree2ws_m3.sh > $L/dyv_wp3_tree2ws.log 2>&1
fresh $EOSMVA/data/mA_M3/ws/run3.root $T || die "data ws mA3 stale"
n=$(find $EOSMVA/sig/mA_M3/ws_Tree2WS -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 10 ] || die "only $n fresh signal ws for mA3"
echo "  data ws + $n signal ws fresh"

echo "=============== [5] envelope change mA24/25 + background fits mA3, mA24, mA25 (local) ==============="
python3 - $SRCF <<'PYS' || die "fTest source edit failed"
import sys
p = sys.argv[1]; s = open(p).read()
a = "case 25: return 6;"
assert s.count(a) == 1, "bestBernOrder(25) line not found"
s = s.replace(a, "case 25: return 5; /* 2026-10-01: Bern6 bias +0.212 after the DY-veto retrain -> 5 (AN 7.4) */")
b = '\tfunctionClasses.push_back("Laurent");\n'
assert s.count(b) == 1, "Laurent push_back not found"
s = s.replace(b, '\t// [PZ 2026-10-01] Lau1 (lowest order) above the 0.2 bias threshold after the 2024 DY-veto retrain\n'
                 '\t// at mA24 (-0.232) and mA25 (-0.215): removed there (AN Sec 7.4; cf. mA8 Pow1 2026-09-25).\n'
                 '\tif (mass_ALP != 24 && mass_ALP != 25) functionClasses.push_back("Laurent");\n')
open(p, "w").write(s); print("  fTest: bestBernOrder(25)=5, Laurent dropped at mA24/25")
PYS
( cd $BK && make -f makefile -B > $L/wp3env_build.log 2>&1 ) || die "background build failed (wp3env_build.log)"
[ "$(stat -c %Y $BK/bin/fTest_ALP_turnOn)" -ge "$(stat -c %Y $SRCF)" ] || die "fTest binary older than its source"
T=$(date +%s)
for m in 3 24 25; do bash $D/fit_bkg_single/fit_bkg_m$m.sh > $L/wp3env_bkg_m$m.log 2>&1 & done
wait
for m in 3 24 25; do
  MP=$FR/$m/CMS-HGG_mva_13p6TeV_multipdf.root
  fresh $MP $T && [ $(stat -c %s $MP) -ge 1024 ] || die "mA$m multipdf missing/stale/truncated"
  echo "  mA$m envelope:"; sed 's/^/    /' $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt
done
grep -qi "Laurent\|Lau" $FR/24/HZAmassInde_fTest/EnvelopeResults.txt && die "mA24 envelope still has Laurent"
grep -qi "Laurent\|Lau" $FR/25/HZAmassInde_fTest/EnvelopeResults.txt && die "mA25 envelope still has Laurent"
grep -q "Bern6\|Bernstein6\|order 6" $FR/25/HZAmassInde_fTest/EnvelopeResults.txt && die "mA25 envelope still has Bern6"

echo "=============== [6] signal model mA3 (fTest -> syst -> fit -> plotter -> effSigma) ==============="
T=$(date +%s)
for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
  sed -e 's/^mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )$/mAs=( 3 )/' -e 's/^mAs_lowMA=( .* )$/mAs_lowMA=( )/' \
      $F/shellScripts/sig_sys/$s.sh > $S1/$s.sh
  grep -q '^mAs=( 3 )$' $S1/$s.sh || die "single-mass copy of $s failed"
  echo "  [$(date '+%H:%M:%S')] $s"
  bash $S1/$s.sh > $L/dyv_wp3_sig_$s.log 2>&1
done
bash $F/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh > $L/dyv_wp3_sig_5.log 2>&1
for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/3_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || die "signal fit mA3 $y $ch stale"
done; done
echo "  signal fits mA3: 10/10 fresh"

echo "=============== [7] all datacards ==============="
T=$(date +%s)
cd $F/Datacard
sh 1_runjob_gen_datacard_makeYields.sh  > $L/dyv_wp3_dc_1.log 2>&1
sh 2_runjob_gen_datacard_makeDatacard.sh > $L/dyv_wp3_dc_2.log 2>&1
sh 3_rysn_datacard.sh                    > $L/dyv_wp3_dc_3.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt $T || bad=1; done
[ $bad -eq 0 ] || die "datacards stale"
for m in $(seq 1 30); do
  n=$(grep -c rateParam $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt)
  case " 1 2 3 4 5 6 7 8 9 10 15 20 25 30 " in *" $m "*) want=0;; *) want=10;; esac
  [ "$n" -eq "$want" ] || die "mA$m has $n rateParam lines (want $want)"
done
echo "  datacards 30/30 fresh, rateParams OK"

echo "=============== [8] t2w + blind limits (all 30), compare ==============="
T=$(date +%s)
cd $C && sh 2_text2ws.sh > $L/dyv_wp3_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/dyv_wp3_limits.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || bad=1; done
[ $bad -eq 0 ] || die "limits stale"
python3 - $C/output_combine_results $V/output_combine_results_$TAG > $L/dyv_wp3_changed.txt <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for m in range(1, 31):
    n, o = med(sys.argv[1], m), med(sys.argv[2], m)
    if abs(n / o - 1) > 1e-4: print(m, "%.5f %.5f %.4f" % (o, n, n / o))
PY
echo "  masses whose median limit moved (mA old new ratio):"; sed 's/^/    /' $L/dyv_wp3_changed.txt
CHANGED=$(awk '{print $1}' $L/dyv_wp3_changed.txt | tr '\n' ' ')
for m in 3 24 25; do echo "$CHANGED" | grep -qw $m || die "mA$m limit did not change -- the change did not propagate"; done

echo "=============== [9] condor: bias + impact for changed masses: $CHANGED ==============="
for m in $CHANGED; do
  cp -p $FZ/${m}_Datacard_leptons.root $FZ/${m}_Datacard_leptons.root.$TAG 2>/dev/null
  cp -p $C/root_t2w/${m}_Datacard_leptons.root $FZ/
  [ -d $BN/bias_outputs_highstat/mA_$m ] && mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG
done
T0=$(date +%s)
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $CHANGED 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $CHANGED 2>&1)
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "bias/impact submission failed"
echo "SUBMITTED $(date '+%F %T') clusters $CL" | tee -a $HB

echo "=============== [10] wait, merge, summarize bias for $CHANGED ==============="
while true; do
  q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  while read -r id mem; do
    [ -z "${id:-}" ] && continue
    if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
    echo "$(date '+%F %T') held $id mem $mem -> $new" >> $HB; condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null
  done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
  echo "$(date '+%F %T') idle=$idle running=$run held=$held" >> $HB
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { echo "only held left" >> $HB; break; }
  sleep 1800
done
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $CHANGED > $L/wp3env_bias_merge.log 2>&1
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/wp3env_bias_collect.log 2>&1
python3 - $BN/bias_outputs_highstat $C/output_impacts $T0 $CHANGED <<'PYB' | tee -a $HB
import json, os, sys
B, I, T0, ms = sys.argv[1], sys.argv[2], float(sys.argv[3]), [int(x) for x in sys.argv[4:]]
bad = 0
for m in ms:
    p = "%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (B, m, m)
    ip = "%s/%d_impacts.json" % (I, m)
    if not os.path.exists(p): print("mA%d bias MISSING" % m); bad += 1; continue
    fr = json.load(open(p))["fit_results"]
    items = sorted(((k, v["mean"]) for k, v in fr.items()), key=lambda x: -abs(x[1]))
    r = None
    if os.path.exists(ip) and os.path.getmtime(ip) >= T0:
        r = [q for q in json.load(open(ip))["POIs"] if q["name"] == "r"][0]["fit"][1]
    over = [k for k, v in items if abs(v) > 0.2]
    bad += bool(over) + (r is None)
    print("mA%-2d bias %s | impact r=%s%s" % (m, " ".join("%s:%+.3f" % x for x in items),
          "%.3f" % r if r is not None else "MISSING", "  ABOVE 0.2: %s" % over if over else ""))
print("VERDICT: %s" % ("PASSED" if bad == 0 else "CHECK (%d issues)" % bad))
PYB
echo "$(date '+%F %H:%M') DONE dyveto mA3 0.988 + env mA24/25 (log dyveto_wp3_env2425.log, heartbeat dyveto_wp3_env2425_heartbeat.txt)" >> $L/PENDING_RESULTS.txt
