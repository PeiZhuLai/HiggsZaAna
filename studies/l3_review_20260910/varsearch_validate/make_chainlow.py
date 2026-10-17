#!/usr/bin/env python3
"""Build a low-mass-only full chain (m_a = 1, 2, 3) for one low-mass candidate model (user 2026-10-15: run the full chain
with low_swapPt_relE2 at the x0.5 cuts; the high-mass model stays production, so m_a >= 4 is untouched).
Derived from chain4/ by explicit, asserted text replacements (every replacement must match exactly the expected
number of times):
  <V>/chainlow_<tag>/build_code_copy.sh  code copy pointing at  scored = sandbox_varsearch_v1/run3_bdt_scored_<cand>,
                                         MVA-cut outputs = sandbox_varsearch_v1/chainlow_<tag>/root_MVAcut
                                         + mass lists of the data/signal/Tree2WS/datacard drivers reduced to 1 2 3
  <V>/chainlow_<tag>/s1/                 S1 scripts on the same paths
  <V>/chainlow_<tag>/chainlow.sh         chain4.sh for 1 2 3: no interpolation JSON stage (no interpolated point),
                                         closure on the candidate's step-3 scans, bias outputs on EOS, output_combine_results
                                         link for the bias script (the 10-15 lesson)
  <V>/chainlow_<tag>/plot_output/MVAcut_points_run3_chain4.json   working points
Usage: make_chainlow.py <tag> <candidate model name> <cut mA1> <cut mA2> <cut mA3>
"""
import json
import os
import re
import shutil
import sys

V = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate"
SBE = "/eos/home-p/pelai/HZa/sandbox_varsearch_v1"
tag, cand, cuts = sys.argv[1], sys.argv[2], [float(x) for x in sys.argv[3:6]]
CH = "%s/chainlow_%s" % (V, tag)
MV = "%s/chainlow_%s" % (SBE, tag)
SC = "%s/run3_bdt_scored_%s" % (SBE, cand)
assert not os.path.exists(CH), CH + " exists"
os.makedirs(CH + "/logs"); os.makedirs(CH + "/plot_output"); os.makedirs(CH + "/s1")


def rep(s, a, b, n=1):
    c = s.count(a)
    assert c == n, "expected %d x %r, found %d" % (n, a[:80], c)
    return s.replace(a, b)


# working points
json.dump({"results": [{"mA": m, "MVAcut": c, "note": "chainlow_%s: %s x0.5 (closure + data GOF pass)" % (tag, cand)}
                       for m, c in zip((1, 2, 3), cuts)]}, open(CH + "/plot_output/MVAcut_points_run3_chain4.json", "w"), indent=1)

# code copy builder
b = open(V + "/chain4/build_code_copy.sh").read()
b = rep(b, "CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chain4", "CH=" + CH)
b = rep(b, "SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1", "SBE=" + MV + "   # MVA-cut outputs of this chain")
b = rep(b, '-e "s#/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix#$SBE/run3_bdt_scored#g"',
        '-e "s#/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix#%s#g"' % SC)
b = rep(b, "set -euo pipefail\n", "set -euo pipefail\nC4=%s/chain4/flashggFinalFit\nV4=$C4/shellScripts/bias/Condor\n" % V)
b = rep(b, "du -sh $S", r'''# low-mass only: mass lists of the drivers reduced to 1 2 3 (the sig_sys scripts keep theirs: the chain runs them per mass)
sed -i 's/^mAs = \[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30\]/mAs = [1,2,3]/' $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
grep -q '^mAs = \[1,2,3\]$' $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
python3 - $S/shellScripts/mva/run_apply_bdt_sig_6jobs.sh <<'PYS'
import re, sys
p = sys.argv[1]; s = open(p).read()
t, n = re.subn(r'samples=\(\n[^)]*\)', 'samples=(\n  "mA_M1" "mA_M2" "mA_M3"\n)', s)
assert n == 1; open(p, "w").write(t)
PYS
sed -i -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(1 2 3)/' -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(1 2 3)/' $S/Trees2WS/run_tree2ws.sh
grep -q '^mAs_sig=(1 2 3)$' $S/Trees2WS/run_tree2ws.sh && grep -q '^mAs_data=(1 2 3)$' $S/Trees2WS/run_tree2ws.sh
for f in 1_runjob_gen_datacard_makeYields.sh 2_runjob_gen_datacard_makeDatacard.sh; do
  sed -i 's/^mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )$/mAs=( 1 2 3 )/' $S/Datacard/$f
  grep -q '^mAs=( 1 2 3 )$' $S/Datacard/$f
done
# bias scripts: the chain4 copies (job work dir in the condor scratch, outputs to HS_OUT_BASE on EOS; 10-15 AFS lesson)
for f in runjob_bias_highstat.sh merge_bias_highstat.sh subjob_bias_highstat.sh; do
  cp $V4/$f $S/shellScripts/bias/Condor/$f
  sed -i "s#$C4#$S#g" $S/shellScripts/bias/Condor/$f
  grep -q 'HS_OUT_BASE' $S/shellScripts/bias/Condor/$f
done
du -sh $S''')
open(CH + "/build_code_copy.sh", "w").write(b)

# S1 scripts
for f in ("eval_alt_reweight.py", "make_s1_lnn_json.py"):
    shutil.copy2(V + "/chain4/s1/" + f, CH + "/s1/" + f)
s = open(CH + "/s1/eval_alt_reweight.py").read()
s = rep(s, '"%s/run3_bdt_scored"' % SBE, '"%s"' % SC)
s = rep(s, '"%s/root_MVAcut/sig"' % SBE, '"%s/root_MVAcut/sig"' % MV)
s = rep(s, V + "/chain4/plot_output/", CH + "/plot_output/")
open(CH + "/s1/eval_alt_reweight.py", "w").write(s)

# chain script
c = open(V + "/chain4/chain4.sh").read()
c = rep(c, "#!/usr/bin/env bash\n", "#!/usr/bin/env bash\n# chainlow_%s (make_chainlow.py, 2026-10-15): chain4.sh reduced to m_a = 1, 2, 3 for the low-mass candidate %s\n# (cuts %s); high-mass model = production, m_a >= 4 untouched. Original header:\n" % (tag, cand, cuts))
c = rep(c, "CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chain4", "CH=" + CH)
c = rep(c, "SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1\n", "SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1\nMV=%s\nSC=%s\n" % (MV, SC))
c = rep(c, 'MASS14="1 2 3 4 5 6 7 8 9 10 15 20 25 30"', 'MASS14="1 2 3"')
c = rep(c, "MASS30=$(seq 1 30 | tr '\\n' ' ')", 'MASS30="1 2 3"')
c = rep(c, 'INTERP="11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29"', 'INTERP=""')
c = rep(c, "$SBE/root_MVAcut", "$MV/root_MVAcut", c.count("$SBE/root_MVAcut"))
c = rep(c, "assert sorted(r) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30], sorted(r)", "assert sorted(r) == [1, 2, 3], sorted(r)")
c = rep(c, '''n=$(find $SBE/run3_bdt_scored -name '*.root' -newer $V/models/model_high.pkl | wc -l)
p=$(find $SBE/run3_bdt_scored -name '*.partial' | wc -l)
[ "$n" -eq 1208 ] && [ "$p" -eq 0 ] || die "scored copy incomplete: $n/1208 files newer than the models, $p partial"''',
        '''lst=$(xrdfs eosuser.cern.ch ls -R $(echo $SC | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#') 2>/dev/null)
n=$(echo "$lst" | grep -c '\\.root$'); p=$(echo "$lst" | grep -c '\\.partial$')
[ "$n" -eq 273 ] && [ "$p" -eq 0 ] || die "scored copy incomplete: $n/273 files (255 signal m_a 1-3 + 18 data/DY), $p partial"''')
c = rep(c, "eval_alt_reweight.py --out yields_chain4.csv", "eval_alt_reweight.py --masses 1,2,3 --out yields_chain4.csv")
c = rep(c, "if run 2; then", "if false; then   # no interpolated point among m_a 1-3: interpolation JSON not needed")
c = rep(c, "for m in range(1, 31):", "for m in (1, 2, 3):", c.count("for m in range(1, 31):"))
c = rep(c, "for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30):", "for m in (1, 2, 3):")
c = rep(c, "python3 - $SBE $WPJ <<'PY' | tee -a $L/chain4_heartbeat.txt || die \"data workspace counts differ from the closure-scan data\"",
        "python3 - $SBE $WPJ $MV <<'PY' | tee -a $L/chain4_heartbeat.txt || die \"data workspace counts differ from the closure-scan data\"")
c = rep(c, 'f = ROOT.TFile.Open("%s/root_MVAcut/data/mA_M%d/ws/run3.root" % (B, m))', 'f = ROOT.TFile.Open("%s/root_MVAcut/data/mA_M%d/ws/run3.root" % (sys.argv[3], m))')
c = rep(c, "closure_scan/new_drop_mA", "closure_scan/cand_%s_drop_mA" % cand, c.count("closure_scan/new_drop_mA"))
c = rep(c, "export -f bkg_one; export F FR SBE", "export -f bkg_one; export F FR SBE MV")
c = rep(c, 'hb "[9] closure with the sandbox xsec card"', 'ln -sfn output_combine_results_xsec $S/Combine/output_combine_results   # the bias script reads the limits there\nhb "[9] closure with the sandbox xsec card"')
c = rep(c, "ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec NCHUNKS=10 bash", "HS_OUT_BASE=$MV/bias_outputs_highstat ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec NCHUNKS=10 bash")
c = rep(c, "xargs -P6 -L1 bash -c 'run_interp $0 $1' < $L/closure_interp_cuts.txt", "xargs -r -P6 -L1 bash -c 'run_interp $0 $1' < $L/closure_interp_cuts.txt")
c = rep(c, '[ $(grep -c "closure Z=" $CH/closure_final_summary.txt) -eq 30 ]', '[ $(grep -c "closure Z=" $CH/closure_final_summary.txt) -eq 3 ]')
assert "$SBE/run3_bdt_scored" not in c
open(CH + "/chainlow.sh", "w").write(c)
os.chmod(CH + "/chainlow.sh", 0o755); os.chmod(CH + "/build_code_copy.sh", 0o755)
print("wrote", CH)
