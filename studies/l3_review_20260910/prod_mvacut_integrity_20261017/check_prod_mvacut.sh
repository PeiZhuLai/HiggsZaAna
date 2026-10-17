#!/usr/bin/env bash
# Read-only integrity check of the production MVA-cut outputs (root_MVAcut/{data,sig}), after a sandbox data file
# written by apply_bdt_data.py through the EOS FUSE mount came out with a damaged TTree (uproot read 15 events, ROOT 2,
# hadd printed the error and returned 0). Production writes the same way.
#   data (30 mA): independent count from the production scored data (MVA_Score_mA_M<m> >= cut of the nearest anchor in
#     Plot/output/MVAcut_points_run3.json, the apply_bdt_data.py rule) == uproot entries of output_<era>.root ==
#     ROOT entries (every entry read, through xrootd) ; run3.root == sum of eras (ROOT and uproot) ; ws data sumEntries
#     == run3 ; no ROOT error printed.
#   signal (14 anchors x 5 eras, every tree): ROOT entries (every entry read, xrootd) == uproot entries ; no ROOT error.
# Output: logs/uproot.json, logs/root_check.log, logs/summary.txt. Nothing is written outside this directory.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
PYU=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
env -i HOME=$HOME PATH=/usr/bin:/bin $PYU - > $D/logs/uproot.json 2> $D/logs/uproot.err <<'PY'
import json, uproot, numpy as np
SC = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/Data/"
MV = "/eos/home-p/pelai/HZa/root_MVAcut/"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
cuts = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
near = {m: cuts[min(cuts, key=lambda k: abs(k - m))] for m in range(1, 31)}
out = {"cuts": near, "indep": {}, "data_up": {}, "sig_up": {}}
for e in ERAS:
    t = uproot.open(SC + e + ".root")["inclusive"]
    a = t.arrays(["MVA_Score_mA_M%d" % m for m in range(1, 31)] + ["run", "luminosityBlock", "event"], library="np")
    ids = set(zip(a["run"], a["luminosityBlock"], a["event"]))
    out.setdefault("dup", {})[e] = len(a["event"]) - len(ids)
    for m in range(1, 31):
        out["indep"]["%d_%s" % (m, e)] = int((a["MVA_Score_mA_M%d" % m] >= near[m]).sum())
for m in range(1, 31):
    for f in ERAS + ["run3"]:
        p = MV + "data/mA_M%d/%s.root" % (m, f if f == "run3" else "output_" + f)
        out["data_up"]["%d_%s" % (m, f)] = uproot.open(p)["DiphotonTree/Data_13p6TeV"].num_entries
for m in sorted(cuts):
    for e in ERAS:
        fi = uproot.open(MV + "sig/mA_M%d/output_%s.root" % (m, e))
        out["sig_up"]["%d_%s" % (m, e)] = {k: v.num_entries for k, v in fi.items(cycle=False) if hasattr(v, "num_entries")}
print(json.dumps(out))
PY
[ -s $D/logs/uproot.json ] || { echo "uproot pass failed: $D/logs/uproot.err"; exit 1; }
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && scramv1 runtime -sh)"; set -u
python3 - $D/logs/uproot.json > $D/logs/root_check.log 2>&1 <<'PY'
import sys, json, ROOT
ROOT.gROOT.SetBatch(True)
U = json.load(open(sys.argv[1])); X = "root://eosuser.cern.ch//eos/user/p/pelai/HZa/root_MVAcut/"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
def nread(f, tn):
    t = f.Get(tn)
    if not t: return -1
    n = t.GetEntries()
    for i in range(n):
        if t.GetEntry(i) <= 0: return -2
    return n
bad = []
print("duplicate run:lumi:event in scored data:", U["dup"])
for m in range(1, 31):
    s_root = s_up = s_ind = 0
    for e in ERAS:
        k = "%d_%s" % (m, e)
        f = ROOT.TFile.Open(X + "data/mA_M%d/output_%s.root" % (m, e)); nr = nread(f, "DiphotonTree/Data_13p6TeV"); f.Close()
        nu, ni = U["data_up"][k], U["indep"][k]
        s_root += nr; s_up += nu; s_ind += ni
        if not (nr == nu == ni): bad.append("data mA%d %s: ROOT %d uproot %d independent %d" % (m, e, nr, nu, ni))
    f = ROOT.TFile.Open(X + "data/mA_M%d/run3.root" % m); nr3 = nread(f, "DiphotonTree/Data_13p6TeV"); f.Close()
    f = ROOT.TFile.Open(X + "data/mA_M%d/ws/run3.root" % m); w = f.Get("CMS_hza_workspace"); nws = w.data("Data_13p6TeV").sumEntries() if w else -1; f.Close()
    nu3 = U["data_up"]["%d_run3" % m]
    ok = nr3 == nu3 == s_ind and int(round(nws)) == s_ind
    print("data mA%-2d cut %.5f  independent %4d  eras ROOT %4d uproot %4d  run3 ROOT %4d uproot %4d  ws %6.1f  %s"
          % (m, U["cuts"][str(m)], s_ind, s_root, s_up, nr3, nu3, nws, "OK" if ok else "BAD"))
    if not ok: bad.append("data mA%d run3/ws" % m)
nt = 0
for k, trees in sorted(U["sig_up"].items(), key=lambda kv: (int(kv[0].split("_")[0]), kv[0])):
    m, e = k.split("_", 1)
    f = ROOT.TFile.Open(X + "sig/mA_M%s/output_%s.root" % (m, e))
    for tn, nu in trees.items():
        nt += 1; nr = nread(f, tn)
        if nr != nu: bad.append("sig mA%s %s %s: ROOT %d uproot %d" % (m, e, tn, nr, nu))
    f.Close()
print("signal trees checked: %d" % nt)
print("BAD:", len(bad)); [print("  " + b) for b in bad]
sys.exit(1 if bad else 0)
PY
rc=$?
nerr=$(grep -c "Error in <" $D/logs/root_check.log)
{ grep -E "^data mA|^signal trees|^BAD|^  |^duplicate" $D/logs/root_check.log; echo "ROOT error lines: $nerr"; echo "rc=$rc"; } > $D/logs/summary.txt
cat $D/logs/summary.txt
