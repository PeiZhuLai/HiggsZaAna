#!/usr/bin/env bash
# Read-only: every branch of every tree of the production signal MVA-cut files (root_MVAcut/sig, 14 anchors x 5 eras)
# read with uproot; a damaged basket raises on read (as the pandas index branch of 3 data files did). Also checks every
# branch reads the tree's number of entries. Output: logs/scan_sig_branches.log (one line per bad branch + summary).
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - > $D/logs/scan_sig_branches.log 2>&1 <<'PY'
import uproot
B = "/eos/home-p/pelai/HZa/root_MVAcut/sig"; nb = nt = nf = 0; bad = []
for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30):
    for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        p = "%s/mA_M%d/output_%s.root" % (B, m, e); nf += 1
        f = uproot.open(p)
        for tn, t in f.items(cycle=False):
            if not hasattr(t, "num_entries"): continue
            nt += 1
            for b in t.keys():
                nb += 1
                try:
                    a = t[b].array(library="np")
                    if len(a) != t.num_entries: bad.append("%s %s %s: %d of %d entries" % (p, tn, b, len(a), t.num_entries))
                except Exception as x:
                    bad.append("%s %s %s: %s" % (p, tn, b, str(x).splitlines()[0][:120]))
        print("done mA%d %s" % (m, e), flush=True)
print("files %d trees %d branches %d bad %d" % (nf, nt, nb, len(bad)))
for b in bad: print("BAD", b)
PY
tail -1 $D/logs/scan_sig_branches.log
