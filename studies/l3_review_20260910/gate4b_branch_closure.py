#!/usr/bin/env python
"""
GATE 4b -- branch-level closure for the hadd stage.

GATE 4 checks that entries are conserved. That catches a whole tree being dropped, but
not a BRANCH being dropped or silently retyped. Both happen: hadd fixes the schema from
the FIRST input file, and this production has a genuine type clash --
    Warning in <TTree::CopyEntries>: the export leaf and the import leaf (year.year)
    do not have the same data type (Long64_t vs Double_t)
because `year` is a string column in the parquet, which converts to NaN/double for
2022-2023 and to int64 for 2024. Nothing downstream reads `year`, but the same
mechanism would hide the loss of a branch that does matter.

So: every hadd output must carry exactly the branches of its inputs, no fewer.

Usage: python gate4b_branch_closure.py
Output: logs_fsrfix/gate4b_report.txt
"""
import os, sys, time
import uproot

BASE   = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
REPORT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gate4b_report.txt"
YEARS  = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = ["M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M15","M20","M25","M30"]


def branches(path):
    if not os.path.exists(path):
        return None
    try:
        with uproot.open(path) as f:
            return set(f["inclusive"].keys())
    except Exception:
        return None


def main():
    cases = [("DYGto2LG/run3", "%s/DYGto2LG/run3.root" % BASE,
              ["%s/DYGto2LG/%s.root" % (BASE, y) for y in YEARS]),
             ("DYJetsToLL/run3", "%s/DYJetsToLL/run3.root" % BASE,
              ["%s/DYJetsToLL/%s.root" % (BASE, y) for y in YEARS]),
             ("All_Bkg/run3", "%s/All_Bkg/run3.root" % BASE,
              ["%s/DYJetsToLL/run3.root" % BASE, "%s/DYGto2LG/run3.root" % BASE]),
             ("Data/run3", "%s/Data/run3.root" % BASE,
              ["%s/Data/%s.root" % (BASE, y) for y in YEARS])]
    for m in MASSES:
        cases.append(("mA_%s/run3" % m, "%s/mA_%s/run3.root" % (BASE, m),
                      ["%s/mA_%s/%s.root" % (BASE, m, y) for y in YEARS]))
    cases.append(("All_Sig/run3", "%s/All_Sig/run3.root" % BASE,
                  ["%s/mA_%s/run3.root" % (BASE, m) for m in MASSES]))

    L = ["GATE 4b -- hadd branch closure (output branches must cover every input's)",
         "generated %s" % time.strftime('%F %T'), ""]
    bad = 0
    for label, outf, inputs in cases:
        ob = branches(outf)
        if ob is None:
            L.append("%-22s MISSING_OUTPUT" % label); bad += 1; continue
        union, per_input_missing = set(), {}
        for p in inputs:
            ib = branches(p)
            if ib is None:
                L.append("%-22s MISSING_INPUT %s" % (label, os.path.basename(p))); bad += 1; continue
            union |= ib
            miss = ib - ob
            if miss:
                per_input_missing[os.path.basename(p)] = sorted(miss)
        lost = union - ob
        if lost or per_input_missing:
            L.append("%-22s LOST_BRANCHES n=%d: %s" % (label, len(lost), sorted(lost)[:8]))
            for k, v in per_input_missing.items():
                L.append("      from %-18s %s" % (k, v[:8]))
            bad += 1
        else:
            L.append("%-22s OK   branches=%d" % (label, len(ob)))
    L += ["", "VERDICT: %s" % ("FAILED" if bad else "PASSED")]
    txt = "\n".join(L)
    open(REPORT, "w").write(txt + "\n")
    print(txt)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
