#!/usr/bin/env python
"""
GATE 4 -- hadd closure for the FSR-fix ROOT files.

hadd takes its branch list from the FIRST input file and silently DROPS any tree whose
branches do not match, emitting only a Warning and exiting 0. The file stays valid, so
nothing downstream complains. Here the run3.root hadds put 2022preEE (NanoAODv12) first
and 2024 (NanoAODv15) last, which is exactly the shape that hides the loss.

The only check that catches it is conservation of entries:
    hadd output inclusive entries == sum of inputs' inclusive entries

Usage:  python gate4_hadd_closure.py
Output: logs_fsrfix/gate4_report.txt, non-zero exit if any hadd lost events.
"""
import os, sys, time
from concurrent.futures import ThreadPoolExecutor
import uproot

BASE   = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
REPORT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gate4_report.txt"
YEARS  = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = ["M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M15","M20","M25","M30"]


def entries(path):
    if not os.path.exists(path):
        return None
    try:
        with uproot.open(path) as f:
            return f["inclusive"].num_entries
    except Exception:
        return None


def checks():
    """(label, output_file, [input_files])"""
    out = []
    out.append(("DYGto2LG/run3", "%s/DYGto2LG/run3.root" % BASE,
                ["%s/DYGto2LG/%s.root" % (BASE, y) for y in YEARS]))
    out.append(("DYJetsToLL/2024", "%s/DYJetsToLL/2024.root" % BASE,
                ["%s/DYJetsTo2%s/2024.root" % (BASE, c) for c in ("E", "Mu", "Tau")]))
    out.append(("DYJetsToLL/run3", "%s/DYJetsToLL/run3.root" % BASE,
                ["%s/DYJetsToLL/%s.root" % (BASE, y) for y in YEARS]))
    out.append(("All_Bkg/run3", "%s/All_Bkg/run3.root" % BASE,
                ["%s/DYJetsToLL/run3.root" % BASE, "%s/DYGto2LG/run3.root" % BASE]))
    out.append(("Data/run3", "%s/Data/run3.root" % BASE,
                ["%s/Data/%s.root" % (BASE, y) for y in YEARS]))
    for m in MASSES:
        out.append(("mA_%s/run3" % m, "%s/mA_%s/run3.root" % (BASE, m),
                    ["%s/mA_%s/%s.root" % (BASE, m, y) for y in YEARS]))
    out.append(("All_Sig/run3", "%s/All_Sig/run3.root" % BASE,
                ["%s/mA_%s/run3.root" % (BASE, m) for m in MASSES]))
    return out


def run_one(c):
    label, outf, inputs = c
    o = entries(outf)
    ins = [entries(p) for p in inputs]
    if o is None:
        return (label, "MISSING_OUTPUT", None, None, os.path.basename(outf))
    missing = [os.path.basename(p) for p, v in zip(inputs, ins) if v is None]
    if missing:
        return (label, "MISSING_INPUT", o, None, ",".join(missing)[:60])
    tot = sum(ins)
    if o != tot:
        return (label, "LOST_EVENTS", o, tot, "diff %+d (%.3f%%)"
                % (o - tot, 100.0 * (o - tot) / tot if tot else 0.0))
    return (label, "OK", o, tot, "")


def main():
    cs = checks()
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = list(ex.map(run_one, cs))
    bad = [r for r in res if r[1] != "OK"]
    L = ["GATE 4 -- hadd closure (output entries == sum of input entries)",
         "generated %s" % time.strftime('%F %T'),
         "checks: %d   OK: %d   FAILED: %d" % (len(res), len(res) - len(bad), len(bad)), ""]
    L.append("%-22s %-14s %12s %12s  %s" % ("hadd", "status", "output", "sum(inputs)", "note"))
    for r in res:
        L.append("%-22s %-14s %12s %12s  %s"
                 % (r[0], r[1], r[2] if r[2] is not None else "-",
                    r[3] if r[3] is not None else "-", r[4]))
    L += ["", "VERDICT: %s" % ("FAILED" if bad else "PASSED")]
    txt = "\n".join(L)
    open(REPORT, "w").write(txt + "\n")
    print(txt)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
