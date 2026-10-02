#!/usr/bin/env python3
"""Merge dataVmc partial ROOT files by summing histograms with the same path.

Replacement for `hadd` in Condor/3_merge_dataVmc_condor.sh (2026-09-29). hadd took ~52 min per
tag on these files (61761 TH1F in the merged file; hadd time grows super-linearly with the number
of keys), this takes ~2 min. Same result: every histogram is Clone() of its first occurrence plus
TH1::Add of the others, in input order -- identical to hadd's TH1::Merge for same-binned histograms
(checked bin by bin, contents and errors, against the 2026-09-26 hadd output).

Directories are reproduced; only the highest cycle of each key is read (as hadd does).
Non-histogram objects are copied from their first occurrence.

Usage: merge_dataVmc_hists.py <target.root> <input1.root> [input2.root ...]
"""
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)


def latest_keys(d):
    best = {}
    for k in d.GetListOfKeys():
        n = k.GetName()
        if n not in best or k.GetCycle() > best[n].GetCycle():
            best[n] = k
    return list(best.values())


def collect(d, path, acc, order):
    for k in latest_keys(d):
        name = k.GetName()
        p = f"{path}/{name}" if path else name
        cls = ROOT.TClass.GetClass(k.GetClassName())
        if cls and cls.InheritsFrom("TDirectory"):
            if p not in acc:
                acc[p] = None
                order.append(p)
            collect(k.ReadObj(), p, acc, order)
            continue
        obj = k.ReadObj()
        if p not in acc:
            if obj.InheritsFrom("TH1"):
                obj = obj.Clone(name)
                obj.SetDirectory(0)
            acc[p] = obj
            order.append(p)
        elif obj.InheritsFrom("TH1"):
            acc[p].Add(obj)


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    target, inputs = argv[0], argv[1:]
    acc, order = {}, []
    for fn in inputs:
        f = ROOT.TFile.Open(fn)
        if not f or f.IsZombie():
            sys.exit(f"cannot open {fn}")
        collect(f, "", acc, order)
        f.Close()
    out = ROOT.TFile.Open(target, "RECREATE")
    nh = 0
    for p in order:
        parent, _, name = p.rpartition("/")
        d = out.GetDirectory(parent) if parent else out
        if acc[p] is None:
            d.mkdir(name)
            continue
        d.cd()
        acc[p].Write(name)
        nh += 1
    out.Close()
    print(f"[merge] {target}: {len(inputs)} inputs, {nh} objects")


if __name__ == "__main__":
    main(sys.argv[1:])
