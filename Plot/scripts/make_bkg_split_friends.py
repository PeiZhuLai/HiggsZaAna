#!/usr/bin/env python3
"""Friend trees that tag every scored background event with its BDT training split.

The BDT trainers use the 'train' tree of run3_bdt_inputs_fsrfix/All_Bkg/run3.root; 'validation'
and 'test' are never trained on. The split was made by row index inside each P2Root chunk, which
is lost once the chunks are merged, so every event of a scored background file is matched to its
split by (event, H_m, Z_m) rounded to 1 MeV; events whose masses moved by a re-production (~1e-4 GeV)
fall back to the same event number with |dH_m|, |dZ_m| < 0.01 GeV and exactly one candidate.
Output, one file per scored file and in the same entry order:
    <out>/<sample dir>/<era>.root : tree 'bkgsplit', branch bkg_split
        0 train, 1 validation, 2 test, -1 not matched, -2 key found in more than one split
1_prepare_dataVmc.py --bkgNonTrainOnly attaches these as friends (Plot_Helper._run3_build_chain).
"""
import argparse, os
import numpy as np
import uproot

SCORED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ALLBKG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix/All_Bkg/run3.root"
OUT = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_bkgsplit"
SOURCES = [("DYJetsToLL", e) for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix")] + \
          [(s, "2024") for s in ("DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau")] + \
          [("DYGto2LG_10to100", e) for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")]
SPLITS = {"train": 0, "validation": 1, "test": 2}


def keys(a):
    return zip(a["event"].astype(np.int64).tolist(), np.round(a["H_m"], 3).tolist(), np.round(a["Z_m"], 3).tolist())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scored", default=SCORED)
    ap.add_argument("--allbkg", default=ALLBKG)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    split_of = {}; by_event = {}
    f = uproot.open(args.allbkg)
    for name, code in SPLITS.items():
        b = f[name].arrays(["event", "H_m", "Z_m"], library="np")
        for k in keys(b):
            split_of[k] = code if split_of.get(k, code) == code else -2
        for ev, hm, zm in zip(b["event"].astype(np.int64).tolist(), b["H_m"].tolist(), b["Z_m"].tolist()):
            by_event.setdefault(ev, []).append((hm, zm, code))
    print("split map: %d keys, %d ambiguous" % (len(split_of), sum(v == -2 for v in split_of.values())))

    tot = {c: 0 for c in (-2, -1, 0, 1, 2)}
    for sample, era in SOURCES:
        src = os.path.join(args.scored, sample, era + ".root")
        a = uproot.open(src)["inclusive"].arrays(["event", "H_m", "Z_m"], library="np")
        code = np.array([split_of.get(k, -1) for k in keys(a)], dtype=np.int32)
        for i in np.where(code == -1)[0]:
            cand = [c for hm, zm, c in by_event.get(int(a["event"][i]), [])
                    if abs(hm - a["H_m"][i]) < 0.01 and abs(zm - a["Z_m"][i]) < 0.01]
            if len(cand) == 1:
                code[i] = cand[0]
            elif len(cand) > 1:
                code[i] = -2
        dst = os.path.join(args.out, sample, era + ".root")
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with uproot.recreate(dst) as g:
            g["bkgsplit"] = {"bkg_split": code}
        n = len(code)
        assert uproot.open(dst)["bkgsplit"].num_entries == n
        c = {k: int((code == k).sum()) for k in tot}
        for k in tot:
            tot[k] += c[k]
        print("%-18s %-13s n=%7d train %.3f  val+test %.3f  unmatched %.4f  ambiguous %d -> %s" % (
            sample, era, n, c[0] / n, (c[1] + c[2]) / n, c[-1] / n, c[-2], dst))
    n = sum(tot.values())
    print("total n=%d train %.4f val+test %.4f unmatched %.4f ambiguous %d" % (
        n, tot[0] / n, (tot[1] + tot[2]) / n, tot[-1] / n, tot[-2]))


if __name__ == "__main__":
    main()
