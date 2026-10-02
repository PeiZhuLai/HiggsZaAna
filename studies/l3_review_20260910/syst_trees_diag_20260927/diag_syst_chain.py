#!/usr/bin/env python3
"""READ-ONLY diagnosis: where do signal systematic-variation events get lost?

For each (mA, era) cell and each syst, report per stage:
  P  parquet: chunk rows (find, excl. merged), merged rows, mtime race, sumw, key overlap w/ nominal,
              row-order alignment (same key at same row index), same-10-bucket fraction for common keys
  S  scored root (run3_bdt_scored_fsrfix): inclusive/test entries, test keys common w/ nominal test
  M  MVAcut output (root_MVAcut/sig): entries & sumw per lepton (what Trees2WS reads)
  F  "fixed" emulation: syst inclusive restricted to events whose (run,lumi,event) is in the NOMINAL test
     split (= a key-stable split), then nominal-pass join + syst MVA cut, as apply_bdt_sig does.
Usage: python3 diag_syst_chain.py mA_M5:2024 mA_M20:2022preEE ...
"""
import os, sys, glob, json
import numpy as np
import pyarrow.parquet as pq
import uproot

PARQ = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
SCORED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
MVAC = "/eos/home-p/pelai/HZa/root_MVAcut/sig"
CUTS = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
SYSTS = ["FNUF", "Material", "Electron_scale", "Electron_smear", "Muon_scale", "Muon_smear",
         "Photon_scale", "Photon_smear"]
ID = ["run", "luminosityBlock", "event"]

cuts = {int(r["mA"]): float(r["MVAcut"]) for r in json.load(open(CUTS))["results"]}


def suffix(corr):
    p = corr.split("_")
    return "_" + p[0] + "".join(x.capitalize() for x in p[1:-1]) + p[-1].capitalize() + "01sigma"


def keys_of(df_or_arrs):
    return list(zip(*[np.asarray(df_or_arrs[c]).astype(np.int64) for c in ID]))


def parquet_stage(sample_dir, corr, nom_keys):
    chunks = []
    for root, _, files in os.walk(sample_dir):
        for fn in files:
            if fn.endswith("_%s.parquet" % corr) and fn.startswith("output_job_"):
                chunks.append(os.path.join(root, fn))
    crow = sum(pq.ParquetFile(c).metadata.num_rows for c in chunks)
    cmax = max(os.path.getmtime(c) for c in chunks) if chunks else 0
    mf = os.path.join(sample_dir, "merged_%s.parquet" % corr)
    t = pq.read_table(mf, columns=ID + ["weight_central"]).to_pandas()
    race = os.path.getmtime(mf) < cmax
    k = keys_of(t)
    out = dict(nchunk=len(chunks), chunk_rows=crow, merged_rows=len(t), race=race,
               sumw=float(t.weight_central.sum()), uniq=len(set(k)))
    if nom_keys is not None:
        nk = nom_keys
        n = min(len(nk), len(k))
        out["common"] = len(set(k) & set(nk))
        out["aligned_rows"] = int(sum(1 for i in range(n) if k[i] == nk[i]))
        # same bucket (idx%10) fraction for common keys
        pos_n = {kk: i for i, kk in enumerate(nk)}
        same = tot = test_same = test_tot = 0
        for i, kk in enumerate(k):
            j = pos_n.get(kk)
            if j is None:
                continue
            tot += 1
            same += (i % 10) == (j % 10)
            if j % 10 >= 7:
                test_tot += 1
                test_same += (i % 10) >= 7
        out["same_bucket_frac"] = same / tot if tot else float("nan")
        out["nomtest_in_systtest_frac"] = test_same / test_tot if test_tot else float("nan")
    return out, k


def main(cells):
    for cell in cells:
        ma, era = cell.split(":")
        m = int(ma.split("_M")[1])
        cut = cuts[m]
        mva = "MVA_Score_%s" % ma
        sdir = os.path.join(PARQ, "%s_%s" % (ma, era))
        print("=" * 100)
        print("CELL %s %s  MVAcut=%.3f" % (ma, era, cut))
        # nominal parquet
        pn, nkeys = parquet_stage(sdir, "nominal", None)
        # nominal scored
        fs = uproot.open(os.path.join(SCORED, ma, "%s.root" % era))
        ntest = fs["test"].arrays(ID + [mva, "weight", "n_electrons", "n_muons"], library="pd")
        ninc_n = fs["inclusive"].num_entries
        nom_test_keys = set(keys_of(ntest))
        passn = ntest[ntest[mva] > cut]
        pass_keys = set(keys_of(passn))
        mv = uproot.open(os.path.join(MVAC, ma, "output_%s.root" % era))

        def mstage(suf):
            r = {}
            for lep in ("ele", "mu"):
                tn = "DiphotonTree/ggh_125_Za_%s_13p6TeV_cat0%s" % (lep, suf)
                if tn in mv:
                    w = mv[tn]["weight"].array(library="np")
                    r[lep] = (len(w), float(w.sum()))
                else:
                    r[lep] = (None, None)
            return r

        mn = mstage("")
        hdr = ("%-22s | %6s %8s %8s %5s %10s %8s %8s %6s %6s | %7s %6s %7s | %-24s %-24s | %-24s %-24s"
               % ("corr", "nchnk", "chunkR", "mergR", "race", "sumw", "common", "alignR", "sameB", "tInT",
                  "incl", "test", "tCommN", "MVAc ele (n,sumw)", "MVAc mu (n,sumw)", "FIX ele", "FIX mu"))
        print(hdr)
        print("%-22s | %6d %8d %8d %5s %10.4f %8s %8s %6s %6s | %7d %6d %7s | %-24s %-24s |"
              % ("nominal", pn["nchunk"], pn["chunk_rows"], pn["merged_rows"], pn["race"], pn["sumw"], "-", "-",
                 "-", "-", ninc_n, len(ntest), "-", "(%d, %.4f)" % mn["ele"], "(%d, %.4f)" % mn["mu"]))
        for syst in SYSTS:
            for ud in ("up", "down"):
                corr = "%s_%s" % (syst, ud)
                p, _ = parquet_stage(sdir, corr, nkeys)
                f2 = uproot.open(os.path.join(SCORED, "%s_%s" % (ma, corr), "%s.root" % era))
                stest = f2["test"].arrays(ID + [mva, "weight", "n_electrons", "n_muons"], library="pd")
                sinc = f2["inclusive"].arrays(ID + [mva, "weight", "n_electrons", "n_muons"], library="pd")
                tcomm = len(set(keys_of(stest)) & nom_test_keys)
                ms = mstage(suffix(corr))
                # FIX emulation: key-stable split (nominal test keys), nominal pass join, syst MVA cut
                kin = keys_of(sinc)
                sel = np.array([(kk in nom_test_keys) and (kk in pass_keys) for kk in kin], dtype=bool)
                fx = sinc[sel & (sinc[mva].to_numpy() > cut)]
                fe = fx[fx.n_electrons == 2]
                fm = fx[fx.n_muons == 2]
                print("%-22s | %6d %8d %8d %5s %10.4f %8d %8d %6.3f %6.3f | %7d %6d %7d | %-24s %-24s | %-24s %-24s"
                      % (corr, p["nchunk"], p["chunk_rows"], p["merged_rows"], p["race"], p["sumw"], p["common"],
                         p["aligned_rows"], p["same_bucket_frac"], p["nomtest_in_systtest_frac"],
                         len(sinc), len(stest), tcomm,
                         "(%s, %.4f)" % (ms["ele"][0], ms["ele"][1] or 0), "(%s, %.4f)" % (ms["mu"][0], ms["mu"][1] or 0),
                         "(%d, %.4f)" % (len(fe), fe.weight.sum()), "(%d, %.4f)" % (len(fm), fm.weight.sum())))
        sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
