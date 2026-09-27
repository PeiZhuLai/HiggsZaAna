"""Join the merged-selection flags onto the resolved dataVmc ROOT trees.

The resolved dataVmc inputs (root_P2Root/run3_bdt_scored_fsrfix) carry none of
the MLPhoton information, so "dataVmc after the merged selection" needs the
merged flags attached first. They are matched on (run, luminosityBlock, event),
which both sides carry, and written to a NEW tree so the resolved inputs stay
untouched -- dataVmc then runs unchanged against the new directory.

Added branches:
    pass_allcuts_merged_ML   uint8   the merged selection flag (0 when the event
                                     is absent from the merged parquet)
    has_merged_info          uint8   1 if the event was found at all -- lets you
                                     tell "failed the merged selection" apart
                                     from "was never processed", which otherwise
                                     look identical and quietly bias any ratio
    MLPhoton_lead_mass       float   reco merged-photon mass (NaN if absent)
    MLPhoton_lead_diphotonScore float
    n_MLPhoton_diphoton      int32

Note on coverage: the merged friend parquet only has 2024 for the backgrounds
(2024 uses the flavour-split DYJetsTo2E/2Mu/2Tau instead of an inclusive
DYJetsToLL -- see Plot/lib/Analyzer_Configs.py), so only 2024 can be compared
data-vs-MC. Other eras are copied through with has_merged_info = 0 rather than
silently dropped.

Input : <resolved>/<sample>/<era>.root  +  <friend>/<tag>*/merged_nominal.parquet
Output: <out>/<sample>/<era>.root

Usage (CERN, env `hza_ana`):
    python add_merged_flag.py --samples Data,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau --eras 2024
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pyarrow.parquet as pq
import uproot

RESOLVED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
FRIEND = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend"
OUT = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_mergedflag"

MERGED_COLS = ["run", "luminosityBlock", "event", "pass_allcuts_merged_ML",
               "MLPhoton_lead_mass", "MLPhoton_lead_diphotonScore",
               "n_MLPhoton_diphoton"]

# resolved sample name -> the friend directory prefix
FRIEND_PREFIX = {
    "Data": "Data",
    "DYGto2LG_10to100": "Bkg_DYGto2LG_10to100",
    "DYGto2LG_10to50": "Bkg_DYGto2LG_10to50",
    "DYGto2LG_50to100": "Bkg_DYGto2LG_50to100",
    "DYJetsTo2E": "Bkg_DYJetsTo2E",
    "DYJetsTo2Mu": "Bkg_DYJetsTo2Mu",
    "DYJetsTo2Tau": "Bkg_DYJetsTo2Tau",
    "DYJetsToLL": "Bkg_DYJetsToLL",
}


def load_merged_index(prefix, era, bases=None):
    """(run, lumi, event) -> merged quantities, for every friend dir of this era.

    `bases` is searched in order and the results concatenated: data and the
    backgrounds live under different trees (the backgrounds had to be reproduced
    with the MLPhoton join actually switched on -- see
    HiggsDNA/scripts/run_merged_bkg_ml_friend.sh), and a run wants both at once.
    """
    bases = bases or [FRIEND]
    files = []
    for base in bases:
        files += sorted(glob.glob(os.path.join(base, f"{prefix}_{era}*", "*",
                                               "merged_nominal.parquet")))
        # signal/data friend dirs are sometimes one level shallower
        files += sorted(glob.glob(os.path.join(base, f"{prefix}_{era}*",
                                               "merged_nominal.parquet")))
    files = sorted(set(files))
    if not files:
        return None, 0

    keys, vals = [], []
    for f in files:
        names = pq.read_schema(f).names
        cols = [c for c in MERGED_COLS if c in names]
        if not all(c in cols for c in ("run", "luminosityBlock", "event")):
            continue
        if "pass_allcuts_merged_ML" not in cols:
            # SKIP, do not fill with zeros. The old parquet_friend/ tree and the
            # reproduced parquet_friend_ML/ one describe the SAME events, so when
            # both are searched an event resolves to whichever copy sorts first --
            # and the old copy has no MLPhoton columns, so its events come back as
            # "did not pass". That is how DYJetsTo2E reported pass_merged 0 on
            # 2026-08-30 while its friend genuinely passes 91.83%.
            print(f"  SKIP: no pass_allcuts_merged_ML in {f} "
                  f"({len(names)} cols) -- MLPhoton join was never run for it",
                  flush=True)
            continue
        t = pq.read_table(f, columns=cols)
        run = t["run"].to_numpy().astype(np.int64)
        lumi = t["luminosityBlock"].to_numpy().astype(np.int64)
        evt = t["event"].to_numpy().astype(np.int64)
        keys.append(np.stack([run, lumi, evt], axis=1))
        v = {}
        for c in ("pass_allcuts_merged_ML", "MLPhoton_lead_mass",
                  "MLPhoton_lead_diphotonScore", "n_MLPhoton_diphoton"):
            v[c] = (t[c].to_numpy().astype(np.float64) if c in cols
                    else np.full(len(run), np.nan))
        vals.append(v)

    if not keys:
        return None, 0
    K = np.concatenate(keys)
    V = {c: np.concatenate([v[c] for v in vals]) for c in vals[0]}
    return (K, V), len(files)


def key_hash(run, lumi, evt):
    """Single int64 key. run/lumi fit in 32/32 bits together; event is large, so
    combine with a mix rather than a naive shift (which would collide)."""
    a = (run.astype(np.uint64) * np.uint64(1000003)) ^ lumi.astype(np.uint64)
    return (a * np.uint64(1000033)) ^ evt.astype(np.uint64)


def process(sample, era, args):
    src = os.path.join(args.resolved, sample, f"{era}.root")
    if not os.path.exists(src):
        return f"{sample}/{era}: no resolved input"
    dst_dir = os.path.join(args.out, sample)
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, f"{era}.root")

    prefix = FRIEND_PREFIX.get(sample)
    bases = [b.strip() for b in args.friend.split(",") if b.strip()]
    idx, n_files = (load_merged_index(prefix, era, bases) if prefix else (None, 0))

    with uproot.open(src) as fin:
        tname = [k for k in fin.keys() if k.rstrip(";1") == "inclusive"]
        tname = tname[0] if tname else list(fin.keys())[0]
        tree = fin[tname]
        data = tree.arrays(library="np")

    n = len(data["event"])
    flag = np.zeros(n, dtype=np.uint8)
    found = np.zeros(n, dtype=np.uint8)
    mmass = np.full(n, np.nan, dtype=np.float32)
    mscore = np.full(n, np.nan, dtype=np.float32)
    ndip = np.zeros(n, dtype=np.int32)

    matched = 0
    if idx is not None:
        K, V = idx
        mk = key_hash(K[:, 0], K[:, 1], K[:, 2])
        order = np.argsort(mk, kind="stable")
        mk_sorted = mk[order]
        rk = key_hash(data["run"].astype(np.int64),
                      data["luminosityBlock"].astype(np.int64),
                      data["event"].astype(np.int64))
        pos = np.searchsorted(mk_sorted, rk)
        pos_c = np.clip(pos, 0, len(mk_sorted) - 1)
        hit = mk_sorted[pos_c] == rk
        src_idx = order[pos_c[hit]]
        dst_idx = np.flatnonzero(hit)
        matched = int(hit.sum())
        found[dst_idx] = 1
        p = V["pass_allcuts_merged_ML"][src_idx]
        flag[dst_idx] = np.nan_to_num(p, nan=0.0).astype(np.uint8)
        mmass[dst_idx] = V["MLPhoton_lead_mass"][src_idx]
        mscore[dst_idx] = V["MLPhoton_lead_diphotonScore"][src_idx]
        ndip[dst_idx] = np.nan_to_num(V["n_MLPhoton_diphoton"], nan=0.0)[src_idx]

    data["pass_allcuts_merged_ML"] = flag
    data["has_merged_info"] = found
    data["MLPhoton_lead_mass"] = mmass
    data["MLPhoton_lead_diphotonScore"] = mscore
    data["n_MLPhoton_diphoton"] = ndip

    with uproot.recreate(dst) as fout:
        fout["inclusive"] = data

    return (f"{sample}/{era}: {n} events, {n_files} friend files, "
            f"matched {matched} ({100.0*matched/max(n,1):.1f}%), "
            f"pass_merged {int(flag.sum())}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolved", default=RESOLVED)
    ap.add_argument("--friend", default=FRIEND,
                    help="comma-separated friend parquet base dirs, searched in order")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--samples", default="Data,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau")
    ap.add_argument("--eras", default="2024")
    args = ap.parse_args()

    for sample in args.samples.split(","):
        for era in args.eras.split(","):
            print(process(sample.strip(), era.strip(), args), flush=True)
    print(f"\noutput: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
