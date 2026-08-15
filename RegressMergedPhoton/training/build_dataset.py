"""Pack the per-file npz into balanced, ready-to-train tensors.

Two datasets come out, because the two networks want different things:

  classifier : the three classes balanced IN ENERGY BINS, following
               AN2020-159 §5.5.1. Without this the classes differ by a factor
               ~2 in median energy (hadronic sits at 17 GeV, diphoton at 35)
               and the network can reach good accuracy by learning an energy
               cut instead of a shower shape.

  regressor  : diphoton clusters only, resampled so that log(m/E) is roughly
               FLAT. The natural sample is not flat -- it is a pile-up of 19
               discrete mass points -- whereas the paper trained on a gun that
               was flat in m/E by construction. Without the flattening the fit
               is dominated by whichever mass points happened to get more
               events, which is exactly how the current model ended up with a
               narrow usable domain.

Why pack at all: the full per-file set is ~15 GB of images, and training runs on
the IHEP GPU while the dumps live on CERN EOS. Packing first means only the
~2-3 GB that actually gets trained on has to be copied.

The scan pass reads only the small arrays (label/energy/moe_true), so choosing
what to keep costs almost nothing; images are read once, for the chosen rows.

Input : <npz_base>/<tag>/*.npz
Output: <out>/classifier_{train,val}.npz , <out>/regressor_{train,val}.npz

Usage (CERN, env `hzgdna`):
    /eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python build_dataset.py \\
        --out /eos/cms/store/group/phys_susy/pelai/HZa_merged/train_packed
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

DEFAULT_BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_npz"
LABEL_MONO, LABEL_DI, LABEL_HAD = 0, 1, 2
NAMES = {LABEL_MONO: "mono", LABEL_DI: "di", LABEL_HAD: "had"}


def scan(base):
    """Cheap metadata pass over every npz. Returns per-row arrays plus the file
    table needed to fetch images later."""
    files = sorted(glob.glob(os.path.join(base, "*", "*.npz")))
    if not files:
        raise SystemExit(f"[FATAL] no npz under {base}")
    fidx, ridx, lab, ene, moe = [], [], [], [], []
    for i, f in enumerate(files):
        with np.load(f) as z:
            n = len(z["label"])
            fidx.append(np.full(n, i, dtype=np.int32))
            ridx.append(np.arange(n, dtype=np.int32))
            lab.append(z["label"])
            ene.append(z["energy"])
            moe.append(z["moe_true"])
    return (files, np.concatenate(fidx), np.concatenate(ridx),
            np.concatenate(lab), np.concatenate(ene), np.concatenate(moe))


def balance_in_energy(lab, ene, per_class, n_bins, rng):
    """AN-style balancing: in each energy bin keep the same number of clusters
    from every class, capped by the rarest class in that bin."""
    edges = np.quantile(ene[np.isfinite(ene)], np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    keep = []
    per_bin_target = max(per_class // n_bins, 1)
    for b in range(n_bins):
        in_bin = (ene >= edges[b]) & (ene < edges[b + 1])
        idx_by_class = {L: np.flatnonzero(in_bin & (lab == L))
                        for L in (LABEL_MONO, LABEL_DI, LABEL_HAD)}
        n_take = min(per_bin_target, min(len(v) for v in idx_by_class.values()))
        if n_take == 0:
            continue
        for L, idx in idx_by_class.items():
            keep.append(rng.choice(idx, n_take, replace=False))
    return np.concatenate(keep) if keep else np.array([], dtype=np.int64)


def flatten_in_logmoe(idx, moe, n_target, n_bins, rng):
    """Resample diphoton rows so log(m/E) is ~flat, mimicking the flat-m/E
    particle gun the original model was trained on."""
    v = moe[idx]
    good = np.isfinite(v) & (v > 0)
    idx, v = idx[good], v[good]
    lv = np.log(v)
    edges = np.linspace(lv.min(), lv.max(), n_bins + 1)
    which = np.clip(np.digitize(lv, edges) - 1, 0, n_bins - 1)
    per_bin = max(n_target // n_bins, 1)
    keep = []
    for b in range(n_bins):
        pool = idx[which == b]
        if pool.size == 0:
            continue
        # sample with replacement only where the bin is genuinely thin, so the
        # sparse high-mass tail still gets represented
        keep.append(rng.choice(pool, min(per_bin, pool.size), replace=False)
                    if pool.size >= per_bin
                    else rng.choice(pool, per_bin, replace=True))
    return np.concatenate(keep)


def gather(files, fidx, ridx, sel, extra_keys=()):
    """Read images (and any extra arrays) for the selected global rows.

    Returns `sel` re-ordered by source file together with the arrays, so each
    npz is opened exactly once. The caller MUST use the returned ordering when
    it indexes the metadata -- reading in file order while labelling in the
    original order would pair every image with the wrong truth.
    """
    sel_sorted = sel[np.argsort(fidx[sel], kind="stable")]
    n = len(sel_sorted)
    img = np.empty((n, 30, 30), dtype=np.float32)
    extra = {k: np.empty(n, dtype=np.float32) for k in extra_keys}
    start = 0
    while start < n:
        f = fidx[sel_sorted[start]]
        stop = start
        while stop < n and fidx[sel_sorted[stop]] == f:
            stop += 1
        rows = ridx[sel_sorted[start:stop]]
        with np.load(files[f]) as z:
            img[start:stop] = z["image"][rows]
            for k in extra_keys:
                extra[k][start:stop] = z[k][rows]
        start = stop
    return img, extra, sel_sorted


def split_and_save(out_dir, name, arrays, group, val_frac, rng, strata=None):
    """Split by SOURCE FILE, STRATIFIED BY TAG.

    By file, because clusters from the same event live in the same npz and a
    row-wise split would put siblings on both sides and flatter the validation
    score.

    Stratified by tag, because each npz comes from exactly ONE dataset. An
    unstratified file split therefore samples whole mass points into val or
    train, not a cross-section of both. Measured on the first attempt: val ended
    up containing no M0p1-M0p6 at all, with a median m/E of 0.058 against the
    training set's 0.0077 -- so every validation number, and the early-stopping
    decision built on it, described a different population from the one being
    fitted. Taking val_frac of the files WITHIN each tag fixes it.
    """
    files_present = np.unique(group)
    if strata is None:
        strata = {int(f): "all" for f in files_present}

    val_files = []
    for tag in sorted({strata[int(f)] for f in files_present}):
        tag_files = np.array([f for f in files_present if strata[int(f)] == tag])
        if tag_files.size == 0:
            continue
        perm = rng.permutation(len(tag_files))
        n_val = max(int(round(len(tag_files) * val_frac)), 1)
        val_files.extend(tag_files[perm[:n_val]].tolist())
    is_val = np.isin(group, val_files)

    os.makedirs(out_dir, exist_ok=True)
    for tag, mask in (("train", ~is_val), ("val", is_val)):
        path = os.path.join(out_dir, f"{name}_{tag}.npz")
        np.savez_compressed(path, **{k: v[mask] for k, v in arrays.items()})
        print(f"  -> {path}  ({int(mask.sum())} rows from "
              f"{len(np.unique(group[mask]))} files, "
              f"{os.path.getsize(path)/1e6:.0f} MB)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-class", type=int, default=200_000,
                    help="classifier: target clusters per class (AN used 200k)")
    ap.add_argument("--n-regressor", type=int, default=500_000,
                    help="regressor: target diphoton clusters (AN used 500k)")
    ap.add_argument("--moe-max", type=float, default=0.075,
                    help="drop diphoton rows with m/E above this. The opening "
                         "angle is theta ~ 2*m/E, and DoPairings only merges "
                         "deposits within MATCH_DR = 0.15, so above m/E = 0.075 "
                         "the two photons are never in the same cluster: such a "
                         "row pairs a SINGLE-photon image with the two-photon "
                         "mass. Measured: the m/E>0.2 band has the photons 57 "
                         "crystals apart while the cluster spans 17 pixels. "
                         "Set to 0 to keep everything (not recommended).")
    ap.add_argument("--energy-bins", type=int, default=20)
    ap.add_argument("--moe-bins", type=int, default=40)
    ap.add_argument("--val-frac", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    # Line-buffer stdout: these runs are long and are watched through a log
    # file, where the default block buffering shows nothing for hours.
    sys.stdout.reconfigure(line_buffering=True)

    rng = np.random.default_rng(args.seed)

    print(f"[build] scanning {args.base} ...")
    files, fidx, ridx, lab, ene, moe = scan(args.base)
    # <base>/<tag>/<file>.npz -> tag, used to stratify the train/val split
    strata = {i: os.path.basename(os.path.dirname(f)) for i, f in enumerate(files)}
    print(f"[build] {len(files)} npz, {len(lab)} clusters "
          f"(mono {np.sum(lab==LABEL_MONO)}, di {np.sum(lab==LABEL_DI)}, "
          f"had {np.sum(lab==LABEL_HAD)})")

    # Restrict the diphoton class to the regime where the clustering can
    # actually produce a merged object. Applies to BOTH networks: an unmerged
    # pair is a single-photon image, so leaving it labelled "diphoton" teaches
    # the classifier that single photons are diphotons just as surely as it
    # gives the regressor an impossible target.
    if args.moe_max > 0:
        bad = (lab == LABEL_DI) & ~(np.isfinite(moe) & (moe > 0) & (moe < args.moe_max))
        lab = lab.copy()
        lab[bad] = -1
        print(f"[build] dropped {int(bad.sum())} diphoton rows with m/E >= "
              f"{args.moe_max} ({100*bad.sum()/max(np.sum(lab==LABEL_DI)+bad.sum(),1):.1f}% "
              f"of the diphoton class); {int(np.sum(lab==LABEL_DI))} remain")

    # ---- classifier ----
    print(f"\n[build] classifier: energy-balanced, target {args.per_class}/class")
    sel = balance_in_energy(lab, ene, args.per_class, args.energy_bins, rng)
    print(f"[build]   selected {len(sel)} rows "
          + ", ".join(f"{NAMES[L]} {int(np.sum(lab[sel]==L))}"
                      for L in (LABEL_MONO, LABEL_DI, LABEL_HAD)))
    e_sel = ene[sel]
    for L in (LABEL_MONO, LABEL_DI, LABEL_HAD):
        m = lab[sel] == L
        print(f"[build]   {NAMES[L]:>4s} energy med {np.median(e_sel[m]):7.2f} GeV")
    img, _, sel = gather(files, fidx, ridx, sel)      # sel now in file order
    split_and_save(args.out, "classifier",
                   {"image": img, "label": lab[sel].astype(np.int64),
                    "energy": ene[sel]},
                   fidx[sel], args.val_frac, rng, strata=strata)
    del img

    # ---- regressor ----
    print(f"\n[build] regressor: diphoton, flat in log(m/E), target {args.n_regressor}")
    di = np.flatnonzero(lab == LABEL_DI)
    sel_r = flatten_in_logmoe(di, moe, args.n_regressor, args.moe_bins, rng)
    print(f"[build]   selected {len(sel_r)} rows; m/E "
          f"[{moe[sel_r].min():.6f}, {moe[sel_r].max():.6f}], "
          f"median {np.median(moe[sel_r]):.5f}")
    img, extra, sel_r = gather(files, fidx, ridx, sel_r,
                               extra_keys=("eta", "energy", "m_gen"))
    split_and_save(args.out, "regressor",
                   {"image": img, "moe_true": moe[sel_r].astype(np.float32),
                    "eta": extra["eta"], "energy": extra["energy"],
                    "m_gen": extra["m_gen"]},
                   fidx[sel_r], args.val_frac, rng, strata=strata)

    print("\n[build] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
