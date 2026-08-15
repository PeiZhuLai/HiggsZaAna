"""Summarise the MLPhoton training set: class balance, energy and m/E coverage.

The numbers that decide the training configuration:

  * class counts per tag and overall -- how much rebalancing is needed
  * the ENERGY distribution per class -- AN2020-159 §5.5.1 balances the three
    classes in energy bins before training, because otherwise the classifier
    just learns an energy cut
  * the m/E coverage of the diphoton (regression) sample, compared against the
    [0.005, 0.07] domain the current model was trained on -- the whole point of
    retraining is to cover the low end this analysis actually lives in

Reads only the small arrays (npz is lazily loaded, so the 30x30 images are
never touched).

Input : <npz_base>/<tag>/*.npz
Output: stdout tables (optional --csv for the per-tag table)

Usage (CERN, env `hzgdna`):
    /eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python summarize_training_set.py
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np

DEFAULT_BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_npz"
LABEL_MONO, LABEL_DI, LABEL_HAD = 0, 1, 2
NAMES = {LABEL_MONO: "mono", LABEL_DI: "di", LABEL_HAD: "had"}

TRAIN_MOE_LO, TRAIN_MOE_HI = 0.005, 0.07   # current model's domain, AN Table 7

KEYS = ["label", "energy", "moe_true", "m_gen", "e_frac", "ml_moe",
        "ml_diphotonScore", "n_crystals"]


def load_tag(tag_dir):
    out = {k: [] for k in KEYS}
    n_dropped = 0
    for f in sorted(glob.glob(os.path.join(tag_dir, "*.npz"))):
        try:
            with np.load(f) as z:
                for k in KEYS:
                    if k in z:
                        out[k].append(z[k])
                if "n_dropped" in z:
                    n_dropped += int(z["n_dropped"])
        except Exception as e:                                  # noqa: BLE001
            print(f"  [warn] unreadable {os.path.basename(f)}: {e}")
    if not out["label"]:
        return None, 0
    return {k: (np.concatenate(v) if v else np.array([])) for k, v in out.items()}, n_dropped


def q(a, p):
    return float(np.quantile(a, p)) if a.size else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    tags = sorted(d for d in os.listdir(args.base)
                  if os.path.isdir(os.path.join(args.base, d)))
    if not tags:
        print(f"[FATAL] no tag dirs under {args.base}")
        return 2

    print("=" * 92)
    print("PER-TAG CLASS COUNTS")
    print("=" * 92)
    hdr = (f"{'tag':<28s} {'total':>9s} {'mono':>9s} {'di':>9s} {'had':>9s} "
           f"{'dropped':>8s} {'drop%':>6s}")
    print(hdr)
    print("-" * len(hdr))

    all_data = {}
    rows = []
    grand = {LABEL_MONO: 0, LABEL_DI: 0, LABEL_HAD: 0}
    grand_drop = 0

    for tag in tags:
        d, n_drop = load_tag(os.path.join(args.base, tag))
        if d is None:
            print(f"{tag:<28s} EMPTY")
            continue
        all_data[tag] = d
        lab = d["label"]
        counts = {L: int(np.sum(lab == L)) for L in (LABEL_MONO, LABEL_DI, LABEL_HAD)}
        for L in grand:
            grand[L] += counts[L]
        grand_drop += n_drop
        tot = lab.size
        dpct = 100.0 * n_drop / (tot + n_drop) if (tot + n_drop) else 0.0
        print(f"{tag:<28s} {tot:>9d} {counts[LABEL_MONO]:>9d} {counts[LABEL_DI]:>9d} "
              f"{counts[LABEL_HAD]:>9d} {n_drop:>8d} {dpct:>5.1f}%")
        rows.append(dict(tag=tag, total=tot, mono=counts[LABEL_MONO],
                         di=counts[LABEL_DI], had=counts[LABEL_HAD], dropped=n_drop))

    gtot = sum(grand.values())
    print("-" * len(hdr))
    print(f"{'TOTAL':<28s} {gtot:>9d} {grand[LABEL_MONO]:>9d} {grand[LABEL_DI]:>9d} "
          f"{grand[LABEL_HAD]:>9d} {grand_drop:>8d}")
    print(f"{'':28s} {'':>9s} {100*grand[LABEL_MONO]/gtot:>8.1f}% "
          f"{100*grand[LABEL_DI]/gtot:>8.1f}% {100*grand[LABEL_HAD]/gtot:>8.1f}%")

    # ---- energy per class (drives the AN-style energy balancing) ----
    print()
    print("=" * 92)
    print("ENERGY PER CLASS (GeV) -- classes must be balanced in energy before training")
    print("=" * 92)
    print(f"{'class':>6s} {'N':>10s} {'q05':>8s} {'q16':>8s} {'med':>8s} {'q84':>8s} {'q95':>8s}")
    for L in (LABEL_MONO, LABEL_DI, LABEL_HAD):
        e = np.concatenate([d["energy"][d["label"] == L] for d in all_data.values()])
        print(f"{NAMES[L]:>6s} {e.size:>10d} {q(e,.05):>8.2f} {q(e,.16):>8.2f} "
              f"{q(e,.5):>8.2f} {q(e,.84):>8.2f} {q(e,.95):>8.2f}")

    # ---- regression target coverage ----
    print()
    print("=" * 92)
    print("REGRESSION TARGET m/E, DIPHOTON SAMPLE (per signal mass point)")
    print(f"current model was trained on m/E in [{TRAIN_MOE_LO}, {TRAIN_MOE_HI}] only")
    print("=" * 92)
    hdr2 = (f"{'tag':<12s} {'N_di':>8s} {'E_med':>7s} {'moe q05':>8s} {'moe med':>8s} "
            f"{'moe q95':>8s} {'<0.005':>7s} {'>0.07':>6s} {'dipho>0.9':>10s}")
    print(hdr2)
    print("-" * len(hdr2))

    di_all = []
    for tag in tags:
        if not tag.startswith("sig_") or tag not in all_data:
            continue
        d = all_data[tag]
        m = d["label"] == LABEL_DI
        if m.sum() == 0:
            continue
        moe = d["moe_true"][m]
        moe = moe[np.isfinite(moe)]
        e = d["energy"][m]
        dip = d["ml_diphotonScore"][m]
        di_all.append(moe)
        print(f"{tag:<12s} {int(m.sum()):>8d} {q(e,.5):>7.1f} {q(moe,.05):>8.5f} "
              f"{q(moe,.5):>8.5f} {q(moe,.95):>8.5f} "
              f"{100*np.mean(moe<TRAIN_MOE_LO):>6.1f}% {100*np.mean(moe>TRAIN_MOE_HI):>5.1f}% "
              f"{100*np.nanmean(dip>0.9):>9.1f}%")

    if di_all:
        moe = np.concatenate(di_all)
        print("-" * len(hdr2))
        print(f"{'ALL SIGNAL':<12s} {moe.size:>8d} {'':>7s} {q(moe,.05):>8.5f} "
              f"{q(moe,.5):>8.5f} {q(moe,.95):>8.5f} "
              f"{100*np.mean(moe<TRAIN_MOE_LO):>6.1f}% {100*np.mean(moe>TRAIN_MOE_HI):>5.1f}%")
        print()
        print(f"  target m/E full range: [{moe.min():.6f}, {moe.max():.6f}]")
        print(f"  fraction below the old training floor: "
              f"{100*np.mean(moe < TRAIN_MOE_LO):.1f}%  <-- what retraining must cover")

    # ---- current classifier efficiency on true diphotons ----
    print()
    print("=" * 92)
    print("CURRENT CLASSIFIER EFFICIENCY ON TRUE DIPHOTONS (the gap to be closed)")
    print("=" * 92)
    print(f"{'tag':<12s} {'N_di':>8s} {'eff(dipho>0.9)':>15s} {'eff(>0.5)':>10s}")
    for tag in tags:
        if not tag.startswith("sig_") or tag not in all_data:
            continue
        d = all_data[tag]
        m = d["label"] == LABEL_DI
        dip = d["ml_diphotonScore"][m]
        dip = dip[np.isfinite(dip)]
        if dip.size == 0:
            continue
        print(f"{tag:<12s} {dip.size:>8d} {100*np.mean(dip>0.9):>14.1f}% "
              f"{100*np.mean(dip>0.5):>9.1f}%")

    if args.csv and rows:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\n[summary] per-tag counts -> {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
