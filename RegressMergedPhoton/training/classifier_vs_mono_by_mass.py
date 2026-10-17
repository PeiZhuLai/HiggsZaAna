"""Diphoton-vs-single-photon separation of a classifier, per generated mass.

The classifier's own validation set has no m_gen, so the train/val FILE split of
the pack it was trained on is reconstructed here with the same code and seed
(build_dataset.py defaults), and checked against the row counts that build
printed. Only clusters from validation files are scored, so nothing here was
seen in training. Within those files ALL clusters are used, not just the
energy-balanced subset.

Single photons are reweighted to each mass point's cluster-energy spectrum,
because the score can depend on energy and GJets photons are softer than the
signal.

    python classifier_vs_mono_by_mass.py --model <classifier_best.pt>
"""
import argparse
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_dataset as bd  # noqa: E402
from models import MLPhotonClassifier  # noqa: E402

EXPECTED_VAL = {bd.LABEL_MONO: 20629, bd.LABEL_DI: 18787, bd.LABEL_HAD: 25697}


def val_files(base, seed=1234, per_class=200_000, n_ebins=20, moe_max=0.075, val_frac=0.1):
    rng = np.random.default_rng(seed)
    files, fidx, ridx, lab, ene, moe = bd.scan(base)
    strata = {i: os.path.basename(os.path.dirname(f)) for i, f in enumerate(files)}
    bad = (lab == bd.LABEL_DI) & ~(np.isfinite(moe) & (moe > 0) & (moe < moe_max))
    lab = lab.copy(); lab[bad] = -1
    sel = bd.balance_in_energy(lab, ene, per_class, n_ebins, rng)
    present = np.unique(fidx[sel])
    vf = []
    for tag in sorted({strata[int(f)] for f in present}):
        tf = np.array([f for f in present if strata[int(f)] == tag])
        perm = rng.permutation(len(tf))
        vf.extend(tf[perm[:max(int(round(len(tf) * val_frac)), 1)]].tolist())
    in_val = np.isin(fidx[sel], vf)
    got = {L: int(np.sum(lab[sel][in_val] == L)) for L in EXPECTED_VAL}
    print(f"[split] reconstructed val rows {got}  expected {EXPECTED_VAL}")
    if got != EXPECTED_VAL:
        raise SystemExit("[FATAL] split reconstruction does not match the training build")
    return [files[i] for i in vf], moe_max


@torch.no_grad()
def di_score(model, img, bs=2048):
    out = [torch.exp(model(torch.from_numpy(img[i:i + bs])))[:, 1].numpy()
           for i in range(0, len(img), bs)]
    return np.concatenate(out).astype(np.float64)


def auc(s_sig, s_bkg, w_bkg):
    """P(score_sig > score_bkg) with background weights."""
    order = np.argsort(s_bkg)
    sb, wb = s_bkg[order], w_bkg[order]
    cw = np.concatenate([[0.0], np.cumsum(wb)])
    k = np.searchsorted(sb, s_sig, side="left")
    return float(np.mean(cw[k]) / cw[-1])


def eff_at_fake(s_sig, s_bkg, w_bkg, fake):
    order = np.argsort(-s_bkg)
    cum = np.cumsum(w_bkg[order]) / w_bkg.sum()
    thr = s_bkg[order][min(np.searchsorted(cum, fake), len(order) - 1)]
    return float(np.mean(s_sig > thr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--base", default=bd.DEFAULT_BASE)
    ap.add_argument("--max-had", type=int, default=60000)
    args = ap.parse_args()

    vfiles, moe_max = val_files(args.base)
    img, lab, mgen, ene = [], [], [], []
    for f in vfiles:
        with np.load(f) as z:
            L = z["label"].astype(np.int64)
            ok = (L == bd.LABEL_MONO) | (L == bd.LABEL_HAD) | (
                (L == bd.LABEL_DI) & (z["moe_true"] > 0) & (z["moe_true"] < moe_max))
            img.append(z["image"][ok].astype(np.float32)); lab.append(L[ok])
            mgen.append(z["m_gen"][ok]); ene.append(z["energy"][ok])
    img = np.concatenate(img).reshape(-1, 1, 30, 30); lab = np.concatenate(lab)
    mgen = np.concatenate(mgen); ene = np.concatenate(ene)
    had = np.flatnonzero(lab == bd.LABEL_HAD)
    if len(had) > args.max_had:
        drop = np.random.default_rng(0).choice(had, len(had) - args.max_had, replace=False)
        keep = np.ones(len(lab), bool); keep[drop] = False
        img, lab, mgen, ene = img[keep], lab[keep], mgen[keep], ene[keep]
    print(f"[data] {len(vfiles)} val files: mono {np.sum(lab==0)} di {np.sum(lab==1)} "
          f"had {np.sum(lab==2)}")

    model = MLPhotonClassifier()
    model.load_state_dict(torch.load(args.model, map_location="cpu")); model.eval()
    s = di_score(model, img)

    mono = lab == bd.LABEL_MONO; hadm = lab == bd.LABEL_HAD
    print(f"\n  {'m_gen':>5s} {'N_sig':>6s} {'E_med':>6s} | {'AUC_mono':>8s} {'AUC_monoEw':>10s} "
          f"{'AUC_had':>7s} | {'sig>0.5':>7s} {'mono>0.5':>8s} {'monoEw>0.5':>10s} {'had>0.5':>7s} | "
          f"{'eff@mono1%':>10s} {'eff@mono10%':>11s}")
    for mg in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        sig = (lab == bd.LABEL_DI) & np.isclose(mgen, mg)
        if sig.sum() < 50:
            print(f"  {mg:>5.1f} {int(sig.sum()):>6d}  (too few)"); continue
        edges = np.quantile(ene[sig], np.linspace(0, 1, 21)); edges[0], edges[-1] = 0, np.inf
        hs, _ = np.histogram(ene[sig], edges); hm, _ = np.histogram(ene[mono], edges)
        wbin = np.divide(hs / hs.sum(), hm / hm.sum(), out=np.zeros(len(hs)), where=hm > 0)
        w_ew = wbin[np.clip(np.digitize(ene[mono], edges) - 1, 0, 19)]
        w1 = np.ones(int(mono.sum()))
        print(f"  {mg:>5.1f} {int(sig.sum()):>6d} {np.median(ene[sig]):>6.1f} | "
              f"{auc(s[sig], s[mono], w1):>8.3f} {auc(s[sig], s[mono], w_ew):>10.3f} "
              f"{auc(s[sig], s[hadm], np.ones(int(hadm.sum()))):>7.3f} | "
              f"{np.mean(s[sig] > 0.5):>7.3f} {np.mean(s[mono] > 0.5):>8.3f} "
              f"{np.average(s[mono] > 0.5, weights=w_ew):>10.3f} {np.mean(s[hadm] > 0.5):>7.3f} | "
              f"{eff_at_fake(s[sig], s[mono], w_ew, 0.01):>10.3f} "
              f"{eff_at_fake(s[sig], s[mono], w_ew, 0.10):>11.3f}")


if __name__ == "__main__":
    main()
