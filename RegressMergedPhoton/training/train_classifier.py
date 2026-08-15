"""Train the MLPhoton 3-class cluster classifier.

Classes (index order fixed by MLPhotonProducer): 0 = mono, 1 = di, 2 = hadronic.
Input is the RAW energy image -- the deployed classifier is fed unnormalised
crystal energies, so training must be too.

What success looks like: the deployed model recovers only 7-27% of true a->gg
clusters depending on mass (see ../README.md §5b), and that efficiency is what
limits the merged analysis. So the metric to watch is the DIPHOTON RECALL per
signal mass, not overall accuracy -- with three classes present in different
proportions, accuracy can look fine while the diphoton class is being thrown
away.

The packed training set is already energy-balanced (build_dataset.py), which
matters: unbalanced, the network can score well by learning "hadronic clusters
are the soft ones" instead of learning shower shapes.

Input : <packed>/classifier_{train,val}.npz
Output: <out>/classifier_best.pt , <out>/classifier_last.pt , <out>/history.csv

Usage:
    python train_classifier.py --packed <dir> --out <dir> --epochs 150
    # add --warm-start <classifier.onnx> to initialise from the current model
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import MLPhotonClassifier, load_onnx_weights  # noqa: E402

NAMES = ["mono", "di", "had"]


def load_split(packed, split):
    path = os.path.join(packed, f"classifier_{split}.npz")
    with np.load(path) as z:
        img = z["image"].astype(np.float32).reshape(-1, 1, 30, 30)
        lab = z["label"].astype(np.int64)
        energy = z["energy"].astype(np.float32)
    counts = {i: int(np.sum(lab == i)) for i in range(3)}
    print(f"[data] {split}: {len(lab)} rows  " +
          "  ".join(f"{NAMES[i]} {counts[i]}" for i in range(3)))
    return (torch.from_numpy(img), torch.from_numpy(lab), torch.from_numpy(energy))


@torch.no_grad()
def evaluate(model, loader, device):
    model.eval()
    logps, labs = [], []
    for img, lab, _ in loader:
        logps.append(model(img.to(device)).cpu())
        labs.append(lab)
    logp = torch.cat(logps)
    lab = torch.cat(labs)
    loss = float(F.nll_loss(logp, lab).item())
    pred = logp.argmax(dim=1).numpy()
    lab = lab.numpy()
    prob_di = torch.exp(logp[:, 1]).numpy()

    conf = np.zeros((3, 3), dtype=np.int64)
    for t in range(3):
        for p in range(3):
            conf[t, p] = int(np.sum((lab == t) & (pred == p)))

    is_di = lab == 1
    # Macro-averaged recall (balanced accuracy). This is the selection metric:
    # plain val loss picked epoch 12, where diphoton recall was still climbing
    # (it kept rising to epoch 35), because NLL is dominated by the classes the
    # network already handles. Macro recall weights all three classes equally,
    # which is what "do not throw away real diphotons" actually means.
    recalls = [conf[t, t] / conf[t].sum() if conf[t].sum() else 0.0 for t in range(3)]
    macro_recall = float(np.mean(recalls))
    return {
        "loss": loss,
        "macro_recall": macro_recall,
        "recalls": recalls,
        "acc": float(np.mean(pred == lab)),
        "conf": conf,
        # the number that matters: fraction of true diphotons scored above 0.9
        "di_recall_p9": float(np.mean(prob_di[is_di] > 0.9)) if is_di.any() else float("nan"),
        "di_recall_p5": float(np.mean(prob_di[is_di] > 0.5)) if is_di.any() else float("nan"),
        # and what it costs: hadronic clusters faking a diphoton
        "had_fake_p9": float(np.mean(prob_di[lab == 2] > 0.9)) if (lab == 2).any() else float("nan"),
    }


def print_confusion(conf):
    header = "true\\pred"
    print(f"    {header:>10s} " + " ".join(f"{n:>9s}" for n in NAMES))
    for t in range(3):
        row = conf[t]
        tot = int(row.sum())
        cells = " ".join(f"{row[p]:>9d}" for p in range(3))
        recall = f"   (recall {row[t]/tot:.3f})" if tot else "   (empty)"
        print(f"    {NAMES[t]:>10s} {cells}{recall}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--packed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=150)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--warm-start", default=None, help="classifier.onnx to start from")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--patience", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    # Line-buffer stdout: these runs are long and are watched through a log
    # file, where the default block buffering shows nothing for hours.
    sys.stdout.reconfigure(line_buffering=True)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)
    device = torch.device(args.device)
    print(f"[train] device={device}")

    tr = load_split(args.packed, "train")
    va = load_split(args.packed, "val")
    train_loader = DataLoader(TensorDataset(*tr), batch_size=args.batch_size,
                              shuffle=True, drop_last=True)
    val_loader = DataLoader(TensorDataset(*va), batch_size=1024)

    model = MLPhotonClassifier().to(device)
    if args.warm_start:
        load_onnx_weights(model, args.warm_start)
        print(f"[train] warm-started from {args.warm_start}")
        base = evaluate(model, val_loader, device)
        print(f"[train] starting point: loss {base['loss']:.4f} acc {base['acc']:.4f} "
              f"di_recall(>0.9) {base['di_recall_p9']:.4f}")
        print_confusion(base["conf"])

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, mode="max", factor=0.5,
                                                       patience=8)

    hist_path = os.path.join(args.out, "history.csv")
    hist = open(hist_path, "w", newline="")
    writer = csv.writer(hist)
    writer.writerow(["epoch", "train_loss", "val_loss", "acc", "macro_recall",
                     "di_recall_argmax", "di_recall_p9", "had_fake_p9", "lr", "sec"])

    best = -float("inf")   # selecting on macro recall now: higher is better
    since_best = 0
    for ep in range(1, args.epochs + 1):
        t0 = time.time()
        model.train()
        tot, nb = 0.0, 0
        for img, lab, _ in train_loader:
            img, lab = img.to(device), lab.to(device)
            opt.zero_grad()
            loss = F.nll_loss(model(img), lab)
            loss.backward()
            opt.step()
            tot += float(loss.item())
            nb += 1
        train_loss = tot / max(nb, 1)

        res = evaluate(model, val_loader, device)
        sched.step(res["macro_recall"])
        lr_now = opt.param_groups[0]["lr"]
        dt = time.time() - t0
        writer.writerow([ep, f"{train_loss:.6f}", f"{res['loss']:.6f}",
                         f"{res['acc']:.4f}", f"{res['macro_recall']:.4f}",
                         f"{res['recalls'][1]:.4f}", f"{res['di_recall_p9']:.4f}",
                         f"{res['had_fake_p9']:.5f}", f"{lr_now:.2e}", f"{dt:.1f}"])
        hist.flush()

        flag = ""
        if res["macro_recall"] > best:
            best = res["macro_recall"]
            since_best = 0
            torch.save(model.state_dict(), os.path.join(args.out, "classifier_best.pt"))
            flag = "  *best"
        else:
            since_best += 1

        if ep % 5 == 0 or flag or ep == 1:
            print(f"[ep {ep:4d}] train {train_loss:.4f}  val {res['loss']:.4f}  "
                  f"acc {res['acc']:.4f}  macro_rec {res['macro_recall']:.4f}  "
                  f"di_rec(argmax) {res['recalls'][1]:.4f}  "
                  f"di_rec(>0.9) {res['di_recall_p9']:.4f}  "
                  f"had_fake {res['had_fake_p9']:.5f}  {dt:.0f}s{flag}")

        if since_best >= args.patience:
            print(f"[train] no improvement for {args.patience} epochs; stopping at {ep}")
            break

    torch.save(model.state_dict(), os.path.join(args.out, "classifier_last.pt"))
    hist.close()

    model.load_state_dict(torch.load(os.path.join(args.out, "classifier_best.pt")))
    final = evaluate(model, val_loader, device)
    print(f"\n[train] best model: loss {final['loss']:.4f}  acc {final['acc']:.4f}  "
          f"macro_recall {final['macro_recall']:.4f}")
    print(f"[train] diphoton recall  >0.9: {final['di_recall_p9']:.4f}   "
          f">0.5: {final['di_recall_p5']:.4f}")
    print(f"[train] hadronic faking diphoton (>0.9): {final['had_fake_p9']:.5f}")
    print_confusion(final["conf"])
    print(f"[train] history -> {hist_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
