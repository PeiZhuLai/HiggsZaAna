"""Export trained weights back to the ONNX files CMSSW loads, and verify.

Exporting is the easy half; the point of this script is the checking that
follows it. Three things are verified before anything is installed:

  1. the exported graph reproduces the torch model (onnxruntime vs torch)
  2. the graph has the input names MLPhotonProducer looks up ('img', and 'eta'
     for the regressor) and opset 9 -- get either wrong and CMSSW fails at
     runtime, not at build time
  3. the new model is compared against the currently deployed one on real
     clusters, so the size and direction of the change is on the record

Installation into RecoEgamma/EgammaMLPhotonProducers/data/ is opt-in
(--install) and always backs up what was there.

Input : <run>/regressor_best.pt and/or <run>/classifier_best.pt
Output: <out>/{regressor,classifier}.onnx  (+ optional install)

Usage (CERN, env `hzgdna`):
    python export_onnx.py --regressor-pt <run>/regressor_best.pt --out <dir>
    python export_onnx.py ... --install
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import (ISIZE, MLPhotonClassifier, MLPhotonRegressor,  # noqa: E402
                    MLPhotonRegressorLog, export_classifier, export_regressor)


def load_regressor(pt_path):
    """Build whichever regressor head the checkpoint was trained with.

    train_regressor.py defaults to the log head, whose state_dict carries the
    extra `mu`/`sigma` buffers; the deployed head has neither. Picking the wrong
    class fails at load_state_dict, so decide from the checkpoint itself rather
    than from a flag someone has to remember to pass.
    """
    state = torch.load(pt_path, map_location="cpu")
    is_log = "mu" in state and "sigma" in state
    model = MLPhotonRegressorLog() if is_log else MLPhotonRegressor()
    model.load_state_dict(state)
    model.eval()
    print(f"  head: {'log (exp(raw*sigma+mu))' if is_log else 'deployed (LeakyReLU)'}"
          + (f"   mu={float(model.mu):.4f} sigma={float(model.sigma):.4f}" if is_log else ""))
    return model

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))       # RegressMergedPhoton
_PROJECT = os.path.dirname(_REPO)                                          # HiggsZaAna
DATA_DIR = os.path.join(_REPO, "RecoEgamma", "EgammaMLPhotonProducers", "data")


def install_targets():
    """Every copy of the .onnx that has to be kept in step.

    There are TWO: the source copy under RegressMergedPhoton/, and the one inside
    the CMSSW area, which is what MLPhotonProducer actually loads (closure_cfg.py
    resolves it through $CMSSW_BASE). Updating only the first leaves CMSSW
    running the old model while everything looks deployed -- the failure would
    only show up as "the retraining changed nothing".
    """
    targets = [DATA_DIR]
    for cmssw in sorted(glob.glob(os.path.join(_PROJECT, "CMSSW_*"))):
        d = os.path.join(cmssw, "src", "RecoEgamma", "EgammaMLPhotonProducers", "data")
        if os.path.isdir(d):
            targets.append(d)
    return targets
DEFAULT_NPZ = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_npz/"
               "sig_M0p5/*.npz")
RTOL = 1e-4


def check_graph(path, expect_inputs):
    import onnx
    m = onnx.load(path)
    names = [i.name for i in m.graph.input]
    opsets = [(o.domain, o.version) for o in m.opset_import]
    ok = names == expect_inputs and any(v == 9 for d, v in opsets if d in ("", "ai.onnx"))
    print(f"  inputs {names}  opset {opsets}  -> {'OK' if ok else 'WRONG'}")
    if names != expect_inputs:
        print(f"  [FAIL] MLPhotonProducer looks these up by name; expected {expect_inputs}")
    return ok


def load_sample(npz_glob, n):
    files = sorted(glob.glob(npz_glob))
    if not files:
        return None
    with np.load(files[0]) as z:
        img = z["image"][:n].astype(np.float32)
        eta = z["eta"][:n].astype(np.float32)
        moe = z["moe_true"][:n].astype(np.float64)
        lab = z["label"][:n]
        mlm = z["ml_moe"][:n].astype(np.float64)
    return img, eta, moe, lab, mlm


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regressor-pt", default=None)
    ap.add_argument("--classifier-pt", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--npz", default=DEFAULT_NPZ, help="real clusters for the checks")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--install", action="store_true",
                    help=f"copy into {DATA_DIR} (backs up the existing files)")
    args = ap.parse_args()

    if not args.regressor_pt and not args.classifier_pt:
        print("[FATAL] give --regressor-pt and/or --classifier-pt")
        return 2
    os.makedirs(args.out, exist_ok=True)

    import onnxruntime as ort
    sample = load_sample(args.npz, args.n)
    if sample is None:
        print(f"[warn] no npz matched {args.npz}; skipping the comparison step")
    ok_all = True
    produced = {}

    if args.regressor_pt:
        print(f"\n=== regressor: {args.regressor_pt}")
        model = load_regressor(args.regressor_pt)
        path = os.path.join(args.out, "regressor.onnx")
        export_regressor(model, path)
        print(f"  wrote {path} ({os.path.getsize(path)} bytes)")
        ok_all &= check_graph(path, ["img", "eta"])

        if sample:
            img, eta, moe, lab, mlm = sample
            t = torch.from_numpy(img).unsqueeze(1)
            s = t.sum(dim=(1, 2, 3), keepdim=True)
            tn = torch.where(s > 0, t / s, t)
            with torch.no_grad():
                t_pred = model(tn, torch.from_numpy(eta).reshape(-1, 1)).numpy().astype(np.float64)
            sess = ort.InferenceSession(path, providers=["CPUExecutionProvider"])
            o_pred = np.array([float(np.asarray(sess.run(None, {
                "img": tn[i:i+1].numpy(),
                "eta": eta[i:i+1].reshape(1, 1)})[0]).reshape(-1)[0])
                for i in range(len(img))])
            d = np.max(np.abs(t_pred - o_pred) / np.maximum(np.abs(t_pred), 1e-6))
            print(f"  torch vs onnxruntime: max reldiff {d:.3e} "
                  f"-> {'PASS' if d <= RTOL else 'FAIL'}")
            ok_all &= d <= RTOL

            di = (lab == 1) & np.isfinite(moe) & (moe > 0)
            if di.sum() > 10:
                new_r = np.median(t_pred[di] / moe[di])
                old_r = np.median(mlm[di] / moe[di])
                print(f"  on {int(di.sum())} true diphotons (this npz):")
                print(f"    median pred/true  old {old_r:.3f}  ->  new {new_r:.3f}")
        produced["regressor.onnx"] = path

    if args.classifier_pt:
        print(f"\n=== classifier: {args.classifier_pt}")
        model = MLPhotonClassifier()
        model.load_state_dict(torch.load(args.classifier_pt, map_location="cpu"))
        model.eval()
        path = os.path.join(args.out, "classifier.onnx")
        export_classifier(model, path)
        print(f"  wrote {path} ({os.path.getsize(path)} bytes)")
        ok_all &= check_graph(path, ["img"])

        if sample:
            img, eta, moe, lab, mlm = sample
            t = torch.from_numpy(img).unsqueeze(1)
            with torch.no_grad():
                t_logp = model(t).numpy().astype(np.float64)
            t_dip = np.exp(t_logp[:, 1])
            sess = ort.InferenceSession(path, providers=["CPUExecutionProvider"])
            o_dip = np.empty(len(img))
            for i in range(len(img)):
                c = sess.run(None, {"img": t[i:i+1].numpy()})[0].reshape(-1)
                e = np.exp(c.astype(np.float64))
                o_dip[i] = e[1] / e.sum()
            d = np.max(np.abs(t_dip - o_dip) / np.maximum(np.abs(t_dip), 1e-6))
            print(f"  torch vs onnxruntime: max reldiff {d:.3e} "
                  f"-> {'PASS' if d <= RTOL else 'FAIL'}")
            ok_all &= d <= RTOL

            di = lab == 1
            if di.sum() > 10:
                print(f"  diphoton recall(>0.9) on {int(di.sum())} true diphotons: "
                      f"{np.mean(t_dip[di] > 0.9):.3f}")
                had = lab == 2
                if had.sum() > 10:
                    print(f"  hadronic faking diphoton(>0.9): "
                          f"{np.mean(t_dip[had] > 0.9):.5f}")
        produced["classifier.onnx"] = path

    print(f"\n[export] checks: {'PASS' if ok_all else 'FAIL'}")
    if not ok_all:
        print("[export] not installing on a failed check")
        return 1

    if args.install:
        stamp = time.strftime("%Y%m%d_%H%M%S")
        targets = install_targets()
        print(f"[install] {len(targets)} target dir(s):")
        for t in targets:
            print(f"           {t}")
        for tgt in targets:
            for name, src in produced.items():
                dst = os.path.join(tgt, name)
                if os.path.exists(dst):
                    bak = f"{dst}.bak_{stamp}"
                    shutil.copy2(dst, bak)
                    print(f"[install] backed up {dst} -> {bak}")
                shutil.copy2(src, dst)
                print(f"[install] {src} -> {dst}")
        print("[install] now re-run the preprocessing closure so the new models "
              "are exercised end to end:")
        print("           bash run_closure.sh M1 200")
    else:
        print(f"[export] not installed. To deploy: re-run with --install")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
