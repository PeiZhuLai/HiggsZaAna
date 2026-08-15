"""Model-level closure: the PyTorch rebuild must reproduce the deployed ONNX.

closure_test.py proved the offline PREPROCESSING matches the online C++. This
script proves the offline MODEL matches: it loads the deployed weights into the
torch modules from models.py and checks, on real cluster images from the
training npz, that

    torch(image, eta)  ==  onnxruntime(image, eta)  ==  ml_moe stored in the npz

If this passes, the torch definition is structurally exact, so a model retrained
in torch and exported back to ONNX will drop into CMSSW unchanged.

Input : one or more training npz (default: a signal tag)
Output: stdout report, exit 0 on PASS

Usage (CERN, env `hzgdna`):
    /eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python verify_models.py
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import (MLPhotonClassifier, MLPhotonRegressor,  # noqa: E402
                    export_classifier, export_regressor, load_onnx_weights)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "RecoEgamma", "EgammaMLPhotonProducers", "data")
DEFAULT_NPZ = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_npz/"
               "sig_M0p5/*.npz")

RTOL = 1e-4


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=DEFAULT_NPZ)
    ap.add_argument("--classifier", default=os.path.join(DATA_DIR, "classifier.onnx"))
    ap.add_argument("--regressor", default=os.path.join(DATA_DIR, "regressor.onnx"))
    ap.add_argument("--n", type=int, default=2000, help="clusters to test")
    ap.add_argument("--check-export", action="store_true", default=True,
                    help="also round-trip torch -> onnx -> onnxruntime")
    args = ap.parse_args()

    files = sorted(glob.glob(args.npz))
    if not files:
        print(f"[FATAL] no npz match {args.npz}")
        return 2

    with np.load(files[0]) as z:
        img = z["image"][:args.n].astype(np.float32)
        eta = z["eta"][:args.n].astype(np.float32)
        ml_moe = z["ml_moe"][:args.n].astype(np.float64)
        ml_dip = z["ml_diphotonScore"][:args.n].astype(np.float64)
    n = len(img)
    print(f"[verify] {os.path.basename(files[0])}: {n} clusters")

    # ---- torch, deployed weights -------------------------------------
    cls = MLPhotonClassifier()
    reg = MLPhotonRegressor()
    print(f"[verify] classifier weights: {len(load_onnx_weights(cls, args.classifier))} tensors")
    print(f"[verify] regressor  weights: {len(load_onnx_weights(reg, args.regressor))} tensors")
    cls.eval()
    reg.eval()

    img_t = torch.from_numpy(img).unsqueeze(1)                 # (N,1,30,30) raw
    sums = img_t.sum(dim=(1, 2, 3), keepdim=True)
    img_norm = torch.where(sums > 0, img_t / sums, img_t)      # regressor input
    eta_t = torch.from_numpy(eta).reshape(-1, 1)

    with torch.no_grad():
        t_moe = reg(img_norm, eta_t).numpy().astype(np.float64)
        t_logp = cls(img_t).numpy().astype(np.float64)
    t_dip = np.exp(t_logp[:, 1])

    # ---- onnxruntime on the same inputs ------------------------------
    import onnxruntime as ort
    s_reg = ort.InferenceSession(args.regressor, providers=["CPUExecutionProvider"])
    s_cls = ort.InferenceSession(args.classifier, providers=["CPUExecutionProvider"])
    o_moe = np.empty(n)
    o_dip = np.empty(n)
    for i in range(n):
        r = s_reg.run(None, {"img": img_norm[i:i+1].numpy(),
                             "eta": eta[i:i+1].reshape(1, 1)})[0]
        o_moe[i] = float(np.asarray(r).reshape(-1)[0])
        c = s_cls.run(None, {"img": img_t[i:i+1].numpy()})[0].reshape(-1)
        e = np.exp(c.astype(np.float64))
        o_dip[i] = e[1] / e.sum()

    def report(name, a, b):
        denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-6)
        d = np.abs(a - b) / denom
        ok = np.nanmax(d) <= RTOL
        print(f"  {name:<34s} max reldiff {np.nanmax(d):.3e}   "
              f"{'PASS' if ok else 'FAIL'}")
        return ok

    print("\n[verify] torch vs onnxruntime (same inputs)")
    ok = report("regressor m/E", t_moe, o_moe)
    ok &= report("classifier diphotonScore", t_dip, o_dip)

    print("\n[verify] torch vs values stored by the C++ producer")
    good = np.isfinite(ml_moe)
    ok &= report("regressor m/E vs ml_moe", t_moe[good], ml_moe[good])
    ok &= report("classifier vs ml_diphotonScore", t_dip[good], ml_dip[good])

    # ---- export round-trip -------------------------------------------
    if args.check_export:
        import tempfile
        print("\n[verify] torch -> onnx (opset 9) -> onnxruntime round-trip")
        with tempfile.TemporaryDirectory() as td:
            rp = os.path.join(td, "reg.onnx")
            cp = os.path.join(td, "cls.onnx")
            export_regressor(reg, rp)
            export_classifier(cls, cp)
            sr = ort.InferenceSession(rp, providers=["CPUExecutionProvider"])
            sc = ort.InferenceSession(cp, providers=["CPUExecutionProvider"])
            e_moe = np.empty(n)
            e_dip = np.empty(n)
            for i in range(n):
                e_moe[i] = float(np.asarray(sr.run(None, {
                    "img": img_norm[i:i+1].numpy(),
                    "eta": eta[i:i+1].reshape(1, 1)})[0]).reshape(-1)[0])
                c = sc.run(None, {"img": img_t[i:i+1].numpy()})[0].reshape(-1)
                ex = np.exp(c.astype(np.float64))
                e_dip[i] = ex[1] / ex.sum()
            ok &= report("re-exported regressor", t_moe, e_moe)
            ok &= report("re-exported classifier", t_dip, e_dip)

    print(f"\n[verify] RESULT: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("[verify] The torch definition does not match the deployed graph. "
              "Fix models.py before training -- otherwise the retrained weights "
              "will not mean the same thing in CMSSW.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
