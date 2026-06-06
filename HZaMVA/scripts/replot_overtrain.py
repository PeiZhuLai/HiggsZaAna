#!/usr/bin/env python3
"""Regenerate the overtraining PDF for a trained NN with the capped-KS p-value
(train_nn.plot_overtraining / ks_capped). Use to fix plots whose KS collapsed to ~0 from the
large sample size. Runs in hza_NN. Rebuilds train + validation sets, scores, replots.
  python3 replot_overtrain.py --model model_Za_NN_run3_permass_lam100.pt \
      --out plots_nn/overtrain_model_Za_NN_run3_permass_lam100.pdf
"""
import argparse
import numpy as np, torch
from train_nn import ParamMLP, plot_overtraining, build_training_set

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-sig-per", type=int, default=30000)
    ap.add_argument("--n-bkg", type=int, default=350000)
    ap.add_argument("--ks-cap", type=int, default=4000)
    a = ap.parse_args()
    ck = torch.load(a.model, map_location="cpu", weights_only=False)
    net = ParamMLP(ck["nfeat"]); net.load_state_dict(ck["state_dict"]); net.eval()
    mu, sd = np.asarray(ck["mu"]), np.asarray(ck["sd"])
    print(f"[cfg] model={a.model} ks_cap={a.ks_cap}")
    X, y, w, *_ = build_training_set(a.n_sig_per, a.n_bkg)
    Xv, yv, wv, *_ = build_training_set(max(a.n_sig_per // 2, 1000), a.n_bkg // 2, tree="validation")
    p_sig, p_bk = plot_overtraining(net, mu, sd, (X, y, w), (Xv, yv, wv), a.out, device="cpu", ks_cap=a.ks_cap)
    print(f"[done] {a.out}  capped-KS p: sig={p_sig:.3g} bkg={p_bk:.3g}")

if __name__ == "__main__":
    main()
