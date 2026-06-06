#!/usr/bin/env python3
"""Re-evaluate trained NN models' 125-peak sculpting R over a fine mass grid (mA=1..5,10,20),
to see which masses per-mass DisCo actually helps. Runs in hza_NN (torch + uproot)."""
import numpy as np, torch
from train_nn import BASE_VARS, ParamMLP, _feat, _load, ROOT_DIR

MASSES = [1, 2, 3, 4, 5, 10, 20]

def load_model(path):
    ck = torch.load(path, map_location="cpu", weights_only=False)
    net = ParamMLP(ck["nfeat"]); net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["mu"]), np.asarray(ck["sd"])

bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
Hm = bk["H_m"]; wB = np.clip(bk["factor"], 0, None)
peak = (Hm > 120) & (Hm < 130); frac = wB[peak].sum() / wB.sum()
def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

def R(net, mu, sd, mh, eff=0.002):
    X = (_feat(bk, mh) - mu) / sd
    with torch.no_grad():
        s = torch.sigmoid(net(torch.tensor(X, dtype=torch.float32))).numpy()
    thr = wq(s, wB, 1 - eff); p = s > thr
    return (wB[p & peak].sum() / wB[p].sum()) / frac if wB[p].sum() > 0 else float("nan")

models = [("lam0 (baseline)", "model_Za_NN_run3_lam0.pt"),
          ("per-mass lam100", "model_Za_NN_run3_permass_lam100.pt")]
res = {}
for tag, path in models:
    net, mu, sd = load_model(path)
    res[tag] = {m: R(net, mu, sd, m) for m in MASSES}
    print(f"{tag:18s} " + "  ".join(f"mA{m}={res[tag][m]:.2f}" for m in MASSES))

print("\nR reduction (baseline - per-mass), >0 means DisCo helped:")
print("        " + "  ".join(f"mA{m}={res['lam0 (baseline)'][m]-res['per-mass lam100'][m]:+.2f}" for m in MASSES))
