#!/usr/bin/env python3
"""Skeleton: apply a trained NN (train_nn.py output) to a ROOT file and write per-mass
NN score branches -- the NN analogue of Parque2Root_BDT.add_mva_scores (MVA_Score_mA_M{ma}).

For each mass hypothesis ma in 1..30: param=(ALP_m-ma)/H_m, forward the net, write a
`NN_Score_mA_M{ma}` branch. Runs in hza_NN (torch + uproot). flashgg (CMSSW, no torch)
then only reads these branches -- nothing downstream needs torch.

SKELETON: single-file path is functional; the per-sample/year batch loop and a condor
driver (mirroring Parquet2Rootfile/Condor) are left as TODOs.
"""
from __future__ import annotations
import argparse
import numpy as np, uproot, awkward as ak, torch
from train_nn import BASE_VARS, ParamMLP, _feat   # reuse the exact feature definitions

READ = BASE_VARS + ["H_m", "ALP_m"]   # branches needed to rebuild the 16 features per mass
SLIM_KEEP = ["H_m", "ALP_m", "factor", "event", "Z_m"]   # kept in --slim outputs (for diagnostics)

def load_model(path, device):
    ck = torch.load(path, map_location=device, weights_only=False)
    net = ParamMLP(ck["nfeat"]).to(device)
    net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["mu"]), np.asarray(ck["sd"])

def score_arrays(arrays, net, mu, sd, device, masses=range(1, 31)):
    """Return {NN_Score_mA_M{ma}: array} for every mass hypothesis."""
    out = {}
    for ma in masses:
        X = _feat(arrays, ma)                      # (N,16); param = (ALP_m-ma)/H_m inside _feat
        with torch.no_grad():
            s = torch.sigmoid(net(torch.tensor((X-mu)/sd, dtype=torch.float32, device=device))).cpu().numpy()
        out[f"NN_Score_mA_M{ma}"] = s.astype(np.float32)
    return out

def apply_to_file(in_path, out_path, net, mu, sd, device,
                  trees=("inclusive", "train", "validation", "test"), slim=False):
    """Write NN_Score_mA_M* per tree. slim=False copies every branch (full analysis file);
    slim=True writes only SLIM_KEEP + the scores (small friend-style file for diagnostics)."""
    fin = uproot.open(in_path)
    present = {k.split(";")[0] for k in fin.keys()}
    fout = uproot.recreate(out_path)
    for tree in trees:
        if tree not in present:
            continue
        if slim:
            # Flat branches only -> numpy round-trip is enough (no jagged collections); fast + small.
            keep = [b for b in dict.fromkeys(READ + SLIM_KEEP) if b in fin[tree].keys()]
            np_in = fin[tree].arrays(keep, library="np")
            missing = [b for b in READ if b not in np_in]
            if missing:
                raise KeyError(f"{in_path}:{tree} missing branches for scoring: {missing}")
            out = {b: np_in[b] for b in dict.fromkeys(SLIM_KEEP) if b in np_in}
            out.update(score_arrays(np_in, net, mu, sd, device))
            fout[tree] = out
            print(f"[apply] {out_path}:{tree}  rows={len(np_in['H_m'])}  slim(+30 NN_Score)")
        else:
            # Read as awkward so jagged branches (photon/jet collections) survive the round-trip;
            # uproot's library="np" turns them into object dtype which awkward cannot write back.
            arr = fin[tree].arrays()
            missing = [b for b in READ if b not in arr.fields]
            if missing:
                raise KeyError(f"{in_path}:{tree} missing branches for scoring: {missing}")
            np_in = {b: ak.to_numpy(arr[b]) for b in READ}          # flat feature branches -> numpy
            for name, vals in score_arrays(np_in, net, mu, sd, device).items():
                arr = ak.with_field(arr, vals, name)                # append NN_Score_mA_M* branches
            fout[tree] = arr                        # keeps inclusive/train/validation/test split
            print(f"[apply] {out_path}:{tree}  rows={len(np_in['H_m'])}  full(+30 NN_Score)")
    # TODO(memory): for very large files, score in chunks (uproot iterate) instead of loading all.

# TODO(batch): wrap this in a driver over every sample/year of run3_bdt_scored_nominal and a
# condor submit (mirror Parquet2Rootfile/Condor/{1_make_joblist.py,2_submit.sub}).

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="model_Za_NN_run3.pt")
    ap.add_argument("--input", required=True, help="input ROOT (with _oHm + H_m/ALP_m branches)")
    ap.add_argument("--output", required=True, help="output ROOT (original + NN_Score_mA_M*)")
    ap.add_argument("--slim", action="store_true",
                    help="write only H_m/ALP_m/factor/event/Z_m + NN_Score (small, for diagnostics)")
    ap.add_argument("--device", default=("cuda" if torch.cuda.is_available() else "cpu"))
    a = ap.parse_args()
    print(f"[cfg] model={a.model} device={a.device} slim={a.slim}")
    net, mu, sd = load_model(a.model, a.device)
    apply_to_file(a.input, a.output, net, mu, sd, a.device, slim=a.slim)
    print("[done]", a.output)

if __name__ == "__main__":
    main()
