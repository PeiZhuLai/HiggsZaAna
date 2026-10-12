#!/usr/bin/env python3
"""m_a = 1 GeV: do overlapping photon showers exist in Run 2 simulation as in Run 3?

Same generator-level definition as Plot/scripts/plot_mA1_doublepeak.py (gen dR(gg) < 0.05,
about three ECAL crystals), applied directly to NanoAOD so that Run 2 (UL NanoAODv9) and
Run 3 (NanoAODv15) are compared with identical code and no analysis selection:
  - each ALP photon is matched to the closest reco photon (pT > 10 GeV, |eta| < 2.5) within dR < 0.1;
  - "resolved": the two ALP photons match two different reco photons;
  - "merged":   both match the same reco photon;
  - for resolved events: reco m_gg and photon-pair pT response (reco / gen).
Optionally both reco photons must pass the MVA ID WP90.

Usage: overlap_showers_nano.py --tag run2_UL18 --out out.json --max-events 150000 FILE [FILE ...]
"""

import argparse
import json

import awkward as ak
import numpy as np
import uproot

ALP_IDS = (36, 9000005)
BRANCHES = ["GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi",
            "Photon_pt", "Photon_eta", "Photon_phi", "Photon_mvaID_WP90"]


def dphi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def process(arr):
    pid, mo = arr.GenPart_pdgId, arr.GenPart_genPartIdxMother
    mo_pid = ak.where(mo >= 0, pid[ak.where(mo >= 0, mo, 0)], 0)
    is_alp_gamma = (pid == 22) & ((mo_pid == ALP_IDS[0]) | (mo_pid == ALP_IDS[1]))
    sel = ak.sum(is_alp_gamma, axis=1) == 2
    arr, is_alp_gamma = arr[sel], is_alp_gamma[sel]
    g_pt = ak.to_numpy(arr.GenPart_pt[is_alp_gamma])
    g_eta = ak.to_numpy(arr.GenPart_eta[is_alp_gamma])
    g_phi = ak.to_numpy(arr.GenPart_phi[is_alp_gamma])
    gen_dr = np.hypot(g_eta[:, 0] - g_eta[:, 1], dphi(g_phi[:, 0], g_phi[:, 1]))
    gpx = (g_pt * np.cos(g_phi)).sum(1); gpy = (g_pt * np.sin(g_phi)).sum(1)
    gen_pair_pt = np.hypot(gpx, gpy)

    ph_ok = (arr.Photon_pt > 10) & (abs(arr.Photon_eta) < 2.5)
    out = {"gen_dr": gen_dr, "gen_pair_pt": gen_pair_pt}
    idx = []
    for k in range(2):
        dr = np.sqrt((arr.Photon_eta - g_eta[:, k]) ** 2 + dphi(arr.Photon_phi, g_phi[:, k]) ** 2)
        dr = ak.where(ph_ok, dr, 999.0)
        best = ak.argmin(dr, axis=1, keepdims=True)
        bdr = ak.fill_none(ak.firsts(dr[best]), 999.0)
        idx.append(np.where(ak.to_numpy(bdr) < 0.1, ak.to_numpy(ak.fill_none(ak.firsts(best), -1)), -1))
    i0, i1 = idx
    both = (i0 >= 0) & (i1 >= 0)
    out["resolved"] = both & (i0 != i1)
    out["merged"] = both & (i0 == i1)
    j0, j1 = np.where(i0 >= 0, i0, 0), np.where(i1 >= 0, i1, 0)

    def take(branch, j):
        return ak.to_numpy(ak.fill_none(ak.firsts(arr[branch][ak.local_index(arr[branch]) == j]), 0))
    pt0, pt1 = take("Photon_pt", j0), take("Photon_pt", j1)
    eta0, eta1 = take("Photon_eta", j0), take("Photon_eta", j1)
    phi0, phi1 = take("Photon_phi", j0), take("Photon_phi", j1)
    id0, id1 = take("Photon_mvaID_WP90", j0), take("Photon_mvaID_WP90", j1)
    m2 = 2 * pt0 * pt1 * (np.cosh(eta0 - eta1) - np.cos(phi0 - phi1))
    out["mgg"] = np.sqrt(np.clip(m2, 0, None))
    rpx = pt0 * np.cos(phi0) + pt1 * np.cos(phi1); rpy = pt0 * np.sin(phi0) + pt1 * np.sin(phi1)
    out["resp"] = np.hypot(rpx, rpy) / gen_pair_pt
    out["id90"] = (id0 > 0) & (id1 > 0)
    return out, int(ak.num(sel, axis=0))


def summarize(d, n_read):
    res = {"n_read": n_read, "n_two_alp_photons": int(len(d["gen_dr"]))}
    close = d["gen_dr"] < 0.05
    res["frac_gen_dR_lt_0p05"] = float(close.mean())
    for name, gsel in (("gen_dR_lt_0p05", close), ("gen_dR_ge_0p05", ~close)):
        r = {"n": int(gsel.sum()),
             "frac_resolved": float(d["resolved"][gsel].mean()),
             "frac_merged_one_photon": float(d["merged"][gsel].mean())}
        for idname, idsel in (("noID", np.ones_like(gsel)), ("WP90", d["id90"])):
            s = gsel & d["resolved"] & idsel
            r[idname] = {"n_resolved": int(s.sum()),
                         "frac_of_class": float(s.sum() / max(gsel.sum(), 1)),
                         "median_resp": float(np.median(d["resp"][s])) if s.any() else None,
                         "frac_mgg_ge_1p15": float((d["mgg"][s] >= 1.15).mean()) if s.any() else None}
        res[name] = r
    # among resolved events (what an analysis can select), the share of overlapping showers
    for idname, idsel in (("noID", np.ones_like(close)), ("WP90", d["id90"])):
        s = d["resolved"] & idsel
        res["resolved_" + idname] = {
            "n": int(s.sum()),
            "share_gen_dR_lt_0p05": float(close[s].mean()),
            "share_mgg_ge_1p15": float((d["mgg"][s] >= 1.15).mean()),
            "median_resp": float(np.median(d["resp"][s]))}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-events", type=int, default=150000)
    a = ap.parse_args()
    parts, n_read = [], 0
    for f in a.files:
        if n_read >= a.max_events:
            break
        with uproot.open(f, timeout=600) as fh:
            t = fh["Events"]
            stop = min(t.num_entries, a.max_events - n_read)
            for arr in t.iterate(BRANCHES, entry_stop=stop, step_size=50000):
                o, n = process(arr)
                parts.append(o); n_read += n
                print(a.tag, "read", n_read, flush=True)
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    res = {"tag": a.tag, "files": a.files, **summarize(d, n_read)}
    np.savez_compressed(a.out.replace(".json", ".npz"), **d)
    with open(a.out, "w") as fo:
        json.dump(res, fo, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
