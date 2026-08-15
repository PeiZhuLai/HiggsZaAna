"""Turn RecHit dumps into an MLPhoton training set.

Reads the ROOT files produced by RecHitDumper/test/closure_cfg.py, re-runs the
offline clustering (cluster_py, verified against the online C++ by
closure_test.py), truth-matches every cluster, and writes one npz per input:

    image      (N, 30, 30) float32   RAW crystal energies, exactly what
                                     Cluster::makeImage() produces
    eta        (N,)        float32   cluster eta (the regressor's 2nd input)
    energy     (N,)        float32   cluster total energy (the m/E denominator)
    label      (N,)        int8      0 = monophoton, 1 = diphoton, 2 = hadronic
                                     -- same order as the classifier's outputs
    moe_true   (N,)        float32   m_gen / E_cluster; NaN for non-diphoton
    m_gen      (N,)        float32   matched gen mass; NaN if unmatched
    dr_match   (N,)        float32   dR to the matched gen object
    n_crystals (N,)        int32
    ml_moe     (N,)        float32   the CURRENT model's prediction, kept so the
                                     retrained model can be compared per cluster

The image is stored UNNORMALISED on purpose: the classifier is trained on raw
energies and the regressor on the sum-normalised image, so both models can be
fed from one file (see cluster_py's module docstring).

Labelling
---------
For each cluster, in order:
  1. nearest gen ALP within --dr-alp            -> diphoton, target m_ALP/E
  2. nearest prompt gen photon NOT from an ALP  -> monophoton, target 0
     within --dr-photon
  3. otherwise                                  -> hadronic

Rule 2 excludes ALP daughters via gp_motherPdgId so a half-reconstructed a->gg
never leaks into the monophoton class.

A geometric match ALONE IS NOT ENOUGH for rules 1 and 2, and getting this wrong
silently destroys the regression target. The ALP deposits its energy over an
area, and the clustering routinely splits off fragments; a fragment sitting
within dR < 0.1 of the ALP passes the geometric cut but carries only part of the
energy. Its label would then pair the FULL gen mass with a PARTIAL cluster
energy, so the target m/E comes out far too large. (Measured on M1: labelling on
dR alone gave a median target m/E of 0.033 against the current model's 0.004 --
a factor 8, entirely an artefact of fragment clusters.)

So a matched cluster must also contain the gen object's energy, within
[--e-frac-lo, --e-frac-hi]. Clusters that match geometrically but fail the
energy containment are AMBIGUOUS: they are neither a clean diphoton/monophoton
nor genuinely hadronic, so they are dropped rather than dumped into the hadronic
class (which would teach the classifier that EM fragments are hadrons).

Input : <dump>/*.root  (TTree 'dumper/rechits')
Output: <outdir>/<stem>.npz  (+ a printed per-class tally)

Usage (CERN, env `hzgdna`):
    /eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python make_training_set.py \\
        --input '/eos/.../train_dump/mA_M0p5/*.root' \\
        --outdir /eos/.../train_npz/mA_M0p5
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_py import build_clusters  # noqa: E402

ALP_PDGID = 9000005
LABEL_MONO, LABEL_DI, LABEL_HAD = 0, 1, 2
LABEL_DROP = -1   # geometric match without energy containment -- unusable
LABEL_NAMES = {LABEL_MONO: "monophoton", LABEL_DI: "diphoton", LABEL_HAD: "hadronic"}

BRANCHES = ["run", "lumi", "event", "rh_ieta", "rh_iphi", "rh_eta", "rh_phi",
            "rh_energy", "pv_z",
            "gp_pdgId", "gp_pt", "gp_eta", "gp_phi", "gp_mass", "gp_energy",
            "gp_motherPdgId", "gp_isPromptFinalState",
            "ml_eta", "ml_phi", "ml_moe", "ml_diphotonScore"]


def delta_r(eta1, phi1, eta2, phi2):
    deta = np.asarray(eta1) - np.asarray(eta2)
    dphi = np.arctan2(np.sin(np.asarray(phi1) - np.asarray(phi2)),
                      np.cos(np.asarray(phi1) - np.asarray(phi2)))
    return np.hypot(deta, dphi)


def process_file(path, args):
    tree = uproot.open(path)[args.tree]
    missing = [b for b in BRANCHES if b not in tree]
    if missing:
        print(f"[skip] {path}: missing {missing}")
        return None

    stop = args.max_events if args.max_events > 0 else None
    data = tree.arrays(BRANCHES, entry_stop=stop, library="np")
    n_ev = len(data["event"])

    images, etas, energies, labels = [], [], [], []
    moes, mgens, drs, nxtals, mlmoes, efracs = [], [], [], [], [], []
    mldiphos = []
    n_dropped = [0]   # geometric match without energy containment

    for iev in range(n_ev):
        rh_e = np.asarray(data["rh_energy"][iev], dtype=np.float64)
        if not np.any(rh_e > 0.0):
            continue
        pv_z = float(data["pv_z"][iev])
        if not np.isfinite(pv_z):
            pv_z = 0.0

        clusters = build_clusters(data["rh_ieta"][iev], data["rh_iphi"][iev],
                                  data["rh_eta"][iev], data["rh_phi"][iev], rh_e)
        if not clusters:
            continue

        # ---- gen objects for this event -------------------------------
        gp_pdg = np.asarray(data["gp_pdgId"][iev])
        gp_eta = np.asarray(data["gp_eta"][iev], dtype=np.float64)
        gp_phi = np.asarray(data["gp_phi"][iev], dtype=np.float64)
        gp_mass = np.asarray(data["gp_mass"][iev], dtype=np.float64)
        gp_mom = np.asarray(data["gp_motherPdgId"][iev])
        gp_prompt = np.asarray(data["gp_isPromptFinalState"][iev])

        gp_energy = np.asarray(data["gp_energy"][iev], dtype=np.float64)

        is_alp = np.abs(gp_pdg) == ALP_PDGID
        # prompt photons that are NOT ALP daughters
        is_pho = (np.abs(gp_pdg) == 22) & (gp_prompt == 1) & (np.abs(gp_mom) != ALP_PDGID)

        alp_eta, alp_phi = gp_eta[is_alp], gp_phi[is_alp]
        alp_mass, alp_energy = gp_mass[is_alp], gp_energy[is_alp]
        pho_eta, pho_phi, pho_energy = gp_eta[is_pho], gp_phi[is_pho], gp_energy[is_pho]

        # Current-model predictions. The online MLPhoton collection is filled in
        # cluster order (closure_test.py verifies this holds event by event), so
        # index i lines up. Do NOT try to match on eta: ml_eta is the vertex-
        # corrected etaprime, not the cluster eta.
        ml_moe = np.asarray(data["ml_moe"][iev], dtype=np.float64)
        ml_dipho = np.asarray(data["ml_diphotonScore"][iev], dtype=np.float64)
        n_online = ml_moe.size

        for icl, c in enumerate(clusters):
            img = c.make_image()
            c_eta, c_phi = c.eta(), c.phi()
            energy = c.total_energy()
            if energy <= 0 or not np.isfinite(c_eta):
                continue
            if args.min_energy > 0 and energy < args.min_energy:
                continue

            label = LABEL_HAD
            moe_true = np.nan
            m_gen = np.nan
            dr_best = np.nan
            e_frac = np.nan

            if alp_eta.size:
                dr = delta_r(alp_eta, alp_phi, c_eta, c_phi)
                k = int(np.argmin(dr))
                if dr[k] < args.dr_alp:
                    dr_best = float(dr[k])
                    e_frac = energy / alp_energy[k] if alp_energy[k] > 0 else np.nan
                    if args.e_frac_lo <= e_frac <= args.e_frac_hi:
                        label = LABEL_DI
                        m_gen = float(alp_mass[k])
                        moe_true = m_gen / energy
                    else:
                        # near the ALP but not containing it: a fragment
                        label = LABEL_DROP

            if label == LABEL_HAD and pho_eta.size:
                dr = delta_r(pho_eta, pho_phi, c_eta, c_phi)
                k = int(np.argmin(dr))
                if dr[k] < args.dr_photon:
                    dr_best = float(dr[k])
                    e_frac = energy / pho_energy[k] if pho_energy[k] > 0 else np.nan
                    if args.e_frac_lo <= e_frac <= args.e_frac_hi:
                        label = LABEL_MONO
                        m_gen = 0.0
                        moe_true = 0.0
                    else:
                        label = LABEL_DROP

            if label == LABEL_DROP:
                n_dropped[0] += 1
                continue

            ml_val = float(ml_moe[icl]) if icl < n_online else np.nan
            ml_dip = float(ml_dipho[icl]) if icl < n_online else np.nan

            images.append(img.astype(np.float32))
            etas.append(np.float32(c_eta))
            energies.append(np.float32(energy))
            labels.append(np.int8(label))
            moes.append(np.float32(moe_true))
            mgens.append(np.float32(m_gen))
            drs.append(np.float32(dr_best))
            nxtals.append(np.int32(len(c.Es)))
            mlmoes.append(np.float32(ml_val))
            mldiphos.append(np.float32(ml_dip))
            efracs.append(np.float32(e_frac))

    if not images:
        print(f"[warn] {path}: no clusters survived")
        return None

    return {
        "image": np.stack(images),
        "eta": np.asarray(etas, dtype=np.float32),
        "energy": np.asarray(energies, dtype=np.float32),
        "label": np.asarray(labels, dtype=np.int8),
        "moe_true": np.asarray(moes, dtype=np.float32),
        "m_gen": np.asarray(mgens, dtype=np.float32),
        "dr_match": np.asarray(drs, dtype=np.float32),
        "e_frac": np.asarray(efracs, dtype=np.float32),
        "n_crystals": np.asarray(nxtals, dtype=np.int32),
        "ml_moe": np.asarray(mlmoes, dtype=np.float32),
        "ml_diphotonScore": np.asarray(mldiphos, dtype=np.float32),
        "n_dropped": np.int64(n_dropped[0]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="glob of dump ROOT files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tree", default="dumper/rechits")
    ap.add_argument("--dr-alp", type=float, default=0.10,
                    help="cluster-to-gen-ALP match radius (diphoton label)")
    ap.add_argument("--dr-photon", type=float, default=0.10,
                    help="cluster-to-prompt-photon match radius (monophoton label)")
    ap.add_argument("--e-frac-lo", type=float, default=0.7,
                    help="min E_cluster / E_gen for a match to count. Without "
                         "this, fragment clusters pair the full gen mass with a "
                         "partial energy and the regression target is wrong.")
    ap.add_argument("--e-frac-hi", type=float, default=1.3,
                    help="max E_cluster / E_gen (above it the cluster has "
                         "swallowed unrelated energy)")
    ap.add_argument("--min-energy", type=float, default=1.0,
                    help="drop clusters below this total energy (GeV); soft "
                         "pileup blobs are not what either network is for")
    ap.add_argument("--max-events", type=int, default=-1, help="per file; -1 = all")
    args = ap.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        print(f"[FATAL] no files match {args.input}")
        return 2
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[make] {len(files)} input files -> {args.outdir}")
    tally = {LABEL_MONO: 0, LABEL_DI: 0, LABEL_HAD: 0}
    n_written = 0
    n_dropped_total = 0

    for i, path in enumerate(files):
        out = os.path.join(args.outdir,
                           os.path.splitext(os.path.basename(path))[0] + ".npz")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            with np.load(out) as z:
                for lab in tally:
                    tally[lab] += int(np.sum(z["label"] == lab))
                if "n_dropped" in z:
                    n_dropped_total += int(z["n_dropped"])
            n_written += 1
            print(f"  [{i+1}/{len(files)}] {os.path.basename(path)} -> exists, skipped")
            continue

        res = process_file(path, args)
        if res is None:
            continue
        np.savez_compressed(out, **res)
        n_written += 1
        for lab in tally:
            tally[lab] += int(np.sum(res["label"] == lab))
        n_dropped_total += int(res["n_dropped"])
        print(f"  [{i+1}/{len(files)}] {os.path.basename(path)} -> "
              f"{len(res['label'])} clusters ({int(res['n_dropped'])} dropped)")

    total = sum(tally.values())
    print()
    print(f"[make] wrote {n_written} npz, {total} clusters total")
    for lab, n in tally.items():
        frac = 100.0 * n / total if total else 0.0
        print(f"        {LABEL_NAMES[lab]:>11s}: {n:>9d}  ({frac:5.1f}%)")
    print(f"        {'dropped':>11s}: {n_dropped_total:>9d}  "
          f"(matched geometrically but failed E containment)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
