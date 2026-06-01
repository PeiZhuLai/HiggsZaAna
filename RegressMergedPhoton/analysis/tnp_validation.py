"""
Z→ee tag-and-probe ONNX validation for the MLPhoton regressor.

Follows AN-20-142 §5.2-5.3:
  1) Open DY→ee MLNanoAOD files produced by Prod_MLNanoAOD_run3_mc.py.
  2) Select tag electrons by Electron_cutBased >= 4 (tight ID).
  3) For each event with a tight tag, look for the *other* electron probe.
  4) Require the dielectron invariant mass in (60, 120) GeV (Z peak).
  5) Truth-match the probe to the nearest MLPhoton (ΔR<0.05).
  6) Extract MLPhoton_mass per |η| bin (central <0.5, middle <1.0, forward <1.44).

For data-vs-MC validation (later, when we have Run3 EGamma data MLNanoAOD):
  7) chi² scan over (s_scale, s_smear) that, when applied to MC m_Γ via
     m_Γ → s_scale × Gauss(m_Γ, s_smear), minimises the χ² versus data.

This script implements steps 1-6 and prints/plots the MC-only spectrum.
The data-vs-MC scan is a small follow-up that re-uses the same loader.

Usage:
    python tnp_validation.py <glob_or_dir> [--output out.png] [--max-events N]
"""

from __future__ import annotations

import argparse
import glob
import os
import awkward as ak
import numpy as np
import uproot

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


ELECTRON_TIGHT_BIT = 4  # Electron_cutBased: 0 fail, 1 veto, 2 loose, 3 medium, 4 tight
ZMASS_LO, ZMASS_HI = 60.0, 120.0
ETA_BINS = [(0.0, 0.5, "central"), (0.5, 1.0, "middle"), (1.0, 1.44, "forward")]
MLPHOTON_MATCH_DR = 0.05  # AN §5.2 uses ΔR<0.04
PROBE_ID_MIN_PT = 7.0


def load_branches(file_paths: list[str], max_events: int | None = None):
    """Stream-read the minimal set of branches needed for TnP."""
    branches = [
        "nElectron",
        "Electron_pt", "Electron_eta", "Electron_phi", "Electron_charge",
        "Electron_cutBased",
        "nMLPhoton",
        "MLPhoton_pt", "MLPhoton_eta", "MLPhoton_phi", "MLPhoton_mass",
        "MLPhoton_diphotonScore",
    ]
    total = 0
    chunks: list[ak.Array] = []
    for fp in file_paths:
        try:
            t = uproot.open(fp)["Events"]
        except Exception as e:
            print(f"[load] cannot open {fp}: {e}")
            continue
        avail = [b for b in branches if b in t]
        missing = set(branches) - set(avail)
        if missing:
            print(f"[load] {fp}: missing {missing}, skipping")
            continue
        arr = t.arrays(avail, library="ak")
        chunks.append(arr)
        total += len(arr)
        print(f"[load] {fp}: {len(arr)} events (running total {total})")
        if max_events is not None and total >= max_events:
            break
    if not chunks:
        raise RuntimeError("no readable input files")
    return ak.concatenate(chunks, axis=0)


def tag_probe(events: ak.Array) -> dict[str, np.ndarray]:
    """Find Z→ee pairs and truth-match the probe to MLPhoton.

    Returns a dict of numpy arrays, one row per accepted probe:
        probe_pt, probe_eta, probe_phi, mll, mlphoton_mass, mlphoton_dr
    """
    ele = ak.zip({
        "pt": events.Electron_pt, "eta": events.Electron_eta,
        "phi": events.Electron_phi, "charge": events.Electron_charge,
        "cb": events.Electron_cutBased,
    })
    ml = ak.zip({
        "pt": events.MLPhoton_pt, "eta": events.MLPhoton_eta,
        "phi": events.MLPhoton_phi, "mass": events.MLPhoton_mass,
        "dipscore": events.MLPhoton_diphotonScore,
    })

    n_pho_per_ev = ak.to_numpy(ak.num(ml))
    ml_eta_flat = ak.to_numpy(ak.flatten(ml.eta))
    ml_phi_flat = ak.to_numpy(ak.flatten(ml.phi))
    ml_mass_flat = ak.to_numpy(ak.flatten(ml.mass))
    ml_dip_flat = ak.to_numpy(ak.flatten(ml.dipscore))
    ml_offsets = np.cumsum(np.concatenate([[0], n_pho_per_ev]))

    # We need event index for each probe row. ak.broadcast_arrays gives
    # an event-indexed structure; redo by counting probes per event.
    # Easier alt: recompute event-by-event using the pre-flatten objects above.
    # We discard the previous flat outs and redo from scratch:
    out2 = {k: [] for k in
            ["probe_pt", "probe_eta", "probe_phi", "mll", "mlphoton_mass", "mlphoton_dr"]}

    # Per-event TnP loop using the ak structures
    n_ev = len(events)
    cb = ak.to_list(ele.cb)
    el_pt = ak.to_list(ele.pt)
    el_eta = ak.to_list(ele.eta)
    el_phi = ak.to_list(ele.phi)
    el_q = ak.to_list(ele.charge)

    for ev in range(n_ev):
        n_e = len(el_pt[ev])
        if n_e < 2:
            continue
        # find tag-probe pairs
        for i in range(n_e):
            for j in range(n_e):
                if i == j:
                    continue
                if cb[ev][i] < ELECTRON_TIGHT_BIT:
                    continue  # i must be tight tag
                if el_q[ev][i] * el_q[ev][j] != -1:
                    continue
                if el_pt[ev][j] < PROBE_ID_MIN_PT:
                    continue
                if abs(el_eta[ev][j]) >= 1.44:
                    continue
                # m_ll
                dphi = el_phi[ev][i] - el_phi[ev][j]
                mll = np.sqrt(
                    2.0 * el_pt[ev][i] * el_pt[ev][j]
                    * (np.cosh(el_eta[ev][i] - el_eta[ev][j]) - np.cos(dphi))
                )
                if not (ZMASS_LO < mll < ZMASS_HI):
                    continue

                # MLPhoton truth-match for probe j
                lo = ml_offsets[ev]; hi = ml_offsets[ev + 1]
                if hi == lo:
                    continue
                deta = ml_eta_flat[lo:hi] - el_eta[ev][j]
                dphi = (ml_phi_flat[lo:hi] - el_phi[ev][j] + np.pi) % (2 * np.pi) - np.pi
                dr = np.sqrt(deta**2 + dphi**2)
                k = int(np.argmin(dr))
                if dr[k] > MLPHOTON_MATCH_DR:
                    continue

                out2["probe_pt"].append(el_pt[ev][j])
                out2["probe_eta"].append(el_eta[ev][j])
                out2["probe_phi"].append(el_phi[ev][j])
                out2["mll"].append(mll)
                out2["mlphoton_mass"].append(ml_mass_flat[lo + k])
                out2["mlphoton_dr"].append(dr[k])

    return {k: np.asarray(v) for k, v in out2.items()}


def summarise_by_eta(probes: dict[str, np.ndarray]) -> None:
    eta = np.abs(probes["probe_eta"])
    mll = probes["mll"]
    m = probes["mlphoton_mass"]
    n = len(m)
    print(f"\nTotal selected probes: {n}")
    if n == 0:
        return
    for lo, hi, name in ETA_BINS:
        sel = (eta >= lo) & (eta < hi)
        n_b = sel.sum()
        if n_b == 0:
            print(f"  {name:8s} |eta| in ({lo},{hi}): 0 probes")
            continue
        mb = m[sel]
        print(f"  {name:8s} |eta| in ({lo},{hi}): N={n_b}, "
              f"mean(m_Γ)={mb.mean():.4f} GeV, "
              f"median={np.median(mb):.4f}, "
              f"std={mb.std():.4f}, "
              f"<mll>={mll[sel].mean():.2f}")


def plot_spectra(probes: dict[str, np.ndarray], outpath: str) -> None:
    if not HAVE_MPL or len(probes["mlphoton_mass"]) == 0:
        return
    fig, (ax_mll, ax_m) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: mll diagnostic — confirm Z peak in TnP selection
    mll_bins = np.linspace(60, 120, 61)
    ax_mll.hist(probes["mll"], bins=mll_bins, histtype="step", lw=1.8, color="steelblue",
                label=f"N={len(probes['mll'])}")
    ax_mll.axvline(91.2, ls="--", lw=1.2, color="gray", label="$m_Z = 91.2$ GeV")
    ax_mll.set_xlabel(r"$m_{\ell\ell}$ [GeV]")
    ax_mll.set_ylabel("tag-probe pairs / 1 GeV")
    ax_mll.set_title("TnP selection: dielectron mass")
    ax_mll.legend(fontsize=9)
    ax_mll.grid(alpha=0.3)

    # Right: probe MLPhoton regressed mass by |η| bin
    eta = np.abs(probes["probe_eta"])
    bins = np.linspace(0.0, 1.2, 121)
    for lo, hi, name in ETA_BINS:
        sel = (eta >= lo) & (eta < hi)
        if sel.sum() == 0:
            continue
        ax_m.hist(probes["mlphoton_mass"][sel], bins=bins, histtype="step", lw=1.6,
                  label=f"|η| {lo}-{hi} ({name}, N={sel.sum()})")
    ax_m.set_xlabel(r"$m_\Gamma^{\rm regressed}$ on Z→ee probe [GeV]")
    ax_m.set_ylabel("probes / 10 MeV")
    ax_m.set_title("MLPhoton regressed mass (Run3 2024)")
    ax_m.legend(fontsize=9)
    ax_m.grid(alpha=0.3)

    fig.suptitle("Z→ee tag-and-probe ONNX validation", fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=130)
    print(f"plot saved: {outpath}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="glob, dir, or single .root file with MLNanoAOD")
    ap.add_argument("--output", default="tnp_mlphoton.png")
    ap.add_argument("--max-events", type=int, default=None)
    args = ap.parse_args()

    if os.path.isdir(args.input):
        files = sorted(glob.glob(os.path.join(args.input, "*.root")))
    else:
        files = sorted(glob.glob(args.input))
    if not files:
        raise SystemExit(f"no input files found at {args.input!r}")
    print(f"input files: {len(files)} (first: {files[0]})")

    events = load_branches(files, max_events=args.max_events)
    print(f"loaded {len(events)} events")

    probes = tag_probe(events)
    summarise_by_eta(probes)
    plot_spectra(probes, args.output)


if __name__ == "__main__":
    main()
