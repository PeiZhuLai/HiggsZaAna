#!/usr/bin/env python3
"""
Derive sideband data/MC shape reweights for the Run 3 ZA BDT inputs.

The correction is shape-only by construction: after every variable update the
background sideband yield is rescaled back to the initial background sideband
yield.  The data/background normalization scale is recorded in the JSON but is
not folded into the per-event reweight.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np

pd = None
uproot = None


REWEIGHT_VARS = [
    # 1. angular variables: high importance and strongly correlated with topology
    "var_dR_Za",
    "var_dR_g1g2",
    "var_dR_g1Z",

    # 2. photon shower-shape variables: very important for BDT
    "pho1IetaIeta55",
    "pho2IetaIeta55",

    # 3. photon isolation variables: important and likely related to fake/photon mismodeling
    "pho1ECALIso",
    "pho2ECALIso",
    "ALP_calculatedPhotonIso",

    # 4. photon R9 variables
    "pho1R9",
    "pho2R9",

    # 5. kinematic variables (photon/Higgs pT normalized by H_m to decorrelate from m_llgammagamma)
    "var_PtaOverMh",
    "pho1Pt_oHm",
    "pho2Pt_oHm",
    "H_pt_oHm",

    # 6. derived BDT feature (computed inside load_frame, not stored in ntuples)
    "pho_pt_asym",

    # 7. invariant masses (sideband selection is on H_m: bins inside the signal
    # window get factor=1 since both data and bkg are empty there, so the
    # correction only acts on the sideband bins themselves)
    "Z_m",
    "H_m",

    # 8. mass-hypothesis variable: most important, but safest to correct last
    "param",
]

# Variables that are computed at load time from other branches (no physical branch).
DERIVED_VARS = {"pho_pt_asym"}

MASS_HYPOTHESES = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0]

BRANCH_ALIASES = {
    "pho1ECALIso": ["pho1ECALIso", "pho1PIso_noCorr", "ALP_lead_photon_ecalPFClusterIso"],
    "pho2ECALIso": ["pho2ECALIso", "pho2PIso_noCorr", "ALP_sublead_photon_ecalPFClusterIso"],
    "Z_m": ["Z_m", "Z_mass"],
}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    parser = argparse.ArgumentParser(
        description="Build iterative sideband data/MC reweight factors for ZA BDT inputs."
    )
    parser.add_argument(
        "--data-root",
        default="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/Data/run3.root",
        help="Input Data ROOT file.",
    )
    parser.add_argument(
        "--bkg-root",
        default="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/All_Bkg/run3.root",
        help="Input background ROOT file.",
    )
    parser.add_argument("--tree", default="inclusive", help="Tree name to read from each ROOT file.")
    parser.add_argument(
        "-o",
        "--output-json",
        default=str(project_dir / "reweights" / "sideband_run3_iterative.json"),
        help="Output JSON path.",
    )
    parser.add_argument(
        "--plot-dir",
        default=str(project_dir / "plots_MVA" / "sideband_reweight_run3"),
        help="Directory for validation plots.",
    )
    parser.add_argument("--iterations", type=int, default=5, help="Number of iterative reweighting rounds.")
    parser.add_argument("--bins", type=int, default=30, help="Number of bins per variable.")
    parser.add_argument(
        "--binning",
        choices=["quantile", "uniform"],
        default="quantile",
        help="Binning strategy. Quantile uses weighted background sideband quantiles.",
    )
    parser.add_argument("--clip-min", type=float, default=0.2, help="Minimum bin correction.")
    parser.add_argument("--clip-max", type=float, default=5.0, help="Maximum bin correction.")
    parser.add_argument(
        "--min-bkg-bin-yield",
        type=float,
        default=1e-8,
        help="Bins with background yield <= this value get factor 1.",
    )
    parser.add_argument("--seed", type=int, default=12345, help="Seed for deterministic param reconstruction.")
    parser.add_argument(
        "--weight-branch",
        default="factor",
        help="Preferred event weight branch. Falls back to 'weight' and then unit weights.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Write JSON/CSV only, without validation PDFs.",
    )
    return parser.parse_args()


def load_runtime_dependencies() -> None:
    global pd, uproot
    import pandas as _pd
    import uproot as _uproot

    pd = _pd
    uproot = _uproot


def json_float(value: float) -> float | None:
    value = float(value)
    if not math.isfinite(value):
        return None
    return value


def list_json_float(values: Sequence[float]) -> List[float | None]:
    return [json_float(v) for v in values]


def sideband_mask(frame: pd.DataFrame) -> np.ndarray:
    h_m = pd.to_numeric(frame["H_m"], errors="coerce").to_numpy(dtype=float)
    return ((h_m > 95.0) & (h_m < 115.0)) | ((h_m > 135.0) & (h_m < 180.0))


def pick_existing_branch(
    available: Sequence[str],
    logical_name: str,
) -> str | None:
    candidates = BRANCH_ALIASES.get(logical_name, [logical_name])
    available_set = set(available)
    for candidate in candidates:
        if candidate in available_set:
            return candidate
    return None


def required_logical_branches(vars_to_read: Sequence[str]) -> List[str]:
    branches = {"H_m", "ALP_m"}
    # Derived features need pho1Pt_oHm/pho2Pt_oHm and var_dR_g1g2 in the underlying frame.
    if any(v in DERIVED_VARS for v in vars_to_read):
        branches.update({"pho1Pt_oHm", "pho2Pt_oHm", "var_dR_g1g2"})
    for var in vars_to_read:
        if var == "param":
            continue
        if var in DERIVED_VARS:
            continue
        branches.add(var)
    return sorted(branches)


def resolve_branches(tree: uproot.behaviors.TTree.TTree, vars_to_read: Sequence[str], weight_branch: str) -> Tuple[List[str], Dict[str, str | None]]:
    available = list(tree.keys())
    branch_map: Dict[str, str | None] = {}
    missing = []

    for logical in required_logical_branches(vars_to_read):
        physical = pick_existing_branch(available, logical)
        branch_map[logical] = physical
        if physical is None:
            missing.append(logical)

    event_branch = pick_existing_branch(available, "event")
    if event_branch is not None:
        branch_map["event"] = event_branch

    physical_weight = pick_existing_branch(available, weight_branch)
    if physical_weight is None and weight_branch != "weight":
        physical_weight = pick_existing_branch(available, "weight")
    branch_map["__weight__"] = physical_weight

    if missing:
        raise KeyError(
            "Missing required branches: {}. Available branches include: {}".format(
                ", ".join(missing), ", ".join(available[:80])
            )
        )

    read_branches = sorted({branch for branch in branch_map.values() if branch is not None})
    return read_branches, branch_map


def deterministic_mass_from_event(events: np.ndarray, seed: int, masses: Sequence[float]) -> np.ndarray:
    masses_array = np.asarray(masses, dtype=float)
    selected = np.empty(events.shape[0], dtype=float)
    for idx, event in enumerate(events):
        payload = f"{seed}:{int(event)}".encode("utf-8", errors="ignore")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        selected[idx] = masses_array[int.from_bytes(digest, "little") % len(masses_array)]
    return selected


def deterministic_mass_from_index(n_events: int, seed: int, masses: Sequence[float]) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.choice(np.asarray(masses, dtype=float), size=n_events, replace=True)


def load_frame(path: str, tree_name: str, vars_to_read: Sequence[str], weight_branch: str, seed: int, label: str) -> Tuple[pd.DataFrame, Dict[str, str | None]]:
    with uproot.open(path) as root_file:
        tree = root_file[tree_name]
        read_branches, branch_map = resolve_branches(tree, vars_to_read, weight_branch)
        raw = tree.arrays(read_branches, library="pd")

    frame = pd.DataFrame(index=raw.index)
    for logical in required_logical_branches(vars_to_read):
        physical = branch_map[logical]
        if physical is None:
            raise KeyError(f"{label}: logical branch {logical} was not resolved")
        frame[logical] = pd.to_numeric(raw[physical], errors="coerce")

    for logical in vars_to_read:
        if logical == "param":
            continue
        if logical in DERIVED_VARS:
            continue
        physical = branch_map.get(logical)
        if physical is None:
            raise KeyError(f"{label}: logical branch {logical} was not resolved")
        frame[logical] = pd.to_numeric(raw[physical], errors="coerce")

    # Compute derived features from already-loaded base columns.
    # pho_pt_asym uses H_m-normalized pT; the H_m factor cancels so it equals the raw-pT asymmetry.
    if any(v in DERIVED_VARS for v in vars_to_read):
        pt1 = frame["pho1Pt_oHm"].to_numpy(dtype=float)
        pt2 = frame["pho2Pt_oHm"].to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            if "pho_pt_asym" in vars_to_read:
                frame["pho_pt_asym"] = (pt1 - pt2) / (pt1 + pt2 + 1e-6)

    physical_weight = branch_map["__weight__"]
    if physical_weight is None:
        frame["__weight__"] = 1.0
    else:
        frame["__weight__"] = pd.to_numeric(raw[physical_weight], errors="coerce").fillna(0.0)

    if "param" in vars_to_read:
        h_m = pd.to_numeric(frame["H_m"], errors="coerce").to_numpy(dtype=float)
        alp_m = pd.to_numeric(frame["ALP_m"], errors="coerce").to_numpy(dtype=float)
        if "event" in branch_map and branch_map["event"] is not None:
            events = pd.to_numeric(raw[branch_map["event"]], errors="coerce").fillna(-1).to_numpy(dtype=np.int64)
            mass_hyp = deterministic_mass_from_event(events, seed, MASS_HYPOTHESES)
            param_method = "event_hash"
        else:
            mass_hyp = deterministic_mass_from_index(frame.shape[0], seed, MASS_HYPOTHESES)
            param_method = "row_order_rng"
        with np.errstate(divide="ignore", invalid="ignore"):
            frame["param"] = (alp_m - mass_hyp) / h_m
        frame["__param_mass_hypothesis__"] = mass_hyp
        branch_map["__param_method__"] = param_method

    return frame, branch_map


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values) & np.isfinite(weights) & (np.abs(weights) > 0.0)
    values = values[finite]
    weights = np.abs(weights[finite])
    if values.size == 0:
        raise ValueError("Cannot build quantile bins from an empty finite sample.")

    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    return np.interp(quantiles, cdf, values)


def unique_monotonic_edges(edges: np.ndarray, fallback_values: np.ndarray, n_bins: int) -> np.ndarray:
    edges = np.asarray(edges, dtype=float)
    finite_edges = edges[np.isfinite(edges)]
    finite_edges = np.unique(finite_edges)
    if finite_edges.size >= 2:
        return finite_edges

    finite_values = fallback_values[np.isfinite(fallback_values)]
    if finite_values.size == 0:
        return np.asarray([0.0, 1.0], dtype=float)

    low = float(np.min(finite_values))
    high = float(np.max(finite_values))
    if high <= low:
        span = max(abs(low), 1.0) * 1e-3
        low -= span
        high += span
    return np.linspace(low, high, n_bins + 1)


def build_edges(
    data_sb: pd.DataFrame,
    bkg_sb: pd.DataFrame,
    bkg_weights: np.ndarray,
    vars_to_read: Sequence[str],
    n_bins: int,
    binning: str,
) -> Dict[str, np.ndarray]:
    edges = {}
    quantiles = np.linspace(0.0, 1.0, n_bins + 1)
    for var in vars_to_read:
        data_values = data_sb[var].to_numpy(dtype=float)
        bkg_values = bkg_sb[var].to_numpy(dtype=float)
        combined = np.concatenate([data_values[np.isfinite(data_values)], bkg_values[np.isfinite(bkg_values)]])
        if combined.size == 0:
            raise ValueError(f"{var}: no finite values in sideband.")

        if binning == "quantile":
            var_edges = weighted_quantile(bkg_values, bkg_weights, quantiles)
        else:
            var_edges = np.linspace(float(np.min(combined)), float(np.max(combined)), n_bins + 1)

        var_edges = unique_monotonic_edges(var_edges, combined, n_bins)
        var_edges[0] = min(var_edges[0], float(np.min(combined)))
        var_edges[-1] = max(var_edges[-1], float(np.max(combined)))
        edges[var] = var_edges
    return edges


def hist(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values) & np.isfinite(weights)
    if not np.any(finite):
        return np.zeros(len(edges) - 1, dtype=float)
    return np.histogram(values[finite], bins=edges, weights=weights[finite])[0].astype(float)


def normalized_hist(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    h = hist(values, weights, edges)
    total = float(np.sum(h))
    if total <= 0.0 or not math.isfinite(total):
        return np.zeros_like(h)
    return h / total


def lookup_factors(values: np.ndarray, edges: np.ndarray, factors: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(edges, values, side="right") - 1
    idx = np.clip(idx, 0, len(factors) - 1)
    out = factors[idx]
    return np.where(np.isfinite(values) & np.isfinite(out), out, 1.0)


def weighted_ks_from_hist(data_density: np.ndarray, bkg_density: np.ndarray) -> float:
    return float(np.max(np.abs(np.cumsum(data_density) - np.cumsum(bkg_density))))


def validation_metrics(
    data_sb: pd.DataFrame,
    bkg_sb: pd.DataFrame,
    data_weights: np.ndarray,
    bkg_weights: np.ndarray,
    edges_map: Mapping[str, np.ndarray],
    iteration: int,
    stage: str,
) -> List[Dict[str, float | int | str | None]]:
    rows: List[Dict[str, float | int | str | None]] = []
    data_yield = float(np.sum(data_weights[np.isfinite(data_weights)]))
    bkg_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))

    for var, edges in edges_map.items():
        data_density = normalized_hist(data_sb[var].to_numpy(dtype=float), data_weights, edges)
        bkg_density = normalized_hist(bkg_sb[var].to_numpy(dtype=float), bkg_weights, edges)
        diff = data_density - bkg_density
        denom = np.where(data_density + bkg_density > 0.0, data_density + bkg_density, np.nan)
        chi2 = np.nansum((diff * diff) / denom)
        ndf = int(np.count_nonzero(np.isfinite(denom)) - 1)
        rows.append(
            {
                "iteration": iteration,
                "stage": stage,
                "var": var,
                "data_yield": json_float(data_yield),
                "bkg_yield": json_float(bkg_yield),
                "bkg_over_data": json_float(bkg_yield / data_yield) if data_yield > 0 else None,
                "shape_l1": json_float(0.5 * np.sum(np.abs(diff))),
                "shape_ks": json_float(weighted_ks_from_hist(data_density, bkg_density)),
                "shape_chi2_ndf": json_float(chi2 / ndf) if ndf > 0 else None,
            }
        )
    return rows


def derive_bin_factors(
    data_values: np.ndarray,
    bkg_values: np.ndarray,
    data_weights: np.ndarray,
    bkg_weights: np.ndarray,
    edges: np.ndarray,
    clip_min: float,
    clip_max: float,
    min_bkg_bin_yield: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    data_shape = normalized_hist(data_values, data_weights, edges)
    bkg_shape = normalized_hist(bkg_values, bkg_weights, edges)

    raw = np.ones_like(bkg_shape)
    valid = np.isfinite(data_shape) & np.isfinite(bkg_shape) & (bkg_shape > min_bkg_bin_yield)
    raw[valid] = data_shape[valid] / bkg_shape[valid]
    raw[~np.isfinite(raw)] = 1.0

    clipped = np.clip(raw, clip_min, clip_max)
    info = {
        "n_bins": int(len(raw)),
        "n_invalid_denominator_bins": int(np.count_nonzero(~valid)),
        "n_clipped_low_bins": int(np.count_nonzero(raw < clip_min)),
        "n_clipped_high_bins": int(np.count_nonzero(raw > clip_max)),
    }
    return raw, clipped, info


def ensure_output_dirs(output_json: Path, plot_dir: Path) -> None:
    output_json.parent.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)


def write_validation_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def draw_validation_pdf(
    out_path: Path,
    data_sb: pd.DataFrame,
    bkg_sb: pd.DataFrame,
    data_weights: np.ndarray,
    bkg_weights: np.ndarray,
    edges_map: Mapping[str, np.ndarray],
    title: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_path) as pdf:
        for var, edges in edges_map.items():
            data_counts = hist(data_sb[var].to_numpy(dtype=float), data_weights, edges)
            bkg_counts = hist(bkg_sb[var].to_numpy(dtype=float), bkg_weights, edges)
            data_total = np.sum(data_counts)
            bkg_total = np.sum(bkg_counts)
            data_density = data_counts / data_total if data_total > 0 else np.zeros_like(data_counts)
            bkg_density = bkg_counts / bkg_total if bkg_total > 0 else np.zeros_like(bkg_counts)

            centers = 0.5 * (edges[:-1] + edges[1:])
            ratio = np.divide(data_density, bkg_density, out=np.full_like(data_density, np.nan), where=bkg_density > 0)

            fig, (ax, rax) = plt.subplots(
                2,
                1,
                figsize=(8, 7),
                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
                sharex=True,
            )
            ax.errorbar(centers, data_density, yerr=np.sqrt(np.maximum(data_counts, 0.0)) / data_total if data_total > 0 else None, fmt="o", color="black", label="Data sideband")
            ax.step(edges[:-1], bkg_density, where="post", color="tab:blue", linewidth=2, label="Bkg sideband")
            ax.set_ylabel("Normalized entries")
            ax.set_title(f"{title}: {var}")
            ax.legend(frameon=False)
            ax.grid(alpha=0.2)

            rax.axhline(1.0, color="gray", linewidth=1)
            rax.errorbar(centers, ratio, fmt="o", color="black")
            rax.set_ylim(0.0, 2.5)
            rax.set_ylabel("Data/Bkg")
            rax.set_xlabel(var)
            rax.grid(alpha=0.2)
            rax.set_xlim(edges[0], edges[-1])

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def make_payload(
    args: argparse.Namespace,
    branch_maps: Mapping[str, Mapping[str, str | None]],
    edges_map: Mapping[str, np.ndarray],
    steps_by_iteration: Sequence[Mapping[str, object]],
    validation_rows: Sequence[Mapping[str, object]],
    data_yield: float,
    bkg_initial_yield: float,
    bkg_final_yield: float,
) -> Dict[str, object]:
    return {
        "schema_version": 1,
        "created": datetime.now().isoformat(timespec="seconds"),
        "description": "Iterative sideband shape reweighting. Apply steps sequentially to background MC only.",
        "inputs": {
            "data_root": args.data_root,
            "bkg_root": args.bkg_root,
            "tree": args.tree,
            "data_branch_map": dict(branch_maps["data"]),
            "bkg_branch_map": dict(branch_maps["bkg"]),
        },
        "selection": {
            "sideband": "(95 < H_m < 115) or (135 < H_m < 180)",
            "signal_window_excluded": "115 < H_m < 135",
        },
        "settings": {
            "reweight_vars": list(REWEIGHT_VARS),
            "iterations": args.iterations,
            "bins": args.bins,
            "binning": args.binning,
            "clip_min": args.clip_min,
            "clip_max": args.clip_max,
            "min_bkg_bin_yield": args.min_bkg_bin_yield,
            "seed": args.seed,
            "weight_branch": args.weight_branch,
            "underflow_overflow": "use nearest edge bin",
        },
        "param": {
            "definition": "(ALP_m - mA_hypothesis) / H_m",
            "mass_hypotheses": MASS_HYPOTHESES,
            "method": branch_maps["bkg"].get("__param_method__", "unknown"),
            "note": "Use the same method and seed when applying param-dependent weights.",
        },
        "normalization": {
            "data_sideband_yield": json_float(data_yield),
            "bkg_initial_sideband_yield": json_float(bkg_initial_yield),
            "bkg_final_sideband_yield": json_float(bkg_final_yield),
            "data_over_bkg_initial": json_float(data_yield / bkg_initial_yield) if bkg_initial_yield > 0 else None,
            "preservation": "background sideband yield is restored after every variable update",
        },
        "binning": {
            var: {"edges": list_json_float(edges), "n_bins": int(len(edges) - 1)}
            for var, edges in edges_map.items()
        },
        "iterations": list(steps_by_iteration),
        "validation": list(validation_rows),
    }


def main() -> None:
    args = parse_args()
    load_runtime_dependencies()

    if args.iterations < 1:
        raise ValueError("--iterations must be >= 1")
    if args.bins < 1:
        raise ValueError("--bins must be >= 1")
    if args.clip_min <= 0.0 or args.clip_max <= 0.0 or args.clip_min > args.clip_max:
        raise ValueError("Require 0 < --clip-min <= --clip-max")

    output_json = Path(args.output_json)
    plot_dir = Path(args.plot_dir)
    ensure_output_dirs(output_json, plot_dir)

    print(f"[load] Data: {args.data_root}")
    data_frame, data_branch_map = load_frame(args.data_root, args.tree, REWEIGHT_VARS, args.weight_branch, args.seed, "Data")
    print(f"[load] Bkg : {args.bkg_root}")
    bkg_frame, bkg_branch_map = load_frame(args.bkg_root, args.tree, REWEIGHT_VARS, args.weight_branch, args.seed, "Bkg")

    data_sb = data_frame.loc[sideband_mask(data_frame)].copy()
    bkg_sb = bkg_frame.loc[sideband_mask(bkg_frame)].copy()
    if data_sb.empty:
        raise ValueError("Data sideband is empty.")
    if bkg_sb.empty:
        raise ValueError("Background sideband is empty.")

    data_weights = data_sb["__weight__"].to_numpy(dtype=float)
    bkg_base_weights = bkg_sb["__weight__"].to_numpy(dtype=float)
    bkg_weights = bkg_base_weights.copy()

    data_yield = float(np.sum(data_weights[np.isfinite(data_weights)]))
    bkg_initial_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))
    if data_yield <= 0.0:
        raise ValueError(f"Data sideband yield is non-positive: {data_yield}")
    if bkg_initial_yield <= 0.0:
        raise ValueError(f"Background sideband yield is non-positive: {bkg_initial_yield}")

    print(f"[sideband] Data entries={len(data_sb)} yield={data_yield:.6g}")
    print(f"[sideband] Bkg  entries={len(bkg_sb)} yield={bkg_initial_yield:.6g}")
    print(f"[sideband] Data/Bkg initial normalization scale={data_yield / bkg_initial_yield:.6g}")

    edges_map = build_edges(data_sb, bkg_sb, bkg_weights, REWEIGHT_VARS, args.bins, args.binning)

    validation_rows: List[Mapping[str, object]] = []
    validation_rows.extend(validation_metrics(data_sb, bkg_sb, data_weights, bkg_weights, edges_map, 0, "initial"))
    if not args.skip_plots:
        draw_validation_pdf(plot_dir / "validation_iter_00_initial.pdf", data_sb, bkg_sb, data_weights, bkg_weights, edges_map, "Initial")

    iterations_payload: List[Mapping[str, object]] = []
    for iteration in range(1, args.iterations + 1):
        print(f"[iter {iteration}] starting")
        steps = []

        for var in REWEIGHT_VARS:
            before_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))
            raw, clipped, info = derive_bin_factors(
                data_sb[var].to_numpy(dtype=float),
                bkg_sb[var].to_numpy(dtype=float),
                data_weights,
                bkg_weights,
                edges_map[var],
                args.clip_min,
                args.clip_max,
                args.min_bkg_bin_yield,
            )
            event_factors = lookup_factors(bkg_sb[var].to_numpy(dtype=float), edges_map[var], clipped)
            bkg_weights = bkg_weights * event_factors
            raw_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))
            norm_scale = bkg_initial_yield / raw_yield if raw_yield > 0.0 else 1.0
            bkg_weights = bkg_weights * norm_scale
            after_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))

            step = {
                "var": var,
                "edges": list_json_float(edges_map[var]),
                "raw_factors": list_json_float(raw),
                "clipped_factors": list_json_float(clipped),
                "normalization_scale": json_float(norm_scale),
                "bkg_yield_before": json_float(before_yield),
                "bkg_yield_after_raw": json_float(raw_yield),
                "bkg_yield_after_norm": json_float(after_yield),
                **info,
            }
            steps.append(step)

        round_validation = validation_metrics(data_sb, bkg_sb, data_weights, bkg_weights, edges_map, iteration, "post_iteration")
        validation_rows.extend(round_validation)
        mean_l1 = np.mean([row["shape_l1"] for row in round_validation if row["shape_l1"] is not None])
        max_ks = np.max([row["shape_ks"] for row in round_validation if row["shape_ks"] is not None])
        print(
            f"[iter {iteration}] done: bkg_yield={np.sum(bkg_weights):.6g}, "
            f"mean_shape_l1={mean_l1:.4f}, max_shape_ks={max_ks:.4f}"
        )

        iterations_payload.append({"iteration": iteration, "steps": steps})
        if not args.skip_plots:
            draw_validation_pdf(
                plot_dir / f"validation_iter_{iteration:02d}.pdf",
                data_sb,
                bkg_sb,
                data_weights,
                bkg_weights,
                edges_map,
                f"After iteration {iteration}",
            )

    bkg_final_yield = float(np.sum(bkg_weights[np.isfinite(bkg_weights)]))
    validation_csv = output_json.with_suffix(".validation.csv")
    write_validation_csv(validation_csv, validation_rows)

    payload = make_payload(
        args,
        {"data": data_branch_map, "bkg": bkg_branch_map},
        edges_map,
        iterations_payload,
        validation_rows,
        data_yield,
        bkg_initial_yield,
        bkg_final_yield,
    )
    with output_json.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)

    print(f"[output] JSON: {output_json}")
    print(f"[output] validation CSV: {validation_csv}")
    if not args.skip_plots:
        print(f"[output] validation PDFs: {plot_dir}")


if __name__ == "__main__":
    main()
