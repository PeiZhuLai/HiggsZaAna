#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from copy import deepcopy
from itertools import product
from pathlib import Path
from typing import Optional


DEFAULT_HIGGSDNA_DIR = Path(
    "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data"
)

ERA_DIRS = {
    "2022preEE": "2022preEE_UL",
    "2022postEE": "2022postEE_UL",
    "2023preBPix": "2023preBPix_UL",
    "2023postBPix": "2023postBPix_UL",
    "2024": "2024_UL",
    "2025": "2025_UL",
}

ERA_TAGS = {
    "2022preEE": "2022",
    "2022postEE": "2022EE",
    "2023preBPix": "2023",
    "2023postBPix": "2023BPix",
    "2024": "2024",
    "2025": "2025",
}

PHOTON_ERAS = ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")
ELECTRON_ERAS = ("2024", "2025")
CORRECTION_NAMES = ("sf_pass", "unc_pass", "sf_fail", "unc_fail")
EFFICIENCY_NAMES = ("effdata", "systdata", "effmc", "systmc")
R9_SPLIT = 0.96
R9_MAX = 99999.0
ELECTRON_TRIGGER_SFS = (
    "dielleg12trigger",
    "dielleg23trigger",
    "sielleg30trigger",
)
ELECTRON_ISO_SFS = (
    ("elminiIso0p1", "eliso0p1"),
    ("elminiIso0p15", "eliso0p15"),
)
MUON_ISO_ERAS = ("2024", "2025")
MUON_ISO_EFFS = ("muiso0p1", "muiso0p15")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge custom HZa SF correctionlib JSONs into hza-style JSONs."
    )
    parser.add_argument(
        "--higgsdna-dir",
        default=os.environ.get("HiggsDNADir", str(DEFAULT_HIGGSDNA_DIR)),
        help="Directory containing <era>_UL SF data directories.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing input JSON: {path}")
    with path.open() as handle:
        return json.load(handle)


def write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"wrote {path}")


def correction_by_name(correction_set: dict) -> dict[str, dict]:
    return {correction["name"]: correction for correction in corrections_from_payload(correction_set)}


def corrections_from_payload(payload: dict) -> list[dict]:
    if "corrections" in payload:
        return payload["corrections"]
    if all(key in payload for key in ("name", "inputs", "output", "data")):
        return [payload]
    raise ValueError(
        "Unsupported correctionlib JSON payload: expected a CorrectionSet with "
        f"'corrections' or a single Correction object; found keys {sorted(payload)}"
    )


def variable_kind(input_name: str) -> str:
    aliases = {
        "pt": "pt",
        "el_et": "pt",
        "eta": "eta",
        "ph_eta": "eta",
        "ph_sc_eta": "eta",
        "el_sc_eta": "eta",
        "ph_abseta": "abseta",
        "ph_sc_abseta": "abseta",
        "el_sc_abseta": "abseta",
        "abseta": "abseta",
        "event_nPV": "npv",
        "nPV": "npv",
        "r9": "r9",
        "R9": "r9",
    }
    if input_name not in aliases:
        raise ValueError(f"Unsupported correction input: {input_name}")
    return aliases[input_name]


def value_for_input(
    input_name: str,
    pt: Optional[float] = None,
    eta: Optional[float] = None,
    r9: Optional[float] = None,
    npv: Optional[float] = None,
) -> float:
    kind = variable_kind(input_name)
    if kind == "pt":
        if pt is None:
            raise ValueError(f"No pt value available for input {input_name}")
        return pt
    if kind == "eta":
        if eta is None:
            raise ValueError(f"No eta value available for input {input_name}")
        return eta
    if kind == "abseta":
        if eta is None:
            raise ValueError(f"No eta value available for input {input_name}")
        return abs(eta)
    if kind == "r9":
        if r9 is None:
            raise ValueError(f"No r9 value available for input {input_name}")
        return r9
    if kind == "npv":
        if npv is None:
            raise ValueError(f"No event_nPV value available for input {input_name}")
        return npv
    raise ValueError(f"Unsupported correction input: {input_name}")


def edges_for_kind(correction: dict, kind: str) -> list[float]:
    data = correction["data"]
    for input_name, edges in zip(data["inputs"], data["edges"]):
        input_kind = variable_kind(input_name)
        if input_kind == kind:
            return edges
        if kind == "eta" and input_kind == "abseta":
            return sorted({-edge for edge in edges if edge != 0.0} | set(edges))
    raise ValueError(f"No {kind} edges found in correction {correction['name']}")


def merged_edges(*edge_lists: list[float]) -> list[float]:
    return sorted({edge for edges in edge_lists for edge in edges})


def midpoint(edges: list[float], index: int) -> float:
    return 0.5 * (edges[index] + edges[index + 1])


def bin_index(edges: list[float], value: float, flow: Optional[str]) -> int:
    if value < edges[0]:
        if flow == "clamp":
            return 0
        raise ValueError(f"value {value} below edge range {edges}")
    if value >= edges[-1]:
        if flow == "clamp":
            return len(edges) - 2
        raise ValueError(f"value {value} above edge range {edges}")

    for index in range(len(edges) - 1):
        if edges[index] <= value < edges[index + 1]:
            return index

    raise ValueError(f"value {value} not found in edge range {edges}")


def evaluate_multibinning(
    correction: dict,
    pt: Optional[float] = None,
    eta: Optional[float] = None,
    r9: Optional[float] = None,
    npv: Optional[float] = None,
) -> float:
    data = correction["data"]
    if data["nodetype"] != "multibinning":
        raise ValueError(f"Only multibinning corrections are supported: {correction['name']}")

    flat_index = 0
    for input_name, edges in zip(data["inputs"], data["edges"]):
        value = value_for_input(input_name, pt, eta, r9, npv)
        index = bin_index(edges, value, data.get("flow"))
        flat_index = flat_index * (len(edges) - 1) + index

    return data["content"][flat_index]


def point_value(point: dict[str, float], kind: str) -> Optional[float]:
    for input_name, value in point.items():
        if variable_kind(input_name) == kind:
            return value
    return None


def hza_correction_template(source: dict, correction_name: str) -> dict:
    output_descriptions = {
        "sf_pass": "data-MC SF",
        "sf_fail": "data-MC SF",
        "unc_pass": "data-MC unc",
        "unc_fail": "data-MC unc",
        "effdata": "data eff",
        "systdata": "data unc",
        "effmc": "MC eff",
        "systmc": "MC unc",
    }
    correction = deepcopy(source)
    correction["name"] = correction_name
    correction["inputs"] = [
        {"name": "pt", "type": "real", "description": "pt"},
        {"name": "eta", "type": "real", "description": "eta"},
    ]
    correction["output"] = {
        "name": "sf",
        "type": "real",
        "description": output_descriptions.get(correction_name, correction["output"]["description"]),
    }
    correction["data"]["inputs"] = ["pt", "eta"]
    return correction


def hza_input(name: str) -> dict:
    return {"name": name, "type": "real", "description": name}


def hza_r9_correction_template(
    source: dict,
    correction_name: str,
    input_names: list[str],
) -> dict:
    correction = deepcopy(source)
    correction["name"] = correction_name
    correction["inputs"] = [hza_input(input_name) for input_name in input_names]
    correction["output"] = {
        "name": "sf",
        "type": "real",
        "description": "data-MC SF" if correction_name.startswith("sf_") else "data-MC unc",
    }
    correction["data"]["inputs"] = input_names
    return correction


def make_merged_correction(
    correction_name: str,
    source_selector,
    target_pt_edges: list[float],
    target_eta_edges: list[float],
    template_source: dict,
) -> dict:
    correction = hza_correction_template(template_source, correction_name)
    content = []

    for pt_bin in range(len(target_pt_edges) - 1):
        pt = midpoint(target_pt_edges, pt_bin)
        for eta_bin in range(len(target_eta_edges) - 1):
            eta = midpoint(target_eta_edges, eta_bin)
            source = source_selector(correction_name, pt, eta)
            content.append(evaluate_multibinning(source, pt, eta))

    correction["data"]["edges"] = [target_pt_edges, target_eta_edges]
    correction["data"]["content"] = content
    correction["data"]["flow"] = "clamp"
    return correction


def make_r9_merged_correction(
    correction_name: str,
    source_selector,
    target_edges: list[list[float]],
    input_names: list[str],
    template_source: dict,
) -> dict:
    correction = hza_r9_correction_template(template_source, correction_name, input_names)
    content = []

    axis_bins = [range(len(edges) - 1) for edges in target_edges]
    for bin_indices in product(*axis_bins):
        point = {
            input_name: midpoint(edges, index)
            for input_name, edges, index in zip(input_names, target_edges, bin_indices)
        }
        source = source_selector(correction_name, point)
        content.append(
            evaluate_multibinning(
                source,
                pt=point_value(point, "pt"),
                eta=point_value(point, "eta"),
                r9=point_value(point, "r9"),
                npv=point_value(point, "npv"),
            )
        )

    correction["data"]["edges"] = target_edges
    correction["data"]["content"] = content
    correction["data"]["flow"] = "clamp"
    return correction


def make_correction_set(corrections: list[dict]) -> dict:
    return {
        "schema_version": 2,
        "corrections": corrections,
    }


def raw_dir(base_dir: Path, era: str) -> Path:
    return base_dir / ERA_DIRS[era] / "custom_SF_raw"


def era_out_dir(base_dir: Path, era: str) -> Path:
    return base_dir / ERA_DIRS[era]


def merge_photon(base_dir: Path, era: str, hole: bool = False) -> Path:
    suffix = f"{era}Hole" if hole else era
    raw = raw_dir(base_dir, era)
    lowpt = correction_by_name(load_json(raw / f"hza_resolve_phid_lowpt_{suffix}_sf.json"))
    highpt = correction_by_name(load_json(raw / f"hza_resolve_phid_{suffix}_sf.json"))

    sample_name = CORRECTION_NAMES[0]
    target_pt_edges = merged_edges(
        edges_for_kind(lowpt[sample_name], "pt"),
        edges_for_kind(highpt[sample_name], "pt"),
    )
    target_eta_edges = merged_edges(
        edges_for_kind(lowpt[sample_name], "eta"),
        edges_for_kind(highpt[sample_name], "eta"),
    )
    highpt_min = edges_for_kind(highpt[sample_name], "pt")[0]

    def selector(correction_name: str, pt: float, eta: float) -> dict:
        return lowpt[correction_name] if pt < highpt_min else highpt[correction_name]

    corrections = [
        make_merged_correction(
            correction_name,
            selector,
            target_pt_edges,
            target_eta_edges,
            highpt[correction_name],
        )
        for correction_name in CORRECTION_NAMES
    ]

    output_tag = f"{ERA_TAGS[era]}Hole" if hole else ERA_TAGS[era]
    output = era_out_dir(base_dir, era) / f"hza_phid_{output_tag}_scalefactors.json"
    write_json(make_correction_set(corrections), output)
    return output


def csev_output_axes(low_r9: dict, high_r9: dict) -> tuple[list[str], list[list[float]]]:
    sample_name = CORRECTION_NAMES[0]
    input_names = []
    edge_lists = []

    for input_name in low_r9[sample_name]["data"]["inputs"]:
        kind = variable_kind(input_name)
        if kind == "r9":
            continue

        if kind == "abseta":
            output_name = "eta"
        elif kind == "npv":
            output_name = "event_nPV"
        else:
            output_name = kind
        if output_name in input_names:
            continue

        input_names.append(output_name)
        edge_lists.append(
            merged_edges(
                edges_for_kind(low_r9[sample_name], kind),
                edges_for_kind(high_r9[sample_name], kind),
            )
        )

    input_names.append("r9")
    edge_lists.append([0.0, R9_SPLIT, R9_MAX])
    return input_names, edge_lists


def merge_photon_csev(base_dir: Path, era: str, hole: bool = False) -> Optional[Path]:
    suffix = f"{era}Hole" if hole else era
    raw = raw_dir(base_dir, era)
    low_path = raw / f"hza_resolve_phcsev_lr9_summary_{suffix}_sf.json"
    high_path = raw / f"hza_resolve_phcsev_hr9_summary_{suffix}_sf.json"

    if not low_path.exists() or not high_path.exists():
        print(f"skip CSEV merge for {suffix}: missing {low_path.name} or {high_path.name}")
        return None

    low_r9 = correction_by_name(load_json(low_path))
    high_r9 = correction_by_name(load_json(high_path))
    input_names, target_edges = csev_output_axes(low_r9, high_r9)

    def selector(correction_name: str, point: dict[str, float]) -> dict:
        return low_r9[correction_name] if point["r9"] < R9_SPLIT else high_r9[correction_name]

    corrections = [
        make_r9_merged_correction(
            correction_name,
            selector,
            target_edges,
            input_names,
            high_r9[correction_name],
        )
        for correction_name in CORRECTION_NAMES
    ]

    output_tag = f"{ERA_TAGS[era]}Hole" if hole else ERA_TAGS[era]
    output = era_out_dir(base_dir, era) / f"hza_phcsev_{output_tag}_scalefactors.json"
    write_json(make_correction_set(corrections), output)
    return output


def is_ecal_gap(eta: float) -> bool:
    return 1.444 <= abs(eta) < 1.566


def merge_electron_id(base_dir: Path, era: str) -> Path:
    raw = raw_dir(base_dir, era)
    sources = {
        "gap": correction_by_name(load_json(raw / f"hza_elid_gap_{era}_sf.json")),
        "nongap": correction_by_name(load_json(raw / f"hza_elid_nongap_{era}_sf.json")),
        "highpt": correction_by_name(load_json(raw / f"hza_elid_nongap_highpT_{era}_sf.json")),
        "lowpt": correction_by_name(load_json(raw / f"hza_elid_nongap_lowpT_{era}_sf.json")),
    }

    sample_name = CORRECTION_NAMES[0]
    target_pt_edges = merged_edges(
        *(edges_for_kind(source[sample_name], "pt") for source in sources.values())
    )
    target_eta_edges = merged_edges(
        *(edges_for_kind(source[sample_name], "eta") for source in sources.values())
    )

    # Priority: high-pT map covers Et >= 100, ECAL gap overrides the mid/low maps
    # below 100, low-pT covers remaining Et < 15, and nongap covers the rest.
    def selector(correction_name: str, pt: float, eta: float) -> dict:
        if pt >= 100.0:
            return sources["highpt"][correction_name]
        if is_ecal_gap(eta):
            return sources["gap"][correction_name]
        if pt < 15.0:
            return sources["lowpt"][correction_name]
        return sources["nongap"][correction_name]

    corrections = [
        make_merged_correction(
            correction_name,
            selector,
            target_pt_edges,
            target_eta_edges,
            sources["nongap"][correction_name],
        )
        for correction_name in CORRECTION_NAMES
    ]

    output = era_out_dir(base_dir, era) / f"hza_elid_{era}_scalefactors.json"
    write_json(make_correction_set(corrections), output)
    return output


def merge_electron_gap_nongap_efficiency(
    base_dir: Path,
    era: str,
    input_sf: str,
    output_sf: str,
    force_efficiency_output: bool = False,
) -> Path:
    raw = raw_dir(base_dir, era)
    sources = {
        "gap": correction_by_name(load_json(raw / f"hza_{input_sf}_gap_{era}_sf.json")),
        "nongap": correction_by_name(load_json(raw / f"hza_{input_sf}_nongap_{era}_sf.json")),
    }

    source_names = set(sources["gap"]) & set(sources["nongap"])
    if all(name in source_names for name in EFFICIENCY_NAMES):
        correction_names = EFFICIENCY_NAMES
        output_name = f"hza_{output_sf}_{era}_efficiencies.json"
    elif all(name in source_names for name in CORRECTION_NAMES):
        correction_names = CORRECTION_NAMES
        output_type = "efficiencies" if force_efficiency_output else "scalefactors"
        output_name = f"hza_{output_sf}_{era}_{output_type}.json"
    else:
        raise ValueError(
            f"Unsupported electron gap/nongap JSON schema for {input_sf}: "
            f"found {sorted(source_names)}"
        )

    sample_name = correction_names[0]
    target_pt_edges = merged_edges(
        edges_for_kind(sources["gap"][sample_name], "pt"),
        edges_for_kind(sources["nongap"][sample_name], "pt"),
    )
    target_eta_edges = merged_edges(
        edges_for_kind(sources["gap"][sample_name], "eta"),
        edges_for_kind(sources["nongap"][sample_name], "eta"),
    )

    def selector(correction_name: str, pt: float, eta: float) -> dict:
        if is_ecal_gap(eta):
            return sources["gap"][correction_name]
        return sources["nongap"][correction_name]

    corrections = [
        make_merged_correction(
            correction_name,
            selector,
            target_pt_edges,
            target_eta_edges,
            sources["nongap"][correction_name],
        )
        for correction_name in correction_names
    ]

    output = era_out_dir(base_dir, era) / output_name
    write_json(make_correction_set(corrections), output)
    return output


def merge_electron_trigger(base_dir: Path, era: str, trigger_sf: str) -> Path:
    return merge_electron_gap_nongap_efficiency(base_dir, era, trigger_sf, trigger_sf)


def merge_electron_iso(base_dir: Path, era: str, input_sf: str, output_sf: str) -> Path:
    return merge_electron_gap_nongap_efficiency(
        base_dir,
        era,
        input_sf,
        output_sf,
        force_efficiency_output=True,
    )


def muon_iso_efficiency(base_dir: Path, era: str, output_sf: str) -> Path:
    era_dir = era_out_dir(base_dir, era)
    source = era_dir / f"hzg_{output_sf}_{era}_efficiencies.json"
    output = era_dir / f"hza_{output_sf}_{era}_efficiencies.json"
    if not source.exists():
        raise FileNotFoundError(
            f"Missing collected muon miniIso efficiency JSON: {source}. "
            "Run 1_collect_custom_sf.sh after producing the hzg_muiso* efficiency JSONs."
        )
    payload = load_json(source)
    write_json(payload, output)
    return output


def validate_with_correctionlib(paths: list[Path]) -> None:
    try:
        from correctionlib import schemav2
    except ImportError:
        print("correctionlib is not available; skipped schema validation")
        return

    for path in paths:
        payload = load_json(path)
        if hasattr(schemav2.CorrectionSet, "model_validate"):
            schemav2.CorrectionSet.model_validate(payload)
        else:
            schemav2.CorrectionSet.parse_obj(payload)


def main() -> None:
    args = parse_args()
    base_dir = Path(args.higgsdna_dir)

    outputs = []
    for era in PHOTON_ERAS:
        outputs.append(merge_photon(base_dir, era))
        csev_output = merge_photon_csev(base_dir, era)
        if csev_output is not None:
            outputs.append(csev_output)
        if era == "2023postBPix":
            outputs.append(merge_photon(base_dir, era, hole=True))
            csev_hole_output = merge_photon_csev(base_dir, era, hole=True)
            if csev_hole_output is not None:
                outputs.append(csev_hole_output)

    for era in ELECTRON_ERAS:
        outputs.append(merge_electron_id(base_dir, era))
        for trigger_sf in ELECTRON_TRIGGER_SFS:
            outputs.append(merge_electron_trigger(base_dir, era, trigger_sf))
        for input_sf, output_sf in ELECTRON_ISO_SFS:
            outputs.append(merge_electron_iso(base_dir, era, input_sf, output_sf))

    for era in MUON_ISO_ERAS:
        for output_sf in MUON_ISO_EFFS:
            outputs.append(muon_iso_efficiency(base_dir, era, output_sf))
    validate_with_correctionlib(outputs)


if __name__ == "__main__":
    main()
