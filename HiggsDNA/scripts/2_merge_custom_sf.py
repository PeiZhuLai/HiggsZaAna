#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from copy import deepcopy
from pathlib import Path


DEFAULT_HIGGSDNA_DIR = Path(
    "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data"
)

ERA_DIRS = {
    "2022preEE": "2022preEE_UL",
    "2022postEE": "2022postEE_UL",
    "2023preBPix": "2023preBPix_UL",
    "2023postBPix": "2023postBPix_UL",
    "2024": "2024_UL",
}

ERA_TAGS = {
    "2022preEE": "2022",
    "2022postEE": "2022EE",
    "2023preBPix": "2023",
    "2023postBPix": "2023BPix",
    "2024": "2024",
}

CORRECTION_NAMES = ("sf_pass", "unc_pass", "sf_fail", "unc_fail")


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
    return {correction["name"]: correction for correction in correction_set["corrections"]}


def variable_kind(input_name: str) -> str:
    aliases = {
        "pt": "pt",
        "el_et": "pt",
        "eta": "eta",
        "el_sc_eta": "eta",
        "el_sc_abseta": "abseta",
    }
    if input_name not in aliases:
        raise ValueError(f"Unsupported correction input: {input_name}")
    return aliases[input_name]


def value_for_input(input_name: str, pt: float, eta: float) -> float:
    kind = variable_kind(input_name)
    if kind == "pt":
        return pt
    if kind == "eta":
        return eta
    if kind == "abseta":
        return abs(eta)
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


def bin_index(edges: list[float], value: float, flow: str | None) -> int:
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


def evaluate_multibinning(correction: dict, pt: float, eta: float) -> float:
    data = correction["data"]
    if data["nodetype"] != "multibinning":
        raise ValueError(f"Only multibinning corrections are supported: {correction['name']}")

    flat_index = 0
    for input_name, edges in zip(data["inputs"], data["edges"]):
        value = value_for_input(input_name, pt, eta)
        index = bin_index(edges, value, data.get("flow"))
        flat_index = flat_index * (len(edges) - 1) + index

    return data["content"][flat_index]


def hza_correction_template(source: dict, correction_name: str) -> dict:
    correction = deepcopy(source)
    correction["name"] = correction_name
    correction["inputs"] = [
        {"name": "pt", "type": "real", "description": "pt"},
        {"name": "eta", "type": "real", "description": "eta"},
    ]
    correction["output"] = {
        "name": "sf",
        "type": "real",
        "description": "data-MC SF" if correction_name.startswith("sf_") else "data-MC unc",
    }
    correction["data"]["inputs"] = ["pt", "eta"]
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

    # Priority: high-pT map covers Et >= 100, ECAL gap overrides the mid/low maps
    # below 100, low-pT covers remaining Et < 15, and nongap covers the rest.
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


def is_ecal_gap(eta: float) -> bool:
    return 1.444 <= abs(eta) < 1.566


def merge_electron_id_2024(base_dir: Path) -> Path:
    era = "2024"
    raw = raw_dir(base_dir, era)
    sources = {
        "gap": correction_by_name(load_json(raw / "hza_elid_gap_2024_sf.json")),
        "nongap": correction_by_name(load_json(raw / "hza_elid_nongap_2024_sf.json")),
        "highpt": correction_by_name(load_json(raw / "hza_elid_nongap_highpT_2024_sf.json")),
        "lowpt": correction_by_name(load_json(raw / "hza_elid_nongap_lowpT_2024_sf.json")),
    }

    sample_name = CORRECTION_NAMES[0]
    target_pt_edges = merged_edges(
        *(edges_for_kind(source[sample_name], "pt") for source in sources.values())
    )
    target_eta_edges = merged_edges(
        *(edges_for_kind(source[sample_name], "eta") for source in sources.values())
    )

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

    output = era_out_dir(base_dir, era) / "hza_elid_2024_scalefactors.json"
    write_json(make_correction_set(corrections), output)
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
    for era in ERA_DIRS:
        outputs.append(merge_photon(base_dir, era))
        if era == "2023postBPix":
            outputs.append(merge_photon(base_dir, era, hole=True))

    outputs.append(merge_electron_id_2024(base_dir))
    validate_with_correctionlib(outputs)


if __name__ == "__main__":
    main()
