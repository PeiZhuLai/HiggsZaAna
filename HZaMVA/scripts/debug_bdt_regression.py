#!/usr/bin/env python3
"""Systematic regression diagnostics for the ZA BDT pipeline.

This script is designed to be safe to run in two modes:

1. Code-only mode:
   When the external CERN `/afs` or `/eos` inputs are not mounted, it still
   produces a useful report by diffing the training/application code,
   extracting model hyperparameters from pickle payloads, and comparing the
   provided reference result tables.

2. Full runtime mode:
   When ROOT inputs and optional Python dependencies are available in the
   execution environment, it also produces event-level summaries, threshold
   scans, feature diagnostics, smoothing comparisons, and score plots.

The script intentionally does not modify the training/application code.  It is
only a diagnostic layer over the existing pipeline.
"""

from __future__ import annotations

import argparse
import ast
import importlib
import json
import math
import pickletools
import re
import sys
import textwrap
import csv
import warnings
from array import array
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None


DEFAULT_MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
APPLICATION_SCAN_MASSES = list(range(1, 31))
KEY_MASSES = [1, 2, 3, 20, 25, 30]
DEFAULT_THRESHOLD_STEP = 0.001
DEFAULT_N_SCORE_BINS = 200
DEFAULT_SMOOTHING_BINS = 200

REFERENCE_OLD_RESULTS = {
    1: {"threshold": 0.955, "signal_eff": 0.492, "bkg_yield": 83.195},
    2: {"threshold": 0.980, "signal_eff": 0.669, "bkg_yield": 26.373},
    3: {"threshold": 0.985, "signal_eff": 0.763, "bkg_yield": 7.888},
    4: {"threshold": 0.980, "signal_eff": 0.842, "bkg_yield": 5.062},
    5: {"threshold": 0.985, "signal_eff": 0.849, "bkg_yield": 5.123},
    6: {"threshold": 0.990, "signal_eff": 0.817, "bkg_yield": 2.494},
    7: {"threshold": 0.985, "signal_eff": 0.860, "bkg_yield": 5.264},
    8: {"threshold": 0.990, "signal_eff": 0.799, "bkg_yield": 11.402},
    9: {"threshold": 0.990, "signal_eff": 0.784, "bkg_yield": 15.526},
    10: {"threshold": 0.990, "signal_eff": 0.771, "bkg_yield": 10.795},
    15: {"threshold": 0.990, "signal_eff": 0.701, "bkg_yield": 13.445},
    20: {"threshold": 0.990, "signal_eff": 0.626, "bkg_yield": 18.443},
    25: {"threshold": 0.985, "signal_eff": 0.644, "bkg_yield": 36.615},
    30: {"threshold": 0.980, "signal_eff": 0.674, "bkg_yield": 43.640},
}

REFERENCE_CURRENT_OBSERVED = {
    1: {"threshold": 0.920, "signal_eff": 0.336, "bkg_yield": 307.845},
    2: {"threshold": 0.965, "signal_eff": 0.511, "bkg_yield": 90.450},
    3: {"threshold": 0.975, "signal_eff": 0.591, "bkg_yield": 28.033},
    4: {"threshold": 0.975, "signal_eff": 0.674, "bkg_yield": 23.586},
    5: {"threshold": 0.985, "signal_eff": 0.650, "bkg_yield": 17.943},
    6: {"threshold": 0.985, "signal_eff": 0.711, "bkg_yield": 26.077},
    7: {"threshold": 0.990, "signal_eff": 0.645, "bkg_yield": 18.693},
    8: {"threshold": 0.990, "signal_eff": 0.648, "bkg_yield": 18.610},
    9: {"threshold": 0.990, "signal_eff": 0.643, "bkg_yield": 20.952},
    10: {"threshold": 0.990, "signal_eff": 0.642, "bkg_yield": 21.750},
    15: {"threshold": 0.990, "signal_eff": 0.572, "bkg_yield": 37.356},
    20: {"threshold": 0.985, "signal_eff": 0.578, "bkg_yield": 80.128},
    25: {"threshold": 0.985, "signal_eff": 0.489, "bkg_yield": 113.778},
    30: {"threshold": 0.980, "signal_eff": 0.503, "bkg_yield": 191.503},
}

FEATURE_BRANCH_ALIASES = {
    "H_m": ("H_m", "H_mass"),
    "ALP_m": ("ALP_m", "ALP_mass"),
    "param": ("param",),
    "event": ("event", "event_number", "EventNumber"),
    "pho1PIso_noCorr": ("pho1PIso_noCorr", "pho1ECALIso", "ALP_lead_photon_ecalPFClusterIso"),
    "pho2PIso_noCorr": ("pho2PIso_noCorr", "pho2ECALIso", "ALP_sublead_photon_ecalPFClusterIso"),
    "n_electrons": ("n_electrons", "nElectrons"),
    "n_muons": ("n_muons", "nMuons"),
    "Z_lead_lepton_id": ("Z_lead_lepton_id",),
    "Z_sublead_lepton_id": ("Z_sublead_lepton_id",),
}

DEFAULT_STATIC_WARNINGS = {
    "bkg_ratio_high": "background_yield_gt_2x_reference",
    "sig_eff_low": "signal_eff_lt_0p8_reference",
    "sample_yield_diff": "sample_weighted_yield_diff_gt_5pct",
    "feature_nan": "feature_has_nan",
    "feature_inf": "feature_has_inf",
    "feature_order": "feature_order_mismatch",
}


@dataclass
class Issue:
    severity: int
    title: str
    evidence: List[str]
    kind: str = "static"
    recommendation: str = ""


@dataclass
class ScriptSnapshot:
    name: str
    path: str
    variables: List[str] = field(default_factory=list)
    mass_variables: List[str] = field(default_factory=list)
    wt_variables: List[str] = field(default_factory=list)
    file_path: Optional[str] = None
    bkg_name: List[str] = field(default_factory=list)
    data_name: List[str] = field(default_factory=list)
    sig_name: List[str] = field(default_factory=list)
    years: List[str] = field(default_factory=list)
    mass_list: List[float] = field(default_factory=list)
    search_masses: List[float] = field(default_factory=list)
    tree_name: Optional[str] = None
    bkg_tree_name: Optional[str] = None
    sig_tree_name: Optional[str] = None
    bkg_selection: Optional[str] = None
    sig_selection: Optional[str] = None
    convert_selection: Optional[str] = None
    convert_selection_mode: Optional[str] = None
    split_test_size: Optional[float] = None
    split_random_state: Optional[int] = None
    split_stratify: bool = False
    objective_search_space: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    ks_pvalue_target: Optional[float] = None
    has_regularized_candidate_selection: bool = False
    background_param_masses: List[float] = field(default_factory=list)


@dataclass
class SampleSpec:
    study: str
    sample: str
    process: str
    year: str
    tree_name: str
    file_path: Path
    selection: Optional[str]
    weight_expr: str
    mass: Optional[int] = None
    split: str = "inclusive"


@dataclass
class RuntimeContext:
    output_dir: Path
    plots_dir: Path
    old_script: ScriptSnapshot
    new_script: ScriptSnapshot
    application_features: List[str]
    application_scan_masses: List[int]
    legacy_model_cfg: Dict[str, Any]
    current_model_cfg: Dict[str, Any]
    warnings: List[str] = field(default_factory=list)
    issues: List[Issue] = field(default_factory=list)
    logs: List[str] = field(default_factory=list)

    def log(self, message: str) -> None:
        print(message)
        self.logs.append(message)

    def warn(self, message: str) -> None:
        tagged = f"WARNING: {message}"
        self.log(tagged)
        self.warnings.append(message)


def safe_import(module_name: str):
    try:
        return importlib.import_module(module_name)
    except Exception:
        return None


def path_exists(path_like: Optional[str]) -> bool:
    if not path_like:
        return False
    return Path(path_like).exists()


def mass_to_label(mass: float | int) -> str:
    value = float(mass)
    if value.is_integer():
        return str(int(value))
    return str(value).replace(".", "p")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def json_dump(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))


def write_csv_rows(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Optional[Sequence[str]] = None) -> None:
    if fieldnames is None:
        keys = []
        seen = set()
        for row in rows:
            for key in row.keys():
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def ensure_csv_with_header(path: Path, fieldnames: Sequence[str]) -> None:
    write_csv_rows(path, [], fieldnames=fieldnames)


def rows_to_markdown(rows: Sequence[Mapping[str, Any]], columns: Optional[Sequence[str]] = None) -> str:
    if not rows:
        return "_No rows available._"
    if pd is not None:
        frame = pd.DataFrame(rows)
        if columns is not None:
            existing = [column for column in columns if column in frame.columns]
            frame = frame.loc[:, existing]
        return frame.to_markdown(index=False)
    if columns is None:
        seen = []
        used = set()
        for row in rows:
            for key in row.keys():
                if key not in used:
                    used.add(key)
                    seen.append(key)
        columns = seen
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        values = [str(row.get(column, "")) for column in columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def require_pandas() -> None:
    if pd is None:
        raise RuntimeError("pandas is not available in this environment")


def flatten_float_list(values: Sequence[Any]) -> List[float]:
    out: List[float] = []
    for value in values:
        try:
            out.append(float(value))
        except Exception:
            continue
    return out


def first_existing_column(frame: pd.DataFrame, logical_name: str) -> Optional[str]:
    for candidate in FEATURE_BRANCH_ALIASES.get(logical_name, (logical_name,)):
        if candidate in frame.columns:
            return candidate
    return None


def normalize_common_aliases(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    for logical_name, candidates in FEATURE_BRANCH_ALIASES.items():
        if logical_name in frame.columns:
            continue
        for candidate in candidates:
            if candidate in frame.columns:
                frame[logical_name] = frame[candidate]
                break
    return frame


def weighted_mean_std(values: np.ndarray, weights: np.ndarray) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights)
    values = values[mask]
    weights = weights[mask]
    if values.size == 0 or np.sum(np.abs(weights)) == 0:
        return float("nan"), float("nan")
    weights = np.abs(weights)
    mean = np.average(values, weights=weights)
    variance = np.average((values - mean) ** 2, weights=weights)
    return float(mean), float(math.sqrt(max(variance, 0.0)))


def weighted_ks_2samp(x1: np.ndarray, w1: np.ndarray, x2: np.ndarray, w2: np.ndarray) -> Tuple[float, float]:
    scipy_stats = safe_import("scipy.stats")
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    w1 = np.asarray(w1, dtype=float)
    w2 = np.asarray(w2, dtype=float)

    mask1 = np.isfinite(x1) & np.isfinite(w1)
    mask2 = np.isfinite(x2) & np.isfinite(w2)
    x1 = x1[mask1]
    x2 = x2[mask2]
    w1 = np.abs(w1[mask1])
    w2 = np.abs(w2[mask2])

    mask1 = w1 > 0
    mask2 = w2 > 0
    x1 = x1[mask1]
    x2 = x2[mask2]
    w1 = w1[mask1]
    w2 = w2[mask2]

    if x1.size == 0 or x2.size == 0:
        return float("nan"), float("nan")

    order1 = np.argsort(x1)
    order2 = np.argsort(x2)
    x1 = x1[order1]
    x2 = x2[order2]
    w1 = w1[order1]
    w2 = w2[order2]

    cdf1 = np.cumsum(w1) / np.sum(w1)
    cdf2 = np.cumsum(w2) / np.sum(w2)

    x_all = np.sort(np.concatenate([x1, x2]))
    idx1 = np.searchsorted(x1, x_all, side="right") - 1
    idx2 = np.searchsorted(x2, x_all, side="right") - 1
    cdf1_all = np.where(idx1 >= 0, cdf1[np.maximum(idx1, 0)], 0.0)
    cdf2_all = np.where(idx2 >= 0, cdf2[np.maximum(idx2, 0)], 0.0)
    stat = float(np.max(np.abs(cdf1_all - cdf2_all)))

    if scipy_stats is None:
        return stat, float("nan")

    neff1 = np.sum(w1) ** 2 / np.sum(w1 ** 2)
    neff2 = np.sum(w2) ** 2 / np.sum(w2 ** 2)
    en = math.sqrt(neff1 * neff2 / (neff1 + neff2))
    if en <= 0:
        return stat, float("nan")
    pvalue = float(scipy_stats.kstwobign.sf((en + 0.12 + 0.11 / en) * stat))
    return stat, pvalue


def asymptotic_significance(s: float, b: float, d_no_smooth: Optional[float] = None) -> float:
    if b <= 0 or s <= 0:
        return float("nan")
    value = (2.0 * (s + b) * math.log(1.0 + (s / b))) - 2.0 * s
    if value <= 0:
        return float("nan")
    if d_no_smooth is not None and d_no_smooth < 10:
        return float("nan")
    return math.sqrt(value)


def simple_s_over_sqrt_b(s: float, b: float) -> float:
    if b <= 0 or s <= 0:
        return float("nan")
    return s / math.sqrt(b)


def very_large_weight_fraction(weights: np.ndarray) -> Tuple[float, float]:
    weights = np.asarray(weights, dtype=float)
    finite = weights[np.isfinite(weights)]
    if finite.size == 0:
        return float("nan"), float("nan")
    abs_weight = np.abs(finite)
    positive = abs_weight[abs_weight > 0]
    if positive.size == 0:
        return 0.0, 0.0
    median = float(np.median(positive))
    std = float(np.std(positive))
    threshold = max(median * 10.0, np.mean(positive) + 5.0 * std)
    fraction = float(np.mean(abs_weight > threshold))
    return threshold, fraction


def weighted_separation(sig_values: np.ndarray, sig_weights: np.ndarray, bkg_values: np.ndarray, bkg_weights: np.ndarray, bins: int = 40) -> float:
    sig_values = np.asarray(sig_values, dtype=float)
    bkg_values = np.asarray(bkg_values, dtype=float)
    sig_weights = np.asarray(sig_weights, dtype=float)
    bkg_weights = np.asarray(bkg_weights, dtype=float)
    mask_sig = np.isfinite(sig_values) & np.isfinite(sig_weights)
    mask_bkg = np.isfinite(bkg_values) & np.isfinite(bkg_weights)
    sig_values = sig_values[mask_sig]
    bkg_values = bkg_values[mask_bkg]
    sig_weights = np.abs(sig_weights[mask_sig])
    bkg_weights = np.abs(bkg_weights[mask_bkg])
    if sig_values.size == 0 or bkg_values.size == 0:
        return float("nan")
    low = min(float(np.min(sig_values)), float(np.min(bkg_values)))
    high = max(float(np.max(sig_values)), float(np.max(bkg_values)))
    if not math.isfinite(low) or not math.isfinite(high) or low == high:
        return float("nan")
    hist_sig, edges = np.histogram(sig_values, bins=bins, range=(low, high), weights=sig_weights, density=True)
    hist_bkg, _ = np.histogram(bkg_values, bins=edges, weights=bkg_weights, density=True)
    denom = hist_sig + hist_bkg
    mask = denom > 0
    separation = 0.5 * np.sum(((hist_sig[mask] - hist_bkg[mask]) ** 2) / denom[mask]) * (edges[1] - edges[0])
    return float(separation)


def evaluate_simple_ast(node: ast.AST, env: Mapping[str, Any]) -> Any:
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.List):
        return [evaluate_simple_ast(element, env) for element in node.elts]
    if isinstance(node, ast.Tuple):
        return tuple(evaluate_simple_ast(element, env) for element in node.elts)
    if isinstance(node, ast.Dict):
        return {
            evaluate_simple_ast(key, env): evaluate_simple_ast(value, env)
            for key, value in zip(node.keys, node.values)
        }
    if isinstance(node, ast.Name):
        if node.id in env:
            return env[node.id]
        raise KeyError(node.id)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -evaluate_simple_ast(node.operand, env)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "format":
        base = evaluate_simple_ast(node.func.value, env)
        args = [evaluate_simple_ast(arg, env) for arg in node.args]
        kwargs = {kw.arg: evaluate_simple_ast(kw.value, env) for kw in node.keywords}
        return base.format(*args, **kwargs)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return evaluate_simple_ast(node.left, env) + evaluate_simple_ast(node.right, env)
    raise ValueError(f"Unsupported simple AST node: {ast.dump(node, include_attributes=False)}")


def extract_function_source(module_text: str, node: ast.FunctionDef) -> str:
    source = ast.get_source_segment(module_text, node)
    return source or ""


def extract_trial_suggestions(function_node: ast.FunctionDef) -> Dict[str, Dict[str, Any]]:
    search_space: Dict[str, Dict[str, Any]] = {}
    for node in ast.walk(function_node):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        call = node.value
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Attribute):
            continue
        if not call.func.attr.startswith("suggest_"):
            continue
        if not isinstance(call.func.value, ast.Name) or call.func.value.id != "trial":
            continue
        parameter_name = target.id
        entry: Dict[str, Any] = {"method": call.func.attr}
        args: List[Any] = []
        for argument in call.args:
            try:
                args.append(ast.literal_eval(argument))
            except Exception:
                args.append(ast.unparse(argument) if hasattr(ast, "unparse") else repr(argument))
        keywords: Dict[str, Any] = {}
        for keyword in call.keywords:
            try:
                keywords[keyword.arg] = ast.literal_eval(keyword.value)
            except Exception:
                keywords[keyword.arg] = ast.unparse(keyword.value) if hasattr(ast, "unparse") else repr(keyword.value)
        if args:
            entry["args"] = args
        if keywords:
            entry["keywords"] = keywords
        search_space[parameter_name] = entry
    return search_space


def extract_train_test_split_config(module_node: ast.Module, env: Mapping[str, Any]) -> Tuple[Optional[float], Optional[int], bool]:
    test_size: Optional[float] = None
    random_state: Optional[int] = None
    stratify = False
    for node in ast.walk(module_node):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not isinstance(func, ast.Name) or func.id != "train_test_split":
            continue
        for keyword in node.keywords:
            if keyword.arg == "test_size":
                try:
                    test_size = float(evaluate_simple_ast(keyword.value, env))
                except Exception:
                    pass
            elif keyword.arg == "random_state":
                try:
                    random_state = int(evaluate_simple_ast(keyword.value, env))
                except Exception:
                    pass
            elif keyword.arg == "stratify":
                stratify = True
        break
    return test_size, random_state, stratify


def extract_convert_selection_info(function_source: str) -> Tuple[Optional[str], Optional[str]]:
    if not function_source:
        return None, None
    uncommented_lines = []
    for line in function_source.splitlines():
        uncommented_lines.append(line.split("#", 1)[0])
    function_source = "\n".join(uncommented_lines)
    if re.search(r"\bselection\s*=\s*selection\b", function_source):
        return None, "forwarded_selection_argument"
    hardcoded = re.search(r"\bselection\s*=\s*['\"]([^'\"]+)['\"]", function_source)
    if hardcoded:
        return hardcoded.group(1), "hardcoded_tree2array_selection"
    return None, None


def parse_python_module(path: Path) -> ast.Module:
    text = path.read_text()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=SyntaxWarning)
        return ast.parse(text)


def build_script_snapshot(name: str, path: Path) -> ScriptSnapshot:
    text = path.read_text()
    module_node = parse_python_module(path)
    env: Dict[str, Any] = {}
    objective_function: Optional[ast.FunctionDef] = None
    convert_function: Optional[ast.FunctionDef] = None

    for node in module_node.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    try:
                        env[target.id] = evaluate_simple_ast(node.value, env)
                    except Exception:
                        continue
        elif isinstance(node, ast.FunctionDef):
            if node.name == "objective":
                objective_function = node
            elif node.name == "convert":
                convert_function = node

    split_test_size, split_random_state, split_stratify = extract_train_test_split_config(module_node, env)
    convert_selection, convert_selection_mode = extract_convert_selection_info(
        extract_function_source(text, convert_function) if convert_function else ""
    )

    snapshot = ScriptSnapshot(
        name=name,
        path=str(path),
        variables=list(env.get("variables", [])),
        mass_variables=list(env.get("mass_variables", [])),
        wt_variables=list(env.get("wt_variables", [])),
        file_path=env.get("file_path"),
        bkg_name=list(env.get("bkg_name", [])),
        data_name=list(env.get("data_name", [])),
        sig_name=list(env.get("sig_name", [])),
        years=[str(item) for item in env.get("years", [])],
        mass_list=flatten_float_list(env.get("mass_list", [])),
        search_masses=flatten_float_list(env.get("search_mA_list", env.get("search_masses", []))),
        tree_name=env.get("tree_name"),
        bkg_tree_name=env.get("bkg_tree_name"),
        sig_tree_name=env.get("sig_tree_name"),
        bkg_selection=env.get("bkg_data_selection"),
        sig_selection=env.get("sig_selection"),
        convert_selection=convert_selection,
        convert_selection_mode=convert_selection_mode,
        split_test_size=split_test_size,
        split_random_state=split_random_state,
        split_stratify=split_stratify,
        objective_search_space=extract_trial_suggestions(objective_function) if objective_function else {},
        ks_pvalue_target=float(env["KS_PVALUE_TARGET"]) if "KS_PVALUE_TARGET" in env else None,
        has_regularized_candidate_selection="regularized_xgb_candidates" in text and "fit_model_with_ks_target" in text,
        background_param_masses=flatten_float_list(env.get("mass_list", [])),
    )

    if not snapshot.search_masses and snapshot.mass_list:
        snapshot.search_masses = list(snapshot.mass_list)
    return snapshot


def extract_application_features(path: Path) -> List[str]:
    module_node = parse_python_module(path)
    env: Dict[str, Any] = {}
    for node in module_node.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    try:
                        env[target.id] = evaluate_simple_ast(node.value, env)
                    except Exception:
                        continue
    return list(env.get("MVA_FEATURE_COLUMNS", []))


def extract_application_scan_masses(path: Path) -> List[int]:
    module_node = parse_python_module(path)
    for node in module_node.body:
        if not isinstance(node, ast.FunctionDef) or node.name != "add_mva_scores":
            continue
        for default in node.args.defaults:
            try:
                value = evaluate_simple_ast(default, {})
            except Exception:
                continue
            if isinstance(value, range):
                return list(value)
            if isinstance(value, (list, tuple)):
                try:
                    return [int(item) for item in value]
                except Exception:
                    continue
    text = path.read_text()
    explicit_range = re.search(r"def\s+add_mva_scores\s*\([^)]*masses\s*=\s*range\((\d+)\s*,\s*(\d+)\)\)", text)
    if explicit_range:
        start = int(explicit_range.group(1))
        stop = int(explicit_range.group(2))
        return list(range(start, stop))
    return []


def extract_prepare_script_metadata(path: Path) -> Dict[str, Any]:
    text = path.read_text()
    info: Dict[str, Any] = {
        "all_bkg_components": [],
        "dygto2lg_years": [],
        "dyjets_years": [],
        "signal_mass_list": [],
        "signal_years": [],
    }

    if "DYJetsToLL_path" in text and "DYGto2LG_path" in text:
        info["all_bkg_components"] = ["DYJetsToLL", "DYGto2LG"]

    dyg_years = re.search(r"prepare_dygto2lg\(\).*?years=\(\s*([^)]+)\)", text, re.DOTALL)
    if dyg_years:
        info["dygto2lg_years"] = dyg_years.group(1).split()

    sig_meta = re.search(r"prepare_sig\(\).*?massList=\(\s*([^)]+)\).*?years=\(\s*([^)]+)\)", text, re.DOTALL)
    if sig_meta:
        info["signal_mass_list"] = sig_meta.group(1).split()
        info["signal_years"] = sig_meta.group(2).split()

    return info


def extract_largest_xgb_config_blob(path: Path) -> Optional[str]:
    biggest = b""
    for opcode, argument, _ in pickletools.genops(path.read_bytes()):
        if opcode.name == "BINBYTES" and isinstance(argument, (bytes, bytearray)) and len(argument) > len(biggest):
            biggest = bytes(argument)
    if not biggest:
        return None
    return biggest.decode("latin1", errors="ignore")


def extract_model_config_from_pickle(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "exists": False}
    config_blob = extract_largest_xgb_config_blob(path)
    if not config_blob:
        return {"path": str(path), "exists": True, "error": "no_xgboost_config_blob"}

    patterns = {
        "num_trees": r"num_treesSL.*?([0-9]+)",
        "max_depth": r"max_depthSL.*?([0-9]+)",
        "learning_rate": r"learning_rateSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "min_child_weight": r"min_child_weightSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "gamma": r"gammaSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "reg_alpha": r"reg_alphaSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "reg_lambda": r"reg_lambdaSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "subsample": r"subsampleSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "colsample_bytree": r"colsample_bytreeSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "scale_pos_weight": r"scale_pos_weightSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "best_score": r"best_scoreSL.*?([0-9]+(?:\.[0-9Ee+-]+)?)",
        "best_iteration": r"best_iterationSL.*?([0-9]+)",
        "num_feature": r"num_featureSL.*?([0-9]+)",
    }
    payload: Dict[str, Any] = {
        "path": str(path),
        "exists": True,
        "objective": "binary:logistic" if "binary:logistic" in config_blob else None,
        "eval_metric": "logloss" if "logloss" in config_blob else None,
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, config_blob, re.DOTALL)
        if not match:
            payload[key] = None
            continue
        value = match.group(1)
        try:
            payload[key] = int(value) if re.fullmatch(r"[0-9]+", value) else float(value)
        except Exception:
            payload[key] = value
    return payload


def root_selection_to_pandas(selection: str) -> str:
    translated = selection.strip()
    translated = translated.replace("&&", " & ").replace("||", " | ")
    translated = re.sub(r"(?<![=!<>])!(?!=)", "~", translated)
    return translated


def selection_tokens(selection: Optional[str]) -> List[str]:
    if not selection:
        return []
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", selection)
    excluded = {"and", "or", "not", "True", "False"}
    return sorted({token for token in tokens if token not in excluded and not token.isupper()})


def expression_tokens(expression: str) -> List[str]:
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression)
    excluded = {"and", "or", "not", "True", "False", "abs", "sqrt"}
    return sorted({token for token in tokens if token not in excluded})


def read_root_frame(path: Path, tree_name: str, branches: Sequence[str], entry_stop: Optional[int] = None) -> pd.DataFrame:
    uproot = safe_import("uproot")
    if uproot is None:
        raise RuntimeError("uproot is not available in this environment")
    with uproot.open(str(path)) as root_file:
        tree = root_file[tree_name]
        available = set(tree.keys())
        read_branches = [branch for branch in branches if branch in available]
        frame = tree.arrays(read_branches, library="pd", entry_stop=entry_stop)
    for branch in branches:
        if branch not in frame.columns:
            frame[branch] = np.nan
    return frame


def evaluate_expression(frame: pd.DataFrame, expression: str) -> pd.Series:
    if expression in frame.columns:
        return pd.to_numeric(frame[expression], errors="coerce")
    return frame.eval(expression, engine="python")


def apply_selection(frame: pd.DataFrame, selection: Optional[str], context: RuntimeContext, label: str) -> pd.DataFrame:
    if not selection:
        return frame
    translated = root_selection_to_pandas(selection)
    try:
        mask = frame.eval(translated, engine="python")
        if isinstance(mask, pd.Series):
            return frame.loc[mask.fillna(False)].copy()
    except Exception as exc:
        context.warn(f"{label}: failed to evaluate selection `{selection}` ({exc}); using unfiltered frame")
    return frame


def import_sideband_reweighter(repo_root: Path):
    scripts_dir = repo_root / "HZaMVA" / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    try:
        module = importlib.import_module("sideband_reweight")
    except Exception:
        return None
    try:
        return module.load_sideband_reweighter()
    except Exception:
        return None


def deterministic_mass_from_row_order(n_rows: int, seed: int, masses: Sequence[float]) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.choice(np.asarray(masses, dtype=float), size=n_rows, replace=True)


def assign_background_param(frame: pd.DataFrame, masses: Sequence[float], seed: int, sideband_reweighter=None) -> pd.DataFrame:
    frame = normalize_common_aliases(frame)
    if "param" in frame.columns and "mass" in frame.columns:
        return frame
    if sideband_reweighter is not None:
        try:
            sideband_reweighter.ensure_param(frame, output_col="param", mass_output_col="mass")
            return frame
        except Exception:
            pass
    if not masses:
        masses = DEFAULT_MASSES
    h_mass = pd.to_numeric(frame.get("H_m"), errors="coerce").to_numpy(dtype=float)
    alp_mass = pd.to_numeric(frame.get("ALP_m"), errors="coerce").to_numpy(dtype=float)
    mass_hyp = deterministic_mass_from_row_order(frame.shape[0], seed, masses)
    with np.errstate(divide="ignore", invalid="ignore"):
        param = (alp_mass - mass_hyp) / h_mass
    frame["mass"] = mass_hyp
    frame["param"] = param
    return frame


def apply_new_background_reweight(frame: pd.DataFrame, sideband_reweighter) -> pd.DataFrame:
    frame = frame.copy()
    if sideband_reweighter is None:
        return frame
    if "factor" not in frame.columns:
        return frame
    frame["factor_nominal"] = pd.to_numeric(frame["factor"], errors="coerce").fillna(0.0)
    try:
        sideband_reweighter.apply_to_dataframe(
            frame,
            base_weight_col="factor_nominal",
            output_weight_col="factor",
            output_reweight_col="sideband_rwgt",
        )
    except Exception:
        return frame
    return frame


def prepare_sample_frame(
    spec: SampleSpec,
    feature_order: Sequence[str],
    context: RuntimeContext,
    entry_stop: Optional[int] = None,
    sideband_reweighter=None,
) -> pd.DataFrame:
    required_branches = set(feature_order)
    required_branches.update({"H_m", "ALP_m", "event", "n_electrons", "n_muons", "Z_lead_lepton_id", "Z_sublead_lepton_id"})
    required_branches.update(expression_tokens(spec.weight_expr))
    required_branches.update(selection_tokens(spec.selection))
    required_branches.discard("param")
    required_branches.discard("mass")
    frame = read_root_frame(spec.file_path, spec.tree_name, sorted(required_branches), entry_stop=entry_stop)
    frame = normalize_common_aliases(frame)
    frame = apply_selection(frame, spec.selection, context, f"{spec.study}:{spec.sample}:{spec.tree_name}")
    if spec.process == "background" and spec.study == "new":
        frame = assign_background_param(
            frame,
            context.new_script.background_param_masses or context.new_script.mass_list,
            seed=12345,
            sideband_reweighter=sideband_reweighter,
        )
        frame = apply_new_background_reweight(frame, sideband_reweighter)
    elif spec.process == "background" and spec.study == "old":
        frame = assign_background_param(
            frame,
            context.old_script.background_param_masses or context.old_script.mass_list,
            seed=123,
            sideband_reweighter=None,
        )
    elif spec.process == "signal" and spec.mass is not None:
        h_mass = pd.to_numeric(frame.get("H_m"), errors="coerce").to_numpy(dtype=float)
        alp_mass = pd.to_numeric(frame.get("ALP_m"), errors="coerce").to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            frame["mass"] = float(spec.mass)
            frame["param"] = (alp_mass - float(spec.mass)) / h_mass
    elif "param" not in frame.columns:
        frame = assign_background_param(frame, DEFAULT_MASSES, seed=12345, sideband_reweighter=None)

    try:
        frame["event_weight"] = pd.to_numeric(evaluate_expression(frame, spec.weight_expr), errors="coerce").fillna(0.0)
    except Exception as exc:
        context.warn(f"{spec.study}:{spec.sample}:{spec.tree_name}: failed to evaluate weight `{spec.weight_expr}` ({exc})")
        frame["event_weight"] = 0.0
    frame["study"] = spec.study
    frame["sample"] = spec.sample
    frame["process"] = spec.process
    frame["year"] = spec.year
    frame["split"] = spec.split
    frame["tree_name"] = spec.tree_name
    if spec.mass is not None:
        frame["generated_mass"] = float(spec.mass)
    return frame


def build_old_specs(snapshot: ScriptSnapshot, masses: Sequence[int]) -> List[SampleSpec]:
    specs: List[SampleSpec] = []
    base = Path(snapshot.file_path or "")
    year = snapshot.years[0] if snapshot.years else "run2"
    tree_name = snapshot.tree_name or "passedEvents"
    weight_expr = snapshot.wt_variables[0] if snapshot.wt_variables else "factor"
    selection = snapshot.convert_selection or snapshot.sig_selection or "H_m>110 && H_m<180"

    for sample in snapshot.bkg_name:
        specs.append(
            SampleSpec(
                study="old",
                sample=sample,
                process="background",
                year=year,
                tree_name=tree_name,
                file_path=base / sample / f"{tree_name}_{year}.root",
                selection=selection,
                weight_expr=weight_expr,
                split="preselection",
            )
        )
    for sample in snapshot.data_name:
        specs.append(
            SampleSpec(
                study="old",
                sample=sample,
                process="data",
                year=year,
                tree_name=tree_name,
                file_path=base / sample / f"{tree_name}_{year}.root",
                selection=selection,
                weight_expr=weight_expr,
                split="preselection",
            )
        )
    signal_dir = snapshot.sig_name[0] if snapshot.sig_name else "sig"
    for mass in masses:
        specs.append(
            SampleSpec(
                study="old",
                sample=f"{signal_dir}_mA_{mass}",
                process="signal",
                year=year,
                tree_name=tree_name,
                file_path=base / signal_dir / f"ALP_M{mass}.root",
                selection=selection,
                weight_expr=weight_expr,
                mass=int(mass),
                split="preselection",
            )
        )
    return specs


def build_new_specs(snapshot: ScriptSnapshot, masses: Sequence[int]) -> List[SampleSpec]:
    specs: List[SampleSpec] = []
    base = Path(snapshot.file_path or "")
    year = snapshot.years[0] if snapshot.years else "run3"
    bkg_tree = snapshot.bkg_tree_name or "inclusive"
    sig_tree = snapshot.sig_tree_name or "train"
    weight_expr = snapshot.wt_variables[0] if snapshot.wt_variables else "factor"
    bkg_selection = snapshot.bkg_selection or "H_m>95 && H_m<180"
    sig_selection = snapshot.sig_selection or "H_m>115 && H_m<135"

    for sample in snapshot.bkg_name:
        specs.append(
            SampleSpec(
                study="new",
                sample=sample,
                process="background",
                year=year,
                tree_name=bkg_tree,
                file_path=base / sample / "run3.root",
                selection=bkg_selection,
                weight_expr=weight_expr,
                split="inclusive",
            )
        )
    for sample in snapshot.data_name:
        specs.append(
            SampleSpec(
                study="new",
                sample=sample,
                process="data",
                year=year,
                tree_name=bkg_tree,
                file_path=base / sample / "run3.root",
                selection=bkg_selection,
                weight_expr=weight_expr,
                split="inclusive",
            )
        )
    for mass in masses:
        specs.append(
            SampleSpec(
                study="new",
                sample=f"mA_M{mass}",
                process="signal",
                year=year,
                tree_name=sig_tree,
                file_path=base / f"mA_M{mass}" / "run3.root",
                selection=sig_selection,
                weight_expr=weight_expr,
                mass=int(mass),
                split="train_tree",
            )
        )
        specs.append(
            SampleSpec(
                study="new",
                sample=f"mA_M{mass}",
                process="signal",
                year=year,
                tree_name="test",
                file_path=base / f"mA_M{mass}" / "run3.root",
                selection=sig_selection,
                weight_expr=weight_expr,
                mass=int(mass),
                split="test_tree",
            )
        )
        specs.append(
            SampleSpec(
                study="new",
                sample=f"mA_M{mass}",
                process="signal",
                year=year,
                tree_name="inclusive",
                file_path=base / f"mA_M{mass}" / "run3.root",
                selection=sig_selection,
                weight_expr=weight_expr,
                mass=int(mass),
                split="inclusive",
            )
        )
    return specs


def summarize_frame(frame: pd.DataFrame) -> Dict[str, Any]:
    weights = pd.to_numeric(frame.get("event_weight"), errors="coerce").to_numpy(dtype=float)
    mean_weight, std_weight = weighted_mean_std(weights, np.ones_like(weights))
    very_large_threshold, very_large_fraction = very_large_weight_fraction(weights)
    unique_event_fraction = float("nan")
    event_col = first_existing_column(frame, "event")
    if event_col is not None and len(frame) > 0:
        try:
            unique_event_fraction = float(frame[event_col].nunique(dropna=True) / len(frame))
        except Exception:
            unique_event_fraction = float("nan")
    return {
        "n_entries": int(len(frame)),
        "sum_weight": float(np.nansum(weights)),
        "mean_weight": float(mean_weight),
        "std_weight": float(std_weight),
        "min_weight": float(np.nanmin(weights)) if len(weights) else float("nan"),
        "max_weight": float(np.nanmax(weights)) if len(weights) else float("nan"),
        "negative_weight_fraction": float(np.mean(weights < 0)) if len(weights) else float("nan"),
        "very_large_weight_threshold": float(very_large_threshold),
        "very_large_weight_fraction": float(very_large_fraction),
        "sum_abs_weight": float(np.nansum(np.abs(weights))),
        "unique_event_fraction": unique_event_fraction,
    }


def feature_summary_rows(
    frame: pd.DataFrame,
    features: Sequence[str],
    sample_type: str,
    study: str,
    mass: Optional[int] = None,
    comparison_frame: Optional[pd.DataFrame] = None,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    weights = pd.to_numeric(frame.get("event_weight"), errors="coerce").to_numpy(dtype=float)
    comparison_weights = None
    if comparison_frame is not None:
        comparison_weights = pd.to_numeric(comparison_frame.get("event_weight"), errors="coerce").to_numpy(dtype=float)

    for feature in features:
        if feature not in frame.columns:
            rows.append(
                {
                    "study": study,
                    "mass": mass,
                    "feature": feature,
                    "sample_type": sample_type,
                    "mean": np.nan,
                    "std": np.nan,
                    "min": np.nan,
                    "max": np.nan,
                    "nan_fraction": np.nan,
                    "inf_fraction": np.nan,
                    "zero_fraction": np.nan,
                    "separation": np.nan,
                }
            )
            continue

        values = pd.to_numeric(frame[feature], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        mean, std = weighted_mean_std(values[finite], weights[finite] if len(weights) == len(values) else np.ones(np.count_nonzero(finite)))

        row = {
            "study": study,
            "mass": mass,
            "feature": feature,
            "sample_type": sample_type,
            "mean": mean,
            "std": std,
            "min": float(np.nanmin(values)) if values.size else np.nan,
            "max": float(np.nanmax(values)) if values.size else np.nan,
            "nan_fraction": float(np.mean(np.isnan(values))) if values.size else np.nan,
            "inf_fraction": float(np.mean(np.isinf(values))) if values.size else np.nan,
            "zero_fraction": float(np.mean(values == 0.0)) if values.size else np.nan,
            "separation": np.nan,
        }
        if comparison_frame is not None and feature in comparison_frame.columns:
            comp_values = pd.to_numeric(comparison_frame[feature], errors="coerce").to_numpy(dtype=float)
            row["separation"] = weighted_separation(values, weights, comp_values, comparison_weights if comparison_weights is not None else np.ones_like(comp_values))
        rows.append(row)
    return rows


def build_sample_manifest(specs: Sequence[SampleSpec]) -> pd.DataFrame:
    rows = []
    for spec in specs:
        rows.append(
            {
                "study": spec.study,
                "sample": spec.sample,
                "year": spec.year,
                "process": spec.process,
                "split": spec.split,
                "tree_name": spec.tree_name,
                "file_path": str(spec.file_path),
                "path_exists": spec.file_path.exists(),
                "selection": spec.selection,
                "weight_expr": spec.weight_expr,
                "mass": spec.mass,
            }
        )
    if pd is not None:
        return pd.DataFrame(rows)
    return rows


def compare_feature_orders(training_features: Sequence[str], application_features: Sequence[str]) -> Tuple[bool, List[str]]:
    training = list(training_features)
    application = list(application_features)
    if training == application:
        return True, []
    messages = []
    if len(training) != len(application):
        messages.append(f"Different lengths: training={len(training)} application={len(application)}")
    max_len = max(len(training), len(application))
    for index in range(max_len):
        train_feature = training[index] if index < len(training) else "<missing>"
        app_feature = application[index] if index < len(application) else "<missing>"
        if train_feature != app_feature:
            messages.append(f"index {index}: training={train_feature}, application={app_feature}")
    return False, messages


def load_pickle_model(path: Path):
    if not path.exists():
        return None
    xgboost = safe_import("xgboost")
    if xgboost is None:
        return None
    import pickle

    with path.open("rb") as handle:
        return pickle.load(handle)


def compute_scores_with_model(model, frame: pd.DataFrame, features: Sequence[str]) -> np.ndarray:
    matrix = frame.loc[:, list(features)].to_numpy(dtype=float)
    return model.predict_proba(matrix)[:, 1]


def compute_scores(
    frame: pd.DataFrame,
    model,
    features: Sequence[str],
    mass: int,
    context: RuntimeContext,
) -> np.ndarray:
    if model is not None:
        missing = [feature for feature in features if feature not in frame.columns]
        if missing:
            raise KeyError(f"Missing model features: {missing}")
        return compute_scores_with_model(model, frame, features)

    score_branch = f"MVA_Score_mA_M{int(mass)}"
    if score_branch in frame.columns:
        return pd.to_numeric(frame[score_branch], errors="coerce").to_numpy(dtype=float)

    raise RuntimeError("No XGBoost model could be loaded and no precomputed score branch is available")


def compute_auc(y_true: np.ndarray, scores: np.ndarray, sample_weight: Optional[np.ndarray] = None) -> float:
    sklearn_metrics = safe_import("sklearn.metrics")
    if sklearn_metrics is None:
        return float("nan")
    fpr, tpr, _ = sklearn_metrics.roc_curve(y_true, scores, sample_weight=sample_weight)
    return float(sklearn_metrics.auc(fpr, tpr))


def find_threshold_bin(thresholds: np.ndarray, value: float) -> int:
    index = int(np.searchsorted(thresholds, value, side="left"))
    return min(max(index, 0), len(thresholds) - 1)


def build_smoothing_histograms(scores: np.ndarray, weights: np.ndarray, bins: int = DEFAULT_SMOOTHING_BINS, mva_low: float = 0.1):
    ROOT = safe_import("ROOT")
    if ROOT is None:
        return None

    hist = ROOT.TH1D("bkg_raw", "bkg_raw", bins, 0.0, 1.0)
    for score, weight in zip(scores, weights):
        if not np.isfinite(score) or not np.isfinite(weight):
            continue
        hist.Fill(float(score), float(weight))

    x_values = []
    y_values = []
    start_bin = hist.FindBin(mva_low)
    xaxis = hist.GetXaxis()
    for bin_idx in range(start_bin, hist.GetNbinsX() + 1):
        x_values.append(xaxis.GetBinCenter(bin_idx))
        y_values.append(hist.GetBinContent(bin_idx))
    if not x_values:
        return None

    graph = ROOT.TGraph(len(x_values), np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float))
    smoother = ROOT.TGraphSmooth()
    g_smooth = smoother.SmoothSuper(graph)

    smooth = ROOT.TH1D("bkg_smooth", "bkg_smooth", bins, 0.0, 1.0)
    smooth_up = ROOT.TH1D("bkg_smooth_up", "bkg_smooth_up", bins, 0.0, 1.0)
    smooth_dn = ROOT.TH1D("bkg_smooth_dn", "bkg_smooth_dn", bins, 0.0, 1.0)

    x = array("d", [0.0])
    y = array("d", [0.0])
    for bin_idx in range(1, hist.GetNbinsX() + 1):
        if bin_idx < start_bin:
            y[0] = hist.GetBinContent(bin_idx)
        else:
            g_smooth.GetPoint(bin_idx - start_bin, x, y)
        value = max(float(y[0]), 0.0)
        smooth.SetBinContent(bin_idx, value)
        smooth_up.SetBinContent(bin_idx, value + math.sqrt(value))
        smooth_dn.SetBinContent(bin_idx, max(value - math.sqrt(value), 0.0))

    return hist, smooth, smooth_up, smooth_dn


def smoothing_yield_above_threshold(smoothing_hists, threshold: float) -> Tuple[float, float]:
    if smoothing_hists is None:
        return float("nan"), float("nan")
    raw, smooth, _, _ = smoothing_hists
    raw_bin = raw.FindBin(float(threshold))
    smooth_bin = smooth.FindBin(float(threshold))
    raw_yield = float(raw.Integral(raw_bin, raw.GetNbinsX() + 1))
    smooth_yield = float(smooth.Integral(smooth_bin, smooth.GetNbinsX() + 1))
    return raw_yield, smooth_yield


def build_threshold_scan(
    mass: int,
    signal_frame: pd.DataFrame,
    background_frame: pd.DataFrame,
    data_frame: Optional[pd.DataFrame],
    signal_scores: np.ndarray,
    background_scores: np.ndarray,
    data_scores: Optional[np.ndarray],
    thresholds: np.ndarray,
    metric: str,
    smoothing_hists=None,
) -> pd.DataFrame:
    signal_weights = pd.to_numeric(signal_frame["event_weight"], errors="coerce").to_numpy(dtype=float)
    background_weights = pd.to_numeric(background_frame["event_weight"], errors="coerce").to_numpy(dtype=float)
    data_weights = (
        pd.to_numeric(data_frame["event_weight"], errors="coerce").to_numpy(dtype=float)
        if data_frame is not None and data_scores is not None
        else None
    )

    signal_denom_weight = float(np.nansum(signal_weights))
    signal_denom_count = int(len(signal_weights))
    rows: List[Dict[str, Any]] = []

    for threshold in thresholds:
        sig_mask = signal_scores >= threshold
        bkg_mask = background_scores >= threshold
        sig_num_weight = float(np.nansum(signal_weights[sig_mask]))
        sig_num_count = int(np.count_nonzero(sig_mask))
        raw_bkg_weight = float(np.nansum(background_weights[bkg_mask]))

        if smoothing_hists is not None:
            raw_smooth_reference, smooth_bkg_weight = smoothing_yield_above_threshold(smoothing_hists, float(threshold))
        else:
            raw_smooth_reference = float("nan")
            smooth_bkg_weight = raw_bkg_weight

        data_sideband_yield = float("nan")
        if data_frame is not None and data_scores is not None and data_weights is not None:
            data_mask = data_scores >= threshold
            data_sideband_yield = float(np.nansum(data_weights[data_mask]))

        if metric == "asimov":
            metric_value = asymptotic_significance(sig_num_weight, smooth_bkg_weight, data_sideband_yield)
        else:
            metric_value = simple_s_over_sqrt_b(sig_num_weight, smooth_bkg_weight)

        rows.append(
            {
                "mA": mass,
                "threshold": float(threshold),
                "signal_efficiency": sig_num_weight / signal_denom_weight if signal_denom_weight > 0 else np.nan,
                "signal_efficiency_unweighted": sig_num_count / signal_denom_count if signal_denom_count > 0 else np.nan,
                "signal_yield": sig_num_weight,
                "signal_yield_unweighted": sig_num_count,
                "background_yield_raw": raw_bkg_weight,
                "background_yield_smooth": smooth_bkg_weight,
                "background_yield_smooth_reference_raw": raw_smooth_reference,
                "data_sideband_yield": data_sideband_yield,
                "metric": metric_value,
            }
        )
    return pd.DataFrame(rows)


def channel_label(frame: pd.DataFrame) -> pd.Series:
    frame = normalize_common_aliases(frame)
    if "n_electrons" in frame.columns and "n_muons" in frame.columns:
        n_electrons = pd.to_numeric(frame["n_electrons"], errors="coerce").fillna(0)
        n_muons = pd.to_numeric(frame["n_muons"], errors="coerce").fillna(0)
        channel = pd.Series(np.where(n_electrons >= 2, "electron", np.where(n_muons >= 2, "muon", "other")), index=frame.index)
        return channel
    if "Z_lead_lepton_id" in frame.columns:
        lep_id = pd.to_numeric(frame["Z_lead_lepton_id"], errors="coerce").abs()
        return pd.Series(np.where(lep_id == 11, "electron", np.where(lep_id == 13, "muon", "other")), index=frame.index)
    return pd.Series(["combined"] * len(frame), index=frame.index)


def build_channel_summary(signal_frame: pd.DataFrame, background_frame: pd.DataFrame, threshold: float, mass: int) -> pd.DataFrame:
    sig_scores = signal_frame["score"].to_numpy(dtype=float)
    bkg_scores = background_frame["score"].to_numpy(dtype=float)
    sig_weights = signal_frame["event_weight"].to_numpy(dtype=float)
    bkg_weights = background_frame["event_weight"].to_numpy(dtype=float)
    sig_channel = channel_label(signal_frame)
    bkg_channel = channel_label(background_frame)

    rows = []
    for channel in ["electron", "muon", "combined"]:
        if channel == "combined":
            sig_mask = np.ones(len(signal_frame), dtype=bool)
            bkg_mask = np.ones(len(background_frame), dtype=bool)
        else:
            sig_mask = sig_channel.to_numpy() == channel
            bkg_mask = bkg_channel.to_numpy() == channel
        sig_denom = np.nansum(sig_weights[sig_mask])
        sig_num = np.nansum(sig_weights[sig_mask & (sig_scores >= threshold)])
        bkg_num = np.nansum(bkg_weights[bkg_mask & (bkg_scores >= threshold)])
        rows.append(
            {
                "mA": mass,
                "channel": channel,
                "signal_efficiency": sig_num / sig_denom if sig_denom > 0 else np.nan,
                "background_yield": bkg_num,
            }
        )
    return pd.DataFrame(rows)


def run_cutflow(frame: pd.DataFrame, selection: Optional[str], label: str, study: str) -> List[Dict[str, Any]]:
    if not selection:
        return []
    steps = [item.strip() for item in selection.split("&&")]
    current = frame.copy()
    rows: List[Dict[str, Any]] = []
    for step_index, step in enumerate(steps, start=1):
        current = apply_selection(current, step, dummy_context_for_cutflow(), f"{label}:{step_index}")
        rows.append(
            {
                "study": study,
                "label": label,
                "step_index": step_index,
                "cut": step,
                "n_entries": int(len(current)),
                "sum_weight": float(np.nansum(pd.to_numeric(current.get("event_weight"), errors="coerce"))),
            }
        )
    return rows


class _DummyCutflowContext:
    def warn(self, *_args, **_kwargs):
        return None


_DUMMY_CUTFLOW_CONTEXT = _DummyCutflowContext()


def dummy_context_for_cutflow():
    return _DUMMY_CUTFLOW_CONTEXT


def reference_table_dataframe(old_results: Mapping[int, Mapping[str, float]], new_results: Mapping[int, Mapping[str, float]]) -> pd.DataFrame:
    rows = []
    for mass in DEFAULT_MASSES:
        old = old_results.get(int(mass), {})
        new = new_results.get(int(mass), {})
        old_threshold = old.get("threshold", np.nan)
        new_threshold = new.get("threshold", np.nan)
        old_eff = old.get("signal_eff", np.nan)
        new_eff = new.get("signal_eff", np.nan)
        old_bkg = old.get("bkg_yield", np.nan)
        new_bkg = new.get("bkg_yield", np.nan)
        warnings = []
        if math.isfinite(old_bkg) and math.isfinite(new_bkg) and old_bkg > 0 and new_bkg / old_bkg > 2.0:
            warnings.append(DEFAULT_STATIC_WARNINGS["bkg_ratio_high"])
        if math.isfinite(old_eff) and math.isfinite(new_eff) and old_eff > 0 and new_eff / old_eff < 0.8:
            warnings.append(DEFAULT_STATIC_WARNINGS["sig_eff_low"])
        rows.append(
            {
                "mA": mass,
                "old_threshold": old_threshold,
                "new_threshold": new_threshold,
                "threshold_diff": new_threshold - old_threshold if np.isfinite(old_threshold) and np.isfinite(new_threshold) else np.nan,
                "old_signal_eff": old_eff,
                "new_signal_eff": new_eff,
                "signal_eff_ratio": new_eff / old_eff if np.isfinite(old_eff) and old_eff > 0 and np.isfinite(new_eff) else np.nan,
                "old_bkg_yield": old_bkg,
                "new_bkg_yield": new_bkg,
                "bkg_yield_ratio": new_bkg / old_bkg if np.isfinite(old_bkg) and old_bkg > 0 and np.isfinite(new_bkg) else np.nan,
                "warning_flag": ";".join(warnings),
            }
        )
    if pd is not None:
        return pd.DataFrame(rows)
    return rows


def add_static_issues(context: RuntimeContext, prepare_meta: Mapping[str, Any], summary_table) -> None:
    old = context.old_script
    new = context.new_script
    training_features = old.variables + ["param"]
    new_training_features = new.variables + ["param"]

    same_feature_order, feature_mismatch = compare_feature_orders(new_training_features, context.application_features)
    if not same_feature_order:
        context.issues.append(
            Issue(
                severity=0,
                title="Training/application feature order mismatch",
                kind="static",
                evidence=feature_mismatch[:10],
                recommendation="Raise immediately before model application if the feature lists differ.",
            )
        )

    if old.convert_selection and (
        old.convert_selection != new.convert_selection or old.convert_selection != new.bkg_selection
    ):
        context.issues.append(
            Issue(
                severity=0,
                title="Selection and control-region definition changed substantially",
                kind="static",
                evidence=[
                    f"Old tree2array selection is hardcoded to `{old.convert_selection}`.",
                    f"New convert() mode is `{new.convert_selection_mode}`.",
                    f"New convert() selection is `{new.convert_selection}`.",
                    f"New background selection is `{new.bkg_selection}`.",
                    f"New signal selection is `{new.sig_selection}`.",
                    "The old script applies passChaHadIso/passNeuHadIso/passdR_gl/passHOverE inside convert(); the new script now forwards externally supplied selections instead.",
                ],
                recommendation="Validate preselection yields before and after each selection component, especially in low-mass and high-mass regions.",
            )
        )

    old_weight = old.wt_variables[0] if old.wt_variables else "<missing>"
    new_weight = new.wt_variables[0] if new.wt_variables else "<missing>"
    context.issues.append(
        Issue(
            severity=0,
            title="Event-weight definition changed",
            kind="static",
            evidence=[
                f"Old training weight expression: `{old_weight}`.",
                f"New training weight expression: `{new_weight}`.",
                "Parque2Root_BDT.py currently sets `factor = weight_central`, while hzgml/skim_ntuples.py defines `weight_corr = weight_central × SFs × PU × trigger × ...`.",
                "run3_Za_BDT.py also converts negative training weights to abs(weight) and rebalances classes before fitting.",
            ],
            recommendation="Dump per-sample weight stats and compare factor, factor_nominal, and weight_corr-style compositions before retraining.",
        )
    )

    if prepare_meta.get("all_bkg_components"):
        context.issues.append(
            Issue(
                severity=1,
                title="Background composition changed",
                kind="static",
                evidence=[
                    f"Old training background samples: {old.bkg_name}.",
                    f"New top-level background sample: {new.bkg_name}, built from {prepare_meta['all_bkg_components']}.",
                    f"Signal preparation years in run3 inputs: {prepare_meta.get('signal_years', [])}.",
                ],
                recommendation="Compare preselection yields separately for DYJets and DYGto2LG-like components, especially in mA=1,2,3 and mA=20,25,30.",
            )
        )

    if new.has_regularized_candidate_selection:
        context.issues.append(
            Issue(
                severity=1,
                title="Current model selection is no longer pure Optuna best-trial",
                kind="static",
                evidence=[
                    f"KS_PVALUE_TARGET = {new.ks_pvalue_target}.",
                    "run3_Za_BDT.py defines `regularized_xgb_candidates()` and selects the first candidate passing a weighted-KS target.",
                    f"Legacy pickle hyperparameters: {context.legacy_model_cfg}.",
                    f"Current pickle hyperparameters: {context.current_model_cfg}.",
                ],
                recommendation="Check whether the KS-safe regularized model shifted the score tail enough to move thresholds down at low mA.",
            )
        )

    if path_exists(str(Path(new.path).parent / "sideband_reweight.py")):
        context.issues.append(
            Issue(
                severity=1,
                title="Sideband reweighting is injected directly into current background training weights",
                kind="static",
                evidence=[
                    "run3_Za_BDT.py loads sideband_run3_iterative.json and replaces `factor` with sideband-reweighted `factor` for background MC.",
                    "The JSON performs 5 iterative reweighting rounds over all 15 BDT variables, including `param`.",
                    "The reweight mass-hypothesis definition uses event-hash assignment over the 14 generated mass points, not the 1-30 scan.",
                ],
                recommendation="Compare score tails and feature distributions with sideband reweight disabled versus enabled.",
            )
        )

    training_mass_points = [int(mass) for mass in (new.background_param_masses or new.mass_list)]
    application_scan_masses = [int(mass) for mass in context.application_scan_masses]
    if application_scan_masses and training_mass_points and training_mass_points != application_scan_masses:
        context.issues.append(
            Issue(
                severity=1,
                title="Training mass hypotheses do not match the application scan",
                kind="static",
                evidence=[
                    f"Current background param assignment uses masses {training_mass_points}.",
                    f"Application computes scores for masses {application_scan_masses}.",
                    "run3_Za_BDT_v1.py previously assigned background masses from 1-30 for the parameterized training.",
                ],
                recommendation="Quantify how much mA interpolation/extrapolation degrades score ranking for 1/2/3 and 20/25/30.",
            )
        )

    if new.sig_tree_name == "train":
        context.issues.append(
            Issue(
                severity=1,
                title="There is no application split in the current run3 signal pipeline",
                kind="static",
                evidence=[
                    "Parque2Root_BDT.py writes only `train` and `test` trees when `--split` is used.",
                    "Background and data are stored only as `inclusive`.",
                    "run3_Za_BDT.py then performs an additional 50/50 train_test_split inside Python.",
                ],
                recommendation="Explicitly report tree-level and in-script splits before comparing any threshold optimization to older studies.",
            )
        )

    if pd is not None:
        degraded_rows = summary_table[
            (summary_table["bkg_yield_ratio"] > 2.0) | (summary_table["signal_eff_ratio"] < 0.8)
        ].sort_values(["bkg_yield_ratio", "signal_eff_ratio"], ascending=[False, True]).to_dict("records")
    else:
        degraded_rows = [
            row for row in summary_table
            if (row.get("bkg_yield_ratio") is not None and row.get("bkg_yield_ratio", float("nan")) > 2.0)
            or (row.get("signal_eff_ratio") is not None and row.get("signal_eff_ratio", float("nan")) < 0.8)
        ]
        degraded_rows = sorted(
            degraded_rows,
            key=lambda row: (
                row.get("bkg_yield_ratio", float("-inf")),
                -row.get("signal_eff_ratio", float("inf")),
            ),
            reverse=True,
        )
    if degraded_rows:
        evidence = []
        for row in degraded_rows:
            evidence.append(
                "mA={mA}: threshold {old_threshold:.3f}->{new_threshold:.3f}, "
                "eff ratio={signal_eff_ratio:.3f}, bkg ratio={bkg_yield_ratio:.3f}".format(**row)
            )
        context.issues.append(
            Issue(
                severity=0,
                title="Observed performance regression is concentrated in low and high mass regions",
                kind="reference_table",
                evidence=evidence[:10],
                recommendation="Prioritize diagnostics for mA=1,2,3,20,25,30 first.",
            )
        )


def write_markdown_report(context: RuntimeContext, runtime_sections: Dict[str, Any], summary_table) -> None:
    report_path = context.output_dir / "debug_report.md"
    issues = sorted(context.issues, key=lambda item: (item.severity, item.title))
    lines = [
        "# BDT Regression Debug Report",
        "",
        "## Executive Summary",
        "",
        "The strongest code-level regression candidates are ranked below.  When event-level inputs were available, the runtime findings are appended afterwards.",
        "",
    ]

    for issue in issues:
        lines.append(f"### [P{issue.severity}] {issue.title}")
        lines.append("")
        lines.append("Evidence:")
        for item in issue.evidence:
            lines.append(f"- {item}")
        if issue.recommendation:
            lines.append(f"Recommendation: {issue.recommendation}")
        lines.append("")

    lines.append("## Summary Comparison")
    lines.append("")
    summary_columns = [
        "mA",
        "old_threshold",
        "new_threshold",
        "threshold_diff",
        "old_signal_eff",
        "new_signal_eff",
        "signal_eff_ratio",
        "old_bkg_yield",
        "new_bkg_yield",
        "bkg_yield_ratio",
        "warning_flag",
    ]
    if pd is not None:
        lines.append(rows_to_markdown(summary_table.to_dict("records"), columns=summary_columns))
    else:
        lines.append(rows_to_markdown(summary_table, columns=summary_columns))
    lines.append("")

    if runtime_sections.get("runtime_notes"):
        lines.append("## Runtime Notes")
        lines.append("")
        for note in runtime_sections["runtime_notes"]:
            lines.append(f"- {note}")
        lines.append("")

    if context.warnings:
        lines.append("## Warnings")
        lines.append("")
        for warning in context.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    report_path.write_text("\n".join(lines))


def save_plot_score_distributions(
    path: Path,
    signal_train_scores: np.ndarray,
    signal_test_scores: np.ndarray,
    signal_train_weights: np.ndarray,
    signal_test_weights: np.ndarray,
    bkg_train_scores: np.ndarray,
    bkg_test_scores: np.ndarray,
    bkg_train_weights: np.ndarray,
    bkg_test_weights: np.ndarray,
    mass: int,
) -> None:
    matplotlib = safe_import("matplotlib")
    if matplotlib is None:
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    bins = np.linspace(0.0, 1.0, 51)
    axes[0].hist(signal_train_scores, bins=bins, weights=np.abs(signal_train_weights), density=True, histtype="step", linewidth=2.0, label="signal train")
    axes[0].hist(signal_test_scores, bins=bins, weights=np.abs(signal_test_weights), density=True, histtype="step", linewidth=2.0, linestyle="--", label="signal test")
    axes[0].set_title(f"Signal train/test, mA={mass}")
    axes[0].set_xlabel("BDT score")
    axes[0].set_ylabel("A.U.")
    axes[0].legend(frameon=False)

    axes[1].hist(bkg_train_scores, bins=bins, weights=np.abs(bkg_train_weights), density=True, histtype="step", linewidth=2.0, label="background train")
    axes[1].hist(bkg_test_scores, bins=bins, weights=np.abs(bkg_test_weights), density=True, histtype="step", linewidth=2.0, linestyle="--", label="background test")
    axes[1].set_title(f"Background train/test, mA={mass}")
    axes[1].set_xlabel("BDT score")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_plot_signal_background(
    path: Path,
    signal_scores: np.ndarray,
    background_scores: np.ndarray,
    signal_weights: np.ndarray,
    background_weights: np.ndarray,
    mass: int,
) -> None:
    matplotlib = safe_import("matplotlib")
    if matplotlib is None:
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(0.0, 1.0, 51)
    ax.hist(signal_scores, bins=bins, weights=np.abs(signal_weights), density=True, histtype="step", linewidth=2.0, label="signal")
    ax.hist(background_scores, bins=bins, weights=np.abs(background_weights), density=True, histtype="step", linewidth=2.0, label="background")
    ax.set_title(f"Signal vs background, mA={mass}")
    ax.set_xlabel("BDT score")
    ax.set_ylabel("A.U.")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_plot_threshold_scan(path: Path, scan_df: pd.DataFrame, mass: int) -> None:
    matplotlib = safe_import("matplotlib")
    if matplotlib is None:
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    axes[0].plot(scan_df["threshold"], scan_df["metric"], color="tab:blue", linewidth=2.0)
    axes[0].set_ylabel("Metric")
    axes[0].set_title(f"Threshold scan, mA={mass}")
    axes[1].plot(scan_df["threshold"], scan_df["signal_efficiency"], color="tab:green", linewidth=2.0, label="signal efficiency")
    axes[1].plot(scan_df["threshold"], scan_df["background_yield_smooth"], color="tab:red", linewidth=2.0, label="smoothed background")
    axes[1].set_xlabel("Threshold")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_plot_smoothing(path: Path, smoothing_hists, mass: int) -> None:
    matplotlib = safe_import("matplotlib")
    ROOT = safe_import("ROOT")
    if matplotlib is None or ROOT is None or smoothing_hists is None:
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    raw, smooth, _, _ = smoothing_hists
    centers = np.array([raw.GetBinCenter(i) for i in range(1, raw.GetNbinsX() + 1)], dtype=float)
    raw_values = np.array([raw.GetBinContent(i) for i in range(1, raw.GetNbinsX() + 1)], dtype=float)
    smooth_values = np.array([smooth.GetBinContent(i) for i in range(1, smooth.GetNbinsX() + 1)], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.step(centers, raw_values, where="mid", linewidth=1.5, label="raw MC")
    ax.step(centers, smooth_values, where="mid", linewidth=2.0, label="SmoothSuper")
    ax.set_yscale("log")
    ax.set_xlabel("BDT score")
    ax.set_ylabel("Weighted yield / bin")
    ax.set_title(f"Background smoothing, mA={mass}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_plot_feature_distributions(path: Path, signal_frame: pd.DataFrame, background_frame: pd.DataFrame, mass: int, features: Sequence[str]) -> None:
    matplotlib = safe_import("matplotlib")
    if matplotlib is None:
        return
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_features = [feature for feature in features if feature in signal_frame.columns and feature in background_frame.columns][:6]
    if not plot_features:
        return

    ncols = 2
    nrows = int(math.ceil(len(plot_features) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    sig_weights = np.abs(signal_frame["event_weight"].to_numpy(dtype=float))
    bkg_weights = np.abs(background_frame["event_weight"].to_numpy(dtype=float))
    for axis, feature in zip(axes, plot_features):
        sig_values = pd.to_numeric(signal_frame[feature], errors="coerce").to_numpy(dtype=float)
        bkg_values = pd.to_numeric(background_frame[feature], errors="coerce").to_numpy(dtype=float)
        low = np.nanmin([np.nanmin(sig_values), np.nanmin(bkg_values)])
        high = np.nanmax([np.nanmax(sig_values), np.nanmax(bkg_values)])
        if not np.isfinite(low) or not np.isfinite(high) or low == high:
            axis.axis("off")
            continue
        bins = np.linspace(low, high, 40)
        axis.hist(sig_values, bins=bins, weights=sig_weights, density=True, histtype="step", linewidth=2.0, label="signal")
        axis.hist(bkg_values, bins=bins, weights=bkg_weights, density=True, histtype="step", linewidth=2.0, label="background")
        axis.set_title(feature)
        axis.legend(frameon=False, fontsize=8)
    for axis in axes[len(plot_features):]:
        axis.axis("off")
    fig.suptitle(f"Key feature distributions, mA={mass}")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def runtime_sample_and_feature_summaries(
    context: RuntimeContext,
    specs: Sequence[SampleSpec],
    masses: Sequence[int],
    entry_stop: Optional[int],
    strict_feature_order: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_order = context.new_script.variables + ["param"]
    ok_order, mismatches = compare_feature_orders(feature_order, context.application_features)
    if strict_feature_order and not ok_order:
        raise RuntimeError("Feature list mismatch detected:\n" + "\n".join(mismatches))

    repo_root = Path(context.new_script.path).resolve().parents[2]
    sideband_reweighter = import_sideband_reweighter(repo_root)
    sample_rows: List[Dict[str, Any]] = []
    feature_rows: List[Dict[str, Any]] = []
    cutflow_rows: List[Dict[str, Any]] = []
    loaded_frames: Dict[Tuple[str, str, Optional[int], str], pd.DataFrame] = {}

    for spec in specs:
        manifest_key = (spec.study, spec.sample, spec.mass, spec.split)
        if not spec.file_path.exists():
            sample_rows.append(
                {
                    "study": spec.study,
                    "sample": spec.sample,
                    "year": spec.year,
                    "process": spec.process,
                    "split": spec.split,
                    "n_entries": np.nan,
                    "sum_weight": np.nan,
                    "mean_weight": np.nan,
                    "negative_weight_fraction": np.nan,
                    "std_weight": np.nan,
                    "min_weight": np.nan,
                    "max_weight": np.nan,
                    "very_large_weight_threshold": np.nan,
                    "very_large_weight_fraction": np.nan,
                    "sum_abs_weight": np.nan,
                    "unique_event_fraction": np.nan,
                    "file_path": str(spec.file_path),
                    "path_exists": False,
                }
            )
            continue

        try:
            frame = prepare_sample_frame(spec, feature_order, context, entry_stop=entry_stop, sideband_reweighter=sideband_reweighter)
        except Exception as exc:
            context.warn(f"Failed to load {spec.study}:{spec.sample}:{spec.tree_name} from {spec.file_path}: {exc}")
            sample_rows.append(
                {
                    "study": spec.study,
                    "sample": spec.sample,
                    "year": spec.year,
                    "process": spec.process,
                    "split": spec.split,
                    "n_entries": np.nan,
                    "sum_weight": np.nan,
                    "mean_weight": np.nan,
                    "negative_weight_fraction": np.nan,
                    "std_weight": np.nan,
                    "min_weight": np.nan,
                    "max_weight": np.nan,
                    "very_large_weight_threshold": np.nan,
                    "very_large_weight_fraction": np.nan,
                    "sum_abs_weight": np.nan,
                    "unique_event_fraction": np.nan,
                    "file_path": str(spec.file_path),
                    "path_exists": True,
                }
            )
            continue

        loaded_frames[manifest_key] = frame
        summary = summarize_frame(frame)
        sample_row = {
            "study": spec.study,
            "sample": spec.sample,
            "year": spec.year,
            "process": spec.process,
            "split": spec.split,
            **summary,
            "file_path": str(spec.file_path),
            "path_exists": True,
        }
        sample_rows.append(sample_row)

        if spec.process in {"signal", "background"}:
            comparison_key = None
            if spec.process == "signal" and spec.mass is not None and spec.study == "new":
                comparison_key = ("new", context.new_script.bkg_name[0] if context.new_script.bkg_name else "All_Bkg", None, "inclusive")
            comparison_frame = loaded_frames.get(comparison_key) if comparison_key is not None else None
            feature_rows.extend(
                feature_summary_rows(
                    frame,
                    feature_order,
                    sample_type=spec.process,
                    study=spec.study,
                    mass=spec.mass,
                    comparison_frame=comparison_frame,
                )
            )

        cutflow_rows.extend(run_cutflow(frame, spec.selection, f"{spec.sample}:{spec.split}", spec.study))

    split_rows: List[Dict[str, Any]] = []
    for spec in specs:
        split_rows.append(
            {
                "study": spec.study,
                "sample": spec.sample,
                "process": spec.process,
                "split": spec.split,
                "tree_name": spec.tree_name,
                "configured_selection": spec.selection,
                "test_size": context.old_script.split_test_size if spec.study == "old" else context.new_script.split_test_size,
                "random_state": context.old_script.split_random_state if spec.study == "old" else context.new_script.split_random_state,
                "stratify": context.old_script.split_stratify if spec.study == "old" else context.new_script.split_stratify,
            }
        )

    sample_df = pd.DataFrame(sample_rows)
    feature_df = pd.DataFrame(feature_rows)
    split_df = pd.DataFrame(split_rows)
    cutflow_df = pd.DataFrame(cutflow_rows)

    for _, row in feature_df[feature_df["nan_fraction"] > 0].iterrows():
        context.warn(f"{row['study']} {row['sample_type']} feature {row['feature']} has NaN fraction {row['nan_fraction']:.4g}")
    for _, row in feature_df[feature_df["inf_fraction"] > 0].iterrows():
        context.warn(f"{row['study']} {row['sample_type']} feature {row['feature']} has inf fraction {row['inf_fraction']:.4g}")

    return sample_df, feature_df, split_df, cutflow_df


def runtime_mass_score_diagnostics(
    context: RuntimeContext,
    masses: Sequence[int],
    entry_stop: Optional[int],
    metric: str,
    threshold_step: float,
    strict_feature_order: bool,
) -> Tuple[Dict[int, Dict[str, Any]], pd.DataFrame, pd.DataFrame]:
    feature_order = context.new_script.variables + ["param"]
    ok_order, mismatches = compare_feature_orders(feature_order, context.application_features)
    if strict_feature_order and not ok_order:
        raise RuntimeError("Feature list mismatch detected:\n" + "\n".join(mismatches))

    model = load_pickle_model(Path(context.current_model_cfg["path"])) if context.current_model_cfg.get("path") else None
    repo_root = Path(context.new_script.path).resolve().parents[2]
    sideband_reweighter = import_sideband_reweighter(repo_root)

    bkg_spec = SampleSpec(
        study="new",
        sample=context.new_script.bkg_name[0] if context.new_script.bkg_name else "All_Bkg",
        process="background",
        year=context.new_script.years[0] if context.new_script.years else "run3",
        tree_name=context.new_script.bkg_tree_name or "inclusive",
        file_path=Path(context.new_script.file_path or "") / (context.new_script.bkg_name[0] if context.new_script.bkg_name else "All_Bkg") / "run3.root",
        selection=context.new_script.bkg_selection or "H_m>95 && H_m<180",
        weight_expr=context.new_script.wt_variables[0] if context.new_script.wt_variables else "factor",
        split="inclusive",
    )
    data_spec = SampleSpec(
        study="new",
        sample=context.new_script.data_name[0] if context.new_script.data_name else "Data",
        process="data",
        year=context.new_script.years[0] if context.new_script.years else "run3",
        tree_name=context.new_script.bkg_tree_name or "inclusive",
        file_path=Path(context.new_script.file_path or "") / (context.new_script.data_name[0] if context.new_script.data_name else "Data") / "run3.root",
        selection=context.new_script.bkg_selection or "H_m>95 && H_m<180",
        weight_expr=context.new_script.wt_variables[0] if context.new_script.wt_variables else "factor",
        split="inclusive",
    )

    if not bkg_spec.file_path.exists() or not data_spec.file_path.exists():
        return {}, pd.DataFrame(), pd.DataFrame()

    background_frame = prepare_sample_frame(bkg_spec, feature_order, context, entry_stop=entry_stop, sideband_reweighter=sideband_reweighter)
    data_frame = prepare_sample_frame(data_spec, feature_order, context, entry_stop=entry_stop, sideband_reweighter=sideband_reweighter)
    background_frame = background_frame.copy()
    data_frame = data_frame.copy()

    threshold_tables: List[pd.DataFrame] = []
    channel_tables: List[pd.DataFrame] = []
    mass_payload: Dict[int, Dict[str, Any]] = {}

    SR_low, SR_high = 115.0, 135.0
    CR_lo_low, CR_lo_high = 95.0, 115.0
    CR_hi_low, CR_hi_high = 135.0, 180.0

    background_frame["is_sr"] = (background_frame["H_m"] > SR_low) & (background_frame["H_m"] < SR_high)
    data_frame["is_sb"] = ((data_frame["H_m"] > CR_lo_low) & (data_frame["H_m"] < CR_lo_high)) | ((data_frame["H_m"] > CR_hi_low) & (data_frame["H_m"] < CR_hi_high))

    if model is not None:
        background_frame["score"] = compute_scores(background_frame, model, feature_order, masses[0], context)
        data_frame["score"] = compute_scores(data_frame, model, feature_order, masses[0], context)
    elif f"MVA_Score_mA_M{masses[0]}" in background_frame.columns:
        background_frame["score"] = background_frame[f"MVA_Score_mA_M{masses[0]}"]
        data_frame["score"] = data_frame[f"MVA_Score_mA_M{masses[0]}"]
    else:
        context.warn("Current model could not be loaded and no precomputed BDT score branch is available; score diagnostics are skipped.")
        return {}, pd.DataFrame(), pd.DataFrame()

    for mass in masses:
        signal_train_spec = SampleSpec(
            study="new",
            sample=f"mA_M{mass}",
            process="signal",
            year=context.new_script.years[0] if context.new_script.years else "run3",
            tree_name=context.new_script.sig_tree_name or "train",
            file_path=Path(context.new_script.file_path or "") / f"mA_M{mass}" / "run3.root",
            selection=context.new_script.sig_selection or "H_m>115 && H_m<135",
            weight_expr=context.new_script.wt_variables[0] if context.new_script.wt_variables else "factor",
            mass=int(mass),
            split="train_tree",
        )
        signal_test_spec = SampleSpec(
            study="new",
            sample=f"mA_M{mass}",
            process="signal",
            year=context.new_script.years[0] if context.new_script.years else "run3",
            tree_name="test",
            file_path=Path(context.new_script.file_path or "") / f"mA_M{mass}" / "run3.root",
            selection=context.new_script.sig_selection or "H_m>115 && H_m<135",
            weight_expr=context.new_script.wt_variables[0] if context.new_script.wt_variables else "factor",
            mass=int(mass),
            split="test_tree",
        )
        if not signal_train_spec.file_path.exists():
            continue

        signal_train = prepare_sample_frame(signal_train_spec, feature_order, context, entry_stop=entry_stop, sideband_reweighter=sideband_reweighter)
        if signal_test_spec.file_path.exists():
            signal_test = prepare_sample_frame(signal_test_spec, feature_order, context, entry_stop=entry_stop, sideband_reweighter=sideband_reweighter)
        else:
            signal_test = signal_train.iloc[0:0].copy()

        if model is not None:
            signal_train["score"] = compute_scores(signal_train, model, feature_order, mass, context)
            signal_test["score"] = compute_scores(signal_test, model, feature_order, mass, context) if len(signal_test) else np.array([], dtype=float)
            background_frame[f"score_m{mass}"] = compute_scores(background_frame, model, feature_order, mass, context)
            data_frame[f"score_m{mass}"] = compute_scores(data_frame, model, feature_order, mass, context)
        else:
            score_branch = f"MVA_Score_mA_M{mass}"
            signal_train["score"] = pd.to_numeric(signal_train[score_branch], errors="coerce")
            signal_test["score"] = pd.to_numeric(signal_test[score_branch], errors="coerce")
            background_frame[f"score_m{mass}"] = pd.to_numeric(background_frame[score_branch], errors="coerce")
            data_frame[f"score_m{mass}"] = pd.to_numeric(data_frame[score_branch], errors="coerce")

        bkg_scores = background_frame[f"score_m{mass}"].to_numpy(dtype=float)
        data_scores = data_frame[f"score_m{mass}"].to_numpy(dtype=float)
        bkg_sr = background_frame.loc[background_frame["is_sr"]].copy()
        bkg_sr_scores = bkg_sr[f"score_m{mass}"].to_numpy(dtype=float)
        bkg_sr_weights = bkg_sr["event_weight"].to_numpy(dtype=float)
        data_sb = data_frame.loc[data_frame["is_sb"]].copy()
        data_sb_scores = data_sb[f"score_m{mass}"].to_numpy(dtype=float)
        signal_train["is_sr"] = (signal_train["H_m"] > SR_low) & (signal_train["H_m"] < SR_high)
        signal_sr = signal_train.loc[signal_train["is_sr"]].copy()
        signal_sr_scores = signal_sr["score"].to_numpy(dtype=float)

        smoothing_hists = build_smoothing_histograms(bkg_sr_scores, bkg_sr_weights, bins=DEFAULT_SMOOTHING_BINS)
        thresholds = np.arange(0.0, 1.0 + 0.5 * threshold_step, threshold_step)
        scan_df = build_threshold_scan(
            mass,
            signal_sr,
            bkg_sr.assign(score=bkg_sr_scores),
            data_sb.assign(score=data_sb_scores),
            signal_sr_scores,
            bkg_sr_scores,
            data_sb_scores,
            thresholds,
            metric=metric,
            smoothing_hists=smoothing_hists,
        )
        threshold_tables.append(scan_df)
        scan_path = context.output_dir / f"threshold_scan_mass{mass}.csv"
        scan_df.to_csv(scan_path, index=False)

        best_row = scan_df.loc[scan_df["metric"].idxmax()] if not scan_df["metric"].dropna().empty else None
        best_threshold = float(best_row["threshold"]) if best_row is not None else float("nan")
        mass_payload[mass] = {
            "best_threshold": best_threshold,
            "best_metric": float(best_row["metric"]) if best_row is not None else float("nan"),
            "signal_efficiency": float(best_row["signal_efficiency"]) if best_row is not None else float("nan"),
            "background_yield": float(best_row["background_yield_smooth"]) if best_row is not None else float("nan"),
        }

        if len(signal_test) and len(background_frame):
            bkg_train = background_frame.iloc[::2].copy()
            bkg_test = background_frame.iloc[1::2].copy()
            bkg_train_scores = bkg_train[f"score_m{mass}"].to_numpy(dtype=float)
            bkg_test_scores = bkg_test[f"score_m{mass}"].to_numpy(dtype=float)
            p_sig = weighted_ks_2samp(
                signal_train["score"].to_numpy(dtype=float),
                signal_train["event_weight"].to_numpy(dtype=float),
                signal_test["score"].to_numpy(dtype=float),
                signal_test["event_weight"].to_numpy(dtype=float),
            )[1]
            p_bkg = weighted_ks_2samp(
                bkg_train_scores,
                bkg_train["event_weight"].to_numpy(dtype=float),
                bkg_test_scores,
                bkg_test["event_weight"].to_numpy(dtype=float),
            )[1]
            test_y = np.concatenate([np.ones(len(signal_test)), np.zeros(len(bkg_test))])
            test_scores = np.concatenate([signal_test["score"].to_numpy(dtype=float), bkg_test_scores])
            test_weights = np.concatenate([np.abs(signal_test["event_weight"].to_numpy(dtype=float)), np.abs(bkg_test["event_weight"].to_numpy(dtype=float))])
            auc_value = compute_auc(test_y, test_scores, sample_weight=test_weights)
        else:
            bkg_train = background_frame.iloc[::2].copy()
            bkg_test = background_frame.iloc[1::2].copy()
            bkg_train_scores = bkg_train[f"score_m{mass}"].to_numpy(dtype=float)
            bkg_test_scores = bkg_test[f"score_m{mass}"].to_numpy(dtype=float)
            p_sig = float("nan")
            p_bkg = float("nan")
            auc_value = float("nan")

        mass_payload[mass].update(
            {
                "auc_test": auc_value,
                "ks_p_signal": p_sig,
                "ks_p_background": p_bkg,
                "signal_score_mean": float(np.nanmean(signal_sr_scores)) if len(signal_sr_scores) else float("nan"),
                "signal_score_std": float(np.nanstd(signal_sr_scores)) if len(signal_sr_scores) else float("nan"),
                "background_score_mean": float(np.nanmean(bkg_sr_scores)) if len(bkg_sr_scores) else float("nan"),
                "background_score_std": float(np.nanstd(bkg_sr_scores)) if len(bkg_sr_scores) else float("nan"),
                "signal_denominator_yield": float(np.nansum(signal_sr["event_weight"])),
                "signal_denominator_entries": int(len(signal_sr)),
            }
        )

        if best_row is not None:
            channel_tables.append(build_channel_summary(signal_sr, bkg_sr.assign(score=bkg_sr_scores), best_threshold, mass))

        save_plot_score_distributions(
            context.plots_dir / f"bdt_score_train_test_mass{mass}.pdf",
            signal_train["score"].to_numpy(dtype=float),
            signal_test["score"].to_numpy(dtype=float),
            signal_train["event_weight"].to_numpy(dtype=float),
            signal_test["event_weight"].to_numpy(dtype=float),
            bkg_train_scores,
            bkg_test_scores,
            bkg_train["event_weight"].to_numpy(dtype=float),
            bkg_test["event_weight"].to_numpy(dtype=float),
            mass,
        )
        save_plot_signal_background(
            context.plots_dir / f"bdt_score_signal_background_mass{mass}.pdf",
            signal_sr_scores,
            bkg_sr_scores,
            signal_sr["event_weight"].to_numpy(dtype=float),
            bkg_sr_weights,
            mass,
        )
        save_plot_threshold_scan(context.plots_dir / f"threshold_scan_mass{mass}.pdf", scan_df, mass)
        save_plot_smoothing(context.plots_dir / f"background_smoothing_before_after_mass{mass}.pdf", smoothing_hists, mass)
        save_plot_feature_distributions(
            context.plots_dir / f"feature_distribution_key_variables_mass{mass}.pdf",
            signal_sr,
            bkg_sr,
            mass,
            [
                "param",
                "ALP_m",
                "H_m",
                "pho1Pt",
                "pho2Pt",
                "pho1PIso_noCorr",
                "pho2PIso_noCorr",
                "ALP_calculatedPhotonIso",
                "var_dR_Za",
                "var_dR_g1g2",
                "var_dR_g1Z",
            ],
        )

    threshold_summary = pd.DataFrame(
        [
            {
                "mA": mass,
                **payload,
            }
            for mass, payload in sorted(mass_payload.items())
        ]
    )
    channel_df = pd.concat(channel_tables, ignore_index=True) if channel_tables else pd.DataFrame()
    return mass_payload, threshold_summary, channel_df


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Debug ZA BDT performance regressions without touching the core training flow.")
    parser.add_argument("--old-script", default=str(repo_root / "HZaMVA" / "scripts" / "run2_Za_BDT.py"))
    parser.add_argument("--new-script", default=str(repo_root / "HZaMVA" / "scripts" / "run3_Za_BDT.py"))
    parser.add_argument("--apply-script", default=str(repo_root / "Parquet2Rootfile" / "Parque2Root_BDT.py"))
    parser.add_argument("--prepare-script", default=str(repo_root / "Parquet2Rootfile" / "2_prepare_rootfile.sh"))
    parser.add_argument("--legacy-model-path", default=str(repo_root / "HZaMVA" / "keep" / "model_Za_BDT_run3.pkl"))
    parser.add_argument("--current-model-path", default=str(repo_root / "HZaMVA" / "using" / "model_Za_BDT_run3.pkl"))
    parser.add_argument("--output-dir", default=str(repo_root / "debug_bdt_regression"))
    parser.add_argument("--masses", nargs="*", type=int, default=DEFAULT_MASSES)
    parser.add_argument("--metric", choices=["asimov", "s_over_sqrt_b"], default="asimov")
    parser.add_argument("--threshold-step", type=float, default=DEFAULT_THRESHOLD_STEP)
    parser.add_argument("--max-events", type=int, default=None, help="Optional entry_stop when reading ROOT trees.")
    parser.add_argument("--strict-feature-order", action="store_true", default=True)
    parser.add_argument("--skip-runtime", action="store_true", help="Write only static/code-level diagnostics.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir))
    plots_dir = ensure_dir(output_dir / "plots")

    old_script = build_script_snapshot("old_run2", Path(args.old_script))
    new_script = build_script_snapshot("current_run3", Path(args.new_script))
    application_features = extract_application_features(Path(args.apply_script))
    application_scan_masses = extract_application_scan_masses(Path(args.apply_script))
    prepare_meta = extract_prepare_script_metadata(Path(args.prepare_script))

    legacy_model_cfg = extract_model_config_from_pickle(Path(args.legacy_model_path))
    current_model_cfg = extract_model_config_from_pickle(Path(args.current_model_path))

    context = RuntimeContext(
        output_dir=output_dir,
        plots_dir=plots_dir,
        old_script=old_script,
        new_script=new_script,
        application_features=application_features,
        application_scan_masses=application_scan_masses,
        legacy_model_cfg=legacy_model_cfg,
        current_model_cfg=current_model_cfg,
    )
    context.log(f"Output directory: {output_dir}")
    context.log(f"Skip runtime mode: {args.skip_runtime}")
    context.log(f"Old script: {args.old_script}")
    context.log(f"New script: {args.new_script}")
    context.log(f"Application script: {args.apply_script}")
    context.log(f"Application scan masses: {application_scan_masses or '<not detected>'}")

    training_features = new_script.variables + ["param"]
    ok_order, mismatch_messages = compare_feature_orders(training_features, application_features)
    if not ok_order:
        for message in mismatch_messages:
            context.warn(message)
        if args.strict_feature_order:
            context.issues.append(
                Issue(
                    severity=0,
                    title="Feature list mismatch would raise at model application time",
                    kind="static",
                    evidence=mismatch_messages,
                )
            )

    summary_table = reference_table_dataframe(REFERENCE_OLD_RESULTS, REFERENCE_CURRENT_OBSERVED)
    if pd is not None:
        summary_rows = summary_table.to_dict("records")
        summary_table.to_csv(output_dir / "summary_comparison.csv", index=False)
    else:
        summary_rows = list(summary_table)
        write_csv_rows(output_dir / "summary_comparison.csv", summary_rows)

    add_static_issues(context, prepare_meta, summary_table)

    specs = build_old_specs(old_script, args.masses) + build_new_specs(new_script, args.masses)
    sample_manifest_table = build_sample_manifest(specs)
    if pd is not None:
        sample_manifest_rows = sample_manifest_table.to_dict("records")
        sample_manifest_table.to_csv(output_dir / "sample_path_comparison.csv", index=False)
    else:
        sample_manifest_rows = list(sample_manifest_table)
        write_csv_rows(output_dir / "sample_path_comparison.csv", sample_manifest_rows)

    runtime_notes: List[str] = []
    sample_df = pd.DataFrame() if pd is not None else None
    feature_df = pd.DataFrame() if pd is not None else None
    split_df = pd.DataFrame() if pd is not None else None
    cutflow_df = pd.DataFrame() if pd is not None else None
    threshold_runtime_df = pd.DataFrame() if pd is not None else None
    channel_df = pd.DataFrame() if pd is not None else None
    threshold_payload: Dict[int, Dict[str, Any]] = {}

    if args.skip_runtime:
        runtime_notes.append("Runtime diagnostics were skipped by command-line option.")
        context.log("Runtime diagnostics skipped by --skip-runtime.")
    elif pd is None:
        runtime_notes.append("Runtime diagnostics were skipped because pandas is not available in the current environment.")
        context.log("Runtime diagnostics skipped because pandas is unavailable.")
    else:
        try:
            sample_df, feature_df, split_df, cutflow_df = runtime_sample_and_feature_summaries(
                context=context,
                specs=specs,
                masses=args.masses,
                entry_stop=args.max_events,
                strict_feature_order=args.strict_feature_order,
            )
            if sample_df.empty or sample_df["path_exists"].fillna(False).sum() == 0:
                runtime_notes.append("No ROOT inputs were available in the current environment; only code-level diagnostics were produced.")
            else:
                runtime_notes.append("Sample and feature summaries were computed from available ROOT inputs.")
        except Exception as exc:
            runtime_notes.append(f"Sample/feature runtime diagnostics failed: {exc}")
            context.warn(f"Sample/feature runtime diagnostics failed: {exc}")

        try:
            threshold_payload, threshold_runtime_df, channel_df = runtime_mass_score_diagnostics(
                context=context,
                masses=args.masses,
                entry_stop=args.max_events,
                metric=args.metric,
                threshold_step=args.threshold_step,
                strict_feature_order=args.strict_feature_order,
            )
            if threshold_payload:
                runtime_notes.append("Mass-by-mass score, smoothing, and threshold diagnostics were produced.")
                for mass, payload in threshold_payload.items():
                    if mass in REFERENCE_OLD_RESULTS:
                        old_ref = REFERENCE_OLD_RESULTS[mass]
                        if math.isfinite(payload.get("background_yield", float("nan"))) and old_ref["bkg_yield"] > 0 and payload["background_yield"] / old_ref["bkg_yield"] > 2.0:
                            context.warn(
                                f"mA={mass}: background yield {payload['background_yield']:.3f} exceeds 2x the old reference {old_ref['bkg_yield']:.3f}"
                            )
                        if math.isfinite(payload.get("signal_efficiency", float("nan"))) and old_ref["signal_eff"] > 0 and payload["signal_efficiency"] / old_ref["signal_eff"] < 0.8:
                            context.warn(
                                f"mA={mass}: signal efficiency {payload['signal_efficiency']:.3f} is below 80% of the old reference {old_ref['signal_eff']:.3f}"
                            )
            else:
                runtime_notes.append("Score-level diagnostics were skipped because the model or scored inputs were unavailable.")
        except Exception as exc:
            runtime_notes.append(f"Score/threshold runtime diagnostics failed: {exc}")
            context.warn(f"Score/threshold runtime diagnostics failed: {exc}")

    if pd is not None and sample_df is not None and not sample_df.empty:
        sample_df.to_csv(output_dir / "sample_yield_summary.csv", index=False)
    else:
        ensure_csv_with_header(
            output_dir / "sample_yield_summary.csv",
            [
                "study",
                "sample",
                "year",
                "process",
                "split",
                "n_entries",
                "sum_weight",
                "mean_weight",
                "negative_weight_fraction",
            ],
        )

    if pd is not None and feature_df is not None and not feature_df.empty:
        feature_df.to_csv(output_dir / "feature_summary.csv", index=False)
    else:
        ensure_csv_with_header(
            output_dir / "feature_summary.csv",
            ["study", "mass", "feature", "sample_type", "mean", "std", "min", "max", "nan_fraction", "inf_fraction", "zero_fraction", "separation"],
        )

    if pd is not None and split_df is not None and not split_df.empty:
        split_df.to_csv(output_dir / "split_summary.csv", index=False)
    else:
        ensure_csv_with_header(
            output_dir / "split_summary.csv",
            ["study", "sample", "process", "split", "tree_name", "configured_selection", "test_size", "random_state", "stratify"],
        )

    if pd is not None and cutflow_df is not None and not cutflow_df.empty:
        cutflow_df.to_csv(output_dir / "cutflow_comparison.csv", index=False)
    else:
        ensure_csv_with_header(
            output_dir / "cutflow_comparison.csv",
            ["study", "label", "step_index", "cut", "n_entries", "sum_weight"],
        )

    if pd is not None and channel_df is not None and not channel_df.empty:
        channel_df.to_csv(output_dir / "channel_summary.csv", index=False)
    else:
        ensure_csv_with_header(
            output_dir / "channel_summary.csv",
            ["mA", "channel", "signal_efficiency", "background_yield"],
        )

    threshold_scan_columns = [
        "mA",
        "threshold",
        "signal_efficiency",
        "signal_efficiency_unweighted",
        "signal_yield",
        "signal_yield_unweighted",
        "background_yield_raw",
        "background_yield_smooth",
        "background_yield_smooth_reference_raw",
        "data_sideband_yield",
        "metric",
    ]
    for mass in args.masses:
        threshold_path = output_dir / f"threshold_scan_mass{mass}.csv"
        if not threshold_path.exists():
            ensure_csv_with_header(threshold_path, threshold_scan_columns)

    hyperparameters_payload = {
        "old_script": asdict(old_script),
        "new_script": asdict(new_script),
        "application_feature_order": application_features,
        "application_scan_masses": application_scan_masses,
        "prepare_script": prepare_meta,
        "legacy_model_pickle": legacy_model_cfg,
        "current_model_pickle": current_model_cfg,
        "feature_order_match": ok_order,
        "feature_order_mismatches": mismatch_messages,
    }
    json_dump(output_dir / "hyperparameters.json", hyperparameters_payload)

    if pd is not None and threshold_runtime_df is not None and not threshold_runtime_df.empty:
        for mass, payload in threshold_payload.items():
            if mass not in summary_table["mA"].values:
                continue
            summary_table.loc[summary_table["mA"] == mass, "new_threshold"] = payload.get("best_threshold")
            summary_table.loc[summary_table["mA"] == mass, "new_signal_eff"] = payload.get("signal_efficiency")
            summary_table.loc[summary_table["mA"] == mass, "new_bkg_yield"] = payload.get("background_yield")
        summary_table["threshold_diff"] = summary_table["new_threshold"] - summary_table["old_threshold"]
        summary_table["signal_eff_ratio"] = summary_table["new_signal_eff"] / summary_table["old_signal_eff"]
        summary_table["bkg_yield_ratio"] = summary_table["new_bkg_yield"] / summary_table["old_bkg_yield"]
        warning_flags = []
        for _, row in summary_table.iterrows():
            row_warnings = []
            if pd.notna(row["bkg_yield_ratio"]) and row["bkg_yield_ratio"] > 2.0:
                row_warnings.append(DEFAULT_STATIC_WARNINGS["bkg_ratio_high"])
            if pd.notna(row["signal_eff_ratio"]) and row["signal_eff_ratio"] < 0.8:
                row_warnings.append(DEFAULT_STATIC_WARNINGS["sig_eff_low"])
            warning_flags.append(";".join(row_warnings))
        summary_table["warning_flag"] = warning_flags
        summary_table.to_csv(output_dir / "summary_comparison.csv", index=False)

    runtime_sections = {"runtime_notes": runtime_notes}
    write_markdown_report(context, runtime_sections, summary_table if pd is not None else summary_rows)
    (output_dir / "diagnostic_log.txt").write_text("\n".join(context.logs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
