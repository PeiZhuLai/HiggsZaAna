#!/usr/bin/env python3
"""Utilities for applying sideband reweight JSON files."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


DEFAULT_REWEIGHT_JSON = str(Path(__file__).resolve().parents[1] / "reweights" / "sideband_run3_iterative.json")
DEFAULT_REWEIGHT_JSON = os.environ.get("HZA_SIDEBAND_REWEIGHT_JSON", DEFAULT_REWEIGHT_JSON)

VAR_ALIASES = {
    "H_m": ("H_m", "H_mass"),
    "ALP_m": ("ALP_m", "ALP_mass"),
    "pho1ECALIso": ("pho1ECALIso", "pho1PIso_noCorr", "ALP_lead_photon_ecalPFClusterIso"),
    "pho2ECALIso": ("pho2ECALIso", "pho2PIso_noCorr", "ALP_sublead_photon_ecalPFClusterIso"),
}


def _finite_float(value, default=np.nan):
    try:
        value = float(value)
    except Exception:
        return default
    if not math.isfinite(value):
        return default
    return value


def _lookup(values, edges, factors):
    values = np.asarray(values, dtype=float)
    edges = np.asarray(edges, dtype=float)
    factors = np.asarray(factors, dtype=float)
    if len(edges) != len(factors) + 1:
        raise ValueError("Invalid reweight binning: len(edges) must equal len(factors) + 1")

    idx = np.searchsorted(edges, values, side="right") - 1
    idx = np.clip(idx, 0, len(factors) - 1)
    out = factors[idx]
    return np.where(np.isfinite(values) & np.isfinite(out), out, 1.0)


def _mass_from_event(events, seed, masses):
    events = np.asarray(events, dtype=np.int64)
    masses = np.asarray(masses, dtype=float)
    selected = np.empty(events.shape[0], dtype=float)
    for idx, event in enumerate(events):
        payload = f"{seed}:{int(event)}".encode("utf-8", errors="ignore")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        selected[idx] = masses[int.from_bytes(digest, "little") % len(masses)]
    return selected


def _mass_from_row_order(n_rows, seed, masses):
    rng = np.random.default_rng(seed)
    return rng.choice(np.asarray(masses, dtype=float), size=n_rows, replace=True)


def _first_present_dataframe_column(frame, logical_name):
    for candidate in VAR_ALIASES.get(logical_name, (logical_name,)):
        if candidate in frame.columns:
            return candidate
    return None


def _first_present_object_attr(obj, logical_name):
    for candidate in VAR_ALIASES.get(logical_name, (logical_name,)):
        if hasattr(obj, candidate):
            return candidate
    return None


class SidebandReweighter:
    def __init__(self, payload: Mapping[str, object], source_path: str | None = None):
        self.payload = dict(payload)
        self.source_path = source_path
        settings = self.payload.get("settings", {})
        param = self.payload.get("param", {})
        self.reweight_vars = list(settings.get("reweight_vars", []))
        self.seed = int(settings.get("seed", 12345))
        self.mass_hypotheses = list(param.get("mass_hypotheses", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0]))
        self.param_method = str(param.get("method", "event_hash"))
        self.iterations = list(self.payload.get("iterations", []))

    @classmethod
    def from_json(cls, path):
        path = Path(path)
        with path.open() as handle:
            payload = json.load(handle)
        return cls(payload, source_path=str(path))

    @classmethod
    def from_json_if_exists(cls, path=DEFAULT_REWEIGHT_JSON):
        path = Path(path)
        if not path.exists():
            return None
        return cls.from_json(path)

    def mass_hypotheses_for_dataframe(self, frame):
        event_col = _first_present_dataframe_column(frame, "event")
        if self.param_method == "event_hash" and event_col is not None:
            events = frame[event_col].fillna(-1).to_numpy(dtype=np.int64)
            return _mass_from_event(events, self.seed, self.mass_hypotheses)
        return _mass_from_row_order(len(frame), self.seed, self.mass_hypotheses)

    def param_values(self, frame, output_col="param", mass_output_col=None):
        h_col = _first_present_dataframe_column(frame, "H_m")
        alp_col = _first_present_dataframe_column(frame, "ALP_m")
        if h_col is None or alp_col is None:
            raise KeyError("Need H_m/H_mass and ALP_m/ALP_mass to build sideband param.")

        h_m = frame[h_col].to_numpy(dtype=float)
        alp_m = frame[alp_col].to_numpy(dtype=float)
        mass_hyp = self.mass_hypotheses_for_dataframe(frame)

        with np.errstate(divide="ignore", invalid="ignore"):
            param = (alp_m - mass_hyp) / h_m
        if output_col is not None:
            frame[output_col] = param
        if mass_output_col is not None:
            frame[mass_output_col] = mass_hyp
        return param

    def ensure_param(self, frame, output_col="param", mass_output_col=None):
        return self.param_values(frame, output_col=output_col, mass_output_col=mass_output_col)

    def _dataframe_values(self, frame, var):
        if var == "param":
            if "param" not in frame.columns:
                return self.param_values(frame, output_col="param")
            return frame["param"].to_numpy(dtype=float)

        column = _first_present_dataframe_column(frame, var)
        if column is None:
            raise KeyError(f"Missing sideband reweight variable '{var}' in DataFrame.")
        return frame[column].to_numpy(dtype=float)

    def weights_for_dataframe(self, frame):
        weights = np.ones(len(frame), dtype=float)
        for iteration in self.iterations:
            for step in iteration.get("steps", []):
                var = step["var"]
                values = self._dataframe_values(frame, var)
                factors = _lookup(values, step["edges"], step["clipped_factors"])
                norm = _finite_float(step.get("normalization_scale", 1.0), default=1.0)
                weights *= factors * norm
        return np.where(np.isfinite(weights), weights, 1.0)

    def apply_to_dataframe(self, frame, base_weight_col="factor", output_weight_col="factor", output_reweight_col="sideband_rwgt"):
        rwgt = self.weights_for_dataframe(frame)
        frame[output_reweight_col] = rwgt
        if base_weight_col is not None and output_weight_col is not None:
            frame[output_weight_col] = frame[base_weight_col].to_numpy(dtype=float) * rwgt
        return frame

    def _object_value(self, obj, var, row_index=None):
        if var == "param":
            attr = _first_present_object_attr(obj, "param")
            if attr is not None:
                return _finite_float(getattr(obj, attr))

            h_attr = _first_present_object_attr(obj, "H_m")
            alp_attr = _first_present_object_attr(obj, "ALP_m")
            if h_attr is None or alp_attr is None:
                return np.nan

            h_m = _finite_float(getattr(obj, h_attr))
            alp_m = _finite_float(getattr(obj, alp_attr))
            if not math.isfinite(h_m) or h_m == 0.0 or not math.isfinite(alp_m):
                return np.nan

            event_attr = _first_present_object_attr(obj, "event")
            if self.param_method == "event_hash" and event_attr is not None:
                mass_hyp = _mass_from_event(np.asarray([getattr(obj, event_attr)], dtype=np.int64), self.seed, self.mass_hypotheses)[0]
            else:
                idx = 0 if row_index is None else int(row_index)
                mass_hyp = _mass_from_row_order(idx + 1, self.seed, self.mass_hypotheses)[idx]
            return (alp_m - mass_hyp) / h_m

        attr = _first_present_object_attr(obj, var)
        if attr is None:
            return np.nan
        return _finite_float(getattr(obj, attr))

    def weight_for_object(self, obj, row_index=None):
        weight = 1.0
        for iteration in self.iterations:
            for step in iteration.get("steps", []):
                value = self._object_value(obj, step["var"], row_index=row_index)
                factor = _lookup(np.asarray([value], dtype=float), step["edges"], step["clipped_factors"])[0]
                norm = _finite_float(step.get("normalization_scale", 1.0), default=1.0)
                weight *= factor * norm
        return float(weight) if math.isfinite(weight) else 1.0


def load_sideband_reweighter(path=DEFAULT_REWEIGHT_JSON, required=False):
    reweighter = SidebandReweighter.from_json_if_exists(path)
    if reweighter is None and required:
        raise FileNotFoundError(f"Sideband reweight JSON not found: {path}")
    return reweighter
