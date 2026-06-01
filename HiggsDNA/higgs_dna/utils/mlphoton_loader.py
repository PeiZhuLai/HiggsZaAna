"""
MLPhoton friend-tree loader for the slim MLNanoAOD production.

The slim MLNanoAOD ROOT files contain only event IDs (run, luminosityBlock,
event, bunchCrossing, orbitNumber) + the MLPhoton flat table. They are keyed
on EOS by the MiniAODv6 UUID, e.g.:

    /eos/.../MLNanoAOD/<tag>/<tag>_<MINI-UUID>.root

Because central NanoAODv15 production is not 1:1 with MiniAODv6, we cannot
pair files by index. Instead we:

1. Read the NanoAODv15 file's parent MiniAOD UUIDs from a pre-built map
   (see HiggsDNA/scripts/build_parent_map.py).
2. Open the corresponding MLNanoAOD files for those UUIDs.
3. JOIN by (run, luminosityBlock, event) and attach MLPhoton_* branches as
   the events.MLPhoton awkward record onto the main events array.

Public API: ``attach_mlphoton(events, nano_url, parent_map, ml_dir, tag)``.
"""
from __future__ import annotations

import glob
import json
import os
from typing import Iterable

import awkward as ak
import numpy as np
import uproot


# Branches we pull out of each MLNanoAOD file. Must exactly match the slim
# producer (Prod_MLNanoAOD_run3_{mc,data}.py).
_ID_BRANCHES = ["run", "luminosityBlock", "event"]
_ML_BRANCHES = [
    "nMLPhoton",
    "MLPhoton_pt", "MLPhoton_eta", "MLPhoton_phi", "MLPhoton_mass",
    "MLPhoton_massEnergyRatio",
    "MLPhoton_diphotonScore", "MLPhoton_monophotonScore", "MLPhoton_hadronScore",
    "MLPhoton_pfIsolation",
    "MLPhoton_r1", "MLPhoton_r2", "MLPhoton_r3",
]


def _normalize_lfn(url_or_lfn: str) -> str:
    """Strip xrootd prefix and store-prefix-only the path so it matches DAS LFNs."""
    s = url_or_lfn
    if "://" in s:
        s = s.split("://", 1)[1]
        # Drop the host: e.g. cms-xrd-global.cern.ch//store/...
        if "/" in s:
            s = "/" + s.split("/", 1)[1]
    # Collapse the double slash that sometimes follows the host
    while "//" in s:
        s = s.replace("//", "/")
    return s


def _ml_path_for(uuid: str, ml_dir: str, tag: str) -> str | None:
    """Resolve <ml_dir>/<tag>_<UUID>.root, return None if missing."""
    candidate = os.path.join(ml_dir, f"{tag}_{uuid}.root")
    if os.path.isfile(candidate):
        return candidate
    # Tolerate slight name variants by glob (e.g. legacy naming)
    matches = glob.glob(os.path.join(ml_dir, f"*_{uuid}*.root"))
    return matches[0] if matches else None


def _load_ml_chunk(ml_path: str) -> ak.Array | None:
    """Open one MLNanoAOD file and return an ak record with event IDs + MLPhoton.

    Returns None on missing file or empty tree."""
    if ml_path is None:
        return None
    try:
        with uproot.open(ml_path) as f:
            tree = f["Events"]
            arrs = tree.arrays(_ID_BRANCHES + _ML_BRANCHES, library="ak")
    except Exception as e:
        # Caller logs; we just bail out
        return None
    if len(arrs) == 0:
        return None
    return arrs


def _build_event_key(run, lumi, event) -> np.ndarray:
    """Build a 64-bit composite key (run, lumi, event). Assumes:
    - run fits in 32 bits (true for CMS data)
    - lumi fits in 24 bits
    - event fits in 64 bits but is well below 2^40 for our use; we pack
      run<<40 | lumi<<24 | (event & 0xFFFFFF), then if collisions matter fall
      back to a tuple-based lookup.

    For safety we go with a structured (run, lumi, event) tuple → numpy
    rec-array, then hash via numpy.lexsort.
    """
    run_a = np.asarray(run, dtype=np.int64)
    lumi_a = np.asarray(lumi, dtype=np.int64)
    evt_a = np.asarray(event, dtype=np.int64)
    return run_a, lumi_a, evt_a


def attach_mlphoton(
    events: ak.Array,
    nano_url: str,
    parent_map: dict[str, list[str]],
    ml_dir: str,
    tag: str,
    logger=None,
) -> ak.Array:
    """Attach an events.MLPhoton record to ``events`` by event-ID JOIN.

    Parameters
    ----------
    events : ak.Array
        The main NanoAODv15 events array (already loaded by HiggsDNA, with
        run / luminosityBlock / event fields).
    nano_url : str
        The NanoAODv15 file URL (xrootd or local path) used to look up parents.
    parent_map : dict[str, list[str]]
        nano_lfn -> [MiniAOD UUID, ...] mapping built by build_parent_map.py.
    ml_dir : str
        Directory holding MLNanoAOD outputs (e.g. /eos/.../MLNanoAOD/<tag>).
    tag : str
        Output prefix used when writing MLNanoAOD (so file names are
        ``<tag>_<UUID>.root``).
    logger : optional logger

    Returns
    -------
    ak.Array
        ``events`` augmented with an ``MLPhoton`` field (record-of-arrays).
        Events without a matching MLPhoton entry get zero-length MLPhoton.
    """
    def _log(msg, level="info"):
        if logger:
            getattr(logger, level, logger.info)(msg)

    nano_lfn = _normalize_lfn(nano_url)
    parents = parent_map.get(nano_lfn)
    if parents is None:
        # Try a looser match by basename — different DAS replicas may have
        # different LFN prefixes.
        base = os.path.basename(nano_lfn)
        for k in parent_map:
            if os.path.basename(k) == base:
                parents = parent_map[k]
                break
    if not parents:
        _log(f"[mlphoton] no parents for {nano_lfn}; injecting empty MLPhoton", "warning")
        return _attach_empty(events)

    ml_chunks = []
    for uuid in parents:
        ml_path = _ml_path_for(uuid, ml_dir, tag)
        if not ml_path:
            _log(f"[mlphoton] MLNanoAOD missing for parent {uuid} in {ml_dir}", "warning")
            continue
        chunk = _load_ml_chunk(ml_path)
        if chunk is None or len(chunk) == 0:
            continue
        ml_chunks.append(chunk)
    if not ml_chunks:
        _log(f"[mlphoton] no MLNanoAOD content for any parent of {nano_lfn}", "warning")
        return _attach_empty(events)

    ml_all = ak.concatenate(ml_chunks, axis=0)

    # Build lookup: (run, lumi, event) -> ML row index
    main_run, main_lumi, main_evt = _build_event_key(events.run, events.luminosityBlock, events.event)
    ml_run, ml_lumi, ml_evt = _build_event_key(ml_all.run, ml_all.luminosityBlock, ml_all.event)

    # numpy lexsort on the ML side, then searchsorted main keys
    ml_order = np.lexsort((ml_evt, ml_lumi, ml_run))
    ml_run_s = ml_run[ml_order]
    ml_lumi_s = ml_lumi[ml_order]
    ml_evt_s = ml_evt[ml_order]

    # For lookup we need exact match. Use 3-key search by composite tuple cmp.
    # Build composite int128-ish via two int64 keys (run<<22 | lumi) and (event):
    key_hi_ml = (ml_run_s.astype(np.int64) << 22) | (ml_lumi_s.astype(np.int64) & ((1 << 22) - 1))
    key_lo_ml = ml_evt_s.astype(np.int64)
    key_hi_main = (main_run.astype(np.int64) << 22) | (main_lumi.astype(np.int64) & ((1 << 22) - 1))
    key_lo_main = main_evt.astype(np.int64)

    # Sort ML side by (key_hi, key_lo) and search.
    combo_sort = np.lexsort((key_lo_ml, key_hi_ml))
    ml_order = ml_order[combo_sort]
    key_hi_ml = key_hi_ml[combo_sort]
    key_lo_ml = key_lo_ml[combo_sort]

    # Binary search per main event. Numpy's searchsorted needs single-key
    # arrays; we do it by treating (hi, lo) as a single 128-bit composite
    # via tuple comparison through a per-event Python loop fallback when N
    # is small, else use numpy's stable searchsorted on hi then refine.
    idx_hi = np.searchsorted(key_hi_ml, key_hi_main, side="left")
    n_ml = len(key_hi_ml)
    matched = np.full(len(main_run), -1, dtype=np.int64)
    for i in range(len(main_run)):
        j = idx_hi[i]
        while j < n_ml and key_hi_ml[j] == key_hi_main[i] and key_lo_ml[j] != key_lo_main[i]:
            j += 1
        if j < n_ml and key_hi_ml[j] == key_hi_main[i] and key_lo_ml[j] == key_lo_main[i]:
            matched[i] = ml_order[j]
    n_unmatched = int((matched < 0).sum())
    if n_unmatched:
        _log(f"[mlphoton] {n_unmatched}/{len(main_run)} main events missing MLPhoton match", "warning")

    # Build the MLPhoton record: variable-length per event
    n_ml_arr = ak.to_numpy(ml_all.nMLPhoton)
    # offsets in the original (unsorted) ml_all
    ml_offsets = np.concatenate([[0], np.cumsum(n_ml_arr)])
    # gather function: for each main event, find its ML row, take the
    # slice [ml_offsets[row]:ml_offsets[row+1]] of each MLPhoton_* branch.
    out_pt   = []
    out_eta  = []
    out_phi  = []
    out_mass = []
    out_mer  = []
    out_dip  = []
    out_mon  = []
    out_had  = []
    out_iso  = []
    out_r1   = []
    out_r2   = []
    out_r3   = []
    out_n    = np.zeros(len(main_run), dtype=np.int32)

    pt   = ak.to_numpy(ak.flatten(ml_all.MLPhoton_pt))
    eta  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_eta))
    phi  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_phi))
    mass = ak.to_numpy(ak.flatten(ml_all.MLPhoton_mass))
    mer  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_massEnergyRatio))
    dip  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_diphotonScore))
    mon  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_monophotonScore))
    had  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_hadronScore))
    iso  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_pfIsolation))
    r1a  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_r1))
    r2a  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_r2))
    r3a  = ak.to_numpy(ak.flatten(ml_all.MLPhoton_r3))

    for i, row in enumerate(matched):
        if row < 0:
            continue
        lo = int(ml_offsets[row])
        hi = int(ml_offsets[row + 1])
        out_n[i] = hi - lo
        out_pt.append(pt[lo:hi]);    out_eta.append(eta[lo:hi]);   out_phi.append(phi[lo:hi])
        out_mass.append(mass[lo:hi]);out_mer.append(mer[lo:hi])
        out_dip.append(dip[lo:hi]);  out_mon.append(mon[lo:hi]);   out_had.append(had[lo:hi])
        out_iso.append(iso[lo:hi])
        out_r1.append(r1a[lo:hi]);   out_r2.append(r2a[lo:hi]);    out_r3.append(r3a[lo:hi])

    # Some main events have no match → emit empty arrays for them
    def _build_jagged(values_list: list, n: np.ndarray):
        # Build a jagged ak.Array of length len(n), in order:
        counts = n.copy()
        # values_list contains arrays only for matched rows; unmatched rows
        # contributed nothing. So we interleave by counts.
        full = np.concatenate(values_list) if values_list else np.array([], dtype=np.float64)
        return ak.unflatten(full, counts)

    events_with_ml = ak.with_field(
        events,
        ak.zip({
            "pt":              _build_jagged(out_pt, out_n),
            "eta":             _build_jagged(out_eta, out_n),
            "phi":             _build_jagged(out_phi, out_n),
            "mass":            _build_jagged(out_mass, out_n),
            "massEnergyRatio": _build_jagged(out_mer, out_n),
            "diphotonScore":   _build_jagged(out_dip, out_n),
            "monophotonScore": _build_jagged(out_mon, out_n),
            "hadronScore":     _build_jagged(out_had, out_n),
            "pfIsolation":     _build_jagged(out_iso, out_n),
            "r1":              _build_jagged(out_r1, out_n),
            "r2":              _build_jagged(out_r2, out_n),
            "r3":              _build_jagged(out_r3, out_n),
        }),
        "MLPhoton",
    )
    _log(f"[mlphoton] attached {len(main_run)} events ({n_unmatched} unmatched)")
    return events_with_ml


def _attach_empty(events: ak.Array) -> ak.Array:
    """Fallback: attach an empty MLPhoton record to every event."""
    n = len(events)
    empty = ak.unflatten(np.array([], dtype=np.float32), np.zeros(n, dtype=np.int32))
    return ak.with_field(
        events,
        ak.zip({
            "pt":              empty, "eta":             empty, "phi":             empty,
            "mass":            empty, "massEnergyRatio": empty,
            "diphotonScore":   empty, "monophotonScore": empty, "hadronScore":     empty,
            "pfIsolation":     empty, "r1": empty, "r2": empty, "r3": empty,
        }),
        "MLPhoton",
    )


def load_parent_map(path: str) -> dict[str, list[str]]:
    """Convenience: load a parent map JSON (and accept either a list or a single string)."""
    with open(path) as f:
        d = json.load(f)
    return {k: ([v] if isinstance(v, str) else list(v)) for k, v in d.items()}
