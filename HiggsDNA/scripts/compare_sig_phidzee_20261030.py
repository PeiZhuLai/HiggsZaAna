#!/usr/bin/env python3
"""Compare a re-produced signal sample with the previous production (parquet_DNA_tmp_sfunc_ord).

The photon-ID (Z->ee) SF re-production must select exactly the same events: the selection did not
change, only the scale-factor JSONs. So for every merged_*.parquet:
  * the (run, luminosityBlock, event) sets must be identical,
  * every non-weight column must be identical after sorting by event id,
  * weight columns may differ; their new/old ratios are summarized.
Exit code 1 if any sample fails the event-identity or non-weight checks.

usage: compare_sig_phidzee_20261030.py <new_Sig_MC_dir> [sample_year ...]
"""
import os
import sys

import numpy as np
import pyarrow.parquet as pq

OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_sfunc_ord/Sig_MC"  # fsrfix deleted 10-29; sfunc_ord = same events, 10-24 SFs
KEY = ["run", "luminosityBlock", "event"]
# Known, intended differences: Hessian PDF band (2026-10-02 attach_lhe_weights fix) and
# process_id (index of the sample in the run's sample_list).
EXPECTED_DIFF = {"LHEPdfWeight_Up", "LHEPdfWeight_Down", "process_id"}
# FSR-dressed lepton "mass": sqrt(E^2 - p^2) of a ~50 GeV object, so float rounding is amplified to
# O(0.01-0.3) GeV (smoke test: muons only, while Z_mass moved by <1.3e-3 GeV). Reported, not failed.
NOISY = {"Z_lead_lepton_mass", "Z_sublead_lepton_mass"}
# Float-level noise seen in the smoke test (relative <= 1.5e-5 on kinematics, CPU rounding):
REL_TOL = 1e-4
ABS_TOL = 5e-3  # floor for small-valued columns (ALP_mass differed by up to 2.6e-3 GeV in the smoke test)


def load(path):
    df = pq.read_table(path).to_pandas()
    return df.sort_values(KEY).reset_index(drop=True)


def compare(new_dir, sy, fname):
    a = load(os.path.join(OLD, sy, fname))
    b = load(os.path.join(new_dir, sy, fname))
    problems = []
    if len(a) != len(b) or not (a[KEY].values == b[KEY].values).all():
        return ["event sets differ: old %d new %d" % (len(a), len(b))], {}
    if set(a.columns) != set(b.columns):
        problems.append("columns differ: only old %s, only new %s"
                        % (sorted(set(a.columns) - set(b.columns))[:5], sorted(set(b.columns) - set(a.columns))[:5]))
    ratios = {}
    for c in sorted(set(a.columns) & set(b.columns)):
        x, y = a[c].values, b[c].values
        if c.startswith("weight"):
            ok = (x != 0) & np.isfinite(x) & np.isfinite(y)
            r = y[ok] / x[ok]
            if r.size and not np.allclose(r, 1.0, rtol=1e-9):
                ratios[c] = (float(np.median(r)), float(r.min()), float(r.max()), int((x == 0).sum()), int((y == 0).sum()))
            continue
        if c in EXPECTED_DIFF:
            continue
        try:
            same = np.array_equal(x, y, equal_nan=True)
            if not same and np.issubdtype(x.dtype, np.number):
                xf, yf = x.astype(float), y.astype(float)
                d = np.abs(yf - xf)
                ok = (d <= np.maximum(REL_TOL * np.abs(xf), ABS_TOL)) | (np.isnan(xf) & np.isnan(yf))
                same = bool(np.all(ok))
        except TypeError:
            same = all(np.array_equal(np.asarray(u), np.asarray(v)) for u, v in zip(x, y))
        if not same and c in NOISY:
            print("    (noise) %s differs in %d events" % (c, int((np.abs(np.asarray(y, float) - np.asarray(x, float)) > ABS_TOL).sum())))
            continue
        if not same:
            xf = np.asarray(x, dtype=float) if np.issubdtype(np.asarray(x).dtype, np.number) else None
            if xf is not None:
                d = np.abs(np.asarray(y, dtype=float) - xf)
                tol = np.maximum(REL_TOL * np.abs(xf), ABS_TOL)
                problems.append("column %s differs (max abs %.3g, n>tol %d)" % (c, np.nanmax(d), int((d > tol).sum())))
            else:
                problems.append("column %s differs" % c)
    return problems, ratios


def main():
    new_dir = sys.argv[1]
    samples = sys.argv[2:] or sorted(d for d in os.listdir(new_dir) if os.path.isdir(os.path.join(new_dir, d)))
    nbad = 0
    for sy in samples:
        files = sorted(f for f in os.listdir(os.path.join(new_dir, sy)) if f.startswith("merged_") and f.endswith(".parquet"))
        if not files:
            print("!! %s: no merged parquet" % sy)
            nbad += 1
            continue
        for fname in files:
            problems, ratios = compare(new_dir, sy, fname)
            tag = "FAIL" if problems else "ok"
            print("[%s] %s/%s" % (tag, sy, fname))
            for p in problems:
                print("    !! " + p)
            if fname == "merged_nominal.parquet":
                for c, (med, lo, hi, z_old, z_new) in ratios.items():
                    print("    %-62s new/old median %.4f range [%.4f, %.4f] zeros old %d new %d" % (c, med, lo, hi, z_old, z_new))
            nbad += bool(problems)
    print("SUMMARY bad=%d" % nbad)
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
