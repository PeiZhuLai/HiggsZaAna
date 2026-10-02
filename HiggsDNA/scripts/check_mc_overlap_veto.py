#!/usr/bin/env python3
"""Post-production check of the MC overlap removal (run after every merge, before p2root).

For each sample dir (the one holding merged_nominal.parquet and job_*/..._config_job*.json):
  * the role is taken from the input files in a job config, with the SAME classifier the tagger
    uses (higgs_dna/taggers/mc_overlap_roles.py) -- an unclassified sample is an error;
  * role 'veto' (the +jets side, e.g. DY+jets, ttbar): merged_nominal.parquet must have
    NO event with n_iso_photons > 0. n_iso_photons (gen photon pT>15, |eta|<2.6) is a subset
    of the tagger's definition (pT>10, any eta), so any such event means the veto was not applied.
  * role 'keep' / 'none': reported only (for 'keep', n_iso_photons==0 is legitimate: 10-15 GeV
    or forward photons).
Why: on 2026-09-27 the 2024 DY+jets samples turned out to carry 21-23% (by weight) of events with
an isolated gen photon -- the veto had silently not been applied because the dataset names were
not in the tagger's list. This check would have failed at the first merge.

Usage: check_mc_overlap_veto.py <sample_dir_or_parent> [...]
  A parent dir is expanded to all its subdirs that contain merged_nominal.parquet.
Exit 0 if every 'veto' sample is clean and every sample is classified, 1 otherwise.
"""
import glob, importlib.util, json, os, sys
import pyarrow.parquet as pq

_R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "higgs_dna", "taggers", "mc_overlap_roles.py")
_spec = importlib.util.spec_from_file_location("mc_overlap_roles", _R)
roles = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(roles)


def sample_dirs(paths):
    out = []
    for p in paths:
        p = p.rstrip("/")
        if os.path.exists(os.path.join(p, "merged_nominal.parquet")):
            out.append(p)
        else:
            out += sorted(os.path.dirname(m) for m in glob.glob(os.path.join(p, "*", "merged_nominal.parquet")))
    return out


def input_files(sd):
    # look at job_1..job_20 directly: globbing thousands of job dirs over EOS FUSE takes minutes
    for i in range(1, 21):
        for cfg in glob.glob(os.path.join(sd, "job_%d" % i, "*_config_job%d.json" % i)):
            files = json.load(open(cfg)).get("files") or []
            if files:
                return files
    return []


def main(argv):
    bad = 0
    for sd in sample_dirs(argv):
        name = os.path.basename(sd)
        if name.startswith("Data"):
            continue
        files = input_files(sd)
        if not files:
            print("NO_CONFIG  %-34s cannot determine inputs (no job config) -> FAIL" % name); bad += 1; continue
        try:
            rs = {roles.classify_overlap_role(f) for f in files}
        except RuntimeError as e:
            print("UNCLASSIFIED %-32s %s -> FAIL" % (name, str(e)[:160])); bad += 1; continue
        if len(rs) != 1:
            print("MIXED      %-34s roles %s -> FAIL" % (name, sorted(rs))); bad += 1; continue
        role = rs.pop()
        n = pq.read_table(os.path.join(sd, "merged_nominal.parquet"), columns=["n_iso_photons"]).column(0).to_numpy()
        npos = int((n > 0).sum())
        if role == "veto":
            ok = npos == 0
            bad += 0 if ok else 1
            print("%-10s %-34s role=veto n_iso_photons>0: %d / %d" % ("OK" if ok else "FAIL", name, npos, len(n)))
        else:
            print("%-10s %-34s role=%s n_iso_photons>0: %d / %d (info)" % ("OK", name, role, npos, len(n)))
    print("VERDICT: %s" % ("PASSED" if bad == 0 else "FAILED (%d)" % bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
