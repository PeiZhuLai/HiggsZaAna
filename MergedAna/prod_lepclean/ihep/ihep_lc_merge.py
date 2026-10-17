"""Merge + reconcile the IHEP lepclean friend chunks, one line per tag.

Layout from _ihep_lc_file_worker.sh:
    <OUT>/<tag>/chunks/f<idx>/<tag>_<era>/merged_nominal.parquet   (one per nano file)
    <OUT>/<tag>/chunks/f<idx>/EMPTY                                 (0 selected events)
Writes <OUT>/<tag>/merged/merged_nominal.parquet (the extra merged/ level matches the
CERN glob Data_202*_*/*/merged_nominal.parquet), streamed via .tmp + rename.

PASS iff every file idx 0..N-1 has a chunk or EMPTY, merged rows == sum(chunk rows),
no duplicated run:lumi:event, and the n_MLPhoton_preclean column exists.

    python ihep_lc_merge.py <catalog.json> <OUT> [tag ...]      (default: all tags)
Exit 0 only if every requested tag passes.
"""
import json
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


def merge_tag(tag, n, out):
    cdir = os.path.join(out, tag, "chunks")
    chunks, empty = [], set()
    for f in sorted(os.listdir(cdir)) if os.path.isdir(cdir) else []:
        if not f.startswith("f"):
            continue
        fd = os.path.join(cdir, f)
        found = [os.path.join(fd, s, "merged_nominal.parquet") for s in os.listdir(fd)
                 if os.path.exists(os.path.join(fd, s, "merged_nominal.parquet"))]
        if found:
            chunks.append((int(f[1:]), found[0]))
        elif os.path.exists(os.path.join(fd, "EMPTY")):
            empty.add(int(f[1:]))
    have = {i for i, _ in chunks} | empty
    missing = sorted(set(range(n)) - have)
    if missing:
        return False, f"FAIL {tag} files={n} chunks={len(chunks)} empty={len(empty)} missing={len(missing)}{missing[:20]}"
    paths = [p for _, p in sorted(chunks)]
    chunk_rows = sum(pq.ParquetFile(p).metadata.num_rows for p in paths)
    schema = pa.unify_schemas([pq.read_schema(p) for p in paths])
    dest_dir = os.path.join(out, tag, "merged")
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, "merged_nominal.parquet")
    with pq.ParquetWriter(dest + ".tmp", schema) as w:
        for p in paths:
            t = pq.read_table(p)
            w.write_table(t.select(schema.names).cast(schema) if t.schema != schema else t)
    os.replace(dest + ".tmp", dest)
    merged_rows = pq.ParquetFile(dest).metadata.num_rows
    t = pq.read_table(dest, columns=[c for c in ("run", "luminosityBlock", "event", "n_MLPhoton_preclean")
                                     if c in schema.names])
    dup = -1
    if all(c in t.column_names for c in ("run", "luminosityBlock", "event")):
        k = np.stack([t[c].to_numpy().astype(np.int64) for c in ("run", "luminosityBlock", "event")], 1)
        dup = len(k) - len(np.unique(k, axis=0))
    nojoin = (f"{np.mean(t['n_MLPhoton_preclean'].to_numpy() == 0):.2%}"
              if "n_MLPhoton_preclean" in t.column_names else "NO_COLUMN")
    ok = merged_rows == chunk_rows and dup == 0 and nojoin != "NO_COLUMN"
    return ok, (f"{'PASS' if ok else 'FAIL'} {tag} files={n} chunks={len(chunks)} empty={len(empty)} "
                f"chunk_rows={chunk_rows} merged_rows={merged_rows} dup={dup} nojoin={nojoin}")


def main():
    cat = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    tags = sys.argv[3:] or sorted(cat)
    nbad = 0
    for tag in tags:
        n = len(list(cat[tag]["files"].values())[0])
        try:
            ok, line = merge_tag(tag, n, out)
        except Exception as e:  # one unreadable chunk must not hide the other tags
            ok, line = False, f"FAIL {tag} exception {type(e).__name__}: {e}"
        nbad += not ok
        print(line, flush=True)
    print(f"SUMMARY tags={len(tags)} fail={nbad}", flush=True)
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
