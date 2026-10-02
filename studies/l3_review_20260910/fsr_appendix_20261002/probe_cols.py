import pyarrow.parquet as pq, sys
for f in sys.argv[1:]:
    pf=pq.ParquetFile(f); names=pf.schema_arrow.names
    print(f, pf.metadata.num_rows)
    print([n for n in names if any(k in n.lower() for k in ("fsr","z_mass","h_mass","z_m","h_m","weight_central","z_mumu","z_ee","mass"))][:80])
