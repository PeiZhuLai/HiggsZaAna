import os
import tempfile
import unittest
from types import SimpleNamespace

import awkward as ak
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from higgs_dna.constants import CENTRAL_WEIGHT
from higgs_dna.job_management.managers import JobsManager
from higgs_dna.job_management.task import (
    _get_merged_output_schema,
    _iter_parquet_row_groups,
    _normalize_table_to_schema,
)


class TestTaskMergeSchema(unittest.TestCase):
    def test_merge_schema_normalizes_missing_fields_and_weight_precision(self):
        with tempfile.TemporaryDirectory(prefix="task-merge-schema-") as tmpdir:
            path_a = os.path.join(tmpdir, "job_a.parquet")
            path_b = os.path.join(tmpdir, "job_b.parquet")

            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([1, 2], type=pa.int64()),
                        CENTRAL_WEIGHT: pa.array([1.0, 2.0], type=pa.float32()),
                        "weight_extra": pa.array([0.9, 1.1], type=pa.float32()),
                        "probe_pt": pa.array([20.0, 30.0], type=pa.float64()),
                    }
                ),
                path_a,
            )
            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([3], type=pa.int64()),
                        CENTRAL_WEIGHT: pa.array([3.0], type=pa.float64()),
                        "probe_pt": pa.array([40.0], type=pa.float64()),
                    }
                ),
                path_b,
            )

            merged_schema = _get_merged_output_schema(
                [path_a, path_b],
                is_data=False,
                task_name="unit",
                syst_tag="Muon_scale_down",
            )

            self.assertEqual(merged_schema.field(CENTRAL_WEIGHT).type, pa.float64())
            self.assertEqual(merged_schema.field("weight_extra").type, pa.float64())
            self.assertGreaterEqual(
                merged_schema.get_field_index(CENTRAL_WEIGHT + "_no_lumi"),
                0,
            )

            merged_path = os.path.join(tmpdir, "merged.parquet")
            writer = pq.ParquetWriter(merged_path, merged_schema)
            try:
                for path in [path_a, path_b]:
                    for table in _iter_parquet_row_groups(path):
                        weight = pc.cast(table.column(CENTRAL_WEIGHT), pa.float64())
                        table = table.set_column(
                            table.schema.get_field_index(CENTRAL_WEIGHT),
                            CENTRAL_WEIGHT,
                            weight,
                        )
                        table = table.append_column(
                            CENTRAL_WEIGHT + "_no_lumi",
                            weight,
                        )
                        table = _normalize_table_to_schema(
                            table,
                            merged_schema,
                            "unit",
                            "Muon_scale_down",
                            path,
                        )
                        writer.write_table(table)
            finally:
                writer.close()

            merged = pq.read_table(merged_path)
            self.assertEqual(merged.num_rows, 3)
            self.assertEqual(
                merged.column("weight_extra").to_pylist(),
                [0.8999999761581421, 1.100000023841858, 1.0],
            )
            self.assertEqual(
                merged.column(CENTRAL_WEIGHT + "_no_lumi").to_pylist(),
                [1.0, 2.0, 3.0],
            )

    def test_merge_schema_reports_incompatible_field_types(self):
        with tempfile.TemporaryDirectory(prefix="task-merge-schema-bad-") as tmpdir:
            path_a = os.path.join(tmpdir, "job_a.parquet")
            path_b = os.path.join(tmpdir, "job_b.parquet")

            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([1], type=pa.int64()),
                        "flag": pa.array([1], type=pa.int64()),
                    }
                ),
                path_a,
            )
            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([2], type=pa.int64()),
                        "flag": pa.array(["bad"], type=pa.string()),
                    }
                ),
                path_b,
            )

            with self.assertRaisesRegex(RuntimeError, "Conflicting fields: flag"):
                _get_merged_output_schema(
                    [path_a, path_b],
                    is_data=True,
                    task_name="unit",
                    syst_tag="nominal",
                )

    def test_merge_schema_normalizes_awkward_extension_types(self):
        with tempfile.TemporaryDirectory(prefix="task-merge-schema-ext-") as tmpdir:
            path_native = os.path.join(tmpdir, "job_native.parquet")
            path_ext = os.path.join(tmpdir, "job_ext.parquet")

            pq.write_table(
                pa.table(
                    {
                        "fixedGridRhoAll": pa.array([10.0], type=pa.float32()),
                        "event_nPV": pa.array([12], type=pa.uint8()),
                    }
                ),
                path_native,
            )
            pq.write_table(
                ak.to_arrow_table(
                    ak.Array(
                        {
                            "fixedGridRhoAll": np.array([11.0], dtype=np.float32),
                            "event_nPV": np.array([13], dtype=np.uint8),
                        }
                    ),
                    extensionarray=True,
                ),
                path_ext,
            )

            merged_schema = _get_merged_output_schema(
                [path_native, path_ext],
                is_data=True,
                task_name="unit",
                syst_tag="Muon_scale_down",
            )

            self.assertEqual(merged_schema.field("fixedGridRhoAll").type, pa.float32())
            self.assertEqual(merged_schema.field("event_nPV").type, pa.uint8())

    def test_jobs_manager_complete_only_depends_on_task_completion(self):
        manager = object.__new__(JobsManager)
        manager.tasks = [
            SimpleNamespace(complete=True, merged_output_files=False),
            SimpleNamespace(complete=True, merged_output_files=True),
        ]
        self.assertTrue(JobsManager.complete(manager))


if __name__ == "__main__":
    unittest.main()
