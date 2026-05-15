import os
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

import awkward as ak
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from higgs_dna.constants import CENTRAL_WEIGHT
from higgs_dna.job_management.jobs import CondorJob
from higgs_dna.job_management.managers import JobsManager
from higgs_dna.job_management.task import (
    Task,
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

    def test_iter_parquet_row_groups_retries_transient_io_error(self):
        with tempfile.TemporaryDirectory(prefix="task-merge-retry-") as tmpdir:
            path = os.path.join(tmpdir, "job.parquet")
            pq.write_table(
                pa.table({"event": pa.array([1], type=pa.int64())}),
                path,
            )

            real_parquet_file = pq.ParquetFile
            failures_remaining = [1]

            class FlakyParquetFile:
                def __init__(self, file_path):
                    self._real = real_parquet_file(file_path)
                    self.num_row_groups = self._real.num_row_groups
                    self.schema_arrow = self._real.schema_arrow

                def read_row_group(self, row_group_idx):
                    if failures_remaining[0] > 0:
                        failures_remaining[0] -= 1
                        raise OSError("[Errno 5] Input/output error")
                    return self._real.read_row_group(row_group_idx)

            with mock.patch.dict(
                os.environ,
                {
                    "HIGGSDNA_PARQUET_READ_RETRIES": "3",
                    "HIGGSDNA_PARQUET_READ_RETRY_DELAY": "0",
                },
            ):
                with mock.patch(
                    "higgs_dna.job_management.task.pq.ParquetFile",
                    side_effect=lambda file_path: FlakyParquetFile(file_path),
                ):
                    with mock.patch("higgs_dna.job_management.task.time.sleep") as sleep:
                        tables = list(_iter_parquet_row_groups(path))

            self.assertEqual(len(tables), 1)
            self.assertEqual(tables[0].column("event").to_pylist(), [1])
            sleep.assert_called_once_with(0.0)

    def test_jobs_manager_complete_only_depends_on_task_completion(self):
        manager = object.__new__(JobsManager)
        manager.tasks = [
            SimpleNamespace(complete=True, merged_output_files=False),
            SimpleNamespace(complete=True, merged_output_files=True),
        ]
        self.assertTrue(JobsManager.complete(manager))

    def test_condor_request_memory_can_be_overridden_from_environment(self):
        with mock.patch.dict(os.environ, {"HIGGSDNA_CONDOR_REQ_MEMORY": "24000"}):
            requests = CondorJob.get_requests()

        self.assertEqual(requests["REQ_MEMORY"], 24000)
        self.assertEqual(CondorJob.REQUESTS["REQ_MEMORY"], 20000)

    def test_task_merge_outputs_skips_existing_merged_files(self):
        with tempfile.TemporaryDirectory(prefix="task-merge-skip-") as tmpdir:
            input_path = os.path.join(tmpdir, "job_1_nominal.parquet")
            merged_path = os.path.join(tmpdir, "merged_nominal.parquet")

            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([1], type=pa.int64()),
                        CENTRAL_WEIGHT: pa.array([1.0], type=pa.float64()),
                    }
                ),
                input_path,
            )
            pq.write_table(
                pa.table(
                    {
                        "event": pa.array([999], type=pa.int64()),
                        CENTRAL_WEIGHT: pa.array([999.0], type=pa.float64()),
                    }
                ),
                merged_path,
            )

            task = object.__new__(Task)
            task.name = "unit"
            task.output_dir = tmpdir
            task.outputs = {"nominal": [input_path]}
            task.merged_outputs = {}
            task.config = {"sample": {"is_data": False}}
            task.scale1fb = 1.0
            task.lumi = 1.0
            task.wrote_process_ids = True
            task.wrote_years = True
            task.merged_output_files = False
            task.remerge = False

            Task.merge_outputs(task)

            merged = pq.read_table(merged_path)
            self.assertEqual(merged.column("event").to_pylist(), [999])
            self.assertTrue(task.merged_output_files)
            self.assertEqual(task.merged_outputs["nominal"], merged_path)


if __name__ == "__main__":
    unittest.main()
