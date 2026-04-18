#!/usr/bin/env python3
"""
Standalone parquet merger for HiggsDNA outputs.

Expected task layout:

    <root>/<task_name>/job_<n>/
        <task_name>_config_job<n>.json
        <task_name>_summary_job<n>.json
        output_job_<n>_<syst>.parquet

This script can:
    1. Merge per-job parquet files inside each task directory into
       <task_dir>/merged_<syst>.parquet.
    2. Optionally merge all task-level merged parquet files into
       <root>/merged_<syst>.parquet.

By default, existing merged files are preserved. Use --force to overwrite them.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from higgs_dna.constants import CENTRAL_WEIGHT
from higgs_dna.job_management.task import (
    _get_merged_output_schema,
    _iter_parquet_row_groups,
    _normalize_table_to_schema,
    _replace_or_add_column,
)


LOGGER = logging.getLogger("merge_parquet_outputs")
OUTPUT_PATTERN = re.compile(r"output_job_\d+_(.+)\.parquet$")
YEAR_PATTERN = re.compile(r"(20\d{2}(?:preEE|postEE|preBPix|postBPix)?)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge HiggsDNA parquet outputs.")
    parser.add_argument(
        "root",
        type=Path,
        help="Task directory or parent directory containing task directories.",
    )
    parser.add_argument(
        "--task",
        dest="tasks",
        action="append",
        default=[],
        help="Only merge the named task. Can be given multiple times.",
    )
    parser.add_argument(
        "--systematic",
        dest="systematics",
        action="append",
        default=[],
        help="Only merge the named systematic. Can be given multiple times.",
    )
    parser.add_argument(
        "--top-level",
        action="store_true",
        help="Also merge task-level merged parquet files into <root>/merged_<syst>.parquet.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing merged parquet files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be merged without writing parquet outputs.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s %(message)s",
    )


def is_task_dir(path: Path) -> bool:
    if not path.is_dir():
        return False
    return any(child.is_dir() and child.name.startswith("job_") for child in path.iterdir())


def discover_task_dirs(root: Path, task_filters: list[str]) -> list[Path]:
    if is_task_dir(root):
        task_dirs = [root]
    else:
        task_dirs = sorted(path for path in root.iterdir() if is_task_dir(path))

    if task_filters:
        allowed = set(task_filters)
        task_dirs = [path for path in task_dirs if path.name in allowed]

    return task_dirs


def load_json(path: Path, *, retries: int = 1, retry_delay: float = 0.2) -> dict[str, Any]:
    last_error = None
    for attempt in range(retries):
        try:
            with path.open("r") as handle:
                return json.load(handle)
        except (json.JSONDecodeError, OSError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(retry_delay)

    assert last_error is not None
    raise last_error


def find_first_json(task_dir: Path, pattern: str) -> Path | None:
    for path in sorted(task_dir.glob(f"job_*/{pattern}")):
        return path
    return None


def infer_year(task_name: str, sample_config: dict[str, Any]) -> str | None:
    year = sample_config.get("year")
    if year is not None:
        return str(year)

    match = YEAR_PATTERN.search(task_name)
    if match is not None:
        return match.group(1)
    return None


def collect_task_metadata(task_dir: Path) -> dict[str, Any]:
    task_name = task_dir.name
    config_path = find_first_json(task_dir, "*_config_job*.json")
    config = load_json(config_path, retries=3) if config_path is not None else {}
    sample = config.get("sample", {})
    is_data = bool(sample.get("is_data", False))
    summaries = sorted(task_dir.glob("job_*/*_summary_job*.json"))

    outputs_by_syst: dict[str, list[Path]] = defaultdict(list)
    sum_weights = 0.0
    bad_summaries: list[Path] = []

    for parquet_path in sorted(task_dir.glob("job_*/output_job_*_*.parquet")):
        match = OUTPUT_PATTERN.match(parquet_path.name)
        if match is None:
            continue
        outputs_by_syst[match.group(1)].append(parquet_path)

    for summary_path in summaries:
        try:
            summary = load_json(summary_path, retries=3)
        except (json.JSONDecodeError, OSError) as exc:
            bad_summaries.append(summary_path)
            LOGGER.warning(
                "Task '%s' has unreadable summary '%s': %s. Continuing without it.",
                task_name,
                summary_path,
                exc,
            )
            continue

        sum_weights += float(summary.get("sum_weights", 0.0))

    if bad_summaries and not is_data:
        LOGGER.warning(
            "Task '%s' skipped %d unreadable summary file(s). "
            "If this task still needs a fresh MC merge, scale1fb may be underestimated.",
            task_name,
            len(bad_summaries),
        )

    metadata = {
        "task_name": task_name,
        "task_dir": task_dir,
        "config": config,
        "sample": sample,
        "is_data": is_data,
        "norm_factor": sample.get("norm_factor"),
        "lumi": sample.get("lumi"),
        "process_id": sample.get("process_id"),
        "year": infer_year(task_name, sample),
        "sum_weights": sum_weights,
        "bad_summaries": [str(path) for path in bad_summaries],
        "outputs_by_syst": dict(outputs_by_syst),
    }

    if not is_data and metadata["norm_factor"] is not None and sum_weights > 0.0:
        metadata["scale1fb"] = (float(metadata["norm_factor"]) * 1000.0) / sum_weights
    else:
        metadata["scale1fb"] = 0.0

    return metadata


def replace_or_add_field(schema: pa.Schema, field: pa.Field) -> pa.Schema:
    field_idx = schema.get_field_index(field.name)
    if field_idx >= 0:
        return schema.set(field_idx, field)
    return schema.append(field)


def build_constant_column(value: Any, length: int) -> pa.Array:
    scalar = pa.scalar(value)
    return pa.array([value] * length, type=scalar.type)


def merge_parquet_files(
    inputs: list[Path],
    output_path: Path,
    task_name: str,
    syst_tag: str,
    *,
    promote_weight_fields: bool,
    scale_weights: bool,
    scale1fb: float,
    lumi: float,
    extra_fields: dict[str, Any] | None,
    force: bool,
    dry_run: bool,
) -> bool:
    if not inputs:
        return False

    if output_path.exists() and not force:
        LOGGER.info(
            "Task '%s' : merged output '%s' already exists for syst '%s', skipping.",
            task_name,
            output_path,
            syst_tag,
        )
        return False

    if output_path.exists() and force:
        LOGGER.info(
            "Task '%s' : merged output '%s' already exists for syst '%s', overwriting.",
            task_name,
            output_path,
            syst_tag,
        )

    LOGGER.info(
        "Task '%s' : merging %d parquet files into '%s' for syst '%s'.",
        task_name,
        len(inputs),
        output_path,
        syst_tag,
    )
    if dry_run:
        return True

    merged_schema = _get_merged_output_schema(
        [str(path) for path in inputs],
        is_data=not promote_weight_fields,
        task_name=task_name,
        syst_tag=syst_tag,
    )

    if extra_fields:
        for name, value in extra_fields.items():
            merged_schema = replace_or_add_field(
                merged_schema,
                pa.field(name, pa.scalar(value).type),
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    writer = None
    try:
        for index, input_path in enumerate(inputs, start=1):
            if index == 1 or index == len(inputs) or index % 25 == 0:
                LOGGER.info(
                    "Task '%s' : processing parquet %d/%d for syst '%s'.",
                    task_name,
                    index,
                    len(inputs),
                    syst_tag,
                )

            for table in _iter_parquet_row_groups(str(input_path)):
                if scale_weights:
                    central_weight = pc.cast(table.column(CENTRAL_WEIGHT), pa.float64())
                    table = _replace_or_add_column(
                        table,
                        CENTRAL_WEIGHT,
                        pc.multiply(
                            central_weight,
                            pa.scalar(scale1fb * lumi, type=pa.float64()),
                        ),
                    )
                    table = _replace_or_add_column(
                        table,
                        CENTRAL_WEIGHT + "_no_lumi",
                        pc.multiply(
                            central_weight,
                            pa.scalar(scale1fb, type=pa.float64()),
                        ),
                    )

                if extra_fields:
                    for name, value in extra_fields.items():
                        table = _replace_or_add_column(
                            table,
                            name,
                            build_constant_column(value, len(table)),
                        )

                table = _normalize_table_to_schema(
                    table,
                    merged_schema,
                    task_name,
                    syst_tag,
                    str(input_path),
                )

                if writer is None:
                    writer = pq.ParquetWriter(str(output_path), merged_schema)
                writer.write_table(table)
    finally:
        if writer is not None:
            LOGGER.info(
                "Task '%s' : closing merged parquet for syst '%s'.",
                task_name,
                syst_tag,
            )
            writer.close()
            LOGGER.info(
                "Task '%s' : finished merged parquet for syst '%s'.",
                task_name,
                syst_tag,
            )

    return True


def merge_task_dir(
    task_dir: Path,
    *,
    systematics_filter: set[str],
    force: bool,
    dry_run: bool,
) -> dict[str, Any]:
    metadata = collect_task_metadata(task_dir)
    outputs_by_syst = metadata["outputs_by_syst"]
    merged_outputs: dict[str, Path] = {}

    if not outputs_by_syst:
        LOGGER.warning("Task '%s' has no parquet outputs to merge.", metadata["task_name"])
        metadata["merged_outputs"] = merged_outputs
        return metadata

    for syst_tag in sorted(outputs_by_syst):
        if systematics_filter and syst_tag not in systematics_filter:
            continue

        merged_output = task_dir / f"merged_{syst_tag}.parquet"
        merge_parquet_files(
            outputs_by_syst[syst_tag],
            merged_output,
            metadata["task_name"],
            syst_tag,
            promote_weight_fields=not metadata["is_data"],
            scale_weights=not metadata["is_data"],
            scale1fb=float(metadata["scale1fb"]),
            lumi=float(metadata["lumi"] or 0.0),
            extra_fields=None,
            force=force,
            dry_run=dry_run,
        )
        merged_outputs[syst_tag] = merged_output

    metadata["merged_outputs"] = merged_outputs
    return metadata


def merge_task_outputs_to_root(
    root: Path,
    task_results: list[dict[str, Any]],
    *,
    systematics_filter: set[str],
    force: bool,
    dry_run: bool,
) -> None:
    grouped: dict[str, list[tuple[Path, dict[str, Any]]]] = defaultdict(list)

    for task in task_results:
        for syst_tag, merged_output in task.get("merged_outputs", {}).items():
            if systematics_filter and syst_tag not in systematics_filter:
                continue
            if merged_output.exists() or dry_run:
                grouped[syst_tag].append((merged_output, task))

    promote_weight_fields = any(not task["is_data"] for task in task_results)

    for syst_tag, entries in sorted(grouped.items()):
        inputs: list[Path] = []
        per_file_fields: dict[Path, dict[str, Any]] = {}
        for merged_output, task in entries:
            inputs.append(merged_output)
            extra_fields = {}
            if task.get("process_id") is not None:
                extra_fields["process_id"] = task["process_id"]
            if task.get("year") is not None:
                extra_fields["year"] = str(task["year"])
            per_file_fields[merged_output] = extra_fields

        merged_path = root / f"merged_{syst_tag}.parquet"
        if merged_path.exists() and not force:
            LOGGER.info(
                "Root merged output '%s' already exists for syst '%s', skipping.",
                merged_path,
                syst_tag,
            )
            continue

        if dry_run:
            LOGGER.info(
                "Would merge %d task-level parquet files into '%s'.",
                len(inputs),
                merged_path,
            )
            continue

        merged_schema = _get_merged_output_schema(
            [str(path) for path in inputs],
            is_data=not promote_weight_fields,
            task_name=root.name,
            syst_tag=syst_tag,
        )
        for extra_fields in per_file_fields.values():
            for name, value in extra_fields.items():
                merged_schema = replace_or_add_field(
                    merged_schema,
                    pa.field(name, pa.scalar(value).type),
                )

        if merged_path.exists():
            merged_path.unlink()

        writer = None
        try:
            for index, input_path in enumerate(inputs, start=1):
                if index == 1 or index == len(inputs) or index % 25 == 0:
                    LOGGER.info(
                        "Root merge: processing parquet %d/%d for syst '%s'.",
                        index,
                        len(inputs),
                        syst_tag,
                    )

                for table in _iter_parquet_row_groups(str(input_path)):
                    for name, value in per_file_fields[input_path].items():
                        table = _replace_or_add_column(
                            table,
                            name,
                            build_constant_column(value, len(table)),
                        )

                    table = _normalize_table_to_schema(
                        table,
                        merged_schema,
                        root.name,
                        syst_tag,
                        str(input_path),
                    )

                    if writer is None:
                        writer = pq.ParquetWriter(str(merged_path), merged_schema)
                    writer.write_table(table)
        finally:
            if writer is not None:
                writer.close()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    root = args.root.resolve()
    if not root.exists():
        raise FileNotFoundError(root)

    task_dirs = discover_task_dirs(root, args.tasks)
    if not task_dirs:
        LOGGER.error("No task directories found under '%s'.", root)
        return 1

    do_top_level_merge = args.top_level
    if is_task_dir(root) and args.top_level:
        LOGGER.warning(
            "'%s' is already a task directory; skipping --top-level merge.",
            root,
        )
        do_top_level_merge = False

    systematics_filter = set(args.systematics)
    task_results = []
    for task_dir in task_dirs:
        task_results.append(
            merge_task_dir(
                task_dir,
                systematics_filter=systematics_filter,
                force=args.force,
                dry_run=args.dry_run,
            )
        )

    if do_top_level_merge:
        merge_task_outputs_to_root(
            root,
            task_results,
            systematics_filter=systematics_filter,
            force=args.force,
            dry_run=args.dry_run,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
