#!/usr/bin/env python3
"""Create HTCondor submit files for split data/MC plotting jobs."""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path


DEFAULT_PROJECT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna"
REGIONS = ("SR", "CR", "mva")
FINAL_TAGS = ("nominal", "sideband_rwgt")
SAMPLES = (
    ("data", "Data"),
    ("dyll", "DYJetsToLL"),
    ("dyg", "DYGto2LG"),
    ("sig_m1", "M1"),
    ("sig_m2", "M2"),
    ("sig_m3", "M3"),
    ("sig_m4", "M4"),
    ("sig_m5", "M5"),
    ("sig_m6", "M6"),
    ("sig_m7", "M7"),
    ("sig_m8", "M8"),
    ("sig_m9", "M9"),
    ("sig_m10", "M10"),
    ("sig_m15", "M15"),
    ("sig_m20", "M20"),
    ("sig_m25", "M25"),
    ("sig_m30", "M30"),
)


@dataclass(frozen=True)
class CondorJob:
    region_key: str
    final_tag: str
    sample_tag: str
    samples: str

    def as_queue_row(self) -> str:
        return f"{self.region_key} {self.final_tag} {self.sample_tag} {self.samples}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate one-CPU HTCondor jobs for Plot/scripts/1_prepare_dataVmc.py. "
            "Each Data/DY/signal mass sample is submitted as a separate job."
        )
    )
    parser.add_argument(
        "--project-dir",
        default=DEFAULT_PROJECT_DIR,
        help="absolute HiggsZaAna path visible to Condor workers",
    )
    parser.add_argument(
        "--condor-dir",
        default=None,
        help="directory where dataVmc.submit, dataVmc_jobs.txt, and logs are written",
    )
    parser.add_argument(
        "--submit-name",
        default="dataVmc.submit",
        help="HTCondor submit file name",
    )
    parser.add_argument(
        "--jobs-name",
        default="dataVmc_jobs.txt",
        help="queue item file name",
    )
    parser.add_argument(
        "--request-memory",
        default="4000MB",
        help="Condor request_memory value",
    )
    parser.add_argument(
        "--job-flavour",
        default="workday",
        help="CERN HTCondor +JobFlavour value",
    )
    parser.add_argument(
        "--max-events",
        default="",
        help="optional MAX_EVENTS value passed to each job",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="run condor_submit after generating files",
    )
    return parser.parse_args()


def build_jobs() -> list[CondorJob]:
    return [
        CondorJob(region_key=region_key, final_tag=final_tag, sample_tag=sample_tag, samples=samples)
        for final_tag in FINAL_TAGS
        for region_key in REGIONS
        for sample_tag, samples in SAMPLES
    ]


def write_jobs_file(path: Path, jobs: list[CondorJob]) -> None:
    path.write_text("\n".join(job.as_queue_row() for job in jobs) + "\n", encoding="utf-8")


def write_submit_file(
    path: Path,
    *,
    project_dir: Path,
    condor_dir: Path,
    jobs_file: Path,
    request_memory: str,
    job_flavour: str,
    max_events: str,
) -> None:
    executable = project_dir / "Plot" / "Condor" / "run_dataVmc_condor_job.sh"
    log_dir = condor_dir / "logs"
    environment_parts = [f"PROJECT_DIR={project_dir}"]
    if max_events:
        environment_parts.append(f"MAX_EVENTS={max_events}")
    environment = " ".join(environment_parts)

    content = f"""universe = vanilla
executable = {executable}
arguments = "$(region_key) $(final_tag) $(sample_tag) $(samples)"
initialdir = {condor_dir}
getenv = True
should_transfer_files = NO

request_cpus = 1
request_memory = {request_memory}
+JobFlavour = "{job_flavour}"

environment = "{environment}"

output = {log_dir}/$(final_tag)_$(region_key)_$(sample_tag).$(ClusterId).$(ProcId).out
error  = {log_dir}/$(final_tag)_$(region_key)_$(sample_tag).$(ClusterId).$(ProcId).err
log    = {log_dir}/$(final_tag)_$(region_key)_$(sample_tag).$(ClusterId).$(ProcId).log

queue region_key, final_tag, sample_tag, samples from {jobs_file}
"""
    path.write_text(content, encoding="utf-8")


def main() -> int:
    args = parse_args()
    project_dir = Path(args.project_dir).resolve()
    condor_dir = Path(args.condor_dir).resolve() if args.condor_dir else project_dir / "Plot" / "Condor"
    jobs_file = condor_dir / args.jobs_name
    submit_file = condor_dir / args.submit_name

    condor_dir.mkdir(parents=True, exist_ok=True)
    (condor_dir / "logs").mkdir(parents=True, exist_ok=True)

    jobs = build_jobs()
    write_jobs_file(jobs_file, jobs)
    write_submit_file(
        submit_file,
        project_dir=project_dir,
        condor_dir=condor_dir,
        jobs_file=jobs_file,
        request_memory=args.request_memory,
        job_flavour=args.job_flavour,
        max_events=str(args.max_events).strip(),
    )

    print(f"Wrote {len(jobs)} jobs to {jobs_file}")
    print(f"Wrote submit file to {submit_file}")
    print("Submit with:")
    print(f"  condor_submit {submit_file}")

    if args.submit:
        subprocess.run(["condor_submit", str(submit_file)], check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
