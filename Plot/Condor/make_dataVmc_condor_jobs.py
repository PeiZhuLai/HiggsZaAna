#!/usr/bin/env python3
"""Create HTCondor submit files for split data/MC plotting jobs."""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path


DEFAULT_PROJECT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna"
REGIONS = ("SR", "CR", "mva")
KNOWN_FINAL_TAGS = ("nominal", "sideband_rwgt")
DEFAULT_FINAL_TAGS = ("sideband_rwgt",)
DYLL_YEARS = (
    "2022preEE",
    "2022postEE",
    "2023preBPix",
    "2023postBPix",
    "2024",
)
DYG_SPLITS = (
    ("DYGto2LG_10to50", ("2022preEE", "2022postEE")),
    ("DYGto2LG_50to100", ("2022preEE", "2022postEE")),
    ("DYGto2LG_10to100", ("2023preBPix", "2023postBPix", "2024")),
)
SIGNAL_MASSES = (
    "M1",
    "M2",
    "M3",
    "M4",
    "M5",
    "M6",
    "M7",
    "M8",
    "M9",
    "M10",
    "M15",
    "M20",
    "M25",
    "M30",
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
            "Data and signal samples stay one job each, while DY backgrounds are "
            "split by year or (subsample, year)."
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
        "--request-disk",
        default="20000MB",
        help="Condor request_disk value; needed when localizing the conda env to worker scratch",
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
        "--final-tags",
        default=",".join(DEFAULT_FINAL_TAGS),
        help=(
            "comma-separated output tags to run. Default is sideband_rwgt "
            "(sideband reweighting with ordinary systematics, but without "
            "sideband reweight uncertainty). Use 'all' or 'nominal,sideband_rwgt' "
            "to also run the unweighted nominal jobs."
        ),
    )
    parser.add_argument(
        "--localize-conda-env",
        default="auto",
        choices=("auto", "1", "0"),
        help="extract a prebuilt conda tarball to worker scratch before running; auto localizes only when CONDA_PREFIX is under /eos",
    )
    parser.add_argument(
        "--setup-conda-env",
        default="auto",
        choices=("auto", "1", "0"),
        help="run use-anaconda, anaconda, and conda activate inside each Condor job",
    )
    parser.add_argument(
        "--conda-env-name",
        default="higgs-alp-ana",
        help="conda environment name activated inside each Condor job",
    )
    parser.add_argument(
        "--anaconda-setup",
        default="/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh",
        help="anaconda setup script sourced before conda activate",
    )
    parser.add_argument(
        "--conda-tarball",
        default=None,
        help="prebuilt conda env tarball used by Condor jobs",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="run condor_submit after generating files",
    )
    return parser.parse_args()


def parse_final_tags(raw_final_tags: str) -> tuple[str, ...]:
    tags: list[str] = []
    for token in str(raw_final_tags).replace(";", ",").split(","):
        tag = token.strip()
        if not tag:
            continue
        if tag.lower() in ("all", "*"):
            tags = list(KNOWN_FINAL_TAGS)
            break
        if tag not in KNOWN_FINAL_TAGS:
            allowed = ", ".join(KNOWN_FINAL_TAGS)
            raise ValueError(f"Unknown final tag '{tag}'. Allowed tags are: {allowed}, all")
        if tag not in tags:
            tags.append(tag)

    if not tags:
        raise ValueError("--final-tags selected no tags")
    return tuple(tags)


def build_sample_specs() -> list[tuple[str, str]]:
    sample_specs: list[tuple[str, str]] = [("data", "Data")]

    for year in DYLL_YEARS:
        sample_specs.append((f"dyll_{year}", f"DYJetsToLL@{year}"))

    for subsample, years in DYG_SPLITS:
        short_tag = subsample.replace("DYGto2LG_", "")
        for year in years:
            sample_specs.append((f"dyg_{short_tag}_{year}", f"DYGto2LG@{subsample}:{year}"))

    for mass in SIGNAL_MASSES:
        sample_specs.append((f"sig_{mass.lower()}", mass))

    return sample_specs


def build_jobs(final_tags: tuple[str, ...]) -> list[CondorJob]:
    return [
        CondorJob(region_key=region_key, final_tag=final_tag, sample_tag=sample_tag, samples=samples)
        for final_tag in final_tags
        for region_key in REGIONS
        for sample_tag, samples in build_sample_specs()
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
    request_disk: str,
    job_flavour: str,
    max_events: str,
    localize_conda_env: str,
    setup_conda_env: str,
    conda_env_name: str,
    anaconda_setup: str,
    conda_tarball: str,
) -> None:
    executable = project_dir / "Plot" / "Condor" / "run_dataVmc_condor_job.sh"
    log_dir = condor_dir / "logs"
    environment_parts = [f"PROJECT_DIR={project_dir}"]
    if max_events:
        environment_parts.append(f"MAX_EVENTS={max_events}")
    if localize_conda_env:
        environment_parts.append(f"LOCALIZE_CONDA_ENV={localize_conda_env}")
    if setup_conda_env:
        environment_parts.append(f"SETUP_CONDA_ENV={setup_conda_env}")
    if conda_env_name:
        environment_parts.append(f"CONDA_ENV_NAME={conda_env_name}")
    if anaconda_setup:
        environment_parts.append(f"ANACONDA_SETUP={anaconda_setup}")
    if conda_tarball:
        environment_parts.append(f"CONDA_TARBALL={conda_tarball}")
    environment = " ".join(environment_parts)

    content = f"""universe = vanilla
executable = {executable}
arguments = "$(region_key) $(final_tag) $(sample_tag) $(samples)"
initialdir = {condor_dir}
getenv = True
should_transfer_files = NO

request_cpus = 1
request_memory = {request_memory}
request_disk = {request_disk}
+JobFlavour = "{job_flavour}"

environment = "{environment}"

output = {log_dir}/$(final_tag)_$(region_key)_$(sample_tag).$(ClusterId).$(ProcId).out
error  = {log_dir}/$(final_tag)_$(region_key)_$(sample_tag).$(ClusterId).$(ProcId).err
log    = {log_dir}/dataVmc.$(ClusterId).log

queue region_key, final_tag, sample_tag, samples from {jobs_file}
"""
    path.write_text(content, encoding="utf-8")


def main() -> int:
    args = parse_args()
    try:
        final_tags = parse_final_tags(args.final_tags)
    except ValueError as exc:
        print(f"[ERROR] {exc}")
        return 2

    project_dir = Path(args.project_dir).resolve()
    condor_dir = Path(args.condor_dir).resolve() if args.condor_dir else project_dir / "Plot" / "Condor"
    jobs_file = condor_dir / args.jobs_name
    submit_file = condor_dir / args.submit_name
    conda_tarball = (
        Path(args.conda_tarball).resolve()
        if args.conda_tarball
        else condor_dir / "env_cache" / f"{args.conda_env_name}.tar.gz"
    )

    condor_dir.mkdir(parents=True, exist_ok=True)
    (condor_dir / "logs").mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(final_tags)
    write_jobs_file(jobs_file, jobs)
    write_submit_file(
        submit_file,
        project_dir=project_dir,
        condor_dir=condor_dir,
        jobs_file=jobs_file,
        request_memory=args.request_memory,
        request_disk=args.request_disk,
        job_flavour=args.job_flavour,
        max_events=str(args.max_events).strip(),
        localize_conda_env=str(args.localize_conda_env).strip(),
        setup_conda_env=str(args.setup_conda_env).strip(),
        conda_env_name=str(args.conda_env_name).strip(),
        anaconda_setup=str(args.anaconda_setup).strip(),
        conda_tarball=str(conda_tarball),
    )

    print(f"Final tags: {', '.join(final_tags)}")
    print(f"Wrote {len(jobs)} jobs to {jobs_file}")
    print(f"Wrote submit file to {submit_file}")
    print("Submit with:")
    print(f"  condor_submit {submit_file}")

    if args.submit:
        subprocess.run(["condor_submit", str(submit_file)], check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
