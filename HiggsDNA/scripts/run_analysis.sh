#!/usr/bin/env bash

merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge

cmd=(
    python
    scripts/run_analysis.py
    --config "metadata/zgamma.json"
    --log-level "DEBUG"
    --n_cores 5
    --output_dir "cutflow/"
    --unretire_jobs
)

if [[ "${merge_outputs}" == "1" ]]; then
    cmd+=(--merge_outputs)
fi

"${cmd[@]}"
