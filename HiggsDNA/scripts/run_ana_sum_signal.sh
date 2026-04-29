#!/usr/bin/env bash

outdir="/eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Sig_MC"
unretire_jobs="${UNRETIRE_JOBS:-1}" # 1 re-run unfinished jobs # 0 no re-run
retire_jobs="${RETIRE_JOBS:-0}" # 1 merged parquet files # 0 no merged until all jobs finished
merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge

# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Sig_MC/mA_M5_2023preBPix
# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Sig_MC/analysis_manager.pkl

cmd=(
    python
    scripts/run_analysis.py
    --config "metadata/za_sum_signal_run3.json"
    --log-level "INFO"
    --n_cores 15
    --output_dir "${outdir}"
    --batch_system "condor"
)

if [[ "${retire_jobs}" == "1" ]]; then
    cmd+=(--retire_jobs)
elif [[ "${unretire_jobs}" == "1" ]]; then
    cmd+=(--unretire_jobs)
fi

if [[ "${merge_outputs}" == "1" ]]; then
    cmd+=(--merge_outputs)
fi

"${cmd[@]}" #--short #--with_skimmed # local condor
