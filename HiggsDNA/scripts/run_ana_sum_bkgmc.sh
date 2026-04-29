#!/usr/bin/env bash

outdir="/eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC"
unretire_jobs="${UNRETIRE_JOBS:-1}" # 1 re-run unfinished jobs # 0 no re-run
retire_jobs="${RETIRE_JOBS:-0}" # 1 merged parquet files # 0 no merged until all jobs finished
merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge

# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/DYGto2LG_10to100_2023preBPix/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/DYJetsToLL_2023preBPix/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/analysis_manager.pkl

cmd=(
    python
    scripts/run_analysis.py
    --config "metadata/za_sum_bkgmc_run3.json"
    --log-level "INFO"
    --n_cores 10
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

"${cmd[@]}" #--short #--with_skimmed  #--batch_system "local" condor
