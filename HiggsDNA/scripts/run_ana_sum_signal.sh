#!/usr/bin/env bash

outdir="/eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Sig_MC"
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
    --unretire_jobs
    --batch_system "condor"
)

if [[ "${merge_outputs}" == "1" ]]; then
    cmd+=(--merge_outputs)
fi

"${cmd[@]}" #--short #--with_skimmed # local condor
