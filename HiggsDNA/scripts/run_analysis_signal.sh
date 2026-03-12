# outdir="/eos/home-p/pelai/HZa/parquet_DNA/Sig_MC"
# For Gen Information
outdir="/eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC"

rm -fr /eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC/mA_M5_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC/mA_M5_2023postBPix/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC/mA_M1_2024
# rm -fr /eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC/mA_M1_2022postEE
rm -fr /eos/home-p/pelai/HZa/parquet_cutflow_DNA/Sig_MC/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_signal_run3.json" --log-level "DEBUG" --n_cores 15 --output_dir $outdir --unretire_jobs --batch_system "local" --short #--with_skimmed # local condor # INFO # DEBUG
