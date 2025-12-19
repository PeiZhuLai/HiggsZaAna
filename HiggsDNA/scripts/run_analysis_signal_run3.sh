outdir="/eos/home-p/pelai/HZa/parquet_DNA/Sig_MC"

# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Sig_MC

python scripts/run_analysis.py --config "metadata/za_signal_run3.json" --log-level "DEBUG" --n_cores 15 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short #--with_skimmed # local condor 