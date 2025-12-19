outdir="/eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC"

# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC

python scripts/run_analysis.py --config "metadata/za_bkgmc_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short #--with_skimmed  #--batch_system "local" condor