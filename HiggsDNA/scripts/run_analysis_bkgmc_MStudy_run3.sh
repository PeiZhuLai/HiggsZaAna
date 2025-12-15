outdir="/eos/home-p/pelai/HZa/parquet_MStudy_DNA/Bkg_MC"

# rm -fr /eos/home-p/pelai/HZa/parquet_MStudy_DNA/Bkg_MC

python scripts/run_analysis.py --config "metadata/za_bkgmc_MStudy_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" --with_skimmed #--short #--batch_system "local" condor