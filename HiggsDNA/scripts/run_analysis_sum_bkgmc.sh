outdir="/eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC"

# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/DYGto2LG_10to100_2023preBPix/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/DYJetsToLL_2023preBPix/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_sumStudy_DNA/Bkg_MC/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_sum_bkgmc_run3.json" --log-level "INFO" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short #--with_skimmed  #--batch_system "local" condor