outdir="/eos/home-p/pelai/HZa/parquet_tnp_zmmg/mc"

# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/mc/job_1
rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/mc/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_mc_tnp_zmmg_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" --short #--with_skimmed  #--batch_system "local" condor
