outdir="/eos/home-p/pelai/HZa/parquet_tnp_zmmg/data"

# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/data/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/data/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_data_tnp_zmmg_run3.json" --log-level "INFO" --n_cores 10 --output_dir $outdir --unretire_jobs  --batch_system "condor" #--short #--batch_system "local" "condor" INFO DEBUG #--with_skimmed