outdir="/eos/home-p/pelai/HZa/parquet_tnp_zmmg/data"
# outdir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/parquet_DNA/Data"

# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Data
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Data/Data_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Data/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_data_tnp_zmmg_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs  --batch_system "local" --short #--batch_system "local" "condor"