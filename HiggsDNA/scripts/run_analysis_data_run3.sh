outdir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/Data"

python scripts/run_analysis.py --config "metadata/za_data_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs  --batch_system "condor" #--short #--batch_system "local" "condor"