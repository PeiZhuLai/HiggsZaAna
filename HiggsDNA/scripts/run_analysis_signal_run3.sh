outdir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/Sig_MC"

rm -fr Parquet/Sig_MC

python scripts/run_analysis.py --config "metadata/za_signal_run3.json" --log-level "DEBUG" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "local" --short #--short local condor 
