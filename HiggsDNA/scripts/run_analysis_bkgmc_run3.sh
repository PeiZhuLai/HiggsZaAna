outdir="/eos/home-p/pelai/HZa/parquet_cutflow_DNA/Bkg_MC"
# outdir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/parquet_DNA/Bkg_MC"

# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC/DYGto2LG_10to100_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC/DYJetsTo2E_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC/DYJetsTo2Mu_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC/DYJetsTo2Tau_2024/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_DNA/Bkg_MC/analysis_manager.pkl

python scripts/run_analysis.py --config "metadata/za_bkgmc_run3.json" --log-level "INFO" --n_cores 10 --output_dir $outdir --unretire_jobs --batch_system "condor" #--short #--with_skimmed  #--batch_system "local" condor