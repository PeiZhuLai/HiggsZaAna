set +u; source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh; conda activate higgs-alp-ana
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/datavmc_opt_20260929
for r in SR CR mva; do
  ls /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc/ALP_plot_run3_UL_${r}_nominal_part_*.root > parts_$r.txt
  echo "$r parts: $(wc -l < parts_$r.txt)"
  /usr/bin/time -f "$r merge wall %e s, maxrss %M kB" python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts/merge_dataVmc_hists.py /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/datavmc_opt_20260929/merged_${r}_nominal.root $(cat parts_$r.txt)
  python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/datavmc_opt_20260929/compare_merged.py /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc/ALP_plot_run3_UL_${r}_nominal.root /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/datavmc_opt_20260929/merged_${r}_nominal.root; echo "$r compare rc=$?"
done
