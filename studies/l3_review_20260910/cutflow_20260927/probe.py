import importlib.util, json, glob, os, sys
spec=importlib.util.spec_from_file_location("cc","/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/5_collect_cutflow.py")
cc=importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
stage,smp,job=sys.argv[1:4]
S="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1" if stage!="Sig_MC" else "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix"
sj=glob.glob(f"{S}/{stage}/{smp}/{job}/*_summary_job*.json")
for f in sj:
    j=json.load(open(f)); print("summary", j['n_events'], j['sum_weights'], j['n_events_selected'].get('nominal'))
for o in sorted(glob.glob(f"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/{stage}/{smp}/{job}/*.out")):
    t=open(o,errors='ignore').read()
    for p in cc.iter_cutflow_payloads_from_text(t, syst="nominal"):
        if p.get('cut_type') in ('zgammas','zgammas_w'):
            print(os.path.basename(o), p['cut_type'], {k:p['cuts'][k] for k in list(p['cuts'])[:2]}, p['cuts'].get('all cuts'))
    print(os.path.basename(o), "genw", cc.parse_generator_weight_sum(t))
