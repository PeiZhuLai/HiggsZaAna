import json,glob,os,sys
base=sys.argv[1]
for d in sorted(glob.glob(base+'/*_20*')):
    s=sorted(glob.glob(d+'/job_*/*_summary_job*.json'))
    jobs=glob.glob(d+'/job_*')
    ne=0;nsel=0;sw=0;ok=0
    for f in s:
        j=json.load(open(f)); ne+=j['n_events']; sw+=j['sum_weights']; nsel+=j['n_events_selected'].get('nominal',0); ok+=bool(j.get('successful'))
    print(os.path.basename(d), len(jobs), len(s), ok, ne, round(sw,1), nsel)
