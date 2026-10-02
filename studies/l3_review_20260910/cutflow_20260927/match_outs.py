"""Match each FSR-fix fpo=1 job (completion marker = summary json in the staging dir) to the
condor .out in HiggsDNA/eos_logs that belongs to THAT run. eos_logs job dirs are reused across
productions (fpo=4 attempt 12-14 Sep, fpo=1 run from 15 Sep), so the collector's own pick
(largest cluster id with a CutFlow) is not a provenance guarantee. Fingerprint used here:
  CutFlow(zgammas)['all'] == summary n_events  AND  ['all cuts'] == summary n_events_selected['nominal'].
Output: one TSV per stage + a symlink mirror <mirror>/<stage>/<sample>/job_N/<cluster>.out."""
import importlib.util, json, glob, os, sys
from concurrent.futures import ProcessPoolExecutor
spec=importlib.util.spec_from_file_location("cc","/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/5_collect_cutflow.py")
cc=importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
LOGS="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs"
HERE=os.path.dirname(os.path.abspath(__file__))
MIRROR=os.path.join(HERE,"logs_mirror")

def zg(o):
    try: t=open(o,errors='ignore').read()
    except Exception: return None
    last=None
    for p in cc.iter_cutflow_payloads_from_text(t, syst="nominal"):
        if p.get('cut_type')=='zgammas': last=p['cuts']
    if last is None: return None
    return float(last.get('all',-1)), float(last.get('all cuts', last.get('allcuts',-1)))

def one(args):
    stage,smp,jd=args
    job=os.path.basename(jd)
    sj=glob.glob(f"{jd}/*_summary_job*.json")
    if not sj: return (stage,smp,job,'NOSUMMARY',None,0,0)
    j=json.load(open(sj[0])); ne=float(j['n_events']); ns=float(j['n_events_selected'].get('nominal',0))
    outs=glob.glob(f"{LOGS}/{stage}/{smp}/{job}/*.out")
    match=[]
    # Data: CutFlow 'all' is counted after the golden-JSON lumi mask, summary n_events before it,
    # so for Data the fingerprint is: same job end time (|mtime(.out)-mtime(summary)|<30 min),
    # 'all' <= n_events ('all cuts' may differ from the parquet row count by ~1 in a few jobs).
    smt=os.path.getmtime(sj[0])
    for o in outs:
        r=zg(o)
        if not r: continue
        if (r[0]==ne and r[1]==ns) or (stage=='Data' and r[0]<=ne and abs(os.path.getmtime(o)-smt)<1800): match.append(o)
    if not match: return (stage,smp,job,'NOMATCH' if outs else 'NOOUT',None,ne,ns)
    best=cc.pick_best_out_file(match)
    return (stage,smp,job,'OK',best,ne,ns)

if __name__=="__main__":
    base,stage=sys.argv[1],sys.argv[2]
    tasks=[]
    for sd in sorted(glob.glob(f"{base}/{stage}/*_20*")):
        smp=os.path.basename(sd)
        for jd in glob.glob(f"{sd}/job_*"):
            tasks.append((stage,smp,jd))
    with ProcessPoolExecutor(6) as ex:
        res=list(ex.map(one,tasks,chunksize=20))
    with open(os.path.join(HERE,f"match_{stage}.tsv"),"w") as f:
        for r in res: f.write("\t".join(map(str,r))+"\n")
    from collections import defaultdict
    agg=defaultdict(lambda: defaultdict(float))
    for st,smp,job,status,best,ne,ns in res:
        a=agg[smp]; a['jobs']+=1; a[status]+=1; a['ev_all']+=ne; a['sel_all']+=ns
        if status=='OK':
            a['ev_ok']+=ne; a['sel_ok']+=ns
            d=os.path.join(MIRROR,st,smp,job); os.makedirs(d,exist_ok=True)
            l=os.path.join(d,os.path.basename(best))
            if not os.path.islink(l): os.symlink(best,l)
    print(f"{'sample':32s} jobs  OK  NOMATCH NOOUT NOSUM   ev_ok/ev_all   sel_ok/sel_all")
    for smp,a in agg.items():
        print(f"{smp:32s} {int(a['jobs']):5d} {int(a['OK']):5d} {int(a['NOMATCH']):5d} {int(a['NOOUT']):5d} {int(a['NOSUMMARY']):4d}  {a['ev_ok']:.0f}/{a['ev_all']:.0f} ({100*a['ev_ok']/max(a['ev_all'],1):.2f}%)  {a['sel_ok']:.0f}/{a['sel_all']:.0f}")
