"""FSR-fix signal: initial (summary n_events), final (summary n_events_selected[nominal]),
chunk rows (find, excluding merged_nominal.parquet), merged_nominal rows; old = Aug-11 cutflow_list JSON
(zgammas all / all cuts) and the AN table currently in merged_cutflow_latex (Jun-11)."""
import json,glob,os,subprocess,re
import pyarrow.parquet as pq
S="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
OLD="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list"
rows=[]
for d in sorted(glob.glob(S+'/mA_M*_20*')):
    smp=os.path.basename(d); m=re.match(r"mA_M(\d+)_(.*)",smp); ma=int(m.group(1)); era=m.group(2)
    ne=ns=0
    for f in glob.glob(d+'/job_*/*_summary_job*.json'):
        j=json.load(open(f)); ne+=j['n_events']; ns+=j['n_events_selected'].get('nominal',0)
    chunks=subprocess.run(["find",d,"-mindepth",2,"-name","*_nominal.parquet"] if False else ["find",d,"-mindepth","2","-name","output_job_*_nominal.parquet"],capture_output=True,text=True).stdout.split()
    chunks=[c for c in chunks if os.path.basename(c)!="merged_nominal.parquet"]
    cr=sum(pq.ParquetFile(c).metadata.num_rows for c in chunks)
    mf=d+'/merged_nominal.parquet'
    mr=pq.ParquetFile(mf).metadata.num_rows if os.path.exists(mf) else -1
    mt_ok = (not chunks) or os.path.getmtime(mf)>=max(os.path.getmtime(c) for c in chunks)
    oa=oc=None
    oj=f"{OLD}/cutflow_Sig_MC_mA_M{ma}_{era}.json"
    if os.path.exists(oj):
        z=json.load(open(oj))['cutflows']['zgammas']; oa=z['all']; oc=z.get('all cuts')
    rows.append((era,ma,ne,ns,cr,mr,mt_ok,oa,oc))
with open("sig_final_fsrfix_vs_old.tsv","w") as f:
    f.write("era\tmA\tnew_all\tnew_final\tchunk_rows\tmerged_rows\tmtime_ok\told_all\told_final\tnew_eff%\told_eff%\trel_diff%\n")
    for era,ma,ne,ns,cr,mr,mt,oa,oc in sorted(rows,key=lambda r:(r[0],r[1])):
        ne_=100*ns/ne; oe=100*oc/oa if oa else float('nan')
        f.write(f"{era}\t{ma}\t{ne}\t{ns}\t{cr}\t{mr}\t{mt}\t{oa}\t{oc}\t{ne_:.3f}\t{oe:.3f}\t{100*(ne_/oe-1):+.2f}\n")
print(open("sig_final_fsrfix_vs_old.tsv").read())
