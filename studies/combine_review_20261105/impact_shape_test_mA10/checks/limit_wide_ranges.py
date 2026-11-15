#!/usr/bin/env python3
"""Expected limits with the production ranges of the background turn-on parameters and with wider ones.
Input : Combine/root_t2w_xsec/<m>_Datacard_leptons.root (read only). Output: limit_wide_ranges.txt here.
Widened: turn-on to 95-125 GeV, width to 0.1-50, smearing to 0.05-20 GeV, the windows the code already uses at
mA24-29. All limits with --run blind, MH frozen, and rAbsAcc = rRelAcc = 1e-5 to resolve small changes."""
import subprocess, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit'); ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
B='/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/root_t2w_xsec/%d_Datacard_leptons.root'
LIM='--run blind --setParameters MH=125.38 --freezeParameters MH --rMax 2 --cminDefaultMinimizerStrategy 0 --rAbsAcc 1e-5 --rRelAcc 1e-5'.split()
def ranges(m):
    f=ROOT.TFile.Open(B%m); w=f.Get('w'); o=[]
    for v in w.allVars():
        n=v.GetName()
        if not n.startswith(('Exp_','Pow_','Lau_','Bern_')) or v.isConstant(): continue
        if 'turnon' in n or n=='Bern_step': o.append('%s=95,125'%n)
        elif 'width' in n or n=='Bern_stepWidth': o.append('%s=0.1,50'%n)
        elif 'sigma' in n: o.append('%s=0.05,20'%n)
    f.Close(); return ':'.join(o)
def lim(m,tag,extra):
    subprocess.run(['combine','-M','AsymptoticLimits',B%m,'-m','125.38','-n','.%s%d'%(tag,m)]+LIM+extra,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    f=ROOT.TFile.Open('higgsCombine.%s%d.AsymptoticLimits.mH125.38.root'%(tag,m)); x=[e.limit for e in f.Get('limit')][:5]; f.Close(); return x
out=['mA  q2.5 q16 q50 q84 q97.5 (production) | (wider) | q50 ratio wider/production']
for m in map(int,sys.argv[1:]):
    a=lim(m,'p',[]); b=lim(m,'w',['--setParameterRanges',ranges(m)])
    out.append('%-3d %s | %s | %.3f'%(m,' '.join('%.5f'%v for v in a),' '.join('%.5f'%v for v in b),b[2]/a[2])); print(out[-1],flush=True)
open('limit_wide_ranges.txt','w').write('\n'.join(out)+'\n')
