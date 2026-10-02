import sys; sys.argv=['x']
sys.path.insert(0,'/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts')
import numpy as np, json
import plot_mA1_doublepeak as P
cut = {e["mA"]: e["MVAcut"] for e in json.load(open(P.WP_JSON))["results"]}[1]
a=P.load(cut); w=a["factor"]; g=a["GenALP_dR_gg"]
reco=a["ALP_lead_photon_pt"]+a["ALP_sublead_photon_pt"]; gen=a["GenALPLeadPho_pt"]+a["GenALPSubleadPho_pt"]
resp=np.where(gen>0,reco/np.where(gen>0,gen,1),np.nan)
d=P.dr
m1=(np.maximum(d(a["ALP_lead_photon_eta"],a["ALP_lead_photon_phi"],a["GenALPLeadPho_eta"],a["GenALPLeadPho_phi"]),d(a["ALP_sublead_photon_eta"],a["ALP_sublead_photon_phi"],a["GenALPSubleadPho_eta"],a["GenALPSubleadPho_phi"]))<0.02)|(np.maximum(d(a["ALP_lead_photon_eta"],a["ALP_lead_photon_phi"],a["GenALPSubleadPho_eta"],a["GenALPSubleadPho_phi"]),d(a["ALP_sublead_photon_eta"],a["ALP_sublead_photon_phi"],a["GenALPLeadPho_eta"],a["GenALPLeadPho_phi"]))<0.02)
ok=np.isfinite(resp)&(g>0); H=a["H_mass"]; M=a["ALP_mass"]
for lab,m in [("core",ok&(H>120)&(H<130)),("tail",ok&(H>134)&(H<137)),("mainpk",ok&(M<1.15)),("2ndpk",ok&(M>=1.15))]:
    print(lab,"dR<0.05 %.4f matched %.4f medresp %.4f"%(w[m&(g<0.05)].sum()/w[m].sum(), w[m&m1].sum()/w[m].sum(), P.wmedian(resp[m],w[m])))
c=ok&(g<0.05); 
for q in (0.16,0.5,0.84):
    o=np.argsort(M[c]); cw=np.cumsum(w[c][o]); print("m_gg quantile",q, M[c][o][np.searchsorted(cw,q*cw[-1])])
for q in (0.25,0.5,0.75):
    o=np.argsort(resp[c]); cw=np.cumsum(w[c][o]); print("resp quantile dR<0.05",q, resp[c][o][np.searchsorted(cw,q*cw[-1])])
print("raw ok", ok.sum(), "w ok", w[ok].sum(), "w all", w.sum())
