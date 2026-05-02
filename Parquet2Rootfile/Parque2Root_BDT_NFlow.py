#!/usr/bin/env python
import os
import math
from argparse import ArgumentParser
import numpy as np
#import time
import pandas as pd
import uproot
# from root_pandas import *
from tqdm import tqdm
import warnings
from ROOT import Math, TVector2, TVector3, TLorentzVector

from xgboost import XGBClassifier
import pickle
import re

warnings.simplefilter(action='ignore', category=FutureWarning)

# 內嵌模型檔路徑與 lazy cache
MODEL_FILE = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_run3_NFlow.pkl"
_MODEL_CACHE = None
def get_model():
    global _MODEL_CACHE
    if _MODEL_CACHE is None:
        with open(MODEL_FILE, 'rb') as f:
            _MODEL_CACHE = pickle.load(f)
    return _MODEL_CACHE

def getArgs():
    parser = ArgumentParser(description="Skim the input ntuples for Hmumu XGBoost analysis.")
    parser.add_argument('-i', '--input', action='store', default='inputs', help='Path to the input ntuple')
    parser.add_argument('-o', '--output', action='store', default='outputs', help='Path to the output ntuple')
    parser.add_argument('-s', '--split', action='store_true', help='Split train and test branch')
    # 修正: 使用正確型別與安全的預設值
    parser.add_argument('--chunksize', type=int, default=500000000, help='size to process at a time') 
    return  parser.parse_args()

def true_delta_phi(delta_phi):
    if delta_phi > math.pi:
        return 2 * math.pi - delta_phi
    return delta_phi

def compute_is_center(x):

    if (x.H_mass >= 115 and x.H_mass <= 135): return 1
    else: return 0

# others

def compute_Z_cosTheta(x):
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    M, mll = x.H_mass, x.Z_mass
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)

    return cosTheta

def compute_l_costheta(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    # H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    Z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    # H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    # Z_BH = Math.VectorUtil.boost(Z, -H_beta)
    l1_BZ = Math.VectorUtil.boost(l1, -Z_beta)
    # l2_BZ = Math.VectorUtil.boost(l2, -Z_beta)
    gamma_BZ = Math.VectorUtil.boost(gamma, -Z_beta)
    
    ## Z and lepton
    #a = l1_BZ + l2_BZ
    #cosTheta = (a.Vect().Unit()).Dot(l1_BZ.Vect().Unit())
    
    ## formula
    #a = l1_BZ.E() - l2_BZ.E()
    #k1 = math.sqrt(Z_BH.Vect().Dot(Z_BH.Vect()))
    #cosTheta = a/k1

    ## photon and lepton
    costheta = - (gamma_BZ.Vect()).Dot(l1_BZ.Vect())/gamma_BZ.Vect().R()/l1_BZ.Vect().R()

    return costheta

def compute_l_phi(x):
    if (x.Z_lead_lepton_charge < 0 and x.Z_sublead_lepton_charge > 0): 
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    elif (x.Z_sublead_lepton_charge < 0 and x.Z_lead_lepton_charge > 0):
        l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
        l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    else: print('leptons have same sign')
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    
    M, mll = x.H_mass, x.Z_mass
    lZ = math.sqrt((H.Dot(Z) / M) ** 2 - mll ** 2)

    H_transverse_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    H_transverse_beta.SetZ(0)
    hH = Math.VectorUtil.boost(H, -H_transverse_beta)
    hPz, hE = hH.Pz(), hH.E()
    q = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz + hE) / 2, (hE + hPz) / 2)
    q = Math.VectorUtil.boost(q, H_transverse_beta)
    qbar = Math.LorentzVector("ROOT::Math::PxPyPzE4D<float>")(0, 0, (hPz - hE) / 2, (hE - hPz) / 2)
    qbar = Math.VectorUtil.boost(qbar, H_transverse_beta)

    cosTheta = (qbar - q).Dot(Z)/(M * lZ)
    if (abs(cosTheta) > 1): 
        print("cosTheta = ", cosTheta)
        cosTheta = 1./cosTheta
    sinTheta = math.sqrt(1 - cosTheta ** 2)
    
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    l1_BH = Math.VectorUtil.boost(l1, -H_beta)
    l2_BH = Math.VectorUtil.boost(l2, -H_beta)
    Z_BH = Math.VectorUtil.boost(Z, -H_beta)
    N1_BH = l1_BH.Vect().Cross(l2_BH.Vect())
    
    q = Math.VectorUtil.boost(q, -H_beta)
    
    # beamAxis = TVector3 (0, 0, 1)
    # Z3_BH = Z_BH.Vect()
    # NSC_BH = (Z3_BH.Cross(beamAxis)).Unit()
    # tmpSgnPhi1 = Z3_BH.Dot(N1_BH.Cross(NSC_BH))
    
    NSC_BH = (q.Vect().Cross(Z_BH.Vect()))
    # tmpSgnPhi1 = - Z3_BH.Dot(N1_BH.Cross(NSC_BH))
    tmpSgnPhi1 = - N1_BH.Dot(q.Vect())/q.Vect().R()/N1_BH.R()/sinTheta
    
    sgnPhi1 = 0.
    if (abs(tmpSgnPhi1)>0.): sgnPhi1 = tmpSgnPhi1/abs(tmpSgnPhi1)
    dot_BH1SC = - N1_BH.Dot(NSC_BH)/NSC_BH.R()/N1_BH.R()
    if (abs(dot_BH1SC)>=1.): 
      print("dot_BH1SC = ", dot_BH1SC)
      dot_BH1SC *= 1./abs(dot_BH1SC)
    Phi1 = sgnPhi1 * math.acos(dot_BH1SC)

    return Phi1

def compute_dR1lg(x):

    if (((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5) > ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5:
        return ((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5
    else:
        return ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5
    
def compute_dR2lg(x):

    if (((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5) < ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5:
        return ((x.Z_lead_lepton_eta - x.gamma_eta)**2 + x.Z_lead_lepton_deltaphi**2)**0.5
    else:
        return ((x.Z_sublead_lepton_eta - x.gamma_eta)**2 + x.Z_sublead_lepton_deltaphi**2)**0.5

# mine

def compute_l_prodAngle(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    if (x.Z_lead_lepton_charge > 0):
        a = Math.VectorUtil.boost(l1, -z_beta).Vect()
    else :
        a = Math.VectorUtil.boost(l2, -z_beta).Vect()
    return a.Unit().Dot(z_beta.Unit())

def compute_Z_prodAngle(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    Z_BH = Math.VectorUtil.boost(Z, -H_beta).Vect()
    gamma_BH = Math.VectorUtil.boost(gamma, -H_beta).Vect()
    return Z_BH.Unit().Dot(H.Vect().Unit())

def compute_gamma_relEerror(x):

    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return x.gamma_energyErr / gamma.E()

def compute_G_ECM(x):

    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(gamma, -H_beta).E()

def compute_Z_ECM(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    return Math.VectorUtil.boost(Z, -H_beta).E()

def compute_l_rapCM(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    z_beta = TLorentzVector(Z.Px(), Z.Py(), Z.Pz(), Z.E()).BoostVector()
    a = Math.VectorUtil.boost(l1, -z_beta)
    return a.Rapidity()

def compute_Z_rapCM(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    H_beta = TLorentzVector(H.Px(), H.Py(), H.Pz(), H.E()).BoostVector()
    a = Math.VectorUtil.boost(Z, -H_beta)
    return a.Rapidity()

def compute_HZ_deltaRap(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    H = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.H_pt, x.H_eta, x.H_phi, x.H_mass)
    return H.Rapidity() - Z.Rapidity()

def compute_ll_deltaR(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    return Math.VectorUtil.DeltaR(l1, l2)

def compute_leadLG_deltaR(x):

    l1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_lead_lepton_pt, x.Z_lead_lepton_eta, x.Z_lead_lepton_phi, x.Z_lead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l1, gamma)

def compute_subleadLG_deltaR(x):

    l2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_sublead_lepton_pt, x.Z_sublead_lepton_eta, x.Z_sublead_lepton_phi, x.Z_sublead_lepton_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    return Math.VectorUtil.DeltaR(l2, gamma)

def compute_ZG_deltaR(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return Math.VectorUtil.DeltaR(Z, gamma)

def compute_H_ptt(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / (Z - gamma).Pt() * 2.0


def compute_H_al(x):
    
    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

    return (Z.Pt() ** 2 - gamma.Pt() ** 2) / (Z - gamma).Pt()

def compute_H_bt(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    r = abs(Z.Pt()) / abs(gamma.Pt())
    dev = ((Z.Px() - r * gamma.Px()) ** 2 + (Z.Py() - r * gamma.Py()) ** 2) ** 0.5

    return r * abs(Z.Px() * gamma.Py() - gamma.Px() * Z.Py()) / dev * 2.0

def compute_jet_ptt(x):

    if x.n_jets < 2:
        return -9999
    else:
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)

        return abs(Jets_1.Px() * Jets_2.Py() - Jets_2.Px() * Jets_1.Py()) / (Jets_1 - Jets_2).Pt() * 2.0

def compute_max_two_jet_btag(x):

    n_jets = x.n_jets
    jet_btag_array = []

    for i in range(1, min(n_jets+1, 5)):
        jet_name = f"jet_{i}"
        jet_btag_array.append(getattr(x, f"{jet_name}_btagDeepFlavB"))

    jet_btag_array = np.array(jet_btag_array)
    jet_btag_array_sorted = np.sort(jet_btag_array)[::-1]

    return np.sum(jet_btag_array_sorted[:2])


def compute_mass_jj(x):

    if x.n_jets < 2:
        return -9999
    else:
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return j1j2.M()

def compute_delta_eta_jj(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.jet_1_eta-x.jet_2_eta)

def compute_delta_phi_jj(x):

    if x.n_jets >= 2:
        return true_delta_phi(abs(x.jet_1_phi-x.jet_2_phi))
    return -9999

def compute_delta_phi_zg_jj(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return true_delta_phi(abs(zg.Phi()-j1j2.Phi()))
    return -9999

def compute_delta_eta_zg_jj(x):

    if x.n_jets >= 2:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        zg = Z+gamma
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jets_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        j1j2 = Jets_1+Jets_2
        return abs(zg.Eta()-j1j2.Eta())
    return -9999

def compute_photon_zeppenfeld(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.gamma_eta-(x.jet_1_eta+x.jet_2_eta)/2)

def compute_H_zeppenfeld(x):

    if x.n_jets < 2:
        return -9999
    else:
        return abs(x.H_eta-(x.jet_1_eta+x.jet_2_eta)/2)

def compute_pt_balance(x):

    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()/(x.Z_pt+x.gamma_pt+x.jet_1_pt+x.jet_2_pt)

def compute_system_pt(x):
    
    if x.n_jets < 2:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        total = Z+gamma+Jet_1+Jet_2
        return total.Pt()

def compute_jet_pair_pt(x):
    
    if x.n_jets < 2:
        return -9999
    else:
        Jet_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        Jet_2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_2_pt, x.jet_2_eta, x.jet_2_phi, x.jet_2_mass)
        Jet_pair = Jet_1+Jet_2
        return Jet_pair.Pt()

def compute_pt_balance_0j(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
    total = Z+gamma
    return total.Pt()/(x.Z_pt+x.gamma_pt)

def compute_pt_balance_1j(x):

    if x.n_jets < 1:
        return -9999
    else:
        Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)
        Jets_1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.jet_1_pt, x.jet_1_eta, x.jet_1_phi, x.jet_1_mass)
        total = Z+gamma+Jets_1
        return total.Pt()/(x.Z_pt+x.gamma_pt+x.jet_1_pt)

def compute_Delta_Phi(x, var = "gamma_phi", min_jet=0):

    if min_jet:
        if x.n_jets < min_jet: return -9999
        return true_delta_phi(abs(getattr(x, "jet_%d_phi" % min_jet) - x.gamma_phi))
    else:
        return true_delta_phi(abs(getattr(x, var) - x.gamma_phi))

def compute_Delta_R(x, min_jet=0):

    if x.n_jets < min_jet: return -9999
    else: 
        jet = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(getattr(x, "jet_%d_pt" % min_jet), getattr(x, "jet_%d_eta" % min_jet), getattr(x, "jet_%d_phi" % min_jet), getattr(x, "jet_%d_mass" % min_jet))
        gamma = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.gamma_pt, x.gamma_eta, x.gamma_phi, x.gamma_mass)

        return Math.VectorUtil.DeltaR(jet, gamma)

def compute_QG(x):

    if x.Jets_jetMultip >= 1 and (abs(x.Jets_Eta_Lead) > 2.1 or x.Jets_PT_Lead < 50):
        Jets_QGscore_Lead, Jets_QGflag_Lead = -1, -1
    else:
        Jets_QGscore_Lead = x.Jets_NTracks_Lead
        Jets_QGflag_Lead = np.heaviside(x.Jets_NTracks_Lead - 11, 0) + np.heaviside(-x.Jets_NTracks_Lead - 9999, -9999)

    if x.Jets_jetMultip >= 2 and (abs(x.Jets_Eta_Sub) > 2.1 or x.Jets_PT_Sub < 50):
        Jets_QGscore_Sub, Jets_QGflag_Sub = -1, -1
    else:
        Jets_QGscore_Sub = x.Jets_NTracks_Sub
        Jets_QGflag_Sub = np.heaviside(x.Jets_NTracks_Sub - 11, 0) + np.heaviside(-x.Jets_NTracks_Sub - 9999, -9999)

    return Jets_QGscore_Lead, Jets_QGflag_Lead, Jets_QGscore_Sub, Jets_QGflag_Sub

def compute_dR_Z_ALP(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    axion = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.ALP_pt, x.ALP_eta, x.ALP_phi, x.ALP_mass)

    return Math.VectorUtil.DeltaR(Z, axion)

def compute_dR_g1_g2(x):

    g1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.ALP_lead_photon_pt, x.ALP_lead_photon_eta, x.ALP_lead_photon_phi, x.ALP_lead_photon_mass)
    g2 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.ALP_sublead_photon_pt, x.ALP_sublead_photon_eta, x.ALP_sublead_photon_phi, x.ALP_sublead_photon_mass)

    return Math.VectorUtil.DeltaR(g1, g2)

def compute_dR_Z_g1(x):

    Z = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.Z_pt, x.Z_eta, x.Z_phi, x.Z_mass)
    g1 = Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<float>")(x.ALP_lead_photon_pt, x.ALP_lead_photon_eta, x.ALP_lead_photon_phi, x.ALP_lead_photon_mass)

    return Math.VectorUtil.DeltaR(Z, g1)

def compute_delta_eta_g1Z(x):

    return abs(x.ALP_lead_photon_eta-x.Z_eta)

def compute_delta_phi_g1Z(x):

    return true_delta_phi(abs(x.ALP_lead_photon_phi-x.Z_phi))

def use_corrected_alp_photon_columns(data):
    corrected_columns = {
        "ALP_sublead_photon_ecalPFClusterIso": "ALP_sublead_photon_ecalPFClusterIso_corr",
        "ALP_sublead_photon_sieie": "ALP_sublead_photon_sieie_corr",
        "ALP_lead_photon_ecalPFClusterIso": "ALP_lead_photon_ecalPFClusterIso_corr",
        "ALP_lead_photon_sieie": "ALP_lead_photon_sieie_corr",
    }

    replaced_columns = []
    missing_columns = []
    for target, source in corrected_columns.items():
        if source in data.columns:
            data[target] = data[source]
            replaced_columns.append(source)
        else:
            missing_columns.append(source)

    if replaced_columns:
        print(
            "Replaced ALP photon isolation/sieie branches with corrected versions: %s"
            % ", ".join(replaced_columns)
        )
    if missing_columns:
        print(
            "Corrected ALP photon branches not found; keeping nominal branches for: %s"
            % ", ".join(missing_columns)
        )

    return data

def ensure_za_compatibility(data):
    """
    Some ZA parquet configs only write ALP photon branches, while this converter
    still decorates legacy single-gamma variables. Preserve existing gamma_*
    branches when present and otherwise alias them from the leading ALP photon.
    """
    fallback_columns = {
        "gamma_pt": "ALP_lead_photon_pt",
        "gamma_eta": "ALP_lead_photon_eta",
        "gamma_phi": "ALP_lead_photon_phi",
        "gamma_mass": "ALP_lead_photon_mass",
        "gamma_mvaID": "ALP_lead_photon_mvaID",
        "gamma_energyErr": "ALP_lead_photon_energyErr",
        "gamma_sieie": "ALP_lead_photon_sieie",
        "gamma_hoe": "ALP_lead_photon_hoe",
        "gamma_r9": "ALP_lead_photon_r9",
        "gamma_chiso": "ALP_lead_photon_chiso",
        "gamma_alliso": "ALP_lead_photon_alliso",
        "gamma_e_veto": "ALP_lead_photon_electronVeto",
    }

    added_columns = []
    for target, source in fallback_columns.items():
        if target not in data.columns and source in data.columns:
            data[target] = data[source]
            added_columns.append(target)

    if added_columns:
        print("Aliased missing gamma branches from ALP_lead_photon: %s" % ", ".join(added_columns))

    required_columns = ["gamma_pt", "gamma_eta", "gamma_phi", "gamma_mass"]
    missing_required = [column for column in required_columns if column not in data.columns]
    if missing_required:
        raise KeyError(
            "Missing required photon branches after compatibility mapping: %s. "
            "Available columns include: %s"
            % (", ".join(missing_required), ", ".join(list(data.columns)[:50]))
        )

    return data

MVA_FEATURE_COLUMNS = [
    "pho1Pt",
    "pho1R9",
    "pho1IetaIeta55",
    "pho1PIso_noCorr",
    "pho2Pt",
    "pho2R9",
    "pho2IetaIeta55",
    "pho2PIso_noCorr",
    "ALP_calculatedPhotonIso",
    "var_dR_Za",
    "var_dR_g1g2",
    "var_dR_g1Z",
    "var_PtaOverMh",
    "H_pt",
    "param",
]

def compute_MVA_Score(x, wanted_ma):
    """
    與 ALP_plot_param.py 一致的輸入變數順序，回傳 predict_proba 的正類機率。
    """
    # 取得模型（快取）
    try:
        model = get_model()
    except Exception as e:
        # 載入失敗時回傳 NaN
        return np.nan

    # 防呆：H_mass 為 0 或非有限
    denom = x.H_mass if x.H_mass != 0 else np.nan
    if not np.isfinite(denom) or denom == 0:
        return np.nan

    # param 定義
    param = (x.ALP_mass - wanted_ma) / denom

    # 特徵順序（v_1）
    MVA_list = [
        x.pho1Pt,
        x.pho1R9,
        x.pho1IetaIeta55,
        x.pho1PIso_noCorr,
        x.pho2Pt,
        x.pho2R9,
        x.pho2IetaIeta55,
        x.pho2PIso_noCorr,
        x.ALP_calculatedPhotonIso,
        x.var_dR_Za,
        x.var_dR_g1g2,
        x.var_dR_g1Z,
        x.var_PtaOverMh,
        x.H_pt,
        param
    ]

    try:
        proba = model.predict_proba([MVA_list])[:, 1]
        return float(proba[0])
    except Exception:
        return np.nan

def add_mva_scores(data, masses=range(1, 31)):
    """
    Add MVA_Score_mA_M* columns using batched XGBoost inference.

    The old implementation called predict_proba once per row per mass. For large
    parquet files that is effectively dominated by Python overhead. Here each
    mass is evaluated with one dense feature matrix.
    """
    try:
        model = get_model()
    except Exception as e:
        print("WARNING: failed to load MVA model: %s" % e)
        for ma in masses:
            data[f"MVA_Score_mA_M{ma}"] = np.nan
        return data

    h_mass = pd.to_numeric(data["H_mass"], errors="coerce").to_numpy(dtype=float)
    alp_mass = pd.to_numeric(data["ALP_mass"], errors="coerce").to_numpy(dtype=float)
    valid_base = np.isfinite(h_mass) & (h_mass != 0)
    base_features = data[MVA_FEATURE_COLUMNS[:-1]].apply(pd.to_numeric, errors="coerce")

    score_columns = {}
    for ma in masses:
        scores = np.full(data.shape[0], np.nan, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            param = (alp_mass - ma) / h_mass

        features = base_features.copy()
        features["param"] = param
        matrix = features[MVA_FEATURE_COLUMNS].to_numpy(dtype=float)
        valid = valid_base & np.isfinite(matrix).all(axis=1)

        if np.any(valid):
            try:
                scores[valid] = model.predict_proba(matrix[valid])[:, 1]
            except Exception as e:
                print("WARNING: failed to evaluate MVA score for mA_M%s: %s" % (ma, e))

        score_columns[f"MVA_Score_mA_M{ma}"] = scores

    return pd.concat([data, pd.DataFrame(score_columns, index=data.index)], axis=1)

def preselect(data):

    #data.query('(Muons_Minv_MuMu_Paper >= 110) | (Event_Paper_Category >= 17)', inplace=True)
    #data.query('Event_Paper_Category > 0', inplace=True)

    return data

def decorate(data):

    if data.shape[0] == 0: return data

    # Za training Variables
    data['pho1Pt'] = data.ALP_lead_photon_pt
    data['pho1R9'] = data.ALP_lead_photon_r9
    data['pho1IetaIeta55'] = data.ALP_lead_photon_sieie
    data['pho1PIso_noCorr'] = data.ALP_lead_photon_ecalPFClusterIso
    data['pho2Pt'] = data.ALP_sublead_photon_pt
    data['pho2R9'] = data.ALP_sublead_photon_r9
    data['pho2IetaIeta55'] = data.ALP_sublead_photon_sieie
    data['pho2PIso_noCorr'] = data.ALP_sublead_photon_ecalPFClusterIso
    data['ALP_calculatedPhotonIso'] = data.ALP_PhotonIso
    data['var_PtaOverMh'] = data.ALP_pt / data.H_mass
    data['var_dR_Za'] = data.apply(lambda x: compute_dR_Z_ALP(x), axis=1) 
    data['var_dR_g1g2'] = data.apply(lambda x: compute_dR_g1_g2(x), axis=1) 
    data['var_dR_g1Z'] = data.apply(lambda x: compute_dR_Z_g1(x), axis=1) 
    
    data = add_mva_scores(data)

    data['H_m'] = data.H_mass
    data['ALP_m'] = data.ALP_mass
    data['Z_m'] = data.Z_mass

    data['var_dEta_g1Z'] = data.apply(lambda x: compute_delta_eta_g1Z(x), axis=1) 
    data['var_dPhi_g1Z'] = data.apply(lambda x: compute_delta_phi_g1Z(x), axis=1) 
    data['var_PtaOverMa'] = data.ALP_pt / data.ALP_mass
    data['var_Pta'] = data.ALP_pt
    data['var_MhMa'] = data.H_mass + data.ALP_mass
    data['var_MhMZ'] = data.H_mass + data.Z_mass

    # data['weight'] = data.weight_central
    data['weight'] = data.weight_central
    data['factor'] = data.weight_central
    data['is_center'] = data.apply(lambda x: compute_is_center(x), axis=1)

    # data['H_ptt'] = data.apply(lambda x: compute_H_ptt(x), axis=1)
    # data['H_al'] = data.apply(lambda x: compute_H_al(x), axis=1)
    # data['H_bt'] = data.apply(lambda x: compute_H_bt(x), axis=1)
    # data['Z_cos_theta'] = data.apply(lambda x:compute_Z_cosTheta(x), axis=1)
    # data['lep_cos_theta'] = data.apply(lambda x: compute_l_costheta(x), axis=1)
    # data['lep_phi'] = data.apply(lambda x: compute_l_phi(x), axis=1)
    # data['l1g_deltaR'] = data.apply(lambda x: compute_dR1lg(x), axis=1) 
    # data['l2g_deltaR'] = data.apply(lambda x: compute_dR2lg(x), axis=1)

    data['HZ_relM'] = data.H_mass / data.Z_mass
    data['H_relpt'] = data.H_pt / data.H_mass
    data['Z_relpt'] = data.Z_pt / data.H_mass
    data['Z_lead_lepton_relpt'] = data.Z_lead_lepton_pt / data.H_mass
    data['Z_sublead_lepton_relpt'] = data.Z_sublead_lepton_pt / data.H_mass
    data['gamma_relpt'] = data.gamma_pt / data.H_mass
    # data['gamma_ptRelErr'] = data.apply(lambda x:compute_gamma_relEerror(x), axis=1)
    data['G_ECM'] = data.apply(lambda x:compute_G_ECM(x), axis=1)
    data['Z_ECM'] = data.apply(lambda x:compute_Z_ECM(x), axis=1)
    data['Z_rapCM'] = data.apply(lambda x:compute_Z_rapCM(x), axis=1)
    data['l_rapCM'] = data.apply(lambda x:compute_l_rapCM(x), axis=1)
    data['HZ_deltaRap'] = data.apply(lambda x:compute_HZ_deltaRap(x), axis=1)
    data['l_cosProdAngle'] = data.apply(lambda x:compute_l_prodAngle(x), axis=1)
    data['Z_cosProdAngle'] = data.apply(lambda x:compute_Z_prodAngle(x), axis=1)
    # delta R (object, gamma)
    data['ll_deltaR'] = data.apply(lambda x:compute_ll_deltaR(x), axis=1)
    data['leadLG_deltaR'] = data.apply(lambda x:compute_leadLG_deltaR(x), axis=1)
    data['ZG_deltaR'] = data.apply(lambda x:compute_ZG_deltaR(x), axis=1)
    data['subleadLG_deltaR'] = data.apply(lambda x:compute_subleadLG_deltaR(x), axis=1)
    # delta Phi (object, gamma)
    data['H_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'H_phi'), axis=1)
    data['Z_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_phi'), axis=1)
    data['Z_lead_lepton_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_lead_lepton_phi'), axis=1)
    data['Z_sublead_lepton_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'Z_sublead_lepton_phi'), axis=1)
    # Jet 
    # data['max_jet_deltaR'] = data[['jet1G_deltaR', 'jet2G_deltaR']].max(axis=1)
    # data['min_jet_deltaR'] = data[['jet1G_deltaR', 'jet2G_deltaR']].min(axis=1)
    # data['MET_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'MET_phi'), axis=1)
    # data['MET_relpt'] = data.MET_pt / data.H_mass
    # data['system_pt'] = data.apply(lambda x: compute_system_pt(x), axis=1)
    # data['jet_pair_pt'] = data.apply(lambda x: compute_jet_pair_pt(x), axis=1)
    # data['jet_1_relpt'] = data.jet_1_pt / data.H_mass
    # data['jet_2_relpt'] = data.jet_2_pt / data.H_mass
    # for i in np.arange(1,5):
    #     data['jet_%d_deltaphi' %i] = data.apply(lambda x: compute_Delta_Phi(x, "jet", min_jet=i), axis=1)
    #     data['jet%dG_deltaR' %i] = data.apply(lambda x: compute_Delta_R(x, min_jet=i), axis=1)
    # data['mass_jj'] = data.apply(lambda x: compute_mass_jj(x), axis=1)
    # data['max_two_jet_btag'] = data.apply(lambda x: compute_max_two_jet_btag(x), axis=1)
    # data['jet_ptt'] = data.apply(lambda x: compute_jet_ptt(x), axis=1)
    # data['delta_eta_jj'] = data.apply(lambda x: compute_delta_eta_jj(x), axis=1)
    # data['delta_phi_jj'] = data.apply(lambda x: compute_delta_phi_jj(x), axis=1)
    # data['delta_phi_zgjj'] = data.apply(lambda x: compute_delta_phi_zg_jj(x), axis=1)
    # data['delta_eta_zgjj'] = data.apply(lambda x: compute_delta_eta_zg_jj(x), axis=1)
    # data['photon_zeppenfeld'] = data.apply(lambda x: compute_photon_zeppenfeld(x), axis=1)
    # data['H_zeppenfeld'] = data.apply(lambda x: compute_H_zeppenfeld(x), axis=1)
    # data['pt_balance'] = data.apply(lambda x: compute_pt_balance(x), axis=1)
    # data['pt_balance_0j'] = data.apply(lambda x: compute_pt_balance_0j(x), axis=1)
    # data['pt_balance_1j'] = data.apply(lambda x: compute_pt_balance_1j(x), axis=1)
    # data[['Jets_QGscore_Lead', 'Jets_QGflag_Lead', 'Jets_QGscore_Sub', 'Jets_QGflag_Sub']] = data.apply(lambda x: compute_QG(x), axis=1, result_type='expand')
    # data.rename(columns={'Muons_Minv_MuMu_Paper': 'm_mumu', 'Muons_Minv_MuMu_VH': 'm_mumu_VH', 'EventInfo_EventNumber': 'eventNumber', 'Jets_jetMultip': 'n_j'}, inplace=True)
    # data.drop(['PassesttHSelection', 'PassesVHSelection', 'GlobalWeight', 'SampleOverlapWeight', 'EventWeight_MCCleaning_5'], axis=1, inplace=True)
    # Additioanl Lepton
    # data['additional_lepton_1_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'additional_lepton_1_phi', min_jet=0), axis=1)
    # data['additional_lepton_2_deltaphi'] = data.apply(lambda x: compute_Delta_Phi(x, 'additional_lepton_2_phi', min_jet=0), axis=1) 
    
    # --- robust dtype handling (avoid "setting an array element with a sequence") ---
    # Some parquet columns can be jagged/list-like => DataFrame-wide astype(float) will crash.
    def _is_sequence_cell(v):
        return isinstance(v, (list, tuple, np.ndarray))

    # Identify object columns that contain sequences (jagged)
    seq_cols = []
    obj_cols = data.select_dtypes(include=["object"]).columns
    for c in obj_cols:
        s = data[c]
        if s.shape[0] == 0:
            continue
        sample = s.dropna().head(50)
        if sample.map(_is_sequence_cell).any():
            seq_cols.append(c)

    # Convert all non-sequence columns to numeric where possible
    cols_to_numeric = [c for c in data.columns if c not in seq_cols]
    data[cols_to_numeric] = data[cols_to_numeric].apply(pd.to_numeric, errors="coerce")

    # Enforce integer-like columns (fill NaN to keep astype stable)
    int_cols = [
        "is_center",
        "Z_lead_lepton_charge", "Z_lead_lepton_id",
        "Z_sublead_lepton_charge", "Z_sublead_lepton_id",
        "n_leptons", "n_electrons", "n_muons",
        "event",
    ]
    present_int_cols = [c for c in int_cols if c in data.columns]
    if present_int_cols:
        data[present_int_cols] = data[present_int_cols].fillna(0).astype("int64")

    return data

def main():
    
    args = getArgs()

    variables = [
        'H_pt', 'H_eta', 'H_phi', 'H_mass',
        'Z_pt', 'Z_eta', 'Z_phi', 'Z_mass',
        'Z_lead_lepton_pt', 'Z_lead_lepton_eta', 'Z_lead_lepton_phi', 'Z_lead_lepton_mass',
        'Z_lead_lepton_charge', 'Z_lead_lepton_id',
        'Z_sublead_lepton_pt', 'Z_sublead_lepton_eta', 'Z_sublead_lepton_phi', 'Z_sublead_lepton_mass',
        'Z_sublead_lepton_charge', 'Z_sublead_lepton_id',
        'gamma_pt', 'gamma_eta', 'gamma_phi', 'gamma_mass',
        'gamma_mvaID', 'gamma_energyErr',
        "gamma_fsr_pt","gamma_fsr_eta","gamma_fsr_phi","gamma_fsr_mass","gamma_fsr_clean",
        'jet_1_pt', 'jet_1_eta', 'jet_1_phi', 'jet_1_mass', 'jet_1_btagDeepFlavB',
        'jet_2_pt', 'jet_2_eta', 'jet_2_phi', 'jet_2_mass', 'jet_2_btagDeepFlavB',
        'n_jets', 'n_leptons', 'n_electrons', 'n_muons', 'n_iso_photons',
        'MET_pt', 'MET_phi', 
        'weight_central',
        'event'
    ]

    print("Now!!! Porcessing {:s} ......".format(args.input))

    if os.path.isfile(args.output): os.remove(args.output)

    initial_events = 0
    final_events = 0

#    for data in tqdm(read_root(args.input, key='DiMuonNtuple', columns=variables, chunksize=args.chunksize), desc='Processing %s' % args.input, bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):

    data = pd.read_parquet(args.input)
    data = use_corrected_alp_photon_columns(data)
    data = ensure_za_compatibility(data)
    initial_events += data.shape[0]
    #data = preprocess(data)
    data = preselect(data) #TODO add cutflow
    data = decorate(data)
    final_events += data.shape[0]

    indices = np.arange(len(data))
    train_mask = indices % 2 == 0
    test_mask = indices % 2 == 1
    data_train = data[train_mask]
    data_test = data[test_mask]
    # data_zero_jet = data.query("n_jets == 0 & n_leptons == 2 & MET_pt < 90")
    # data_one_jet = data.query("n_jets == 1 & n_leptons == 2 & MET_pt < 90")
    # data_two_jet = data.query("n_jets >= 2 & n_leptons == 2 & n_b_jets == 0")
    # data_zero_to_one_jet = data.query("n_jets <= 1 & n_leptons == 2 & MET_pt < 90")
    # data_VH_ttH = data[(data.n_leptons > 2) & (data.n_b_jets == 0)]
    # data_VH =  data.query("n_leptons >= 3 & n_b_jets == 0 & max_I_mini < 0.15 & H_relpt > 0.3 & MET_pt > 30 & Z_mass > 85 & Z_mass < 95")
    # data_ZH = data.query("n_leptons == 2 & n_jets <= 1 & MET_pt > 90 & H_relpt > 0.4 & Z_mass > 85 & Z_mass < 95") 
    # data_ttH_had = data.query("n_leptons == 2 & n_jets >= 5 & n_b_jets >= 1 & Z_mass > 85 & Z_mass < 95")
    # data_ttH_lep = data.query("((n_leptons == 3 & n_jets >= 3 & n_b_jets >= 1) | (n_leptons >= 4 & n_jets >= 1 & n_b_jets >= 1)) & (max_I_mini < 0.1 & Z_mass > 85 & Z_mass < 95)")
    with uproot.recreate(args.output) as f:
        f['inclusive'] = data
        if args.split:
            f['train'] = data_train
            f['test'] = data_test
        # f['zero_jet'] = data_zero_jet
        # f['one_jet'] = data_one_jet
        # f['zero_to_one_jet'] = data_zero_to_one_jet
        # f['two_jet'] = data_two_jet
        # f['VH_ttH'] = data_VH_ttH
        # f['VH'] = data_VH
        # f['ZH'] = data_ZH
        # f['ttH_had'] = data_ttH_had
        # f['ttH_lep'] = data_ttH_lep
    # data.to_root(args.output, key='inclusive', mode='a', index=False)
    # data_zero_jet.to_root(args.output, key='zero_jet', mode='a', index=False)
    # data_one_jet.to_root(args.output, key='one_jet', mode='a', index=False)
    # data_zero_to_one_jet.to_root(args.output, key='zero_to_one_jet', mode='a', index=False)
    # data_two_jet.to_root(args.output, key='two_jet', mode='a', index=False)
    # data_VH_ttH.to_root(args.output, key='VH_ttH', mode='a', index=False)

    # meta_data = pd.DataFrame({'initial_events': [initial_events], 'final_events': [final_events]})
    # meta_data.to_root(args.output, key='MetaData', mode='a', index=False)

    print("Finished!!! Have gotten the skimmed data in {:s}".format(args.output))

if __name__ == '__main__':
    main()
