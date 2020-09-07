#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from array import array
import argparse
import sys
import ROOT
import json
import re
import os

from DataFormats.FWLite import Events, Handle
from EgammaUser.EgammaDAS2020.EvtData import EvtData, EvtHandles, add_product
import EgammaUser.EgammaDAS2020.CoreTools as CoreTools
import EgammaUser.EgammaDAS2020.GenTools as GenTools
import EgammaUser.EgammaDAS2020.MathTools as MathTools

""" 
This is an example of how to access E/gamma objects in the MiniAOD in fwlite. 
It uses a simple package developed to ease the usage of fwlite
Its actually for HLT analysis but works well in general

"""
def get_trig_indx(selected_name,trig_names):
    """
    returns the index of the trigger (if it exists, otherwise returns the size of trig_names, yes this is a bit c++ like)
    the annoying thing here is we want to match all HLT_Blah_v*
    there are several ways to do it, the quickest is see which starts with HLT_Blah_v
    
    note one can be smarter and cache the result and update when trig_names.parameterSetID() changes 
    """
    for idx,name in enumerate(trig_names.triggerNames()):
        if name.startswith(selected_name):
            return idx
    return len(trig_names)


def match_trig_objs(eta,phi,trig_objs,max_dr=0.1):    
    max_dr2 = max_dr*max_dr
    matched_objs = [obj for obj in trig_objs if MathTools.cal_delta_r2(eta,phi,obj.eta(),obj.phi()) < max_dr2]
    return matched_objs

if __name__ == "__main__":

    #first we load in the FWLite libaries and enable them 
    CoreTools.load_fwlitelibs()

    #now read the cmd line to find our input filename
    #local file "file:<filename>
    #remote file "root:// 
    #note unlike CMSSW proper its not smart enough to resolve to fall over to xrootd automatically
    #so it has to be specified
    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('in_filenames',nargs="+",help='input filenames')
    parser.add_argument('--prefix','-p',default='',help='file prefix')
    parser.add_argument('--out','-o',default="output.root",help='output filename')
    parser.add_argument('--report','-r',default=5000,type=int,help='report every x events')
    args = parser.parse_args()
    
    in_filenames_with_prefix = ['{}{}'.format(args.prefix,x) for x in args.in_filenames]

    #this is a handy class we have made which manages all the handles and labels needed to access products
    #it gets the products on demand so its fine in this case to declare products you wont use
    products = [] 
    add_product(products,"eles","std::vector<pat::Electron>","slimmedElectrons")
    add_product(products,"phos","std::vector<pat::Photon>","slimmedPhotons")
    add_product(products,"gen","std::vector<pat::PackedGenParticle>","packedGenParticles")
    add_product(products,"trig_res","edm::TriggerResults","TriggerResults::HLT")
    add_product(products,"trig_objs","std::vector< pat::TriggerObjectStandAlone >","slimmedPatTrigger")
    evtdata = EvtData(products,verbose=True)

    #now open the output root file
    out_file = ROOT.TFile(args.out,"RECREATE")
    #create a histogram
    #we will do this for signal and bkg seperately
    #for data we will just fill the signal    
    sigmaIEtaIEtaSigHist = ROOT.TH1D("sigmaIEtaIEtaSigHist",";#sigma_{i#etai#eta};#entries",90,0,0.03)
    sigmaIEtaIEtaBkgHist = ROOT.TH1D("sigmaIEtaIEtaBkgHist",";#sigma_{i#etai#eta};#entries",90,0,0.03)
    dEtaInSeedSigHist = ROOT.TH1D("dEtaInSeedSigHist",";#Delta#eta_{in}^{seed};#entries",120,-0.03,0.03)
    dEtaInSeedBkgHist = ROOT.TH1D("dEtaInSeedBkgHist",";#Delta#eta_{in}^{seed};#entries",120,-0.03,0.03)
    chargedIsoSigHist = ROOT.TH1D("chargedIsoSigHist",";charged isol [GeV];#entries",160,0,40)
    chargedIsoBkgHist = ROOT.TH1D("chargedIsoBkgHist",";charged isol [GeV];#entries",160,0,40)
    etSigHist = ROOT.TH1D("etSigHist",";E_{T} [GeV];#entries",50,25,75)
    etBkgHist = ROOT.TH1D("etBkgHist",";E_{T} [GeV];#entries",50,25,75)
    
    trig_psetid = None
    trig_indx = None
    trig_name = "HLT_Ele32_WPTight_Gsf_v"

    events = Events(in_filenames_with_prefix)
    nr_events = events.size()

    for event_nr,event in enumerate(events):
        if event_nr%args.report==0:
            print("processing event {} / {}".format(event_nr,nr_events))
        evtdata.get_handles(event)
        trig_res = evtdata.get("trig_res")
        trig_names = events.object().triggerNames(trig_res)
        if trig_psetid != trig_names.parameterSetID():
            trig_psetid = trig_names.parameterSetID()
            trig_indx = get_trig_indx(trig_name,trig_names)
        trig_pass = trig_res[trig_indx].accept() if trig_indx < trig_res.size() else False
        if not trig_pass:
            continue
            
        trig_objs = evtdata.get("trig_objs")
        eg_trig_objs = [x for x in trig_objs if x.collection().startswith("hltEgammaCandidates")]
        for trig_obj in eg_trig_objs:
            trig_obj.unpackFilterLabels(event.object(),trig_res)

        for ele in evtdata.get("eles"):
            #we will plot this for 25 GeV barrel electrons in our example
            #note we use supercluster eta to for geometric cuts as its from 0,0,0 and therefore represents the detector coordinates
            #electron eta is used for physics quantities as it is w.r.t to the vtx position of the electron and just represents its true eta
            #summary: SC eta = detector eta, electron eta = physics eta

            matched_objs = match_trig_objs(ele.superCluster().eta(),ele.superCluster().phi(),eg_trig_objs)
            pass_trig_obj = False
            for obj in matched_objs:
                if obj.filter("hltEle32WPTightGsfTrackIsoFilter"):
                    pass_trig_obj = True
                    break

            if ele.et()>25 and abs(ele.superCluster().eta())<1.4442 and pass_trig_obj:
                is_data = event.object().event().isRealData()
                #first fill for data & gen matched MC
                if is_data or GenTools.match_to_gen(ele.eta(),ele.phi(),evtdata.get("gen"))[0]:
                    sigmaIEtaIEtaSigHist.Fill(ele.full5x5_sigmaIetaIeta())
                    etSigHist.Fill(ele.et())
                    dEtaInSeedSigHist.Fill(ele.deltaEtaSeedClusterTrackAtVtx())
                    chargedIsoSigHist.Fill(ele.pfIsolationVariables().sumChargedHadronPt)
                #now fill for MC which has not been gen matched 
                elif not is_data: 
                    sigmaIEtaIEtaBkgHist.Fill(ele.full5x5_sigmaIetaIeta())
                    etBkgHist.Fill(ele.et())
                    dEtaInSeedBkgHist.Fill(ele.deltaEtaSeedClusterTrackAtVtx())
                    chargedIsoBkgHist.Fill(ele.pfIsolationVariables().sumChargedHadronPt)

    out_file.Write()
    
    



