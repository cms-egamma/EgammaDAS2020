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
import itertools
import math
import time

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
def get_bincontent_2D(hist,x,y):
    x_bin = hist.GetXaxis().FindBin(x)
    y_bin = hist.GetYaxis().FindBin(y)
    if x_bin<0 or y_bin<0: 
        return None
    else:
        return hist.GetBinContent(x_bin,y_bin)

def rm_hlt_version(hltname):
    """
    removes the version number from an HLT name: HLT_Blah_v4 -> HLT_Blah_v
    """
    match = re.search(r"([\w]*_v)([\d]*$)")
    if match:
        return str(match.group(1))
    else:
        return str(hltname)

def get_trig_indx(selected_name,trig_names):
    """
    returns the index of the trigger (if it exists, otherwise returns the size of trig_names, yes this is a bit c++)
    the annoying thing here is we want to match all HLT_Blah_v*
    there are several ways to do it, the quickest is see which starts with HLT_Blah_v
    
    note one can be smarter and cache the result and update when trig_names.parameterSetID() changes 
    """
    for idx,name in enumerate(trig_names.triggerNames()):
        if name.startswith(selected_name):
            return idx
    return len(trig_names)

def match_trig_objs(eta,phi,max_dr,trig_objs):    
    max_dr2 = max_dr*max_dr
    matched_objs = [x for x in trig_objs if MathTools.cal_delta_r2(eta,phi,obj.eta(),obj.phi()) < max_dr2]
    return matched_objs


def pass_filter(eta,phi,max_dr,filter_name,trig_objs):
    """
    so here we loop through the trigger objects and for each one matching our object in dR
    we check its filters and see if it passed the filter we are interested in
    an object can correspond to many trigger objects (seeded e/gamma, unseeded e/gamma, jet, taus, etc)
    and each trigger object can fire many filters

    this function could be much more efficient but its fine for now
    """
    
    max_dr2 = max_dr*max_dr
    for obj in trig_objs:
        if MathTools.cal_delta_r2(eta,phi,obj.eta(),obj.phi()) < max_dr2:
#           if ROOT.reco.deltaR2(eta,phi,obj.eta(),obj.phi()) < max_dr2:
            if obj.filter(filter_name):
                return True
    return False
        


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
    parser.add_argument('--evtweight','-w',default=1.0,type=float,help='weight for histograms, 1.0 for data, ')
    parser.add_argument('--maxevents','-n',type=int,default=-1,help='max events <0, process all')
    parser.add_argument('--idsf_file','-s',type=str,default=None,help='id sf file')
    parser.add_argument('--recosf_file','-r',type=str,default=None,help='reco sf file')
    parser.add_argument('--out','-o',default="output.root",help='output filename')
    args = parser.parse_args()
    
    in_filenames_with_prefix = CoreTools.get_filenames(args.in_filenames,args.prefix)

    #this is a handy class we have made which manages all the handles and labels needed to access products
    #it gets the products on demand so its fine in this case to declare products you wont use
    products = [] 
    add_product(products,"eles","std::vector<pat::Electron>","slimmedElectrons")
    add_product(products,"trig_res","edm::TriggerResults","TriggerResults::HLT")
    add_product(products,"trig_objs","std::vector<pat::TriggerObjectStandAlone>","slimmedPatTrigger")
    
    evtdata = EvtData(products,verbose=True)

    #ID scale factor histogram
    if args.idsf_file:
        idsf_file  = ROOT.TFile(args.idsf_file,"READ") 
        idsf_hist = idsf_file.EGamma_SF2D 
    else:
        idsf_hist = None

    #reco scale factor histogram
    if args.recosf_file:
        recosf_file  = ROOT.TFile(args.recosf_file,"READ") 
        recosf_hist = recosf_file.EGamma_SF2D 
    else:
        recosf_hist = None


        

    #now open the output root file
    out_file = ROOT.TFile(args.out,"RECREATE")
    #create a histogram
    ROOT.massAllHist = ROOT.TH1D("massAllHist",";M_{ee} (GeV);#entries;",110,40,150)
    ROOT.massTrigHist = ROOT.TH1D("massTrigHist",";M_{ee} (GeV);#entries;",110,40,150)
    ROOT.massIDSelHist = ROOT.TH1D("massIDSelHist",";M_{ee} (GeV);#entries;",110,40,150)
    ROOT.massIDSFCorrHist = ROOT.TH1D("massIDSFCorrHist",";M_{ee} (GeV);#entries;",110,40,150)
    ROOT.massESSCorrHist =  ROOT.TH1D("massESSCorrHist",";M_{ee} (GeV);#entries;",110,40,150)

    id_name = "cutBasedElectronID-Fall17-94X-V2-tight"
    trig_name = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"
    trig_filter_1st = "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"
    trig_filter_2nd = "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"
    
    events = Events(in_filenames_with_prefix)
    nr_events = events.size()
    trig_psetid = None
    trig_indx = None
    for event_nr,event in enumerate(events):
        if args.maxevents >=0 and event_nr>=args.maxevents:
            break
        if event_nr%10000==0:
            print("processing event {} / {} {}".format(event_nr,nr_events,time.ctime()))
        evtdata.get_handles(event)
        eles = evtdata.get("eles")

        trig_res = evtdata.get("trig_res")
        trig_names = events.object().triggerNames(trig_res)
        if trig_psetid != trig_names.parameterSetID():
            trig_psetid = trig_names.parameterSetID()
            trig_indx = get_trig_indx(trig_name,trig_names)       
        
        
        trig_pass = trig_res[trig_indx].accept() if trig_indx<trig_res.size() else False
        
        eg_trig_objs = []

        for ele1,ele2 in itertools.combinations(eles,2):
            if (ele1.et()>25 and abs(ele1.superCluster().eta())<1.4442 and
                ele2.et()>25 and abs(ele2.superCluster().eta())<1.4442 ):
                #note we have a problem with p4() due to how its declared in the c++
                #easier to just use polarP4 
                mass = (ele1.polarP4()+ele2.polarP4()).mag()
                ROOT.massAllHist.Fill(mass,args.evtweight)
                
                #this is bad practice but otherwise I'm going to get a very nested code
                #so we continue when we fail a selection
                #there are better ways but I wanted to keep this simple
                if not (ele1.electronID(id_name) and ele2.electronID(id_name) ):
                    continue
                ROOT.massIDSelHist.Fill(mass,args.evtweight) 
                #normally we would have the trigger selection outside the electron loop 
                #but here I want to show its effect on the mass
            
                if not eg_trig_objs:
                     trig_objs = evtdata.get("trig_objs")
                     #now we need to unpack the labels so we can determine which filters an object passed
                     #note, when doing c++ you would need to make a copy of the trigger object
                     #as products from the event are constant (for very good reason!)
                     #python has no concept of this...
                     #to save time we wil only do this for e/gamma trigger objects which are all from 
                     #the collections hltEgammaCandidates or hltEgammaCandidatesUnseeded (unseeded = no L1 match requirement)
                     eg_trig_objs = [x for x in trig_objs if x.collection().startswith("hltEgammaCandidates")]
                     for trig_obj in eg_trig_objs:
                         trig_obj.unpackFilterLabels(event.object(),trig_res)
                
 

                ele1_pass_leading = pass_filter(ele1.superCluster().eta(),ele1.superCluster().phi(),0.1,trig_filter_1st,eg_trig_objs)           
                ele2_pass_leading = pass_filter(ele2.superCluster().eta(),ele2.superCluster().phi(),0.1,trig_filter_1st,eg_trig_objs)
                ele1_pass_subleading = pass_filter(ele1.superCluster().eta(),ele1.superCluster().phi(),0.1,trig_filter_2nd,eg_trig_objs)
                ele2_pass_subleading = pass_filter(ele2.superCluster().eta(),ele2.superCluster().phi(),0.1,trig_filter_2nd,eg_trig_objs)
                trig_pass_pair = True if (ele1_pass_leading and ele2_pass_subleading) or (ele2_pass_leading and ele1_pass_subleading) else False
               
                    
                if not (trig_pass and trig_pass_pair):
                    continue
                ROOT.massTrigHist.Fill(mass,args.evtweight) 

                id_sf_weight = 1.
                reco_sf_weight = 1.
                #now do ID scale factors
                if idsf_hist:
                    ele1_idsf = get_bincontent_2D(idsf_hist,ele1.superCluster().eta(),ele1.et())
                    ele2_idsf = get_bincontent_2D(idsf_hist,ele2.superCluster().eta(),ele2.et())
                    id_sf_weight = ele1_idsf*ele2_idsf
                if recosf_hist:
                    ele1_recosf = get_bincontent_2D(recosf_hist,ele1.superCluster().eta(),ele1.et())
                    ele2_recosf = get_bincontent_2D(recosf_hist,ele2.superCluster().eta(),ele2.et())
                    reco_sf_weight = ele1_recosf*ele2_recosf

                ROOT.massIDSFCorrHist.Fill(mass,args.evtweight*id_sf_weight*reco_sf_weight)

                #now do energy corrections
                massCorr = mass*math.sqrt(ele1.userFloat("ecalTrkEnergyPostCorr")/ele1.energy()*
                                          ele2.userFloat("ecalTrkEnergyPostCorr")/ele2.energy())

                ROOT.massESSCorrHist.Fill(massCorr,args.evtweight*id_sf_weight)

                            
    out_file.Write()
    
    
