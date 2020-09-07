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

import EgammaUser.EgammaDAS2020.CoreTools as CoreTools
import EgammaUser.EgammaDAS2020.MathTools as MathTools

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection

""" 
This is an example of how to make a z mass in nano

"""
def get_bincontent_2D(hist,x,y):
    x_bin = hist.GetXaxis().FindBin(x)
    y_bin = hist.GetYaxis().FindBin(y)
    if x_bin<0 or y_bin<0: 
        return None
    else:
        return hist.GetBinContent(x_bin,y_bin)

def pass_filter(ele,trig_objs,filter_nr,min_pt=0,max_dr=0.1):
    max_dr2 = max_dr*max_dr
    for trig_obj in trig_objs:
        if (trig_obj.id ==11 and 
            trig_obj.pt >= min_pt and
            MathTools.cal_delta_r2(ele.eta+ele.deltaEtaSC,ele.phi,trig_obj.eta,trig_obj.phi) < max_dr2 and
            (trig_obj.filterBits & (0x1 << filter_nr) ) !=0 ):
            return True
    return False

def is_good_lumi(grl_dict,run,lumi):
    if str(run) in grl_dict:
        lumi_ranges = grl_dict[str(run)]
        for lumi_range in lumi_ranges:
            if lumi>=lumi_range[0] and lumi<=lumi_range[1]:
                return True
    return False


if __name__ == "__main__":


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
    parser.add_argument('--grl','-g',default=None,help='good lumi json')
    args = parser.parse_args()

    max_events = args.maxevents if args.maxevents>0 else ROOT.TTree.kMaxEntries

    in_filenames_with_prefix = CoreTools.get_filenames(args.in_filenames,args.prefix)

    Events = ROOT.TChain("Events","chain");
    for x in in_filenames_with_prefix:
        print("adding {}".format(x))
        Events.Add(x)
    
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

    id_nr = 4
    trig_name = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"
        
    trig_psetid = None
    trig_indx = None

    Events.Draw(">>elist","Sum$(Electron_pt>25 && abs(Electron_eta+Electron_deltaEtaSC)<1.4442)>=2","entrylist goff",max_events)
    elist = ROOT.gDirectory.Get('elist')
    elist.SetDirectory(0) #removing it from the file not to write it out
    Events.SetEntryList(elist)

    branchsel = BranchSelection("EgammaUser/EgammaDAS2020/data/nano_electron_trig_branches.txt")
    branchsel.selectBranches(Events)

    grl_dict = None
    if args.grl:
        with open(args.grl,'r') as f:
            grl_dict = json.load(f)

    
    nr_events = elist.GetN()
    for event_nr in range(0,nr_events):

        if args.maxevents >=0 and event_nr>=args.maxevents:
            break
        if event_nr%10000==0:
            print("processing event {} / {} {}".format(event_nr,nr_events,time.ctime()))
        
        entry_nr = Events.GetEntryNumber(event_nr)
        Events.GetEntry(entry_nr)
        event = Events

        good_lumi = is_good_lumi(grl_dict,event.run,event.luminosityBlock) if grl_dict else True
                
        eles = Collection(event,"Electron")
        trig_objs = Collection(Events,"TrigObj")

        trig_pass = getattr(event,trig_name)
        
        for ele1,ele2 in itertools.combinations(eles,2):
            if (ele1.pt>25 and abs(ele1.eta+ele1.deltaEtaSC)<1.4442 and
                ele2.pt>25 and abs(ele2.eta+ele2.deltaEtaSC)<1.4442 ):
                #note we have a problem with p4() due to how its declared in the c++
                #easier to just use polarP4 
                mass = math.sqrt(2*ele1.pt*ele2.pt*(math.cosh(ele1.eta-ele2.eta) - math.cos(ele1.phi-ele2.phi)))
                ROOT.massAllHist.Fill(mass,args.evtweight)
                
                #this is bad practice but otherwise I'm going to get a very nested code
                #so we continue when we fail a selection
                #there are better ways but I wanted to keep this simple
                if not (ele1.cutBased>=id_nr and ele2.cutBased>=id_nr):
                    continue
                ROOT.massIDSelHist.Fill(mass,args.evtweight) 
                #normally we would have the trigger selection outside the electron loop 
                #but here I want to show its effect on the mass
             

                ele1_pass_leading = pass_filter(ele1,trig_objs,0,min_pt=23)
                ele2_pass_leading = pass_filter(ele2,trig_objs,0,min_pt=23)
                ele1_pass_subleading = pass_filter(ele1,trig_objs,0,min_pt=12)
                ele2_pass_subleading = pass_filter(ele2,trig_objs,0,min_pt=12)
                trig_pass_pair = True if (ele1_pass_leading and ele2_pass_subleading) or (ele2_pass_leading and ele1_pass_subleading) else False
                                   
                if not (trig_pass and trig_pass_pair):
                    continue
                ROOT.massTrigHist.Fill(mass,args.evtweight) 

                id_sf_weight = 1.
                reco_sf_weight = 1.
                #now do ID scale factors
                if idsf_hist:
                    ele1_idsf = get_bincontent_2D(idsf_hist,ele1.eta+ele1.deltaEtaSC,ele1.pt)
                    ele2_idsf = get_bincontent_2D(idsf_hist,ele2.eta+ele2.deltaEtaSC,ele2.pt)
                    id_sf_weight = ele1_idsf*ele2_idsf
                if recosf_hist:
                    ele1_recosf = get_bincontent_2D(recosf_hist,ele1.eta+ele1.deltaEtaSC,ele1.pt)
                    ele2_recosf = get_bincontent_2D(recosf_hist,ele2.eta+ele2.deltaEtaSC,ele2.pt)
                    reco_sf_weight = ele1_recosf*ele2_recosf

                ROOT.massIDSFCorrHist.Fill(mass,args.evtweight*id_sf_weight*reco_sf_weight)

                #now energy corrections are already applied
                ROOT.massESSCorrHist.Fill(mass,args.evtweight*id_sf_weight)

                            
    out_file.Write()
    
    
