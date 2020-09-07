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
import time

from DataFormats.FWLite import Events, Handle
from EgammaUser.EgammaDAS2020.EvtData import EvtData, EvtHandles, add_product

import EgammaUser.EgammaDAS2020.CoreTools as CoreTools
import EgammaUser.EgammaDAS2020.GenTools as GenTools
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.branchselection import BranchSelection

""" 
This is a simple example of nanoAOD analysis (or similar tree)

"""


if __name__ == "__main__":

    
    #now read the cmd line to find our input filename
    #local file "file:<filename>
    #remote file "root:// 
    #note unlike CMSSW proper its not smart enough to resolve to fall over to xrootd automatically
    #so it has to be specified
    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('in_filenames',nargs="+",help='input filenames')
    parser.add_argument('--prefix','-p',default='',help='file prefix')
    parser.add_argument('--out','-o',default='output.root',help='output filename')
    parser.add_argument('--report','-r',default=25000,type=int,help='report every x events')
    parser.add_argument('--max_events','-n',default=-1,type=int,help='max number of events to run over')
    args = parser.parse_args()

    max_events = args.max_events if args.max_events>0 else ROOT.TTree.kMaxEntries

    in_filenames_with_prefix = CoreTools.get_filenames(args.in_filenames,args.prefix)

    Events = ROOT.TChain("Events","chain");
    for x in in_filenames_with_prefix:
        print("adding {}".format(x))
        Events.Add(x)
    
    

    #now open the output root file
    out_file = ROOT.TFile(args.out,"RECREATE")
    #create a histogram, we will make 3 just to demonstrate 3 ways you can do this, normally you would just do one. The best way is probably to follow the root exercise with RDataFrame
    sigmaIEtaIEtaSigHist = ROOT.TH1D("sigmaIEtaIEtaSigHist",";#sigma_{i#etai#eta};#entries;",90,0,0.03)
    sigmaIEtaIEtaBkgHist = ROOT.TH1D("sigmaIEtaIEtaBkgHist",";#sigma_{i#etai#eta};#entries;",90,0,0.03)
    chargedIsoSigHist = ROOT.TH1D("chargedIsoSigHist",";charged isol [GeV];#entries",160,0,40)
    chargedIsoBkgHist = ROOT.TH1D("chargedIsoBkgHist",";charged isol [GeV];#entries",160,0,40)
    etSigHist = ROOT.TH1D("etSigHist",";E_{T} [GeV];#entries",50,25,75)
    etBkgHist = ROOT.TH1D("etBkgHist",";E_{T} [GeV];#entries",50,25,75)

    

    #Nano is just a flat tree, therefore this ways work for any similar flat tree and 
    #are not specific to nano

    nr_events = Events.GetEntries()

    print("start time",time.ctime())
    #now we will demonstrate running outside the Draw command which is more flexable
    #however this is slooow so we do things to speed it up like 
    #pre skim things as it is slower with an entry list
    #and dropping branches we dont need so we dont waste time reading them
    #note I highly suggest you look into RDataFrames as that should be much faster
    Events.Draw(">>elist","Sum$(Electron_pt>25 && abs(Electron_eta+Electron_deltaEtaSC)<1.4442)>=1","entrylist goff",max_events)
    elist = ROOT.gDirectory.Get('elist')
    elist.SetDirectory(0) #removing it from the file not to write it out
    Events.SetEntryList(elist)

    branchsel = BranchSelection("EgammaUser/EgammaDAS2020/data/nano_electron_branches.txt")
    branchsel.selectBranches(Events)
    
    nr_events = elist.GetN()
    for event_nr in range(0,nr_events):
        entry_nr = Events.GetEntryNumber(event_nr)
        Events.GetEntry(entry_nr)
        event = Events
        if event_nr%args.report==0:
            print("processing event {} / {} {}".format(event_nr,nr_events,time.ctime()))

        eles = Collection(event,"Electron")
        for ele in eles:
            if ele.pt>25 and abs(ele.eta+ele.deltaEtaSC)<1.4442:
                if ele.genPartIdx>=0:
                    sigmaIEtaIEtaSigHist.Fill(ele.sieie)
                    etSigHist.Fill(ele.pt)
                    chargedIsoSigHist.Fill(ele.pfRelIso03_chg*ele.pt)
                else:
                    sigmaIEtaIEtaBkgHist.Fill(ele.sieie)
                    etBkgHist.Fill(ele.pt)
                    chargedIsoBkgHist.Fill(ele.pfRelIso03_chg*ele.pt)

    
    
    out_file.Write()
