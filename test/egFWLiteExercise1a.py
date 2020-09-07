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
    
    events = Events(in_filenames_with_prefix)
    nr_events = events.size()
    for event_nr,event in enumerate(events):
        if event_nr%args.report==0:
            print("processing event {} / {}".format(event_nr,nr_events))
        evtdata.get_handles(event)
        for ele in evtdata.get("eles"):
            #we will plot this for 25 GeV barrel electrons in our example
            #note we use supercluster eta to for geometric cuts as its from 0,0,0 and therefore represents the detector coordinates
            #electron eta is used for physics quantities as it is w.r.t to the vtx position of the electron and just represents its true eta
            #summary: SC eta = detector eta, electron eta = physics eta
            if ele.et()>25 and abs(ele.superCluster().eta())<1.4442:
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
    
    



