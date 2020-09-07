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

def set_style_att(hist,color=None,line_width=None,marker_style=None,line_style=None,marker_size=None):
    if color!=None:
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
    if line_width!=None:
        hist.SetLineWidth(line_width)
    if line_style!=None:
        hist.SetLineWidth(line_style)
    if marker_style!=None:
        hist.SetMarkerStyle(marker_style)
    if marker_size!=None:
        hist.SetMarkerSize(marker_size)


def adjust_yaxis_dist(hist1,hist2):
    max_val = max(hist1.GetMaximum(),hist2.GetMaximum())
    hist1.GetYaxis().SetRangeUser(0.1,max_val*1.3)
    
def adjust_yaxis_eff(hist1,hist2):
    max_val = max(hist1.GetMaximum(),hist2.GetMaximum())
    min_val = min(hist1.GetMinimum(),hist2.GetMinimum())
    print(min_val,max_val)

    if min_val>0.95:
        y_min = 0.9
        y_max = 1.05
    elif min_val>0.85:
        y_min = 0.8
        y_max = 1.1
    elif min_val>0.75:
        y_min = 0.7
        y_max = 1.2
    elif min_val>0.6:
        y_min = 0.5
        y_max = 1.2
    else:
        y_min = 0
        y_max = 1.3
    hist1.GetYaxis().SetRangeUser(y_min,y_max)
    
def adjust_yaxis(hist1,hist2,is_eff=False):
    if is_eff:
        adjust_yaxis_eff(hist1,hist2)
    else:
        adjust_yaxis_dist(hist1,hist2)

def plot_with_ratio(numer,numer_label,denom,denom_label,div_opt="",labels_text=[]):

    set_style_att(numer,color=1,marker_style=8)
    set_style_att(denom,color=4,line_width=2)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",900,750)
    c1.cd()
    spectrum_pad = ROOT.TPad("spectrumPad","newpad",0.01,0.30,0.99,0.99)
    spectrum_pad.Draw() 
    spectrum_pad.cd()
    xaxis_title = denom.GetXaxis().GetTitle()
    denom.GetXaxis().SetTitle()
    denom.Draw("HISTE")
    numer.Draw("EP SAME")
    leg = ROOT.TLegend(0.115,0.766,0.415,0.888)
    leg.AddEntry(numer,numer_label,"LP")
    leg.AddEntry(denom,denom_label,"LP")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()

    labels = []    
    for cut_nr,cut_label in labels_text:
        y_max = 0.882-0.06*cut_nr
        y_min = 0.882-0.06*(cut_nr+1)
        label = ROOT.TPaveLabel(0.681,y_min,0.894,y_max,cut_label,"brNDC")
        label.SetFillStyle(0)
        label.SetBorderSize(0)
        label.SetTextFont(42)
        label.SetTextAlign(12);
        label.SetTextSize(0.657895);
        label.Draw()
        labels.append(label)

    
    c1.cd()
    ratio_pad = ROOT.TPad("ratioPad", "newpad",0.01,0.01,0.99,0.33)
    ratio_pad.Draw()
    ratio_pad.cd()
    ratio_pad.SetGridy()
    ratio_pad.SetTopMargin(0.05)
    ratio_pad.SetBottomMargin(0.3)
    ratio_pad.SetFillStyle(0)
    ratio_hist = numer.Clone("ratioHist")
    ratio_hist.SetDirectory(0)
    ratio_hist.Sumw2()
    ratio_hist.Divide(numer,denom,1,1,div_opt)

    set_style_att(ratio_hist,color=1,marker_style=8)
    ratio_hist.SetTitle("")
    #  ratio_hist.GetXaxis().SetLabelSize(ratio_hist.GetXaxis().GetLabelSize()*(0.99-0.33)/0.33)
    ratio_hist.GetXaxis().SetLabelSize(0.1)
    ratio_hist.GetXaxis().SetTitleSize(0.1)
    ratio_hist.GetXaxis().SetTitle(xaxis_title)
    ratio_hist.GetYaxis().SetLabelSize(0.1)
    ratio_hist.GetYaxis().SetTitleSize(0.1)
    ratio_hist.GetYaxis().SetTitleOffset(0.3) 
    ratio_hist.GetYaxis().SetTitle("ratio")   
    ratio_hist.GetYaxis().SetRangeUser(0.5,1.5)
    ratio_hist.GetYaxis().SetNdivisions(505)
    
    
    ratio_hist.Draw("EP")
    spectrum_pad.cd()
    return c1,spectrum_pad,ratio_pad,ratio_hist,leg,labels

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='prints E/gamma pat::Electrons/Photons')
    parser.add_argument('--tar','-t',type=str,required=True,help='target filename')
    parser.add_argument('--tar_name',type=str,default="target",help='target leg entry')
    parser.add_argument('--ref','-r',type=str,required=True,help='ref filename')
    parser.add_argument('--ref_name',type=str,default="reference",help='reference leg entry')
    args = parser.parse_args()

    if hasattr(args,"tar"):
        tar = ROOT.TFile(args.tar,"READ")

    if hasattr(args,"ref"):
        ref = ROOT.TFile(args.ref,"READ")
    
