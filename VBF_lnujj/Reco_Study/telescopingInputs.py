import ROOT
import copy
import math
#from math import *
from array import array

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetPalette(1);


############################################################
############################################
# Job steering #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--reco', action='store_true', dest='reco', default=False, help='no X11 windows')
parser.add_option('--gen', action='store_true', dest='gen', default=False, help='no X11 windows')
(options, args) = parser.parse_args()


def makeROCFromHisto(hsig,hbkg,LtoR):
    
    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);
        
    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();
    
    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    for i in range(1,nbins+1):
        
        effBkg = 0;
        effSig = 0;
        
        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;
        
        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;
        
        print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig , "i  ", i;
        
        xval.append( effSig );
        yval.append( 1-effBkg );
        ctr = ctr + 1;
    
    tg = ROOT.TGraph( nbins, xval, yval );
    return tg;

if __name__ == '__main__':
    
    
    fns = "";
    fnb = "";
    if options.reco:
        fns = "Tree_TTbar_vbf.root";
        fnb = "Tree_TTbar_vbf.root";
    if options.gen:
        fns = "Tree_TTbar_vbf.root";
        fnb = "Tree_TTbar_vbf.root";
        

    fileS = ROOT.TFile(fns);
    fileB = ROOT.TFile(fnb);
    tS = fileS.Get("QGl");
    tB = fileB.Get("QGl");


    xs_bkg = 34.1*1.5*1000.; #fb, k-factor of 1.5
    xs_sig = 0.415*0.577*1000.; #fb
    
    nSig = 99529;
    nBkg = 2626868;
    
    LUMI = 20;
    
    # jet of interest (e.g. ak5);
    joi = 19;
    
    cc_joi_sig = 0;
    cc_joi_bkg = 0;
    
    mass_lo = 105;
    mass_hi = 135;
    ptcut = 25;
    
    mpeak = 115;
    
    
    h_sig_absm = ROOT.TH1F("h_sig_absm",";QGl (dis); N", 100,0,1);
    h_bkg_absm = ROOT.TH1F("h_bkg_absm",";QGl (dis); N", 100,0,1);
                        
    ###################################
    print "tS.GetEntries() = ", tS.GetEntries()
    for i in range(tS.GetEntries()):
        tS.GetEntry(i);
        #if i%10000 == 0: print "i = ", i;
        #if i == 5000: break;

        nInterp = 1;
        nPassed = 0;
        z = nPassed/float(nInterp);
        #print nPassed, nInterp, z
        #h_sig_absm.Fill( math.fabs( 126 ) );
	tS.Draw("h_quarklike_QGl>>h_sig_absm","(h_quarklike_pt[0]>30 && h_quarklike_pt[0]>-1)","goff" );

    print "tB.GetEntries() = ", tB.GetEntries()
    for i in range(tB.GetEntries()):
        tB.GetEntry(i);
        #if i%10000 == 0: print "i = ", i;
        #if i == 5000: break;
        
        nInterp = 1;
        nPassed = 0;
        z = nPassed/float(nInterp);
        #h_bkg_absm.Fill( math.fabs( 126 ) );
        tB.Draw("h_gluonlike_QGl>>h_bkg_absm","(h_gluonlike_pt[0]>30 && h_gluonlike_pt[0]>-1)","goff" );

        
    #########################################################
    # save histos
    odir = "";
    if options.reco: odir = "RECO";
    if options.gen: odir = "GEN";
    h_sig_absm.SaveAs("histsForQGl/"+odir+"/h_sig_absm.root");
    h_bkg_absm.SaveAs("histsForQGl/"+odir+"/h_bkg_absm.root");

    # draw histos
    sigHistos = [h_sig_absm];
    bkgHistos = [h_bkg_absm];
        
    for i in range(len(bkgHistos)):
        bkgHistos[i].SetLineColor(2);
        sigHistos[i].Scale(1./sigHistos[i].Integral());
        if bkgHistos[i].Integral() > 0: bkgHistos[i].Scale(1./bkgHistos[i].Integral());
        
        
    c_absm = ROOT.TCanvas("c_absm","c_absm",1000,800);
    h_sig_absm.Draw();
    h_bkg_absm.Draw("sames");
    c_absm.SaveAs("histsForQGl/"+odir+"/c_absm.eps");
        
        
    ## make some ROCs
    roc_B = makeROCFromHisto(h_sig_absm,h_bkg_absm,True);
    #roc_A = makeROCFromHisto(h_wgt_sig,h_wgt_bkg,True);

    roc_B.SetLineColor(2);
        
    canRoc = ROOT.TCanvas("canRoc","canRoc",800,800);
    canRoc.cd();
    hrl = canRoc.DrawFrame(0,0,1.0,1.0);
    hrl.GetXaxis().SetTitle("#epsilon_{sig}");
    hrl.GetYaxis().SetTitle("1 - #epsilon_{bkg}");
    roc_B.Draw();
    #roc_A.Draw();
    canRoc.SaveAs("histsForQGl/"+odir+"/canRoc.eps");
