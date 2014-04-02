/*//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Channel: VBF WW semileptonic Muon channel (W -> uv, W -> jj, gamma)

Usage: Reads in multiple ROOT files (RD Trees) representing data, MC
       signal, and background.  The background normalization is done
       using theoretical cross sections, observed luminosity, and
       sample sizes (user can also add K-factors, etc.).  User
       specifies kinematical distribution to observe, as well as
       histogram range and binning.  User can also specify the
       Systematic/Statistical uncertainty.

Output: ROOT file containing 1D histograms for:
	(1) Observed data (data_obs)
	(2) SM backgrounds + SM signal (background)
	(3) Background + uncertainty (background_mudijet_backshapeUp)
	(4) Background - uncertainty (background_mudijet_backshapeDown)
	(5-?) Anomalous signal - SM signal (signal_??)

Syntax: root -l -b -q mkROOTaqgc.C

///////////////////////////////////////////////////////////////////////
*//////////////////////////////////////////////////////////////////////


/////////////////////
// Load header files:
#include <iostream>
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"


///////////////////////////////////////////
// Define Type to store histogram settings:
struct plotVar_t {
  char* plotvar;
  double MINRange;
  double MAXRange;
  int    NBINS;
  int    slog;
  char* xlabel;
};


///////////////////////
// Begin Main function:
void InputPre(){


  //////////////////////////////////////
  // Set MC normalization scale factors:
  const double intLUMI = 19297;
/*
// WA+jets norm - data driven : 
//  const double WAJets_scale  = 1.171 * 9.37246 * intLUMI/1049473;
//  const double WAJetsPT100_scale   = 1.171* 0.423028445 *intLUMI/369309; // V
  const double WAJets_scale  = 9.37246 * intLUMI/1049473;
  const double WAJetsPT100_scale   =  0.423028445 *intLUMI/369309; // V
  const double ZAJets_scale  = 0.63196 * intLUMI/979207;
  const double ZZ_scale      = 8.059 * intLUMI/9702850;

  const double ttbarA_scale  = 1.44 * intLUMI/604113;
  const double SToppS_scale  = 1.75776  * intLUMI/139974;
  const double SToppT_scale  = 30.0042  * intLUMI/1935066;
  const double SToppTW_scale = 11.1773  * intLUMI/493458;
  const double STopS_scale   = 3.89394  * intLUMI/259960;
  const double STopT_scale   = 55.531  * intLUMI/3758221;
  const double STopTW_scale  = 11.1773  * intLUMI/497657;

  // PT dependent k-factor included (same for WWA,WZA aQGC)
  const double WWA_scale     =  2.1* 0.01362 * intLUMI/198777;
  const double WWA2_scale    =  2.1* 0.01409 * intLUMI/199183;
  const double WZA_scale     =  2.1* 0.00578008 * intLUMI/497450;

  const double KOG_m5m5_scale     = 2.1* 0.0470065 * intLUMI/194029;
  const double KOG_m2m5_scale     = 2.1* 0.0249555 * intLUMI/195494;
  const double KOG_m3m5_scale     = 2.1* 0.0302425 * intLUMI/196004;
  const double KOG_p5m5_scale     = 2.1* 0.046106 * intLUMI/184331;
  const double KOG_p2m5_scale     = 2.1* 0.0245765 * intLUMI/197493;
  const double KOG_p3m5_scale     = 2.1* 0.029646 * intLUMI/139920;
*/
                                const double wpj_kfactor = 1.16*1.16; //actual K factor
                        //        const double wpj_kfactor = 1.0;
                               //const double wpj_kfactor = 1.27809;
                              //  const double ttbar_kfactor = 0.955077;
                                const double ttbar_kfactor = 1.0;
                                const double mh126_scale   = 0.0776*(3.0/2.0)* intLUMI/333061; // V
                                const double madww2jets_scale   = 0.03907/2.0* intLUMI/414963; // V
                                const double W4Jets_scale  = 214.0 * intLUMI*wpj_kfactor/4369420;
                                const double W3Jets_scale  = 519.0 * intLUMI*wpj_kfactor/15059503;
                                const double W2Jets_scale  = 1750.0 * intLUMI*wpj_kfactor/33004921;
                                const double W1Jets_scale  = 5400.0 * intLUMI*wpj_kfactor/19871598;
                                const double WW_scale      = 54.838   * intLUMI/9450414; // V
                                const double WZ_scale      = 32.3161   * intLUMI/10000267; // V
                                const double ZZ_scale      = 8.059   * intLUMI/9702850; // V
                                const double QCD_scale     = 364000000 *    0.00037 * intLUMI/7529312 ; // V
                                const double ZJets_scale   = 3503.71  * intLUMI/30209426; // V
                                const double ttbar_scale   = 225.197 * intLUMI*ttbar_kfactor/6893735; // V
                                const double SToppS_scale  = 1.75776  * intLUMI/139974; // V anti-top
                                const double SToppT_scale  = 30.0042  * intLUMI/1935066; // V
                                const double SToppTW_scale = 11.1773  * intLUMI/493458; // V
                                const double STopS_scale   = 3.89394  * intLUMI/259960; // top
                                const double STopT_scale   = 55.531  * intLUMI/3758221; //V
                                const double STopTW_scale  = 11.1773  * intLUMI/497657; // V

  ///////////////////////////////////////////////////////////////////////////////////
  // Specify what kinematical distribution to observe, as well as histogram settings:
  // 
//  plotVar_t pv = {"Photon_Et[iPhoton12]",30,1000,7,3,"Photon E_{T} (GeV)"};
//  plotVar_t pv = {"Photon_Et[iPhoton12]",45,325,7,3,"Photon E_{T} (GeV)"};
//  plotVar_t pv = {"masslvjja",0,1200,12,3,"M_{l\nujj\gamma} (GeV)"};
//  plotVar_t pv = {"masslvjja",800,6800,15,3,"M_{l\nujj\gamma} (GeV)"};
//  plotVar_t pv = {"c2jMass12",10, 160, 10,3,"M_{jj} (GeV)"};
//  plotVar_t pv = {"hvbf_wjj_m",30, 600, 50,3,"M_{jj} (GeV)"};
  plotVar_t pv = {"hvbf_topWm",30, 212, 14,3,"Top M_{jj} (GeV)"};

  if ( !strlen(pv.plotvar) ) break;
  std::cout << TString(pv.plotvar) << "\t"<<pv.MINRange<<"\t" << pv.MAXRange<<"\t" << pv.NBINS<<"\tTHE CUT " << endl;


  ////////////////////////////////
  // Specify event selection cuts:
//  the_cut = TCut("effwt*puwt*(hvbf_event && ((hvbf_wjj_m > 40.0 && hvbf_wjj_m < 65.0) || (hvbf_wjj_m > 95.0 && hvbf_wjj_m < 120.0))&& hvbf_wjet1_btagCSV < 0.679 && hvbf_wjet2_btagCSV < 0.679 )");

//  TCut the_cutE("effwt*puwt*(hvbf_event && ((hvbf_wjj_m > 40.0 && hvbf_wjj_m < 65.0) || (hvbf_wjj_m > 95.0 && hvbf_wjj_m < 120.0))&& hvbf_wjet1_btagCSV < 0.679 && hvbf_wjet2_btagCSV < 0.679 )");

  the_cut = TCut("effwt*puwt*(hvbf_topWm>30 &&hvbf_topWm<212 &&  abs(W_muon_eta)<2.4)");
  TCut the_cutE("effwt*puwt*(hvbf_topWm>30 &&hvbf_topWm<212 &&  abs(W_muon_eta)<2.4)");

//  the_cut = TCut("(hvbf_topWm>0&& abs(W_muon_eta)<2.4)");
//  TCut the_cutE("(hvbf_topWm>0&& abs(W_muon_eta)<2.4)");


 // ! no need for extra PUwt stat err. sumw2 is enough
  //#TCut WA23Jcut(TString(the_cut)+TString("*((W_Photon_pt_gen>115)? 0. : 1.) "));
  //#TCut WA23JPT100cut(TString(the_cut)+TString("*((W_Photon_pt_gen>115)? 1. : 0.) "));

  //////////////////////////////////////////////////////////////////
  // Specify Fake-Photon contribution function (fake rate function):
  //#TCut fkPhoton_cut(TString(the_cutPlj)+TString("*(0.01373795 + 308.9628/(Photon_Et[(iPhoton12plj>-1)? iPhoton12plj : 0]^2.29711))"));

  ///////////////////////////
  // Create output ROOT file:
  TFile f("mu_VBFWW_TTbarControlRegion_corr.root", "RECREATE");

  //////////////////////////////////////////////////
  // Create file pointers for each sample ROOT file:
 TFile *fin2,*madww2jets_file,*mh126_file,*wwShape_file,*wzShape_file,*zzShape_file,*w4jetsShape_file,*w3jetShape_file,*w2jetShape_file,*w1jetShape_file,*ttbar_file,*qcd_file1,*zjets_file,*stops_file,*stopt_file,*stoptW_file;


  //////////////////////////////
  // Open each sample ROOT file:
        fin2            = new TFile("InData/RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root", "read");
        madww2jets_file       = new TFile("InData/RD_mu_mh126_CMSSW532.root", "READ");
        mh126_file       = new TFile("InData/RD_mu_phantom_CMSSW532.root", "READ");
        wwShape_file    = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
        wzShape_file    = new TFile("InData/RD_mu_WZ_CMSSW532.root", "READ");
        zzShape_file    = new TFile("InData/RD_mu_ZZ_CMSSW532.root", "READ");
        w4jetsShape_file = new TFile("InData/RD_mu_W4Jets_old_CMSSW532.root","READ");
        w3jetShape_file = new TFile("InData/RD_mu_W3Jets_CMSSW532.root","READ");
        w2jetShape_file = new TFile("InData/RD_mu_W2Jets_CMSSW532.root","READ");
        w1jetShape_file = new TFile("InData/RD_mu_W1Jets_CMSSW532.root","READ");
        ttbar_file      = new TFile("InData/RD_mu_TTbar_CMSSW532.root", "READ");
        zjets_file      = new TFile("InData/RD_mu_ZpJ_CMSSW532.root", "READ");
        stops_file      = new TFile("InData/RD_mu_STopS_T_CMSSW532.root", "READ");
        stopt_file      = new TFile("InData/RD_mu_STopT_T_CMSSW532.root", "READ");
        stoptW_file     = new TFile("InData/RD_mu_STopTW_T_CMSSW532.root", "READ");
        stopps_file     =  new TFile("InData/RD_mu_STopS_Tbar_CMSSW532.root", "READ");
        stopptW_file    =  new TFile("InData/RD_mu_STopTW_Tbar_CMSSW532.root", "READ");
        stoppt_file     =  new TFile("InData/RD_mu_STopT_Tbar_CMSSW532.root", "READ");


  ///////////////////////////////////////////////////
  // Retrieve ROOT tree with kinematic distributions:
      TTree* treedata = (TTree*) fin2->Get("WJet");
        double nData = treedata->GetEntries();
        std::cout << "ndata =" << nData <<std::endl;
        TTree* treevhmadww2jets  = (TTree*)       madww2jets_file->Get("WJet");
        TTree* treemh126  = (TTree*)       mh126_file->Get("WJet");
        TTree* treeww    = (TTree*)    wwShape_file->Get("WJet");
        TTree* treewz    = (TTree*)    wzShape_file->Get("WJet");
        TTree* treezz    = (TTree*)    zzShape_file->Get("WJet");
        TTree* treew4j    = (TTree*) w4jetsShape_file->Get("WJet");
        TTree* treew3j    = (TTree*) w3jetShape_file->Get("WJet");
        TTree* treew2j    = (TTree*) w2jetShape_file->Get("WJet");
        TTree* treew1j    = (TTree*) w1jetShape_file->Get("WJet");
        TTree* treettb   = (TTree*)      ttbar_file->Get("WJet");
        TTree* treezj    = (TTree*)      zjets_file->Get("WJet");

        TTree* treests   = (TTree*)      stops_file->Get("WJet");
        TTree* treestt   = (TTree*)      stopt_file->Get("WJet");
        TTree* treestw   = (TTree*)     stoptW_file->Get("WJet");
        TTree* tree64 = (TTree*) 	stopps_file->Get("WJet");
        TTree* tree65 = (TTree*) 	stoppt_file->Get("WJet");
        TTree* tree66 = (TTree*) 	stopptW_file->Get("WJet");



  ////////////////////////////////////////////////////////////
  // Create kinematic-distribution histograms for each sample:
        TH1* data_obs  = new TH1D("data_obs",  "data_obs",  pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1madww2jets = new TH1D("th1madww2jets", "th1madww2jets", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1mh126 = new TH1D("th1mh126", "th1mh126", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1ww = new TH1D("th1ww", "th1ww", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1wz = new TH1D("th1wz", "th1wz", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1zz = new TH1D("th1zz", "th1zz", pv.NBINS, pv.MINRange, pv.MAXRange);

        TH1* th1w4jets  = new TH1D("th1w4jets",  "th1w4jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
        TH1* th1w3jets  = new TH1D("th1w3jets",  "th1w3jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
        TH1* th1w2jets  = new TH1D("th1w2jets",  "th1w2jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
        TH1* th1w1jets  = new TH1D("th1w1jets",  "th1w1jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);

        TH1* th1Top = new TH1D("th1Top", "th1Top", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1zjets = new TH1D("th1zjets", "th1zjets", pv.NBINS, pv.MINRange, pv.MAXRange);

        TH1* th1stops = new TH1D("th1stops", "th1stops", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1stopt = new TH1D("th1stopt", "th1stopt", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1stoptw = new TH1D("th1stoptw", "th1stoptw", pv.NBINS, pv.MINRange, pv.MAXRange);

        TH1* th1stopps = new TH1D("th1stopps", "th1stopps", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1stoppt = new TH1D("th1stoppt", "th1stoppt", pv.NBINS, pv.MINRange, pv.MAXRange);
        TH1* th1stopptw = new TH1D("th1stopptw", "th1stopptw", pv.NBINS, pv.MINRange, pv.MAXRange);


  TH1F* result;
  /////////////////////////////////////////////////////////////////////////
  // Specify histograms to store Sum of Squares of Weights for each sample:
  data_obs->Sumw2();
  th1madww2jets->Sumw2();
  th1mh126->Sumw2(); 
  th1ww->Sumw2();
  th1wz->Sumw2();
  th1zz->Sumw2();
  th1w4jets->Sumw2();
  th1w3jets->Sumw2();
  th1w2jets->Sumw2();
  th1w1jets->Sumw2();

  th1Top->Sumw2();
  th1stops->Sumw2();
  th1stopt->Sumw2();
  th1stoptw->Sumw2();
  th1stopps->Sumw2();
  th1stoppt->Sumw2();
  th1stopptw->Sumw2();


  ///////////////////////////////////////////////////////////////////////////////////
  // Fill kinematical distribution for each sample according to event selection cuts:


  std::cout<<"\nFill Data Histogram..."<<std::endl;
  treedata->Draw(TString(pv.plotvar)+TString(">>data_obs"), the_cut, "goff");
  data_obs->AddBinContent(pv.NBINS,data_obs->GetBinContent(pv.NBINS+1));data_obs->SetBinContent(pv.NBINS+1,0.);

//  std::cout<<"Fill Fake Photon Histogram..."<<std::endl;
//  treedata->Draw(TString(pv.plotvar)+TString(">>th1fkdata"), fkPhoton_cut, "goff");
//  treedata->Draw(TString("c2jMass12plj>>th1fkdata"), fkPhoton_cut, "goff");
//  th1fkdata->AddBinContent(pv.NBINS,th1fkdata->GetBinContent(pv.NBINS+1));th1fkdata->SetBinContent(pv.NBINS+1,0.);

  treevhmadww2jets->Draw(TString(pv.plotvar)+TString(">>th1madww2jets"), the_cut, "goff");
  th1madww2jets->AddBinContent(pv.NBINS,th1madww2jets->GetBinContent(pv.NBINS+1));th1madww2jets->SetBinContent(pv.NBINS+1,0.);

  treemh126->Draw(TString(pv.plotvar)+TString(">>th1mh126"), the_cut, "goff");
  th1mh126->AddBinContent(pv.NBINS,th1mh126->GetBinContent(pv.NBINS+1));th1mh126->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill SM Dibosn  Histogram..."<<std::endl;
  treeww->Draw(TString(pv.plotvar)+TString(">>th1ww"), the_cut, "goff");
  th1ww->AddBinContent(pv.NBINS,th1ww->GetBinContent(pv.NBINS+1));th1ww->SetBinContent(pv.NBINS+1,0.);

  treewz->Draw(TString(pv.plotvar)+TString(">>th1wz"), the_cut, "goff");
  th1wz->AddBinContent(pv.NBINS,th1wz->GetBinContent(pv.NBINS+1));th1wz->SetBinContent(pv.NBINS+1,0.);

  treezz->Draw(TString(pv.plotvar)+TString(">>th1zz"), the_cut, "goff");
  th1zz->AddBinContent(pv.NBINS,th1zz->GetBinContent(pv.NBINS+1));th1zz->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill W+Jets Histogram..."<<std::endl;
  treew4j->Draw(TString(pv.plotvar)+TString(">>th1w4jets"), the_cut, "goff");
  th1w4jets->AddBinContent(pv.NBINS,th1w4jets->GetBinContent(pv.NBINS+1));th1w4jets->SetBinContent(pv.NBINS+1,0.);

  treew3j->Draw(TString(pv.plotvar)+TString(">>th1w3jets"), the_cut, "goff");
  th1w3jets->AddBinContent(pv.NBINS,th1w3jets->GetBinContent(pv.NBINS+1));th1w3jets->SetBinContent(pv.NBINS+1,0.);

  treew2j->Draw(TString(pv.plotvar)+TString(">>th1w2jets"), the_cut, "goff");
  th1w2jets->AddBinContent(pv.NBINS,th1w2jets->GetBinContent(pv.NBINS+1));th1w2jets->SetBinContent(pv.NBINS+1,0.);

  treew1j->Draw(TString(pv.plotvar)+TString(">>th1w1jets"), the_cut, "goff");
  th1w1jets->AddBinContent(pv.NBINS,th1w1jets->GetBinContent(pv.NBINS+1));th1w1jets->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill ttbar Histogram..."<<std::endl;
  treettb->Draw(TString(pv.plotvar)+TString(">>th1Top"), the_cut, "goff");
  th1Top->AddBinContent(pv.NBINS,th1Top->GetBinContent(pv.NBINS+1));th1Top->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill Z+Jets Histogram..."<<std::endl;
  treezj->Draw(TString(pv.plotvar)+TString(">>th1zjets"), the_cut, "goff");
  th1zjets->AddBinContent(pv.NBINS,th1zjets->GetBinContent(pv.NBINS+1));th1zjets->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopS Histogram..."<<std::endl;
  treests->Draw(TString(pv.plotvar)+TString(">>th1stops"), the_cut, "goff");
  th1stops->AddBinContent(pv.NBINS,th1stops->GetBinContent(pv.NBINS+1));th1stops->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopT Histogram..."<<std::endl;
  treestt->Draw(TString(pv.plotvar)+TString(">>th1stopt"), the_cut, "goff");
  th1stopt->AddBinContent(pv.NBINS,th1stopt->GetBinContent(pv.NBINS+1));th1stopt->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopTW Histogram..."<<std::endl;
  treestw->Draw(TString(pv.plotvar)+TString(">>th1stoptw"), the_cut, "goff");
  th1stoptw->AddBinContent(pv.NBINS,th1stoptw->GetBinContent(pv.NBINS+1));th1stoptw->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopPS Histogram..."<<std::endl;
  tree64->Draw(TString(pv.plotvar)+TString(">>th1stopps"), the_cut, "goff");
  th1stopps->AddBinContent(pv.NBINS,th1stopps->GetBinContent(pv.NBINS+1));th1stopps->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopPT Histogram..."<<std::endl;
  tree65->Draw(TString(pv.plotvar)+TString(">>th1stoppt"), the_cut, "goff");
  th1stoppt->AddBinContent(pv.NBINS,th1stoppt->GetBinContent(pv.NBINS+1));th1stoppt->SetBinContent(pv.NBINS+1,0.);

  std::cout<<"Fill STopPTW Histogram..."<<std::endl;
  tree66->Draw(TString(pv.plotvar)+TString(">>th1stopptw"), the_cut, "goff");
  th1stopptw->AddBinContent(pv.NBINS,th1stopptw->GetBinContent(pv.NBINS+1));th1stopptw->SetBinContent(pv.NBINS+1,0.);


  ///////////////////////////////////////////////////////////////////////////
/*  // Set uncertainty in each bin of kinematical distribution for each sample:
    float fperr = 0.;
    // adding 20% sys on photon fakerate 
    for(int hi=1;hi<=pv.NBINS;hi++) {
        fperr = (th1fkdata->GetBinError(hi))*(th1fkdata->GetBinError(hi));
        fperr + = 0.04*(th1fkdata->GetBinContent(hi))*(th1fkdata->GetBinContent(hi));
        th1fkdata->SetBinError(hi,sqrt(fperr));
    }

//  Adding Normalization uncertainty to WA+jets
    for(int hi=1;hi<=pv.NBINS;hi++) {
        fperr = (th1wajets->GetBinError(hi))*(th1wajets->GetBinError(hi));
        fperr + = 0.078*0.078*(th1wajets->GetBinContent(hi))*(th1wajets->GetBinContent(hi));
        th1wajets->SetBinError(hi,sqrt(fperr));
        fperr = (th1wajetsPT100->GetBinError(hi))*(th1wajetsPT100->GetBinError(hi));
        fperr + = 0.078*0.078*(th1wajetsPT100->GetBinContent(hi))*(th1wajetsPT100->GetBinContent(hi));
        th1wajetsPT100->SetBinError(hi,sqrt(fperr));
    }
*/

  /////////////////////////
  // Normalize each sample:
  std::cout<<"\nScale Histograms..."<<std::endl;
  th1Top->Scale(ttbar_scale);
  th1stops->Scale(STopS_scale);
  th1stopt->Scale(STopT_scale);
  th1stoptw->Scale(STopTW_scale);
  th1stopps->Scale(SToppS_scale);
  th1stoppt->Scale(SToppT_scale);
  th1stopptw->Scale(SToppTW_scale);
  th1w4jets->Scale(W4Jets_scale );
  th1w3jets->Scale(W3Jets_scale);
  th1w2jets->Scale(W2Jets_scale);
  th1w1jets->Scale(W1Jets_scale);
  th1ww->Scale(WW_scale);
  th1wz->Scale(WZ_scale);
  th1zz->Scale(ZZ_scale);
  th1zjets->Scale(ZJets_scale);

   th1madww2jets->Scale(madww2jets_scale);
   th1mh126->Scale(mh126_scale);

  ///////////////////////////
  // Combine certain samples:

  TH1D *th1sumwjets = (TH1D*)th1w4jets->Clone("th1sumwjets");
  th1sumwjets->Reset();
  th1sumwjets->Add(th1w4jets,1);
  th1sumwjets->Add(th1w3jets,1);
  th1sumwjets->Add(th1w2jets,1);
  th1sumwjets->Add(th1w1jets,1);

  TH1D *th1sumstop = (TH1D*)th1stopt->Clone("th1sumstop");
  th1sumstop->Reset();
  th1sumstop->Add(th1stopptw,1);
  th1sumstop->Add(th1stoppt,1);
  th1sumstop->Add(th1stopps,1);
  th1sumstop->Add(th1stoptw,1);
  th1sumstop->Add(th1stopt,1);
  th1sumstop->Add(th1stops,1);

  TH1D *th1sumdiboson = (TH1D*)th1ww->Clone("th1sumdiboson");
  th1sumdiboson->Reset();
  th1sumdiboson->Add(th1ww,1);
  th1sumdiboson->Add(th1wz,1);
  th1sumdiboson->Add(th1zz,1);

        // Sum all the backgrounds
  TH1D *th1sumbkg = (TH1D*)th1w4jets->Clone("th1sumbkg");
  th1sumbkg->Reset();
  th1sumbkg->Add(th1sumwjets,1);
  th1sumbkg->Add(th1sumdiboson,1);
  th1sumbkg->Add(th1Top,1);
  th1sumbkg->Add(th1sumstop,1);
  th1sumbkg->Add(th1zjets,1);

/*
  ///////////////////////////////////////////
  // Combine all background + SM WW samples:
  ////////////////////////////////////////////////
  // Create Background +/- Uncertainty histograms:
  std::cout<<"Set error..."<<std::endl;
  TH1D *background_mudijet_backshapeUp = (TH1D*)background->Clone("background_mudijet_backshapeUp");
  TH1D *background_mudijet_backshapeDown = (TH1D*)background->Clone("background_mudijet_backshapeDown");
  background_mudijet_backshapeUp->Reset();
  background_mudijet_backshapeDown->Reset();
//  Double_t change;
  for(int ibin=1;ibin<=pv.NBINS+1;ibin++){
     //////////////////////////////////////////////////////
     // Instead use pre-existing Control Plot uncertainties
     background_mudijet_backshapeUp->SetBinContent(ibin,background->GetBinContent(ibin)+background->GetBinError(ibin));
     background_mudijet_backshapeDown->SetBinContent(ibin,background->GetBinContent(ibin)-background->GetBinError(ibin));
  }

*/
  ////////////////////////
  // Set histogram labels:
  const double BINWIDTH = ((pv.MAXRange-pv.MINRange)/pv.NBINS);
  char tmpc[100];    sprintf(tmpc,"Events / %.1f GeV",BINWIDTH);
  if (pv.slog==1)    sprintf(tmpc,"Events / %.1f",BINWIDTH);
  if (pv.slog==2)    sprintf(tmpc,"Events / %.2f",BINWIDTH);
  if (pv.slog==3)    sprintf(tmpc,"Events / %.0f GeV",BINWIDTH);
  if (pv.slog==6)    sprintf(tmpc,"Events / %.2f rad",BINWIDTH);

  data_obs->SetYTitle(tmpc);
  data_obs->GetXaxis()->SetTitle(pv.xlabel);
  data_obs->GetYaxis()->CenterTitle(true);

  //////////////////////////////////////////////////////////
  // Save Observed Data, Background (+ SM WWA), Background +
  // Uncertainty, Background - Uncertainty, Anomalous Signal
  // histograms to output ROOT file:
  std::cout<<"Save Histograms..."<<std::endl;
  f.cd();
  data_obs->Write();
  //background->Write();
  //  th1fkdata->Write();
  th1sumwjets->Write();
  th1sumstop->Write();
  th1sumdiboson->Write();
  th1sumbkg->Write();
  th1ww->Write();
  th1wz->Write();
  th1zz->Write();

  th1w4jets->Write();
  th1w3jets->Write();
  th1w2jets->Write();
  th1w1jets->Write();

  th1zjets->Write();
  th1Top->Write();
  th1stops->Write();
  th1stopt->Write();
  th1stoptw->Write();

  th1stopps->Write();
  th1stoppt->Write();
  th1stopptw->Write();

  th1madww2jets->Write();
  th1mh126->Write();
/*
  TObjArray *mc = new TObjArray(7);
  mc->Add(th1wajets);
  mc->Add(th1wwa);
  mc->Add(th1wza);
  mc->Add(th1zz);
  mc->Add(th1fkdata);
  mc->Add(th1zajets);
  mc->Add(th1Top);
  mc->Add(th1stops);
  cout<<"wajets"<<"\t"<<th1wajets->Integral()/background->Integral()<<endl;
  cout<<"wwa"<<"\t"<<th1wwa->Integral()/background->Integral()<<endl;
  cout<<"wza"<<"\t"<<th1wza->Integral()/background->Integral()<<endl;
  cout<<"zz"<<"\t"<<th1zz->Integral()/background->Integral()<<endl;
  cout<<"fake"<<"\t"<<th1fkdata->Integral()/background->Integral()<<endl;
  cout<<"Zajets"<<"\t"<<th1zajets->Integral()/background->Integral()<<endl;
  cout<<"Top"<<"\t"<<th1Top->Integral()/background->Integral()<<endl;
  cout<<"single top"<<"\t"<<th1stops->Integral()/background->Integral()<<endl;

  TFractionFitter* fit = new TFractionFitter(data_obs, mc);
//  fit->Constrain(1,0.9*th1wajets->Integral()/background->Integral(),1.1*th1wajets->Integral()/background->Integral()); 
  fit->Constrain(2,0.7178*th1wwa->Integral()/background->Integral(),1.2821*th1wwa->Integral()/background->Integral()); 
  fit->Constrain(3,0.7178*th1wza->Integral()/background->Integral(),1.2821*th1wza->Integral()/background->Integral()); 
  fit->Constrain(4,0.9*th1zz->Integral()/background->Integral(),1.1*th1zz->Integral()/background->Integral()); 
  fit->Constrain(5,0.8*th1fkdata->Integral()/background->Integral(),1.2*th1fkdata->Integral()/background->Integral()); 
  fit->Constrain(6,0.776*th1zajets->Integral()/background->Integral(),1.2245*th1zajets->Integral()/background->Integral()); 
  fit->Constrain(7,0.78725*th1Top->Integral()/background->Integral(),1.2127*th1Top->Integral()/background->Integral()); 
  fit->Constrain(8,0.7*th1stops->Integral()/background->Integral(),1.3*th1stops->Integral()/background->Integral()); 

 fit->SetRangeX(1,10);

  Int_t status = fit->Fit();
  cout << "fit status: " << status << endl;
   if (status == 0) {                       // check on fit status
     result = (TH1F*) fit->GetPlot();
     result->SetMarkerStyle(25);
//     data_obs->Draw("Ep");
     result->Draw();
     result->SetName("result");
     result->Write();
   }

  Double_t value,error;
  fit->GetResult(0, value, error);
   cout<<"wajets  "<<value/th1wajets->Integral()*background->Integral()<<"\t"<<th1wajets->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(1, value, error);
   cout<<"wwa  "<<value/th1wwa->Integral()*background->Integral()<<"\t"<<th1wwa->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(2, value, error);
   cout<<"wza  "<<value/th1wza->Integral()*background->Integral()<<"\t"<<th1wza->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(3, value, error);
   cout<<"zz  "<<value/th1zz->Integral()*background->Integral()<<"\t"<<th1zz->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(4, value, error);
   cout<<"fake  "<<value/th1fkdata->Integral()*background->Integral()<<"\t"<<th1fkdata->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(5, value, error);
   cout<<"zajets  "<<value/th1zajets->Integral()*background->Integral()<<"\t"<<th1zajets->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(6, value, error);
   cout<<"Top  "<<value/th1Top->Integral()*background->Integral()<<"\t"<<th1Top->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;
  fit->GetResult(7, value, error);
   cout<<"Single Top  "<<value/th1stops->Integral()*background->Integral()<<"\t"<<th1stops->Integral()/background->Integral()<<"\t"<<value<<"\t"<<error/value<<endl;

cout<<background->Integral()<<"  "<<result->Integral()<<endl;
  background_mudijet_backshapeUp->Write();
  background_mudijet_backshapeDown->Write();
  signal_a0w_5->Write();
  signal_a0w_3->Write();
  signal_a0w_2->Write();
  signal_a0w_m2->Write();
  signal_a0w_m3->Write();
  signal_a0w_m5->Write();
*/
}// End Main function
