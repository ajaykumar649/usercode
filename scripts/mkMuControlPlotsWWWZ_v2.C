#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TMath.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
//=====================================================================================
// SYNOPSIS:
//   1. Prepare "InData" and "OutDir" directories; e.g., "ln -s . OutDir" to go to current dir
//   2. root [0] .x mkControlPlots.C(0) for electron data, or
//      root [0] .x mkControlPlots.C(1) for muon data
//
// ====================================================================================
// Self Function
// ====================================================================================
/*
 */
void cmspre(double intlumifbinv)
{
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);

	latex.SetTextAlign(31); // align right
	latex.DrawLatex(0.85,0.93,"#sqrt{s} = 8 TeV");
	latex.SetTextAlign(31); // align right
	latex.DrawLatex(0.65,0.93,Form("#int #font[12]{L} dt = %.1f fb^{-1}", (float)intlumifbinv));
	//  latex.DrawLatex(0.65,0.93,Form("2012B", (float)intlumifbinv));

	latex.SetTextAlign(11); // align left
	//  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");
	latex.DrawLatex(0.15,0.96,"CMS preliminary");

}
			struct plotVar_t {
			char* plotvar;
			double MINRange;
			double MAXRange;
			int    NBINS;
			int    slog;
			char* xlabel;
			char* outfile;
			double AMINRange;
			double AMAXRange;
			int    ANBINS;
			int    mva_in;
			int    hplot;
			int  drawleg;
			};
			void mkMuControlPlotsWWWZ_v2(bool domu=true,bool domva=false,
				bool dovbf=false)
			{
				gROOT->ProcessLine(".L tdrstyle.C");
				if (domu)
				{
					const double intLUMI = 19300.;
				}
				else
				{
					const double intLUMI = 19200.;
				}
				const double wpj_kfactor = 1.16;
				double wpj_fitscale =1.01; 
                                double ttbar_fitscale = 1.0;
				double diboson_fitscale = 1.43;	
			//	const double mh126_scale   = 0.0776*(3.0/2.0)* intLUMI/333061; // V
                        //        const double madww2jets_scale   = 0.03907/2.0* intLUMI/414963; // V
				const double WJets_scale   = 36257.2* intLUMI/18353019; // V
				//const double W4Jets_scale  = 214.0 * intLUMI*wpj_kfactor/4369420;
                                const double W4Jets_scale  = 214.0 * intLUMI*wpj_fitscale*wpj_kfactor/12842803;
				const double W3Jets_scale  = 519.0 * intLUMI*wpj_fitscale*wpj_kfactor/15059503;
				const double W2Jets_scale  = 1750.0 * intLUMI*wpj_fitscale*wpj_kfactor/33004921;
                                const double W1Jets_scale  = 5400.0 * intLUMI*wpj_fitscale*wpj_kfactor/19871598;
				const double WW_scale      = 54.838  * intLUMI*diboson_fitscale/9450414; // V
				const double WZ_scale      = 32.3161 * intLUMI*diboson_fitscale/10000267; // V
				//const double ZZ_scale      = 8.059   * intLUMI/9702850; // V
				const double QCD_scale     = 364000000 *    0.00037 * intLUMI/7529312 ; // V
				const double ZJets_scale   = 3503.71  * intLUMI*wpj_fitscale/30209426; // V
				const double ttbar_scale   = 225.197 * intLUMI*ttbar_fitscale/6893735; // V
				const double SToppS_scale  = 1.75776  * intLUMI*ttbar_fitscale/139974; // V anti-top
				const double SToppT_scale  = 30.0042  * intLUMI*ttbar_fitscale/1935066; // V
				const double SToppTW_scale = 11.1773  * intLUMI*ttbar_fitscale/493458; // V
				const double STopS_scale   = 3.89394  * intLUMI*ttbar_fitscale/259960; // top
				const double STopT_scale   = 55.531  * intLUMI*ttbar_fitscale/3758221; //V
				const double STopTW_scale  = 11.1773  * intLUMI*ttbar_fitscale/497657; // V

	const plotVar_t plotvars[] =
       	{
//    plotvar	MINRange  MAXRange  NBINS  slog xlabel outfile AMINRange  AMAXRange ANBINS mva_in hplot drawleg

	{ "Mass2j_PFCor",     20, 400, 18, 3,   "m_{jj} (GeV)",            "mjj_initial",       20,  400, 19, 0, 1, 1 },
	{ "JetPFCor_Pt[0]",   15,   160, 17, 3, "Leading Jet  p_{T}", "jet1_pt",      15,  160, 18, 0, 0, 1 },
	{ "JetPFCor_Eta[0]", -2.6 , 2.6, 12, 1, "Leading Jet  #eta",  "jet1_eta",   -2.6 , 2.6, 13, 0, 0, 0 },
	{ "JetPFCor_Pt[1]",   20,   150, 12, 3, "Jet 2 p_{T}",        "jet2_pt",      20,  150, 13, 0, 0, 0 },
	{ "JetPFCor_Eta[1]", -2.6 , 2.6, 12, 1, "Jet 2 #eta",         "jet2_eta",   -2.6 , 2.6, 13, 0, 0, 0 },
	{ "W_muon_pt",         15, 155, 14, 3,  "Muon p_{T} (GeV)",     "W_muon_pt",       20,  155, 13, 0, 0, 0 },
	{ "W_muon_eta",      -2.5, 2.5, 18, 1,  "Muon #eta",            "W_muon_eta",    -2.5,  2.5, 16, 0, 0, 0 },
	{ "event_met_pfmet",  20, 155, 13, 3,   "pf MET (GeV)",  "event_met_pfmet",     15, 155, 14, 0, 0, 0 },
	{ "W_mt",             20, 140, 12, 3,   "W Transverse Mass (GeV)", "W_mt",      20,  140, 12, 0, 1, 1 },
	{ "JetPFCor_dphiMET[0]", -3.14, 3.14, 16, 1, "#Delta #phi (Leading Jet, MET)", "deltaphi_jetldmet", -3.4, 3.4, 17, 0, 0, 0 },
	{ "JetPFCor_Pt[2]",   15,   160, 12, 3, "Jet 3 p_{T}",        "jet3_pt",      15,  160, 13, 0, 0, 0 },
	{ "(JetPFCor_Eta[0]-JetPFCor_Eta[1])",
	       -3.0, 3.0, 15, 1,  "#Delta #eta (j,j)",    "deltaeta_jj",   -3.2,  3.2, 16, 0, 0, 0 },
	{ "sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",
	            10, 250, 16, 3,   "Dijet System p_{T} (GeV)", "dijet_pt", 10, 250, 17, 0, 0, 0},  
        { "Mass2j_PFCor",     55, 95, 8, 3,   "m_{jj} (GeV)",            "mjj_final",       55,  95, 8, 0, 1, 1 },
	{ "", 0.0,0.0,0,0,"","",0.,0.,0.,0,0,0 },
	};

	//  const char* cutsVector[cutno].c_str() = "1";
	//  double BINWIDTH = ((MAXRange-MINRange)/NBINS);

//       	TCut cutsVector[cutno].c_str()E("effwt*puwt*puwt*(((sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&(abs(JetPFCor_Eta[0])<2.4)&&(abs(JetPFCor_Eta[1])<2.4)&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.))&&(Mass2j_PFCor>60.0)&&(Mass2j_PFCor<90.0))");

//	cutsVector[cutno].c_str()= TCut("effwt*puwt*(((sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&(abs(JetPFCor_Eta[0])<2.4)&&(abs(JetPFCor_Eta[1])<2.4)&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.))&&(Mass2j_PFCor>60.0)&&(Mass2j_PFCor<90.0))");

	std::vector<std::string> cutsVector;
	cutsVector.push_back("effwt*puwt*( vbf_event==0)");
	cutsVector.push_back("effwt*puwt*((vbf_event==0) && ( Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90. ))");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ))");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4))");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.)   )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ))");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.)   )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1)  )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) &&(W_mt>30.)  )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) &&(W_mt>30.) &&(abs(JetPFCor_dphiMET[0])>0.4) )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) &&(W_mt>30.) &&(abs(JetPFCor_dphiMET[0])>0.4) && (JetPFCor_Pt[2]<30.)  )");
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) &&(W_mt>30.) &&(abs(JetPFCor_dphiMET[0])>0.4) && (JetPFCor_Pt[2]<30.)&&(  abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)  )"); 
        cutsVector.push_back("effwt*puwt*((vbf_event==0) && (Mass2j_PFCor > 60. &&  Mass2j_PFCor< 90.  ) && (JetPFCor_Pt[0]>40. ) && ((abs(JetPFCor_Eta[0]))<2.4) && (JetPFCor_Pt[1]>35.) && ((abs(JetPFCor_Eta[1])) <2.4 ) &&(W_muon_pt>25.) && (abs(W_muon_eta)<2.1) && (event_met_pfmet>25) &&(W_mt>30.) &&(abs(JetPFCor_dphiMET[0])>0.4) && (JetPFCor_Pt[2]<30.)&&(  abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5) && ( sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.) )");    


	TFile f("mu1_plotvar2_histo.root", "RECREATE");
	// Get the input trees:
	// Data
	TFile *fin2,*wwShape_file,*wzShape_file,*w4jetsShape_file,*w3jetShape_file,*w2jetShape_file,*w1jetShape_file,*ttbar_file,*qcd_file1,*zjets_file,*stops_file,*stopt_file,*stoptW_file;
	if (domu) {
	fin2            = new TFile("InData/RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root", "read");
	//    H190_file       = new TFile("InData/RD_mu_HWWMH190_CMSSW532_private.root", "READ");
	//madww2jets_file       = new TFile("InData/RD_mu_mh126_CMSSW532.root", "READ");
	//mh126_file       = new TFile("InData/RD_mu_phantom_CMSSW532.root", "READ");
	//    wwShape_file    = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
	wwShape_file    = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
	wzShape_file    = new TFile("InData/RD_mu_WZ_CMSSW532.root", "READ");
//	zzShape_file    = new TFile("InData/RD_mu_ZZ_CMSSW532.root", "READ");
	//if (dovbf)
	//{
//	w4jetsShape_file = new TFile("InData/RD_mu_W4Jets_old_CMSSW532.root","READ");
      w4jetsShape_file = new TFile("InData/RD_mu_W4Jets_CMSSW532_new.root","READ");
//       w4jetsShape_file = new TFile("InData/RD_mu_W4Jets_CMSSW532.root","READ");
	w3jetShape_file = new TFile("InData/RD_mu_W3Jets_CMSSW532.root","READ");
	w2jetShape_file = new TFile("InData/RD_mu_W2Jets_CMSSW532.root","READ");
        w1jetShape_file = new TFile("InData/RD_mu_W1Jets_CMSSW532.root","READ");
	//}
	//else
	//{	
	//w4jetsShape_file = new TFile("InData/RD_mu_WpJ_CMSSW532.root", "READ");
	//}
	ttbar_file      = new TFile("InData/RD_mu_TTbar_CMSSW532.root", "READ");
//    qcd_file1       = new TFile("InData/RDQCD_WmunuJets_DataAll_GoldenJSON_2p1invfb.root", "READ");
//    qcd_file1       = new TFile("InData/RD_mu_QCDMu_CMSSW532.root", "READ");
	zjets_file      = new TFile("InData/RD_mu_ZpJ_CMSSW532.root", "READ");
	stops_file      = new TFile("InData/RD_mu_STopS_T_CMSSW532.root", "READ");
	stopt_file      = new TFile("InData/RD_mu_STopT_T_CMSSW532.root", "READ");
	stoptW_file     = new TFile("InData/RD_mu_STopTW_T_CMSSW532.root", "READ");
	} 
	else 
	{ // electrons
	fin2            = new TFile("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root", "READ");
	//    H190_file       = new TFile("InData/RD_el_HWWMH190_CMSSW428.root", "READ");
	//    H500_file       = new TFile("InData/RD_el_HWWMH500_CMSSW428.root", "READ");
	madww2jets_file       = new TFile("InData/RD_el_VBFHWWMH500_CMSSW532_private.root", "READ");
	wwShape_file    = new TFile("InData/RD_el_WW_CMSSW532.root", "READ");
	wzShape_file    = new TFile("InData/RD_el_WZ_CMSSW532.root", "READ");
	//zzShape_file    = new TFile("InData/RD_el_ZZ_CMSSW532.root", "READ");
	//if (dovbf)
	//{
	  w4jetsShape_file = new TFile("InData/RD_el_W4Jets_CMSSW532_old.root","READ");
	w3jetsShape_file = new TFile("InData/RD_el_W3Jets_CMSSW532.root","READ");
	w2jetsShape_file = new TFile("InData/RD_el_W2Jets_CMSSW532.root","READ");
        w1jetsShape_file = new TFile("InData/RD_el_W1Jets_CMSSW532.root","READ");
//	}
//	else
//	w4jetsShape_file = new TFile("InData/RD_el_WpJ_CMSSW532.root", "READ");
	ttbar_file      = new TFile("InData/RD_el_TTbar_CMSSW532.root", "READ");
//    qcd_file1       = new TFile("InData/RDQCD_WenuJets_Isog0p3NoElMVA_1p6invfb.root", "READ");
//    qcd_file1       = new TFile("InData/RDQCD_WenuJets_DataAll_GoldenJSON_2p1invfb2011.root", "READ");
	zjets_file      = new TFile("InData/RD_el_ZpJ_CMSSW532.root", "READ");
	stops_file      = new TFile("InData/RD_el_STopS_T_CMSSW532.root", "READ");
	stopt_file      = new TFile("InData/RD_el_STopT_T_CMSSW532.root", "READ");
	stoptW_file     = new TFile("InData/RD_el_STopTW_T_CMSSW532.root", "READ");
	}
	TTree* treedata = (TTree*) fin2->Get("WJet");
	double nData = treedata->GetEntries();
	std::cout << "ndata =" << nData <<std::endl;
	//  TTree* treeh190  = (TTree*)       H190_file->Get("WJet");
	//TTree* treeh500  = (TTree*)       H500_file->Get("WJet");
	//	TTree* treevhmadww2jets  = (TTree*)       madww2jets_file->Get("WJet");
	//	TTree* treemh126  = (TTree*)       mh126_file->Get("WJet");
	//  TTree* treeh300  = (TTree*)       H300_file->Get("WJet");
	TTree* treeww    = (TTree*)    wwShape_file->Get("WJet");
	TTree* treewz    = (TTree*)    wzShape_file->Get("WJet");
	//TTree* treezz    = (TTree*)    zzShape_file->Get("WJet");
	TTree* treewj    = (TTree*) w4jetsShape_file->Get("WJet");
	TTree* treew3j    = (TTree*) w3jetShape_file->Get("WJet");
	TTree* treew2j    = (TTree*) w2jetShape_file->Get("WJet");
        TTree* treew1j    = (TTree*) w1jetShape_file->Get("WJet");
	TTree* treettb   = (TTree*)      ttbar_file->Get("WJet");
	//  TTree* treeqcd   = (TTree*)       qcd_file1->Get("WJet");
	TTree* treezj    = (TTree*)      zjets_file->Get("WJet");
	TTree* treests   = (TTree*)      stops_file->Get("WJet");
	TTree* treestt   = (TTree*)      stopt_file->Get("WJet");
	TTree* treestw   = (TTree*)     stoptW_file->Get("WJet");

	//loop over cuts
	for (int cutno=0;cutno<cutsVector.size();cutno++)
	{
// loop over varialbes
	//for (int ivar=0; ; ivar++) 
	//{
	int ivar=cutno;
	plotVar_t pv;
	if (dovbf)
	pv = vbfplotvars[ivar];
	else
	if(domva) 
	pv = higgsplotvars[ivar];
	else
	pv = plotvars[ivar];
	if ( !strlen(pv.plotvar) ) break;
	std::cout << TString(pv.plotvar) << "\t"<<pv.MINRange<<"\t" << pv.MAXRange<<"\t" << pv.NBINS<<"\tTHE CUT " << endl;


	/*    if (domu) {
      if (strstr(pv.plotvar,"el")) continue;
      } else {
      if (strstr(pv.plotvar,"mu")) continue;
      }
	 */
	cout<<"pv.plotvar   "<<pv.plotvar<<endl;

	//    if (dovbf && pv.mva_in)
	//    cutsVector[cutno].c_str() = TCut("effwt*puwt*(hvbf_event && hvbf_wjj_m > 65.0 && hvbf_wjj_m < 95.0)"); // plot only events in the signal region

	const double BINWIDTH = ((pv.MAXRange-pv.MINRange)/pv.NBINS);

	//cout<<"pv.plotvar   "<<pv.plotvar<<endl;

	TH1* th1data  = new TH1D("th1data",  "th1data",  pv.NBINS, pv.MINRange, pv.MAXRange);
	TH1* th1data1 = new TH1D("th1data1", "th1data1", pv.NBINS, pv.MINRange, pv.MAXRange);
	//	TBox *errbox = new TBox(pv.AMINRange,0.95,pv.AMAXRange,1.05);
        TBox *errbox = new TBox(pv.AMINRange,0.974,pv.AMAXRange,1.026);

	errbox->SetFillColor(kYellow);
	treedata->Draw(TString(pv.plotvar)+TString(">>th1data"), cutsVector[cutno].c_str(), "goff");
        treedata->Draw(TString(pv.plotvar)+TString(">>th1data1"), cutsVector[cutno].c_str(), "goff");
	//treedata->Draw(TString(pv.plotvar)+TString(">>th1data1"), cutsVector[cutno].c_str()E, "goff");


	/*	// Get Signal MC

	TH1* th1madww2jets = new TH1D("th1madww2jets", "th1madww2jets", pv.NBINS, pv.MINRange, pv.MAXRange);
	treevhmadww2jets->Draw(TString(pv.plotvar)+TString(">>th1madww2jets"), cutsVector[cutno].c_str(), "goff");
	th1madww2jets->Scale(madww2jets_scale);

	TH1* th1mh126 = new TH1D("th1mh126", "th1mh126", pv.NBINS, pv.MINRange, pv.MAXRange);
	treemh126->Draw(TString(pv.plotvar)+TString(">>th1mh126"), cutsVector[cutno].c_str(), "goff");
	th1mh126->Scale(mh126_scale);

*/		// Get WW/WZ/ZZ

		TH1* th1ww = new TH1D("th1ww", "th1ww", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1wz = new TH1D("th1wz", "th1wz", pv.NBINS, pv.MINRange, pv.MAXRange);
		//TH1* th1zz = new TH1D("th1zz", "th1zz", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1wwSQ = new TH1D("th1wwSQ", "th1wwSQ", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1wzSQ = new TH1D("th1wzSQ", "th1wzSQ", pv.NBINS, pv.MINRange, pv.MAXRange);
		//TH1* th1zzSQ = new TH1D("th1zzSQ", "th1zzSQ", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1ww->Sumw2();
		th1wz->Sumw2();
		//th1zz->Sumw2();

		treeww->Draw(TString(pv.plotvar)+TString(">>th1ww"), cutsVector[cutno].c_str(), "goff");
		treewz->Draw(TString(pv.plotvar)+TString(">>th1wz"), cutsVector[cutno].c_str(), "goff");
		//treezz->Draw(TString(pv.plotvar)+TString(">>th1zz"), cutsVector[cutno].c_str(), "goff");
		treeww->Draw(TString(pv.plotvar)+TString(">>th1wwSQ"), cutsVector[cutno].c_str(), "goff");
		treewz->Draw(TString(pv.plotvar)+TString(">>th1wzSQ"), cutsVector[cutno].c_str(), "goff");
		//treezz->Draw(TString(pv.plotvar)+TString(">>th1zzSQ"), cutsVector[cutno].c_str()E, "goff");
		for(int hi=1;hi<=pv.NBINS;hi++)th1ww->SetBinError(hi,sqrt(th1wwSQ->GetBinContent(hi)));
		for(int hi=1;hi<=pv.NBINS;hi++)th1wz->SetBinError(hi,sqrt(th1wzSQ->GetBinContent(hi)));
		//for(int hi=1;hi<=pv.NBINS;hi++)th1zz->SetBinError(hi,sqrt(th1zzSQ->GetBinContent(hi)));
		// Get WJets
		TH1* th1wjets  = new TH1D("th1wjets",  "th1wjets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
		TH1* th1w3jets  = new TH1D("th1w3jets",  "th1w3jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
		TH1* th1w2jets  = new TH1D("th1w2jets",  "th1w2jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);
                TH1* th1w1jets  = new TH1D("th1w1jets",  "th1w1jets",  pv.NBINS ,pv.MINRange,pv.MAXRange);

		TH1* th1wjetsSQ  = new TH1D("th1wjetsSQ",  "th1wjetsSQ",  pv.NBINS ,pv.MINRange,pv.MAXRange);
		TH1* th1w3jetsSQ  = new TH1D("th1w3jetsSQ",  "th1w3jetsSQ",  pv.NBINS ,pv.MINRange,pv.MAXRange);
		TH1* th1w2jetsSQ  = new TH1D("th1w2jetsSQ",  "th1w2jetsSQ",  pv.NBINS ,pv.MINRange,pv.MAXRange);
                TH1* th1w1jetsSQ  = new TH1D("th1w1jetsSQ",  "th1w1jetsSQ",  pv.NBINS ,pv.MINRange,pv.MAXRange);

                th1wjets->Sumw2();
                th1w3jets->Sumw2();
                th1w2jets->Sumw2();
                th1w1jets->Sumw2();

		treewj->Draw(TString(pv.plotvar)+TString(">>th1wjets"), cutsVector[cutno].c_str(), "goff");
		treew3j->Draw(TString(pv.plotvar)+TString(">>th1w3jets"), cutsVector[cutno].c_str(), "goff");
		treew2j->Draw(TString(pv.plotvar)+TString(">>th1w2jets"), cutsVector[cutno].c_str(), "goff");
                treew1j->Draw(TString(pv.plotvar)+TString(">>th1w1jets"), cutsVector[cutno].c_str(), "goff");

		treewj->Draw(TString(pv.plotvar)+TString(">>th1wjetsSQ"), cutsVector[cutno].c_str(), "goff");
		treew3j->Draw(TString(pv.plotvar)+TString(">>th1w3jetsSQ"), cutsVector[cutno].c_str(), "goff");
		treew2j->Draw(TString(pv.plotvar)+TString(">>th1w2jetsSQ"), cutsVector[cutno].c_str(), "goff");
                treew1j->Draw(TString(pv.plotvar)+TString(">>th1w1jetsSQ"), cutsVector[cutno].c_str(), "goff");

		for(int hi=1;hi<=pv.NBINS;hi++) th1wjets->SetBinError(hi,sqrt(th1wjetsSQ->GetBinContent(hi)));
		for(int hi=1;hi<=pv.NBINS;hi++) th1w3jets->SetBinError(hi,sqrt(th1w3jetsSQ->GetBinContent(hi)));
		for(int hi=1;hi<=pv.NBINS;hi++) th1w2jets->SetBinError(hi,sqrt(th1w2jetsSQ->GetBinContent(hi)));
                for(int hi=1;hi<=pv.NBINS;hi++) th1w1jets->SetBinError(hi,sqrt(th1w1jetsSQ->GetBinContent(hi)));

		// Get ttbar

		TH1* th1Top = new TH1D("th1Top", "th1Top", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1TopSQ = new TH1D("th1TopSQ", "th1TopSQ", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1Top->Sumw2();
		// cross section: 157.5 pb, events_gen = 3701947 (These are summer11 TTJets sample

		treettb->Draw(TString(pv.plotvar)+TString(">>th1Top"), cutsVector[cutno].c_str(), "goff");
		treettb->Draw(TString(pv.plotvar)+TString(">>th1TopSQ"), cutsVector[cutno].c_str(), "goff");
		for(int hi=1;hi<=pv.NBINS;hi++)th1Top->SetBinError(hi,sqrt(th1TopSQ->GetBinContent(hi)));
		/*
		// Get QCD
		TH1* th1qcd = new TH1D("th1qcd", "th1qcd", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1qcd->Sumw2();
		treeqcd->Draw(TString(pv.plotvar)+TString(">>th1qcd"), cutsVector[cutno].c_str(), "goff");
		int n2 = treeqcd->Draw(TString(pv.plotvar),  cutsVector[cutno].c_str()2, "goff");
		int n3 = treeqcd->Draw(TString(pv.plotvar),  cutsVector[cutno].c_str()3 , "goff");
		std::cout << "got qcd " << " n2 " << n2 <<  " n3  " << n3 <<std::endl; 
		 */
		// Get Z+Jets

		TH1* th1zjets = new TH1D("th1zjets", "th1zjets", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1zjetsSQ = new TH1D("th1zjetsSQ", "th1zjetsSQ", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1zjets->Sumw2();
		treezj->Draw(TString(pv.plotvar)+TString(">>th1zjets"), cutsVector[cutno].c_str(), "goff");
		treezj->Draw(TString(pv.plotvar)+TString(">>th1zjetsSQ"), cutsVector[cutno].c_str(), "goff");
		for(int hi=1;hi<=pv.NBINS;hi++)th1zjets->SetBinError(hi,sqrt(th1zjetsSQ->GetBinContent(hi)));
                //treeww->Draw(TString("numPFCorJets")+TString(">>tmpHist"), cutsVector[iExtraStep].c_str(), "goff");

		// Get Single top

		TH1* th1stops = new TH1D("th1stops", "th1stops", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1stopt = new TH1D("th1stopt", "th1stopt", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1stoptw = new TH1D("th1stoptw", "th1stoptw", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1stops->Sumw2();
		th1stopt->Sumw2();
		th1stoptw->Sumw2();

		treests->Draw(TString(pv.plotvar)+TString(">>th1stops"), cutsVector[cutno].c_str(), "goff");
		treestt->Draw(TString(pv.plotvar)+TString(">>th1stopt"), cutsVector[cutno].c_str(), "goff");
		treestw->Draw(TString(pv.plotvar)+TString(">>th1stoptw"), cutsVector[cutno].c_str(), "goff");

		if (domu)
		{
		TFile* stopps_file =  new TFile("InData/RD_mu_STopS_Tbar_CMSSW532.root", "READ");
		TTree* tree64 = (TTree*) stopps_file->Get("WJet");
		TFile* stoppt_file =  new TFile("InData/RD_mu_STopT_Tbar_CMSSW532.root", "READ");
		TTree* tree65 = (TTree*) stoppt_file->Get("WJet");
		TFile* stopptW_file =  new TFile("InData/RD_mu_STopTW_Tbar_CMSSW532.root", "READ");
		TTree* tree66 = (TTree*) stopptW_file->Get("WJet");
		}
		else
		{
		TFile* stopps_file =  new TFile("InData/RD_el_STopS_Tbar_CMSSW532.root", "READ");
		TTree* tree64 = (TTree*) stopps_file->Get("WJet");
		TFile* stoppt_file =  new TFile("InData/RD_el_STopT_Tbar_CMSSW532.root", "READ");
		TTree* tree65 = (TTree*) stoppt_file->Get("WJet");
		TFile* stopptW_file =  new TFile("InData/RD_el_STopTW_Tbar_CMSSW532.root", "READ");
		TTree* tree66 = (TTree*) stopptW_file->Get("WJet");
		}
		TH1* th1stopps = new TH1D("th1stopps", "th1stopps", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1stoppt = new TH1D("th1stoppt", "th1stoppt", pv.NBINS, pv.MINRange, pv.MAXRange);
		TH1* th1stopptw = new TH1D("th1stopptw", "th1stopptw", pv.NBINS, pv.MINRange, pv.MAXRange);
		th1stopps->Sumw2();
		th1stoppt->Sumw2();
		th1stopptw->Sumw2();
		tree64->Draw(TString(pv.plotvar)+TString(">>th1stopps"), cutsVector[cutno].c_str(), "goff");
		tree65->Draw(TString(pv.plotvar)+TString(">>th1stoppt"), cutsVector[cutno].c_str(), "goff");
		tree66->Draw(TString(pv.plotvar)+TString(">>th1stopptw"), cutsVector[cutno].c_str(), "goff");

		// Setup the canvas
		//    gROOT->ProcessLine(".L tdrstyle.C");
		setTDRStyle();
		tdrStyle->SetErrorX(0.5);
		tdrStyle->SetPadRightMargin(0.05);

		tdrStyle->SetLegendBorderSize(0);

		th1data->Sumw2();

		TCanvas* c1 = new TCanvas(pv.plotvar,pv.plotvar,10,10, 800, 800);
		TPad *d1, *d2;

		c1->Divide(1,2,0,0);
		d1 = (TPad*)c1->GetPad(1);
		d1->SetPad(0.01,0.30,0.95,0.99);
		d2 = (TPad*)c1->GetPad(2);
		d2->SetPad(0.01,0.02,0.95,0.30);

		// Setup the stack, scale the histos

		THStack* hs = new THStack("hs","MC contribution");
		th1Top->Scale(ttbar_scale);
		th1Top->SetFillColor(kGreen+2);
		th1Top->SetLineWidth(0);

		th1stops->Scale(STopS_scale);
		th1stops->SetFillColor(7);
		th1stops->SetLineWidth(2);

		th1stopt->Scale(STopT_scale);
		th1stopt->SetFillColor(13);
		th1stopt->SetLineWidth(2);

		th1stoptw->Scale(STopTW_scale);
		th1stoptw->SetFillColor(9);
		th1stoptw->SetLineWidth(2);

		th1stopps->Scale(SToppS_scale);
		th1stopps->SetFillColor(7);
		th1stopps->SetLineWidth(2);

		th1stoppt->Scale(SToppT_scale);
		th1stoppt->SetFillColor(13);
		th1stoppt->SetLineWidth(2);

		th1stopptw->Scale(SToppTW_scale);
		th1stopptw->SetFillColor(9);
		th1stopptw->SetLineWidth(2);

		th1wjets->Scale(W4Jets_scale);
		th1wjets->SetFillColor(kRed);
		th1wjets->SetLineWidth(0);

		th1w3jets->Scale(W3Jets_scale);
		th1w3jets->SetFillColor(kRed+1);
		th1w3jets->SetLineWidth(0);

		th1w2jets->Scale(W2Jets_scale);
		th1w2jets->SetFillColor(kRed+2);
		th1w2jets->SetLineWidth(0);

                th1w1jets->Scale(W1Jets_scale);
                th1w1jets->SetFillColor(kRed+3);
                th1w1jets->SetLineWidth(0);

		th1ww->Scale(WW_scale);
		th1ww->SetFillColor(kAzure+8);
		th1ww->SetLineWidth(0);

		th1wz->Scale(WZ_scale);
		th1wz->SetFillColor(11);
		th1wz->SetLineWidth(0);
/*
		th1zz->Scale(ZZ_scale);
		th1zz->SetFillColor(11);
		th1zz->SetLineWidth(0);
*/
		// th1qcd->Scale(QCD_scale);

		//    th1qcd->SetFillColor(kGray+1);
		//    th1qcd->SetLineColor(kGray+1);
		//    th1qcd->SetLineWidth(0);
		th1zjets->Scale(ZJets_scale);
		th1zjets->SetFillColor(kYellow);
		th1zjets->SetLineWidth(0);
		//    std::cout << " qcd " << th1qcd->Integral()   << std::endl;
//		std::cout << "mad ww2jets   "   <<th1madww2jets->Integral()   << std::endl;
//		std::cout << "mh 126 "   <<th1mh126->Integral()   << std::endl;
		std::cout << "ttbar "   << th1Top->Integral()   << std::endl;
		std::cout << "s top s chann "   << th1stops->Integral()   << std::endl;
		std::cout << "s top t chann "   << th1stopt->Integral()   << std::endl;
		std::cout << "s top tw chann "   << th1stoptw->Integral()   << std::endl;
		std::cout << "s top s chann bar "   << th1stopps->Integral()    << std::endl;
		std::cout << "s top t chann bar "   << th1stoppt->Integral()    << std::endl;
		std::cout << "s top tw  chann bar "   << th1stopptw->Integral()    << std::endl;
		std::cout << "ww "   << th1ww->Integral()    << std::endl;
		std::cout << "wz "   << th1wz->Integral()    << std::endl;
		//std::cout << "zz "   << th1zz->Integral()    << std::endl;
		std::cout << "W4jets "   << th1wjets->Integral()    << std::endl;
		std::cout << "W3jets "   << th1w3jets->Integral()    << std::endl;
		std::cout << "W2jets "    << th1w2jets->Integral() << std::endl;
                std::cout << "W1jets "    << th1w1jets->Integral() << std::endl;
		std::cout << "Zjets "    << th1zjets->Integral() << std::endl;

		double den_qcd = 
		th1Top->Integral()+
		th1stops->Integral()+
		th1stopt->Integral()+
		th1stoptw->Integral()+
		th1stopps->Integral()+
		th1stoppt->Integral()+
		th1stopptw->Integral()+
		th1wjets->Integral()+
		th1w3jets->Integral()+
		th1w2jets->Integral()+
                th1w1jets->Integral()+
		th1ww->Integral()+
		th1wz->Integral()+
		//th1zz->Integral()+
		th1zjets->Integral();
		/*
	   double qcd_scale;
	   if (domu)
	   qcd_scale = (n2*0.002 + n3*0.000) / (n2+n3) ;//muon
	   else
	   qcd_scale = (n2*0.0637 + n3*0.02) / (n2+n3); //electron
	 */
	//    std::cout << " qcd_scale  " << qcd_scale <<std::endl;
	//    th1qcd->Scale(qcd_scale*den_qcd/th1qcd->Integral()); 
		double den = 
		th1Top->Integral()+
		th1stops->Integral()+
		th1stopt->Integral()+
		th1stoptw->Integral()+
		th1stopps->Integral()+
		th1stoppt->Integral()+
		th1stopptw->Integral()+
		th1wjets->Integral()+
		th1w3jets->Integral()+
		th1w2jets->Integral()+
                th1w1jets->Integral()+
		th1ww->Integral()+
		th1wz->Integral()+
		//th1zz->Integral()+
		th1zjets->Integral();
		//      th1qcd->Integral();
                std::cout <<"  " <<std::endl<<std::endl;
                std::cout << "den_qcd = " <<den_qcd <<std::endl;
	//	std::cout << "den = " <<den <<std::endl;
		std::cout <<" data " <<  th1data->Integral() << std::endl;
                std::cout <<" scale factor " <<  th1data->Integral()/den_qcd  << std::endl;
  /*              std::cout << "S/B 1 = " <<(th1mh126->Integral()/den) <<std::endl;
                std::cout << "S/B 2 = " <<(th1madww2jets->Integral()/den) <<std::endl;
                std::cout << "S/sqrt(B) 1 = " <<(th1mh126->Integral()/sqrt(den)) <<std::endl;
                std::cout << "S/sqrt(B )2 = " <<(th1madww2jets->Integral()/sqrt(den)) <<std::endl;
*/
                std::cout <<"  " <<std::endl<<std::endl;
//    th1qcd->Scale   (th1data->Integral()/den); std::cout <<"qcd "   << th1qcd->Integral()   << std::endl;
      th1Top->Scale   (th1data->Integral()/den); std::cout <<"tt "    << th1Top->Integral()   << std::endl;
      th1stops->Scale (th1data->Integral()/den); std::cout <<"stops " << th1stops->Integral() << std::endl;
      th1stopt->Scale (th1data->Integral()/den); std::cout <<"stopt " << th1stopt->Integral() << std::endl;
      th1stoptw->Scale(th1data->Integral()/den); std::cout <<"stoptw "<< th1stoptw->Integral()<< std::endl;
    th1stopps->Scale (th1data->Integral()/den); std::cout <<"stops bar " << th1stopps->Integral() << std::endl;
   th1stoppt->Scale (th1data->Integral()/den); std::cout <<"stopt bar " << th1stoppt->Integral() << std::endl;
   th1stopptw->Scale(th1data->Integral()/den); std::cout <<"stoptw bar "<< th1stopptw->Integral()<< std::endl;
     th1wjets->Scale (th1data->Integral()/den); std::cout <<"wjets " << th1wjets->Integral() << std::endl;
     th1w3jets->Scale (th1data->Integral()/den); std::cout <<"w3jets " << th1w3jets->Integral() << std::endl;
      th1w2jets->Scale (th1data->Integral()/den); std::cout <<"w2jets " << th1w2jets->Integral() << std::endl;
      th1w1jets->Scale (th1data->Integral()/den); std::cout <<"w1jets " << th1w1jets->Integral() << std::endl;
      th1ww->Scale    (th1data->Integral()/den); std::cout <<"ww "    << th1ww->Integral()    << std::endl;
      th1wz->Scale    (th1data->Integral()/den); std::cout << "wz "   << th1wz->Integral()    << std::endl;
//      th1zz->Scale    (th1data->Integral()/den); std::cout << "zz "   << th1zz->Integral()    << std::endl;
     th1zjets->Scale (th1data->Integral()/den); std::cout << "zjets "    << th1zjets->Integral() << std::endl;
     
		//th1madww2jets->Scale(th1data->Integral()/den);
		//th1mh126->Scale(th1data->Integral()/den);
		//    th1H300->Scale(th1data->Integral()/den);
		//    th1H190->Scale(th1data->Integral()/den);
		//  cout<<"(th1data->Integral()/den) = "<< (th1data->Integral()/den) <<endl;
	double den2 =
	th1Top->Integral()+
	th1stops->Integral()+
	th1stopt->Integral()+
	th1stoptw->Integral()+
	th1stopps->Integral()+
	th1stoppt->Integral()+
	th1stopptw->Integral()+
	th1wjets->Integral()+
	th1w3jets->Integral()+
	th1w2jets->Integral()+
        th1w1jets->Integral()+
	th1ww->Integral()+
	th1wz->Integral()+
	//th1zz->Integral()+
	th1zjets->Integral();
	//      th1qcd->Integral();
//	std::cout << "den2 " << den2 << std::endl;
	th1Top->Add(th1stopptw,1);
	th1Top->Add(th1stoppt,1);
	th1Top->Add(th1stopps,1);
	th1Top->Add(th1stoptw,1);
	th1Top->Add(th1stopt,1);
	th1Top->Add(th1stops,1);
	th1ww->Add(th1wz,1);
	//th1ww->Add(th1zz,1);
	th1wjets->Add(th1w3jets,1);
	th1wjets->Add(th1w2jets,1);
        th1wjets->Add(th1w1jets,1);

	// Sum all the backgrounds
	TH1D *th1tot = (TH1D*)th1wjets->Clone("th1tot");
	th1tot->Reset();
	th1tot->Add(th1ww,1);
	//    th1tot->Add(th1qcd,1);
	th1tot->Add(th1Top,1);
	th1tot->Add(th1wjets,1);
	th1tot->Add(th1zjets,1);

	TH1D* th1totClone = ( TH1D*) th1tot->Clone("th1totClone");
	th1totClone->SetMarkerStyle(0);
	th1totClone->SetFillStyle(3003);
	th1totClone->SetFillColor(11);
	th1totClone->SetLineColor(0);
	double binErr(0.0);
	for(int i=0; i<th1totClone->GetNbinsX(); ++i) {
	binErr = sqrt(
	(th1ww->GetBinError(i))**2 +
	// (th1qcd->GetBinError(i))**2 +
	(th1Top->GetBinError(i))**2 +
	(th1wjets->GetBinError(i))**2 +
	(th1zjets->GetBinError(i))**2);
	th1totClone->SetBinError(i, binErr);
	}
	// Compose the stack
	hs->Add(th1zjets);
	d1->cd();
	gPad->SetBottomMargin(0.0);
	gPad->SetTopMargin(0.1);
	gPad->SetRightMargin(0.05);
	gPad->SetLeftMargin(0.14);
	//    hs->Add(th1qcd);
	hs->Add(th1Top);
	hs->Add(th1ww);
	hs->Add(th1wjets);

// Set up the legend

float  legX0=0.65, legX1=0.99, legY0=0.4, legY1=0.88;
// float  legX0=0.35, legX1=0.85, legY0=0.4, legY1=0.88;
// float  legX0=0.18, legX1=0.52, legY0=0.4, legY1=0.88;
TLegend * Leg = new TLegend( legX0, legY0, legX1, legY1);
Leg->SetFillColor(0);
Leg->SetFillStyle(0);
Leg->SetTextSize(0.04);
if (domu)
Leg->AddEntry(th1data,  "Muon Data",  "PLE");
else
Leg->AddEntry(th1data,  "Electron Data",  "PLE");
Leg->AddEntry(th1wjets,  "W+Jets",  "f");
Leg->AddEntry(th1ww,  "Diboson ",  "f");
Leg->AddEntry(th1Top,  "top",  "f");
//    Leg->AddEntry(th1qcd,  "QCD",  "f");
Leg->AddEntry(th1zjets,  "Z+Jets",  "f");
Leg->AddEntry(th1tot,  "MC Uncertainty",  "f");
//Leg->AddEntry(th1madww2jets,  "mad VBFWW X100  ",  "L");
//Leg->AddEntry(th1mh126,  "phan VBFWW X100 ",  "L");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")&&!strstr(pv.plotvar,"mva2j500mu")) Leg->AddEntry(th1H300,  "H (300) x 200",  "L");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j500mu")) Leg->AddEntry(th1H190,  "H (190) x 100",  "L");
Leg->SetFillColor(0);
/*th1madww2jets->SetLineColor(kBlue);
th1madww2jets->SetLineWidth(3);
th1madww2jets->Scale(100);
//    th1madww2jets->Draw("sames");
th1mh126->SetLineColor(kBlack);
th1mh126->SetLineWidth(3);
th1mh126->SetLineStyle(2);
th1mh126->Scale(100);
*/
TH1* th1totempty = new TH1D("th1totempty", "th1totempty", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
th1data->SetMarkerStyle(20);
th1data->SetMarkerSize(1.25);
th1data->SetLineWidth(2);

th1tot->SetFillStyle(3001);
th1tot->SetFillColor(1);
th1tot->SetLineColor(1);
th1tot->SetMarkerStyle(0);

char tmpc[100];    sprintf(tmpc,"Events / %.1f GeV",BINWIDTH);
if (pv.slog==1)    sprintf(tmpc,"Events/ %.1f",BINWIDTH);
if (pv.slog==2)    sprintf(tmpc,"Events/ %.2f",BINWIDTH);
if (pv.slog==3)    sprintf(tmpc,"Events/ %.0f GeV",BINWIDTH);
if (pv.slog==6)    sprintf(tmpc,"Events/ %.1f rad",BINWIDTH);
th1totempty->SetYTitle(tmpc);
//  th1totempty->GetYaxis()->SetTitleSize(0.1);
th1totempty->GetYaxis()->SetTitleOffset(1.2);
th1totempty->GetYaxis()->SetLabelOffset(0.01);
//  th1totempty->GetYaxis()->CenterTitle(true);
th1totempty->GetYaxis()->SetLabelSize(0.04);
// th1totClone->Draw("e3");   

th1tot->SetMinimum(0.0);
int maxbin = th1data->GetMaximumBin();
float maxval = th1data->GetBinContent(maxbin);
std::cout << "maxval " <<maxval <<std::endl;
//    th1totempty->SetMaximum(2.5*maxval);
th1totempty->SetMaximum(1.2*maxval);
th1totempty->SetMinimum(0.0);
if(pv.slog==3) th1totempty->SetMaximum(1.6*maxval);
th1data->SetMinimum(0.0);
// Draw it all
th1totempty->Draw();
//th1tot->Draw("e2same");
th1data->Draw("esame");
hs->Draw("samehist");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")) 
//th1madww2jets->Draw("same");
//th1mh126->Draw("same");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")&&!strstr(pv.plotvar,"mva2j500mu")) th1H300->Draw("same");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j500mu")) th1H190->Draw("same");
th1tot->Draw("e2same");

th1data->Draw("esame");
cmspre(intLUMI/1000.0);    
if (pv.drawleg ==1)  Leg->Draw();  
// th1data->Draw("Axissame");
gPad->RedrawAxis();
d2->cd();

TH1F    * hhratio    = (TH1F*) th1data->Clone("hhratio")  ;
hhratio->Sumw2();
hhratio->SetStats(0);

gPad->SetLeftMargin(0.14);
gPad->SetTopMargin(0);
gPad->SetRightMargin(0.05);
gPad->SetFrameBorderSize(0);
gPad->SetBottomMargin(0.3);
gPad->SetTickx();

hhratio->SetMarkerSize(1.25);
//  hhratio->GetYaxis()->SetRangeUser(0.48,1.52);
hhratio->GetYaxis()->SetRangeUser(0.3,1.7);
hhratio->GetXaxis()->SetTitle(pv.xlabel);
hhratio->GetXaxis()->SetTitleOffset(0.9);
hhratio->GetXaxis()->SetTitleSize(0.15);
hhratio->GetXaxis()->SetLabelSize(0.15);
hhratio->SetYTitle("Ratio Data/MC");
hhratio->GetYaxis()->SetTitleSize(0.1);
hhratio->GetYaxis()->SetTitleOffset(0.5);
hhratio->GetYaxis()->CenterTitle(true);
hhratio->GetYaxis()->SetLabelSize(0.1);
std::cout << hhratio->GetNbinsX() << std::endl;
std::cout << th1tot->GetNbinsX() << std::endl;
hhratio->Divide(th1tot);
double binError(0.0), mcbinentry(0.0), mcerror(0.0);
for(int i=0; i<hhratio->GetNbinsX(); ++i) {
binError = hhratio->GetBinError(i);
mcerror = th1tot->GetBinError(i);
mcbinentry = th1tot->GetBinContent(i);
if(mcbinentry>0.) mcerror /= mcbinentry;
else mcerror = 0.0;
binError = sqrt(binError**2 + mcerror**2);
hhratio->SetBinError(i, binError);
}
TH1D *th1emptyclone = new TH1D("th1emptyclone", "th1emptyclone", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
th1emptyclone->GetYaxis()->SetRangeUser(0.6,1.3999);
th1emptyclone->GetXaxis()->SetTitle(pv.xlabel);
th1emptyclone->GetXaxis()->SetTitleOffset(0.9);
th1emptyclone->GetXaxis()->SetTitleSize(0.15);
th1emptyclone->GetXaxis()->SetLabelSize(0.15);
th1emptyclone->SetYTitle("Ratio Data/MC");
th1emptyclone->GetYaxis()->SetTitleSize(0.1);
th1emptyclone->GetXaxis()->SetNdivisions(505);
th1emptyclone->GetYaxis()->SetNdivisions(505);
th1emptyclone->GetYaxis()->SetTitleOffset(0.5);
th1emptyclone->GetYaxis()->CenterTitle(true);
th1emptyclone->GetYaxis()->SetLabelSize(0.1);
th1emptyclone->Draw();
errbox->Draw();
hhratio->Draw("esame");
TLine *line; line = new TLine(pv.AMINRange,1.0,pv.AMAXRange,1.0);
line->SetLineStyle(1);
line->SetLineWidth(1);
line->SetLineColor(1);
line->Draw();

TString outfile = TString("OutDir/")+(domu?TString("mu_"):TString("el_"))+TString(pv.outfile);
c1->Print(outfile+".png");
c1->Print(outfile+".C");
//gPad->WaitPrimitive();
c1->Modified();
c1->Update();
//c1->SaveAs(outfile+".pdf"); 
//} // var loop
}// cutno loop
// f.Write();
}

