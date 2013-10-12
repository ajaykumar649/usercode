// c++ -o read_007 `root-config --glibs --cflags` -lm read_007.cpp
//  /home/ajay/Work/LPC/Morion13/VBF/RD_mu_TTbar_CMSSW532.root
//c++ -o test0010_method5 `root-config --glibs --cflags` -lm test0010_method5.cpp
/*
third method
A) consider only jets with |eta|<2.5, chose the 2 with mjj more similar
to mW
as the ones from W
B) consider ALL the remaining jets, search the pair of jets that satisfy
the default tag selection,
if more then one then keep the one with largest Mjj
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>
#include "TROOT.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"

#include "conversions.h"
#include "vbf.h"
#include "vbf_histos.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
//#include "lorentzVector.h"
struct dist_sort: public std::binary_function<pair <lorentzVector *, lorentzVector *>, pair <lorentzVector *, lorentzVector *>, bool>
{
  bool operator() (const pair <lorentzVector *, lorentzVector *> & x, const pair <lorentzVector *, lorentzVector *> & y)
    {
      return  (deltaR2 (x.first->Eta (), x.first->Phi (), x.second->Eta (), x.second->Phi()) < 
               deltaR2 (y.first->Eta (), y.first->Phi (), y.second->Eta (), y.second->Phi())) ;
    }
} ;


struct local_TLVP_ptSort: public std::binary_function<const lorentzVector&,const lorentzVector& , bool>
{
  bool operator() (const lorentzVector& x, const lorentzVector& y)
    {
      return x.Pt () < y.Pt () ;

    }
} ;

float my_deltaR_function(float eta1, float phi1, float eta2, float phi2)
{
float delta_phi(float phi1, float phi2);
        const float PI=2.0*acos(0.);
        const float TWOPI=2.0*PI;
        float deta = eta1-eta2;
        float dphi = delta_phi(phi1,phi2);
        return sqrt(deta*deta + dphi*dphi);
}

float delta_phi(float phi1, float phi2)
{
        const float PI=2.0*acos(0.);
        const float TWOPI=2.0*PI;
        float PHI=fabs(phi1-phi2);
        return (PHI<=PI)? PHI : TWOPI-PHI;
}

int main (int argc, char** argv)
{
  if (argc != 2)
    {
      std::cout << ">>> Usage:   " << argv[0] << "   treeFile.root maxNumber deltaR" << std::endl;
      return -1;
    }

//  int maxNumber = atoi (argv[1]) ;
//  double maxDR = atof (argv[2]) ;
  double maxDR = atof (argv[1]) ;
	char method5;

   gROOT->Reset();
 TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RD_mu_TTbar_CMSSW532.root");
   if (!f) {
     f = new TFile("RD_mu_TTbar_CMSSW532.root");
   }
	TTree * WJet;   
 f->GetObject("WJet",WJet);

   Int_t           numPFCorJets;
   Int_t           numPFCorJetBTags;
   Float_t         JetPFCor_Px[8];
   Float_t         JetPFCor_Py[8];
   Float_t         JetPFCor_Pz[8];
   Float_t         JetPFCor_E[8];
   Float_t         JetPFCorVBFTag_Pt[8];
   Float_t         JetPFCorVBFTag_Eta[8];
   Float_t         JetPFCorVBFTag_Phi[8];
   Float_t         JetPFCorVBFTag_E[8];

   Float_t         W_Parton_px[2];
   Float_t         W_Parton_py[2];
   Float_t         W_Parton_pz[2];
   Float_t         W_Parton_E[2];
   Float_t         W_Parton_pt[2];

   Float_t         W_TagQuark_px[2];
   Float_t         W_TagQuark_py[2];
   Float_t         W_TagQuark_pz[2];
   Float_t         W_TagQuark_E[2];
   Float_t         W_TagQuark_pt[2];

  Float_t         W_H_mass_gen;
  Float_t         W_H_px_gen;
  Float_t         W_H_py_gen;
  Float_t         W_H_pz_gen;
  Float_t         W_H_e_gen;

   Float_t         W_muon_pt;

   Float_t         W_muon_px;
   Float_t         W_muon_py;
   Float_t         W_muon_pz;
   Float_t         W_muon_e;
   Float_t         event_met_pfmet;
   Float_t         event_met_pfmetPhi;
   Float_t         W_muon_dz000;
   Float_t         W_muon_dzPV;
   Float_t         W_muon_et;
   Float_t         W_muon_eta;
   Float_t         W_muon_theta;
   Float_t         W_muon_phi;




  TBranch        *b_numPFCorJets ;   //!
  TBranch        *b_numPFCorJetBTags ;   //!
  TBranch        *b_JetPFCor_Px ;   //!
  TBranch        *b_JetPFCor_Py ;   //!
  TBranch        *b_JetPFCor_Pz ;   //!
  TBranch        *b_JetPFCor_E ;   //!
  TBranch        *b_JetPFCorVBFTag_Pt;   //!
  TBranch        *b_JetPFCorVBFTag_Eta;   //!
  TBranch        *b_JetPFCorVBFTag_Phi;   //!
  TBranch        *b_JetPFCorVBFTag_E;   //!

   TBranch        *b_W_muon_pt;   //!
   TBranch        *b_W_muon_eta;   //!

   TBranch        *b_W_muon_px;   //!
   TBranch        *b_W_muon_py;   //!
   TBranch        *b_W_muon_pz;   //!
   TBranch        *b_W_muon_e;   //!
   TBranch        *b_event_met_pfmet;   //!
   TBranch        *b_event_met_pfmetPhi;   //!
   TBranch        *b_W_muon_dz000;   //!
   TBranch        *b_W_muon_dzPV;   //!


 // TBranch        *b_GroomedJet_numberbjets ;
  TBranch        *b_W_H_mass_gen ;
  TBranch        *b_W_H_px_gen ;
  TBranch        *b_W_H_py_gen ;
  TBranch        *b_W_H_pz_gen ;
  TBranch        *b_W_H_e_gen ;
//  TBranch        *b_GroomedJet_CA[6]_tau2tau1;   //!

   TBranch        *b_W_Parton_px;   //!
   TBranch        *b_W_Parton_py;   //!
   TBranch        *b_W_Parton_pz;   //!
   TBranch        *b_W_Parton_E;   //!
   TBranch        *b_W_TagQuark_px;   //!
   TBranch        *b_W_TagQuark_py;   //!
   TBranch        *b_W_TagQuark_pz;   //!
   TBranch        *b_W_TagQuark_E;   //!
   TBranch        *b_W_Parton_pt;   //!
   TBranch        *b_W_TagQuark_pt;   //!



  WJet->SetBranchAddress ("numPFCorJets",       &numPFCorJets,      &b_numPFCorJets      ) ;
  WJet->SetBranchAddress ("numPFCorJetBTags",   &numPFCorJetBTags,  &b_numPFCorJetBTags  ) ;
  WJet->SetBranchAddress ("JetPFCor_Px",        JetPFCor_Px,        &b_JetPFCor_Px       ) ;
  WJet->SetBranchAddress ("JetPFCor_Py",        JetPFCor_Py,        &b_JetPFCor_Py       ) ;
  WJet->SetBranchAddress ("JetPFCor_Pz",        JetPFCor_Pz,        &b_JetPFCor_Pz       ) ;
  WJet->SetBranchAddress ("JetPFCor_E",         JetPFCor_E,         &b_JetPFCor_E        ) ;
  WJet->SetBranchAddress ("JetPFCorVBFTag_Pt",  JetPFCorVBFTag_Pt,  &b_JetPFCorVBFTag_Pt ) ;
  WJet->SetBranchAddress ("JetPFCorVBFTag_Eta", JetPFCorVBFTag_Eta, &b_JetPFCorVBFTag_Eta) ;
  WJet->SetBranchAddress ("JetPFCorVBFTag_Phi", JetPFCorVBFTag_Phi, &b_JetPFCorVBFTag_Phi) ;
  WJet->SetBranchAddress ("JetPFCorVBFTag_E",   JetPFCorVBFTag_E,   &b_JetPFCorVBFTag_E  ) ;


  WJet->SetBranchAddress ("W_H_mass_gen", &W_H_mass_gen, &b_W_H_mass_gen) ;
  WJet->SetBranchAddress ("W_H_px_gen",   &W_H_px_gen,   &b_W_H_px_gen  ) ;
  WJet->SetBranchAddress ("W_H_py_gen",   &W_H_py_gen,   &b_W_H_py_gen  ) ;
  WJet->SetBranchAddress ("W_H_pz_gen",   &W_H_pz_gen,   &b_W_H_pz_gen  ) ;
  WJet->SetBranchAddress ("W_H_e_gen",    &W_H_e_gen,    &b_W_H_e_gen   ) ;

   WJet->SetBranchAddress("W_Parton_px[2]", W_Parton_px, &b_W_Parton_px);
   WJet->SetBranchAddress("W_Parton_py[2]", W_Parton_py, &b_W_Parton_py);
   WJet->SetBranchAddress("W_Parton_pz[2]", W_Parton_pz, &b_W_Parton_pz);
   WJet->SetBranchAddress("W_Parton_E[2]", W_Parton_E, &b_W_Parton_E);

   WJet->SetBranchAddress("W_TagQuark_px[2]", W_TagQuark_px, &b_W_TagQuark_px);
   WJet->SetBranchAddress("W_TagQuark_py[2]", W_TagQuark_py, &b_W_TagQuark_py);
   WJet->SetBranchAddress("W_TagQuark_pz[2]", W_TagQuark_pz, &b_W_TagQuark_pz);
   WJet->SetBranchAddress("W_TagQuark_E[2]", W_TagQuark_E, &b_W_TagQuark_E);

   WJet->SetBranchAddress("W_Parton_pt[2]", W_Parton_pt, &b_W_Parton_pt);
   WJet->SetBranchAddress("W_TagQuark_pt[2]", W_TagQuark_pt, &b_W_TagQuark_pt);

   WJet->SetBranchAddress("W_muon_px", &W_muon_px, &b_W_muon_px);
   WJet->SetBranchAddress("W_muon_py", &W_muon_py, &b_W_muon_py);
   WJet->SetBranchAddress("W_muon_pz", &W_muon_pz, &b_W_muon_pz);
   WJet->SetBranchAddress("W_muon_e", &W_muon_e, &b_W_muon_e);
   WJet->SetBranchAddress("W_muon_pt", &W_muon_pt, &b_W_muon_pt);
   WJet->SetBranchAddress("W_muon_eta", &W_muon_eta, &b_W_muon_eta);

   WJet->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_event_met_pfmet);
   WJet->SetBranchAddress("event_met_pfmetPhi", &event_met_pfmetPhi, &b_event_met_pfmetPhi);
   WJet->SetBranchAddress("W_muon_dz000", &W_muon_dz000, &b_W_muon_dz000);
   WJet->SetBranchAddress("W_muon_dzPV", &W_muon_dzPV, &b_W_muon_dzPV);




   //PG histograms
//method 1 hist
// tag jets
   TH1F h_vbf_jj_pt ("h_vbf_jj_pt", "h_vbf_jj_pt", 100, 0, 500) ;
   TH1F h_vbf_jj_eta ("h_vbf_jj_eta", "h_vbf_jj_eta", 100, -5, 5) ;
   TH1F h_vbf_jj_e ("h_vbf_jj_e", "h_vbf_jj_e", 100, 0, 500) ;
   TH1F h_vbf_jj_phi ("h_vbf_jj_phi", "h_vbf_jj_phi", 100, -3, 3) ;

   TH1F h_vbf_aj_pt ("h_vbf_aj_pt", "h_vbf_aj_pt", 100, 0, 500) ;
   TH1F h_vbf_aj_eta ("h_vbf_aj_eta", "h_vbf_aj_eta", 100, -5, 5) ;
   TH1F h_vbf_aj_e ("h_vbf_aj_e", "h_vbf_aj_e", 100, 0, 500) ;
   TH1F h_vbf_aj_phi ("h_vbf_aj_phi", "h_vbf_aj_phi", 100, -3, 3) ;

   TH1F h_vbf_bj_pt ("h_vbf_bj_pt", "h_vbf_bj_pt", 100, 0, 500) ;
   TH1F h_vbf_bj_eta ("h_vbf_bj_eta", "h_vbf_bj_eta", 100, -5, 5) ;
   TH1F h_vbf_bj_e ("h_vbf_bj_e", "h_vbf_bj_e", 100, 0, 500) ;
   TH1F h_vbf_bj_phi ("h_vbf_bj_phi", "h_vbf_bj_phi", 100, -3, 3) ;

   TH1F h_vbf_jj_deta ("h_vbf_jj_deta", "h_vbf_jj_deta", 100, -5, 5) ;
   TH1F h_vbf_jj_dphi ("h_vbf_jj_dphi", "h_vbf_jj_dphi", 100, -3, 3) ;
// W JETS
   TH1F h_vbf_wjj_pt ("h_vbf_wjj_pt", "h_vbf_wjj_pt", 100, 0, 500) ;
   TH1F h_vbf_wjj_eta ("h_vbf_wjj_eta", "h_vbf_wjj_eta", 100, -5, 5) ;
   TH1F h_vbf_wjj_e ("h_vbf_wjj_e", "h_vbf_wjj_e", 100, 0, 500) ;
   TH1F h_vbf_wjj_phi ("h_vbf_wjj_phi", "h_vbf_wjj_phi", 100, -3, 3) ;
   TH1F h_vbf_wjj_m ("h_vbf_wjj_m", "h_vbf_wjj_m", 100, 0, 500) ;

   TH1F h_vbf_waj_pt ("h_vbf_waj_pt", "h_vbf_waj_pt", 100, 0, 500) ;
   TH1F h_vbf_waj_eta ("h_vbf_waj_eta", "h_vbf_waj_eta", 100, -5, 5) ;
   TH1F h_vbf_waj_e ("h_vbf_waj_e", "h_vbf_waj_e", 100, 0, 500) ;
   TH1F h_vbf_waj_phi ("h_vbf_waj_phi", "h_vbf_waj_phi", 100, -3, 3) ;

   TH1F h_vbf_wbj_pt ("h_vbf_wbj_pt", "h_vbf_wbj_pt", 100, 0, 500) ;
   TH1F h_vbf_wbj_eta ("h_vbf_wbj_eta", "h_vbf_wbj_eta", 100, -5, 5) ;
   TH1F h_vbf_wbj_e ("h_vbf_wbj_e", "h_vbf_wbj_e", 100, 0, 500) ;
   TH1F h_vbf_wbj_phi ("h_vbf_wbj_phi", "h_vbf_wbj_phi", 100, -3, 3) ;
   
   TH1F h_vbf_lvjj_pt ("h_vbf_lvjj_pt", "h_vbf_lvjj_pt", 100, 0, 500) ;
   TH1F h_vbf_lvjj_eta ("h_vbf_lvjj_eta", "h_vbf_lvjj_eta", 100, -5, 5) ;
   TH1F h_vbf_lvjj_e ("h_vbf_lvjj_e", "h_vbf_lvjj_e", 100, 0, 500) ;
   TH1F h_vbf_lvjj_phi ("h_vbf_lvjj_phi", "h_vbf_lvjj_phi", 100, -3, 3) ;
   TH1F h_vbf_lvjj_m ("h_vbf_lvjj_m", "h_vbf_lvjj_m", 100, 0, 1000) ;
   TH1F h_vbf_lvj_m ("h_vbf_lvj_m", "h_vbf_lvj_m", 100, 0, 500) ;

   TH1F h_vbf_wjj_deta ("h_vbf_wjj_deta", "h_vbf_wjj_deta", 100, -5, 5) ;
   TH1F h_vbf_wjj_dphi ("h_vbf_wjj_dphi", "h_vbf_wjj_dphi", 100, -3, 3) ;
   TH1F h_vbf_jj_m ("h_vbf_jj_m", "h_vbf_jj_m", 100, 0, 2000) ;
	
/*	tmva_t->Branch("vbf_event",&vbf_event, "vbf_event/I");
	tmva_t->Branch("vbf_aj_id",&vbf_aj_id, "vbf_aj_id/I");
	tmva_t->Branch("vbf_bj_id",&vbf_bj_id, "vbf_bj_id/I");
	tmva_t->Branch("vbf_waj_id",&vbf_waj_id, "vbf_waj_id/I");
	tmva_t->Branch("vbf_wbj_id",&vbf_wbj_id, "vbf_wbj_id/I");		
*/


   int countevt=0;
   Long64_t nentries = WJet->GetEntries();
   cout << nentries << endl ;
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += WJet->GetEntry(i);
        // if (jentry % 1000 == 0) cout << "evnt " << jentry << endl ;
        lorentzVector W_H_gen ;
        vector <lorentzVector> jets ;
        vector <lorentzVector> genjets ;
	vector <lorentzVector> genjets_v;
//        vector<int> genUsed (genjets.size (), 0) ; //PG vector of flags to know if a jet is already used

        int numFwdJets = 0 ;

        for (int i = 0 ; i < 8 ; ++i) if (JetPFCorVBFTag_Pt[i] > 20.) ++numFwdJets ;
	int numGFwdJets = 0 ;

        for (int i = 0 ; i < 2 ; ++i) if (W_TagQuark_pt[i] > 20.) ++numGFwdJets ;

        int numGPFJets = 0 ;

        for (int i = 0 ; i < 2 ; ++i) if (W_Parton_pt[i] > 20.) ++numGPFJets ;

        int jetsNum = numPFCorJets + numFwdJets ;
	//cout<<"before continue "<<jetsNum<<endl;
        if (jetsNum < 3) continue ;
	//cout<<"after continue "<<jetsNum<<endl;
        int GjetsNum = numGFwdJets + numGPFJets ;
       // if (GjetsNum < 3) continue ;
       //PG read jets from the event and merge the collections
        for (int i = 0 ; i < numPFCorJets ; ++i)
          {
            lorentzVector dummy ;
            dummy.SetPxPyPzE (
                JetPFCor_Px[i],
                JetPFCor_Py[i],
                JetPFCor_Pz[i],
                JetPFCor_E[i]
              ) ;
            jets.push_back (dummy) ;  
          }
  
        for (int i = 0 ; i < numFwdJets ; ++i)
          {
            lorentzVector dummy ;
            dummy.SetPxPyPzE (
                hep_to_polar_x (JetPFCorVBFTag_Pt[i], JetPFCorVBFTag_Eta[i], JetPFCorVBFTag_Phi[i]),
                hep_to_polar_y (JetPFCorVBFTag_Pt[i], JetPFCorVBFTag_Eta[i], JetPFCorVBFTag_Phi[i]),
                hep_to_polar_z (JetPFCorVBFTag_Pt[i], JetPFCorVBFTag_Eta[i], JetPFCorVBFTag_Phi[i]),
                JetPFCorVBFTag_E[i]
              ) ;
            jets.push_back (dummy) ;  
          }
        sort (jets.rbegin (), jets.rend (), local_TLVP_ptSort ()) ;
//	cout<<" jets size  after sorting"<<jets.size()<<endl; 
	for (int i =0; i<jets.size();i++)
	{
	//cout<<"sorted jets pt "<<jets.at (i).Pt()<<endl; //to check jets are sorted 
	}
       /* for (int i = 0 ; i < numGPFJets ; ++i)
          {
	//if(W_Parton_Pt[i]<20)continue;	
            lorentzVector dummy1 ;
            dummy1.SetPxPyPzE (
                W_Parton_px[i],
                W_Parton_py[i],
                W_Parton_pz[i],
                W_Parton_E[i]
              ) ;
            genjets.push_back (dummy1) ;
          }
    
        for (int i = 0 ; i < numGFwdJets ; ++i)
          {
	 //if(W_TagQuark_Pt[i]<20)continue;
            lorentzVector dummy1 ;
            dummy1.SetPxPyPzE (
                W_TagQuark_px[i],
                W_TagQuark_py[i],
                W_TagQuark_pz[i],
                W_TagQuark_E[i]
              ) ;
            genjets.push_back (dummy1) ;
          }
      sort (genjets.rbegin (), genjets.rend (), local_TLVP_ptSort ()) ;
*/


	// method 1

	 lorentzVector vbf_ajp(0,0,0,0), vbf_bjp(0,0,0,0);
         lorentzVector wjj_ajp(0,0,0,0), wjj_bjp(0,0,0,0); 
         float best_detatagjj = 0; // float best_mtagjj =0;
         float best_mjj = 0; // float best_mjj =0;
         float c_mjj =1000;
 	 float Jpt = 25.0;

   float vbf_jj_e =-999,   vbf_jj_pt =-999,   vbf_jj_eta=-999,  vbf_jj_phi =-999, vbf_jj_m=-999;
   float vbf_aj_e =-999,   vbf_aj_pt =-999,   vbf_aj_eta=-999,  vbf_aj_phi =-999;   
   float vbf_bj_e =-999,   vbf_bj_pt =-999,   vbf_bj_eta=-999,  vbf_bj_phi =-999;   
   float vbf_jj_deta=-999; Float_t vbf_jj_dphi=-999; 
    int   vbf_jj_type=0,   vbf_n_excj=0,   vbf_n_exfj=0,   vbf_n_gdjj=0;

   float vbf_wjj_e =-999,   vbf_wjj_pt =-999,   vbf_wjj_eta=-999,  vbf_wjj_phi =-999, vbf_wjj_m=-999;   
   float vbf_waj_e =-999,   vbf_waj_pt =-999,   vbf_waj_eta=-999,  vbf_waj_phi =-999;   
   float vbf_wbj_e =-999,   vbf_wbj_pt =-999,   vbf_wbj_eta=-999,  vbf_wbj_phi =-999;   
   float vbf_wjj_deta=-999; Float_t vbf_wjj_dphi=-999;
   float vbf_lvjj_e =-999,   vbf_lvjj_pt =-999,   vbf_lvjj_eta=-999,  vbf_lvjj_phi =-999, vbf_lvjj_m =-999, vbf_lvj_m=-999 ;


// 	 float jess 1.0; 
         int   n_excj =0, n_exfj = 0, n_gdjj = 0, jj_type = 0, tag_i_id = -1, tag_j_id = -1, wjj_a_id = -1, wjj_b_id = -1;

	// method 3 A
        for ( int k=0; k < (int) jets.size(); ++k) // loop over jets
                {
                  if (fabs(jets.at(k).Eta()) > 2.5) continue;    //choose w jets from central region closest to W 
                        lorentzVector k_p;
                        k_p.SetPxPyPzE (
                        hep_to_polar_x (jets.at(k).Pt(), jets.at(k).Eta(), jets.at(k).Phi()),
                        hep_to_polar_y (jets.at(k).Pt(), jets.at(k).Eta(), jets.at(k).Phi()),
                        hep_to_polar_z (jets.at(k).Pt(), jets.at(k).Eta(), jets.at(k).Phi()),
                        jets.at(k).E()
                        );
                   for ( int m=k+1; m < (int) jets.size(); ++m) // loop over  jets
                    {
                    if (fabs(jets.at(m).Eta()) > 2.5) continue;   //cancel events with jet pt<25
                        lorentzVector m_p;
                       m_p.SetPxPyPzE (
                        hep_to_polar_x (jets.at(m).Pt(), jets.at(m).Eta(), jets.at(m).Phi()),
                        hep_to_polar_y (jets.at(m).Pt(), jets.at(m).Eta(), jets.at(m).Phi()),
                        hep_to_polar_z (jets.at(m).Pt(), jets.at(m).Eta(), jets.at(m).Phi()),
                        jets.at(m).E()
                        );
                        if ( fabs(sqrt((k_p+m_p).M2())-80.0) < c_mjj )
                        {
                              c_mjj =fabs(sqrt((k_p+m_p).M2())-80.0);
                                //cout<<" c_mjj  "<<c_mjj<<endl;        
                                wjj_a_id = k;
                                wjj_b_id = m;
                                wjj_ajp=k_p;
                                wjj_bjp=m_p;
                              //cout<<"  "<<ie<<"  "<<wjj_a_id<<"   "<<wjj_b_id<<"   "<<c_mjj<<endl;

                        } //loop to choose best mjj jets (jets closest to mass of W)
                   }//internal loop over  jets 
                }  //loop over  jets
        if (wjj_a_id !=-1 && wjj_b_id != -1)
        	{       // Find two  W jets closeset to W
                        vbf_wjj_e      = (wjj_ajp+wjj_bjp).E();
                        vbf_wjj_pt     = (wjj_ajp+wjj_bjp).Pt();
                        vbf_wjj_eta    = (wjj_ajp+wjj_bjp).Eta();
                        vbf_wjj_phi    = (wjj_ajp+wjj_bjp).Phi();
                        vbf_wjj_m      = sqrt((wjj_ajp+wjj_bjp).M2());
                        vbf_waj_e      = (wjj_ajp).E();
                        vbf_waj_pt     = (wjj_ajp).Pt();
                        vbf_waj_eta    = (wjj_ajp).Eta();
                        vbf_waj_phi    = (wjj_ajp).Phi();
                        //vbf_waj_m      = (wjj_ajp).M();
                        vbf_wbj_e      = (wjj_bjp).E();
                        vbf_wbj_pt     = (wjj_bjp).Pt();
                        vbf_wbj_eta    = (wjj_bjp).Eta();
                        vbf_wbj_phi    = (wjj_bjp).Phi();
                        vbf_wjj_deta=vbf_waj_eta-vbf_wbj_eta;
                        vbf_wjj_dphi= vbf_waj_phi-vbf_wbj_phi;

		}

               // vbf tag jet pair method5 B
	if (wjj_a_id!=-1 && wjj_b_id!=-1)
	{
         for ( size_t i=0; i < jets.size(); ++i) 
         {
	if (fabs(jets.at(i).Eta())>4.7) continue;
	  if (i == wjj_a_id || i== wjj_b_id ) continue;
         	lorentzVector i_p;
            	i_p.SetPxPyPzE (
                hep_to_polar_x (jets.at(i).Pt(), jets.at(i).Eta(), jets.at(i).Phi()),
                hep_to_polar_y (jets.at(i).Pt(), jets.at(i).Eta(), jets.at(i).Phi()),
                hep_to_polar_z (jets.at(i).Pt(), jets.at(i).Eta(), jets.at(i).Phi()),
                jets.at(i).E()
              );

            for (size_t j=i+1; j <jets.size(); ++j) 
            {
            if (fabs(jets.at(j).Eta()) > 4.7) continue;
          if (j == wjj_a_id || j== wjj_b_id ) continue;
                lorentzVector j_p;
            	j_p.SetPxPyPzE (
                hep_to_polar_x (jets.at(j).Pt(), jets.at(j).Eta(), jets.at(j).Phi()),
                hep_to_polar_y (jets.at(j).Pt(), jets.at(j).Eta(), jets.at(j).Phi()),
                hep_to_polar_z (jets.at(j).Pt(), jets.at(j).Eta(), jets.at(j).Phi()),
                jets.at(j).E()
              ) ;
  //             if ( (jets.at(i).Eta()*jets.at(j).Eta())>0 )  continue;     // 1.  have to be one forward, one backward
  //             if ( (fabs(jets.at(i).Eta()-jets.at(j).Eta())<3.5) || (sqrt ((i_p+j_p).M2()) <500)) continue;// 2.Tag pair delta eta>3.5, Mjj>500
		//cout<<"  "<<sqrt ((i_p+j_p).M2())<<endl;
               // if find more than one combinations
               //if ( (fas(i_Eta-j_Eta)>best_detatagjj) )      // 3   Select best combination with maximum deta Eta
               if ( sqrt((i_p+j_p).M2()) > best_mjj )
               {                          // 3   Select best combination with maximum Mjj because of the bad angular resolution in the HF
                  best_detatagjj = fabs(jets.at(i).Eta()-jets.at(j).Eta()); 
		  n_gdjj++;
                  best_mjj = sqrt ((i_p+j_p).M2());
                  tag_i_id = i; 
		  tag_j_id = j; 
		  vbf_ajp = i_p; 
		  vbf_bjp = j_p;

               } //loop to find two tag jets with highest mjj
           } // over over loop over Nmax-1 reco jets inside 
        } //loop over Nmax reco jets 
	}
	        // method5 B
	
		if (tag_i_id !=-1 && tag_j_id != -1)
		{
                  vbf_jj_e      = (vbf_ajp+ vbf_bjp).E();
                  vbf_jj_pt     = (vbf_ajp+ vbf_bjp).Pt();
                  vbf_jj_eta    = (vbf_ajp+ vbf_bjp).Eta();
                  vbf_jj_phi    = (vbf_ajp+ vbf_bjp).Phi();
                  vbf_jj_m      = sqrt ((vbf_ajp+vbf_bjp).M2());
                  vbf_aj_e      = (vbf_ajp).E();
                  vbf_aj_pt     = (vbf_ajp).Pt();
                  vbf_aj_eta    = (vbf_ajp).Eta();
                  vbf_aj_phi    = (vbf_ajp).Phi();
                  //vbf_aj_m      = (i_p).M();
                  vbf_bj_e      = (vbf_bjp).E();
                  vbf_bj_pt     = (vbf_bjp).Pt();
                  vbf_bj_eta    = (vbf_bjp).Eta();
                  vbf_bj_phi    = (vbf_bjp).Phi();
                 // vbf_bj_m      = (j_p).M();

                  vbf_jj_deta   =vbf_aj_eta-vbf_bj_eta;
                  vbf_jj_dphi   = vbf_aj_phi-vbf_bj_phi;
                  //cout<<"  "<<vbf_jj_dphi<<endl;
        } //loop  
	// method5 B

	lorentzVector  lepton;
        if (W_muon_pt < 25. || fabs(W_muon_dz000)>0.02 || fabs(W_muon_dzPV)>0.5 || fabs(W_muon_eta)>2.5) continue;
	lepton.SetPxPyPzE(W_muon_px,W_muon_py,W_muon_pz,W_muon_e);
	lorentzVector nutrino;
        if (event_met_pfmet <25) continue;
                nutrino.SetPxPyPzE (
                hep_to_polar_x (event_met_pfmet,0.,event_met_pfmetPhi),
                hep_to_polar_y (event_met_pfmet,0.,event_met_pfmetPhi),
                hep_to_polar_z (event_met_pfmet,0.,event_met_pfmetPhi),
//                0.,
                sqrt( pow((hep_to_polar_x (event_met_pfmet,0.,event_met_pfmetPhi)),2)+pow((hep_to_polar_y (event_met_pfmet,0.,event_met_pfmetPhi)),2)));
		//cout<<tag_i_id<<"   "<<tag_j_id<<endl;
//        if (wjj_a_id !=-1 && wjj_b_id != -1 && tag_i_id !=-1 && tag_j_id != -1)
//		{
      		    vbf_lvjj_e      = (lepton+nutrino+wjj_ajp+wjj_bjp).E();
		    vbf_lvjj_pt     = (lepton+nutrino+wjj_ajp+wjj_bjp).Pt();
        	    vbf_lvjj_eta    = (lepton+nutrino+wjj_ajp+wjj_bjp).Eta();
		    vbf_lvjj_phi    = (lepton+nutrino+wjj_ajp+wjj_bjp).Phi();
		    vbf_lvjj_m      = sqrt ((lepton+nutrino+wjj_ajp+wjj_bjp).M2());
                    vbf_lvj_m      = sqrt((lepton+nutrino).M2());

	//			cout<<" four body mass "<<vbf_lvjj_m<<endl;
		           // vbf_lvjj_y      = (lepton+nutrino+wjj_ajp+wjj_bjp).Rapidity();
			//cout<<ie<<"  "<<tag_i_id<<"  "<<tag_j_id<<"   "<<wjj_a_id<<"  "<<wjj_b_id<<"   "<<vbf_waj_eta<<endl;
        if (wjj_a_id !=-1 && wjj_b_id != -1 && tag_i_id !=-1 && tag_j_id != -1 && vbf_lvjj_m>0 && vbf_lvj_m >30. && vbf_wjj_m >30. )
                {
                if (vbf_aj_eta*vbf_bj_eta>0) continue;
                if (( (fabs(vbf_jj_deta))<3.5) || vbf_jj_m < 500) continue;
                                cout<<" four body mass "<<vbf_lvjj_m<<endl;

	// fill all histo 
	h_vbf_jj_pt.Fill(vbf_jj_pt);
        h_vbf_jj_eta.Fill(vbf_jj_eta);
        h_vbf_jj_phi.Fill(vbf_jj_phi);
        h_vbf_jj_e.Fill(vbf_jj_e);
        h_vbf_jj_m.Fill(vbf_jj_m);

        h_vbf_aj_pt.Fill(vbf_aj_pt);
        h_vbf_aj_eta.Fill(vbf_aj_eta);
        h_vbf_aj_phi.Fill(vbf_aj_phi);
        h_vbf_aj_e.Fill(vbf_aj_e);

        h_vbf_bj_pt.Fill(vbf_bj_pt);
        h_vbf_bj_eta.Fill(vbf_bj_eta);
        h_vbf_bj_phi.Fill(vbf_bj_phi);
        h_vbf_bj_e.Fill(vbf_bj_e);

        h_vbf_jj_deta.Fill(vbf_jj_deta);
        h_vbf_jj_dphi.Fill(vbf_jj_dphi);

        h_vbf_wjj_pt.Fill(vbf_wjj_pt);
        h_vbf_wjj_eta.Fill(vbf_wjj_eta);
        h_vbf_wjj_phi.Fill(vbf_wjj_phi);
        h_vbf_wjj_e.Fill(vbf_wjj_e);
        h_vbf_wjj_m.Fill(vbf_wjj_m);

        h_vbf_waj_pt.Fill(vbf_waj_pt);
        h_vbf_waj_eta.Fill(vbf_waj_eta);
        h_vbf_waj_phi.Fill(vbf_waj_phi);
        h_vbf_waj_e.Fill(vbf_waj_e);

        h_vbf_wbj_pt.Fill(vbf_wbj_pt);
        h_vbf_wbj_eta.Fill(vbf_wbj_eta);
        h_vbf_wbj_phi.Fill(vbf_wbj_phi);
        h_vbf_wbj_e.Fill(vbf_wbj_e);

        h_vbf_wjj_deta.Fill(vbf_wjj_deta);
        h_vbf_wjj_dphi.Fill(vbf_wjj_dphi);

        h_vbf_lvjj_pt.Fill(vbf_lvjj_pt);
        h_vbf_lvjj_eta.Fill(vbf_lvjj_eta);
        h_vbf_lvjj_phi.Fill(vbf_lvjj_phi);
        h_vbf_lvjj_e.Fill(vbf_lvjj_e);
        h_vbf_lvjj_m.Fill(vbf_lvjj_m);
        h_vbf_lvj_m.Fill(vbf_lvj_m);
	
       }  //loop closing for having four jets
//method 1 ends
                               // cout<<" four body mass "<<vbf_lvjj_m<<endl;

    }   //PG loop over events

  TString outFileName = "TTbar_output_" ;
//  outFileName += maxNumber;
  outFileName += "method5";
  outFileName += "_";
  outFileName += int (maxDR * 100);
  outFileName += ".root" ;
  TFile outFile (outFileName, "recreate") ;

  h_vbf_jj_pt.Write ();
  h_vbf_jj_eta.Write ();
  h_vbf_jj_phi.Write ();
  h_vbf_jj_e.Write ();
  h_vbf_jj_m.Write ();

  h_vbf_aj_pt.Write ();
  h_vbf_aj_eta.Write ();
  h_vbf_aj_phi.Write ();
  h_vbf_aj_e.Write ();

  h_vbf_bj_pt.Write ();
  h_vbf_bj_eta.Write ();
  h_vbf_bj_phi.Write ();
  h_vbf_bj_e.Write ();

  h_vbf_jj_deta.Write ();
  h_vbf_jj_dphi.Write ();

  h_vbf_wjj_pt.Write ();
  h_vbf_wjj_eta.Write ();
  h_vbf_wjj_phi.Write ();
  h_vbf_wjj_e.Write ();
  h_vbf_wjj_m.Write ();

  h_vbf_waj_pt.Write ();
  h_vbf_waj_eta.Write ();
  h_vbf_waj_phi.Write ();
  h_vbf_waj_e.Write ();

  h_vbf_wbj_pt.Write ();
  h_vbf_wbj_eta.Write ();
  h_vbf_wbj_phi.Write ();
  h_vbf_wbj_e.Write ();

  h_vbf_lvjj_pt.Write ();
  h_vbf_lvjj_eta.Write ();
  h_vbf_lvjj_phi.Write ();
  h_vbf_lvjj_e.Write ();
  h_vbf_lvjj_m.Write ();
  h_vbf_lvj_m.Write ();

  h_vbf_wjj_deta.Write ();
  h_vbf_wjj_dphi.Write ();

  outFile.Close () ;
 return 0 ;
}

