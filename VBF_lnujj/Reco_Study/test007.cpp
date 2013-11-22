// c++ -o read_007 `root-config --glibs --cflags` -lm read_007.cpp
//  ./read_007 /Users/govoni/data/ReducedTrees/RD_el_VBFHWWMH900_CMSSW532_private.root
//  ./read_007 /Users/govoni/data/ReducedTrees/RD_mu_VBFHWWMH300_CMSSW532_private.root
//  /home/ajay/Work/LPC/Morion13/VBF/RD_mu_VBFHWWMH300_CMSSW532_private.root

/*
- choose the groomed jet belonging to the W boson in the signal case
- clean AK5 jets from CA[6] ones and search for VBF topology
- choose the VBF jets as the two with the highest pTs
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
//using namespace ROOT::Math;
//using namespace std ;
//typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;



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


// function defined..................................below.................
float delta_phi(float phi1, float phi2)
{
        const float PI=2.0*acos(0.);
        const float TWOPI=2.0*PI;
        float PHI=fabs(phi1-phi2);
        return (PHI<=PI)? PHI : TWOPI-PHI;
}



/*
struct TLVP_ptSort: public std::binary_function<int, int, bool>
{
  bool operator() (const lorentzVector * x, const lorentzVector * y)
    {
      return x->Pt () < y->Pt () ;
    }
} ;

*/

int main (int argc, char** argv)
{
  if (argc != 3)
    {
      std::cout << ">>> Usage:   " << argv[0] << "   treeFile.root maxNumber deltaR" << std::endl;
      return -1;
    }

  int maxNumber = atoi (argv[1]) ;
  double maxDR = atof (argv[2]) ;

  //PG input tree, variables, branches
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

//  TString inputFile = argv[1] ;
//  TChain * WJet = new TChain ("WJet") ;
//  WJet->Add ("/home/ajay/Work/LPC/Morion13/VBF/RD_mu_VBFHWWMH300_CMSSW532_private.root") ;
//  WJet->Add (inputFile) ;
//  if (WJet == 0) return 1 ;

   gROOT->Reset();
 TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RD_mu_VBFHWWMH300_CMSSW532_private.root");
   if (!f) {
     f = new TFile("RD_mu_VBFHWWMH300_CMSSW532_private.root");
   }
	TTree * WJet;   
 f->GetObject("WJet",WJet);



// Declaration of leaf types
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


  //PG the Higgs boson MC truth
  Float_t         W_H_mass_gen;
  Float_t         W_H_px_gen;
  Float_t         W_H_py_gen;
  Float_t         W_H_pz_gen;
  Float_t         W_H_e_gen;

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

   //PG histograms
   //PG ----------




   TH1F h_deltaR ("h_deltaR", "h_deltaR", 10, -5.0, 5.0) ;
   TH1F h_matchedSize ("h_matchedSize", "h_matchedSize", 10, 0, 10) ;
   TH1F h_parton1_pt ("h_parton1_pt", "h_parton1_pt", 100, 0, 500) ;
   TH1F h_parton1_eta ("h_parton1_eta", "h_parton1_eta", 100, -5, 5) ;
   TH1F h_parton2_pt ("h_parton2_pt", "h_parton2_pt", 100, 0, 500) ;
   TH1F h_parton2_eta ("h_parton2_eta", "h_parton2_eta", 100, -5, 5) ;
   TH1F h_parton3_pt ("h_parton3_pt", "h_parton3_pt", 100, 0, 500) ;
   TH1F h_parton3_eta ("h_parton3_eta", "h_parton3_eta", 100, -5, 5) ;
   TH1F h_parton4_pt ("h_parton4_pt", "h_parton4_pt", 100, 0, 500) ;
   TH1F h_parton4_eta ("h_parton4_eta", "h_parton4_eta", 100, -5, 5) ;




   //PG loop over events
/*   if (WJet == 0) return 0;
   Long64_t nentries = WJet->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int JetIndex[8] = {0,0,0,0,0,0,0,0};

   cout << nentries << endl ;
   nentries = 500000 ;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
//  for (Long64_t jentry=0; jentry<10;jentry++) 
      // if (Cut(ientry) < 0) continue;
*/
int countevt=0;

   Long64_t nentries = WJet->GetEntries();
   cout << nentries << endl ;
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += WJet->GetEntry(i);
//        if (jentry % 1000 == 0) cout << "evnt " << jentry << endl ;
        lorentzVector W_H_gen ;
        vector <lorentzVector> jets ;
        vector <lorentzVector> genjets ;
//        vector <lorentzVector> genjets_v;

      //  nb = WJet->GetEntry(jentry);   nbytes += nb;
  
        //PG build the total default jet collection
        //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        int numFwdJets = 0 ;
        for (int i = 0 ; i < 8 ; ++i) if (JetPFCorVBFTag_Pt[i] > 1.) ++numFwdJets ;
  
        int numGFwdJets = 0 ;
        for (int i = 0 ; i < 2 ; ++i) if (W_TagQuark_pt[i] > 20.) ++numGFwdJets ;
  
        int numGPFJets = 0 ;
        for (int i = 0 ; i < 2 ; ++i) if (W_Parton_pt[i] > 20.) ++numGPFJets ;
  
        int jetsNum = numPFCorJets + numFwdJets ;
        if (jetsNum < 3) continue ;
  
        int GjetsNum = numGFwdJets + numGPFJets ;
        if (GjetsNum < 3) continue ;
//int countevt=0;
if(GjetsNum==4)countevt++;
//cout<<countevt<<endl; 
 
cout<<countevt<<endl;
 
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
//	 cout<<" jets size  after sorting"<<jets.size()<<endl; 
        for (int i = 0 ; i < numGPFJets ; ++i)
          {
//	if(W_Parton_Pt[i]<20)continue;	
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
 //       if(W_TagQuark_Pt[i]<20)continue;
            lorentzVector dummy1 ;
            dummy1.SetPxPyPzE (
                W_TagQuark_px[i],
                W_TagQuark_py[i],
                W_TagQuark_pz[i],
                W_TagQuark_E[i]
              ) ;
            genjets.push_back (dummy1) ;
          }
        int genjetpt20 =0;
	for( int i=0;i<genjets.size();i++)
	{if(genjets[i].Pt()<20) continue;
	//	cout<<genjets[i].Pt()<<endl;
	genjetpt20=1;
	}	 
      sort (genjets.rbegin (), genjets.rend (), local_TLVP_ptSort ()) ;
         //cout<<" genjets size  after sorting"<<genjets.size()<<endl;

        //PG take only the first maxNumber jets of the reco jets collection
        vector<lorentzVector> firstN ;
        for (int i = 0 ; i < maxNumber ; ++i)
          {
            if (jets.size () == i) break ;
            firstN.push_back (jets.at (i)) ;
          }
  
        //PG      deltaR      (recoIndex , genIndex)
        std::map <float, std::pair<int, int> > matchingsMap ;
    
        for (int iReco = 0 ; iReco < firstN.size (); ++iReco) // loop over reco jets
          {
            for (int iGen = 0 ; iGen < genjets.size (); ++iGen) // loop over gen particles
              {
                float deltaR = my_deltaR_function (
                           firstN.at (iReco).Eta (),firstN.at (iReco).Phi (), 
                           genjets.at (iGen).Eta (), genjets.at (iGen).Phi ()) ;
	        h_deltaR.Fill (deltaR) ;
                matchingsMap[deltaR] = std::pair<int, int> (iReco, iGen) ;
              }
          }
    
        vector<int> recoUsed (firstN.size (), 0) ;   //PG vector of flags to know if a jet is already used
        vector<int> genUsed (genjets.size (), 0) ; //PG vector of flags to know if a jet is already used
        std::map<float, std::pair<int, int> > cleanedMatchings ;
        // remove elements with jets that have been already used
        for (std::map<float, std::pair<int, int> >::const_iterator iMap = matchingsMap.begin () ;
             iMap != matchingsMap.end () ;
             ++iMap)
          {
            if (iMap->first > maxDR) break ;
            if (recoUsed.at (iMap->second.first) == 1 ||
                genUsed.at (iMap->second.second) == 1)
              continue ; 
            cleanedMatchings[iMap->first] = iMap->second ;
            recoUsed.at (iMap->second.first) = 1 ;
            genUsed.at (iMap->second.second) = 1 ;

          }  
          h_matchedSize.Fill (cleanedMatchings.size ()) ;
	if(genjets.size()<4)
	{	
//	cout<<"genjets size is less than 4"<<endl;
	}
	else
	{ 
        h_parton1_pt.Fill (genjets.at (0).Pt());
        h_parton1_eta.Fill (genjets.at (0).Eta());
        h_parton2_pt.Fill (genjets.at (1).Pt());
        h_parton2_eta.Fill (genjets.at (1).Eta());
        h_parton3_pt.Fill (genjets.at (2).Pt());
        h_parton3_eta.Fill (genjets.at (2).Eta());
        h_parton4_pt.Fill (genjets.at (3).Pt());
        h_parton4_eta.Fill (genjets.at (3).Eta());

	}
    }   //PG loop over events



  TString outFileName = "output_" ;
  outFileName += maxNumber;
  outFileName += "_";
  outFileName += int (maxDR * 100);
  outFileName += ".root" ;
  TFile outFile (outFileName, "recreate") ;
  h_matchedSize.Write () ;
  h_deltaR.Write () ;
  h_parton1_pt.Write ();
  h_parton1_eta.Write ();
  h_parton2_pt.Write ();
  h_parton2_eta.Write ();
  h_parton3_pt.Write ();
  h_parton3_eta.Write ();
  h_parton4_pt.Write ();
  h_parton4_eta.Write ();


  outFile.Close () ;

 return 0 ;
}

