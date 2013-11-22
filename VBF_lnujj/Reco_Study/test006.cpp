// c++ -o read_007 `root-config --glibs --cflags` -lm read_007.cpp


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
	//cout<<"   "<<maxDR<<endl;
	//PG input tree, variables, branches
	//PG ---- ---- ---- ---- ---- ---- ---- ---- ---- 

	//  TString inputFile = argv[1] ;
	//  TChain * WJet = new TChain ("WJet") ;
	//  WJet->Add ("/home/ajay/Work/LPC/Morion13/VBF/RD_mu_VBFHWWMH300_CMSSW532_private.root") ;
	//  WJet->Add (inputFile) ;
	//  if (WJet == 0) return 1 ;

	gROOT->Reset();
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("WmunuJetAnalysisntuple_1_1_gqg.root");
	if (!f) {
		f = new TFile("WmunuJetAnalysisntuple_1_1_gqg.root");
        //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RD_mu_W4Jets_CMSSW532_new.root");
        //if (!f) {
         //       f = new TFile("RD_mu_W4Jets_CMSSW532_new.root");
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

        Float_t         JetPFCor_QGLikelihood[8];
	Float_t         JetPFCorVBFTag_QGLikelihood[8];


	Float_t         W_gqmass_gen[20];
	Float_t         W_gq_px_gen[20];
	Float_t         W_gq_py_gen[20];
	Float_t         W_gq_pz_gen[20];
	Float_t         W_gq_e_gen[20];
	Float_t         W_gq_pt_gen[20];
	Float_t         W_gq_et_gen[20];
	Float_t         W_gq_eta_gen[20];
	Float_t         W_gq_phi_gen[20];
	Float_t         W_gq_vx_gen[20];
	Float_t         W_gq_vy_gen[20];
	Float_t         W_gq_vz_gen[20];
	Float_t         W_gq_y_gen[20];
	Int_t           W_gq_Id_gen[20];




	//PG the Higgs boson MC truth

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
        TBranch        *b_JetPFCor_QGLikelihood;   //!
        TBranch        *b_JetPFCorVBFTag_QGLikelihood;   //!


	TBranch        *b_W_gq_px_gen;   //!
	TBranch        *b_W_gq_py_gen;   //!
	TBranch        *b_W_gq_pz_gen;   //!
	TBranch        *b_W_gq_e_gen;   //!
	TBranch        *b_W_gq_pt_gen;   //!
	TBranch        *b_W_gq_et_gen;   //!
	TBranch        *b_W_gq_eta_gen;   //!
	TBranch        *b_W_gq_phi_gen;   //!
	TBranch        *b_W_gq_vx_gen;   //!
	TBranch        *b_W_gq_vy_gen;   //!
	TBranch        *b_W_gq_vz_gen;   //!
	TBranch        *b_W_gq_y_gen;   //!
	TBranch        *b_W_gq_Id_gen;   //!



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
        WJet->SetBranchAddress("JetPFCor_QGLikelihood",JetPFCor_QGLikelihood);
        WJet->SetBranchAddress("JetPFCorVBFTag_QGLikelihood",JetPFCorVBFTag_QGLikelihood);


	WJet->SetBranchAddress("W_gqmass_gen[20]",W_gqmass_gen);
	WJet->SetBranchAddress("W_gq_px_gen[20]",W_gq_px_gen);
	WJet->SetBranchAddress("W_gq_py_gen[20]",W_gq_py_gen);
	WJet->SetBranchAddress("W_gq_pz_gen[20]",W_gq_pz_gen);
	WJet->SetBranchAddress("W_gq_e_gen[20]",W_gq_e_gen);
	WJet->SetBranchAddress("W_gq_pt_gen[20]",W_gq_pt_gen);
	WJet->SetBranchAddress("W_gq_et_gen[20]",W_gq_et_gen);
	WJet->SetBranchAddress("W_gq_eta_gen[20]",W_gq_eta_gen);
	WJet->SetBranchAddress("W_gq_phi_gen[20]",W_gq_phi_gen);
	WJet->SetBranchAddress("W_gq_vx_gen[20]",W_gq_vx_gen);
	WJet->SetBranchAddress("W_gq_vy_gen[20]",W_gq_vy_gen);
	WJet->SetBranchAddress("W_gq_vz_gen[20]",W_gq_vz_gen);
	WJet->SetBranchAddress("W_gq_y_gen[20]",W_gq_y_gen);
	WJet->SetBranchAddress("W_gq_Id_gen[20]",W_gq_Id_gen);




	//PG histograms
	//PG ----------




	TH1F h_qdeltaR ("h_qdeltaR", "h_qdeltaR", 10, -5.0, 15.0) ;
	TH1F h_gdeltaR ("h_gdeltaR", "h_gdeltaR", 10, -5.0, 15.0) ;

	TH1F h_qmatchedSize ("h_qmatchedSize", "h_qmatchedSize", 10, 0, 10) ;
	TH1F h_gmatchedSize ("h_gmatchedSize", "h_gmatchedSize", 10, 0, 10) ;



	TH1F h_quark1_pt ("h_quark1_pt", "h_quark1_pt", 100, 0, 500) ;
	TH1F h_quark1_eta ("h_quark1_eta", "h_quark1_eta", 100, -5, 5) ;
	TH1F h_quark2_pt ("h_quark2_pt", "h_quark2_pt", 100, 0, 500) ;
	TH1F h_quark2_eta ("h_quark2_eta", "h_quark2_eta", 100, -5, 5) ;
	TH1F h_quark3_pt ("h_quark3_pt", "h_quark3_pt", 100, 0, 500) ;
	TH1F h_quark3_eta ("h_quark3_eta", "h_quark3_eta", 100, -5, 5) ;
	TH1F h_quark4_pt ("h_quark4_pt", "h_quark4_pt", 100, 0, 500) ;
	TH1F h_quark4_eta ("h_quark4_eta", "h_quark4_eta", 100, -5, 5) ;

	TH1F h_gluon1_pt ("h_gluon1_pt", "h_gluon1_pt", 100, 0, 500) ;
	TH1F h_gluon1_eta ("h_gluon1_eta", "h_gluon1_eta", 100, -5, 5) ;
	TH1F h_gluon2_pt ("h_gluon2_pt", "h_gluon2_pt", 100, 0, 500) ;
	TH1F h_gluon2_eta ("h_gluon2_eta", "h_gluon2_eta", 100, -5, 5) ;
	TH1F h_gluon3_pt ("h_gluon3_pt", "h_gluon3_pt", 100, 0, 500) ;
	TH1F h_gluon3_eta ("h_gluon3_eta", "h_gluon3_eta", 100, -5, 5) ;  
	TH1F h_gluon4_pt ("h_gluon4_pt", "h_gluon4_pt", 100, 0, 500) ;
	TH1F h_gluon4_eta ("h_gluon4_eta", "h_gluon4_eta", 100, -5, 5) ;



	TH1F h_quarklike_jet1_pt ("h_quarklike_jet1_pt", "h_quarklike_jet1_pt", 100, 0, 500) ;
	TH1F h_quarklike_jet1_eta ("h_quarklike_jet1_eta", "h_quarklike_jet1_eta", 100, -5, 5) ;
	TH1F h_quarklike_jet2_pt ("h_quarklike_jet2_pt", "h_quarklike_jet2_pt", 100, 0, 500) ;
	TH1F h_quarklike_jet2_eta ("h_quarklike_jet2_eta", "h_quarklike_jet2_eta", 100, -5, 5) ;
	TH1F h_quarklike_jet3_pt ("h_quarklike_jet3_pt", "h_quarklike_jet3_pt", 100, 0, 500) ;
	TH1F h_quarklike_jet3_eta ("h_quarklike_jet3_eta", "h_quarklike_jet3_eta", 100, -5, 5) ;
	TH1F h_quarklike_jet4_pt ("h_quarklike_jet4_pt", "h_quarklike_jet4_pt", 100, 0, 500) ;
	TH1F h_quarklike_jet4_eta ("h_quarklike_jet4_eta", "h_quarklike_jet4_eta", 100, -5, 5) ;



	TH1F h_gluonlike_jet1_pt ("h_gluonlike_jet1_pt", "h_gluonlike_jet1_pt", 100, 0, 500) ;
	TH1F h_gluonlike_jet1_eta ("h_gluonlike_jet1_eta", "h_gluonlike_jet1_eta", 100, -5, 5) ;
	TH1F h_gluonlike_jet2_pt ("h_gluonlike_jet2_pt", "h_gluonlike_jet2_pt", 100, 0, 500) ;
	TH1F h_gluonlike_jet2_eta ("h_gluonlike_jet2_eta", "h_gluonlike_jet2_eta", 100, -5, 5) ;
	TH1F h_gluonlike_jet3_pt ("h_gluonlike_jet3_pt", "h_gluonlike_jet3_pt", 100, 0, 500) ;
	TH1F h_gluonlike_jet3_eta ("h_gluonlike_jet3_eta", "h_gluonlike_jet3_eta", 100, -5, 5) ;
	TH1F h_gluonlike_jet4_pt ("h_gluonlike_jet4_pt", "h_gluonlike_jet4_pt", 100, 0, 500) ;
	TH1F h_gluonlike_jet4_eta ("h_gluonlike_jet4_eta", "h_gluonlike_jet4_eta", 100, -5, 5) ;





	//PG loop over events
	Long64_t nentries = WJet->GetEntries();
	cout << nentries << endl ;
	Long64_t nbytes = 0;
	for (Long64_t ent=0; ent<nentries;ent++) {
		nbytes += WJet->GetEntry(ent);
		//        if (jentry % 1000 == 0) cout << "evnt " << jentry << endl ;
		lorentzVector W_H_gen ;
		vector <lorentzVector> jets ;
                vector <lorentzVector> QGd ;

		vector <lorentzVector> quarks ;
		vector <lorentzVector> gluons ;
		vector <lorentzVector> quarklike_jets;
		vector <lorentzVector> gluonlike_jets;


		//        vector <lorentzVector> quarks_v;

		//PG build the total default jet collection
		//PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		int numFwdJets = 0 ;
		for (int i = 0 ; i < 8 ; ++i) if (JetPFCorVBFTag_Pt[i] > 1.) ++numFwdJets ;

		int numGPFJets = 0 ;
		//for (int i = 0 ; i < 2 ; ++i) if (W_Parton_pt[i] > 1.) ++numGPFJets ;
		for (int i = 0 ; i < 20 ; ++i) if ( W_gq_Id_gen[i]!=0 && W_gq_pt_gen[i] > 20.) ++numGPFJets ;

		int jetsNum = numPFCorJets + numFwdJets ;
		if (jetsNum < 1) continue ;

		int GjetsNum =  numGPFJets ;
		if (GjetsNum < 1) continue ;

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
		//		 cout<<" jets size  after sorting"<<jets.size()<<endl; 
		for (int i = 0 ; i < numGPFJets ; ++i)
		{
			if((abs(W_gq_Id_gen[i]))<7 && (W_gq_Id_gen[i]!=0))
			{
				lorentzVector dummy1 ;
				dummy1.SetPxPyPzE (
						W_gq_px_gen[i],
						W_gq_py_gen[i],
						W_gq_pz_gen[i],
						W_gq_e_gen[i]
						) ;
				quarks.push_back (dummy1) ;
			}
		}
		for (int i = 0 ; i < numGPFJets ; ++i)
		{
			if((abs(W_gq_Id_gen[i]))==21)
			{
				lorentzVector dummy1 ;
				dummy1.SetPxPyPzE (
						W_gq_px_gen[i],
						W_gq_py_gen[i],
						W_gq_pz_gen[i],
						W_gq_e_gen[i]
						) ;
				gluons.push_back (dummy1) ;
			}
		}

		cout<<" quarks size    "<<quarks.size()<<endl;
		cout<<" gluons size    "<<gluons.size()<<endl;
                cout<<" jets size    "<<jets.size()<<endl;

		sort (quarks.rbegin (), quarks.rend (), local_TLVP_ptSort ()) ;
		//cout<<" quarks size  after sorting   "<<quarks.size()<<endl;
		sort (gluons.rbegin (), gluons.rend (), local_TLVP_ptSort ()) ;
		//cout<<" gluons size  after sorting   "<<gluons.size()<<endl;
		//PG      deltaR      (recoIndex , genIndex)
		std::map <float, std::pair<int, int> >q_matchingsMap ;
		for (int iReco = 0 ; iReco < jets.size (); ++iReco) // loop over reco jets
		{
			for (int iGen = 0 ; iGen < quarks.size (); ++iGen) // loop over quarks
			{
				float deltaR = my_deltaR_function (
						jets.at (iReco).Eta (),jets.at (iReco).Phi (), 
						quarks.at (iGen).Eta (), quarks.at (iGen).Phi ()) ;
				h_qdeltaR.Fill (deltaR) ;
				q_matchingsMap[deltaR] = std::pair<int, int> (iReco, iGen) ;
			}
		}

				for (std::map<float, std::pair<int, int> >::const_iterator ipair = q_matchingsMap.begin () ;
				ipair != q_matchingsMap.end () ;
				++ipair)
				{
				//cout <<"delta R between jets and quarks : "<< ipair->first <<" reco jet index  "<<ipair->second.first<<" quark index  "<<ipair->second.second << endl;
				cout<<"......  "<<endl;
				}

		 

		vector<int> q_recoUsed (jets.size (), 0) ;   //PG vector of flags to know if a jet is already used
		vector<int> q_genUsed (quarks.size (), 0) ; //PG vector of flags to know if a jet is already used
		std::map<float, std::pair<int, int> > q_cleanedMatchings ;
		// remove elements with jets that have been already used
		for (std::map<float, std::pair<int, int> >::const_iterator iMap = q_matchingsMap.begin () ;
				iMap != q_matchingsMap.end () ;
				++iMap)
		{
			if (iMap->first > maxDR) break ;
			if (q_recoUsed.at (iMap->second.first) == 1 ||
					q_genUsed.at (iMap->second.second) == 1)
				continue ; 
			q_cleanedMatchings[iMap->first] = iMap->second ;
			q_recoUsed.at (iMap->second.first) = 1 ;
			q_genUsed.at (iMap->second.second) = 1 ;

		}  
		for (std::map<float, std::pair<int, int> >::const_iterator ipair1 = q_cleanedMatchings.begin () ;
				ipair1 != q_cleanedMatchings.end () ;
				++ipair1)
		{
			cout <<"delta R between clean jets and quarks : "<< ipair1->first <<" reco jet index  "<<ipair1->second.first<<" quark index  "<<ipair1->second.second << endl;
			//		h_quarklikejet1_pt.Fill (jets.at (ipair1->second.first ).Pt());
			// cout<<"......  "<<endl;
		}
		h_qmatchedSize.Fill (q_cleanedMatchings.size ()) ;

		//gluons matching
		std::map <float, std::pair<int, int> > g_matchingsMap ;
		for (int iReco = 0 ; iReco < jets.size (); ++iReco) // loop over reco jets
		{
			for (int iGen = 0 ; iGen < gluons.size (); ++iGen) // loop over gluons
			{
				float deltaR = my_deltaR_function (
						jets.at (iReco).Eta (),jets.at (iReco).Phi (),
						gluons.at (iGen).Eta (), gluons.at (iGen).Phi ()) ;
				h_gdeltaR.Fill (deltaR) ;
				g_matchingsMap[deltaR] = std::pair<int, int> (iReco, iGen) ;
			}
		}

		vector<int> g_recoUsed (jets.size (), 0) ;   //PG vector of flags to know if a jet is already used
		vector<int> g_genUsed (gluons.size (), 0) ; //PG vector of flags to know if a jet is already used
		std::map<float, std::pair<int, int> > g_cleanedMatchings ;
		// remove elements with jets that have been already used
		for (std::map<float, std::pair<int, int> >::const_iterator iMap = g_matchingsMap.begin () ;
				iMap != g_matchingsMap.end () ;
				++iMap)
		{
			if (iMap->first > maxDR) break ;
			if (g_recoUsed.at (iMap->second.first) == 1 ||
					g_genUsed.at (iMap->second.second) == 1)
				continue ;
			g_cleanedMatchings[iMap->first] = iMap->second ;
			g_recoUsed.at (iMap->second.first) = 1 ;
			g_genUsed.at (iMap->second.second) = 1 ;

		}
		h_gmatchedSize.Fill (g_cleanedMatchings.size ()) ;


		for (std::map<float, std::pair<int, int> >::const_iterator ipair2 = g_cleanedMatchings.begin () ;
				ipair2 != g_cleanedMatchings.end () ;
				++ipair2)
		{
			cout <<"delta R between clean jets and gluons : "<< ipair2->first <<" reco jet index  "<<ipair2->second.first<<" gluon index  "<<ipair2->second.second << endl;
			//		h_gluonlikejet1_pt.Fill (jets.at (ipair2->second.first ).Pt());
			// cout<<"......  "<<endl;
		}


		quarklike_jets.clear();

		for (std::map<float, std::pair<int, int> >::const_iterator iMap = q_cleanedMatchings.begin () ;
				iMap != q_cleanedMatchings.end () ;
				++iMap)
		{
			lorentzVector dummy3 ;
			dummy3.SetPxPyPzE (
					jets.at (iMap->second.first).Px (),
					jets.at (iMap->second.first).Py (),
					jets.at (iMap->second.first).Pz (),
					jets.at (iMap->second.first).E ()
					) ;
			//      quarklike_jets.clear();
			quarklike_jets.push_back(dummy3) ;
		}

		gluonlike_jets.clear();


		for (std::map<float, std::pair<int, int> >::const_iterator iMap = g_cleanedMatchings.begin () ;
				iMap != g_cleanedMatchings.end () ;
				++iMap)
		{
			lorentzVector dummy4 ;
			dummy4.SetPxPyPzE (
					jets.at (iMap->second.first).Px (),
					jets.at (iMap->second.first).Py (),
					jets.at (iMap->second.first).Pz (),
					jets.at (iMap->second.first).E ()
					) ;
			gluonlike_jets.push_back(dummy4) ;
		}

		sort (gluonlike_jets.rbegin (), gluonlike_jets.rend (), local_TLVP_ptSort ()) ;
		sort (quarklike_jets.rbegin (), quarklike_jets.rend (), local_TLVP_ptSort ()) ;

		cout<<"gluonlike_jets size "<<gluonlike_jets.size()<<endl;
		cout<<"quarklike_jets size "<<quarklike_jets.size()<<endl;

		if(quarks.size()==1)
		{
			h_quark1_pt.Fill (quarks.at (0).Pt());
			h_quark1_eta.Fill (quarks.at (0).Eta());
		}
		if(quarks.size()==2)
		{        h_quark1_pt.Fill (quarks.at (0).Pt());
			h_quark1_eta.Fill (quarks.at (0).Eta());
			h_quark2_pt.Fill (quarks.at (1).Pt());
			h_quark2_eta.Fill (quarks.at (1).Eta());
		}

		if(quarks.size()==3)
		{
			h_quark1_pt.Fill (quarks.at (0).Pt());
			h_quark1_eta.Fill (quarks.at (0).Eta());
			h_quark2_pt.Fill (quarks.at (1).Pt());
			h_quark2_eta.Fill (quarks.at (1).Eta());
			h_quark3_pt.Fill (quarks.at (2).Pt());
			h_quark3_eta.Fill (quarks.at (2).Eta());
		}

		if(quarks.size()==4)
		{
			h_quark1_pt.Fill (quarks.at (0).Pt());
			h_quark1_eta.Fill (quarks.at (0).Eta());
			h_quark2_pt.Fill (quarks.at (1).Pt());
			h_quark2_eta.Fill (quarks.at (1).Eta());
			h_quark3_pt.Fill (quarks.at (2).Pt());
			h_quark3_eta.Fill (quarks.at (2).Eta());
			h_quark4_pt.Fill (quarks.at (3).Pt());
			h_quark4_eta.Fill (quarks.at (3).Eta());

		}


		if(gluons.size()==1)
		{
			h_gluon1_pt.Fill (gluons.at (0).Pt());
			h_gluon1_eta.Fill (gluons.at (0).Eta());
		}
		if(gluons.size()==2)
		{        h_gluon1_pt.Fill (gluons.at (0).Pt());
			h_gluon1_eta.Fill (gluons.at (0).Eta());
			h_gluon2_pt.Fill (gluons.at (1).Pt());
			h_gluon2_eta.Fill (gluons.at (1).Eta());
		}

		if(gluons.size()==3)
		{
			h_gluon1_pt.Fill (gluons.at (0).Pt());
			h_gluon1_eta.Fill (gluons.at (0).Eta());
			h_gluon2_pt.Fill (gluons.at (1).Pt());
			h_gluon2_eta.Fill (gluons.at (1).Eta());
			h_gluon3_pt.Fill (gluons.at (2).Pt());
			h_gluon3_eta.Fill (gluons.at (2).Eta());
		}

		if(gluons.size()==4)
		{
			h_gluon1_pt.Fill (gluons.at (0).Pt());
			h_gluon1_eta.Fill (gluons.at (0).Eta());
			h_gluon2_pt.Fill (gluons.at (1).Pt());
			h_gluon2_eta.Fill (gluons.at (1).Eta());
			h_gluon3_pt.Fill (gluons.at (2).Pt());
			h_gluon3_eta.Fill (gluons.at (2).Eta());
			h_gluon4_pt.Fill (gluons.at (3).Pt());
			h_gluon4_eta.Fill (gluons.at (3).Eta());

		}





		if(quarklike_jets.size()==1)
		{
			h_quarklike_jet1_pt.Fill (quarklike_jets.at (0).Pt());
			h_quarklike_jet1_eta.Fill (quarklike_jets.at (0).Eta());
		}
		if(quarklike_jets.size()==2)
		{
			h_quarklike_jet1_pt.Fill (quarklike_jets.at (0).Pt());
			h_quarklike_jet1_eta.Fill (quarklike_jets.at (0).Eta());
			h_quarklike_jet2_pt.Fill (quarklike_jets.at (1).Pt());
			h_quarklike_jet2_eta.Fill (quarklike_jets.at (1).Eta());
		}
		if(quarklike_jets.size()==3)
		{
			h_quarklike_jet1_pt.Fill (quarklike_jets.at (0).Pt());
			h_quarklike_jet1_eta.Fill (quarklike_jets.at (0).Eta());
			h_quarklike_jet2_pt.Fill (quarklike_jets.at (1).Pt());
			h_quarklike_jet2_eta.Fill (quarklike_jets.at (1).Eta());
			h_quarklike_jet3_pt.Fill (quarklike_jets.at (2).Pt());
			h_quarklike_jet3_eta.Fill (quarklike_jets.at (2).Eta());
		}
		if(quarklike_jets.size()==4)
		{
			h_quarklike_jet1_pt.Fill (quarklike_jets.at (0).Pt());
			h_quarklike_jet1_eta.Fill (quarklike_jets.at (0).Eta());
			h_quarklike_jet2_pt.Fill (quarklike_jets.at (1).Pt());
			h_quarklike_jet2_eta.Fill (quarklike_jets.at (1).Eta());
			h_quarklike_jet3_pt.Fill (quarklike_jets.at (2).Pt());
			h_quarklike_jet3_eta.Fill (quarklike_jets.at (2).Eta());
			h_quarklike_jet4_pt.Fill (quarklike_jets.at (3).Pt());
			h_quarklike_jet4_eta.Fill (quarklike_jets.at (3).Eta());
		}


		//..................................
		if(gluonlike_jets.size()==1)
		{
			h_gluonlike_jet1_pt.Fill (gluonlike_jets.at (0).Pt());
			h_gluonlike_jet1_eta.Fill (gluonlike_jets.at (0).Eta());
		}
		if(gluonlike_jets.size()==2)
		{
			h_gluonlike_jet1_pt.Fill (gluonlike_jets.at (0).Pt());
			h_gluonlike_jet1_eta.Fill (gluonlike_jets.at (0).Eta());
			h_gluonlike_jet2_pt.Fill (gluonlike_jets.at (1).Pt());
			h_gluonlike_jet2_eta.Fill (gluonlike_jets.at (1).Eta());
		}
		if(gluonlike_jets.size()==3)
		{
			h_gluonlike_jet1_pt.Fill (gluonlike_jets.at (0).Pt());
			h_gluonlike_jet1_eta.Fill (gluonlike_jets.at (0).Eta());
			h_gluonlike_jet2_pt.Fill (gluonlike_jets.at (1).Pt());
			h_gluonlike_jet2_eta.Fill (gluonlike_jets.at (1).Eta());
			h_gluonlike_jet3_pt.Fill (gluonlike_jets.at (2).Pt());
			h_gluonlike_jet3_eta.Fill (gluonlike_jets.at (2).Eta());
		}
		if(gluonlike_jets.size()==4)
		{
			h_gluonlike_jet1_pt.Fill (gluonlike_jets.at (0).Pt());
			h_gluonlike_jet1_eta.Fill (gluonlike_jets.at (0).Eta());
			h_gluonlike_jet2_pt.Fill (gluonlike_jets.at (1).Pt());
			h_gluonlike_jet2_eta.Fill (gluonlike_jets.at (1).Eta());
			h_gluonlike_jet3_pt.Fill (gluonlike_jets.at (2).Pt());
			h_gluonlike_jet3_eta.Fill (gluonlike_jets.at (2).Eta());
			h_gluonlike_jet4_pt.Fill (gluonlike_jets.at (3).Pt());
			h_gluonlike_jet4_eta.Fill (gluonlike_jets.at (3).Eta());
		}





	}   //PG loop over events

	TString outFileName = "output_" ;
	outFileName += maxNumber;
	outFileName += "_";
	outFileName += int (maxDR * 100);
	outFileName += ".root" ;
	TFile outFile (outFileName, "recreate") ;
	h_qmatchedSize.Write () ;
	h_gmatchedSize.Write () ;

	h_qdeltaR.Write () ;
	h_gdeltaR.Write () ;


	h_quarklike_jet1_pt.Write ();
	h_quarklike_jet1_eta.Write ();
	h_quarklike_jet2_pt.Write ();
	h_quarklike_jet2_eta.Write ();
	h_quarklike_jet3_pt.Write ();
	h_quarklike_jet3_eta.Write ();
	h_quarklike_jet4_pt.Write ();
	h_quarklike_jet4_eta.Write ();

	h_gluonlike_jet1_pt.Write ();
	h_gluonlike_jet1_eta.Write ();
	h_gluonlike_jet2_pt.Write ();
	h_gluonlike_jet2_eta.Write ();
	h_gluonlike_jet3_pt.Write ();
	h_gluonlike_jet3_eta.Write ();
	h_gluonlike_jet4_pt.Write ();
	h_gluonlike_jet4_eta.Write ();


	h_quark1_pt.Write ();
	h_quark1_eta.Write ();
	h_quark2_pt.Write ();
	h_quark2_eta.Write ();
	h_quark3_pt.Write ();
	h_quark3_eta.Write ();
	h_quark4_pt.Write ();
	h_quark4_eta.Write ();

	h_gluon1_pt.Write ();
	h_gluon1_eta.Write ();
	h_gluon2_pt.Write ();
	h_gluon2_eta.Write ();
	h_gluon3_pt.Write ();
	h_gluon3_eta.Write ();
	h_gluon4_pt.Write ();
	h_gluon4_eta.Write ();



	outFile.Close () ;

	return 0 ;
}

