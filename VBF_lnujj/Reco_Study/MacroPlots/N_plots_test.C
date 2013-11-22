void N_plots_test(){

  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",800,800); 

 // c1->Divide(2,2);


 TFile f1 ("Tree_W4Jets_vbf.root"); 
        //f1.ls();

TTree *T1 = (TTree*)f1.Get("QGl");

TH1F *h1_1 = new TH1F("h1_1","qgl",100,0,600);
 T1->Draw("h_quarklike_pt[0]>>h1_1","(h_quarklike_pt[0]>30 && h_quarklike_pt[0]>-1)" );

  //TAxis *x = h1_1->GetXaxis();
  //TAxis *y = h1_1->GetYaxis();
  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(2);
  h1_1->SetLineStyle(1);
  h1_1->SetStats(kFALSE);
  Double_t scaleh11 = 1./h1_1->Integral();
  h1_1->Scale(scaleh11);
 // h1_1->Draw();

TH1F *h2_1 = new TH1F("h2_1","qgl",100,0,600);
 T1->Draw("h_gluonlike_pt[0]>>h2_1","(h_gluonlike_pt[0]>30 && h_gluonlike_pt[0]>-1)" );
  h2_1->SetLineColor(kRed);
  h2_1->SetLineWidth(2);
  h2_1->SetLineStyle(2);
  h2_1->SetStats(kFALSE);
  Double_t scaleh12 = 1./h2_1->Integral();
  h2_1->Scale(scaleh12);
  h2_1->Draw();
  h1_1->Draw("sames");

  

gPad->Modified(); gPad->Update();

/*
 TLegend *tleg4_1 = new TLegend(0.2305481,0.665852,0.5550148,0.8811685,NULL,"brNDC");
tleg4_1 ->SetTextSize(0.043);
tleg4_1->SetBorderSize(1);
tleg4_1->SetLineColor(1);
tleg4_1->SetLineStyle(1);
tleg4_1->SetLineWidth(2);
tleg4_1->SetFillColor(19);
tleg4_1->SetFillStyle(1001);
  tleg4_1->AddEntry(h4_1,"Madgraph");
  tleg4_1->AddEntry(h4_2,"vbf@nlo");
  tleg4_1->AddEntry(h4_3,"phantom");

   tleg4_1->Draw();
*/
//******************************
 c1->SaveAs("W4Jets_newaNormalizedtoArea_pt_dis[0]30.png");
}

