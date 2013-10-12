void plots_mmmm4_eta_wjets_tagjets(){

  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",800,800); 

  c1->Divide(2,2);


 TFile *f1 = new TFile("output_method1_50.root"); 
 TFile *f2 = new TFile("output_method2_50.root");
 TFile *f3 = new TFile("output_method3_50.root");
 TFile *f4 = new TFile("output_method4_50.root");
        f1.ls();
        f2.ls();
        f3.ls();
        f4.ls();
   TH1F *h1_1 = (TH1F*)f1->Get("h_vbf_aj_eta");
   TH1F *h2_1 = (TH1F*)f1->Get("h_vbf_bj_eta");
   TH1F *h3_1 = (TH1F*)f1->Get("h_vbf_waj_eta");
   TH1F *h4_1 = (TH1F*)f1->Get("h_vbf_wbj_eta");
  
  
  TH1F *h1_2 = (TH1F*)f2->Get("h_vbf_aj_eta");
  TH1F *h2_2 = (TH1F*)f2->Get("h_vbf_bj_eta");
  TH1F *h3_2 = (TH1F*)f2->Get("h_vbf_waj_eta");
  TH1F *h4_2 = (TH1F*)f2->Get("h_vbf_wbj_eta");
  
  TH1F *h1_3 = (TH1F*)f3->Get("h_vbf_aj_eta");
  TH1F *h2_3 = (TH1F*)f3->Get("h_vbf_bj_eta");
  TH1F *h3_3 = (TH1F*)f3->Get("h_vbf_waj_eta");
  TH1F *h4_3 = (TH1F*)f3->Get("h_vbf_wbj_eta");
  
   TH1F *h1_4 = (TH1F*)f4->Get("h_vbf_aj_eta");
  TH1F *h2_4 = (TH1F*)f4->Get("h_vbf_bj_eta");
  TH1F *h3_4 = (TH1F*)f4->Get("h_vbf_waj_eta");
  TH1F *h4_4 = (TH1F*)f4->Get("h_vbf_wbj_eta");
 

 
 c1->cd(1); 

//  TLegend *tleg1_1 = new TLegend(0.30,0.65,0.55,0.85);

  TAxis *x = h1_1->GetXaxis();
 // x->SetTitle("Number of tracks");
 // x->SetTitleColor(kBlue);
//  x->CenterTitle(true);
  TAxis *y = h1_1->GetYaxis();
 // y->SetTitle("Number of events");
//  y->SetTitleColor(kRed);
//  y->CenterTitle(true);
//  y->SetTitleOffset(1.2);
 //  h1_1->SetMinimum(0);
 // h1_1->SetLineColor(1);
// h1_1->SetMarkerColor(1);
//  h1_1->SetMarkerSize(1);
 // h1_1->SetMarkerStyle(29);


  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(1);
  h1_1->SetStats(kFALSE);
  h1_1->Draw("");
 
  h1_2->SetLineColor(kRed);
  h1_2->SetLineWidth(1);
  h1_2->SetLineStyle(2);
  h1_2->SetStats(kFALSE);
  h1_2->Draw("sames");

  h1_3->SetLineColor(kGreen);
  h1_3->SetLineWidth(1);
  h1_3->SetLineStyle(2);
  h1_3->SetStats(kFALSE);
  h1_3->Draw("sames");

  h1_4->SetLineColor(kBlue);
  h1_4->SetLineWidth(1);
  h1_4->SetLineStyle(2);
  h1_4->SetStats(kFALSE);
  h1_4->Draw("sames");
 

gPad->Modified(); gPad->Update();

 
TLegend *tleg1_1 = new TLegend(0.1108809,0.6774907,0.4610479,0.8898976,NULL,"brNDC");
tleg1_1 ->SetTextSize(0.043);
tleg1_1->SetBorderSize(1);
tleg1_1->SetLineColor(1);
tleg1_1->SetLineStyle(1);
tleg1_1->SetLineWidth(1);
tleg1_1->SetFillColor(19);
tleg1_1->SetFillStyle(1001);
  tleg1_1->AddEntry(h1_1,"method1");
  tleg1_1->AddEntry(h1_2,"method2");
  tleg1_1->AddEntry(h1_3,"method3 ");
  tleg1_1->AddEntry(h1_4,"method4 ");
  tleg1_1->Draw();

 
 c1->cd(2); 
  
  TAxis *x = h2_1->GetXaxis();
 // x->SetTitle("Number of tracks");
 // x->SetTitleColor(kBlue);
//  x->CenterTitle(true);
  TAxis *y = h2_1->GetYaxis();
 // y->SetTitle("Number of events");
//  y->SetTitleColor(kRed);
//  y->CenterTitle(true);
 // y->SetTitleOffset(1.2);
 // h2_1->SetMinimum(0);
 // h2_1->SetMarkerColor(1);
 // h2_1->SetMarkerSize(1);
 // h2_1->SetMarkerStyle(29);

  h2_1->SetLineWidth(1);
  h2_1->SetLineColor(kBlack);
  h2_1->SetStats(kFALSE);
  h2_1->Draw("");

  h2_2->SetLineColor(kRed);
  h2_2->SetLineWidth(1);
  h2_2->SetLineStyle(2);
  h2_2->SetStats(kFALSE);
  h2_2->Draw("sames");

  h2_3->SetLineColor(kGreen);
  h2_3->SetLineWidth(1);
  h2_3->SetLineStyle(2);
  h2_3->SetStats(kFALSE);
  h2_3->Draw("sames");

  h2_4->SetLineColor(kBlue);
  h2_4->SetLineWidth(1);
  h2_4->SetLineStyle(2);
  h2_4->SetStats(kFALSE);
  h2_4->Draw("sames");

gPad->Modified(); gPad->Update();



TLegend *tleg2_1 = new TLegend(0.5550148,0.2003026,0.8923317,0.4243482,NULL,"brNDC");
  tleg2_1 ->SetTextSize(0.043);
 tleg2_1->SetBorderSize(1);
tleg2_1->SetLineColor(1);
tleg2_1->SetLineStyle(1);
tleg2_1->SetLineWidth(1);
tleg2_1->SetFillColor(19);
tleg2_1->SetFillStyle(1001);
  tleg2_1->AddEntry(h2_1,"method1");
  tleg2_1->AddEntry(h2_2,"method2");
    tleg2_1->AddEntry(h2_3,"method3");
    tleg2_1->AddEntry(h2_4,"method4");
  tleg2_1->Draw();

//******************************
//***************************
 c1->cd(3); 
  
  TAxis *x = h3_1->GetXaxis();
 // x->SetTitle("Number of tracks");
 // x->SetTitleColor(kBlue);
 // x->CenterTitle(true);
  TAxis *y = h3_1->GetYaxis();
 // y->SetTitle("Number of events");
//  y->SetTitleColor(kRed);
 // y->CenterTitle(true);
 // y->SetTitleOffset(1.2);
 // h2_1->SetMinimum(0);
 // h2_1->SetMarkerColor(1);
 // h2_1->SetMarkerSize(1);
 // h2_1->SetMarkerStyle(29);

  h3_1->SetLineWidth(1);
  h3_1->SetLineColor(kBlack);
  h3_1->SetStats(kFALSE);
  h3_1->Draw("");

  h3_2->SetLineColor(kRed);
  h3_2->SetLineWidth(1);
  h3_2->SetLineStyle(2);
  h3_2->SetStats(kFALSE);
  h3_2->Draw("sames");

  h3_3->SetLineColor(kGreen);
  h3_3->SetLineWidth(1);
  h3_3->SetLineStyle(2);
  h3_3->SetStats(kFALSE);
  h3_3->Draw("sames");

  h3_4->SetLineColor(kBlue);
  h3_4->SetLineWidth(1);
  h3_4->SetLineStyle(2);
  h3_4->SetStats(kFALSE);
  h3_4->Draw("sames");

gPad->Modified(); gPad->Update();



TLegend *tleg3_1 = new TLegend(0.2458076,0.6571229,0.5831245,0.8811685,NULL,"brNDC");
 tleg3_1 ->SetTextSize(0.043);
tleg3_1->SetBorderSize(1);
tleg3_1->SetLineColor(1);
tleg3_1->SetLineStyle(1);
tleg3_1->SetLineWidth(1);
tleg3_1->SetFillColor(19);
tleg3_1->SetFillStyle(1001);
 tleg3_1->AddEntry(h2_1,"method1");
  tleg3_1->AddEntry(h2_2,"method2");
  tleg3_1->AddEntry(h2_3,"method3");
  tleg3_1->AddEntry(h2_4,"method4");
 // tleg2_1->AddEntry(h2_2,"May16ReReco-42X_158028-158383");
  tleg2_1->Draw();

 
//******************************
//***************************
 
 c1->cd(4); 
 
  TAxis *x = h4_1->GetXaxis();
 // x->SetTitle("Number of tracks");
 // x->SetTitleColor(kBlue);
 // x->CenterTitle(true);
  TAxis *y = h4_1->GetYaxis();
 // y->SetTitle("Number of events");
//  y->SetTitleColor(kRed);
 // y->CenterTitle(true);
 // y->SetTitleOffset(1.2);
 // h4_1->SetMinimum(0);
 // h4_1->SetMarkerColor(1);
 // h4_1->SetMarkerSize(1);
 // h4_1->SetMarkerStyle(29);

  h4_1->SetLineWidth(1);
  h4_1->SetLineColor(kBlack);
  h4_1->SetStats(kFALSE);
  h4_1->Draw("");  
  
  h4_2->SetLineWidth(1);
  h4_2->SetLineStyle(2);
  h4_2->SetLineColor(kRed);
  h4_2->SetStats(kFALSE);
  h4_2->Draw("sames"); 

  h4_3->SetLineWidth(1);
  h4_3->SetLineStyle(2);
  h4_3->SetLineColor(kGreen);
  h4_3->SetStats(kFALSE);
  h4_3->Draw("sames"); 
   
   h4_4->SetLineWidth(1);
  h4_4->SetLineStyle(2);
  h4_4->SetLineColor(kBlue);
  h4_4->SetStats(kFALSE);
  h4_4->Draw("sames"); 
   
 
gPad->Modified(); gPad->Update();


 TLegend *tleg4_1 = new TLegend(0.2305481,0.665852,0.5550148,0.8811685,NULL,"brNDC");
tleg4_1 ->SetTextSize(0.043);
tleg4_1->SetBorderSize(1);
tleg4_1->SetLineColor(1);
tleg4_1->SetLineStyle(1);
tleg4_1->SetLineWidth(1);
tleg4_1->SetFillColor(19);
tleg4_1->SetFillStyle(1001);
  tleg4_1->AddEntry(h4_1,"method1");
  tleg4_1->AddEntry(h4_2,"method2");
  tleg4_1->AddEntry(h4_3,"method3");
  tleg4_1->AddEntry(h4_4,"method4");
   tleg4_1->Draw();

//******************************
 c1->SaveAs("mmmm4_etas.png");
}

