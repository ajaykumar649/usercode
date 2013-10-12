void N_plots_mmmm4_eta_wjets_tagjets(){

  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",800,800); 

  c1->Divide(2,2);


 TFile *f1 = new TFile("VBF300_output_method1_50.root"); 
 TFile *f2 = new TFile("VBF300_output_method2_50.root");
 TFile *f3 = new TFile("VBF300_output_method3_50.root");
 TFile *f4 = new TFile("VBF300_output_method4_50.root");
 TFile *f5 = new TFile("VBF300_output_method5_50.root");
 TFile *f6 = new TFile("VBF300_output_method6_50.root");
 TFile *f7 = new TFile("VBF300_output_method7_50.root");
 TFile *f8 = new TFile("VBF300_output_method8_50.root");

        f1.ls();
        f2.ls();
        f3.ls();
        f4.ls();
        f5.ls();
        f6.ls();
        f7.ls();
        f8.ls();

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
 
  TH1F *h1_5 = (TH1F*)f5->Get("h_vbf_aj_eta");
  TH1F *h2_5 = (TH1F*)f5->Get("h_vbf_bj_eta");
  TH1F *h3_5 = (TH1F*)f5->Get("h_vbf_waj_eta");
  TH1F *h4_5 = (TH1F*)f5->Get("h_vbf_wbj_eta");


  TH1F *h1_6 = (TH1F*)f6->Get("h_vbf_aj_eta");
  TH1F *h2_6 = (TH1F*)f6->Get("h_vbf_bj_eta");
  TH1F *h3_6 = (TH1F*)f6->Get("h_vbf_waj_eta");
  TH1F *h4_6 = (TH1F*)f6->Get("h_vbf_wbj_eta");

  TH1F *h1_7 = (TH1F*)f7->Get("h_vbf_aj_eta");
  TH1F *h2_7 = (TH1F*)f7->Get("h_vbf_bj_eta");
  TH1F *h3_7 = (TH1F*)f7->Get("h_vbf_waj_eta");
  TH1F *h4_7 = (TH1F*)f7->Get("h_vbf_wbj_eta");

  TH1F *h1_8 = (TH1F*)f8->Get("h_vbf_aj_eta");
  TH1F *h2_8 = (TH1F*)f8->Get("h_vbf_bj_eta");
  TH1F *h3_8 = (TH1F*)f8->Get("h_vbf_waj_eta");
  TH1F *h4_8 = (TH1F*)f8->Get("h_vbf_wbj_eta");


 
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

  h1_6->SetLineColor(kOrange);
  h1_6->SetLineWidth(2);
  h1_6->SetLineStyle(2);
  h1_6->SetStats(kFALSE);
  Double_t scaleh16 = 1./h1_6->Integral();
  h1_6->Scale(scaleh16);
  h1_6->Draw("");


  h1_8->SetLineColor(kViolet);
  h1_8->SetLineWidth(2);
  h1_8->SetLineStyle(2);
  h1_8->SetStats(kFALSE);
  Double_t scaleh18 = 1./h1_8->Integral();
  h1_8->Scale(scaleh18);
  h1_8->Draw("sames");


  h1_3->SetLineColor(kGreen);
  h1_3->SetLineWidth(2);
  h1_3->SetLineStyle(2);
  h1_3->SetStats(kFALSE);
  Double_t scaleh13 = 1./h1_3->Integral();
  h1_3->Scale(scaleh13);
  h1_3->Draw("sames");

  h1_4->SetLineColor(kBlue);
  h1_4->SetLineWidth(2);
  h1_4->SetLineStyle(2);
  h1_4->SetStats(kFALSE);
  Double_t scaleh14 = 1./h1_4->Integral();
  h1_4->Scale(scaleh14);
  h1_4->Draw("sames");

  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(2);
  h1_1->SetStats(kFALSE);
  Double_t scaleh11 = 1./h1_1->Integral();
  h1_1->Scale(scaleh11);
  h1_1->Draw("sames");
 
  h1_2->SetLineColor(kRed);
  h1_2->SetLineWidth(2);
  h1_2->SetLineStyle(2);
  h1_2->SetStats(kFALSE);
  Double_t scaleh12 = 1./h1_2->Integral();
  h1_2->Scale(scaleh12);
  h1_2->Draw("sames");

  h1_5->SetLineColor(kYellow);
  h1_5->SetLineWidth(2);
  h1_5->SetLineStyle(2);
  h1_5->SetStats(kFALSE);
  Double_t scaleh15 = 1./h1_5->Integral();
  h1_5->Scale(scaleh15);
  h1_5->Draw("sames");
/*
  h1_6->SetLineColor(kOrange);
  h1_6->SetLineWidth(2);
  h1_6->SetLineStyle(2);
  h1_6->SetStats(kFALSE);
  Double_t scaleh16 = 1./h1_6->Integral();
  h1_6->Scale(scaleh16);
  h1_6->Draw("sames");
*/
  h1_7->SetLineColor(kCyan);
  h1_7->SetLineWidth(2);
  h1_7->SetLineStyle(2);
  h1_7->SetStats(kFALSE);
  Double_t scaleh17 = 1./h1_7->Integral();
  h1_7->Scale(scaleh17);
  h1_7->Draw("sames");
/*
  h1_8->SetLineColor(kViolet);
  h1_8->SetLineWidth(2);
  h1_8->SetLineStyle(2);
  h1_8->SetStats(kFALSE);
  Double_t scaleh18 = 1./h1_8->Integral();
  h1_8->Scale(scaleh18);
  h1_8->Draw("sames");
*/
/*
  h1_3->SetLineColor(kGreen);
  h1_3->SetLineWidth(2);
  h1_3->SetLineStyle(2);
  h1_3->SetStats(kFALSE);
  Double_t scaleh13 = 1./h1_3->Integral();
  h1_3->Scale(scaleh13);
  h1_3->Draw("sames");
*/
/*
  h1_4->SetLineColor(kBlue);
  h1_4->SetLineWidth(2);
  h1_4->SetLineStyle(2);
  h1_4->SetStats(kFALSE);
  Double_t scaleh14 = 1./h1_4->Integral();
  h1_4->Scale(scaleh14);
  h1_4->Draw("sames");
 */

gPad->Modified(); gPad->Update();

 
TLegend *tleg1_1 = new TLegend(0.1108809,0.6774907,0.4610479,0.8898976,NULL,"brNDC");
//TLegend *tleg1_1 = new TLegend(0.01,.01,0.91,0.91,NULL,"brNDC");
tleg1_1 ->SetTextSize(0.043);
tleg1_1->SetBorderSize(1);
tleg1_1->SetLineColor(1);
tleg1_1->SetLineStyle(1);
tleg1_1->SetLineWidth(2);
tleg1_1->SetFillColor(19);
tleg1_1->SetFillStyle(1001);
  tleg1_1->AddEntry(h1_1,"method1");
  tleg1_1->AddEntry(h1_2,"method2");
  tleg1_1->AddEntry(h1_3,"method3 ");
  tleg1_1->AddEntry(h1_4,"method4 ");
  tleg1_1->AddEntry(h1_5,"method5 ");
  tleg1_1->AddEntry(h1_6,"method6 ");
  tleg1_1->AddEntry(h1_7,"method7 ");
  tleg1_1->AddEntry(h1_8,"method8 ");

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

  h2_6->SetLineColor(kOrange);
  h2_6->SetLineWidth(2);
  h2_6->SetLineStyle(2);
  h2_6->SetStats(kFALSE);
  Double_t scaleh26 = 1./h2_6->Integral();
  h2_6->Scale(scaleh26);
  h2_6->Draw("");


  h2_8->SetLineColor(kViolet);
  h2_8->SetLineWidth(2);
  h2_8->SetLineStyle(2);
  h2_8->SetStats(kFALSE);
  Double_t scaleh28 = 1./h2_8->Integral();
  h2_8->Scale(scaleh28);
  h2_8->Draw("sames");


  h2_3->SetLineColor(kGreen);
  h2_3->SetLineWidth(2);
  h2_3->SetLineStyle(2);
  h2_3->SetStats(kFALSE);
  Double_t scaleh23 = 1./h2_3->Integral();
  h2_3->Scale(scaleh23);
  h2_3->Draw("sames");


  h2_4->SetLineColor(kBlue);
  h2_4->SetLineWidth(2);
  h2_4->SetLineStyle(2);
  h2_4->SetStats(kFALSE);
  Double_t scaleh24 = 1./h2_4->Integral();
  h2_4->Scale(scaleh24);
  h2_4->Draw("sames");


  h2_1->SetLineWidth(2);
  h2_1->SetLineColor(kBlack);
  h2_1->SetStats(kFALSE);
  Double_t scaleh21 = 1./h2_1->Integral();
  h2_1->Scale(scaleh21);
  h2_1->Draw("sames");

  h2_2->SetLineColor(kRed);
  h2_2->SetLineWidth(2);
  h2_2->SetLineStyle(2);
  h2_2->SetStats(kFALSE);
  Double_t scaleh22 = 1./h2_2->Integral();
  h2_2->Scale(scaleh22);
  h2_2->Draw("sames");

  h2_3->SetLineColor(kGreen);
  h2_3->SetLineWidth(2);
  h2_3->SetLineStyle(2);
  h2_3->SetStats(kFALSE);
  Double_t scaleh23 = 1./h2_3->Integral();
  h2_3->Scale(scaleh23);
  h2_3->Draw("sames");

  h2_5->SetLineColor(kYellow);
  h2_5->SetLineWidth(2);
  h2_5->SetLineStyle(2);
  h2_5->SetStats(kFALSE);
  Double_t scaleh25 = 1./h2_5->Integral();
  h2_5->Scale(scaleh25);
  h2_5->Draw("sames");
/*
  h2_6->SetLineColor(kOrange);
  h2_6->SetLineWidth(2);
  h2_6->SetLineStyle(2);
  h2_6->SetStats(kFALSE);
  Double_t scaleh26 = 1./h2_6->Integral();
  h2_6->Scale(scaleh26);
  h2_6->Draw("sames");
*/
  h2_7->SetLineColor(kCyan);
  h2_7->SetLineWidth(2);
  h2_7->SetLineStyle(2);
  h2_7->SetStats(kFALSE);
  Double_t scaleh27 = 1./h2_7->Integral();
  h2_7->Scale(scaleh27);
  h2_7->Draw("sames");
/*
  h2_8->SetLineColor(kViolet);
  h2_8->SetLineWidth(2);
  h2_8->SetLineStyle(2);
  h2_8->SetStats(kFALSE);
  Double_t scaleh28 = 1./h2_8->Integral();
  h2_8->Scale(scaleh28);
  h2_8->Draw("sames");
*/

gPad->Modified(); gPad->Update();



TLegend *tleg2_1 = new TLegend(0.5550148,0.2003026,0.8923317,0.4243482,NULL,"brNDC");
  tleg2_1 ->SetTextSize(0.043);
 tleg2_1->SetBorderSize(1);
tleg2_1->SetLineColor(1);
tleg2_1->SetLineStyle(1);
tleg2_1->SetLineWidth(2);
tleg2_1->SetFillColor(19);
tleg2_1->SetFillStyle(1001);
    tleg2_1->AddEntry(h2_1,"method1");
    tleg2_1->AddEntry(h2_2,"method2");
    tleg2_1->AddEntry(h2_3,"method3");
    tleg2_1->AddEntry(h2_4,"method4");
    tleg2_1->AddEntry(h2_5,"method5");
    tleg2_1->AddEntry(h2_6,"method6");
    tleg2_1->AddEntry(h2_7,"method7");
    tleg2_1->AddEntry(h2_8,"method8");

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

  h3_3->SetLineColor(kGreen);
  h3_3->SetLineWidth(2);
  h3_3->SetLineStyle(2);
  h3_3->SetStats(kFALSE);
  Double_t scaleh33 = 1./h3_3->Integral();
  h3_3->Scale(scaleh33);
  h3_3->Draw("");


  h3_8->SetLineColor(kViolet);
  h3_8->SetLineWidth(2);
  h3_8->SetLineStyle(2);
  h3_8->SetStats(kFALSE);
  Double_t scaleh38 = 1./h3_8->Integral();
  h3_8->Scale(scaleh38);
  h3_8->Draw("sames");

  h3_4->SetLineColor(kBlue);
  h3_4->SetLineWidth(2);
  h3_4->SetLineStyle(2);
  h3_4->SetStats(kFALSE);
  Double_t scaleh34 = 1./h3_4->Integral();
  h3_4->Scale(scaleh34);
  h3_4->Draw("sames");
/*
  h3_3->SetLineColor(kGreen);
  h3_3->SetLineWidth(2);
  h3_3->SetLineStyle(2);
  h3_3->SetStats(kFALSE);
  Double_t scaleh33 = 1./h3_3->Integral();
  h3_3->Scale(scaleh33);
  h3_3->Draw("sames");
*/

  h3_1->SetLineWidth(2);
  h3_1->SetLineColor(kBlack);
  h3_1->SetStats(kFALSE);
  Double_t scaleh31 = 1./h3_1->Integral();
  h3_1->Scale(scaleh31);
  h3_1->Draw("sames");

  h3_2->SetLineColor(kRed);
  h3_2->SetLineWidth(2);
  h3_2->SetLineStyle(2);
  h3_2->SetStats(kFALSE);
  Double_t scaleh32 = 1./h3_2->Integral();
  h3_2->Scale(scaleh32);
  h3_2->Draw("sames");

  h3_5->SetLineColor(kYellow);
  h3_5->SetLineWidth(2);
  h3_5->SetLineStyle(2);
  h3_5->SetStats(kFALSE);
  Double_t scaleh35 = 1./h3_5->Integral();
  h3_5->Scale(scaleh35);
  h3_5->Draw("sames");

  h3_6->SetLineColor(kOrange);
  h3_6->SetLineWidth(2);
  h3_6->SetLineStyle(2);
  h3_6->SetStats(kFALSE);
  Double_t scaleh36 = 1./h3_6->Integral();
  h3_6->Scale(scaleh36);
  h3_6->Draw("sames");

  h3_7->SetLineColor(kCyan);
  h3_7->SetLineWidth(2);
  h3_7->SetLineStyle(2);
  h3_7->SetStats(kFALSE);
  Double_t scaleh37 = 1./h3_7->Integral();
  h3_7->Scale(scaleh37);
  h3_7->Draw("sames");
/*
  h3_8->SetLineColor(kViolet);
  h3_8->SetLineWidth(2);
  h3_8->SetLineStyle(2);
  h3_8->SetStats(kFALSE);
  Double_t scaleh38 = 1./h3_8->Integral();
  h3_8->Scale(scaleh38);
  h3_8->Draw("sames");
*/
/*
  h3_3->SetLineColor(kGreen);
  h3_3->SetLineWidth(2);
  h3_3->SetLineStyle(2);
  h3_3->SetStats(kFALSE);
  Double_t scaleh33 = 1./h3_3->Integral();
  h3_3->Scale(scaleh33);
  h3_3->Draw("sames");
*/
/*  h3_4->SetLineColor(kBlue);
  h3_4->SetLineWidth(2);
  h3_4->SetLineStyle(2);
  h3_4->SetStats(kFALSE);
  Double_t scaleh34 = 1./h3_4->Integral();
  h3_4->Scale(scaleh34);
  h3_4->Draw("sames");
*/
gPad->Modified(); gPad->Update();



TLegend *tleg3_1 = new TLegend(0.2458076,0.6571229,0.5831245,0.8811685,NULL,"brNDC");
 tleg3_1 ->SetTextSize(0.043);
tleg3_1->SetBorderSize(1);
tleg3_1->SetLineColor(1);
tleg3_1->SetLineStyle(1);
tleg3_1->SetLineWidth(2);
tleg3_1->SetFillColor(19);
tleg3_1->SetFillStyle(1001);
  tleg3_1->AddEntry(h2_1,"method1");
  tleg3_1->AddEntry(h2_2,"method2");
  tleg3_1->AddEntry(h2_3,"method3");
  tleg3_1->AddEntry(h2_4,"method4");
  tleg3_1->AddEntry(h2_5,"method5");
  tleg3_1->AddEntry(h2_6,"method6");
  tleg3_1->AddEntry(h2_7,"method7");
  tleg3_1->AddEntry(h2_8,"method8");

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

  h4_3->SetLineWidth(2);
  h4_3->SetLineStyle(2);
  h4_3->SetLineColor(kGreen);
  h4_3->SetStats(kFALSE);
  Double_t scaleh43 = 1./h4_3->Integral();
  h4_3->Scale(scaleh43);
  h4_3->Draw("");


  h4_2->SetLineWidth(2);
  h4_2->SetLineStyle(2);
  h4_2->SetLineColor(kRed);
  h4_2->SetStats(kFALSE);
  Double_t scaleh42 = 1./h4_2->Integral();
  h4_2->Scale(scaleh42);
  h4_2->Draw("sames");

/*
  h4_3->SetLineWidth(2);
  h4_3->SetLineStyle(2);
  h4_3->SetLineColor(kGreen);
  h4_3->SetStats(kFALSE);
  Double_t scaleh43 = 1./h4_3->Integral();
  h4_3->Scale(scaleh43);
  h4_3->Draw("sames"); 
*/

  h4_1->SetLineWidth(2);
  h4_1->SetLineColor(kBlack);
  h4_1->SetStats(kFALSE);
  Double_t scaleh41 = 1./h4_1->Integral();
  h4_1->Scale(scaleh41);
  h4_1->Draw("sames");  
/*  
  h4_2->SetLineWidth(2);
  h4_2->SetLineStyle(2);
  h4_2->SetLineColor(kRed);
  h4_2->SetStats(kFALSE);
  Double_t scaleh42 = 1./h4_2->Integral();
  h4_2->Scale(scaleh42);
  h4_2->Draw("sames"); 
*/
/*
  h4_3->SetLineWidth(2);
  h4_3->SetLineStyle(2);
  h4_3->SetLineColor(kGreen);
  h4_3->SetStats(kFALSE);
  Double_t scaleh43 = 1./h4_3->Integral();
  h4_3->Scale(scaleh43);
  h4_3->Draw("sames"); 
*/   
  h4_4->SetLineWidth(2);
  h4_4->SetLineStyle(2);
  h4_4->SetLineColor(kBlue);
  h4_4->SetStats(kFALSE);
  Double_t scaleh44 = 1./h4_4->Integral();
  h4_4->Scale(scaleh44);
  h4_4->Draw("sames"); 
   
  h4_5->SetLineWidth(2);
  h4_5->SetLineStyle(2);
  h4_5->SetLineColor(kYellow);
  h4_5->SetStats(kFALSE);
  Double_t scaleh45 = 1./h4_5->Integral();
  h4_5->Scale(scaleh45);
  h4_5->Draw("sames");

  h4_6->SetLineWidth(2);
  h4_6->SetLineStyle(2);
  h4_6->SetLineColor(kOrange);
  h4_6->SetStats(kFALSE);
  Double_t scaleh46 = 1./h4_6->Integral();
  h4_6->Scale(scaleh46);
  h4_6->Draw("sames");

  h4_7->SetLineWidth(2);
  h4_7->SetLineStyle(2);
  h4_7->SetLineColor(kCyan);
  h4_7->SetStats(kFALSE);
  Double_t scaleh47 = 1./h4_7->Integral();
  h4_7->Scale(scaleh47);
  h4_7->Draw("sames");

  h4_8->SetLineWidth(2);
  h4_8->SetLineStyle(2);
  h4_8->SetLineColor(kViolet);
  h4_8->SetStats(kFALSE);
  Double_t scaleh48 = 1./h4_8->Integral();
  h4_8->Scale(scaleh48);
  h4_8->Draw("sames");
 
gPad->Modified(); gPad->Update();


 TLegend *tleg4_1 = new TLegend(0.2305481,0.665852,0.5550148,0.8811685,NULL,"brNDC");
tleg4_1 ->SetTextSize(0.043);
tleg4_1->SetBorderSize(1);
tleg4_1->SetLineColor(1);
tleg4_1->SetLineStyle(1);
tleg4_1->SetLineWidth(2);
tleg4_1->SetFillColor(19);
tleg4_1->SetFillStyle(1001);
  tleg4_1->AddEntry(h4_1,"method1");
  tleg4_1->AddEntry(h4_2,"method2");
  tleg4_1->AddEntry(h4_3,"method3");
  tleg4_1->AddEntry(h4_4,"method4");
  tleg4_1->AddEntry(h4_5,"method5");
  tleg4_1->AddEntry(h4_6,"method6");
  tleg4_1->AddEntry(h4_7,"method7");
  tleg4_1->AddEntry(h4_8,"method8");

   tleg4_1->Draw();

//******************************
 c1->SaveAs("testN_mmmm4.png");
}

