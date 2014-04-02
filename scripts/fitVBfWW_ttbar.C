const float IntLUMI = 19.3;

RooRealVar *mjj_;
using namespace RooFit;


void fitVBfWW_ttbar(){

  TFile* fin
    = new TFile("mu_VBFWW_TTbarControlRegion_corr.root");

// Read Histograms 
  TH1D* data_obs = (TH1D*) fin->Get("data_obs");
  TH1D* th1ww = (TH1D*) fin->Get("th1ww");
  TH1D* th1wz = (TH1D*) fin->Get("th1wz");
  TH1D* th1zz = (TH1D*) fin->Get("th1zz");
  TH1D* th1sumdiboson = (TH1D*) fin->Get("th1sumdiboson");
  TH1D* th1w4jets = (TH1D*) fin->Get("th1w4jets");
  TH1D* th1w3jets = (TH1D*) fin->Get("th1w3jets");
  TH1D* th1w2jets = (TH1D*) fin->Get("th1w2jets");
  TH1D* th1w1jets = (TH1D*) fin->Get("th1w1jets");
  TH1D* th1sumwjets = (TH1D*) fin->Get("th1sumwjets");
  TH1D* th1zjets = (TH1D*) fin->Get("th1zjets");
  TH1D* th1Top = (TH1D*) fin->Get("th1Top");
  TH1D* th1stops = (TH1D*) fin->Get("th1stops");
  TH1D* th1stopt = (TH1D*) fin->Get("th1stopt");
  TH1D* th1stoptw = (TH1D*) fin->Get("th1stoptw");
  TH1D* th1stopps = (TH1D*) fin->Get("th1stopps");
  TH1D* th1stoppt = (TH1D*) fin->Get("th1stoppt");
  TH1D* th1stopptw = (TH1D*) fin->Get("th1stopptw");
  TH1D* th1sumstop = (TH1D*) fin->Get("th1sumstop");
  TH1D* th1madww2jets = (TH1D*) fin->Get("th1madww2jets");
  TH1D* th1mh126 = (TH1D*) fin->Get("th1mh126");

// Get Bin nos from data hist
  int NbinsX = data_obs->GetNbinsX();
  double xmin = data_obs->GetXaxis()->GetBinLowEdge(1);
  double xmax = data_obs->GetXaxis()->GetBinLowEdge(NbinsX+1);

//choose bin range to use
  data_obs->GetXaxis()->SetRangeUser(xmin, xmax);
  th1w4jets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1w3jets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1w2jets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1w1jets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1sumwjets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1ww->GetXaxis()->SetRangeUser(xmin, xmax);
  th1wz->GetXaxis()->SetRangeUser(xmin, xmax);
  th1zz->GetXaxis()->SetRangeUser(xmin, xmax);
  th1sumdiboson->GetXaxis()->SetRangeUser(xmin,xmax);
  th1zjets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1Top->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stops->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stopt->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stoptw->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stopps->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stoppt->GetXaxis()->SetRangeUser(xmin, xmax);
  th1stopptw->GetXaxis()->SetRangeUser(xmin, xmax);
  th1sumstop->GetXaxis()->SetRangeUser(xmin, xmax);

  th1madww2jets->GetXaxis()->SetRangeUser(xmin, xmax);
  th1mh126->GetXaxis()->SetRangeUser(xmin, xmax);

  mjj_ = new RooRealVar( "Mjj", "m_{jj}", xmin, xmax, "GeV");
  RooRealVar Mass = *mjj_;
// Define Roopdf for each contributing component hist
  RooDataHist* data = new RooDataHist("data","data", *mjj_, data_obs);
  RooHistPdf* pdfw4jets = makePdf(th1w4jets, "pdfw4jets");
  RooHistPdf* pdfw3jets = makePdf(th1w3jets, "pdfw3jets");
  RooHistPdf* pdfw2jets = makePdf(th1w2jets, "pdfw2jets");
  RooHistPdf* pdfw1jets = makePdf(th1w1jets, "pdfw1jets");
  RooHistPdf* pdfsumwjets = makePdf(th1sumwjets, "pdfsumwjets");

  RooHistPdf* pdfww = makePdf(th1ww, "pdfww");
  RooHistPdf* pdfwz = makePdf(th1wz, "pdfwz");
  RooHistPdf* pdfzz = makePdf(th1zz, "pdfzz");
  RooHistPdf* pdfsumdiboson = makePdf(th1sumdiboson, "pdfsumdiboson");
  RooHistPdf* pdfzjets = makePdf(th1zjets, "pdfzjets");
  RooHistPdf* pdfTop = makePdf(th1Top, "pdfTop");
  RooHistPdf* pdfstops = makePdf(th1stops, "pdfstops");
  RooHistPdf* pdfstopt = makePdf(th1stopt, "pdfstopt");
  RooHistPdf* pdfstoptw = makePdf(th1stoptw, "pdfstoptw");
  RooHistPdf* pdfstopps = makePdf(th1stopps, "pdfstopps");
  RooHistPdf* pdfstoppt = makePdf(th1stoppt, "pdfstoppt");
  RooHistPdf* pdfstopptw = makePdf(th1stopptw, "pdfstopptw");
  RooHistPdf* pdfsumstop = makePdf(th1sumstop, "pdfsumstop");

  RooHistPdf* pdfmadww2jets = makePdf(th1madww2jets, "pdfmadww2jets");
  RooHistPdf* pdfmh126 = makePdf(th1mh126, "pdfmh126");

// find Norm 
  double w4jetsNorm = th1w4jets->Integral();
  double w3jetsNorm = th1w3jets->Integral();
  double w2jetsNorm = th1w2jets->Integral();
  double w1jetsNorm = th1w1jets->Integral();
  double sumwjetsNorm = th1sumwjets->Integral();

  double wwNorm = th1ww->Integral();
  double wzNorm = th1wz->Integral();
  double zzNorm = th1zz->Integral();
  double sumdibosonNorm = th1sumdiboson->Integral();

  double zjetsNorm = th1zjets->Integral();
  double TopNorm = th1Top->Integral();
  double stoptNorm = th1stopt->Integral();
  double stopsNorm = th1stops->Integral();
  double stoptwNorm = th1stoptw->Integral();
  double stopptNorm = th1stoppt->Integral();
  double stoppsNorm = th1stopps->Integral();
  double stopptwNorm = th1stopptw->Integral();
  double sumstopNorm = th1sumstop->Integral();

  double madww2jetsNorm = th1madww2jets->Integral();
  double mh126Norm = th1mh126->Integral();
//  RooRealVar nw3jets("nw3jets","nw3jets", fkNorm, 0.0, 1000.);
//  RooRealVar nw4jets("nw4jets","nw4jets", 400.0, 0.0, 10000.);
//  RooRealVar nsumwjets("nsumwjets","nsumwjets", 400.0, 0.0, 35000.);
  RooRealVar nsumwjets("nsumwjets","nsumwjets", sumwjetsNorm);
  RooRealVar nww("nww","nww", wwNorm);
  RooRealVar nwz("nwz","nwz", wzNorm);
  RooRealVar nzz("nzz","nzz", zzNorm);
  RooRealVar nsumdiboson("nsumdiboson","nsumdiboson", sumdibosonNorm);
  RooRealVar nzjets("nzjets","nzjets", zjetsNorm);

//  RooRealVar nTop("nTop","nTop", TopNorm);
  RooRealVar nTop("nTop","nTop", TopNorm, 0.0, 55000.);
//  RooRealVar nsumstop("nsumstop","nsumstop", sumstopNorm, 0.0, 25000.);
  RooRealVar nsumstop("nsumstop","nsumstop", sumstopNorm);
//  RooArgList* components
//    = new RooArgList(*pdfsumwjets, *pdfww, *pdfwz, *pdfzz, *pdfzjets, *pdfTop, *pdfsumstop);
//  RooArgList* yields = new RooArgList(nsumwjets, nww, nwz, nzz, nzjets, nTop,nsumstop);

  RooArgList* components
    = new RooArgList(*pdfsumwjets, *pdfsumdiboson, *pdfzjets, *pdfTop, *pdfsumstop);
  RooArgList* yields = new RooArgList(nsumwjets, nsumdiboson, nzjets, nTop,nsumstop);


  RooAddPdf totalPdf("totalPdf","extended sum pdf", *components, *yields);
//  RooGaussian consNwjets("consNw3jets","", nw3jets, RooConst(fkNorm),RooConst(0.2*w3jetsNorm)) ;
//  RooGaussian consNww("consNww","", nww, RooConst(wwNorm),RooConst(0.1*wwNorm)) ;
//  RooGaussian consNwz("consNwz","", nwz, RooConst(wzNorm),RooConst(0.1*wzNorm)) ;
//  RooGaussian consNzz("consNzz","", nzz, RooConst(zzNorm),RooConst(0.1*zzNorm)) ;
  RooGaussian consNsumwjets("consNsumwjets","", nsumwjets, RooConst(sumwjetsNorm),RooConst(0.2*sumwjetsNorm));
  RooGaussian consNsumdiboson("consNsumdiboson","", nsumdiboson, RooConst(sumdibosonNorm),RooConst(0.1*sumdibosonNorm)) ;
  RooGaussian consNzjets("consNzjets","", nzjets, RooConst(zjetsNorm),RooConst(0.10*zjetsNorm)) ;
  RooGaussian consNsumstop("consNsumstop","", nsumstop, RooConst(sumstopNorm),RooConst(0.5*sumstopNorm)) ;
//  RooGaussian consNTop("consNTop","", nTop, RooConst(TopNorm),RooConst(0.063*TopNorm)) ;


  RooFitResult *fitResult
    = totalPdf.fitTo(*data, Save(true),
//ExternalConstraints(consNfk),
//ExternalConstraints(consNww),
//ExternalConstraints(consNwz),
//ExternalConstraints(consNzz),
ExternalConstraints(consNsumwjets),
ExternalConstraints(consNsumdiboson),
ExternalConstraints(consNzjets),
ExternalConstraints(consNsumstop),
//ExternalConstraints(consNTop),
RooFit::Extended(true),
//RooFit::Minos(true),
//RooFit::Hesse(false),
//PrintEvalErrors(-1),
// RooFit::Range(rangeString),
Warnings(false)
);

   fitResult->Print("v");


   std::cout << "===================== TTBar k-factor = " <<
//     nsumwjets.getVal() / sumwjetsNorm << " +- " << nsumwjets.getError() / sumwjetsNorm << std::endl;
     nTop.getVal() / TopNorm << " +- " << nTop.getError() / TopNorm << std::endl;


   // ********** Make and save Canvas for the plots ********** //
   gROOT->ProcessLine(".L ~kalanand/tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.19);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);
  RooAbsData::ErrorType errorType = RooAbsData::SumW2;



   TCanvas* c = new TCanvas("fit","",500,500);
   RooPlot* frame1 = Mass.frame();
   data->plotOn(frame1,RooFit::DataError(errorType), Name("h_data"));
   totalPdf.plotOn(frame1,ProjWData(*data), Name("h_total"));
   totalPdf.plotOn(frame1,ProjWData(*data),Components("pdfTop"),
LineColor(kRed), LineStyle(2), Name("h_ttbar"));

   totalPdf.plotOn(frame1,ProjWData(*data),Components("pdfsumdiboson,pdfzjets,pdfsumswjets,pdfsumstop"),
LineColor(kBlack), LineStyle(2), Name("h_others"));
   totalPdf.plotOn(frame1,ProjWData(*data));

   frame1->SetMinimum(0);
   frame1->SetMaximum(1.35* frame1->GetMaximum());
   frame1->Draw("e0");


   std::cout << "===================== chi2/ dof = " << frame1->chiSquare() << std::endl;

   TPaveText *plotlabel4 = new TPaveText(0.25,0.66,0.5,0.81,"NDC");
   plotlabel4->SetTextColor(kBlack);
   plotlabel4->SetFillColor(kWhite);
   plotlabel4->SetBorderSize(0);
   plotlabel4->SetTextAlign(12);
   plotlabel4->SetTextSize(0.04);
   //char temp[50];
   //sprintf(temp, "#chi^{2} / dof = %.2f", frame1->chiSquare());
   //plotlabel4->AddText(temp);
   //plotlabel4->Draw();

   cmsPrelim2();
   TLegend* legend = new TLegend(0.55,0.72,0.88,0.91);
   RooHist* datahist = frame1->getHist("h_data");
   RooCurve* totalhist = frame1->getCurve("h_total");
   RooCurve* ttbarhist = frame1->getCurve("h_ttbar");
   RooCurve* otherhist = frame1->getCurve("h_others");

   legend->AddEntry( datahist, "Data", "PE");
   legend->AddEntry( totalhist, "Fit", "L");
   legend->AddEntry( ttbarhist, "TTBar", "L");
   legend->AddEntry( otherhist, "Other processes", "L");
   legend->SetFillColor(0);
   legend->Draw();
   c->SaveAs( "WW_ttbarKfactorFit.png");
   c->SaveAs( "WW_ttbarKfactorFit.pdf");
}


RooHistPdf* makePdf(TH1D* hist, char* name) {

  RooDataHist* rdh = new RooDataHist("rdh","", *mjj_, hist);
  RooHistPdf* pdf = new RooHistPdf(name,name,*mjj_,*rdh);
  return pdf;
}




////CMS Preliminary label and lumu -- upper left corner
void cmsPrelim2()
{
   const float LUMINOSITY = IntLUMI;
   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);

   latex.SetTextAlign(31); // align right
   latex.DrawLatex(0.90,0.96,"#sqrt{s} = 8 TeV");
   if (LUMINOSITY > 0.) {
      latex.SetTextAlign(11); // align left
      latex.DrawLatex(0.21,0.85,Form("#int #font[12]{L} dt = %.1f fb^{-1}", LUMINOSITY));
   }
   latex.SetTextAlign(11); // align left
   latex.DrawLatex(0.18,0.96,"CMS preliminary");
}
