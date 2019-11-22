#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<map>
#include<set>
#include "stdlib.h"
using namespace std;

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"

void paper_plot_01_bias()
{
  gROOT->ProcessLine(".x ./lhcbStyle_edit.C");
  TString roostr = "";
  
  TGaxis *yaxis = new TGaxis();
  yaxis->SetMaxDigits(3);

  TFile *input_file = new TFile("example1_RelBias_Mu15_Ndet10.root", "read");
  TH1D *h1_rel_Poisson = (TH1D*)input_file->Get("h1_Poisson");
  TH1D *h1_rel_Neyman  = (TH1D*)input_file->Get("h1_Neyman");
  TH1D *h1_rel_Pearson = (TH1D*)input_file->Get("h1_Pearson");
  TH1D *h1_rel_CNP     = (TH1D*)input_file->Get("h1_CNP");
  TH1D *h1_rel_Gauss     = (TH1D*)input_file->Get("h1_Gauss");
  
  int color_Neyman = kRed;
  int color_Pearson = kBlue;
  int color_Poisson = kCyan;
  int color_Gauss = kOrange;  
  int color_CNP = kBlack;
 
  // color_Neyman = 1;
  // color_Pearson = 1;
  // color_Poisson = 1;
  // color_Gauss = 1;
  // color_CNP = 1;


  roostr = "canv_simple_example_bias";
  TCanvas *canv_simple_example_bias = new TCanvas(roostr, roostr, 900, 650);

  canv_simple_example_bias->SetRightMargin(0.05);
  canv_simple_example_bias->SetTopMargin(0.09);
  canv_simple_example_bias->SetBottomMargin(0.18);
  canv_simple_example_bias->SetLeftMargin(0.12);

  TH1D *h1_basic = new TH1D("h1_basic","h1_basic",200,-40,40);
  for(int ibin=1; ibin<=200; ibin++) h1_basic->SetBinContent(ibin, -200);
  h1_basic->SetLineWidth(0.1);
    
  h1_basic->Draw();
  h1_basic->SetMinimum( 0 );
  h1_basic->SetMaximum( 2.8 * h1_rel_Poisson->GetMaximum() );
  h1_basic->SetLineColor(kBlack);
  h1_basic->SetXTitle("(#hat{#mu}-#mu)/#mu [%]");
  h1_basic->GetXaxis()->CenterTitle();
  h1_basic->GetYaxis()->CenterTitle();
  h1_basic->SetYTitle("Entries");
  h1_basic->GetXaxis()->SetLabelSize(0.07);
  h1_basic->GetXaxis()->SetTitleSize(0.07);
  h1_basic->GetYaxis()->SetLabelSize(0.07);
  h1_basic->GetYaxis()->SetTitleSize(0.07);  
  h1_basic->GetXaxis()->SetTitleOffset(1.18);
  h1_basic->GetYaxis()->SetTitleOffset(0.85);

  h1_rel_Poisson->Draw("same");
  h1_rel_Poisson->SetLineColorAlpha(color_Poisson, 0);
  h1_rel_Poisson->SetFillColorAlpha(color_Poisson, 0.5);
  h1_rel_Poisson->SetFillStyle(1001);

  h1_rel_Gauss->Draw("same");
  h1_rel_Gauss->SetLineColorAlpha(color_Gauss, 0);
  h1_rel_Gauss->SetFillColorAlpha(color_Gauss, 1);
  h1_rel_Gauss->SetFillStyle(3021);

  h1_rel_Neyman->Draw("same");
  h1_rel_Neyman->SetLineColorAlpha(color_Neyman, 1);
  h1_rel_Neyman->SetLineStyle(7);

  h1_rel_Pearson->Draw("same");
  h1_rel_Pearson->SetLineColorAlpha(color_Pearson, 1);
  h1_rel_Pearson->SetLineStyle(5);

  h1_rel_CNP->Draw("same");
  h1_rel_CNP->SetLineColorAlpha(color_CNP, 1);
  h1_rel_CNP->SetLineStyle(1);

  h1_basic->Draw("same axis");

  TLegend *lg_rel = new TLegend(0.154+0.2, 0.580-0.1, 0.6247+0.28, 0.881+0.005);
  lg_rel->SetBorderSize(1);
  lg_rel->SetTextFont(42);
  lg_rel->SetHeader("                    Rel.Bias, RMS");
  lg_rel->AddEntry(h1_rel_Poisson,  TString::Format("Poisson        0, %4.2f", h1_rel_Poisson->GetRMS() ), "F" );
  lg_rel->AddEntry(h1_rel_Gauss,    TString::Format("Gauss    %5.2f, %4.2f", h1_rel_Gauss->GetMean(), h1_rel_Gauss->GetRMS() ), "F" );
  lg_rel->AddEntry(h1_rel_Neyman,   TString::Format("Neyman %5.2f, %4.2f", h1_rel_Neyman->GetMean(), h1_rel_Neyman->GetRMS() ), "l" );
  lg_rel->AddEntry(h1_rel_Pearson,  TString::Format("Pearson %5.2f, %4.2f", h1_rel_Pearson->GetMean(), h1_rel_Pearson->GetRMS() ), "l" );
  lg_rel->AddEntry(h1_rel_CNP,      TString::Format("CNP       %5.2f, %4.2f", h1_rel_CNP->GetMean(), h1_rel_CNP->GetRMS() ), "l" );
  lg_rel->Draw("same");
  lg_rel->SetTextFont(42);
  lg_rel->SetTextSize(0.065);
  
  canv_simple_example_bias->SaveAs("canv_simple_example_bias.pdf");
    
  ////////////////////////////////////////////
  /*
  roostr = "canv_rel";
  TCanvas *canv_rel = new TCanvas(roostr, roostr, 900, 650);
  canv_rel->SetRightMargin(0.05);
  canv_rel->SetTopMargin(0.09);
  canv_rel->SetBottomMargin(0.18);
  canv_rel->SetLeftMargin(0.12);

  TH1D *h1_basic = new TH1D("h1_basic","h1_basic",200,-40,40);
  h1_basic->Draw();
  h1_basic->SetMaximum( 2.8 * h1_rel_Poisson->GetMaximum() );
  h1_basic->SetLineColor(kBlack);
  h1_basic->SetXTitle("(#hat{#mu}-#mu)/#mu [%]");
  h1_basic->GetXaxis()->CenterTitle();
  h1_basic->GetYaxis()->CenterTitle();
  h1_basic->SetYTitle("Entries");
  h1_basic->GetXaxis()->SetLabelSize(0.07);
  h1_basic->GetXaxis()->SetTitleSize(0.07);
  h1_basic->GetYaxis()->SetLabelSize(0.07);
  h1_basic->GetYaxis()->SetTitleSize(0.07);  
  h1_basic->GetXaxis()->SetTitleOffset(1.18);
  h1_basic->GetYaxis()->SetTitleOffset(0.85);

  h1_rel_Poisson->Draw("same");
  h1_rel_Neyman->Draw("same");  
  h1_rel_Pearson->Draw("same");
  h1_rel_Gauss->Draw("same");
  h1_rel_CNP->Draw("same");

  h1_rel_Neyman->SetLineColor(color_Neyman);
  h1_rel_Pearson->SetLineColor(color_Pearson);
  h1_rel_Poisson->SetLineColor(color_Poisson);
  h1_rel_Gauss->SetLineColor(color_Gauss);
  h1_rel_CNP->SetLineColor(color_CNP);


  h1_rel_Neyman->SetLineStyle(7);
  h1_rel_Pearson->SetLineStyle(5);
  
  // h1_rel_Poisson->SetLineStyle(3);
  // h1_rel_Neyman->SetLineStyle(3);
  // h1_rel_Pearson->SetLineStyle(5);
  // h1_rel_CNP->SetLineStyle(1);  
  // h1_rel_Gauss->SetLineStyle(1); 
  

 
  // h1_rel_Poisson->SetFillColor(color_Poisson);
  // h1_rel_Poisson->SetFillStyle(1001);

  // h1_rel_Poisson->SetLineStyle(1);
  // h1_rel_Poisson->SetLineColor(color_Poisson);

  // h1_rel_Gauss->SetFillStyle(3021);
  // h1_rel_Gauss->SetFillColor(color_Gauss);
  
  
  
  
  h1_basic->Draw("same axis");

  TLegend *lg_rel = new TLegend(0.154+0.2, 0.580-0.1, 0.6247+0.28, 0.881+0.005);
  lg_rel->SetBorderSize(1);
  lg_rel->SetTextFont(42);
  lg_rel->SetHeader("                    Rel.Bias, RMS");
  lg_rel->AddEntry(h1_rel_Poisson,  TString::Format("Poisson        0, %4.2f", h1_rel_Poisson->GetRMS() ), "F" );
  lg_rel->AddEntry(h1_rel_Gauss,    TString::Format("Gauss    %5.2f, %4.2f", h1_rel_Gauss->GetMean(), h1_rel_Gauss->GetRMS() ), "F" );
  lg_rel->AddEntry(h1_rel_Neyman,   TString::Format("Neyman %5.2f, %4.2f", h1_rel_Neyman->GetMean(), h1_rel_Neyman->GetRMS() ), "l" );
  lg_rel->AddEntry(h1_rel_Pearson,  TString::Format("Pearson %5.2f, %4.2f", h1_rel_Pearson->GetMean(), h1_rel_Pearson->GetRMS() ), "l" );
  lg_rel->AddEntry(h1_rel_CNP,      TString::Format("CNP       %5.2f, %4.2f", h1_rel_CNP->GetMean(), h1_rel_CNP->GetRMS() ), "l" );
  lg_rel->Draw("same");
  lg_rel->SetTextFont(42);
  lg_rel->SetTextSize(0.065);
  
  canv_rel->SaveAs("canv_simple_example_bias.png");
  canv_rel->SaveAs("canv_simple_example_bias.pdf");
  */
}
