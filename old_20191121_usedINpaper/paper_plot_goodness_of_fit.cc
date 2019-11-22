#include<iostream>
#include<fstream>
#include<cmath>
#include "stdlib.h" 
using namespace std;

#include<map>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
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

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "Math/Polynomial.h"

#include "format.icc"

void paper_plot_goodness_of_fit()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");

  int color_Poisson = kCyan;
  int color_Neyman = kRed;
  int color_Pearson = kBlue;
  int color_Gauss = kOrange;
  int color_CNP = kBlack;
  
  // color_Neyman = 1;
  // color_Pearson = 1;
  // color_Poisson = 1;
  // color_Gauss = 1;
  // color_CNP = 1;
  

  TString roostr = "";





  
  TFile *roofile = new TFile("./file_goodness_of_fit_Ndet10_v2.root", "read");
  TTree *tree_bias = (TTree*)roofile->Get("tree_bias");
  
  // Declaration of leaf types
  Int_t           ndf_val;
  Double_t        mutrue;
  Double_t        chi2hat_Poisson;
  Double_t        chi2hat_Neyman;
  Double_t        chi2hat_Pearson;
  Double_t        chi2hat_Gauss;
  Double_t        chi2hat_CNP;

  // List of branches
  TBranch        *b_ndf_val;   //!
  TBranch        *b_mutrue;   //!
  TBranch        *b_chi2hat_Poisson;   //!
  TBranch        *b_chi2hat_Neyman;   //!
  TBranch        *b_chi2hat_Pearson;   //!
  TBranch        *b_chi2hat_Gauss;   //!
  TBranch        *b_chi2hat_CNP;   //!

  // Set branch addresses and branch pointers
  tree_bias->SetBranchAddress("ndf_val", &ndf_val, &b_ndf_val);
  tree_bias->SetBranchAddress("mutrue", &mutrue, &b_mutrue);
  tree_bias->SetBranchAddress("chi2hat_Poisson", &chi2hat_Poisson, &b_chi2hat_Poisson);
  tree_bias->SetBranchAddress("chi2hat_Neyman", &chi2hat_Neyman, &b_chi2hat_Neyman);
  tree_bias->SetBranchAddress("chi2hat_Pearson", &chi2hat_Pearson, &b_chi2hat_Pearson);
  tree_bias->SetBranchAddress("chi2hat_Gauss", &chi2hat_Gauss, &b_chi2hat_Gauss);
  tree_bias->SetBranchAddress("chi2hat_CNP", &chi2hat_CNP, &b_chi2hat_CNP);

  int entries = tree_bias->GetEntries();

  //////////////////////////

  TGraph *gh15_rel_Poisson = new TGraph(); roostr = "gh15_rel_Poisson"; gh15_rel_Poisson->SetName(roostr);
  TGraph *gh15_rel_Neyman = new TGraph(); roostr = "gh15_rel_Neyman"; gh15_rel_Neyman->SetName(roostr);
  TGraph *gh15_rel_Pearson = new TGraph(); roostr = "gh15_rel_Pearson"; gh15_rel_Pearson->SetName(roostr);
  TGraph *gh15_rel_Gauss = new TGraph(); roostr = "gh15_rel_Gauss"; gh15_rel_Gauss->SetName(roostr);
  TGraph *gh15_rel_CNP = new TGraph(); roostr = "gh15_rel_CNP"; gh15_rel_CNP->SetName(roostr);

  TGraph *gh15_abs_Poisson = new TGraph(); roostr = "gh15_abs_Poisson"; gh15_abs_Poisson->SetName(roostr);
  TGraph *gh15_abs_Neyman = new TGraph(); roostr = "gh15_abs_Neyman"; gh15_abs_Neyman->SetName(roostr);
  TGraph *gh15_abs_Pearson = new TGraph(); roostr = "gh15_abs_Pearson"; gh15_abs_Pearson->SetName(roostr);
  TGraph *gh15_abs_Gauss = new TGraph(); roostr = "gh15_abs_Gauss"; gh15_abs_Gauss->SetName(roostr);
  TGraph *gh15_abs_CNP = new TGraph(); roostr = "gh15_abs_CNP"; gh15_abs_CNP->SetName(roostr);

  int line15 = 0;

  //////////////////////////
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_bias->GetEntry(ientry);

      line15++;      
      gh15_abs_Poisson->SetPoint(line15-1, mutrue, chi2hat_Poisson/ndf_val - 1);
      gh15_abs_Neyman->SetPoint(line15-1, mutrue, chi2hat_Neyman/ndf_val - 1);
      gh15_abs_Pearson->SetPoint(line15-1, mutrue, chi2hat_Pearson/ndf_val - 1);
      gh15_abs_Gauss->SetPoint(line15-1, mutrue, chi2hat_Gauss/ndf_val - 1);
      gh15_abs_CNP->SetPoint(line15-1, mutrue, chi2hat_CNP/ndf_val - 1);

      gh15_rel_Poisson->SetPoint(line15-1, mutrue, chi2hat_Poisson/ndf_val );
      gh15_rel_Neyman->SetPoint(line15-1, mutrue, chi2hat_Neyman/ndf_val );
      gh15_rel_Pearson->SetPoint(line15-1, mutrue, chi2hat_Pearson/ndf_val );
      gh15_rel_Gauss->SetPoint(line15-1, mutrue, chi2hat_Gauss/ndf_val );
      gh15_rel_CNP->SetPoint(line15-1, mutrue, chi2hat_CNP/ndf_val );
    
  }

  ///////////////////////////////////////////////////////////////////////

  roostr = "canv_abs_AA";
  TCanvas *canv_abs_AA = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_abs_AA, 0.18, 0.09, 0.2, 0.05);
  canv_abs_AA->SetLogx();
  
  roostr = "h1_basic_AA";
  TH1D *h1_basic_AA = new TH1D(roostr, roostr, 100, 0.5, 2000);
  for(int ibin=1; ibin<=100; ibin++) h1_basic_AA->SetBinContent(ibin, -200);
  
  h1_basic_AA->GetYaxis()->SetRangeUser(-0.1, 0.6);
  h1_basic_AA->Draw();
  h1_basic_AA->SetLineWidth(0.1);
  roo_title(h1_basic_AA, "#mu", "(<#chi^{2}> - ndf) / ndf");
  roo_center_title(h1_basic_AA);
  roo_offset(h1_basic_AA, 1.2, 1.1);
  h1_basic_AA->GetXaxis()->SetLabelSize(0.07);
  h1_basic_AA->GetXaxis()->SetTitleSize(0.07);
  h1_basic_AA->GetYaxis()->SetLabelSize(0.07);
  h1_basic_AA->GetYaxis()->SetTitleSize(0.07);

  double markerSize_AA = 1.4;
  
  //////
  gh15_abs_Poisson->Draw("same pl");
  gh15_abs_Poisson->SetLineColor(color_Poisson);
  gh15_abs_Poisson->SetLineStyle(1);
  gh15_abs_Poisson->SetMarkerStyle(26);
  gh15_abs_Poisson->SetMarkerSize(markerSize_AA);
  gh15_abs_Poisson->SetMarkerColor(color_Poisson);
  
  gh15_abs_Gauss->Draw("same pl");
  gh15_abs_Gauss->SetLineColor(color_Gauss);
  gh15_abs_Gauss->SetLineStyle(1);
  gh15_abs_Gauss->SetMarkerStyle(30);
  gh15_abs_Gauss->SetMarkerSize(markerSize_AA+0.2);
  gh15_abs_Gauss->SetMarkerColor(color_Gauss);
  
  gh15_abs_Pearson->Draw("same pl");
  gh15_abs_Pearson->SetLineColor(color_Pearson);
  gh15_abs_Pearson->SetLineStyle(1);
  gh15_abs_Pearson->SetMarkerStyle(24);
  gh15_abs_Pearson->SetMarkerSize(markerSize_AA);
  gh15_abs_Pearson->SetMarkerColor(color_Pearson);
  
  gh15_abs_Neyman->Draw("same pl");
  gh15_abs_Neyman->SetLineColor(color_Neyman);
  gh15_abs_Neyman->SetLineStyle(1);
  gh15_abs_Neyman->SetMarkerStyle(25);
  gh15_abs_Neyman->SetMarkerSize(markerSize_AA);
  gh15_abs_Neyman->SetMarkerColor(color_Neyman);
  
  gh15_abs_CNP->Draw("same pl");
  gh15_abs_CNP->SetLineColor(color_CNP);
  gh15_abs_CNP->SetLineStyle(1);
  gh15_abs_CNP->SetMarkerStyle(27);
  gh15_abs_CNP->SetMarkerSize(markerSize_AA+0.2);
  gh15_abs_CNP->SetMarkerColor(color_CNP);

  h1_basic_AA->Draw("same axis");

  TLegend *lg_basic_AA = new TLegend(0.5991091+0.03,0.4807692,0.8964365,0.8766026);
  lg_basic_AA->SetBorderSize(1);
  lg_basic_AA->SetTextFont(42);
  // lg_basic_AA->AddEntry("", "ndf = 10", "" );
  lg_basic_AA->AddEntry(gh15_abs_Poisson, "Poisson", "pl" );
  lg_basic_AA->AddEntry(gh15_abs_Gauss, "Gauss", "pl" );
  lg_basic_AA->AddEntry(gh15_abs_Neyman, "Neyman", "pl" );
  lg_basic_AA->AddEntry(gh15_abs_Pearson, "Pearson", "pl" );
  lg_basic_AA->AddEntry(gh15_abs_CNP, "CNP", "pl" );
  lg_basic_AA->Draw("same");
  lg_basic_AA->SetTextFont(42);
  lg_basic_AA->SetTextSize(0.065);

  canv_abs_AA->SaveAs("canv_goodness_of_fit_scan_Ndet10_rel.pdf");


  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
  {
    roostr = "../data_goodness_of_fit_Ndet010_50M/out_sum_chi2_mu0015.root";
    TFile *file_ = new TFile(roostr, "read");
    TTree *tree_chi2 = (TTree*)file_->Get("tree_chi2");
    
    // Declaration of leaf types
    Double_t        Nmu;
    Int_t           Ndet;
    Double_t        mu_bestFit_Poisson;
    Double_t        mu_bestFit_Neyman;
    Double_t        mu_bestFit_Pearson;
    Double_t        mu_bestFit_CNP;
    Double_t        mu_bestFit_Gauss;
    Double_t        chi2true_Poisson;
    Double_t        chi2true_Neyman;
    Double_t        chi2true_Pearson;
    Double_t        chi2true_CNP;
    Double_t        chi2true_Gauss;
    Double_t        chi2bestFit_Poisson;
    Double_t        chi2bestFit_Neyman;
    Double_t        chi2bestFit_Pearson;
    Double_t        chi2bestFit_CNP;
    Double_t        chi2bestFit_Gauss;

    // List of branches
    TBranch        *b_Nmu;   //!
    TBranch        *b_Ndet;   //!
    TBranch        *b_mu_bestFit_Poisson;   //!
    TBranch        *b_mu_bestFit_Neyman;   //!
    TBranch        *b_mu_bestFit_Pearson;   //!
    TBranch        *b_mu_bestFit_CNP;   //!
    TBranch        *b_mu_bestFit_Gauss;   //!
    TBranch        *b_chi2true_Poisson;   //!
    TBranch        *b_chi2true_Neyman;   //!
    TBranch        *b_chi2true_Pearson;   //!
    TBranch        *b_chi2true_CNP;   //!
    TBranch        *b_chi2true_Gauss;   //!
    TBranch        *b_chi2bestFit_Poisson;   //!
    TBranch        *b_chi2bestFit_Neyman;   //!
    TBranch        *b_chi2bestFit_Pearson;   //!
    TBranch        *b_chi2bestFit_CNP;   //!
    TBranch        *b_chi2bestFit_Gauss;   //!
 
    // Set branch addresses and branch pointers
    tree_chi2->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
    tree_chi2->SetBranchAddress("Ndet", &Ndet, &b_Ndet);
    tree_chi2->SetBranchAddress("mu_bestFit_Poisson", &mu_bestFit_Poisson, &b_mu_bestFit_Poisson);
    tree_chi2->SetBranchAddress("mu_bestFit_Neyman", &mu_bestFit_Neyman, &b_mu_bestFit_Neyman);
    tree_chi2->SetBranchAddress("mu_bestFit_Pearson", &mu_bestFit_Pearson, &b_mu_bestFit_Pearson);
    tree_chi2->SetBranchAddress("mu_bestFit_CNP", &mu_bestFit_CNP, &b_mu_bestFit_CNP);
    tree_chi2->SetBranchAddress("mu_bestFit_Gauss", &mu_bestFit_Gauss, &b_mu_bestFit_Gauss);
    tree_chi2->SetBranchAddress("chi2true_Poisson", &chi2true_Poisson, &b_chi2true_Poisson);
    tree_chi2->SetBranchAddress("chi2true_Neyman", &chi2true_Neyman, &b_chi2true_Neyman);
    tree_chi2->SetBranchAddress("chi2true_Pearson", &chi2true_Pearson, &b_chi2true_Pearson);
    tree_chi2->SetBranchAddress("chi2true_CNP", &chi2true_CNP, &b_chi2true_CNP);
    tree_chi2->SetBranchAddress("chi2true_Gauss", &chi2true_Gauss, &b_chi2true_Gauss);
    tree_chi2->SetBranchAddress("chi2bestFit_Poisson", &chi2bestFit_Poisson, &b_chi2bestFit_Poisson);
    tree_chi2->SetBranchAddress("chi2bestFit_Neyman", &chi2bestFit_Neyman, &b_chi2bestFit_Neyman);
    tree_chi2->SetBranchAddress("chi2bestFit_Pearson", &chi2bestFit_Pearson, &b_chi2bestFit_Pearson);
    tree_chi2->SetBranchAddress("chi2bestFit_CNP", &chi2bestFit_CNP, &b_chi2bestFit_CNP);
    tree_chi2->SetBranchAddress("chi2bestFit_Gauss", &chi2bestFit_Gauss, &b_chi2bestFit_Gauss);

    int entries = tree_chi2->GetEntries();

    cout<<endl<<" ---> Entries "<<entries<<endl<<endl;

    //entries = 1000000;
    
    // int nbins_x = 58;
    // double xlow = 0;
    // double xhgh = 50;

    int nbins_x = 58;
    double xlow = 0;
    double xhgh = 50;
    
    roostr = "h1_chi2_Poisson"; TH1D *h1_chi2_Poisson = new TH1D(roostr, roostr, nbins_x, xlow, xhgh);
    roostr = "h1_chi2_Gauss"; TH1D *h1_chi2_Gauss = new TH1D(roostr, roostr, nbins_x, xlow, xhgh);
    roostr = "h1_chi2_Pearson"; TH1D *h1_chi2_Pearson = new TH1D(roostr, roostr, nbins_x, xlow, xhgh);
    roostr = "h1_chi2_Neyman"; TH1D *h1_chi2_Neyman = new TH1D(roostr, roostr, nbins_x, xlow, xhgh);
    roostr = "h1_chi2_CNP"; TH1D *h1_chi2_CNP = new TH1D(roostr, roostr, nbins_x, xlow, xhgh);

    double mean_Poisson = 0;
    double mean_Gauss = 0;
    double mean_Pearson = 0;
    double mean_Neyman = 0;
    double mean_CNP = 0;
    
    /////////////////////////////////
    
    int nbins_xD = 32;
    double xDlow = 0;
    double xDhgh = 16;

    
    // int nbins_xD = 80;
    // double xDlow = 0;
    // double xDhgh = 16;
    

    roostr = "h1_Dchi2_Poisson"; TH1D *h1_Dchi2_Poisson = new TH1D(roostr, roostr, nbins_xD, xDlow, xDhgh);
    roostr = "h1_Dchi2_Gauss"; TH1D *h1_Dchi2_Gauss = new TH1D(roostr, roostr, nbins_xD, xDlow, xDhgh);
    roostr = "h1_Dchi2_Pearson"; TH1D *h1_Dchi2_Pearson = new TH1D(roostr, roostr, nbins_xD, xDlow, xDhgh);
    roostr = "h1_Dchi2_Neyman"; TH1D *h1_Dchi2_Neyman = new TH1D(roostr, roostr, nbins_xD, xDlow, xDhgh);
    roostr = "h1_Dchi2_CNP"; TH1D *h1_Dchi2_CNP = new TH1D(roostr, roostr, nbins_xD, xDlow, xDhgh);

    ////////////////////////////////
    
    for(int ientry=0; ientry<entries; ientry++) {
      tree_chi2->GetEntry(ientry);
      if(ientry%500000==0) cout<<TString::Format(" ---> procssing %8.4f", ientry*100./entries )<<endl;

      h1_chi2_Poisson->Fill( chi2true_Poisson );
      h1_chi2_Gauss->Fill( chi2true_Gauss );
      h1_chi2_Pearson->Fill( chi2true_Pearson );
      h1_chi2_Neyman->Fill( chi2true_Neyman );
      h1_chi2_CNP->Fill( chi2true_CNP );

      mean_Poisson += chi2true_Poisson/entries;
      mean_Gauss += chi2true_Gauss/entries;
      mean_Pearson += chi2true_Pearson/entries;
      mean_Neyman += chi2true_Neyman/entries;
      mean_CNP += chi2true_CNP/entries;

      h1_Dchi2_Poisson->Fill( chi2true_Poisson - chi2bestFit_Poisson );
      h1_Dchi2_Gauss->Fill( chi2true_Gauss - chi2bestFit_Gauss );
      h1_Dchi2_Pearson->Fill( chi2true_Pearson - chi2bestFit_Pearson );
      h1_Dchi2_Neyman->Fill( chi2true_Neyman - chi2bestFit_Neyman );
      h1_Dchi2_CNP->Fill( chi2true_CNP - chi2bestFit_CNP );

    }

    h1_chi2_Poisson->Scale(1./entries);
    h1_chi2_Gauss->Scale(1./entries);
    h1_chi2_Pearson->Scale(1./entries);
    h1_chi2_Neyman->Scale(1./entries);
    h1_chi2_CNP->Scale(1./entries);
    
    h1_Dchi2_Poisson->Scale(1./entries);
    h1_Dchi2_Gauss->Scale(1./entries);
    h1_Dchi2_Pearson->Scale(1./entries);
    h1_Dchi2_Neyman->Scale(1./entries);
    h1_Dchi2_CNP->Scale(1./entries);

    /////////////////////////////////////////////////////////////////////////////
    
    roostr = "canv_chi2";
    TCanvas *canv_chi2 = new TCanvas(roostr, roostr, 900, 650);
    roo_canv_margin(canv_chi2, 0.18, 0.09, 0.2, 0.05);
 
    h1_chi2_Poisson->GetYaxis()->SetRangeUser(0, 0.1);
    roo_title(h1_chi2_Poisson, "#chi^{2}", "Probability");
    roo_center_title(h1_chi2_Poisson);
    roo_offset(h1_chi2_Poisson, 1.2, 1.1);
    h1_chi2_Poisson->GetXaxis()->SetLabelSize(0.07);
    h1_chi2_Poisson->GetXaxis()->SetTitleSize(0.07);
    h1_chi2_Poisson->GetYaxis()->SetLabelSize(0.07);
    h1_chi2_Poisson->GetYaxis()->SetTitleSize(0.07);
    h1_chi2_Poisson->GetYaxis()->SetNdivisions(509);
    h1_chi2_Poisson->Draw();
    h1_chi2_Poisson->SetLineColorAlpha(color_Poisson, 0);
    h1_chi2_Poisson->SetLineStyle(1);
    h1_chi2_Poisson->SetFillColorAlpha(color_Poisson, 0.5);
    h1_chi2_Poisson->SetFillStyle(1001);

    h1_chi2_Gauss->Draw("same");
    h1_chi2_Gauss->SetLineColorAlpha(color_Gauss, 0.5);
    h1_chi2_Gauss->SetLineStyle(1);
    h1_chi2_Gauss->SetFillColor(color_Gauss);
    h1_chi2_Gauss->SetFillStyle(3021);

    h1_chi2_Neyman->Draw("same");
    h1_chi2_Neyman->SetLineColor(color_Neyman);
    h1_chi2_Neyman->SetLineStyle(7);

    h1_chi2_Pearson->Draw("same");
    h1_chi2_Pearson->SetLineColor(color_Pearson);
    h1_chi2_Pearson->SetLineStyle(5);
    
    h1_chi2_CNP->Draw("same");
    h1_chi2_CNP->SetLineColor(color_CNP);
    h1_chi2_CNP->SetLineStyle(1);
    
    double width_chi2_10 = (xhgh-xlow)/nbins_x;
    TF1 *f_chi2_10 = new TF1("f_chi2_10",
			     TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", 10, width_chi2_10),
			     xlow, xhgh);
    f_chi2_10->Draw("same");
    f_chi2_10->SetLineColorAlpha(kRed, 0.5);


    h1_chi2_Poisson->Draw("same axis");
    
    TLegend *lg_basic_chi2 = new TLegend(0.518931,0.4791667-0.08,0.9187082,0.875);
    lg_basic_chi2->SetBorderSize(1);
    lg_basic_chi2->SetTextFont(42);
    lg_basic_chi2->AddEntry(f_chi2_10, "#chi^{2}(10)", "l" );
    lg_basic_chi2->AddEntry(h1_chi2_Poisson, TString::Format("Poisson, %5.2f", mean_Poisson), "F" );
    lg_basic_chi2->AddEntry(h1_chi2_Gauss,   TString::Format("Gauss,   %5.2f", mean_Gauss), "F" );
    lg_basic_chi2->AddEntry(h1_chi2_Neyman,  TString::Format("Neyman, %5.2f", mean_Neyman), "l" );
    lg_basic_chi2->AddEntry(h1_chi2_Pearson, TString::Format("Pearson, %5.2f", mean_Pearson), "l" );
    lg_basic_chi2->AddEntry(h1_chi2_CNP,     TString::Format("CNP,       %5.2f", mean_CNP), "l" );

    lg_basic_chi2->Draw("same");
    lg_basic_chi2->SetTextFont(42);
    lg_basic_chi2->SetTextSize(0.065);

    canv_chi2->SaveAs("canv_goodness_of_fit_mu15_Ndet10.pdf");
    
    //////////////////////////////////////////////////////////////////////////////////////////////


    
    /////////////////////////////////////////////////////////////////////////////
    
    roostr = "canv_Dchi2";
    TCanvas *canv_Dchi2 = new TCanvas(roostr, roostr, 900, 650);
    roo_canv_margin(canv_Dchi2, 0.18, 0.09, 0.2, 0.05);
    canv_Dchi2->SetLogy();
      
    h1_Dchi2_Poisson->GetYaxis()->SetRangeUser(1.2e-5, 9);
    roo_title(h1_Dchi2_Poisson, "#Delta#chi^{2}", "Probability");
    roo_center_title(h1_Dchi2_Poisson);
    roo_offset(h1_Dchi2_Poisson, 1.2, 1.1);
    h1_Dchi2_Poisson->GetXaxis()->SetLabelSize(0.07);
    h1_Dchi2_Poisson->GetXaxis()->SetTitleSize(0.07);
    h1_Dchi2_Poisson->GetYaxis()->SetLabelSize(0.07);
    h1_Dchi2_Poisson->GetYaxis()->SetTitleSize(0.07);
    h1_Dchi2_Poisson->GetYaxis()->SetNdivisions(509);
    h1_Dchi2_Poisson->Draw();
    h1_Dchi2_Poisson->SetLineColor(color_Poisson);

    h1_Dchi2_Poisson->SetLineColorAlpha(color_Poisson, 0);
    h1_Dchi2_Poisson->SetLineStyle(1);
    h1_Dchi2_Poisson->SetFillColorAlpha(color_Poisson, 0.5);
    h1_Dchi2_Poisson->SetFillStyle(1001);

    h1_Dchi2_Gauss->Draw("same");
    h1_Dchi2_Gauss->SetLineColorAlpha(color_Gauss, 0.5);
    h1_Dchi2_Gauss->SetLineStyle(1);
    h1_Dchi2_Gauss->SetFillColor(color_Gauss);
    h1_Dchi2_Gauss->SetFillStyle(3021);

    h1_Dchi2_Neyman->Draw("same");
    h1_Dchi2_Neyman->SetLineColor(color_Neyman);
    h1_Dchi2_Neyman->SetLineStyle(7);

    h1_Dchi2_Pearson->Draw("same");
    h1_Dchi2_Pearson->SetLineColor(color_Pearson);
    h1_Dchi2_Pearson->SetLineStyle(5);
    
    h1_Dchi2_CNP->Draw("same");
    h1_Dchi2_CNP->SetLineColor(color_CNP);
    h1_Dchi2_CNP->SetLineStyle(1);
    

    
    double width_Dchi2_10 = (xDhgh-xDlow)/nbins_xD;
    TF1 *f_Dchi2_10 = new TF1("f_Dchi2_10",
			     TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", 1, width_Dchi2_10),
			     xDlow, xDhgh);
    f_Dchi2_10->Draw("same");
    f_Dchi2_10->SetLineColorAlpha(kRed, 0.5);

    h1_Dchi2_Poisson->Draw("same axis");

    TLegend *lg_basic_Dchi2 = new TLegend(0.4231626,0.6394231,0.9187082,0.875);
    lg_basic_Dchi2->SetBorderSize(1);
    lg_basic_Dchi2->SetTextFont(42);
    lg_basic_Dchi2->SetNColumns(2);
    lg_basic_Dchi2->AddEntry(f_Dchi2_10, "#chi^{2}(1)", "l" );
    lg_basic_Dchi2->AddEntry(h1_Dchi2_Neyman, TString::Format("Neyman"), "l" ); 
    lg_basic_Dchi2->AddEntry(h1_Dchi2_Poisson, TString::Format("Poisson"), "F" );
    lg_basic_Dchi2->AddEntry(h1_Dchi2_Pearson, TString::Format("Pearson"), "l" );    
    lg_basic_Dchi2->AddEntry(h1_Dchi2_Gauss, TString::Format("Gauss"), "F" );
    lg_basic_Dchi2->AddEntry(h1_Dchi2_CNP, TString::Format("CNP"), "l" );

    lg_basic_Dchi2->Draw("same");
    lg_basic_Dchi2->SetTextFont(42);
    lg_basic_Dchi2->SetTextSize(0.065);

    canv_Dchi2->SaveAs("canv_DeltaChi2_mu15_Ndet10.pdf");


    // TFile *roofile_out = new TFile("file_paper_plot_goodness_of_fit.root", "recreate");
    
    // h1_Dchi2_Poisson->Write();
    // h1_Dchi2_Gauss->Write();
    // h1_Dchi2_Neyman->Write();
    // h1_Dchi2_Pearson->Write();
    // h1_Dchi2_CNP->Write();

    // h1_chi2_Poisson->Write();
    // h1_chi2_Gauss->Write();
    // h1_chi2_Neyman->Write();
    // h1_chi2_Pearson->Write();
    // h1_chi2_CNP->Write();

    // roofile_out->Close();
  }

  
}
