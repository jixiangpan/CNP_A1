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

void plot_check_deltaChi2_dist()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";
  
  int color_Poisson = kBlack;
  int color_Neyman  = kCyan;
  int color_Pearson = kBlue;
  int color_CNP     = kGreen;
  int color_Gauss   = kMagenta;
  
  /////////////////////////////////////////////////////////////////////////////
  
  roostr = "./data_01_check_deltaChi2_dist/out_chi2_mu0015.000_det100_run001000_event010000.root";
  TFile *file_input = new TFile(roostr, "read");
  TTree *tree_chi2 = (TTree*)file_input->Get("tree_chi2");

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

  roostr = "h1_mu_bestFit_Pearson";
  TH1D *h1_mu_bestFit_Pearson = new TH1D(roostr, roostr, 59, 12, 18);

  roostr = "h1_deltaChi2_Pearson";
  TH1D *h1_deltaChi2_Pearson = new TH1D(roostr, roostr, 100, 0, 10);
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_chi2->GetEntry( ientry );
    h1_mu_bestFit_Pearson->Fill( mu_bestFit_Pearson );

    h1_deltaChi2_Pearson->Fill( chi2true_Pearson - chi2bestFit_Pearson );
  }
  
  h1_deltaChi2_Pearson->Scale(1./entries);
    
  ////////////////////////////////////////
  
  roostr = "canv_h1_mu_bestFit_Pearson";
  TCanvas *canv_h1_mu_bestFit_Pearson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h1_mu_bestFit_Pearson, 0.18, 0.1, 0.18, 0.06);  
  h1_mu_bestFit_Pearson->Draw();  
  h1_mu_bestFit_Pearson->SetLineColor(kBlack);
  roo_title(h1_mu_bestFit_Pearson, "#hat{#mu}", "Entries");
  roo_center_title(h1_mu_bestFit_Pearson);
  roo_offset(h1_mu_bestFit_Pearson, 1.1, 1.3);

  ////////////////////////////////////////
  
  double val_mean = h1_mu_bestFit_Pearson->GetMean();
  double val_sigma = h1_mu_bestFit_Pearson->GetRMS();

  roostr = "h1_testChi2_Pearson";
  TH1D *h1_testChi2_Pearson = new TH1D(roostr, roostr, 100, 0, 10);
  int cycle_Pearson = 100000;
  
  roostr = "h1_testChi2_Pearson2";
  TH1D *h1_testChi2_Pearson2 = new TH1D(roostr, roostr, 100, 0, 10);
  
  for(int ientry=1; ientry<=cycle_Pearson; ientry++) {
    double val_x = h1_mu_bestFit_Pearson->GetRandom();
    double normal_x = (val_x-val_mean)/val_sigma;
    double val_testChi2 = normal_x * normal_x;
    h1_testChi2_Pearson->Fill( val_testChi2 );

    double user_val = pow( (val_x-15)/val_sigma, 2 );
    h1_testChi2_Pearson2->Fill( user_val );
  }
  h1_testChi2_Pearson->Scale(1./cycle_Pearson);
  h1_testChi2_Pearson2->Scale(1./cycle_Pearson);




  


  
  roostr = "canv_h1_testChi2_Pearson";
  TCanvas *canv_h1_testChi2_Pearson = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_testChi2_Pearson->SetLogy();
  roo_canv_margin(canv_h1_testChi2_Pearson, 0.18, 0.1, 0.18, 0.06);
  // h1_testChi2_Pearson->Draw();  
  // h1_testChi2_Pearson->SetLineColor(kBlue);
  // roo_title(h1_testChi2_Pearson, "#chi^{2}_{#mu} - #chi^{2}_{#hat{#mu}}", "Probability");
  // roo_center_title(h1_testChi2_Pearson);
  // roo_offset(h1_testChi2_Pearson, 1.1, 1.3);
  
  h1_deltaChi2_Pearson->Draw();
  h1_deltaChi2_Pearson->SetLineColor(kBlack);
  roo_title(h1_deltaChi2_Pearson, "#chi^{2}_{#mu} - #chi^{2}_{#hat{#mu}}", "Probability");
  roo_center_title(h1_deltaChi2_Pearson);
  roo_offset(h1_deltaChi2_Pearson, 1.1, 1.3);
  
  h1_testChi2_Pearson2->Draw("same");
  h1_testChi2_Pearson2->SetLineColor(kOrange);
  
  double low_chi2 = 0;
  double hgh_chi2 = 10;
  int nbins_chi2 = 100;
  double width_chi2_150 = (hgh_chi2-low_chi2)/nbins_chi2;
  TF1 *f_chi2_150 = new TF1("f_chi2_150",
			    TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", 1, width_chi2_150),
			    low_chi2, hgh_chi2);
  f_chi2_150->Draw("same");
  f_chi2_150->SetLineColor(kRed);
  
  h1_testChi2_Pearson->Draw("same axis");


  canv_h1_mu_bestFit_Pearson->SaveAs("canv_h1_mu_bestFit_Pearson.png");
  canv_h1_testChi2_Pearson->SaveAs("canv_h1_testChi2_Pearson.png");
}
