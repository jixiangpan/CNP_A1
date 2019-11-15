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


double func_deltaChi2_Pearson(double *x, double *par)
{
  double xcur = x[0];
  double mu_true = par[0];
  double n_dets = par[1];

  double result = n_dets * (mu_true-xcur) + n_dets*xcur*xcur*( 1./mu_true - 1./xcur );
  return result;
}



void plot_check_deltaChi2_vs_muhat()
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

  roostr = "h1_deltaChi2_Pearson";
  TH1D *h1_deltaChi2_Pearson = new TH1D(roostr, roostr, 100, 0, 10);
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_chi2->GetEntry( ientry );
    h1_deltaChi2_Pearson->Fill( chi2true_Pearson - chi2bestFit_Pearson );
  }  
  h1_deltaChi2_Pearson->Scale(1./entries);
    
  ////////////////////////////////////////
  
  roostr = "canv_h1_deltaChi2_Pearson";
  TCanvas *canv_h1_deltaChi2_Pearson = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_deltaChi2_Pearson->SetLogy();
  roo_canv_margin(canv_h1_deltaChi2_Pearson, 0.18, 0.1, 0.18, 0.06);
  h1_deltaChi2_Pearson->Draw();  
  h1_deltaChi2_Pearson->SetLineColor(kBlack);
  roo_title(h1_deltaChi2_Pearson, "#chi^{2}_{#mu} - #chi^{2}_{#hat{#mu}}", "Probability");
  roo_center_title(h1_deltaChi2_Pearson);
  roo_offset(h1_deltaChi2_Pearson, 1.1, 1.3);
  

  roostr = "roofunc_deltaChi2_Pearson";
  double f_xmin = 5;
  double f_xmax = 25;
  double f_nbins = 2000000;
  TF1 *roofunc_deltaChi2_Pearson = new TF1(roostr, func_deltaChi2_Pearson, f_xmin, f_xmax, 2);
  roofunc_deltaChi2_Pearson->SetParameters(15, 100);

  roostr = "h1_testChi2_Pearson";
  TH1D *h1_testChi2_Pearson = new TH1D(roostr, roostr, 100, 0, 10);
  
  for(int ibin=1; ibin<=f_nbins; ibin++) {
    double val_x = f_xmin + (ibin-1)* (f_xmax-f_xmin)/f_nbins;
    double val_chi2 = roofunc_deltaChi2_Pearson->Eval(val_x);
    h1_testChi2_Pearson->Fill( val_chi2 );
  }
  h1_testChi2_Pearson->Scale(1./f_nbins);

  
  h1_testChi2_Pearson->Draw("same");
  h1_testChi2_Pearson->SetLineColor(kRed);
}
