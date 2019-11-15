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

void plot_FC_LC()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";
  
  int color_Poisson = kBlack;
  int color_Neyman  = kMagenta;
  int color_Pearson = kBlue;
  int color_CNP     = kGreen;
  int color_Gauss   = kOrange;
  
  /////////////////////////////////////////////////////////////////////////////

  roostr = "./data_03_FC_Ndet10/out_FC_CL_det10_intervals.root";
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
  Double_t        CL68_low_Poisson;
  Double_t        CL68_hgh_Poisson;
  Double_t        CL68_low_Neyman;
  Double_t        CL68_hgh_Neyman;
  Double_t        CL68_low_Pearson;
  Double_t        CL68_hgh_Pearson;
  Double_t        CL68_low_CNP;
  Double_t        CL68_hgh_CNP;
  Double_t        CL68_low_Gauss;
  Double_t        CL68_hgh_Gauss;
  Double_t        CL90_low_Poisson;
  Double_t        CL90_hgh_Poisson;
  Double_t        CL90_low_Neyman;
  Double_t        CL90_hgh_Neyman;
  Double_t        CL90_low_Pearson;
  Double_t        CL90_hgh_Pearson;
  Double_t        CL90_low_CNP;
  Double_t        CL90_hgh_CNP;
  Double_t        CL90_low_Gauss;
  Double_t        CL90_hgh_Gauss;
  Double_t        CL95_low_Poisson;
  Double_t        CL95_hgh_Poisson;
  Double_t        CL95_low_Neyman;
  Double_t        CL95_hgh_Neyman;
  Double_t        CL95_low_Pearson;
  Double_t        CL95_hgh_Pearson;
  Double_t        CL95_low_CNP;
  Double_t        CL95_hgh_CNP;
  Double_t        CL95_low_Gauss;
  Double_t        CL95_hgh_Gauss;

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
  TBranch        *b_CL68_low_Poisson;   //!
  TBranch        *b_CL68_hgh_Poisson;   //!
  TBranch        *b_CL68_low_Neyman;   //!
  TBranch        *b_CL68_hgh_Neyman;   //!
  TBranch        *b_CL68_low_Pearson;   //!
  TBranch        *b_CL68_hgh_Pearson;   //!
  TBranch        *b_CL68_low_CNP;   //!
  TBranch        *b_CL68_hgh_CNP;   //!
  TBranch        *b_CL68_low_Gauss;   //!
  TBranch        *b_CL68_hgh_Gauss;   //!
  TBranch        *b_CL90_low_Poisson;   //!
  TBranch        *b_CL90_hgh_Poisson;   //!
  TBranch        *b_CL90_low_Neyman;   //!
  TBranch        *b_CL90_hgh_Neyman;   //!
  TBranch        *b_CL90_low_Pearson;   //!
  TBranch        *b_CL90_hgh_Pearson;   //!
  TBranch        *b_CL90_low_CNP;   //!
  TBranch        *b_CL90_hgh_CNP;   //!
  TBranch        *b_CL90_low_Gauss;   //!
  TBranch        *b_CL90_hgh_Gauss;   //!
  TBranch        *b_CL95_low_Poisson;   //!
  TBranch        *b_CL95_hgh_Poisson;   //!
  TBranch        *b_CL95_low_Neyman;   //!
  TBranch        *b_CL95_hgh_Neyman;   //!
  TBranch        *b_CL95_low_Pearson;   //!
  TBranch        *b_CL95_hgh_Pearson;   //!
  TBranch        *b_CL95_low_CNP;   //!
  TBranch        *b_CL95_hgh_CNP;   //!
  TBranch        *b_CL95_low_Gauss;   //!
  TBranch        *b_CL95_hgh_Gauss;   //!
  
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
  tree_chi2->SetBranchAddress("CL68_low_Poisson", &CL68_low_Poisson, &b_CL68_low_Poisson);
  tree_chi2->SetBranchAddress("CL68_hgh_Poisson", &CL68_hgh_Poisson, &b_CL68_hgh_Poisson);
  tree_chi2->SetBranchAddress("CL68_low_Neyman", &CL68_low_Neyman, &b_CL68_low_Neyman);
  tree_chi2->SetBranchAddress("CL68_hgh_Neyman", &CL68_hgh_Neyman, &b_CL68_hgh_Neyman);
  tree_chi2->SetBranchAddress("CL68_low_Pearson", &CL68_low_Pearson, &b_CL68_low_Pearson);
  tree_chi2->SetBranchAddress("CL68_hgh_Pearson", &CL68_hgh_Pearson, &b_CL68_hgh_Pearson);
  tree_chi2->SetBranchAddress("CL68_low_CNP", &CL68_low_CNP, &b_CL68_low_CNP);
  tree_chi2->SetBranchAddress("CL68_hgh_CNP", &CL68_hgh_CNP, &b_CL68_hgh_CNP);
  tree_chi2->SetBranchAddress("CL68_low_Gauss", &CL68_low_Gauss, &b_CL68_low_Gauss);
  tree_chi2->SetBranchAddress("CL68_hgh_Gauss", &CL68_hgh_Gauss, &b_CL68_hgh_Gauss);
  tree_chi2->SetBranchAddress("CL90_low_Poisson", &CL90_low_Poisson, &b_CL90_low_Poisson);
  tree_chi2->SetBranchAddress("CL90_hgh_Poisson", &CL90_hgh_Poisson, &b_CL90_hgh_Poisson);
  tree_chi2->SetBranchAddress("CL90_low_Neyman", &CL90_low_Neyman, &b_CL90_low_Neyman);
  tree_chi2->SetBranchAddress("CL90_hgh_Neyman", &CL90_hgh_Neyman, &b_CL90_hgh_Neyman);
  tree_chi2->SetBranchAddress("CL90_low_Pearson", &CL90_low_Pearson, &b_CL90_low_Pearson);
  tree_chi2->SetBranchAddress("CL90_hgh_Pearson", &CL90_hgh_Pearson, &b_CL90_hgh_Pearson);
  tree_chi2->SetBranchAddress("CL90_low_CNP", &CL90_low_CNP, &b_CL90_low_CNP);
  tree_chi2->SetBranchAddress("CL90_hgh_CNP", &CL90_hgh_CNP, &b_CL90_hgh_CNP);
  tree_chi2->SetBranchAddress("CL90_low_Gauss", &CL90_low_Gauss, &b_CL90_low_Gauss);
  tree_chi2->SetBranchAddress("CL90_hgh_Gauss", &CL90_hgh_Gauss, &b_CL90_hgh_Gauss);
  tree_chi2->SetBranchAddress("CL95_low_Poisson", &CL95_low_Poisson, &b_CL95_low_Poisson);
  tree_chi2->SetBranchAddress("CL95_hgh_Poisson", &CL95_hgh_Poisson, &b_CL95_hgh_Poisson);
  tree_chi2->SetBranchAddress("CL95_low_Neyman", &CL95_low_Neyman, &b_CL95_low_Neyman);
  tree_chi2->SetBranchAddress("CL95_hgh_Neyman", &CL95_hgh_Neyman, &b_CL95_hgh_Neyman);
  tree_chi2->SetBranchAddress("CL95_low_Pearson", &CL95_low_Pearson, &b_CL95_low_Pearson);
  tree_chi2->SetBranchAddress("CL95_hgh_Pearson", &CL95_hgh_Pearson, &b_CL95_hgh_Pearson);
  tree_chi2->SetBranchAddress("CL95_low_CNP", &CL95_low_CNP, &b_CL95_low_CNP);
  tree_chi2->SetBranchAddress("CL95_hgh_CNP", &CL95_hgh_CNP, &b_CL95_hgh_CNP);
  tree_chi2->SetBranchAddress("CL95_low_Gauss", &CL95_low_Gauss, &b_CL95_low_Gauss);
  tree_chi2->SetBranchAddress("CL95_hgh_Gauss", &CL95_hgh_Gauss, &b_CL95_hgh_Gauss);

  int entries = tree_chi2->GetEntries();
  entries = 500000;
  cout<<endl<<" ---> Entries "<<entries<<endl<<endl;

  int bins_CL68 = 100;
  double low_CL68 = 2;
  double hgh_CL68 = 4;
  double average_CL68_Poisson = 0;
  roostr = "h1_CL68_Poisson"; TH1D *h1_CL68_Poisson = new TH1D(roostr, roostr, bins_CL68, low_CL68, hgh_CL68);
  double average_CL68_Neyman = 0;
  roostr = "h1_CL68_Neyman"; TH1D *h1_CL68_Neyman = new TH1D(roostr, roostr, bins_CL68, low_CL68, hgh_CL68);
  double average_CL68_Pearson = 0;
  roostr = "h1_CL68_Pearson"; TH1D *h1_CL68_Pearson = new TH1D(roostr, roostr, bins_CL68, low_CL68, hgh_CL68);
  double average_CL68_CNP = 0;
  roostr = "h1_CL68_CNP"; TH1D *h1_CL68_CNP = new TH1D(roostr, roostr, bins_CL68, low_CL68, hgh_CL68);
  double average_CL68_Gauss = 0;
  roostr = "h1_CL68_Gauss"; TH1D *h1_CL68_Gauss = new TH1D(roostr, roostr, bins_CL68, low_CL68, hgh_CL68);
  
  for(int ientry=0; ientry<entries; ientry++) {
    if(ientry%100000==0) cout<<TString::Format(" ---> Processing %6.4f", ientry*1./entries)<<endl;    
    tree_chi2->GetEntry(ientry);

    double CL68_Poisson = CL68_hgh_Poisson - CL68_low_Poisson;
    h1_CL68_Poisson->Fill( CL68_Poisson );
    average_CL68_Poisson += CL68_Poisson/entries;
    double CL68_Neyman = CL68_hgh_Neyman - CL68_low_Neyman;
    h1_CL68_Neyman->Fill( CL68_Neyman );
    average_CL68_Neyman += CL68_Neyman/entries;
    double CL68_Pearson = CL68_hgh_Pearson - CL68_low_Pearson;
    h1_CL68_Pearson->Fill( CL68_Pearson );
    average_CL68_Pearson += CL68_Pearson/entries;
    double CL68_CNP = CL68_hgh_CNP - CL68_low_CNP;
    h1_CL68_CNP->Fill( CL68_CNP );
    average_CL68_CNP += CL68_CNP/entries;
    double CL68_Gauss = CL68_hgh_Gauss - CL68_low_Gauss;
    h1_CL68_Gauss->Fill( CL68_Gauss );
    average_CL68_Gauss += CL68_Gauss/entries;
    
  }

  cout<<endl;
  cout<<TString::Format(" Poisson 68C.L. %6.3f", average_CL68_Poisson)<<endl;
  cout<<TString::Format("     CNP 68C.L. %6.3f", average_CL68_CNP)<<endl;
  cout<<TString::Format("   Gauss 68C.L. %6.3f", average_CL68_Gauss)<<endl;
  cout<<TString::Format("  Neyman 68C.L. %6.3f", average_CL68_Neyman)<<endl;
  cout<<TString::Format(" Pearson 68C.L. %6.3f", average_CL68_Pearson)<<endl;
  cout<<endl;

  double ymax = 0;
  double ymax_Poisson = h1_CL68_Poisson->GetMaximum();
  if( ymax_Poisson>ymax ) ymax = ymax_Poisson;
  double ymax_Neyman = h1_CL68_Neyman->GetMaximum();
  if( ymax_Neyman>ymax ) ymax = ymax_Neyman;
  double ymax_Pearson = h1_CL68_Pearson->GetMaximum();
  if( ymax_Pearson>ymax ) ymax = ymax_Pearson;
  double ymax_CNP = h1_CL68_CNP->GetMaximum();
  if( ymax_CNP>ymax ) ymax = ymax_CNP;
  double ymax_Gauss = h1_CL68_Gauss->GetMaximum();
  if( ymax_Gauss>ymax ) ymax = ymax_Gauss;
  
  roostr = "canv_CL68";
  TCanvas *canv_CL68 = new TCanvas(roostr, roostr, 680, 650);
  roo_canv_margin(canv_CL68, 0.18, 0.1, 0.18, 0.06);  
  h1_CL68_Poisson->Draw();  
  h1_CL68_Poisson->SetLineColor(color_Poisson);
  roo_title(h1_CL68_Poisson, "68% C.L. interval", "Entries");
  roo_center_title(h1_CL68_Poisson);
  roo_offset(h1_CL68_Poisson, 1.1, 1.3);
  
  // h1_CL68_Neyman->Draw("same");  
  // h1_CL68_Neyman->SetLineColor(color_Neyman);

  // h1_CL68_Pearson->Draw("same");  
  // h1_CL68_Pearson->SetLineColor(color_Pearson);

  // h1_CL68_CNP->Draw("same");  
  // h1_CL68_CNP->SetLineColor(color_CNP);

  // h1_CL68_Gauss->Draw("same");  
  // h1_CL68_Gauss->SetLineColor(color_Gauss);

  h1_CL68_Poisson->Draw("same axis");  


  h1_CL68_Poisson->GetYaxis()->SetRangeUser(0.5, ymax * 1.4);
  
  // TLegend *lg_CL68 = new TLegend(0.55-0.01,0.56,0.83+0.07-0.01,0.88);
  // lg_CL68->SetBorderSize(0);
  // lg_CL68->AddEntry(  h1_CL68_Poisson, TString::Format("Poisson  %5.3f",  average_CL68_Poisson), "l" );
  // lg_CL68->AddEntry(  h1_CL68_CNP,     TString::Format("CNP       %5.3f",  average_CL68_CNP), "l" );
  // lg_CL68->AddEntry(  h1_CL68_Gauss,   TString::Format("Gauss    %5.3f",  average_CL68_Gauss), "l" );
  // lg_CL68->AddEntry(  h1_CL68_Pearson, TString::Format("Pearson %5.3f",  average_CL68_Pearson), "l" );    
  // lg_CL68->AddEntry(  h1_CL68_Neyman,  TString::Format("Neyman %5.3f",  average_CL68_Neyman), "l" );
		    
  // lg_CL68->Draw("same");
  // lg_CL68->SetTextSize(0.05);
  // lg_CL68->SetTextFont(42);

}
