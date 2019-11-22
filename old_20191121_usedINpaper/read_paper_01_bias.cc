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

//#include "format.icc"

void read_paper_01_bias()
{
  //gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";
  
  int color_Poisson = kBlack;
  int color_Neyman  = kCyan;
  int color_Pearson = kBlue;
  int color_CNP     = kGreen;
  int color_Gauss   = kMagenta;
  
  /////////////////////////////////////////////////////////////////////////////
  
  roostr = "../data_11_paper_example1/out_chi2_mu0015.000_det010_run001000_event010000.root";
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

  double xlow = -40-0.6666666666666666/2;
  double xhgh = 40-0.6666666666666666/2;
  int nbins = (xhgh-xlow)/0.6666666666666666/2;

  roostr = "h1_Poisson"; TH1D *h1_Poisson = new TH1D(roostr, roostr, nbins, xlow, xhgh);
  roostr = "h1_Neyman"; TH1D *h1_Neyman = new TH1D(roostr, roostr, nbins, xlow, xhgh);
  roostr = "h1_Pearson"; TH1D *h1_Pearson = new TH1D(roostr, roostr, nbins, xlow, xhgh);
  roostr = "h1_Gauss"; TH1D *h1_Gauss = new TH1D(roostr, roostr, nbins, xlow, xhgh);
  roostr = "h1_CNP"; TH1D *h1_CNP = new TH1D(roostr, roostr, nbins, xlow, xhgh);

  //entries = 1000000;

  double bias_Poisson = 0;
  double bias_Neyman = 0;
  double bias_Pearson = 0;
  double bias_Gauss = 0;
  double bias_CNP = 0;
  
  for(int ientry=0; ientry<entries; ientry++) {
    if(ientry%100000==0) cout<<TString::Format(" ---> processing %6.4f", ientry*1./entries )<<endl;
    
    tree_chi2->GetEntry( ientry );

    double relbais_Poisson = (mu_bestFit_Poisson-Nmu)*100./Nmu;
    double relbais_Neyman = (mu_bestFit_Neyman-Nmu)*100./Nmu;
    double relbais_Pearson = (mu_bestFit_Pearson-Nmu)*100./Nmu;
    double relbais_Gauss = (mu_bestFit_Gauss-Nmu)*100./Nmu;
    double relbais_CNP = (mu_bestFit_CNP-Nmu)*100./Nmu;

    h1_Poisson->Fill( relbais_Poisson );
    h1_Neyman->Fill( relbais_Neyman );
    h1_Pearson->Fill( relbais_Pearson );
    h1_Gauss->Fill( relbais_Gauss );
    h1_CNP->Fill( relbais_CNP );

    bias_Poisson += relbais_Poisson/entries;
    bias_Neyman += relbais_Neyman/entries;
    bias_Pearson += relbais_Pearson/entries;
    bias_Gauss += relbais_Gauss/entries;
    bias_CNP += relbais_CNP/entries;
    
  }

  cout<<endl;
  cout<<TString::Format(" bias Poisson %7.2f", bias_Poisson)<<endl;
  cout<<TString::Format(" bias Neyman  %7.2f", bias_Neyman)<<endl;
  cout<<TString::Format(" bias Pearson %7.2f", bias_Pearson)<<endl;
  cout<<TString::Format(" bias Gauss   %7.2f", bias_Gauss)<<endl;
  cout<<TString::Format(" bias CNP     %7.2f", bias_CNP)<<endl;
  cout<<endl;
  
  roostr = "canv_test";
  TCanvas *canv_test = new TCanvas(roostr, roostr, 900, 650);
  
  h1_Poisson->Draw();
  //h1_Poisson->SetLineColor(color_Poisson);
  
  h1_Neyman->Draw("same");
  //h1_Neyman->SetLineColor(color_Neyman);
  
  h1_Pearson->Draw("same");
  //h1_Pearson->SetLineColor(color_Pearson);
  
  h1_Gauss->Draw("same");
  //h1_Gauss->SetLineColor(color_Gauss);
  
  h1_CNP->Draw("same");
  //h1_CNP->SetLineColor(color_CNP);

  TFile *outfile = new TFile("example1_RelBias_Mu15_Ndet10.root", "recreate");
  h1_Poisson->Write();
  h1_Neyman->Write();
  h1_Pearson->Write();
  h1_Gauss->Write();
  h1_CNP->Write();
  outfile->Close();
  
}
