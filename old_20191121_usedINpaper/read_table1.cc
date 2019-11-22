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

void read_table1()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";

  const int user_Ndets = 18;
  int array_Ndet[user_Ndets] = {3,4,5,6,8,10,15,20,25,30,35,40,50,60,70,80,90,100};

  const int user_Nmus = 3;
  int array_Nmu[user_Nmus] = {15, 70, 150};
  
  //for(int idx=0; idx<user_Ndets; idx++) cout<<TString::Format(" %2d, %3d", idx, array_Ndet[idx])<<endl;

  
  
  TFile *roofile = new TFile("file_table1.root", "recreate");  
  TTree *tree_bias = new TTree("tree_bias", "tree_bias");

  int ndf_val = 0;
  double mutrue = 0;
  double muhat_Poisson = 0;
  double muhat_Neyman = 0;
  double muhat_Pearson = 0;
  double muhat_Gauss = 0;
  double muhat_CNP = 0;
  
  tree_bias->Branch( "ndf_val",  &ndf_val,  "ndf_val/I" );
  tree_bias->Branch( "mutrue",  &mutrue,  "mutrue/D" );
  tree_bias->Branch( "muhat_Poisson",  &muhat_Poisson,  "muhat_Poisson/D" );
  tree_bias->Branch( "muhat_Neyman",  &muhat_Neyman,  "muhat_Neyman/D" );
  tree_bias->Branch( "muhat_Pearson",  &muhat_Pearson,  "muhat_Pearson/D" );
  tree_bias->Branch( "muhat_Gauss",  &muhat_Gauss,  "muhat_Gauss/D" );
  tree_bias->Branch( "muhat_CNP",  &muhat_CNP,  "muhat_CNP/D" );


  ///////////////////////////////////////////////////////////////////////////////////

  
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
 
  
  for( int imu=0; imu<user_Nmus; imu++ ) {
    for( int idet=0; idet<user_Ndets; idet++ ) {
      
      mutrue = array_Nmu[imu];
      ndf_val = array_Ndet[idet];

      roostr = TString::Format("../data_hep/data_talbe1/out_chi2_mu%08.3f_det%03d_run%06d_event%06d.root", mutrue, ndf_val, 1000, 1000);
      TFile *file_input = new TFile(roostr, "read");
      TTree *tree_chi2 = (TTree*)file_input->Get("tree_chi2");

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

      cout<<endl<<roostr<<endl;
      cout<<" ---> Entries "<<entries<<endl<<endl;
      
      muhat_Poisson = 0;
      muhat_Neyman = 0;
      muhat_Pearson = 0;
      muhat_Gauss = 0;
      muhat_CNP = 0;
  
      for(int ientry=0; ientry<entries; ientry++) {
	tree_chi2->GetEntry(ientry);

	muhat_Poisson += mu_bestFit_Poisson/entries;
	muhat_Neyman += mu_bestFit_Neyman/entries;
	muhat_Pearson += mu_bestFit_Pearson/entries;
	muhat_Gauss += mu_bestFit_Gauss/entries;
	muhat_CNP += mu_bestFit_CNP/entries;
	
      }
      
      file_input->Close();

      tree_bias->Fill();
      
    }
  }

  roofile->cd();
  tree_bias->Write();
  roofile->Close();
  
}
  
