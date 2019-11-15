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

//void plot_FC_Coverage_Interval()
int main(int argc, char** argv)
{

  /////////////////////////////////////////////////////////////////

  double Nnu0 = 12; /// setting
  
  int Ndet = 10;
  int Nrun = 1000;
  int Nevent = 1000;
  
  /////////////////////////////////////////////////////////////////
  
  int from = 0;
  int to = 0;

  for(int i=1; i<argc; i++) {
  
    if( strcmp(argv[i],"-from")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>from ) ) { cerr<<" ---> Error from !"<<endl; exit(1); }
    }
  
    if( strcmp(argv[i],"-to")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>to ) ) { cerr<<" ---> Error to !"<<endl; exit(1); }
    }
  
    if( from+to ==0 ) {
      cerr<<" ERROR from+to=0"<<endl;
    }

  }
  
  cout<<endl<<" --- from/to "<<from<<"\t"<<to<<endl<<endl;
  
  
  //gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";
  
  int color_Poisson = kBlack;
  int color_Neyman  = kRed;
  int color_Pearson = kBlue;
  int color_CNP     = kGreen;
  int color_Gauss   = kOrange;

  ///////////////////////////////////////////// --->

  const double sigma68 = 0.6827;
  const double sigma90 = 0.9;
  const double sigma95 = 0.95;

  
  int nbins_deltaChi2 = 6000;
  double low_deltaChi2 = 0;
  double hgh_deltaChi2 = 24;
  
  roostr = "h1_deltaChi2_Poisson"; TH1D *h1_deltaChi2_Poisson = new TH1D(roostr, roostr, nbins_deltaChi2, low_deltaChi2, hgh_deltaChi2);
  roostr = "h1_deltaChi2_Neyman"; TH1D *h1_deltaChi2_Neyman = new TH1D(roostr, roostr, nbins_deltaChi2, low_deltaChi2, hgh_deltaChi2);
  roostr = "h1_deltaChi2_Pearson"; TH1D *h1_deltaChi2_Pearson = new TH1D(roostr, roostr, nbins_deltaChi2, low_deltaChi2, hgh_deltaChi2);
  roostr = "h1_deltaChi2_CNP"; TH1D *h1_deltaChi2_CNP = new TH1D(roostr, roostr, nbins_deltaChi2, low_deltaChi2, hgh_deltaChi2);
  roostr = "h1_deltaChi2_Gauss"; TH1D *h1_deltaChi2_Gauss = new TH1D(roostr, roostr, nbins_deltaChi2, low_deltaChi2, hgh_deltaChi2);
  
  ///////////////////////////////////////////// <---

  double mu_true = 0;

  double chi2_68_Poisson = 0;
  double chi2_68_Neyman = 0;
  double chi2_68_Pearson = 0;
  double chi2_68_CNP = 0;
  double chi2_68_Gauss = 0;

  double chi2_90_Poisson = 0;
  double chi2_90_Neyman = 0;
  double chi2_90_Pearson = 0;
  double chi2_90_CNP = 0;
  double chi2_90_Gauss = 0;
  
  double chi2_95_Poisson = 0;
  double chi2_95_Neyman = 0;
  double chi2_95_Pearson = 0;
  double chi2_95_CNP = 0;
  double chi2_95_Gauss = 0;
  
  TTree *tree_fc = new TTree("tree_fc", "tree_fc");

  tree_fc->Branch( "mu_true", &mu_true, "mu_true/D" );

  tree_fc->Branch( "chi2_68_Poisson", &chi2_68_Poisson, "chi2_68_Poisson/D" );
  tree_fc->Branch( "chi2_68_Neyman", &chi2_68_Neyman, "chi2_68_Neyman/D" );
  tree_fc->Branch( "chi2_68_Pearson", &chi2_68_Pearson, "chi2_68_Pearson/D" );
  tree_fc->Branch( "chi2_68_CNP", &chi2_68_CNP, "chi2_68_CNP/D" );
  tree_fc->Branch( "chi2_68_Gauss", &chi2_68_Gauss, "chi2_68_Gauss/D" );
  
  tree_fc->Branch( "chi2_90_Poisson", &chi2_90_Poisson, "chi2_90_Poisson/D" );
  tree_fc->Branch( "chi2_90_Neyman", &chi2_90_Neyman, "chi2_90_Neyman/D" );
  tree_fc->Branch( "chi2_90_Pearson", &chi2_90_Pearson, "chi2_90_Pearson/D" );
  tree_fc->Branch( "chi2_90_CNP", &chi2_90_CNP, "chi2_90_CNP/D" );
  tree_fc->Branch( "chi2_90_Gauss", &chi2_90_Gauss, "chi2_90_Gauss/D" );
  
  tree_fc->Branch( "chi2_95_Poisson", &chi2_95_Poisson, "chi2_95_Poisson/D" );
  tree_fc->Branch( "chi2_95_Neyman", &chi2_95_Neyman, "chi2_95_Neyman/D" );
  tree_fc->Branch( "chi2_95_Pearson", &chi2_95_Pearson, "chi2_95_Pearson/D" );
  tree_fc->Branch( "chi2_95_CNP", &chi2_95_CNP, "chi2_95_CNP/D" );
  tree_fc->Branch( "chi2_95_Gauss", &chi2_95_Gauss, "chi2_95_Gauss/D" );
  
  ///////////////////////////////////////////// --->

  for(int idx=from; idx<=to; idx++) {

    h1_deltaChi2_Poisson->Reset();
    h1_deltaChi2_Neyman->Reset();
    h1_deltaChi2_Pearson->Reset();
    h1_deltaChi2_CNP->Reset();
    h1_deltaChi2_Gauss->Reset();
    
    ///////////////////
    
    double Nnu_user = Nnu0 + (idx-1)*0.005;
    
    roostr = TString::Format("./data_03_FC/out_chi2_mu%08.3f_det%03d_run%06d_event%06d.root",
			     Nnu_user, Ndet, Nrun, Nevent);
    
    ifstream check_file(roostr);

    if(!check_file) {
      cerr<<" ---> Warning, no file: "<<roostr<<endl;
      continue;
    }

    cout<<TString::Format(" ---> processing %4d, ", idx)<<roostr<<endl;

    ////////////////////////////////////////////////////////
    
    TFile *roofile = new TFile(roostr, "read");
    TTree *tree_chi2 = (TTree*)roofile->Get("tree_chi2");

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

    for(int ientry=0; ientry<=entries; ientry++) {
      tree_chi2->GetEntry(ientry);
      
      h1_deltaChi2_Poisson->Fill( chi2true_Poisson-chi2bestFit_Poisson );
      h1_deltaChi2_Neyman->Fill( chi2true_Neyman-chi2bestFit_Neyman );
      h1_deltaChi2_Pearson->Fill( chi2true_Pearson-chi2bestFit_Pearson );
      h1_deltaChi2_CNP->Fill( chi2true_CNP-chi2bestFit_CNP );
      h1_deltaChi2_Gauss->Fill( chi2true_Gauss-chi2bestFit_Gauss );
    }// ientry

    h1_deltaChi2_Poisson->Scale(1./entries);
    h1_deltaChi2_Neyman->Scale(1./entries);
    h1_deltaChi2_Pearson->Scale(1./entries);
    h1_deltaChi2_CNP->Scale(1./entries);
    h1_deltaChi2_Gauss->Scale(1./entries);

    //////////////////////////////////////////////////////
    
    chi2_68_Poisson = 0;
    chi2_90_Poisson = 0;
    chi2_68_Neyman = 0;
    chi2_90_Neyman = 0;
    chi2_68_Pearson = 0;
    chi2_90_Pearson = 0;
    chi2_68_CNP = 0;
    chi2_90_CNP = 0;
    chi2_68_Gauss = 0;
    chi2_90_Gauss = 0;
  
    chi2_95_Poisson = 0;
    chi2_95_Neyman = 0;
    chi2_95_Pearson = 0;
    chi2_95_CNP = 0;
    chi2_95_Gauss = 0;
    
    mu_true = Nmu;

    int flag_68 = 1;
    int flag_90 = 1;
    
    //////////////////////////// Poisson
    flag_68 = 1;
    flag_90 = 1;
    
    for(int ibin=1; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Poisson->Integral(1, ibin);
      if( content>sigma68 ) {
	chi2_68_Poisson = h1_deltaChi2_Poisson->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_68 = ibin-1;
	break;
      }      
    }   
    for(int ibin=flag_68; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Poisson->Integral(1, ibin);
      if( content>sigma90 ) {
	chi2_90_Poisson = h1_deltaChi2_Poisson->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_90 = ibin-1;
	break;
      }      
    }
    for(int ibin=flag_90; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Poisson->Integral(1, ibin);
      if( content>sigma95 ) {
	chi2_95_Poisson = h1_deltaChi2_Poisson->GetXaxis()->GetBinUpEdge( ibin-1 );
	break;
      }      
    }

    //////////////////////////// Neyman
    flag_68 = 1;
    flag_90 = 1;
    
    for(int ibin=1; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Neyman->Integral(1, ibin);
      if( content>sigma68 ) {
	chi2_68_Neyman = h1_deltaChi2_Neyman->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_68 = ibin-1;
	break;
      }      
    }   
    for(int ibin=flag_68; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Neyman->Integral(1, ibin);
      if( content>sigma90 ) {
	chi2_90_Neyman = h1_deltaChi2_Neyman->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_90 = ibin-1;
	break;
      }      
    }
    for(int ibin=flag_90; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Neyman->Integral(1, ibin);
      if( content>sigma95 ) {
	chi2_95_Neyman = h1_deltaChi2_Neyman->GetXaxis()->GetBinUpEdge( ibin-1 );
	break;
      }      
    }

    //////////////////////////// Pearson
    flag_68 = 1;
    flag_90 = 1;
    
    for(int ibin=1; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Pearson->Integral(1, ibin);
      if( content>sigma68 ) {
	chi2_68_Pearson = h1_deltaChi2_Pearson->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_68 = ibin-1;
	break;
      }      
    }   
    for(int ibin=flag_68; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Pearson->Integral(1, ibin);
      if( content>sigma90 ) {
	chi2_90_Pearson = h1_deltaChi2_Pearson->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_90 = ibin-1;
	break;
      }      
    }
    for(int ibin=flag_90; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Pearson->Integral(1, ibin);
      if( content>sigma95 ) {
	chi2_95_Pearson = h1_deltaChi2_Pearson->GetXaxis()->GetBinUpEdge( ibin-1 );
	break;
      }      
    }

    //////////////////////////// CNP
    flag_68 = 1;
    flag_90 = 1;
    
    for(int ibin=1; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_CNP->Integral(1, ibin);
      if( content>sigma68 ) {
	chi2_68_CNP = h1_deltaChi2_CNP->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_68 = ibin-1;
	break;
      }      
    }   
    for(int ibin=flag_68; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_CNP->Integral(1, ibin);
      if( content>sigma90 ) {
	chi2_90_CNP = h1_deltaChi2_CNP->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_90 = ibin-1;
	break;
      }      
    }
    for(int ibin=flag_90; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_CNP->Integral(1, ibin);
      if( content>sigma95 ) {
	chi2_95_CNP = h1_deltaChi2_CNP->GetXaxis()->GetBinUpEdge( ibin-1 );
	break;
      }      
    }

    //////////////////////////// Gauss
    flag_68 = 1;
    flag_90 = 1;
    
    for(int ibin=1; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Gauss->Integral(1, ibin);
      if( content>sigma68 ) {
	chi2_68_Gauss = h1_deltaChi2_Gauss->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_68 = ibin-1;
	break;
      }      
    }   
    for(int ibin=flag_68; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Gauss->Integral(1, ibin);
      if( content>sigma90 ) {
	chi2_90_Gauss = h1_deltaChi2_Gauss->GetXaxis()->GetBinUpEdge( ibin-1 );
	flag_90 = ibin-1;
	break;
      }      
    }
    for(int ibin=flag_90; ibin<=nbins_deltaChi2; ibin++) {
      double content = h1_deltaChi2_Gauss->Integral(1, ibin);
      if( content>sigma95 ) {
	chi2_95_Gauss = h1_deltaChi2_Gauss->GetXaxis()->GetBinUpEdge( ibin-1 );
	break;
      }      
    }




    
    tree_fc->Fill();
      
  }// idx for file


  roostr = TString::Format("out_FC_%06d_%06d.root", from, to);
  TFile *roofile = new TFile(roostr, "recreate");
  tree_fc->Write();
  roofile->Close();
    
  return 0;
}
