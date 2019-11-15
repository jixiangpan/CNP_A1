#include<iostream>
#include<fstream>
#include<cmath>
#include "stdlib.h" 
using namespace std;

#include<map>
#include<set>

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

double chi2_Poisson(double mu_pred, map<int,double> array_val, int ndet)
{
  double chi2 = 0;

  for(int idet=1; idet<=ndet; idet++) {
    double val_m = array_val[idet];

    if( val_m==0 ) chi2 += 2*mu_pred;
    else chi2 += 2*( mu_pred -val_m +val_m*log(val_m/mu_pred) );
  }
  
  return chi2;
}

double chi2_Neyman(double mu_pred, map<int,double> array_val, int ndet)
{
  double chi2 = 0;

  for(int idet=1; idet<=ndet; idet++) {
    double val_m = array_val[idet];

    if( val_m==0 ) chi2 += 2*mu_pred;
    else chi2 += pow(mu_pred-val_m, 2)/val_m;
  }
  
  return chi2;
}

double chi2_Pearson(double mu_pred, map<int,double> array_val, int ndet)
{
  double chi2 = 0;

  for(int idet=1; idet<=ndet; idet++) {
    double val_m = array_val[idet];

    chi2 += pow(mu_pred-val_m, 2)/mu_pred;
  }
  
  return chi2;
}

double chi2_CNP(double mu_pred, map<int,double> array_val, int ndet)
{
  double chi2 = 0;

  for(int idet=1; idet<=ndet; idet++) {
    double val_m = array_val[idet];

    if( val_m==0 ) chi2 += 2*mu_pred;
    else chi2 += pow( mu_pred-val_m, 2 )/3 * ( 1./val_m + 2./mu_pred );
  }
  
  return chi2;
}

double chi2_Gauss(double mu_pred, map<int,double> array_val, int ndet)
{
  double chi2 = 0;

  for(int idet=1; idet<=ndet; idet++) {
    double val_m = array_val[idet];

    double m_saturated = sqrt(4*val_m*val_m + 1)/2 - 0.5;
    
    if( val_m==0 ) chi2 += 2*mu_pred;
    else chi2 += ( pow(val_m-mu_pred,2)/mu_pred -pow(val_m-m_saturated,2)/m_saturated + log(mu_pred/m_saturated) );;
  }
  
  return chi2;
}




void read_FC_CL(double Nmu, int Ndet, int Nrun, int Nevent, int fflag)
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  TString roostr = "";
  
  if( Nmu==0 || Ndet==0 || Nrun==0 || Nevent==0 ) {
    cout<<endl<<" ---> ERROR , 0 in inputs "<<endl<<endl;
    exit(1);
  }

  cout<<endl;
  cout<<TString::Format(" ---> Nmu    %8.3f", Nmu)<<endl;
  cout<<TString::Format(" ---> Ndet   %8d", Ndet)<<endl;
  cout<<TString::Format(" ---> Nrun   %8d", Nrun)<<endl;
  cout<<TString::Format(" ---> Nevent %8d", Nevent)<<endl;
  cout<<endl;

  //////////////////////////////////////////////////////////////////////////////////////

  int line_outScaned_range = 0;
  
  TFile *file_fc = new TFile("./data_03_FC_Ndet10/out_FC_Ndet10.root", "read");
  TTree *tree_fc = (TTree*)file_fc->Get("tree_fc");  

  // Declaration of leaf types
  Double_t        mu_true;
  Double_t        chi2_68_Poisson;
  Double_t        chi2_68_Neyman;
  Double_t        chi2_68_Pearson;
  Double_t        chi2_68_CNP;
  Double_t        chi2_68_Gauss;
  Double_t        chi2_90_Poisson;
  Double_t        chi2_90_Neyman;
  Double_t        chi2_90_Pearson;
  Double_t        chi2_90_CNP;
  Double_t        chi2_90_Gauss;
  Double_t        chi2_95_Poisson;
  Double_t        chi2_95_Neyman;
  Double_t        chi2_95_Pearson;
  Double_t        chi2_95_CNP;
  Double_t        chi2_95_Gauss;

  // List of branches
  TBranch        *b_mu_true;   //!
  TBranch        *b_chi2_68_Poisson;   //!
  TBranch        *b_chi2_68_Neyman;   //!
  TBranch        *b_chi2_68_Pearson;   //!
  TBranch        *b_chi2_68_CNP;   //!
  TBranch        *b_chi2_68_Gauss;   //!
  TBranch        *b_chi2_90_Poisson;   //!
  TBranch        *b_chi2_90_Neyman;   //!
  TBranch        *b_chi2_90_Pearson;   //!
  TBranch        *b_chi2_90_CNP;   //!
  TBranch        *b_chi2_90_Gauss;   //!
  TBranch        *b_chi2_95_Poisson;   //!
  TBranch        *b_chi2_95_Neyman;   //!
  TBranch        *b_chi2_95_Pearson;   //!
  TBranch        *b_chi2_95_CNP;   //!
  TBranch        *b_chi2_95_Gauss;   //!
 
  // Set branch addresses and branch pointers
  tree_fc->SetBranchAddress("mu_true", &mu_true, &b_mu_true);
  tree_fc->SetBranchAddress("chi2_68_Poisson", &chi2_68_Poisson, &b_chi2_68_Poisson);
  tree_fc->SetBranchAddress("chi2_68_Neyman", &chi2_68_Neyman, &b_chi2_68_Neyman);
  tree_fc->SetBranchAddress("chi2_68_Pearson", &chi2_68_Pearson, &b_chi2_68_Pearson);
  tree_fc->SetBranchAddress("chi2_68_CNP", &chi2_68_CNP, &b_chi2_68_CNP);
  tree_fc->SetBranchAddress("chi2_68_Gauss", &chi2_68_Gauss, &b_chi2_68_Gauss);
  tree_fc->SetBranchAddress("chi2_90_Poisson", &chi2_90_Poisson, &b_chi2_90_Poisson);
  tree_fc->SetBranchAddress("chi2_90_Neyman", &chi2_90_Neyman, &b_chi2_90_Neyman);
  tree_fc->SetBranchAddress("chi2_90_Pearson", &chi2_90_Pearson, &b_chi2_90_Pearson);
  tree_fc->SetBranchAddress("chi2_90_CNP", &chi2_90_CNP, &b_chi2_90_CNP);
  tree_fc->SetBranchAddress("chi2_90_Gauss", &chi2_90_Gauss, &b_chi2_90_Gauss);
  tree_fc->SetBranchAddress("chi2_95_Poisson", &chi2_95_Poisson, &b_chi2_95_Poisson);
  tree_fc->SetBranchAddress("chi2_95_Neyman", &chi2_95_Neyman, &b_chi2_95_Neyman);
  tree_fc->SetBranchAddress("chi2_95_Pearson", &chi2_95_Pearson, &b_chi2_95_Pearson);
  tree_fc->SetBranchAddress("chi2_95_CNP", &chi2_95_CNP, &b_chi2_95_CNP);
  tree_fc->SetBranchAddress("chi2_95_Gauss", &chi2_95_Gauss, &b_chi2_95_Gauss); 

  int entries_fc = tree_fc->GetEntries();
  map<int, double>fc_mu_true;
  
  map<int, double>fc_sigma68_Poisson;
  map<int, double>fc_sigma68_Neyman;
  map<int, double>fc_sigma68_Pearson;
  map<int, double>fc_sigma68_CNP;
  map<int, double>fc_sigma68_Gauss;

  map<int, double>fc_sigma90_Poisson;
  map<int, double>fc_sigma90_Neyman;
  map<int, double>fc_sigma90_Pearson;
  map<int, double>fc_sigma90_CNP;
  map<int, double>fc_sigma90_Gauss;

  map<int, double>fc_sigma95_Poisson;
  map<int, double>fc_sigma95_Neyman;
  map<int, double>fc_sigma95_Pearson;
  map<int, double>fc_sigma95_CNP;
  map<int, double>fc_sigma95_Gauss;

  for(int ientry=0; ientry<entries_fc; ientry++) {
    tree_fc->GetEntry(ientry);

    fc_mu_true[ientry] = mu_true;

    fc_sigma68_Poisson[ientry] = chi2_68_Poisson;
    fc_sigma68_Neyman[ientry]  = chi2_68_Neyman;
    fc_sigma68_Pearson[ientry] = chi2_68_Pearson;
    fc_sigma68_CNP[ientry]     = chi2_68_CNP;
    fc_sigma68_Gauss[ientry]   = chi2_68_Gauss;

    fc_sigma90_Poisson[ientry] = chi2_90_Poisson;
    fc_sigma90_Neyman[ientry]  = chi2_90_Neyman;
    fc_sigma90_Pearson[ientry] = chi2_90_Pearson;
    fc_sigma90_CNP[ientry]     = chi2_90_CNP;
    fc_sigma90_Gauss[ientry]   = chi2_90_Gauss;

    fc_sigma95_Poisson[ientry] = chi2_95_Poisson;
    fc_sigma95_Neyman[ientry]  = chi2_95_Neyman;
    fc_sigma95_Pearson[ientry] = chi2_95_Pearson;
    fc_sigma95_CNP[ientry]     = chi2_95_CNP;
    fc_sigma95_Gauss[ientry]   = chi2_95_Gauss;

  }
  
  file_fc->Close();

  
  //////////////////////////////////////// Tree
  
  double mu_bestFit_Poisson = 0;
  double mu_bestFit_Neyman  = 0;
  double mu_bestFit_Pearson = 0;
  double mu_bestFit_CNP     = 0;
  double mu_bestFit_Gauss   = 0;

  double chi2true_Poisson = 0;
  double chi2true_Neyman  = 0;
  double chi2true_Pearson = 0;
  double chi2true_CNP     = 0;
  double chi2true_Gauss   = 0;

  double chi2bestFit_Poisson = 0;
  double chi2bestFit_Neyman  = 0;
  double chi2bestFit_Pearson = 0;
  double chi2bestFit_CNP     = 0;
  double chi2bestFit_Gauss   = 0;

  double CL68_low_Poisson = 0;
  double CL68_hgh_Poisson = 0;
  double CL68_low_Neyman = 0;
  double CL68_hgh_Neyman = 0;
  double CL68_low_Pearson = 0;
  double CL68_hgh_Pearson = 0;
  double CL68_low_CNP = 0;
  double CL68_hgh_CNP = 0;
  double CL68_low_Gauss = 0;
  double CL68_hgh_Gauss = 0;

  double CL90_low_Poisson = 0;
  double CL90_hgh_Poisson = 0;
  double CL90_low_Neyman = 0;
  double CL90_hgh_Neyman = 0;
  double CL90_low_Pearson = 0;
  double CL90_hgh_Pearson = 0;
  double CL90_low_CNP = 0;
  double CL90_hgh_CNP = 0;
  double CL90_low_Gauss = 0;
  double CL90_hgh_Gauss = 0;
  
  double CL95_low_Poisson = 0;
  double CL95_hgh_Poisson = 0;
  double CL95_low_Neyman = 0;
  double CL95_hgh_Neyman = 0;
  double CL95_low_Pearson = 0;
  double CL95_hgh_Pearson = 0;
  double CL95_low_CNP = 0;
  double CL95_hgh_CNP = 0;
  double CL95_low_Gauss = 0;
  double CL95_hgh_Gauss = 0;
  
  roostr = TString::Format("out_FC_CL_mu%08.3f_det%03d_run%06d_event%06d_%02d.root", Nmu, Ndet, Nrun, Nevent, fflag);
  TFile *out_file = new TFile(roostr, "recreate");
  out_file->cd();
    
  TTree *tree_chi2 = new TTree("tree_chi2", "tree_chi2");
  
  tree_chi2->Branch( "Nmu",  &Nmu,  "Nmu/D" );
  tree_chi2->Branch( "Ndet", &Ndet, "Ndet/I" );

  tree_chi2->Branch( "mu_bestFit_Poisson", &mu_bestFit_Poisson, "mu_bestFit_Poisson/D" );
  tree_chi2->Branch( "mu_bestFit_Neyman",  &mu_bestFit_Neyman,  "mu_bestFit_Neyman/D" );
  tree_chi2->Branch( "mu_bestFit_Pearson", &mu_bestFit_Pearson, "mu_bestFit_Pearson/D" );
  tree_chi2->Branch( "mu_bestFit_CNP",     &mu_bestFit_CNP,     "mu_bestFit_CNP/D" );
  tree_chi2->Branch( "mu_bestFit_Gauss",   &mu_bestFit_Gauss,   "mu_bestFit_Gauss/D" );
  
  tree_chi2->Branch( "chi2true_Poisson", &chi2true_Poisson, "chi2true_Poisson/D" );
  tree_chi2->Branch( "chi2true_Neyman",  &chi2true_Neyman,  "chi2true_Neyman/D" );
  tree_chi2->Branch( "chi2true_Pearson", &chi2true_Pearson, "chi2true_Pearson/D" );
  tree_chi2->Branch( "chi2true_CNP",     &chi2true_CNP,     "chi2true_CNP/D" );
  tree_chi2->Branch( "chi2true_Gauss",   &chi2true_Gauss,   "chi2true_Gauss/D" );
  
  tree_chi2->Branch( "chi2bestFit_Poisson", &chi2bestFit_Poisson, "chi2bestFit_Poisson/D" );
  tree_chi2->Branch( "chi2bestFit_Neyman",  &chi2bestFit_Neyman,  "chi2bestFit_Neyman/D" );
  tree_chi2->Branch( "chi2bestFit_Pearson", &chi2bestFit_Pearson, "chi2bestFit_Pearson/D" );
  tree_chi2->Branch( "chi2bestFit_CNP",     &chi2bestFit_CNP,     "chi2bestFit_CNP/D" );
  tree_chi2->Branch( "chi2bestFit_Gauss",   &chi2bestFit_Gauss,   "chi2bestFit_Gauss/D" );

  tree_chi2->Branch( "CL68_low_Poisson", &CL68_low_Poisson, "CL68_low_Poisson/D" );
  tree_chi2->Branch( "CL68_hgh_Poisson", &CL68_hgh_Poisson, "CL68_hgh_Poisson/D" );
  tree_chi2->Branch( "CL68_low_Neyman",  &CL68_low_Neyman,  "CL68_low_Neyman/D" );
  tree_chi2->Branch( "CL68_hgh_Neyman",  &CL68_hgh_Neyman,  "CL68_hgh_Neyman/D" );
  tree_chi2->Branch( "CL68_low_Pearson", &CL68_low_Pearson, "CL68_low_Pearson/D" );
  tree_chi2->Branch( "CL68_hgh_Pearson", &CL68_hgh_Pearson, "CL68_hgh_Pearson/D" );
  tree_chi2->Branch( "CL68_low_CNP",     &CL68_low_CNP,     "CL68_low_CNP/D" );
  tree_chi2->Branch( "CL68_hgh_CNP",     &CL68_hgh_CNP,     "CL68_hgh_CNP/D" );
  tree_chi2->Branch( "CL68_low_Gauss",   &CL68_low_Gauss,   "CL68_low_Gauss/D" );
  tree_chi2->Branch( "CL68_hgh_Gauss",   &CL68_hgh_Gauss,   "CL68_hgh_Gauss/D" );
  
  tree_chi2->Branch( "CL90_low_Poisson", &CL90_low_Poisson, "CL90_low_Poisson/D" );
  tree_chi2->Branch( "CL90_hgh_Poisson", &CL90_hgh_Poisson, "CL90_hgh_Poisson/D" );
  tree_chi2->Branch( "CL90_low_Neyman",  &CL90_low_Neyman,  "CL90_low_Neyman/D" );
  tree_chi2->Branch( "CL90_hgh_Neyman",  &CL90_hgh_Neyman,  "CL90_hgh_Neyman/D" );
  tree_chi2->Branch( "CL90_low_Pearson", &CL90_low_Pearson, "CL90_low_Pearson/D" );
  tree_chi2->Branch( "CL90_hgh_Pearson", &CL90_hgh_Pearson, "CL90_hgh_Pearson/D" );
  tree_chi2->Branch( "CL90_low_CNP",     &CL90_low_CNP,     "CL90_low_CNP/D" );
  tree_chi2->Branch( "CL90_hgh_CNP",     &CL90_hgh_CNP,     "CL90_hgh_CNP/D" );
  tree_chi2->Branch( "CL90_low_Gauss",   &CL90_low_Gauss,   "CL90_low_Gauss/D" );
  tree_chi2->Branch( "CL90_hgh_Gauss",   &CL90_hgh_Gauss,   "CL90_hgh_Gauss/D" );
  
  tree_chi2->Branch( "CL95_low_Poisson", &CL95_low_Poisson, "CL95_low_Poisson/D" );
  tree_chi2->Branch( "CL95_hgh_Poisson", &CL95_hgh_Poisson, "CL95_hgh_Poisson/D" );
  tree_chi2->Branch( "CL95_low_Neyman",  &CL95_low_Neyman,  "CL95_low_Neyman/D" );
  tree_chi2->Branch( "CL95_hgh_Neyman",  &CL95_hgh_Neyman,  "CL95_hgh_Neyman/D" );
  tree_chi2->Branch( "CL95_low_Pearson", &CL95_low_Pearson, "CL95_low_Pearson/D" );
  tree_chi2->Branch( "CL95_hgh_Pearson", &CL95_hgh_Pearson, "CL95_hgh_Pearson/D" );
  tree_chi2->Branch( "CL95_low_CNP",     &CL95_low_CNP,     "CL95_low_CNP/D" );
  tree_chi2->Branch( "CL95_hgh_CNP",     &CL95_hgh_CNP,     "CL95_hgh_CNP/D" );
  tree_chi2->Branch( "CL95_low_Gauss",   &CL95_low_Gauss,   "CL95_low_Gauss/D" );
  tree_chi2->Branch( "CL95_hgh_Gauss",   &CL95_hgh_Gauss,   "CL95_hgh_Gauss/D" );
  
  //////////////////////////////////////// Toy
  
  map<int, double>array_val_m;
  
  for(int irun=1; irun<=Nrun; irun++) {

    if(irun%1==0) cout<<TString::Format(" ---> processing %6.4f", irun*1./Nrun)<<endl;
    
    TRandom3 *roo_random = new TRandom3(0);

    for(int ievent=1; ievent<=Nevent; ievent++) {
      array_val_m.clear();
        
      double val_sum_m = 0;
      double val_sum_inverse_m = 0;
      double val_sum_m2 = 0;
      int N_zero_m = 0;

      for(int idet=1; idet<=Ndet; idet++) {
	double val_m = 1* roo_random->Poisson( Nmu );
	array_val_m[idet] = val_m;
	  
	if( val_m==0 ) {
	  N_zero_m++;
	}
	else {
	  val_sum_m += val_m;
	  val_sum_inverse_m += 1./val_m;
	  val_sum_m2 += val_m*val_m;
	}
	
      }// idet

      //////
      chi2true_Poisson = chi2_Poisson( Nmu, array_val_m, Ndet );
      chi2true_Neyman  = chi2_Neyman( Nmu, array_val_m, Ndet );      
      chi2true_Pearson = chi2_Pearson( Nmu, array_val_m, Ndet );
      chi2true_CNP     = chi2_CNP( Nmu, array_val_m, Ndet );
      chi2true_Gauss   = chi2_Gauss( Nmu, array_val_m, Ndet );

      mu_bestFit_Poisson = val_sum_m/Ndet;
      mu_bestFit_Neyman = ( 2*(Ndet-N_zero_m) - Ndet )/val_sum_inverse_m;
      mu_bestFit_Pearson = sqrt( val_sum_m2/Ndet );
      mu_bestFit_Gauss = ( sqrt(4*val_sum_m2/Ndet +1) -1)/2;
      
      chi2bestFit_Poisson = chi2_Poisson( mu_bestFit_Poisson, array_val_m, Ndet );
      chi2bestFit_Neyman  = chi2_Neyman( mu_bestFit_Neyman, array_val_m, Ndet );      
      chi2bestFit_Pearson = chi2_Pearson( mu_bestFit_Pearson, array_val_m, Ndet );
      chi2bestFit_Gauss   = chi2_Gauss( mu_bestFit_Gauss, array_val_m, Ndet );

      ////// CNP
      mu_bestFit_CNP  = 0;
      chi2bestFit_CNP = 1e6;
	
      // P(x) = p[0] + p[1]*x + p[2] *x 2 + ... + p[n] x n
      double p0 = -1.*val_sum_m2;
      double p1 = 0;
      double p2 = 3*N_zero_m;
      double p3 = val_sum_inverse_m;
      double pars_polynomial[] = {p0,p1,p2,p3};
      
      int order_polynomial = 3;
      ROOT::Math::Polynomial *func_polynomial = new ROOT::Math::Polynomial(order_polynomial);
      func_polynomial->SetParameters( pars_polynomial );
      vector<double>solution_func_polynomial = func_polynomial->FindRealRoots();
      int size_solution_polynomial = solution_func_polynomial.size();
      delete func_polynomial;

      double chi2_min = 1e6;
      for(int idx=0; idx<size_solution_polynomial; idx++) {
	double b = solution_func_polynomial[idx];
	if(b<1e-4) b = 1e-4;
	double chi2_test = chi2_CNP(b, array_val_m, Ndet);
	if( chi2_test<0 ) chi2_test = 1e7;
	if( chi2_min>chi2_test ) chi2_min = chi2_test;
      }
      
      for(int idx=0; idx<size_solution_polynomial; idx++) {
	double b = solution_func_polynomial[idx];
	if(b<1e-4) b = 1e-4;
	double chi2_test = chi2_CNP(b, array_val_m, Ndet);
	if( chi2_min==chi2_test ) {	    
	  mu_bestFit_CNP  = b;
	  chi2bestFit_CNP = chi2_test;
	  break;
	}
      }

      // if( mu_bestFit_Poisson<Nmu-0.5 ) continue;
      // if( mu_bestFit_Poisson>Nmu+0.5 ) continue;

      double low_mu_bestFit = 4;
      double hgh_mu_bestFit = 25;
      
      if( mu_bestFit_Poisson<low_mu_bestFit ||
	  mu_bestFit_Poisson>hgh_mu_bestFit ||
	  mu_bestFit_Neyman<low_mu_bestFit ||
	  mu_bestFit_Neyman>hgh_mu_bestFit ||
	  mu_bestFit_Pearson<low_mu_bestFit ||
	  mu_bestFit_Pearson>hgh_mu_bestFit ||
	  mu_bestFit_CNP<low_mu_bestFit ||
	  mu_bestFit_CNP>hgh_mu_bestFit ||
	  mu_bestFit_Gauss<low_mu_bestFit ||
	  mu_bestFit_Gauss>hgh_mu_bestFit	  
	  ) {
	cout<<" -------------> WARNING: mu_bestFit out of the scaned range"<<endl;
	line_outScaned_range++;
	continue;
      }
	
      
      
      // cout<<endl;
      // cout<<TString::Format(" --> Poisson bestFit %8.3f", mu_bestFit_Poisson )<<endl;
      // cout<<TString::Format(" --> Neyman  bestFit %8.3f", mu_bestFit_Neyman )<<endl;
      // cout<<TString::Format(" --> Pearson bestFit %8.3f", mu_bestFit_Pearson )<<endl;
      // cout<<TString::Format(" --> CNP     bestFit %8.3f", mu_bestFit_CNP )<<endl;
      // cout<<TString::Format(" --> Gauss   bestFit %8.3f", mu_bestFit_Gauss )<<endl;
      
      ///////////////////////////////////////////////////// FC intervals

      CL68_low_Poisson = 0;
      CL68_hgh_Poisson = 0;
      CL68_low_Neyman = 0;
      CL68_hgh_Neyman = 0;
      CL68_low_Pearson = 0;
      CL68_hgh_Pearson = 0;
      CL68_low_CNP = 0;
      CL68_hgh_CNP = 0;
      CL68_low_Gauss = 0;
      CL68_hgh_Gauss = 0;

      CL90_low_Poisson = 0;
      CL90_hgh_Poisson = 0;
      CL90_low_Neyman = 0;
      CL90_hgh_Neyman = 0;
      CL90_low_Pearson = 0;
      CL90_hgh_Pearson = 0;
      CL90_low_CNP = 0;
      CL90_hgh_CNP = 0;
      CL90_low_Gauss = 0;
      CL90_hgh_Gauss = 0;
  
      CL95_low_Poisson = 0;
      CL95_hgh_Poisson = 0;
      CL95_low_Neyman = 0;
      CL95_hgh_Neyman = 0;
      CL95_low_Pearson = 0;
      CL95_hgh_Pearson = 0;
      CL95_low_CNP = 0;
      CL95_hgh_CNP = 0;
      CL95_low_Gauss = 0;
      CL95_hgh_Gauss = 0;
      
      ///////////////////////////////////////////////////////////////////////// CL68
      
      ///////////////////////////////// Poisson
      set<double>cl_sigma68_Poisson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Poisson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Poisson;
	if( deltaChi2<fc_sigma68_Poisson[idx] ) cl_sigma68_Poisson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma68_Poisson;
      for( it_cl_sigma68_Poisson=cl_sigma68_Poisson.begin(); it_cl_sigma68_Poisson!=cl_sigma68_Poisson.end(); it_cl_sigma68_Poisson++ ) {
	double result = *it_cl_sigma68_Poisson;
      }
      CL68_low_Poisson = *(cl_sigma68_Poisson.begin());
      CL68_hgh_Poisson = *(cl_sigma68_Poisson.rbegin());
      double width_CL68_Poisson = CL68_hgh_Poisson - CL68_low_Poisson;

      ///////////////////////////////// Neyman
      set<double>cl_sigma68_Neyman;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Neyman( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Neyman;
	if( deltaChi2<fc_sigma68_Neyman[idx] ) cl_sigma68_Neyman.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma68_Neyman;
      for( it_cl_sigma68_Neyman=cl_sigma68_Neyman.begin(); it_cl_sigma68_Neyman!=cl_sigma68_Neyman.end(); it_cl_sigma68_Neyman++ ) {
	double result = *it_cl_sigma68_Neyman;
      }
      CL68_low_Neyman = *(cl_sigma68_Neyman.begin());
      CL68_hgh_Neyman = *(cl_sigma68_Neyman.rbegin());
      double width_CL68_Neyman = CL68_hgh_Neyman - CL68_low_Neyman;

      ///////////////////////////////// Pearson
      set<double>cl_sigma68_Pearson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Pearson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Pearson;
	if( deltaChi2<fc_sigma68_Pearson[idx] ) cl_sigma68_Pearson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma68_Pearson;
      for( it_cl_sigma68_Pearson=cl_sigma68_Pearson.begin(); it_cl_sigma68_Pearson!=cl_sigma68_Pearson.end(); it_cl_sigma68_Pearson++ ) {
	double result = *it_cl_sigma68_Pearson;
      }
      CL68_low_Pearson = *(cl_sigma68_Pearson.begin());
      CL68_hgh_Pearson = *(cl_sigma68_Pearson.rbegin());
      double width_CL68_Pearson = CL68_hgh_Pearson - CL68_low_Pearson;

      ///////////////////////////////// CNP
      set<double>cl_sigma68_CNP;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_CNP( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_CNP;
	if( deltaChi2<fc_sigma68_CNP[idx] ) cl_sigma68_CNP.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma68_CNP;
      for( it_cl_sigma68_CNP=cl_sigma68_CNP.begin(); it_cl_sigma68_CNP!=cl_sigma68_CNP.end(); it_cl_sigma68_CNP++ ) {
	double result = *it_cl_sigma68_CNP;
      }
      CL68_low_CNP = *(cl_sigma68_CNP.begin());
      CL68_hgh_CNP = *(cl_sigma68_CNP.rbegin());
      double width_CL68_CNP = CL68_hgh_CNP - CL68_low_CNP;

      ///////////////////////////////// Gauss
      set<double>cl_sigma68_Gauss;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Gauss( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Gauss;
	if( deltaChi2<fc_sigma68_Gauss[idx] ) cl_sigma68_Gauss.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma68_Gauss;
      for( it_cl_sigma68_Gauss=cl_sigma68_Gauss.begin(); it_cl_sigma68_Gauss!=cl_sigma68_Gauss.end(); it_cl_sigma68_Gauss++ ) {
	double result = *it_cl_sigma68_Gauss;
      }
      CL68_low_Gauss = *(cl_sigma68_Gauss.begin());
      CL68_hgh_Gauss = *(cl_sigma68_Gauss.rbegin());
      double width_CL68_Gauss = CL68_hgh_Gauss - CL68_low_Gauss;

      // cout<<endl;
      // cout<<TString::Format( " ---> 68 C.L. Poisson  %8.3f %8.3f, width %8.3f", CL68_low_Poisson, CL68_hgh_Poisson, width_CL68_Poisson)<<endl;      
      // cout<<TString::Format( " ---> 68 C.L. Neyman   %8.3f %8.3f, width %8.3f", CL68_low_Neyman, CL68_hgh_Neyman, width_CL68_Neyman)<<endl;      
      // cout<<TString::Format( " ---> 68 C.L. Pearson  %8.3f %8.3f, width %8.3f", CL68_low_Pearson, CL68_hgh_Pearson, width_CL68_Pearson)<<endl;      
      // cout<<TString::Format( " ---> 68 C.L. CNP      %8.3f %8.3f, width %8.3f", CL68_low_CNP, CL68_hgh_CNP, width_CL68_CNP)<<endl;      
      // cout<<TString::Format( " ---> 68 C.L. Gauss    %8.3f %8.3f, width %8.3f", CL68_low_Gauss, CL68_hgh_Gauss, width_CL68_Gauss)<<endl;      
      // cout<<endl;
     
      ///////////////////////////////////////////////////////////////////////// CL90
      
      ///////////////////////////////// Poisson
      set<double>cl_sigma90_Poisson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Poisson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Poisson;
	if( deltaChi2<fc_sigma90_Poisson[idx] ) cl_sigma90_Poisson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma90_Poisson;
      for( it_cl_sigma90_Poisson=cl_sigma90_Poisson.begin(); it_cl_sigma90_Poisson!=cl_sigma90_Poisson.end(); it_cl_sigma90_Poisson++ ) {
	double result = *it_cl_sigma90_Poisson;
      }
      CL90_low_Poisson = *(cl_sigma90_Poisson.begin());
      CL90_hgh_Poisson = *(cl_sigma90_Poisson.rbegin());
      double width_CL90_Poisson = CL90_hgh_Poisson - CL90_low_Poisson;

      ///////////////////////////////// Neyman
      set<double>cl_sigma90_Neyman;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Neyman( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Neyman;
	if( deltaChi2<fc_sigma90_Neyman[idx] ) cl_sigma90_Neyman.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma90_Neyman;
      for( it_cl_sigma90_Neyman=cl_sigma90_Neyman.begin(); it_cl_sigma90_Neyman!=cl_sigma90_Neyman.end(); it_cl_sigma90_Neyman++ ) {
	double result = *it_cl_sigma90_Neyman;
      }
      CL90_low_Neyman = *(cl_sigma90_Neyman.begin());
      CL90_hgh_Neyman = *(cl_sigma90_Neyman.rbegin());
      double width_CL90_Neyman = CL90_hgh_Neyman - CL90_low_Neyman;

      ///////////////////////////////// Pearson
      set<double>cl_sigma90_Pearson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Pearson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Pearson;
	if( deltaChi2<fc_sigma90_Pearson[idx] ) cl_sigma90_Pearson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma90_Pearson;
      for( it_cl_sigma90_Pearson=cl_sigma90_Pearson.begin(); it_cl_sigma90_Pearson!=cl_sigma90_Pearson.end(); it_cl_sigma90_Pearson++ ) {
	double result = *it_cl_sigma90_Pearson;
      }
      CL90_low_Pearson = *(cl_sigma90_Pearson.begin());
      CL90_hgh_Pearson = *(cl_sigma90_Pearson.rbegin());
      double width_CL90_Pearson = CL90_hgh_Pearson - CL90_low_Pearson;

      ///////////////////////////////// CNP
      set<double>cl_sigma90_CNP;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_CNP( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_CNP;
	if( deltaChi2<fc_sigma90_CNP[idx] ) cl_sigma90_CNP.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma90_CNP;
      for( it_cl_sigma90_CNP=cl_sigma90_CNP.begin(); it_cl_sigma90_CNP!=cl_sigma90_CNP.end(); it_cl_sigma90_CNP++ ) {
	double result = *it_cl_sigma90_CNP;
      }
      CL90_low_CNP = *(cl_sigma90_CNP.begin());
      CL90_hgh_CNP = *(cl_sigma90_CNP.rbegin());
      double width_CL90_CNP = CL90_hgh_CNP - CL90_low_CNP;

      ///////////////////////////////// Gauss
      set<double>cl_sigma90_Gauss;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Gauss( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Gauss;
	if( deltaChi2<fc_sigma90_Gauss[idx] ) cl_sigma90_Gauss.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma90_Gauss;
      for( it_cl_sigma90_Gauss=cl_sigma90_Gauss.begin(); it_cl_sigma90_Gauss!=cl_sigma90_Gauss.end(); it_cl_sigma90_Gauss++ ) {
	double result = *it_cl_sigma90_Gauss;
      }
      CL90_low_Gauss = *(cl_sigma90_Gauss.begin());
      CL90_hgh_Gauss = *(cl_sigma90_Gauss.rbegin());
      double width_CL90_Gauss = CL90_hgh_Gauss - CL90_low_Gauss;

      // cout<<endl;
      // cout<<TString::Format( " ---> 90 C.L. Poisson  %8.3f %8.3f, width %8.3f", CL90_low_Poisson, CL90_hgh_Poisson, width_CL90_Poisson)<<endl;      
      // cout<<TString::Format( " ---> 90 C.L. Neyman   %8.3f %8.3f, width %8.3f", CL90_low_Neyman, CL90_hgh_Neyman, width_CL90_Neyman)<<endl;      
      // cout<<TString::Format( " ---> 90 C.L. Pearson  %8.3f %8.3f, width %8.3f", CL90_low_Pearson, CL90_hgh_Pearson, width_CL90_Pearson)<<endl;      
      // cout<<TString::Format( " ---> 90 C.L. CNP      %8.3f %8.3f, width %8.3f", CL90_low_CNP, CL90_hgh_CNP, width_CL90_CNP)<<endl;      
      // cout<<TString::Format( " ---> 90 C.L. Gauss    %8.3f %8.3f, width %8.3f", CL90_low_Gauss, CL90_hgh_Gauss, width_CL90_Gauss)<<endl;      
      // cout<<endl;
     
      ///////////////////////////////////////////////////////////////////////// CL95
      
      ///////////////////////////////// Poisson
      set<double>cl_sigma95_Poisson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Poisson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Poisson;
	if( deltaChi2<fc_sigma95_Poisson[idx] ) cl_sigma95_Poisson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma95_Poisson;
      for( it_cl_sigma95_Poisson=cl_sigma95_Poisson.begin(); it_cl_sigma95_Poisson!=cl_sigma95_Poisson.end(); it_cl_sigma95_Poisson++ ) {
	double result = *it_cl_sigma95_Poisson;
      }
      CL95_low_Poisson = *(cl_sigma95_Poisson.begin());
      CL95_hgh_Poisson = *(cl_sigma95_Poisson.rbegin());
      double width_CL95_Poisson = CL95_hgh_Poisson - CL95_low_Poisson;

      ///////////////////////////////// Neyman
      set<double>cl_sigma95_Neyman;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Neyman( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Neyman;
	if( deltaChi2<fc_sigma95_Neyman[idx] ) cl_sigma95_Neyman.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma95_Neyman;
      for( it_cl_sigma95_Neyman=cl_sigma95_Neyman.begin(); it_cl_sigma95_Neyman!=cl_sigma95_Neyman.end(); it_cl_sigma95_Neyman++ ) {
	double result = *it_cl_sigma95_Neyman;
      }
      CL95_low_Neyman = *(cl_sigma95_Neyman.begin());
      CL95_hgh_Neyman = *(cl_sigma95_Neyman.rbegin());
      double width_CL95_Neyman = CL95_hgh_Neyman - CL95_low_Neyman;

      ///////////////////////////////// Pearson
      set<double>cl_sigma95_Pearson;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Pearson( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Pearson;
	if( deltaChi2<fc_sigma95_Pearson[idx] ) cl_sigma95_Pearson.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma95_Pearson;
      for( it_cl_sigma95_Pearson=cl_sigma95_Pearson.begin(); it_cl_sigma95_Pearson!=cl_sigma95_Pearson.end(); it_cl_sigma95_Pearson++ ) {
	double result = *it_cl_sigma95_Pearson;
      }
      CL95_low_Pearson = *(cl_sigma95_Pearson.begin());
      CL95_hgh_Pearson = *(cl_sigma95_Pearson.rbegin());
      double width_CL95_Pearson = CL95_hgh_Pearson - CL95_low_Pearson;

      ///////////////////////////////// CNP
      set<double>cl_sigma95_CNP;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_CNP( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_CNP;
	if( deltaChi2<fc_sigma95_CNP[idx] ) cl_sigma95_CNP.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma95_CNP;
      for( it_cl_sigma95_CNP=cl_sigma95_CNP.begin(); it_cl_sigma95_CNP!=cl_sigma95_CNP.end(); it_cl_sigma95_CNP++ ) {
	double result = *it_cl_sigma95_CNP;
      }
      CL95_low_CNP = *(cl_sigma95_CNP.begin());
      CL95_hgh_CNP = *(cl_sigma95_CNP.rbegin());
      double width_CL95_CNP = CL95_hgh_CNP - CL95_low_CNP;

      ///////////////////////////////// Gauss
      set<double>cl_sigma95_Gauss;
      for(int idx=0; idx<entries_fc; idx++) {
	double chi2true =  chi2_Gauss( fc_mu_true[idx], array_val_m, Ndet );
	double deltaChi2 = chi2true - chi2bestFit_Gauss;
	if( deltaChi2<fc_sigma95_Gauss[idx] ) cl_sigma95_Gauss.insert( fc_mu_true[idx] );
      }
      set<double>::iterator it_cl_sigma95_Gauss;
      for( it_cl_sigma95_Gauss=cl_sigma95_Gauss.begin(); it_cl_sigma95_Gauss!=cl_sigma95_Gauss.end(); it_cl_sigma95_Gauss++ ) {
	double result = *it_cl_sigma95_Gauss;
      }
      CL95_low_Gauss = *(cl_sigma95_Gauss.begin());
      CL95_hgh_Gauss = *(cl_sigma95_Gauss.rbegin());
      double width_CL95_Gauss = CL95_hgh_Gauss - CL95_low_Gauss;

      // cout<<endl;
      // cout<<TString::Format( " ---> 95 C.L. Poisson  %8.3f %8.3f, width %8.3f", CL95_low_Poisson, CL95_hgh_Poisson, width_CL95_Poisson)<<endl;      
      // cout<<TString::Format( " ---> 95 C.L. Neyman   %8.3f %8.3f, width %8.3f", CL95_low_Neyman, CL95_hgh_Neyman, width_CL95_Neyman)<<endl;      
      // cout<<TString::Format( " ---> 95 C.L. Pearson  %8.3f %8.3f, width %8.3f", CL95_low_Pearson, CL95_hgh_Pearson, width_CL95_Pearson)<<endl;      
      // cout<<TString::Format( " ---> 95 C.L. CNP      %8.3f %8.3f, width %8.3f", CL95_low_CNP, CL95_hgh_CNP, width_CL95_CNP)<<endl;      
      // cout<<TString::Format( " ---> 95 C.L. Gauss    %8.3f %8.3f, width %8.3f", CL95_low_Gauss, CL95_hgh_Gauss, width_CL95_Gauss)<<endl;      
      // cout<<endl;



      

      
      tree_chi2->Fill();
      
    }// ievent    
  }// irun
  cout<<endl;

  
  cout<<TString::Format(" mu_bestFit out of the scaned range: %d", line_outScaned_range)<<endl<<endl;

  
  tree_chi2->Write();
  out_file->Close();

}

