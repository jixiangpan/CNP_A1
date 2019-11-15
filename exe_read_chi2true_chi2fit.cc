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




//void read_read_chi2true_chi2fit(double Nmu, int Ndet, int Nrun, int Nevent)
int main(int argc, char** argv)
{
  //gROOT->ProcessLine(".x ./DrawOption.cc");
  TString roostr = "";
  
  double Nmu    = 0;
  int    Ndet   = 0;
  int    Nrun   = 0;
  int    Nevent = 0;
 
  for(int i=1; i<argc; i++) {
    
    if( strcmp(argv[i],"-Nmu")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>Nmu ) ) { cerr<<" ---> Error Nmu !"<<endl; exit(1); }
    }
    
    if( strcmp(argv[i],"-Ndet")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>Ndet ) ) { cerr<<" ---> Error Ndet !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-Nrun")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>Nrun ) ) { cerr<<" ---> Error Nrun !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-Nevent")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>Nevent ) ) { cerr<<" ---> Error Nevent !"<<endl; exit(1); }
    }
    
  }

  if( Nmu==0 || Ndet==0 || Nrun==0 || Nevent==0 ) {
    cerr<<endl<<" ---> ERROR , 0 in inputs "<<endl<<endl;
    exit(1);
  }

  cout<<endl;
  cout<<TString::Format(" ---> Nmu    %8.3f", Nmu)<<endl;
  cout<<TString::Format(" ---> Ndet   %8d", Ndet)<<endl;
  cout<<TString::Format(" ---> Nrun   %8d", Nrun)<<endl;
  cout<<TString::Format(" ---> Nevent %8d", Nevent)<<endl;
  cout<<endl;

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
  
  //////////////////////////////////////// Toy
  
  map<int, double>array_val_m;
  
  for(int irun=1; irun<=Nrun; irun++) {

    if(irun%100==0) cout<<TString::Format(" ---> processing %6.4f", irun*1./Nrun)<<endl;
    
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
      
      tree_chi2->Fill();
      
    }// ievent    
  }// irun
  cout<<endl;

  roostr = TString::Format("out_chi2_mu%08.3f_det%03d_run%06d_event%06d.root", Nmu, Ndet, Nrun, Nevent);
  TFile *out_file = new TFile(roostr, "recreate");
  tree_chi2->Write();
  out_file->Close();

  
  return 0;
}

