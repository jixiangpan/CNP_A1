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
#include "TGraphErrors.h"
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

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"


//
// y = p0 + p1*x
//
// x = 2,3,4,5,6,7,8,9,10
// Poisson(y)
//

////////////////////////////////////////////// Global

map<int, double>array_val_meas;
const int xmin = 1;
const int xmax = 100;

double func_linear(double x, double p0, double p1)
{
  return p0 + p1*x;
}

////////////////////////////////////////////// FCN

double fcn_Poisson(const double *par)
{
  double chi2 = 0;
  
  double p0 = par[0];
  double p1 = par[1];

  for(int idx=xmin; idx<=xmax; idx++) {
    double xcur = idx*0.1;
    double val_pred = func_linear(xcur, p0, p1);
    double val_meas = array_val_meas[idx];

    if( val_meas==0 ) chi2 += 2*val_pred;
    else chi2 += 2*( val_pred - val_meas + val_meas*log(val_meas/val_pred) );
  }

  return chi2;
}

double fcn_Neyman(const double *par)
{
  double chi2 = 0;
  
  double p0 = par[0];
  double p1 = par[1];

  for(int idx=xmin; idx<=xmax; idx++) {
    double xcur = idx*0.1;
    double val_pred = func_linear(xcur, p0, p1);
    double val_meas = array_val_meas[idx];

    if( val_meas==0 ) chi2 += 2*val_pred;
    else chi2 += pow( val_pred-val_meas,2 )/val_meas;
  }

  return chi2;
}

double fcn_Pearson(const double *par)
{
  double chi2 = 0;
  
  double p0 = par[0];
  double p1 = par[1];

  for(int idx=xmin; idx<=xmax; idx++) {
    double xcur = idx*0.1;
    double val_pred = func_linear(xcur, p0, p1);
    double val_meas = array_val_meas[idx];

    chi2 += pow( val_pred-val_meas,2 )/val_pred;
  }

  return chi2;
}

double fcn_CNP(const double *par)
{
  double chi2 = 0;
  
  double p0 = par[0];
  double p1 = par[1];

  for(int idx=xmin; idx<=xmax; idx++) {
    double xcur = idx*0.1;
    double val_pred = func_linear(xcur, p0, p1);
    double val_meas = array_val_meas[idx];

    if( val_meas==0 ) chi2 += 2*val_pred;
    else chi2 += pow( val_pred-val_meas, 2 )/3 * ( 1./val_meas + 2./val_pred );
  }

  return chi2;
}

double fcn_Gauss(const double *par)
{
  double chi2 = 0;
  
  double p0 = par[0];
  double p1 = par[1];

  for(int idx=xmin; idx<=xmax; idx++) {
    double xcur = idx*0.1;
    double val_pred = func_linear(xcur, p0, p1);
    double val_meas = array_val_meas[idx];

    double m_saturated = sqrt(4*val_meas*val_meas + 1)/2 - 0.5;
    
    if( val_meas==0 ) chi2 += 2*val_pred;
    else chi2 += ( pow(val_meas-val_pred,2)/val_pred -pow(val_meas-m_saturated,2)/m_saturated + log(val_pred/m_saturated) );
  }

  return chi2;
}



////////////////////////////////////////////// MAIN

void read_linearFit(double p0_true, double p1_true, int Nrun, int Nevent)
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  TString roostr = "";

  const double par_true[] = {p0_true, p1_true};
  
  //////////////////////////////////////// Tree

  double p0_bestFit_Poisson = 0;
  double p0_bestFit_Neyman  = 0;
  double p0_bestFit_Pearson = 0;
  double p0_bestFit_CNP     = 0;
  double p0_bestFit_Gauss   = 0;

  double p1_bestFit_Poisson = 0;
  double p1_bestFit_Neyman  = 0;
  double p1_bestFit_Pearson = 0;
  double p1_bestFit_CNP     = 0;
  double p1_bestFit_Gauss   = 0;

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
  
  tree_chi2->Branch( "p0_bestFit_Poisson", &p0_bestFit_Poisson, "p0_bestFit_Poisson/D" );
  tree_chi2->Branch( "p0_bestFit_Neyman",  &p0_bestFit_Neyman,  "p0_bestFit_Neyman/D" );
  tree_chi2->Branch( "p0_bestFit_Pearson", &p0_bestFit_Pearson, "p0_bestFit_Pearson/D" );
  tree_chi2->Branch( "p0_bestFit_CNP",     &p0_bestFit_CNP,     "p0_bestFit_CNP/D" );
  tree_chi2->Branch( "p0_bestFit_Gauss",   &p0_bestFit_Gauss,   "p0_bestFit_Gauss/D" );
  
  tree_chi2->Branch( "p1_bestFit_Poisson", &p1_bestFit_Poisson, "p1_bestFit_Poisson/D" );
  tree_chi2->Branch( "p1_bestFit_Neyman",  &p1_bestFit_Neyman,  "p1_bestFit_Neyman/D" );
  tree_chi2->Branch( "p1_bestFit_Pearson", &p1_bestFit_Pearson, "p1_bestFit_Pearson/D" );
  tree_chi2->Branch( "p1_bestFit_CNP",     &p1_bestFit_CNP,     "p1_bestFit_CNP/D" );
  tree_chi2->Branch( "p1_bestFit_Gauss",   &p1_bestFit_Gauss,   "p1_bestFit_Gauss/D" );
  
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
    
  for(int irun=1; irun<=Nrun; irun++) {

    if(irun%100==0) cout<<TString::Format(" -------------------> processing %6.4f", irun*1./Nrun)<<endl;
    
    TRandom3 *roo_random = new TRandom3(0);

    for(int ievent=1; ievent<=Nevent; ievent++) {
      array_val_meas.clear();

      // TGraphErrors *gh_test = new TGraphErrors(); gh_test->SetName("gh_test");
      // int line_gh_test = 0;
      
      for(int idx=xmin; idx<=xmax; idx++) {
	double xcur = idx*0.1;
	double meas_true = func_linear(xcur, p0_true, p1_true);
	double meas_rand = roo_random->Poisson( meas_true ) *1.;
	//meas_rand = meas_true;
	array_val_meas[idx] = meas_rand;

	// line_gh_test++;
	// gh_test->SetPoint(line_gh_test-1, xcur, meas_rand);
	// gh_test->SetPointError(line_gh_test-1, 0, sqrt(meas_rand) );
	
      }

      // gh_test->Draw("apE");
      // gh_test->GetXaxis()->SetTitle("x");
      // gh_test->GetYaxis()->SetTitle("n");
      

      ////////////////////////////////////////////////////////////// Poisson
            
      ROOT::Minuit2::Minuit2Minimizer min_Poisson( ROOT::Minuit2::kMigrad );
      min_Poisson.SetPrintLevel(0);
      min_Poisson.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
      min_Poisson.SetMaxFunctionCalls(500000);
      min_Poisson.SetMaxIterations(500000);
      min_Poisson.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
      min_Poisson.SetPrecision(1e-18); //precision in the target function
  
      /// set fitting parameters
      ROOT::Math::Functor Chi2Functor_Poisson( &fcn_Poisson, 2 );
      min_Poisson.SetFunction(Chi2Functor_Poisson);

      min_Poisson.SetVariable( 0, "par_p0", p0_true, 1e-2);
      min_Poisson.SetVariable( 1, "par_p1", p1_true, 1e-2);
      //min_Poisson.SetFixedVariable( 1, "par_p1", 0);
      
      /// do the minimization
      min_Poisson.Minimize();
      int status_Poisson = min_Poisson.Status();
      const double *par_Poisson = min_Poisson.X();

      p0_bestFit_Poisson = par_Poisson[0];
      p1_bestFit_Poisson = par_Poisson[1];
      chi2bestFit_Poisson = min_Poisson.MinValue();
      chi2true_Poisson = fcn_Poisson(par_true);
      
      /// https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
      // min_Poisson.SetErrorDef( 2.2977 );// 1 - 0.317 with ndf = 2
      // //min.SetErrorDef( 6.2021 );// 1 - 0.045 with ndf = 2
      // //min.SetErrorDef(11.6182 );// 1 - 0.003 with ndf = 2      
      // double array_p0_Poisson[80] = {0};
      // double array_p1_Poisson[80] = {0};
      // unsigned int idx_p0_Poisson = 0;
      // unsigned int idx_p1_Poisson = 1;
      // unsigned int npoints_Poisson = 80;
      // min_Poisson.Contour(idx_p0_Poisson, idx_p1_Poisson, npoints_Poisson, array_p0_Poisson, array_p1_Poisson);
      // TGraph *gh_contour_Poisson = new TGraph(); gh_contour_Poisson->SetName("gh_contour_Poisson");
      // for(int i=1; i<=80; i++) {    
      // 	// cout<<TString::Format(" ---> %3d %8.4f %8.4f",
      // 	// 		      i, array_p0_Poisson[i-1], array_p1_Poisson[i-1])<<endl;
      // 	gh_contour_Poisson->SetPoint(i-1, array_p0_Poisson[i-1], array_p1_Poisson[i-1]);
      // }    
      // gh_contour_Poisson->Draw("apl");

      ////////////////////////////////////////////////////////////// Neyman
            
      ROOT::Minuit2::Minuit2Minimizer min_Neyman( ROOT::Minuit2::kMigrad );
      min_Neyman.SetPrintLevel(0);
      min_Neyman.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
      min_Neyman.SetMaxFunctionCalls(500000);
      min_Neyman.SetMaxIterations(500000);
      min_Neyman.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
      min_Neyman.SetPrecision(1e-18); //precision in the target function
  
      /// set fitting parameters
      ROOT::Math::Functor Chi2Functor_Neyman( &fcn_Neyman, 2 );
      min_Neyman.SetFunction(Chi2Functor_Neyman);

      min_Neyman.SetVariable( 0, "par_p0", p0_true, 1e-2);
      min_Neyman.SetVariable( 1, "par_p1", p1_true, 1e-2);
      //min_Neyman.SetFixedVariable( 1, "par_p1", 0);
      
      /// do the minimization
      min_Neyman.Minimize();
      int status_Neyman = min_Neyman.Status();
      const double *par_Neyman = min_Neyman.X();

      p0_bestFit_Neyman = par_Neyman[0];
      p1_bestFit_Neyman = par_Neyman[1];
      chi2bestFit_Neyman = min_Neyman.MinValue();
      chi2true_Neyman = fcn_Neyman(par_true);
      
      ////////////////////////////////////////////////////////////// Pearson
            
      ROOT::Minuit2::Minuit2Minimizer min_Pearson( ROOT::Minuit2::kMigrad );
      min_Pearson.SetPrintLevel(0);
      min_Pearson.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
      min_Pearson.SetMaxFunctionCalls(500000);
      min_Pearson.SetMaxIterations(500000);
      min_Pearson.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
      min_Pearson.SetPrecision(1e-18); //precision in the target function
  
      /// set fitting parameters
      ROOT::Math::Functor Chi2Functor_Pearson( &fcn_Pearson, 2 );
      min_Pearson.SetFunction(Chi2Functor_Pearson);

      min_Pearson.SetVariable( 0, "par_p0", p0_true, 1e-2);
      min_Pearson.SetVariable( 1, "par_p1", p1_true, 1e-2);
      //min_Pearson.SetFixedVariable( 1, "par_p1", 0);
      
      /// do the minimization
      min_Pearson.Minimize();
      int status_Pearson = min_Pearson.Status();
      const double *par_Pearson = min_Pearson.X();

      p0_bestFit_Pearson = par_Pearson[0];
      p1_bestFit_Pearson = par_Pearson[1];
      chi2bestFit_Pearson = min_Pearson.MinValue();
      chi2true_Pearson = fcn_Pearson(par_true);
      
      ////////////////////////////////////////////////////////////// CNP
            
      ROOT::Minuit2::Minuit2Minimizer min_CNP( ROOT::Minuit2::kMigrad );
      min_CNP.SetPrintLevel(0);
      min_CNP.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
      min_CNP.SetMaxFunctionCalls(500000);
      min_CNP.SetMaxIterations(500000);
      min_CNP.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
      min_CNP.SetPrecision(1e-18); //precision in the target function
  
      /// set fitting parameters
      ROOT::Math::Functor Chi2Functor_CNP( &fcn_CNP, 2 );
      min_CNP.SetFunction(Chi2Functor_CNP);

      min_CNP.SetVariable( 0, "par_p0", p0_true, 1e-2);
      min_CNP.SetVariable( 1, "par_p1", p1_true, 1e-2);
      //min_CNP.SetFixedVariable( 1, "par_p1", 0);
      
      /// do the minimization
      min_CNP.Minimize();
      int status_CNP = min_CNP.Status();
      const double *par_CNP = min_CNP.X();

      p0_bestFit_CNP = par_CNP[0];
      p1_bestFit_CNP = par_CNP[1];
      chi2bestFit_CNP = min_CNP.MinValue();
      chi2true_CNP = fcn_CNP(par_true);
      
      ////////////////////////////////////////////////////////////// Gauss
            
      ROOT::Minuit2::Minuit2Minimizer min_Gauss( ROOT::Minuit2::kMigrad );
      min_Gauss.SetPrintLevel(0);
      min_Gauss.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
      min_Gauss.SetMaxFunctionCalls(500000);
      min_Gauss.SetMaxIterations(500000);
      min_Gauss.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
      min_Gauss.SetPrecision(1e-18); //precision in the target function
  
      /// set fitting parameters
      ROOT::Math::Functor Chi2Functor_Gauss( &fcn_Gauss, 2 );
      min_Gauss.SetFunction(Chi2Functor_Gauss);

      min_Gauss.SetVariable( 0, "par_p0", p0_true, 1e-2);
      min_Gauss.SetVariable( 1, "par_p1", p1_true, 1e-2);
      //min_Gauss.SetFixedVariable( 1, "par_p1", 0);
      
      /// do the minimization
      min_Gauss.Minimize();
      int status_Gauss = min_Gauss.Status();
      const double *par_Gauss = min_Gauss.X();

      p0_bestFit_Gauss = par_Gauss[0];
      p1_bestFit_Gauss = par_Gauss[1];
      chi2bestFit_Gauss = min_Gauss.MinValue();
      chi2true_Gauss = fcn_Gauss(par_true);
      
      ////////////////////////////////////////////////////////////// Fill Tree

      if( status_Poisson+status_Neyman+status_Pearson+status_CNP+status_Gauss==0 ) {	
	tree_chi2->Fill();
      }
      else {
	int check = 0;
	// cout<<" ---> Bad fitting"<<endl;
	// cout<<" ---> Poisson "<<status_Poisson<<endl;
	// cout<<" ---> Neyman  "<<status_Neyman<<endl;
	// cout<<" ---> Pearson "<<status_Pearson<<endl;
	// cout<<" ---> CNP     "<<status_CNP<<endl;
	// cout<<" ---> Gauss   "<<status_Gauss<<endl;
	// cout<<endl;
      }
      
    }// ievent

  }// irun

  roostr = TString::Format("out_linear_p0_%03d_p1_%03d_Nrun%06d_Nevent%06d.root", (int)p0_true, (int)p1_true, Nrun, Nevent);
  TFile *outfile = new TFile(roostr, "recreate");
  tree_chi2->Write();
  outfile->Close();
  
}

