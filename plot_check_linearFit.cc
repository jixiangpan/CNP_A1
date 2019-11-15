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

void plot_check_linearFit()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");
  
  TString roostr = "";
  
  int color_Poisson = kBlack;
  int color_Neyman  = kMagenta;
  int color_Pearson = kBlue;
  int color_CNP     = kGreen;
  int color_Gauss   = kOrange;
  
  /////////////////////////////////////////////////////////////////////////////

  double p0_true = 5;
  double p1_true = 1;
  
  //roostr = "./data_02_check_linearFit/out_linear_p0_008_p1_020_Nrun001000_Nevent001000.root";
  //roostr = "./data_02_check_linearFit/out_linear_p0_004_p1_020_Nrun000100_Nevent001000.root";
  //roostr = "./data_02_check_linearFit/out_linear_p0_040_p1_100_Nrun001000_Nevent001000.root";

  //roostr = "./data_02_check_linearFit/out_linear_p0_015_p1_000_Nrun001000_Nevent001000.root";
  //roostr = "./data_02_check_linearFit/out_linear_p0_015_p1_000_Nrun001000_Nevent001000_Ndet35.root";

  //roostr = "./data_02_check_linearFit/out_linear_p0_008_p1_020_Nrun001000_Nevent001000_Ndet35.root";

  //roostr = "./data_02_check_linearFit/out_linear_p0_015_p1_000_Nrun001000_Nevent001000.root";

  roostr = "./data_02_check_linearFit/out_linear_p0_005_p1_001_Nrun001000_Nevent001000_Ndet100.root";
  
  
  TFile *file_input = new TFile(roostr, "read");
  TTree *tree_chi2 = (TTree*)file_input->Get("tree_chi2");

  // Declaration of leaf types
  Double_t        p0_bestFit_Poisson;
  Double_t        p0_bestFit_Neyman;
  Double_t        p0_bestFit_Pearson;
  Double_t        p0_bestFit_CNP;
  Double_t        p0_bestFit_Gauss;
  Double_t        p1_bestFit_Poisson;
  Double_t        p1_bestFit_Neyman;
  Double_t        p1_bestFit_Pearson;
  Double_t        p1_bestFit_CNP;
  Double_t        p1_bestFit_Gauss;  
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
  TBranch        *b_p0_bestFit_Poisson;   //!
  TBranch        *b_p0_bestFit_Neyman;   //!
  TBranch        *b_p0_bestFit_Pearson;   //!
  TBranch        *b_p0_bestFit_CNP;   //!
  TBranch        *b_p0_bestFit_Gauss;   //!
  TBranch        *b_p1_bestFit_Poisson;   //!
  TBranch        *b_p1_bestFit_Neyman;   //!
  TBranch        *b_p1_bestFit_Pearson;   //!
  TBranch        *b_p1_bestFit_CNP;   //!
  TBranch        *b_p1_bestFit_Gauss;   //!  
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
  tree_chi2->SetBranchAddress("p0_bestFit_Poisson", &p0_bestFit_Poisson, &b_p0_bestFit_Poisson);
  tree_chi2->SetBranchAddress("p0_bestFit_Neyman", &p0_bestFit_Neyman, &b_p0_bestFit_Neyman);
  tree_chi2->SetBranchAddress("p0_bestFit_Pearson", &p0_bestFit_Pearson, &b_p0_bestFit_Pearson);
  tree_chi2->SetBranchAddress("p0_bestFit_CNP", &p0_bestFit_CNP, &b_p0_bestFit_CNP);
  tree_chi2->SetBranchAddress("p0_bestFit_Gauss", &p0_bestFit_Gauss, &b_p0_bestFit_Gauss);
  tree_chi2->SetBranchAddress("p1_bestFit_Poisson", &p1_bestFit_Poisson, &b_p1_bestFit_Poisson);
  tree_chi2->SetBranchAddress("p1_bestFit_Neyman", &p1_bestFit_Neyman, &b_p1_bestFit_Neyman);
  tree_chi2->SetBranchAddress("p1_bestFit_Pearson", &p1_bestFit_Pearson, &b_p1_bestFit_Pearson);
  tree_chi2->SetBranchAddress("p1_bestFit_CNP", &p1_bestFit_CNP, &b_p1_bestFit_CNP);
  tree_chi2->SetBranchAddress("p1_bestFit_Gauss", &p1_bestFit_Gauss, &b_p1_bestFit_Gauss);  
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

  ////////////////////////// Goodness of fit
  int nbins_chi2 = 200;
  double chi2_low = 0;
  double chi2_hgh = 40;
  
  roostr = "h1_chi2true_Poisson"; TH1D *h1_chi2true_Poisson = new TH1D(roostr, roostr, nbins_chi2, chi2_low, chi2_hgh);
  roostr = "h1_chi2true_Neyman"; TH1D *h1_chi2true_Neyman = new TH1D(roostr, roostr, nbins_chi2, chi2_low, chi2_hgh);
  roostr = "h1_chi2true_Pearson"; TH1D *h1_chi2true_Pearson = new TH1D(roostr, roostr, nbins_chi2, chi2_low, chi2_hgh);
  roostr = "h1_chi2true_CNP"; TH1D *h1_chi2true_CNP = new TH1D(roostr, roostr, nbins_chi2, chi2_low, chi2_hgh);
  roostr = "h1_chi2true_Gauss"; TH1D *h1_chi2true_Gauss = new TH1D(roostr, roostr, nbins_chi2, chi2_low, chi2_hgh);
  
  ////////////////////////// deltaChi2
  int nbins_deltaChi2 = 149;
  double deltaChi2_low = 0;
  double deltaChi2_hgh = 8;
  
  roostr = "h1_deltaChi2_Poisson"; TH1D *h1_deltaChi2_Poisson = new TH1D(roostr, roostr, nbins_deltaChi2, deltaChi2_low, deltaChi2_hgh);
  roostr = "h1_deltaChi2_Neyman"; TH1D *h1_deltaChi2_Neyman = new TH1D(roostr, roostr, nbins_deltaChi2, deltaChi2_low, deltaChi2_hgh);
  roostr = "h1_deltaChi2_Pearson"; TH1D *h1_deltaChi2_Pearson = new TH1D(roostr, roostr, nbins_deltaChi2, deltaChi2_low, deltaChi2_hgh);
  roostr = "h1_deltaChi2_CNP"; TH1D *h1_deltaChi2_CNP = new TH1D(roostr, roostr, nbins_deltaChi2, deltaChi2_low, deltaChi2_hgh);
  roostr = "h1_deltaChi2_Gauss"; TH1D *h1_deltaChi2_Gauss = new TH1D(roostr, roostr, nbins_deltaChi2, deltaChi2_low, deltaChi2_hgh);

  ////////////////////////// pars
  int nbins_p0 = 151;
  double p0_low = 0;
  double p0_hgh = 15;
  // int nbins_p0 = 100;
  // double p0_low = 10;
  // double p0_hgh = 100;
  // int nbins_p0 = 100;
  // double p0_low = 5;
  // double p0_hgh = 35;

  
  roostr = "h1_p0_Poisson"; TH1D *h1_p0_Poisson = new TH1D(roostr, roostr, nbins_p0, p0_low, p0_hgh);
  roostr = "h1_p0_Neyman"; TH1D *h1_p0_Neyman = new TH1D(roostr, roostr, nbins_p0, p0_low, p0_hgh);
  roostr = "h1_p0_Pearson"; TH1D *h1_p0_Pearson = new TH1D(roostr, roostr, nbins_p0, p0_low, p0_hgh);
  roostr = "h1_p0_CNP"; TH1D *h1_p0_CNP = new TH1D(roostr, roostr, nbins_p0, p0_low, p0_hgh);
  roostr = "h1_p0_Gauss"; TH1D *h1_p0_Gauss = new TH1D(roostr, roostr, nbins_p0, p0_low, p0_hgh);
  double bias_p0_Poisson = 0;
  double bias_p0_Neyman = 0;
  double bias_p0_Pearson = 0;
  double bias_p0_CNP = 0;
  double bias_p0_Gauss = 0;
  
  int nbins_p1 = 111;
  double p1_low = 0;
  double p1_hgh = 3;
  // int nbins_p1 = 100;
  // double p1_low = 50;
  // double p1_hgh = 200; 
  roostr = "h1_p1_Poisson"; TH1D *h1_p1_Poisson = new TH1D(roostr, roostr, nbins_p1, p1_low, p1_hgh);
  roostr = "h1_p1_Neyman"; TH1D *h1_p1_Neyman = new TH1D(roostr, roostr, nbins_p1, p1_low, p1_hgh);
  roostr = "h1_p1_Pearson"; TH1D *h1_p1_Pearson = new TH1D(roostr, roostr, nbins_p1, p1_low, p1_hgh);
  roostr = "h1_p1_CNP"; TH1D *h1_p1_CNP = new TH1D(roostr, roostr, nbins_p1, p1_low, p1_hgh);
  roostr = "h1_p1_Gauss"; TH1D *h1_p1_Gauss = new TH1D(roostr, roostr, nbins_p1, p1_low, p1_hgh);
  double bias_p1_Poisson = 0;
  double bias_p1_Neyman = 0;
  double bias_p1_Pearson = 0;
  double bias_p1_CNP = 0;
  double bias_p1_Gauss = 0;
  
  // int nbins_2d_p0 = 80;
  // double low_2d_p0 = -50;
  // double hgh_2d_p0 = 50;
  // int nbins_2d_p1 = 80;
  // double low_2d_p1 = -50;
  // double hgh_2d_p1 = 50;
  
  int nbins_2d_p0 = 80;
  double low_2d_p0 = -60;
  double hgh_2d_p0 = 60;
  int nbins_2d_p1 = 80;
  double low_2d_p1 = -60;
  double hgh_2d_p1 = 60;
  
  roostr = "h2_p1vsp0_Poisson"; TH2D *h2_p1vsp0_Poisson = new TH2D(roostr, roostr, nbins_2d_p0, low_2d_p0, hgh_2d_p0, nbins_2d_p1, low_2d_p1, hgh_2d_p1);
  roostr = "h2_p1vsp0_Neyman"; TH2D *h2_p1vsp0_Neyman = new TH2D(roostr, roostr, nbins_2d_p0, low_2d_p0, hgh_2d_p0, nbins_2d_p1, low_2d_p1, hgh_2d_p1);
  roostr = "h2_p1vsp0_Pearson"; TH2D *h2_p1vsp0_Pearson = new TH2D(roostr, roostr, nbins_2d_p0, low_2d_p0, hgh_2d_p0, nbins_2d_p1, low_2d_p1, hgh_2d_p1);
  roostr = "h2_p1vsp0_CNP"; TH2D *h2_p1vsp0_CNP = new TH2D(roostr, roostr, nbins_2d_p0, low_2d_p0, hgh_2d_p0, nbins_2d_p1, low_2d_p1, hgh_2d_p1);
  roostr = "h2_p1vsp0_Gauss"; TH2D *h2_p1vsp0_Gauss = new TH2D(roostr, roostr, nbins_2d_p0, low_2d_p0, hgh_2d_p0, nbins_2d_p1, low_2d_p1, hgh_2d_p1);
  
  ///////////////////////////////////////////////////////////////////
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_chi2->GetEntry( ientry );

    ////////////////////////// Goodness of fit
    h1_chi2true_Poisson->Fill( chi2true_Poisson );
    h1_chi2true_Neyman->Fill( chi2true_Neyman );
    h1_chi2true_Pearson->Fill( chi2true_Pearson );
    h1_chi2true_CNP->Fill( chi2true_CNP );
    h1_chi2true_Gauss->Fill( chi2true_Gauss );

    ////////////////////////// deltaChi2
    h1_deltaChi2_Poisson->Fill( chi2true_Poisson - chi2bestFit_Poisson );
    h1_deltaChi2_Neyman->Fill( chi2true_Neyman - chi2bestFit_Neyman );
    h1_deltaChi2_Pearson->Fill( chi2true_Pearson - chi2bestFit_Pearson );
    h1_deltaChi2_CNP->Fill( chi2true_CNP - chi2bestFit_CNP );
    h1_deltaChi2_Gauss->Fill( chi2true_Gauss - chi2bestFit_Gauss );

    ////////////////////////// pars
    h1_p0_Poisson->Fill( p0_bestFit_Poisson );
    h1_p0_Neyman->Fill( p0_bestFit_Neyman );
    h1_p0_Pearson->Fill( p0_bestFit_Pearson );
    h1_p0_CNP->Fill( p0_bestFit_CNP );
    h1_p0_Gauss->Fill( p0_bestFit_Gauss );
    bias_p0_Poisson += (p0_bestFit_Poisson-p0_true)/p0_true/entries;
    bias_p0_Neyman += (p0_bestFit_Neyman-p0_true)/p0_true/entries;
    bias_p0_Pearson += (p0_bestFit_Pearson-p0_true)/p0_true/entries;
    bias_p0_CNP += (p0_bestFit_CNP-p0_true)/p0_true/entries;
    bias_p0_Gauss += (p0_bestFit_Gauss-p0_true)/p0_true/entries;
    
    h1_p1_Poisson->Fill( p1_bestFit_Poisson );
    h1_p1_Neyman->Fill( p1_bestFit_Neyman );
    h1_p1_Pearson->Fill( p1_bestFit_Pearson );
    h1_p1_CNP->Fill( p1_bestFit_CNP );
    h1_p1_Gauss->Fill( p1_bestFit_Gauss );
    bias_p1_Poisson += (p1_bestFit_Poisson-p1_true)/p1_true/entries;
    bias_p1_Neyman += (p1_bestFit_Neyman-p1_true)/p1_true/entries;
    bias_p1_Pearson += (p1_bestFit_Pearson-p1_true)/p1_true/entries;
    bias_p1_CNP += (p1_bestFit_CNP-p1_true)/p1_true/entries;
    bias_p1_Gauss += (p1_bestFit_Gauss-p1_true)/p1_true/entries;

    h2_p1vsp0_Poisson->Fill( (p0_bestFit_Poisson-p0_true)*100/p0_true, (p1_bestFit_Poisson-p1_true)*100/p1_true );
    h2_p1vsp0_Neyman->Fill( (p0_bestFit_Neyman-p0_true)*100/p0_true, (p1_bestFit_Neyman-p1_true)*100/p1_true );
    h2_p1vsp0_Pearson->Fill( (p0_bestFit_Pearson-p0_true)*100/p0_true, (p1_bestFit_Pearson-p1_true)*100/p1_true );
    h2_p1vsp0_CNP->Fill( (p0_bestFit_CNP-p0_true)*100/p0_true, (p1_bestFit_CNP-p1_true)*100/p1_true );
    h2_p1vsp0_Gauss->Fill( (p0_bestFit_Gauss-p0_true)*100/p0_true, (p1_bestFit_Gauss-p1_true)*100/p1_true );
    
  }

  h1_chi2true_Poisson->Scale( 1./entries );
  h1_chi2true_Neyman->Scale( 1./entries );
  h1_chi2true_Pearson->Scale( 1./entries );
  h1_chi2true_CNP->Scale( 1./entries );
  h1_chi2true_Gauss->Scale( 1./entries );

  h1_deltaChi2_Poisson->Scale( 1./entries );
  h1_deltaChi2_Neyman->Scale( 1./entries );
  h1_deltaChi2_Pearson->Scale( 1./entries );
  h1_deltaChi2_CNP->Scale( 1./entries );
  h1_deltaChi2_Gauss->Scale( 1./entries );

  
  //////////////////////////////////////////////////////////////////////////////////////// h1_chi2true_Poisson
  
  roostr = "canv_h1_chi2true_Poisson";
  TCanvas *canv_h1_chi2true_Poisson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h1_chi2true_Poisson, 0.18, 0.1, 0.18, 0.06);  
  h1_chi2true_Poisson->Draw();  
  h1_chi2true_Poisson->SetLineColor(color_Poisson);
  roo_title(h1_chi2true_Poisson, "#chi^{2}_{true}", "Probability");
  roo_center_title(h1_chi2true_Poisson);
  roo_offset(h1_chi2true_Poisson, 1.1, 1.3);
  
  h1_chi2true_Neyman->Draw("same");  
  h1_chi2true_Neyman->SetLineColor(color_Neyman);
      
  h1_chi2true_Pearson->Draw("same");  
  h1_chi2true_Pearson->SetLineColor(color_Pearson);
      
  h1_chi2true_CNP->Draw("same");  
  h1_chi2true_CNP->SetLineColor(color_CNP);
      
  h1_chi2true_Gauss->Draw("same");  
  h1_chi2true_Gauss->SetLineColor(color_Gauss);

  int ndf_chi2true = 10;
  TF1 *f_chi2true = new TF1("f_chi2true",
			    TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", ndf_chi2true, (chi2_hgh-chi2_low)/nbins_chi2  ),
			    chi2_low, chi2_hgh);
  f_chi2true->Draw("same");
  f_chi2true->SetLineColor(kRed);
  
  TLegend *lg_chi2true = new TLegend(0.55+0.07,0.56-0.1,0.83+0.07,0.88);
  lg_chi2true->SetBorderSize(0);
  lg_chi2true->AddEntry(  f_chi2true, TString::Format("chi2, ndf=%d", ndf_chi2true), "l");
  lg_chi2true->AddEntry(  h1_chi2true_Poisson, TString::Format("Poisson"), "l" );
  lg_chi2true->AddEntry(  h1_chi2true_Neyman,  TString::Format("Neyman"), "l" );
  lg_chi2true->AddEntry(  h1_chi2true_Pearson, TString::Format("Pearson"), "l" );
  lg_chi2true->AddEntry(  h1_chi2true_CNP,     TString::Format("CNP"), "l" );
  lg_chi2true->AddEntry(  h1_chi2true_Gauss,   TString::Format("Gauss"), "l" );			    
  lg_chi2true->Draw("same");
  lg_chi2true->SetTextSize(0.055);
  lg_chi2true->SetTextFont(42);

  //////////////////////////////////////////////////////////////////////////////////////// h1_deltaChi2_Poisson
  
  roostr = "canv_h1_deltaChi2_Poisson";
  TCanvas *canv_h1_deltaChi2_Poisson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h1_deltaChi2_Poisson, 0.18, 0.1, 0.18, 0.06);  
  h1_deltaChi2_Poisson->Draw();  
  h1_deltaChi2_Poisson->SetLineColor(color_Poisson);
  roo_title(h1_deltaChi2_Poisson, "#chi^{2}_{true}", "Probability");
  roo_center_title(h1_deltaChi2_Poisson);
  roo_offset(h1_deltaChi2_Poisson, 1.1, 1.3);
  
  h1_deltaChi2_Neyman->Draw("same");  
  h1_deltaChi2_Neyman->SetLineColor(color_Neyman);
      
  h1_deltaChi2_Pearson->Draw("same");  
  h1_deltaChi2_Pearson->SetLineColor(color_Pearson);
      
  h1_deltaChi2_CNP->Draw("same");  
  h1_deltaChi2_CNP->SetLineColor(color_CNP);
      
  h1_deltaChi2_Gauss->Draw("same");  
  h1_deltaChi2_Gauss->SetLineColor(color_Gauss);

  h1_deltaChi2_Poisson->Draw("same axis");
  
  int ndf_deltaChi2 = 2;
  TF1 *f_deltaChi2 = new TF1("f_deltaChi2",
			    TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", ndf_deltaChi2, (deltaChi2_hgh-deltaChi2_low)/nbins_deltaChi2  ),
			    chi2_low, chi2_hgh);
  f_deltaChi2->Draw("same");
  f_deltaChi2->SetLineColor(kRed);
  
  TLegend *lg_deltaChi2 = new TLegend(0.55+0.07,0.56-0.1,0.83+0.07,0.88);
  lg_deltaChi2->SetBorderSize(0);
  lg_deltaChi2->AddEntry(  f_deltaChi2, TString::Format("chi2, ndf=%d", ndf_deltaChi2), "l");
  lg_deltaChi2->AddEntry(  h1_deltaChi2_Poisson, TString::Format("Poisson"), "l" );
  lg_deltaChi2->AddEntry(  h1_deltaChi2_Neyman,  TString::Format("Neyman"), "l" );
  lg_deltaChi2->AddEntry(  h1_deltaChi2_Pearson, TString::Format("Pearson"), "l" );
  lg_deltaChi2->AddEntry(  h1_deltaChi2_CNP,     TString::Format("CNP"), "l" );
  lg_deltaChi2->AddEntry(  h1_deltaChi2_Gauss,   TString::Format("Gauss"), "l" );			    
  lg_deltaChi2->Draw("same");
  lg_deltaChi2->SetTextSize(0.055);
  lg_deltaChi2->SetTextFont(42);

  //////////////////////////////////////////////////////////////////////////////////////// h1_p0_Poisson
  
  roostr = "canv_h1_p0_Poisson";
  TCanvas *canv_h1_p0_Poisson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h1_p0_Poisson, 0.18, 0.1, 0.18, 0.06);  
  h1_p0_Poisson->Draw();  
  h1_p0_Poisson->SetLineColor(color_Poisson);
  roo_title(h1_p0_Poisson, "#hat{p}_{0}", "Entries");
  roo_center_title(h1_p0_Poisson);
  roo_offset(h1_p0_Poisson, 1.1, 1.3);
  
  h1_p0_Neyman->Draw("same");  
  h1_p0_Neyman->SetLineColor(color_Neyman);
      
  h1_p0_Pearson->Draw("same");  
  h1_p0_Pearson->SetLineColor(color_Pearson);
      
  h1_p0_CNP->Draw("same");  
  h1_p0_CNP->SetLineColor(color_CNP);
      
  h1_p0_Gauss->Draw("same");  
  h1_p0_Gauss->SetLineColor(color_Gauss);

  h1_p0_Poisson->Draw("same axis");
  
  // TLegend *lg_p0 = new TLegend(0.55-0.01,0.56,0.83+0.07-0.01,0.88);
  // lg_p0->SetBorderSize(0);
  // lg_p0->AddEntry(  h1_p0_Poisson, TString::Format("Poisson,  Bias %5.2f", bias_p0_Poisson*100)+"%", "l" );
  // lg_p0->AddEntry(  h1_p0_Neyman, TString::Format("Neyman, Bias %5.2f", bias_p0_Neyman*100)+"%", "l" );
  // lg_p0->AddEntry(  h1_p0_Pearson, TString::Format("Pearson, Bias %5.2f", bias_p0_Pearson*100)+"%", "l" );
  // lg_p0->AddEntry(  h1_p0_CNP, TString::Format("CNP,       Bias %5.2f", bias_p0_CNP*100)+"%", "l" );
  // lg_p0->AddEntry(  h1_p0_Gauss, TString::Format("Gauss,    Bias %5.2f", bias_p0_Gauss*100)+"%", "l" );
		    
  // lg_p0->Draw("same");
  // lg_p0->SetTextSize(0.05);
  // lg_p0->SetTextFont(42);
  
  TLegend *lg_p0 = new TLegend(0.55-0.01,0.56,0.83+0.07-0.01,0.88);
  lg_p0->SetBorderSize(0);
  lg_p0->AddEntry(  h1_p0_Poisson, TString::Format("Poisson,  Abs Bias %5.2f", bias_p0_Poisson*p0_true), "l" );
  lg_p0->AddEntry(  h1_p0_Neyman, TString::Format("Neyman, Abs Bias %5.2f", bias_p0_Neyman*p0_true), "l" );
  lg_p0->AddEntry(  h1_p0_Pearson, TString::Format("Pearson, Abs Bias %5.2f", bias_p0_Pearson*p0_true), "l" );
  lg_p0->AddEntry(  h1_p0_CNP, TString::Format("CNP,       Abs Bias %5.2f", bias_p0_CNP*p0_true), "l" );
  lg_p0->AddEntry(  h1_p0_Gauss, TString::Format("Gauss,    Abs Bias %5.2f", bias_p0_Gauss*p0_true), "l" );
		    
  lg_p0->Draw("same");
  lg_p0->SetTextSize(0.05);
  lg_p0->SetTextFont(42);

  //////////////////////////////////////////////////////////////////////////////////////// h1_p1_Poisson
 
  roostr = "canv_h1_p1_Poisson";
  TCanvas *canv_h1_p1_Poisson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h1_p1_Poisson, 0.18, 0.1, 0.18, 0.06);  
  h1_p1_Poisson->Draw();  
  h1_p1_Poisson->SetLineColor(color_Poisson);
  roo_title(h1_p1_Poisson, "#hat{p}_{1}", "Entries");
  roo_center_title(h1_p1_Poisson);
  roo_offset(h1_p1_Poisson, 1.1, 1.3);
  
  h1_p1_Neyman->Draw("same");  
  h1_p1_Neyman->SetLineColor(color_Neyman);
      
  h1_p1_Pearson->Draw("same");  
  h1_p1_Pearson->SetLineColor(color_Pearson);
      
  h1_p1_CNP->Draw("same");  
  h1_p1_CNP->SetLineColor(color_CNP);
      
  h1_p1_Gauss->Draw("same");  
  h1_p1_Gauss->SetLineColor(color_Gauss);

  h1_p1_Poisson->Draw("same axis");
  
  // TLegend *lg_p1 = new TLegend(0.55-0.01,0.56,0.83+0.07-0.01,0.88);
  // lg_p1->SetBorderSize(0);
  // lg_p1->AddEntry(  h1_p1_Poisson, TString::Format("Poisson,  Bias %5.2f", bias_p1_Poisson*100)+"%", "l" );
  // lg_p1->AddEntry(  h1_p1_Neyman, TString::Format("Neyman, Bias %5.2f", bias_p1_Neyman*100)+"%", "l" );
  // lg_p1->AddEntry(  h1_p1_Pearson, TString::Format("Pearson, Bias %5.2f", bias_p1_Pearson*100)+"%", "l" );
  // lg_p1->AddEntry(  h1_p1_CNP, TString::Format("CNP,       Bias %5.2f", bias_p1_CNP*100)+"%", "l" );
  // lg_p1->AddEntry(  h1_p1_Gauss, TString::Format("Gauss,    Bias %5.2f", bias_p1_Gauss*100)+"%", "l" );
		    
  // lg_p1->Draw("same");
  // lg_p1->SetTextSize(0.05);
  // lg_p1->SetTextFont(42);

  TLegend *lg_p1 = new TLegend(0.55-0.01,0.56,0.83+0.07-0.01,0.88);
  lg_p1->SetBorderSize(0);
  lg_p1->AddEntry(  h1_p1_Poisson, TString::Format("Poisson,  Abs Bias %5.2f", bias_p1_Poisson*p1_true ), "l" );
  lg_p1->AddEntry(  h1_p1_Neyman, TString::Format("Neyman, Abs Bias %5.2f", bias_p1_Neyman*p1_true ), "l" );
  lg_p1->AddEntry(  h1_p1_Pearson, TString::Format("Pearson, Abs Bias %5.2f", bias_p1_Pearson*p1_true ), "l" );
  lg_p1->AddEntry(  h1_p1_CNP, TString::Format("CNP,       Abs Bias %5.2f", bias_p1_CNP*p1_true ), "l" );
  lg_p1->AddEntry(  h1_p1_Gauss, TString::Format("Gauss,    Abs Bias %5.2f", bias_p1_Gauss*p1_true ), "l" );
		    
  lg_p1->Draw("same");
  lg_p1->SetTextSize(0.05);
  lg_p1->SetTextFont(42);

  //////////////////////////////////////////////////////////////////////////////////////// h2_p1vsp0_Poisson

  TLine *line_p0 = new TLine( 0,low_2d_p1, 0, hgh_2d_p1 );
  TLine *line_p1 = new TLine( low_2d_p0, 0, hgh_2d_p0, 0);
  line_p0->SetLineStyle(7);
  line_p1->SetLineStyle(7);
      
  roostr = "canv_h2_p1vsp0_Poisson";
  TCanvas *canv_h2_p1vsp0_Poisson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h2_p1vsp0_Poisson, 0.18, 0.1, 0.18, 0.18);
  canv_h2_p1vsp0_Poisson->SetLogz();
  h2_p1vsp0_Poisson->Draw("colz");  
  h2_p1vsp0_Poisson->SetLineColor(color_Poisson);
  roo_title(h2_p1vsp0_Poisson, "(#hat{p}_{0} - p_{0}) / p_{0} [%]", "(#hat{p}_{1} - p_{1}) / p_{1} [%]");
  roo_center_title(h2_p1vsp0_Poisson);
  roo_offset(h2_p1vsp0_Poisson, 1.1, 1.3);
  line_p0->Draw("same");
  line_p1->Draw("same");
  
  roostr = "canv_h2_p1vsp0_Neyman";
  TCanvas *canv_h2_p1vsp0_Neyman = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h2_p1vsp0_Neyman, 0.18, 0.1, 0.18, 0.18);
  canv_h2_p1vsp0_Neyman->SetLogz();
  h2_p1vsp0_Neyman->Draw("colz");  
  h2_p1vsp0_Neyman->SetLineColor(color_Neyman);
  roo_title(h2_p1vsp0_Neyman, "(#hat{p}_{0} - p_{0}) / p_{0} [%]", "(#hat{p}_{1} - p_{1}) / p_{1} [%]");
  roo_center_title(h2_p1vsp0_Neyman);
  roo_offset(h2_p1vsp0_Neyman, 1.1, 1.3);
  line_p0->Draw("same");
  line_p1->Draw("same");
  
  roostr = "canv_h2_p1vsp0_Pearson";
  TCanvas *canv_h2_p1vsp0_Pearson = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h2_p1vsp0_Pearson, 0.18, 0.1, 0.18, 0.18);
  canv_h2_p1vsp0_Pearson->SetLogz();
  h2_p1vsp0_Pearson->Draw("colz");  
  h2_p1vsp0_Pearson->SetLineColor(color_Pearson);
  roo_title(h2_p1vsp0_Pearson, "(#hat{p}_{0} - p_{0}) / p_{0} [%]", "(#hat{p}_{1} - p_{1}) / p_{1} [%]");
  roo_center_title(h2_p1vsp0_Pearson);
  roo_offset(h2_p1vsp0_Pearson, 1.1, 1.3);
  line_p0->Draw("same");
  line_p1->Draw("same");
  
  roostr = "canv_h2_p1vsp0_CNP";
  TCanvas *canv_h2_p1vsp0_CNP = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h2_p1vsp0_CNP, 0.18, 0.1, 0.18, 0.18);
  canv_h2_p1vsp0_CNP->SetLogz();
  h2_p1vsp0_CNP->Draw("colz");  
  h2_p1vsp0_CNP->SetLineColor(color_CNP);
  roo_title(h2_p1vsp0_CNP, "(#hat{p}_{0} - p_{0}) / p_{0} [%]", "(#hat{p}_{1} - p_{1}) / p_{1} [%]");
  roo_center_title(h2_p1vsp0_CNP);
  roo_offset(h2_p1vsp0_CNP, 1.1, 1.3);
  line_p0->Draw("same");
  line_p1->Draw("same");
  
  roostr = "canv_h2_p1vsp0_Gauss";
  TCanvas *canv_h2_p1vsp0_Gauss = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_h2_p1vsp0_Gauss, 0.18, 0.1, 0.18, 0.18);
  canv_h2_p1vsp0_Gauss->SetLogz();
  h2_p1vsp0_Gauss->Draw("colz");  
  h2_p1vsp0_Gauss->SetLineColor(color_Gauss);
  roo_title(h2_p1vsp0_Gauss, "(#hat{p}_{0} - p_{0}) / p_{0} [%]", "(#hat{p}_{1} - p_{1}) / p_{1} [%]");
  roo_center_title(h2_p1vsp0_Gauss);
  roo_offset(h2_p1vsp0_Gauss, 1.1, 1.3);
  line_p0->Draw("same");
  line_p1->Draw("same");
  
  ///////////////////////////////////////////

  canv_h1_chi2true_Poisson->SaveAs("canv_h1_chi2true_Poisson.png");
  
  canv_h1_deltaChi2_Poisson->SaveAs("canv_h1_deltaChi2_Poisson.png");
  
  canv_h1_p0_Poisson->SaveAs("canv_h1_p0_Poisson.png");
  canv_h1_p1_Poisson->SaveAs("canv_h1_p1_Poisson.png");

  canv_h2_p1vsp0_Poisson->SaveAs("canv_h2_p1vsp0_Poisson.png");
  canv_h2_p1vsp0_Neyman->SaveAs("canv_h2_p1vsp0_Neyman.png");
  canv_h2_p1vsp0_Pearson->SaveAs("canv_h2_p1vsp0_Pearson.png");
  canv_h2_p1vsp0_CNP->SaveAs("canv_h2_p1vsp0_CNP.png");
  canv_h2_p1vsp0_Gauss->SaveAs("canv_h2_p1vsp0_Gauss.png");
  
}
