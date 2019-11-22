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

void paper_plot_table1()
{
  gROOT->ProcessLine(".x ./DrawOption.cc");

  int color_Poisson = kCyan;
  int color_Neyman = kBlue;
  int color_Pearson = kRed;
  int color_Gauss = kOrange;
  int color_CNP = kBlack;
  
  // color_Neyman = 1;
  // color_Pearson = 1;
  // color_Poisson = 1;
  // color_Gauss = 1;
  // color_CNP = 1;

  TString roostr = "";

  TFile *roofile = new TFile("file_table1_50M.root", "read");
  TTree *tree_bias = (TTree*)roofile->Get("tree_bias");
  
  // Declaration of leaf types
  Int_t           ndf_val;
  Double_t        mutrue;
  Double_t        muhat_Poisson;
  Double_t        muhat_Neyman;
  Double_t        muhat_Pearson;
  Double_t        muhat_Gauss;
  Double_t        muhat_CNP;

  // List of branches
  TBranch        *b_ndf_val;   //!
  TBranch        *b_mutrue;   //!
  TBranch        *b_muhat_Poisson;   //!
  TBranch        *b_muhat_Neyman;   //!
  TBranch        *b_muhat_Pearson;   //!
  TBranch        *b_muhat_Gauss;   //!
  TBranch        *b_muhat_CNP;   //!

  // Set branch addresses and branch pointers
  tree_bias->SetBranchAddress("ndf_val", &ndf_val, &b_ndf_val);
  tree_bias->SetBranchAddress("mutrue", &mutrue, &b_mutrue);
  tree_bias->SetBranchAddress("muhat_Poisson", &muhat_Poisson, &b_muhat_Poisson);
  tree_bias->SetBranchAddress("muhat_Neyman", &muhat_Neyman, &b_muhat_Neyman);
  tree_bias->SetBranchAddress("muhat_Pearson", &muhat_Pearson, &b_muhat_Pearson);
  tree_bias->SetBranchAddress("muhat_Gauss", &muhat_Gauss, &b_muhat_Gauss);
  tree_bias->SetBranchAddress("muhat_CNP", &muhat_CNP, &b_muhat_CNP);

  int entries = tree_bias->GetEntries();

  //////////////////////////

  TGraph *gh15_rel_Poisson = new TGraph(); roostr = "gh15_rel_Poisson"; gh15_rel_Poisson->SetName(roostr);
  TGraph *gh15_rel_Neyman = new TGraph(); roostr = "gh15_rel_Neyman"; gh15_rel_Neyman->SetName(roostr);
  TGraph *gh15_rel_Pearson = new TGraph(); roostr = "gh15_rel_Pearson"; gh15_rel_Pearson->SetName(roostr);
  TGraph *gh15_rel_Gauss = new TGraph(); roostr = "gh15_rel_Gauss"; gh15_rel_Gauss->SetName(roostr);
  TGraph *gh15_rel_CNP = new TGraph(); roostr = "gh15_rel_CNP"; gh15_rel_CNP->SetName(roostr);

  TGraph *gh15_abs_Poisson = new TGraph(); roostr = "gh15_abs_Poisson"; gh15_abs_Poisson->SetName(roostr);
  TGraph *gh15_abs_Neyman = new TGraph(); roostr = "gh15_abs_Neyman"; gh15_abs_Neyman->SetName(roostr);
  TGraph *gh15_abs_Pearson = new TGraph(); roostr = "gh15_abs_Pearson"; gh15_abs_Pearson->SetName(roostr);
  TGraph *gh15_abs_Gauss = new TGraph(); roostr = "gh15_abs_Gauss"; gh15_abs_Gauss->SetName(roostr);
  TGraph *gh15_abs_CNP = new TGraph(); roostr = "gh15_abs_CNP"; gh15_abs_CNP->SetName(roostr);

  int line15 = 0;

  //////////////////////////

  TGraph *gh150_rel_Poisson = new TGraph(); roostr = "gh150_rel_Poisson"; gh150_rel_Poisson->SetName(roostr);
  TGraph *gh150_rel_Neyman = new TGraph(); roostr = "gh150_rel_Neyman"; gh150_rel_Neyman->SetName(roostr);
  TGraph *gh150_rel_Pearson = new TGraph(); roostr = "gh150_rel_Pearson"; gh150_rel_Pearson->SetName(roostr);
  TGraph *gh150_rel_Gauss = new TGraph(); roostr = "gh150_rel_Gauss"; gh150_rel_Gauss->SetName(roostr);
  TGraph *gh150_rel_CNP = new TGraph(); roostr = "gh150_rel_CNP"; gh150_rel_CNP->SetName(roostr);

  TGraph *gh150_abs_Poisson = new TGraph(); roostr = "gh150_abs_Poisson"; gh150_abs_Poisson->SetName(roostr);
  TGraph *gh150_abs_Neyman = new TGraph(); roostr = "gh150_abs_Neyman"; gh150_abs_Neyman->SetName(roostr);
  TGraph *gh150_abs_Pearson = new TGraph(); roostr = "gh150_abs_Pearson"; gh150_abs_Pearson->SetName(roostr);
  TGraph *gh150_abs_Gauss = new TGraph(); roostr = "gh150_abs_Gauss"; gh150_abs_Gauss->SetName(roostr);
  TGraph *gh150_abs_CNP = new TGraph(); roostr = "gh150_abs_CNP"; gh150_abs_CNP->SetName(roostr);

  int line150 = 0;

  //////////////////////////
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree_bias->GetEntry(ientry);

    if( ndf_val==4 || ndf_val==6 ) continue;
    
    if( mutrue==15 ) {
      line15++;      
      gh15_abs_Poisson->SetPoint(line15-1, ndf_val, muhat_Poisson-mutrue);
      gh15_abs_Neyman->SetPoint(line15-1, ndf_val, muhat_Neyman-mutrue);
      gh15_abs_Pearson->SetPoint(line15-1, ndf_val, muhat_Pearson-mutrue);
      gh15_abs_Gauss->SetPoint(line15-1, ndf_val, muhat_Gauss-mutrue);
      gh15_abs_CNP->SetPoint(line15-1, ndf_val, muhat_CNP-mutrue);

      gh15_rel_Poisson->SetPoint(line15-1, ndf_val, muhat_Poisson/mutrue*100 - 100);
      gh15_rel_Neyman->SetPoint(line15-1, ndf_val, muhat_Neyman/mutrue*100 - 100);
      gh15_rel_Pearson->SetPoint(line15-1, ndf_val, muhat_Pearson/mutrue*100 - 100);
      gh15_rel_Gauss->SetPoint(line15-1, ndf_val, muhat_Gauss/mutrue*100 - 100);
      gh15_rel_CNP->SetPoint(line15-1, ndf_val, muhat_CNP/mutrue*100 - 100);
    }

    
    if( mutrue==150 ) {
      line150++;      
      gh150_abs_Poisson->SetPoint(line150-1, ndf_val, muhat_Poisson-mutrue);
      gh150_abs_Neyman->SetPoint(line150-1, ndf_val, muhat_Neyman-mutrue);
      gh150_abs_Pearson->SetPoint(line150-1, ndf_val, muhat_Pearson-mutrue);
      gh150_abs_Gauss->SetPoint(line150-1, ndf_val, muhat_Gauss-mutrue);
      gh150_abs_CNP->SetPoint(line150-1, ndf_val, muhat_CNP-mutrue);

      gh150_rel_Poisson->SetPoint(line150-1, ndf_val, muhat_Poisson/mutrue*100 - 100);
      gh150_rel_Neyman->SetPoint(line150-1, ndf_val, muhat_Neyman/mutrue*100 - 100);
      gh150_rel_Pearson->SetPoint(line150-1, ndf_val, muhat_Pearson/mutrue*100 - 100);
      gh150_rel_Gauss->SetPoint(line150-1, ndf_val, muhat_Gauss/mutrue*100 - 100);
      gh150_rel_CNP->SetPoint(line150-1, ndf_val, muhat_CNP/mutrue*100 - 100);
    }
      
  }

  ///////////////////////////////////////////////////////////////////////

  TF1 *f1_zero = new TF1("f1_zero", "0", 0, 1000);
  f1_zero->SetLineStyle(1);
  f1_zero->SetLineColor(kGray);
    
  roostr = "canv_abs_AA";
  TCanvas *canv_abs_AA = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_abs_AA, 0.18, 0.09, 0.2, 0.05);

  roostr = "h1_basic_AA";
  TH1D *h1_basic_AA = new TH1D(roostr, roostr, 100, 0, 105);
  for(int ibin=1; ibin<=100; ibin++) h1_basic_AA->SetBinContent(ibin, -200);
  
  h1_basic_AA->GetYaxis()->SetRangeUser(-0.2, 0.05);
  h1_basic_AA->Draw();
  h1_basic_AA->SetLineWidth(0.1);
  roo_title(h1_basic_AA, "Number of measurements", "#hat{#mu} - #mu");
  roo_center_title(h1_basic_AA);
  roo_offset(h1_basic_AA, 1.1, 1.4);
  h1_basic_AA->GetXaxis()->SetLabelSize(0.07);
  h1_basic_AA->GetXaxis()->SetTitleSize(0.07);
  h1_basic_AA->GetYaxis()->SetLabelSize(0.07);
  h1_basic_AA->GetYaxis()->SetTitleSize(0.07);

  double markerSize_AA = 1.4;
  
  //////
  gh150_abs_Gauss->Draw("same pl");
  gh150_abs_Gauss->SetLineColor(color_Gauss);
  gh150_abs_Gauss->SetLineStyle(1);
  gh150_abs_Gauss->SetMarkerStyle(25);
  gh150_abs_Gauss->SetMarkerSize(markerSize_AA);
  gh150_abs_Gauss->SetMarkerColor(color_Gauss);
  
  gh15_abs_Gauss->Draw("same pl");
  gh15_abs_Gauss->SetLineColor(color_Gauss);
  gh15_abs_Gauss->SetLineStyle(7);
  gh15_abs_Gauss->SetMarkerStyle(5);
  gh15_abs_Gauss->SetMarkerSize(markerSize_AA);
  gh15_abs_Gauss->SetMarkerColor(color_Gauss);

  //////
  gh150_abs_CNP->Draw("same pl");
  gh150_abs_CNP->SetLineColor(color_CNP);
  gh150_abs_CNP->SetLineStyle(1);
  gh150_abs_CNP->SetMarkerStyle(24);
  gh150_abs_CNP->SetMarkerSize(markerSize_AA);
  gh150_abs_CNP->SetMarkerColor(color_CNP);
  
  gh15_abs_CNP->Draw("same pl");
  gh15_abs_CNP->SetLineColor(color_CNP);
  gh15_abs_CNP->SetLineStyle(7);
  gh15_abs_CNP->SetMarkerStyle(26);
  gh15_abs_CNP->SetMarkerSize(markerSize_AA);
  gh15_abs_CNP->SetMarkerColor(color_CNP);

  f1_zero->Draw("same");
  h1_basic_AA->Draw("same axis");
    
  TLegend *lg_basic_AA = new TLegend(0.4922049,0.2467949,0.9209354,0.5480769);
  lg_basic_AA->SetBorderSize(1);
  lg_basic_AA->SetTextFont(42);
  lg_basic_AA->AddEntry(gh15_abs_Gauss,   TString::Format("Gauss, #mu = 15"), "pl" );
  lg_basic_AA->AddEntry(gh150_abs_Gauss,  TString::Format("Gauss, #mu = 150"), "pl" );    
  lg_basic_AA->AddEntry(gh15_abs_CNP,     TString::Format("CNP,   #mu = 15"), "pl" );
  lg_basic_AA->AddEntry(gh150_abs_CNP,    TString::Format("CNP,   #mu = 150"), "pl" );

  lg_basic_AA->Draw("same");
  lg_basic_AA->SetTextFont(42);
  lg_basic_AA->SetTextSize(0.068);

  canv_abs_AA->SaveAs("canv_table_CNP_Gauss.pdf");
  
  ///////////////////////////////////////////////////////////////////////
  
  roostr = "canv_abs_BB";
  TCanvas *canv_abs_BB = new TCanvas(roostr, roostr, 900, 650);
  roo_canv_margin(canv_abs_BB, 0.18, 0.09, 0.2, 0.05);

  roostr = "h1_basic_BB";
  TH1D *h1_basic_BB = new TH1D(roostr, roostr, 100, 0, 105);
  for(int ibin=1; ibin<=100; ibin++) h1_basic_BB->SetBinContent(ibin, -200);
  
  h1_basic_BB->GetYaxis()->SetRangeUser(-1.2, 0.6);
  h1_basic_BB->Draw();
  h1_basic_BB->SetLineWidth(0.1);
  roo_title(h1_basic_BB, "Number of measurements", "#hat{#mu} - #mu");
  roo_center_title(h1_basic_BB);
  roo_offset(h1_basic_BB, 1.1, 1.2);
  h1_basic_BB->GetXaxis()->SetLabelSize(0.07);
  h1_basic_BB->GetXaxis()->SetTitleSize(0.07);
  h1_basic_BB->GetYaxis()->SetLabelSize(0.07);
  h1_basic_BB->GetYaxis()->SetTitleSize(0.07);
  h1_basic_BB->GetYaxis()->SetNdivisions(508);
  
  //////
  gh150_abs_Neyman->Draw("same pl");
  gh150_abs_Neyman->SetLineColor(color_Neyman);
  gh150_abs_Neyman->SetLineStyle(1);
  gh150_abs_Neyman->SetMarkerStyle(24);
  gh150_abs_Neyman->SetMarkerSize(markerSize_AA);
  gh150_abs_Neyman->SetMarkerColor(color_Neyman);
  
  gh15_abs_Neyman->Draw("same pl");
  gh15_abs_Neyman->SetLineColor(color_Neyman);
  gh15_abs_Neyman->SetLineStyle(7);
  gh15_abs_Neyman->SetMarkerStyle(26);
  gh15_abs_Neyman->SetMarkerSize(markerSize_AA);
  gh15_abs_Neyman->SetMarkerColor(color_Neyman);

  //////
  gh150_abs_Pearson->Draw("same pl");
  gh150_abs_Pearson->SetLineColor(color_Pearson);
  gh150_abs_Pearson->SetLineStyle(1);
  gh150_abs_Pearson->SetMarkerStyle(25);
  gh150_abs_Pearson->SetMarkerSize(markerSize_AA);
  gh150_abs_Pearson->SetMarkerColor(color_Pearson);
  
  gh15_abs_Pearson->Draw("same pl");
  gh15_abs_Pearson->SetLineColor(color_Pearson);
  gh15_abs_Pearson->SetLineStyle(7);
  gh15_abs_Pearson->SetMarkerStyle(5);
  gh15_abs_Pearson->SetMarkerSize(markerSize_AA);
  gh15_abs_Pearson->SetMarkerColor(color_Pearson);

  f1_zero->Draw("same");
  
  h1_basic_BB->Draw("same axis");
  
  TLegend *lg_basic_BB = new TLegend(0.4365256,0.3253205,0.9242762,0.6266026);
  lg_basic_BB->SetBorderSize(1);
  lg_basic_BB->SetTextFont(42);
  lg_basic_BB->AddEntry(gh15_abs_Neyman,    TString::Format("Neyman, #mu = 15"), "pl" );
  lg_basic_BB->AddEntry(gh150_abs_Neyman,   TString::Format("Neyman, #mu = 150"), "pl" );  
  lg_basic_BB->AddEntry(gh15_abs_Pearson,   TString::Format("Pearson, #mu = 15"), "pl" );
  lg_basic_BB->AddEntry(gh150_abs_Pearson,  TString::Format("Pearson, #mu = 150"), "pl" );  

  lg_basic_BB->Draw("same");
  lg_basic_BB->SetTextFont(42);
  lg_basic_BB->SetTextSize(0.068);

  canv_abs_BB->SaveAs("canv_table_Pearson_Neyman.pdf");
  
  
}
