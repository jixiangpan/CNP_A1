void paper_plot_01_chi2scan(){

  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.07,"Y");
  gStyle->SetTitleSize(0.07,"Y");
  gStyle->SetLabelSize(0.07,"X");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetNdivisions(506,"X");
  gStyle->SetNdivisions(506,"Y");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  
  int color_Neyman = kRed;
  int color_Pearson = kBlue;
  int color_Poisson = kCyan;//;kGreen;
  int color_Gauss = kOrange;  
  int color_CNP = kBlack;

  const int Ndet = 10;

  int val_true = 15;

  TRandom3 *random = new TRandom3(1815896008+39); // 1815896008+9, 1815896008+12, 1815896008+31, 1815896008+37
  cout<<"seed "<<random->GetSeed()<<endl;

  double muhat = 0;
  
  Double_t num[Ndet];
  for (int i=0;i!=Ndet;i++){
    num[i] = random->Poisson(val_true);

    muhat += num[i]/Ndet;
  }

  cout<<endl<<" ---> Poisson muhat "<<muhat<<endl<<endl;


  TGraph *g_poisson = new TGraph();
  TGraph *g_neyman = new TGraph();
  TGraph *g_pearson = new TGraph();
  TGraph *g_cnp = new TGraph();
  TGraph *g_gauss = new TGraph();

  Double_t chi2_min[4]={1e9,1e9,1e9,1e9};
  Double_t mu_min[4]={1e9,1e9,1e9,1e9};
  
  double ylow = 0;
  double yhgh = 50;

  double xlow = 5;
  double xhgh = 30;  
  double nbins_x = 2000;
  double step = (xhgh-xlow)/nbins_x;
  
  for (Int_t i=0;i!=nbins_x;i++){
    //double mu = 8+(i+0.5)*30./1000.;
    double mu = xlow + i*step;
    
    double chi2 = 0;
    for (Int_t j=0;j!=Ndet;j++){
      chi2 += pow(mu-num[j],2)/mu;
    }
    if (chi2 < chi2_min[1] ){
      chi2_min[1] = chi2;
      mu_min[1] = mu;
    }
    g_pearson->SetPoint(i,mu,chi2);
    
    
    chi2 = 0;
    for (Int_t j=0;j!=Ndet;j++){
      chi2 += pow(mu-num[j],2)/num[j];
    }
    if (chi2 < chi2_min[2] ) {
      chi2_min[2] = chi2;
      mu_min[2] = mu;
    }
    g_neyman->SetPoint(i,mu,chi2);

    
    chi2 =0;
    for (Int_t j=0;j!=Ndet;j++){
      chi2 += (pow(mu-num[j],2)/mu*2+pow(mu-num[j],2)/num[j])/3.;
    }
    if (chi2 < chi2_min[3] ) {
      chi2_min[3] = chi2;
      mu_min[3] = mu;
    }
    g_cnp->SetPoint(i,mu,chi2);

    
    chi2 = 0;
    for (Int_t j=0;j!=Ndet;j++){
      chi2 += 2*(mu - num[j] + num[j] * log(num[j]/mu));

      if( num[j]==0 ) {
	cout<<endl<<" Meas = 0"<<endl<<endl;
      }
    }
    if (chi2 < chi2_min[0] ) {
      chi2_min[0] = chi2;
      mu_min[0] = mu;
    }
    g_poisson->SetPoint(i,mu,chi2);

    
    chi2 = 0;
    for (Int_t j=0;j!=Ndet;j++){
      double mi = num[j];
      double mp = 0.5*( sqrt(4*mi*mi+1) -1 );

      double p1 = pow( mu-mi,2 )/mu;
      double p2 = pow( mp-mi,2 )/mp;
      double p3 = log(mu/mp);
      
      chi2 += (p1-p2+p3);
    }
    g_gauss->SetPoint(i,mu,chi2);
    
  }

  TCanvas *canv_simple_example_1toy = new TCanvas("canv_simple_example_1toy", "canv_simple_example_1toy", 900, 650);
  canv_simple_example_1toy->SetRightMargin(0.05);
  canv_simple_example_1toy->SetTopMargin(0.09);
  canv_simple_example_1toy->SetBottomMargin(0.18);
  canv_simple_example_1toy->SetLeftMargin(0.12);

  g_poisson->Draw("AL");
  //g_poisson->GetXaxis()->SetRangeUser(9.1,21.9);
  g_poisson->GetXaxis()->SetRangeUser(5,27);
  g_poisson->GetYaxis()->SetRangeUser(ylow, yhgh);
  g_pearson->Draw("Lsame");
  g_neyman->Draw("Lsame");
  g_gauss->Draw("Lsame");
  g_cnp->Draw("Lsame");

  g_poisson->SetLineStyle(3);
  g_neyman->SetLineStyle(7);
  g_pearson->SetLineStyle(5);
  g_gauss->SetLineStyle(2);
  g_cnp->SetLineStyle(1);

  g_poisson->SetLineColor(color_Poisson);
  g_neyman->SetLineColor(color_Neyman);
  g_pearson->SetLineColor(color_Pearson);
  g_gauss->SetLineColor(color_Gauss);
  g_cnp->SetLineColor(color_CNP);
  
  g_poisson->SetLineWidth(2);
  g_neyman->SetLineWidth(2);
  g_pearson->SetLineWidth(2);
  g_gauss->SetLineWidth(2);
  g_cnp->SetLineWidth(2);


  //////////////////////// vertical line
  
  TLine *l1 = new TLine(mu_min[0],ylow,mu_min[0],yhgh);
  TLine *l2 = new TLine(mu_min[1],ylow,mu_min[1],yhgh);
  TLine *l3 = new TLine(mu_min[2],ylow,mu_min[2],yhgh);
  TLine *l4 = new TLine(mu_min[3],ylow,mu_min[3],yhgh);

  // l1->Draw("same");
  // l2->Draw("same");
  // l3->Draw("same");
  // l4->Draw("same");
  
  l1->SetLineColor(color_Poisson);
  l2->SetLineColor(color_Pearson);
  l3->SetLineColor(color_Neyman);
  l4->SetLineColor(color_CNP);
  
  l2->SetLineStyle(2);
  l3->SetLineStyle(3);
  l4->SetLineStyle(7);

  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l4->SetLineWidth(2);

  g_poisson->GetXaxis()->SetTitle("Test #mu");
  g_poisson->GetYaxis()->SetTitle("#chi^{2}");
  g_poisson->GetXaxis()->CenterTitle();
  g_poisson->GetYaxis()->CenterTitle();
  // g_poisson->GetXaxis()->SetLabelSize(0.06);
  // g_poisson->GetXaxis()->SetTitleSize(0.06);
  // g_poisson->GetYaxis()->SetLabelSize(0.06);
  // g_poisson->GetYaxis()->SetTitleSize(0.06);
  g_poisson->GetYaxis()->SetTitleOffset(0.85);
  g_poisson->GetXaxis()->SetTitleOffset(1.18);
  
  TLegend *le1 = new TLegend(0.3797327-0.01,0.5865385-0.075,0.6202673+0.01,0.8894231);
  le1->SetTextSize(0.07);
  le1->SetFillColor(10);
  le1->AddEntry(g_poisson,"Poisson","l");
  le1->AddEntry(g_gauss,"Gauss","l");
  le1->AddEntry(g_neyman,"Neyman","l");
  le1->AddEntry(g_pearson,"Pearson","l");  
  le1->AddEntry(g_cnp,"CNP","l");
  le1->Draw();
  
  //canv_simple_example_1toy->SaveAs("canv_simple_example_1toy.png");
  canv_simple_example_1toy->SaveAs("canv_simple_example_1toy.pdf");
  
}
