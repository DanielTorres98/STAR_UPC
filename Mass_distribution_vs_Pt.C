void Mass_distribution_vs_Pt(){

    TFile * InputFile = new TFile( "UpcOutput_Rho.root", "READ" );

    TH1D* H1D_MpiMass_0_20MeV = (TH1D*) InputFile->Get("hRhoMass_0_20");
    TH1D* H1D_MpiMass_20_40MeV = (TH1D*) InputFile->Get("hRhoMass_20_40");
    TH1D* H1D_MpiMass_40_60MeV = (TH1D*) InputFile->Get("hRhoMass_40_60");
    TH1D* H1D_MpiMass_60_80MeV = (TH1D*) InputFile->Get("hRhoMass_60_80");
    TH1D* H1D_MpiMass_80_100MeV = (TH1D*) InputFile->Get("hRhoMass_80_100");
    TH1D* H1D_MpiMass_100_150MeV = (TH1D*) InputFile->Get("hRhoMass_100_150");
    TH1D* H1D_MpiMass_150_200MeV = (TH1D*) InputFile->Get("hRhoMass_150_200");
    TH1D* H1D_MpiMass_200_1000MeV = (TH1D*) InputFile->Get("hRhoMass_200_1000");

    H1D_MpiMass_0_20MeV->SetTitle("#pi^{+}#pi^{-} mass distribution");
    H1D_MpiMass_0_20MeV->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
    H1D_MpiMass_0_20MeV->GetYaxis()->SetTitle("dN/dM_{#pi#pi} [1/(MeV/c^{2})]");
    H1D_MpiMass_0_20MeV->SetMarkerStyle(7);
    H1D_MpiMass_0_20MeV->GetYaxis()->SetTitleOffset(1.2);
    H1D_MpiMass_0_20MeV->GetYaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_20MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_20MeV->GetXaxis()->SetTitleOffset(0.8);


    // Canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLineWidth(2);

    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.11);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.05);
    gPad->SetFrameLineWidth(2);
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    TLatex tl;
    tl.SetTextSize(0.06);
    tl.SetNDC();
    c1->cd();

    H1D_MpiMass_0_20MeV->SetLineColor(kRed); H1D_MpiMass_0_20MeV->Draw();
    H1D_MpiMass_20_40MeV->SetLineColor(kOrange); H1D_MpiMass_20_40MeV->Draw("same");
    H1D_MpiMass_40_60MeV->SetLineColor(kYellow); H1D_MpiMass_40_60MeV->Draw("same");
    H1D_MpiMass_60_80MeV->SetLineColor(kGreen); H1D_MpiMass_60_80MeV->Draw("same");
    H1D_MpiMass_80_100MeV->SetLineColor(kCyan); H1D_MpiMass_80_100MeV->Draw("same");
    H1D_MpiMass_100_150MeV->SetLineColor(kBlue); H1D_MpiMass_100_150MeV->Draw("same");
    H1D_MpiMass_150_200MeV->SetLineColor(kAzure); H1D_MpiMass_150_200MeV->Draw("same");
    H1D_MpiMass_200_1000MeV->SetLineColor(kBlack); H1D_MpiMass_200_1000MeV->Draw("same");

    TLegend *legend = new TLegend(0.7, 0.7, 0.99, 0.99);
    legend->AddEntry(H1D_MpiMass_0_20MeV, "0.00 GeV < P_{T} < 0.02 GeV", "l");
    legend->AddEntry(H1D_MpiMass_20_40MeV, "0.02 GeV < P_{T} < 0.04 GeV", "l");
    legend->AddEntry(H1D_MpiMass_40_60MeV, "0.04 GeV < P_{T} < 0.06 GeV", "l");
    legend->AddEntry(H1D_MpiMass_60_80MeV, "0.06 GeV < P_{T} < 0.08 GeV", "l");
    legend->AddEntry(H1D_MpiMass_80_100MeV, "0.08 GeV < P_{T} < 0.10 GeV", "l");
    legend->AddEntry(H1D_MpiMass_100_150MeV, "0.10 GeV < P_{T} < 0.15 GeV", "l");
    legend->AddEntry(H1D_MpiMass_150_200MeV, "0.15 GeV < P_{T} < 0.20 GeV", "l");
    legend->AddEntry(H1D_MpiMass_200_1000MeV, "0.20 GeV < P_{T} < 1.00 GeV", "l");
    legend->SetTextSize(0.025);
    legend->Draw();

    c1->SetLeftMargin(0.15); 
    c1->SaveAs("Images/MpiMass_vs_Pt.pdf");
    c1->SaveAs("Images/MpiMass_vs_Pt.png");
    delete c1;


}