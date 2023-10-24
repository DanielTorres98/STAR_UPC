// Code to generate plots from the Analysis output file.

#include <TMath.h>

#include "utils/Fit_functions.h"

void Phi_distribution_fit(){

    TFile * InputFile = new TFile( "UpcOutput_Rho.root", "READ" );

    TH1D* H1D_phi_0_60MeV = (TH1D*) InputFile->Get("hphi");
    H1D_phi_0_60MeV->SetTitle("#phi, p_{T}<60 (MeV/c^{2})");
    H1D_phi_0_60MeV->SetMarkerStyle(5);
    H1D_phi_0_60MeV->GetYaxis()->SetTitle("dN/d#phi (norm. to unity)");
    H1D_phi_0_60MeV->GetYaxis()->SetTitleSize(0.05);
    H1D_phi_0_60MeV->GetXaxis()->SetTitle("#phi");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(1, "-#pi");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(25, "-#pi/2");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(50, "0");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(75, "#pi/2");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(100, "#pi");
    H1D_phi_0_60MeV->LabelsOption("h"); // Change x axis to horizontal orientation
    H1D_phi_0_60MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_phi_0_60MeV->GetXaxis()->SetLabelSize(0.06);

    TF1* fit_phi = new TF1("fit_phi", Interference_phi_fit, -TMath::Pi(), TMath::Pi(), 1); // One parameter (A)

    // Optionally, set some initial parameter values to guide the fit. Initial guess for A (amplitude)
    fit_phi->SetParameter(0, 1.0);

    // Perform the fit
    H1D_phi_0_60MeV->Fit("fit_phi", "ERS0", "", -3.14, 3.14);
    double fit_phi_params = fit_phi->GetParameter(0);
    TF1 *Interference_phi = new TF1("Interference_phi", Interference_phi_fit, -3.14, 3.14, 1);
    Interference_phi->SetParameters(fit_phi_params);
    Interference_phi->SetLineColor(kRed);
    Interference_phi->SetLineStyle(1);

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
    c1->SetLeftMargin(0.15);
 
    H1D_phi_0_60MeV->Draw("E1");
    Interference_phi->Draw("Same");
    c1->SaveAs("Images/phi_interference.pdf");
    c1->SaveAs("Images/phi_interference.png");

    delete c1;
}
