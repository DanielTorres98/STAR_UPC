// Code to generate plots from the Analysis output file.

#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770
#define RhoLifeTime        0.148
#define OmegaPdgMass       0.7824
#define OmegaLifeTime      0.017
#define speedoflight       29.9792

#include <complex>
#include <cmath>
#include <TMath.h>
#include </Users/dt36/Documents/Research/BNL/UPC/STAR_UPC/utils/Fit_functions.h>

void AddEmptyLegendEntry(TLegend *legend, const char *text = "", int color = 0);

void Make_Plots(){

    TFile * InputFile = new TFile( "UpcOutput_Rho.root", "READ" );

    TH1D* H1D_MpiMass_0_100MeV = (TH1D*) InputFile->Get("hMpiMass_0_100MeV");
    H1D_MpiMass_0_100MeV->SetTitle("#pi^{+}#pi^{-} mass distribution p_{T} < 100 MeV");
    H1D_MpiMass_0_100MeV->GetYaxis()->SetRangeUser(-0.005, 0.04);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitle("dN/dM_{#pi#pi} (1/MeV/c^{2})");
    H1D_MpiMass_0_100MeV->SetMarkerStyle(24);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleOffset(1.2);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleOffset(0.8);

    TH1D* H1D_phi_0_60MeV = (TH1D*) InputFile->Get("hphi");
    H1D_phi_0_60MeV->SetMarkerStyle(24);
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(1, "-#pi");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(25, "-#pi/2");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(50, "0");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(75, "#pi/2");
    H1D_phi_0_60MeV->GetXaxis()->SetBinLabel(100, "#pi");
    H1D_phi_0_60MeV->LabelsOption("h");
    H1D_phi_0_60MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_phi_0_60MeV->GetXaxis()->SetLabelSize(0.06);

    // Fit
    TF1 *Soding_Param = new TF1("Soding_Param", SodingEqn, RhoPdgMass-RhoLifeTime, RhoPdgMass+RhoLifeTime, 4);

    // Set parameter names
    Soding_Param->SetParNames("A", "B", "M_#rho", "#Gamma_#rho", "C", "D", "E");

    // Set initial parameter values for A, B RhoMass, GammaRho respectively.
    Soding_Param->SetParameters(1.0, 1.0, RhoPdgMass, RhoLifeTime);

    // Perform the first fit and set its line and marker color
    H1D_MpiMass_0_100MeV->Fit("Soding_Param", "ERS0", "", 0.6, 0.9);

    // Perform the second fit and set its line and marker color
    // fitParams = (A, B, MRho, Gammarho)
    double fit_Soding_Params[4]; Soding_Param->GetParameters(fit_Soding_Params);
    double fit_Soding_Errors[4]; for (int i = 0; i < 4; i++) fit_Soding_Errors[i] = Soding_Param->GetParError(i);

    TF1 *Soding_Param_all = new TF1("Soding_Param_all", SodingEqn, 0.44, 1.3, 4);
    Soding_Param_all->SetLineColor(kBlack);
    Soding_Param_all->SetLineStyle(1);
    Soding_Param_all->SetParameters(fit_Soding_Params[0], fit_Soding_Params[1], 
                                    fit_Soding_Params[2], fit_Soding_Params[3]);

    TF1 *Soding_Param_term1 = new TF1("Soding_Param_term1", SodingEqn_addition_term, 0.44, 1.3, 4);
    Soding_Param_term1->SetLineColor(kBlack);
    Soding_Param_term1->SetLineStyle(2);
    Soding_Param_term1->SetParameters(fit_Soding_Params[0], fit_Soding_Params[1],
                                      fit_Soding_Params[2], fit_Soding_Params[3]);

    TF1 *Soding_Param_term2 = new TF1("SodingEqn_addition_term", SodingEqn_interference_term, 0.44, 1.3, 4);
    Soding_Param_term2->SetLineColor(kBlack);
    Soding_Param_term2->SetLineStyle(3);
    Soding_Param_term2->SetParameters(fit_Soding_Params[0], fit_Soding_Params[1],
                                      fit_Soding_Params[2], fit_Soding_Params[3]);

    // Fit
    TF1 *BW_rho_omega_photoproduction_param = new TF1("BW_rho_omega_photoproduction_param", 
                            BW_rho_omega_photoproduction, RhoPdgMass-RhoLifeTime, RhoPdgMass+RhoLifeTime, 8);

    // Set parameter names
    BW_rho_omega_photoproduction_param->SetParNames("M_#rho", "#Gamma_#rho", "A_#rho", "B_{#pi#pi}", "C_#omega",
                                                  "M_#omega", "#Gamma_#omega", "#phi_#omega");

    // Set initial parameter values for A, B RhoMass, GammaRho respectively.
    BW_rho_omega_photoproduction_param->SetParameters(RhoPdgMass, RhoLifeTime, 1.0, 1.0, 1.0, OmegaPdgMass, OmegaLifeTime,
                                                     0.0);

    // Perform the first fit and set its line and marker color
    H1D_MpiMass_0_100MeV->Fit("BW_rho_omega_photoproduction_param", "ERS0", "", 0.6, 0.9);

    // Perform the second fit and set its line and marker color
    // fitParams = (A, B, MRho, Gammarho)
    double BW_rho_omega_photoproduction_fit[8];
    BW_rho_omega_photoproduction_param->GetParameters(BW_rho_omega_photoproduction_fit);

    // double fit_Soding_Errors[7]; for (int i = 0; i < 4; i++) fit_Soding_Errors[i] = Soding_Param->GetParError(i);

    TF1 *BW_rho_omega_fit = new TF1("BW_rho_omega_fit", BW_rho_omega_photoproduction, 0.44, 1.3, 4);
    BW_rho_omega_fit->SetLineColor(kBlack);
    BW_rho_omega_fit->SetLineStyle(1);
    BW_rho_omega_fit->SetParameters(BW_rho_omega_photoproduction_fit[0], BW_rho_omega_photoproduction_fit[1],
                                    BW_rho_omega_photoproduction_fit[2], BW_rho_omega_photoproduction_fit[3],
                                    BW_rho_omega_photoproduction_fit[4], BW_rho_omega_photoproduction_fit[5],
                                    BW_rho_omega_photoproduction_fit[6], BW_rho_omega_photoproduction_fit[7]);

    TF1* fit_phi = new TF1("fit_phi", Interference_phi_fit, 0, 10, 1); // One parameter (A)


    // Optionally, set some initial parameter values to guide the fit. Initial guess for A (amplitude)
    fit_phi->SetParameter(0, 1.0);

    // Perform the fit
    H1D_phi_0_60MeV->Fit("fit_phi", "ERS", "", -3.14, 3.14);
    // fitParams = (A, B, MRho, Gammarho)
    double fit_phi_params; fit_phi->GetParameter(fit_phi_params);
    TF1 *Interference_phi = new TF1("Interference_phi", Interference_phi_fit, -3.14, 3.14, 1);
    Interference_phi->SetParameters(fit_phi_params);
    Interference_phi->SetLineColor(kBlue);
    Interference_phi->SetLineStyle(2);

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

    H1D_MpiMass_0_100MeV->Draw();

    Soding_Param_all->Draw("Same");

    Soding_Param_term1->Draw("Same");

    Soding_Param_term2->Draw("Same");

    TLegend *legend = new TLegend(0.65, 0.6, 0.9, 0.9);
    legend->AddEntry(H1D_MpiMass_0_100MeV, "Data points", "p");
    legend->AddEntry(Soding_Param_all, "Full Sodding", "l");
    AddEmptyLegendEntry(legend, "Parametrization");
    legend->AddEntry(Soding_Param_term1, "Pole term", "l");
    legend->AddEntry(Soding_Param_term2, "Interference term", "l");
    legend->SetTextSize(0.035);
    // Remove the legend box
    legend->SetBorderSize(0);
    legend->Draw();

    c1->SetLeftMargin(0.15); 
    c1->SaveAs("Images/MpiMass_0_100MeV.pdf");
    c1->SaveAs("Images/MpiMass_0_100MeV.png");
    delete c1;

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);

    H1D_phi_0_60MeV->Draw();
    Interference_phi->Draw("Same");
    c2->SaveAs("Images/phi_interference.pdf");
    c2->SaveAs("Images/phi_interference.png");
    delete c2;

    // Canvas
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);

    H1D_MpiMass_0_100MeV->Draw();

    BW_rho_omega_fit->Draw("Same");
    c3->SaveAs("Images/Rho_omega_photoproduction.pdf");
    c3->SaveAs("Images/Rho_omega_photoproduction.png");
    delete c3;


}


void AddEmptyLegendEntry(TLegend *legend, const char *text = "", int color = 0) {
    // Create a "dummy" histogram with no contents and assign it a color
    TH1F *dummy = new TH1F(Form("dummy%d", color), "", 0, 0, 0);
    dummy->SetLineColor(color);

    // Add the dummy entry to the legend with the specified text
    legend->AddEntry(dummy, text, "l");

    // Make sure the dummy histogram doesn't appear in the plot
    dummy->SetFillColor(0);
    dummy->SetLineColor(0);
    dummy->SetLineWidth(0);
    dummy->SetMarkerSize(0);
}