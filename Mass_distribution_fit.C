// Code to generate Mass distribution fits from the Analysis output file.

#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770
#define RhoLifeTime        0.148
#define OmegaPdgMass       0.782
#define OmegaLifeTime      0.017
#define speedoflight       29.9792

#include <complex>
#include <cmath>
#include "TLatex.h"
#include <TMath.h>

#include "utils/Fit_functions.h"
#include "utils/utils.h"

void Mass_distribution_fit(){

    TFile * InputFile = new TFile( "UpcOutput_Rho.root", "READ" );

    // TH1D* H1D_MpiMass_0_100MeV = (TH1D*) InputFile->Get("hRhoMass");
    TH1D* H1D_MpiMass_0_100MeV = (TH1D*) InputFile->Get("hMpiMass_0_100MeV");
    H1D_MpiMass_0_100MeV->SetTitle("#pi^{+}#pi^{-} mass distribution p_{T} < 100 MeV/c^2");
    H1D_MpiMass_0_100MeV->Scale(1/1074.60);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitle("d#sigma/dM_{#pi#pi} [#mub/(GeV/c^{2})]");
    H1D_MpiMass_0_100MeV->SetMarkerStyle(7);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleOffset(1.2);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleOffset(0.8);

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
    Soding_Param_all->SetLineColor(kBlack); Soding_Param_all->SetLineStyle(1);
    Soding_Param_all->SetParameters(fit_Soding_Params);

    TF1 *Soding_Param_term1 = new TF1("Soding_Param_term1", Pole_term_squared, 0.44, 1.3, 4);
    Soding_Param_term1->SetLineColor(kBlack); Soding_Param_term1->SetLineStyle(2);
    Soding_Param_term1->SetParameter(0, fit_Soding_Params[0]);
    Soding_Param_term1->SetParameter(1, fit_Soding_Params[2]);
    Soding_Param_term1->SetParameter(2, fit_Soding_Params[3]);

    TF1 *Soding_Param_term2 = new TF1("SodingEqn_addition_term", Interference_term, 0.44, 1.3, 4);
    Soding_Param_term2->SetLineColor(kBlack); Soding_Param_term2->SetLineStyle(3);
    Soding_Param_term2->SetParameters(fit_Soding_Params);

    // Fit
    TF1 *BW_rho_omega_photoproduction_param = new TF1("BW_rho_omega_photoproduction_param",
                            BW_rho_omega_photoproduction, RhoPdgMass-2*RhoLifeTime, RhoPdgMass+2*RhoLifeTime, 10);

    // Set parameter names
    BW_rho_omega_photoproduction_param->SetParNames("M_#rho", "#Gamma_#rho", "A_#rho", "B_{#pi#pi}", "C_#omega",
                                                  "M_#omega", "#Gamma_#omega", "#phi_#omega", "p0", "p1");

    // Set initial parameter values for A, B RhoMass, GammaRho respectively.
    BW_rho_omega_photoproduction_param->SetParameters(RhoPdgMass, RhoLifeTime, 0.8, -1.21, 0.5, OmegaPdgMass, OmegaLifeTime,
                                                     1.46, 1.0, -1.0);

    BW_rho_omega_photoproduction_param->SetParLimits(0, 0.750, 0.790); // M_rho
    BW_rho_omega_photoproduction_param->SetParLimits(1, 0.1, 0.2);     // Gamma_rho
    BW_rho_omega_photoproduction_param->SetParLimits(2, 0.0, 5.0);     // A
    BW_rho_omega_photoproduction_param->SetParLimits(3, -10.0, 0.0);   // B
    BW_rho_omega_photoproduction_param->SetParLimits(4, 0.0, 1.0);     // C
    BW_rho_omega_photoproduction_param->SetParLimits(5, 0.770, 0.790); // M_omega
    BW_rho_omega_photoproduction_param->SetParLimits(6, 0.01, 0.1);    // Gamma_omega
    BW_rho_omega_photoproduction_param->SetParLimits(7, -3.14, 3.14);  // phi
    BW_rho_omega_photoproduction_param->SetParLimits(8, 0.0, 3.0);     // p0
    BW_rho_omega_photoproduction_param->SetParLimits(9, -2.0, 0.0);    // p1


    // Perform the first fit and set its line and marker color
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Genetic");
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.000001);
    H1D_MpiMass_0_100MeV->Fit("BW_rho_omega_photoproduction_param", "ERS0", "", 0.6, 1.3);
    double chi2 = BW_rho_omega_photoproduction_param->GetChisquare();
    int ndf = BW_rho_omega_photoproduction_param->GetNDF();
    double chi2PerNdf = chi2 / ndf;

    // Perform the second fit and set its line and marker color
    // fitParams = (A, B, MRho, Gammarho)
    double BW_rho_omega_photoproduction_param_values[8];
    BW_rho_omega_photoproduction_param->GetParameters(BW_rho_omega_photoproduction_param_values);
    double M_rho = BW_rho_omega_photoproduction_param_values[0];
    double Gamma_rho = BW_rho_omega_photoproduction_param_values[1];
    double A_rho = BW_rho_omega_photoproduction_param_values[2];
    double B_pipi = BW_rho_omega_photoproduction_param_values[3];
    double C_omega = BW_rho_omega_photoproduction_param_values[4];
    double M_omega = BW_rho_omega_photoproduction_param_values[5];
    double Gamma_omega = BW_rho_omega_photoproduction_param_values[6];
    double phi_omega = BW_rho_omega_photoproduction_param_values[7];

    TF1 *BW_rho_omega_photoproduction_fit = new TF1("BW_rho_omega_photoproduction_fit", BW_rho_omega_photoproduction, 0.6, 1.3, 10);
    BW_rho_omega_photoproduction_fit->SetLineColor(kBlack);
    BW_rho_omega_photoproduction_fit->SetLineStyle(1);
    BW_rho_omega_photoproduction_fit->SetParameters(BW_rho_omega_photoproduction_param_values);

    TF1 *BW_rho_omega_photoproduction_rho_pole_term = new TF1("BW_rho_omega_photoproduction_rho_pole_term", 
                                                              Pole_term_squared, 0.6, 1.3, 3);
    BW_rho_omega_photoproduction_rho_pole_term->SetLineColor(kRed);
    BW_rho_omega_photoproduction_rho_pole_term->SetLineStyle(1);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(0, A_rho);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(1, M_rho);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(2, Gamma_rho);


    double Branching_ratio_omega_2_pipi = 0.0153;
    TF1 *BW_rho_omega_photoproduction_omega_pole_term = new TF1("BW_rho_omega_photoproduction_omega_pole_term", 
                                                              Pole_term_squared, 0.6, 1.3, 3);
    BW_rho_omega_photoproduction_omega_pole_term->SetLineColor(kCyan);
    BW_rho_omega_photoproduction_omega_pole_term->SetLineStyle(1);
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(0, C_omega*sqrt(Branching_ratio_omega_2_pipi));
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(1, M_omega);
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(2, Gamma_rho);

    TF1 *B_pipi_line = new TF1("B_pipi_line", "[0]", 0.6, 1.3);

    B_pipi_line->SetParameter(0, B_pipi*B_pipi);
    // Set line properties if needed (e.g., color, width...)ß
    B_pipi_line->SetLineColor(kBlue);   // Set color of the line to blue
    B_pipi_line->SetLineWidth(2);      // Set line width

    // Multiply
    TF1 *rho_interference_term = new TF1("rho_interference_term", Interference_term, 0.6, 1.3, 4);
    rho_interference_term->SetLineStyle(1);
    rho_interference_term->SetParameter(0, A_rho);
    rho_interference_term->SetParameter(1, B_pipi);
    rho_interference_term->SetParameter(2, M_rho);
    rho_interference_term->SetParameter(3, Gamma_rho);
    rho_interference_term->SetLineColor(kGreen);

    TF1 *omega_interference_term = new TF1("omega_interference_term", Interference_term, 0.6, 1.3, 4);
    omega_interference_term->SetLineStyle(2);
    omega_interference_term->SetParameter(0, C_omega*sqrt(Branching_ratio_omega_2_pipi));
    omega_interference_term->SetParameter(1, B_pipi);
    omega_interference_term->SetParameter(2, M_omega);
    omega_interference_term->SetParameter(3, Gamma_omega);
    omega_interference_term->SetLineColor(kSpring);

    TF1 *BW_interference_term = new TF1("BW_interference_term", BW_Interference_term, 0.6, 1.3, 7);
    BW_interference_term->SetLineColor(kMagenta);
    BW_interference_term->SetLineStyle(2);
    BW_interference_term->SetParameter(0, A_rho);
    BW_interference_term->SetParameter(1, M_rho);
    BW_interference_term->SetParameter(2, Gamma_rho);
    BW_interference_term->SetParameter(3, C_omega*sqrt(Branching_ratio_omega_2_pipi));
    BW_interference_term->SetParameter(4, M_omega);
    BW_interference_term->SetParameter(5, Gamma_omega);
    BW_interference_term->SetParameter(6, phi_omega);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    SetCanvasStyle(c1);
    c1->cd();

    H1D_MpiMass_0_100MeV->Draw("E1");

    Soding_Param_all->Draw("Same");

    Soding_Param_term1->Draw("Same");

    Soding_Param_term2->Draw("Same");

    TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->AddEntry(H1D_MpiMass_0_100MeV, "Data points", "p");
    legend->AddEntry(Soding_Param_all, "Modified S\\ddot{o}ding", "l");
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

    // Canvas
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    SetCanvasStyle(c3);
    H1D_MpiMass_0_100MeV->Draw("E1");
    double hist_maxima = H1D_MpiMass_0_100MeV->GetMaximum();
    H1D_MpiMass_0_100MeV->GetYaxis()->SetRangeUser(-1.0, hist_maxima*1.2);
    // gStyle->SetOptFit(1111);
    BW_rho_omega_photoproduction_fit->Draw("Same");

    BW_rho_omega_photoproduction_rho_pole_term->Draw("Same");

    BW_rho_omega_photoproduction_omega_pole_term->Draw("Same");

    B_pipi_line->Draw("Same");

    rho_interference_term->Draw("Same");

    omega_interference_term->Draw("Same");

    BW_interference_term->Draw("Same");

    TLegend *legend2 = new TLegend(0.6, 0.5, 0.88, 0.8);
    legend2->AddEntry(H1D_MpiMass_0_100MeV, "Data points", "p");
    legend2->AddEntry(BW_rho_omega_photoproduction_fit, "#rho-#omega full fit", "l");
    legend2->AddEntry(BW_rho_omega_photoproduction_rho_pole_term, "#rho BW fit", "l");
    legend2->AddEntry(BW_rho_omega_photoproduction_omega_pole_term, "#omega BW fit", "l");
    legend2->AddEntry(B_pipi_line, "B_{#pi#pi}", "l");
    legend2->AddEntry(rho_interference_term, "#rho interference #pi#pi", "l");
    legend2->AddEntry(omega_interference_term, "#omega interference #pi#pi", "l");
    legend2->AddEntry(BW_interference_term, "#omega-#rho interference", "l");

    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.035);
    legend2->Draw();

    TPaveText *pt = new TPaveText(0.6, 0.3, 0.85, 0.4, "NDC"); // NDC sets coordinates relative to pad dimensions
    pt->AddText(Form("#chi^{2}/DOF = %.2f", chi2PerNdf));
    pt->SetFillColor(0); // Transparent background
    pt->SetTextSize(0.05);
    pt->Draw();
    c3->SetLeftMargin(0.15);
    c3->SaveAs("Images/Rho_omega_photoproduction.pdf");
    c3->SaveAs("Images/Rho_omega_photoproduction.png");
    delete c3;
}
