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

void Mass_distribution_fit_xy_scan(){

   const int N = 52;
   double x[N] = {0.6036217304,
                    0.6084507042,
                    0.63138833,
                    0.6410462777,
                    0.6434607646,
                    0.6555331992,
                    0.6639839034,
                    0.6688128773,
                    0.6784708249,
                    0.6881287726,
                    0.6929577465,
                    0.6965794769,
                    0.7062374245,
                    0.7110663984,
                    0.723138833,
                    0.7315895372,
                    0.7400402414,
                    0.7460764588,
                    0.7533199195,
                    0.7557344064,
                    0.7617706237,
                    0.7665995976,
                    0.7690140845,
                    0.7714285714,
                    0.7762575453,
                    0.7786720322,
                    0.7822937626,
                    0.783501006,
                    0.785915493,
                    0.7883299799,
                    0.7919517103,
                    0.7943661972,
                    0.8040241449,
                    0.8064386318,
                    0.8088531187,
                    0.8136820926,
                    0.8197183099,
                    0.8269617706,
                    0.8293762575,
                    0.8438631791,
                    0.8559356137,
                    0.8668008048,
                    0.8716297787,
                    0.890945674,
                    0.909054326,
                    0.9718309859,
                    1.010462777,
                    1.050301811,
                    1.091348089,
                    1.13722334,
                    1.171026157,
                    1.278470825}; 
    double y[N] = {
                    6.245847176,
                    6.710963455,
                    7.574750831,
                    8.438538206,
                    8.837209302,
                    9.568106312,
                    9.966777409,
                    10.7641196,
                    11.29568106,
                    12.62458472,
                    13.48837209,
                    14.08637874,
                    14.88372093,
                    15.94684385,
                    17.6744186,
                    18.87043189,
                    19.46843854,
                    19.86710963,
                    19.53488372,
                    20,
                    20.06644518,
                    19.80066445,
                    19.40199336,
                    19.6013289,
                    18.87043189,
                    18.00664452,
                    16.41196013,
                    14.88372093,
                    14.01993355,
                    13.15614618,
                    12.62458472,
                    12.02657807,
                    11.1627907,
                    10.56478405,
                    10.29900332,
                    9.833887043,
                    8.77076412,
                    7.707641196,
                    7.043189369,
                    5.714285714,
                    4.850498339,
                    4.053156146,
                    3.720930233,
                    2.657807309,
                    2.126245847,
                    0.9966777409,
                    0.6644518272,
                    0.5315614618,
                    0.3986710963,
                    0.1993355482,
                    0.1993355482,
                    0.06644518272};

    TGraph *H1D_MpiMass_0_100MeV = new TGraph(N, x, y);

    H1D_MpiMass_0_100MeV->SetTitle("#pi^{+}#pi^{-} mass distribution p_{T} < 100 MeV/c^2");
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^s{2})");
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitle("d#sigma/dM_{#pi#pi} [#mub/(GeV/c^{2})]");
    H1D_MpiMass_0_100MeV->SetMarkerStyle(7);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleOffset(1.2);
    H1D_MpiMass_0_100MeV->GetYaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetRangeUser(0.5, 1.3);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleSize(0.05);
    H1D_MpiMass_0_100MeV->GetXaxis()->SetTitleOffset(0.8);

    // Fit
    TF1 *BW_rho_omega_photoproduction_param = new TF1("BW_rho_omega_photoproduction_param",
                            BW_rho_omega_photoproduction, RhoPdgMass-2*RhoLifeTime, RhoPdgMass+2*RhoLifeTime, 8);

    // Set parameter names
    BW_rho_omega_photoproduction_param->SetParNames("M_#rho", "#Gamma_#rho", "A_#rho", "B_{#pi#pi}", "C_#omega",
                                                  "M_#omega", "#Gamma_#omega", "#phi_#omega");

    // Set initial parameter values for A, B RhoMass, GammaRho respectively.
    BW_rho_omega_photoproduction_param->SetParameters(RhoPdgMass, RhoLifeTime, 1.5, -1.0, 0.55, OmegaPdgMass, OmegaLifeTime,
                                                     1.0);

    BW_rho_omega_photoproduction_param->SetParLimits(0, 0.6, 1.0);
    BW_rho_omega_photoproduction_param->SetParLimits(1, 0.0, 0.4);
    BW_rho_omega_photoproduction_param->SetParLimits(2, 0.0, 5.0); // A
    BW_rho_omega_photoproduction_param->SetParLimits(3, -10, 0); // B
    BW_rho_omega_photoproduction_param->SetParLimits(4, 0.0, 5.0); // C
    BW_rho_omega_photoproduction_param->SetParLimits(5, 0.6, 1.0);
    BW_rho_omega_photoproduction_param->SetParLimits(6, 0.0, 0.4);
    BW_rho_omega_photoproduction_param->SetParLimits(7, -3.14, 3.14);


    // Perform the first fit and set its line and marker color
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Genetic");
    H1D_MpiMass_0_100MeV->Fit("BW_rho_omega_photoproduction_param", "ERS0", "", 0.6, 1.1);
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

    TF1 *BW_rho_omega_photoproduction_fit = new TF1("BW_rho_omega_photoproduction_fit", BW_rho_omega_photoproduction, 0.44, 1.3, 8);
    BW_rho_omega_photoproduction_fit->SetLineColor(kBlack);
    BW_rho_omega_photoproduction_fit->SetLineStyle(1);
    BW_rho_omega_photoproduction_fit->SetParameters(BW_rho_omega_photoproduction_param_values);

    TF1 *BW_rho_omega_photoproduction_rho_pole_term = new TF1("BW_rho_omega_photoproduction_rho_pole_term", 
                                                              Pole_term_squared, 0.44, 1.3, 3);
    BW_rho_omega_photoproduction_rho_pole_term->SetLineColor(kRed);
    BW_rho_omega_photoproduction_rho_pole_term->SetLineStyle(1);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(0, A_rho);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(1, M_rho);
    BW_rho_omega_photoproduction_rho_pole_term->SetParameter(2, Gamma_rho);


    double Branching_ratio_omega_2_pipi = 0.0153;
    TF1 *BW_rho_omega_photoproduction_omega_pole_term = new TF1("BW_rho_omega_photoproduction_omega_pole_term", 
                                                              Pole_term_squared, 0.44, 1.3, 3);
    BW_rho_omega_photoproduction_omega_pole_term->SetLineColor(kCyan);
    BW_rho_omega_photoproduction_omega_pole_term->SetLineStyle(1);
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(0, C_omega*sqrt(Branching_ratio_omega_2_pipi));
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(1, M_omega);
    BW_rho_omega_photoproduction_omega_pole_term->SetParameter(2, Gamma_rho);

    TF1 *B_pipi_line = new TF1("B_pipi_line", "[0]", 0.44, 1.3);

    B_pipi_line->SetParameter(0, B_pipi*B_pipi);
    // Set line properties if needed (e.g., color, width...)ß
    B_pipi_line->SetLineColor(kBlue);   // Set color of the line to blue
    B_pipi_line->SetLineWidth(2);      // Set line width

    // Multiply
    TF1 *rho_interference_term = new TF1("rho_interference_term", Interference_term, 0.44, 1.3, 4);
    rho_interference_term->SetLineStyle(1);
    rho_interference_term->SetParameter(0, A_rho);
    rho_interference_term->SetParameter(1, B_pipi);
    rho_interference_term->SetParameter(2, M_rho);
    rho_interference_term->SetParameter(3, Gamma_rho);
    rho_interference_term->SetLineColor(kGreen);

    TF1 *omega_interference_term = new TF1("omega_interference_term", Interference_term, 0.44, 1.3, 4);
    omega_interference_term->SetLineStyle(2);
    omega_interference_term->SetParameter(0, C_omega*sqrt(Branching_ratio_omega_2_pipi));
    omega_interference_term->SetParameter(1, B_pipi);
    omega_interference_term->SetParameter(2, M_omega);
    omega_interference_term->SetParameter(3, Gamma_omega);
    omega_interference_term->SetLineColor(kSpring);

    TF1 *BW_interference_term = new TF1("BW_interference_term", BW_Interference_term, 0.44, 1.3, 7);
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
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    SetCanvasStyle(c3);
    H1D_MpiMass_0_100MeV->Draw();
    double hist_maxima = H1D_MpiMass_0_100MeV->GetMaximum();
    cout << hist_maxima << endl;
    H1D_MpiMass_0_100MeV->GetYaxis()->SetRangeUser(-5, 30);
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
    c3->SaveAs("Images/Rho_omega_photoproduction_xy_scan.pdf");
    c3->SaveAs("Images/Rho_omega_photoproduction_xy_scan.png");
    delete c3;
}
