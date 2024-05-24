#include <complex>
#include <cmath>
#include "TLatex.h"
#include <TMath.h>
#include <locale>

void vertical_line(double x, int color = kRed);

void Invariant_mass(double ptLow=0.0, double ptHigh=0.099){

    TFile * InputFile = new TFile(
    "/Users/daniel/Documents/Research/BNL/STAR_UPC/output_files/UPC_output.root",
    "READ");

    TFile * InputFile_v2 = new TFile(
    "/Users/daniel/Documents/Research/BNL/STAR_UPC/output_files/UpcOutput_test.root",
    "READ");

    TH2F *H2D_MassVsPt = (TH2F*)InputFile->Get("hMassVsPt");
    TH1F *H1D_MassVsPt_BL = (TH1F*)InputFile_v2->Get("hMpiMass_0_100MeV");

    // Projecting a Pt range from the 2D histogram
    // Parameters: name of the projection, first bin, last bin, option
    int binLow = H2D_MassVsPt->GetXaxis()->FindBin(ptLow);
    int binUp  = H2D_MassVsPt->GetXaxis()->FindBin(ptHigh);
    TH1D *H1D_Mass = H2D_MassVsPt ->ProjectionY("H1D_Mass", binLow, binUp, "e"); // "e" option is for error calculation

    std::ostringstream titleStream;
    titleStream << "#pi^{+}#pi^{-} mass distribution (" << ptLow << " < P_{T}(#pi^{+}#pi^{-}) < " << ptHigh << " GeV/c)";
    H1D_Mass->SetTitle("");
    H1D_Mass->GetYaxis()->SetTitle("dN/dM_{#pi#pi} [1/(GeV/c^{2})]");
    H1D_Mass->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
    // H1D_Mass->SetMarkerStyle(31);
    H1D_Mass->GetYaxis()->SetTitleOffset(0.6);
    H1D_Mass->GetYaxis()->SetTitleSize(0.08);
    H1D_Mass->GetXaxis()->SetTitleSize(0.08);
    H1D_Mass->GetXaxis()->SetTitleOffset(0.6);

    // H1D_Mass->Scale(1/1000.0);

    // H1D_Mass->Rebin(2);
    // H1D_Mass->Scale(2);

    // Hide the statistics box
    H1D_Mass ->SetStats(0);
    H1D_Mass->GetXaxis()->SetRangeUser(0.4, 1.3);

    // Setting the maximum number of digits to force scientific notation
    TGaxis::SetMaxDigits(3);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Histograms", 800, 600);

    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.12);
    
    H1D_Mass->SetLineColor(kBlue);
    H1D_Mass->Draw("E1");
    H1D_MassVsPt_BL->SetMarkerStyle(32);
    H1D_MassVsPt_BL->Draw("SAME P");

    double data_points = H1D_MassVsPt_BL->GetEntries();
    double new_data_points = H1D_Mass->GetEntries();

    cout << "ratio: " << new_data_points/data_points << endl;

    vertical_line(0.770);

    // pt->Draw();

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.9);
    leg->AddEntry(H1D_Mass,  "UPC Data");
    leg->AddEntry(H1D_MassVsPt_BL, "UPC Prev Code");
    // leg->AddEntry(bw_only, "BW only");
    // leg->AddEntry(poly_only, "BKG only");
    // leg->SetTextSize(0.08);
    leg->Draw();

    // Double_t mass = bw->GetParameter("Mass");
    // Double_t width = bw->GetParameter("Width");
    // cout << "Fitted Mass: " << mass << " GeV/c^2" << endl;
    // cout << "Fitted Width: " << width << " GeV/c^2" << endl;

    // Save the canvas
    std::ostringstream Image_name;
    Image_name << "Images/Mpipi_distribution_pt_" << ptLow << "_" << ptHigh << "_GeV.pdf";
    c1->SaveAs(Image_name.str().c_str());
}

void vertical_line(double x, int color = kBlack) {

    double y_min = 0.0;
    double y_max = 24.;

    // Create and draw a line at x with given color
    TLine *line = new TLine(x, y_min, x, y_max);
    line->SetLineColor(color);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("SAME");
}