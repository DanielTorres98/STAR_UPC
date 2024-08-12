#include <complex>
#include <cmath>
#include <cstring>  
#include "TLatex.h"
#include <TMath.h>
#include <locale>

void normalizeHistogram(TH1* histogram);

double ComputeChi2(TH1 *hist, TF1 *func);

void Phi_Distribution(const char* physics="pion"){

    float ptlow = 0.0;
    float pthigh = 0.10;
    float masslow = 0.65;
    float masshigh = 0.9;
    // float ptlow = 0.5;
    // float pthigh = 0.6;
    // float masslow = 0.49;
    // float masshigh = 0.505;

    TFile * InputFile = new TFile(
     "/Users/daniel/Research/STAR_UPC/output_files/UPC_output.root",
    "READ");

    TH3F *H3D_PhiVsPtVsMass = (TH3F*)InputFile->Get("H3D_PhiVsPtVsMass");

    int bin_pt_low = H3D_PhiVsPtVsMass->GetYaxis()->FindBin(ptlow);
    int bin_pt_high = H3D_PhiVsPtVsMass->GetYaxis()->FindBin(pthigh);
    int bin_mass_low = H3D_PhiVsPtVsMass->GetZaxis()->FindBin(masslow);
    int bin_mass_high = H3D_PhiVsPtVsMass->GetZaxis()->FindBin(masshigh);
    TH1D* H1D_Phi = H3D_PhiVsPtVsMass->ProjectionX("H1D_Phi", bin_pt_low, bin_pt_high, bin_mass_low, bin_mass_high);

    std::ostringstream titleStream;
    titleStream << ptlow << " < P_{T} (GeV/c) < " << pthigh << ", "
                << masslow << " < M_{#pi#pi} (GeV/c^2) < " << masshigh;
    H1D_Phi->SetTitle("");
    H1D_Phi->GetYaxis()->SetTitle("US - LS (norm. to unity)");
    H1D_Phi->GetXaxis()->SetTitle("#phi([#pi^{+}+#pi^{+}], [#pi^{+}-#pi^{+}])");
    H1D_Phi->SetMarkerStyle(31);
    H1D_Phi->GetYaxis()->SetTitleSize(0.08);
    H1D_Phi->GetYaxis()->SetTitleOffset(0.6);
    H1D_Phi->GetXaxis()->SetTitleSize(0.08);
    H1D_Phi->GetXaxis()->SetTitleOffset(0.6);
    H1D_Phi->GetYaxis()->SetRangeUser(0.0, 2.);

    H1D_Phi->Rebin(2);
    normalizeHistogram(H1D_Phi);

    // Define the fitting function
    TF1 *fitFunc = new TF1("fitFunc", "1+[0]*cos(2*x)", -TMath::Pi(), TMath::Pi());
    H1D_Phi->Fit(fitFunc, "R"); // "R" option for fitting within the function range

    double chi2_me = ComputeChi2(H1D_Phi, fitFunc);
    // Extract the fit parameters and Chi^2/DOF
    double a = fitFunc->GetParameter(0);
    double aError = fitFunc->GetParError(0);
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double chi2ndf = (ndf != 0) ? chi2 / ndf : 0;  // Check to avoid division by zero

    cout << "chi2 = " << chi2_me << endl;
    cout << "ndf" << ndf << endl;
    cout << endl;
    // Display fit results and Chi^2/DOF using TPaveText
    TPaveText *pt = new TPaveText(0.2, 0.2, 0.7, 0.4, "NDC"); // Adjusted coordinates for bottom right
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001); // Solid fill
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12); // Left-aligned text
    pt->SetTextFont(42); // Helvetica
    pt->SetTextSize(0.06);
    pt->AddText("Fit Results:");
    pt->AddText(Form("A = %f #pm %f", a, aError));
    pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf)); 


    H1D_Phi->SetStats(0);
    TCanvas *c1 = new TCanvas("c1", "Histograms", 800, 600);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.12);

    H1D_Phi->SetLineColor(kBlack);
    H1D_Phi->Draw();
    fitFunc->Draw("SAME");
    pt->Draw();

    // Add a legend
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(H1D_Phi,  "UPC data");
    leg->AddEntry(fitFunc, "1+ACos(2#phi)");
    leg->SetTextSize(0.06); 
    leg->Draw();
    c1->SaveAs("Images/Phi_Distribution.pdf");
}

void normalizeHistogram(TH1* histogram) {
    /**
     * @brief Normalize histagram so that the average is the unity. It also scale the error
     * accordingly.
     *
     * @param histogram: histogram to normalize.
     */
    // Calculate the sum of bin contents and the number of bins
    double sum = 0.0;
    int numBins = histogram->GetNbinsX();

    for (int i = 1; i <= numBins; ++i) {
        sum += histogram->GetBinContent(i);
    }

    // Calculate the average value of the histogram
    double average = sum / numBins;

    // Normalize the histogram by dividing each bin content by the average value
    for (int i = 1; i <= numBins; ++i) {
        double binContent = histogram->GetBinContent(i);
        histogram->SetBinContent(i, binContent / average);
        float normalizedError = histogram->GetBinError(i) / average;
        histogram->SetBinError(i, normalizedError);
    }
}

// Function to compute chi-squared
double ComputeChi2(TH1 *hist, TF1 *func) {
    if (!hist || !func) {
        std::cerr << "Histogram or function not provided." << std::endl;
        return -1;
    }

    double chi2 = 0.0;
    int nBins = hist->GetNbinsX();

    for (int i = 1; i <= nBins; ++i) {
        double observed = hist->GetBinContent(i);
        double error = hist->GetBinError(i);
        double expected = func->Eval(hist->GetBinCenter(i));

        if (error != 0) {
            chi2 += TMath::Power((observed - expected) / error, 2);
        }
    }

    return chi2;
}