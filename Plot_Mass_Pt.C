#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TColor.h>

void Plot_Mass_Pt() {
    // Open the ROOT file containing the histograms
    TFile *file = TFile::Open("UpcOutput.root");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    // Create a TCanvas to hold the plots
    TCanvas *canvas = new TCanvas("canvas", "Gradient Plots", 800, 600);

    // Define custom colors for the gradient
    const Int_t numColors = 10;
    Int_t colors[numColors] = {
        kBlue, kAzure - 4, kCyan - 6, kTeal - 8, kGreen + 2,
        kGreen + 3, kYellow - 7, kOrange - 3, kRed - 4, kRed
    };

    // Array of histogram names
    const char *histogramNames[] = {
        "hRhoMass_0_20", "hRhoMass_20_40", "hRhoMass_40_60", "hRhoMass_60_80",
        "hRhoMass_80_100", "hRhoMass_100_150", "hRhoMass_150_200", "hRhoMass_200_1000"
    };

    // Create and draw histograms using the custom gradient colors
    for (Int_t i = 0; i < 8; ++i) {
        TH1F *histogram = dynamic_cast<TH1F*>(file->Get(histogramNames[i]));

        if (!histogram) {
            std::cerr << "Error: Could not retrieve histogram " << histogramNames[i] << std::endl;
            continue;
        }

        histogram->SetLineColor(colors[i % numColors]);  // Set custom color
        if (i == 0) {
            histogram->Draw();
        } else {
            histogram->Draw("same");
        }
    }

    // Update the canvas
    canvas->Update();

    // Close the ROOT file
    file->Close();
}