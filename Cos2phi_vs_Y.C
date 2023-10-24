
#include <TCanvas.h>
#include <TH1F.h>
#include <TAxis.h>
#include "utils/utils.h"

void Cos2phi_vs_Y() {

    TFile * InputFile = new TFile( "UpcOutput_Rho.root", "READ" );

    TH2D* H2D_2Cos2phi_vs_Y = (TH2D*) InputFile->Get("hCos2phivsY");
    // Set the y-axis title of the original 2D histogram with LaTeX-like symbols

    // Create a TProfile to store the 1D projection of the mean value of Y for each bin in X
    //
    TProfile* profileX = H2D_2Cos2phi_vs_Y->ProfileX();

    // Use TLatex to manually add LaTeX-like symbols to the title
    // profileX->GetYaxis()->SetTitle("");
    profileX->GetYaxis()->SetTitle("<2cos(2#phi)>");
    profileX->SetMarkerStyle(7);
    profileX->GetYaxis()->SetTitleOffset(1.2);
    profileX->GetYaxis()->SetTitleSize(0.05);
    profileX->GetXaxis()->SetTitle("y");
    profileX->GetXaxis()->SetTitleSize(0.05);
    profileX->GetXaxis()->SetTitleOffset(0.8);

    TCanvas* c1 = new TCanvas("c1", "Canvas", 800, 600);
    
    // Call the function to set the canvas style
    SetCanvasStyle(c1);

    c1->cd();
    profileX -> Draw("E1");

    c1->SetLeftMargin(0.15); 
    c1->SaveAs("Images/Cos2phi_vs_Y.pdf");
    c1->SaveAs("Images/Cos2phi_vs_Y.png");

    // Don't forget to delete the canvas when you're done
    delete c1;

}