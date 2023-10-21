void pions()
{
    TCanvas *c1 = new TCanvas();

    // TFile *input = new TFile('FemtoDst_Run10AuAu_tofMult_allTracks.root', 'read');
    TFile *input = new TFile("./Data/run11upc_pions.femtodst.root", "read");

    TTree *tree = (TTree*)input->Get("mEvent;1");

    float x, y;

    tree->SetBranchAddress("mvx", &x);
    tree->SetBranchAddress("mvy", &y);

    // takes the number of entries in the tree
    int entries = tree->GetEntries();

    TH1F *histx = new TH1F("hist", "Transverse Momentum", 100, -2, 2);
    // TH1F *histy = new TH1F("hist", "Histogram", 100, -2, 2);
    for(int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        auto pt = pow(x,2)+pow(y,2);
        histx->Fill(pt);
        // histy->Fill(y);
    }

    histx->GetXaxis()->SetTitle("P_T (GeV)");
    histx->GetYaxis()->SetTitle("Counts (N)");

    // Double_t factor = 1.;
    // histx->Scale(factor/entries);
    histx->Draw("hist");

    // histy->Draw("SAME");
    // input->Close();
}