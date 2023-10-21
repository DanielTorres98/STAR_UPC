void open_UpcOutput() {

    TFile *_file0 = TFile::Open("UpcOutput.root");
    TH1F *hdEdxVsP = (TTree*)input->Get("hdEdxVsP;1");
    hdEdxVsP->Draw("colz");

}