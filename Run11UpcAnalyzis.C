#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770
#define c                  29.9792

#include <cmath>
#include <TMath.h>

//I want to define a function. I have to define it before the main function "FindRhos". I'll write
// what this function does *after* FindRhos. I would encourage you to define functions for any
// repeated code, e.g. track cuts
double GetMassSqr(double p, double beta);

double GetDeltaTOF(double s, double p, double mass);

// Function to fit: 1 - A*cos(2*phi)
Double_t myFitFunction(Double_t* x, Double_t* par);

void normalizeHistogram(TH1* histogram);

void normalizeHistogramProbability(TH1F* histogram);

void Run11UpcAnalyzis() {


    TFile *myFile = TFile::Open("./Data/FemtoDst_Run10AuAu_wZDC.root");
    TTreeReader myReader("FemtoDst", myFile);

    // Output file where histograms go. Histograms defined below a TFile are automatically
    // associated with that TFile.
    TFile * OutputFile = new TFile( "UpcOutput.root", "RECREATE" );

    // 2D histograms. Defined by (name, title, #binsx, lower range x, upper range x, binsy, low y,
    // upper y)
    //
    TH2F * hdEdxVsP = new TH2F("hdEdxVsP", "Track energy loss (dE/dx) vs q*momentum (GeV/c)", 
                               1000, -3., 3., 1000, 0., 10.);
    TH2F * hMassVsPt = new TH2F("hMassVsPt", "Invariant Mass vs Transverse momentum", 
                               1000, 0, 0.1, 1000, 0.3, 1);
    hMassVsPt->GetXaxis()->SetTitle("P_{T}"); hMassVsPt->GetYaxis()->SetTitle("M_{#rho}");

    TH2F * hPxPy = new TH2F("hPxPy", "Counts per momentum tranversal direction", 
                               100, -0.1, 0.1, 100, -0.1, 0.1);
    TH2F * hTofVsP = new TH2F("hTofVsP", "TOF m^2 (GeV/c) vs momentum (GeV/c)", 1000, 0, 3.5,
                              1000, 0.0, 1.2);
    TH2F * hXpionXkaon = new TH2F("hXpionXkaon", "#Chi_{#pi#pi}^{2} vs #Chi_{kk}^{2}", 
                                  1000, 0, 1000, 1000, 0, 200);
    TH2F * hXpionXee = new TH2F("hXpionXee", "#Chi_{#pi#pi}^{2} vs #Chi_{ee}^{2}", 
                                1000, 0, 200, 1000, 0, 200);
    TH2F * hXpionXpp = new TH2F("hXpionXpp", "#Chi_{#pi#pi}^{2} vs #Chi_{pp}^{2}", 
                                1000, 0, 1000, 1000, 0, 200);
    int numBinsX = 29;
    Double_t xBins[numBinsX+1];
    double x = 0.0;
    for (Int_t i = 0; i <= numBinsX; ++i) {
        xBins[i] = x; // Bin edges grow quadratically
        if (i<=15) {
            x = x+0.1/15;
        }
        else {
            x = x+0.15/15;
        }
    }
    // Create a 2D histogram with custom bin edges along the x-axis
    TH2F* hCos2phivsPT = new TH2F("hCos2phivsPT", "2cos(2#phi) vs P_{T}",
                            numBinsX, xBins, 200, -2, 2);
    TH2F* hCos2phivsY = new TH2F("hCos2phivsY", "2cos(2#phi) vs y",
                            50, 0, 1, 200, -2, 2);
    
    //1D histogram. Defined by (name, title, #binsx, lower range x, upper range x)
    //
    TH1F * hRhoMass = new TH1F("hRhoMass", "#rho mass distribution", 100, 0.6, 0.9);
    hRhoMass->GetXaxis()->SetTitle("m_{#rho}(GeV)"); hRhoMass->GetYaxis()->SetTitle("Counts");
    
    // Mass plots
    //
    TH1F * hMpiMass_0_100MeV = new TH1F("hMpiMass_0_100MeV", 
                                        "#pi^{+}#pi^{-} mass distribution P_{T} < 100 MeV",
                                        100, 0.46, 1.1);
    hMpiMass_0_100MeV->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV)");
    hMpiMass_0_100MeV->GetYaxis()->SetTitle("dN/dM_{#pi^{+}#pi^{-}}(GeV)");

    TH1F * hKaonMass = new TH1F("hKaonMass", "#pi^{+} #pi^{-} pairs", 100, 0.44, 0.56);
    hKaonMass->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}(GeV)}");
    hKaonMass->GetYaxis()->SetTitle("dN/dM_{#pi^{+}#pi^{-}}(GeV)");

    TH1F * hRhoTP = new TH1F("hRhoTP", "#rho Transverse Momentum", 80, 0., 0.35);
    TH1F * hDDTOF = new TH1F("hDDTOF", "#Delta#Delta TOF", 100, 0, 100);
    TH1F * hDeltaPhi = new TH1F("hDeltaPhi", "#Delta#phi", 80, -3.14, 3.14);
    TH1F * hcosphi = new TH1F("hcosphi", "cos#phi", 80, -1, 1);
    hDeltaPhi->GetXaxis()->SetTitle("#Delta#phi"); hDeltaPhi->GetYaxis()->SetTitle("Counts");
    TH1F * hPhipPPhim = new TH1F("hPhipPPhim", "#phi(#pi^{+} + #pi^{-})", 80, -3.14, 3.14);
    hPhipPPhim->GetXaxis()->SetTitle("#phi(#pi^{+} + #pi^{-})"); hPhipPPhim->GetYaxis()->SetTitle("Counts");
    TH1F * hPhipMPhim = new TH1F("hPhipMPhim", "#phi(#pi^{+} - #pi^{-})", 80, -3.14, 3.14);
    hPhipMPhim->GetXaxis()->SetTitle("#phi(#pi^{+} - #pi^{-})"); hPhipMPhim->GetYaxis()->SetTitle("Counts");
    TH1F * hPhip = new TH1F("hPhip", "#pi^{+}", 80, -3.14, 3.14);
    hPhip->GetXaxis()->SetTitle("#phi"); hPhip->GetYaxis()->SetTitle("Counts");
    TH1F * hPhim = new TH1F("hPhim", "#pi^{-}", 80, -3.14, 3.14);
    TH1F * hphi = new TH1F("hphi", "#phi", 120, -3.14, 3.14);
    hphi->GetXaxis()->SetTitle("#phi"); hphi->GetYaxis()->SetTitle("Counts");
    TH1F * hCos2phi = new TH1F("hCos2phi", "cos2#phi", 80, -1, 1);
    hCos2phi->GetXaxis()->SetTitle("#cos2phi"); hCos2phi->GetYaxis()->SetTitle("Counts");
    hPhim->GetXaxis()->SetTitle("#phi"); hPhim->GetYaxis()->SetTitle("Counts");

    TTreeReaderArray<Float_t> trackPt(myReader, "Tracks.mPt"); //Transverse momentum in GeV
    TTreeReaderArray<Float_t> trackEta(myReader, "Tracks.mEta"); //Pseudorapidity
    // Azimuthal angle in STAR coordinate system
    TTreeReaderArray<Float_t> trackPhi(myReader, "Tracks.mPhi");
    TTreeReaderArray<Char_t> trackQ(myReader, "Tracks.mQ"); //Track charge +/- 1.
    TTreeReaderArray<Float_t> trackdEdx(myReader, "Tracks.mDedx"); //energy lost in TPC
    // Distance of closest approach to primary vertex
    TTreeReaderArray<Float_t> trackDCA(myReader, "Tracks.mDCA");
    TTreeReaderArray<Float_t> trackBeta(myReader, "BTofPidTraits.mBTofBeta"); //Velocity fraction of track as measured by TOF
    TTreeReaderArray<Float_t> trackLength(myReader, "BTofPidTraits.mLength"); 
    TTreeReaderArray<Float_t> trackTOF(myReader, "BTofPidTraits.mBTof"); 

    // Number of standard deviations from dEdx fit for a given particle.
    TTreeReaderArray<Float_t> trackNSigmaPion(myReader, "Tracks.mNSigmaPion"); //Sigma Pion
    TTreeReaderArray<Float_t> trackNSigmaKaon(myReader, "Tracks.mNSigmaKaon"); //Sigma Kaon
    TTreeReaderArray<Float_t> trackNSigmaElectron(myReader, "Tracks.mNSigmaElectron"); //Sigma Electron
    TTreeReaderArray<Float_t> trackNSigmaProton(myReader, "Tracks.mNSigmaProton"); //Sigma Proton

    vector<TLorentzVector> uusiglv; // uu pairs with daughters from k0s mass range
    vector<TLorentzVector> uunotsiglv; // uu pairs with daughters not from k0s mass range


    // Defining Variables
    unsigned int nPos, nNeg;
    unsigned int iEvent = 0; // event loop iterator index
    unsigned int TotTracks; //# of tracks in given event

    int prev_frac = -1;
    int fraction;
    int EventsInFile = myReader.GetEntries();
    cout << "There are " << EventsInFile << " events in this file" << endl;
    // Set smaller if you want a small fraction of the data. If you want all just make this a big
    // number
    //
    int EventLimit = EventsInFile;
    cout << "The event limit is " << EventLimit << " events" << endl; 
    int totalevents = std::min(EventsInFile, EventLimit);

    //TLorentzVector is a 4-vector class (https://root.cern.ch/doc/master/classTLorentzVector.html)
    //I'll declare the + and - track vectors here and then reset the values for every track.
    //
    TLorentzVector pTrackMomentum;
    TLorentzVector nTrackMomentum;
    TLorentzVector TrackMomentum;
    TLorentzVector ParentMomentum;
    TLorentzVector Polarization;

    float pTrackMSqr, nTrackMSqr, TrackMSqr, RhoMsqr, cosphi, sinphi;
    float NSigmaPion1, NSigmaPion2, NSigmaKaon1, NSigmaKaon2, NSigmaKaon, NSigmaPion;
    float NSigmaElectron1, NSigmaElectron2, NSigmaProton1, NSigmaProton2, NSigmaElectron;
    float NSigmaProton, t2, t1, DeltaTOF_expected, DeltaTOF, DeltaDeltaTOF, SurvivingFraction;


    //"while" loop to go over all of the events. This could also be written as a "for" loop
    //
    int number_of_zeros = 0;
    SurvivingFraction = 0;
    
    while (myReader.Next() && iEvent < EventLimit) {
        iEvent ++;
        TotTracks = trackPt.GetSize();

        //get out of loop if the event has anything other than 2 tracks. You will probably want to
        //uncomment this cut later, but if you want to look at more kinds of tracks it might be
        //good to comment out now. 
        if(TotTracks != 2) continue; 


        //little block of code to keep track of % done running code. Can be deleted.
        fraction = iEvent * 100.0 / totalevents;
            if (fraction % 1 == 0 && fraction != prev_frac) {
            printf("\b\b\b%2d%%", fraction);
            std::cout << std::flush;
            prev_frac = fraction;
        }

        nPos = 0;
        nNeg = 0;
        //I'll make nested loops here. So a loop of negative tracks inside of a loop of positive tracks. In principle you only care about events with two tracks (I made that cut earlier), but I want to keep the code general.
        //start first track loop, over + tracks.
        bool found_zero = false;
        int count;
        count = 1;
        for (unsigned int iTp = 0; iTp < TotTracks; iTp++ ){
            //4 momentum is uniquely defined by pT, eta, phi, and mass. We assume the tracks are
            //pions and later evaluate the truth of this idea and correct accordingly.
            TrackMomentum.SetPtEtaPhiM(trackPt[iTp], trackEta[iTp], trackPhi[iTp], PionPdgMass);
            TrackMSqr = GetMassSqr(TrackMomentum.P(), trackBeta[iTp]);
            t1 = GetDeltaTOF(trackLength[iTp], nTrackMomentum.P(), PionPdgMass);
            NSigmaPion1 = trackNSigmaPion[iTp];
            NSigmaKaon1 = trackNSigmaKaon[iTp];
            NSigmaElectron1 = trackNSigmaElectron[iTp];
            NSigmaProton1 = trackNSigmaProton[iTp];
            // cout << "GetMassSqr->" << TrackMSqr << " m^2->" << TrackMomentum.M2() <<endl;
            // cout << "P->" << TrackMomentum.P()<< " beta->" << trackBeta[iTp] <<endl;
            hTofVsP->Fill(TrackMomentum.P(), trackBeta[iTp]);

            if (!found_zero && trackBeta[iTp] == 0) {
                found_zero = true;
                number_of_zeros += 1;
            }

            //this loop is for positive pions, but before I make cuts you can fill general track 
            //information here fill 2d histogram with x value, yvalue
            //it's convention to separate negative and positive particles. There is no deep meaning.
            pTrackMomentum.SetPtEtaPhiM(trackPt[iTp], trackEta[iTp], trackPhi[iTp], PionPdgMass);
            if(trackQ[iTp] > 0) hdEdxVsP->Fill(pTrackMomentum.P(), trackdEdx[iTp]);

            else hdEdxVsP->Fill(-1.*pTrackMomentum.P(), trackdEdx[iTp]);
            //get out of this loop if the charge is < 0, so beyond
            //this all tracks indexed by iTp should be positive
            if(trackQ[iTp] < 0) continue; 

            //call function defined above and detailed below
            pTrackMSqr = GetMassSqr(pTrackMomentum.P(), trackBeta[iTp]); 
            // add cut on nSigma and (later) mass squared here.
            //negative track loop, nested
            for (unsigned int iTn = 0; iTn < TotTracks; iTn++ ){

                if(trackQ[iTn] > 0) continue; //negative track cut

                // 
                NSigmaPion2 = trackNSigmaPion[iTn];
                NSigmaKaon2 = trackNSigmaKaon[iTn];
                NSigmaElectron2 = trackNSigmaElectron[iTn];
                NSigmaProton2 = trackNSigmaProton[iTn];

                NSigmaPion = NSigmaPion1*NSigmaPion1 + NSigmaPion2*NSigmaPion2;
                NSigmaKaon = NSigmaKaon1*NSigmaKaon1 + NSigmaKaon2*NSigmaKaon2;
                NSigmaElectron = NSigmaElectron1*NSigmaElectron1 + NSigmaElectron2*NSigmaElectron2;
                NSigmaProton = NSigmaProton1*NSigmaProton1 + NSigmaProton2*NSigmaProton2;

                nTrackMomentum.SetPtEtaPhiM( trackPt[iTn], trackEta[iTn], trackPhi[iTn], PionPdgMass);

                nTrackMSqr = GetMassSqr(nTrackMomentum.P(), trackBeta[iTn]);
                t2 = GetDeltaTOF(trackLength[iTn], nTrackMomentum.P(), PionPdgMass);

                DeltaTOF = abs(trackTOF[iTn] - trackTOF[iTp]);
                DeltaTOF_expected = abs(t2-t1);

                DeltaDeltaTOF = abs(DeltaTOF - DeltaTOF_expected);

                ParentMomentum = nTrackMomentum + pTrackMomentum;
                Polarization = nTrackMomentum - pTrackMomentum;
                float DeltaPhi = ParentMomentum.Phi() - Polarization.Phi();
                float RhoMass = ParentMomentum.M();
                RhoMsqr = GetMassSqr(ParentMomentum.P(), trackBeta[iTn]);
                DeltaTOF_expected = GetMassSqr(ParentMomentum.P(), trackBeta[iTn]);
                hXpionXkaon->Fill(NSigmaKaon, NSigmaPion);
                hXpionXee->Fill(NSigmaElectron, NSigmaPion);
                hXpionXpp->Fill(NSigmaProton, NSigmaPion);
                hMassVsPt->Fill(ParentMomentum.Pt(), RhoMass);

                if (ParentMomentum.Pt()<0.150){
                    hMpiMass_0_100MeV->Fill(RhoMass);
                }
                if (ParentMomentum.Pt()>0.150){
                    hKaonMass->Fill(RhoMass);
                }
                if ((DeltaDeltaTOF < 0.750 && DeltaDeltaTOF > 0) && NSigmaPion < 8) {
                    hRhoMass->Fill(RhoMass);
                    if (RhoMass > 0.65 && RhoMass < 0.90 
                        // && fabs(ParentMomentum.Rapidity()) >= 0.6
                        // && fabs(ParentMomentum.Rapidity()) <= 1.0
                        ) {
                        // Defining new reference frame variables.
                        cosphi = (ParentMomentum.Px()*Polarization.Px() + ParentMomentum.Py()*
                               Polarization.Py())/(ParentMomentum.Pt()*Polarization.Pt());
                        sinphi = (ParentMomentum.Px()*Polarization.Py() - ParentMomentum.Py()*
                               Polarization.Px())/(ParentMomentum.Pt()*Polarization.Pt());
                        float phi = TMath::ACos(cosphi);
                        if (ParentMomentum.Pt()*sinphi< 0){
                            phi = -1*phi;
                        }

                        count ++;
                        SurvivingFraction ++;

                        // Histograms
                        // 2D histogram
                        hPxPy->Fill(ParentMomentum.Pt()*cos(phi), ParentMomentum.Pt()*sin(phi));
                        hCos2phivsPT->Fill(ParentMomentum.Pt(), 2*cos(2*phi));

                        // 1D Histograms 
                        hDeltaPhi->Fill(DeltaPhi);
                        hPhipPPhim->Fill(ParentMomentum.Phi());
                        hPhipMPhim->Fill(Polarization.Phi());
                        hcosphi->Fill(cosphi);
                        hPhip->Fill(pTrackMomentum.Phi());
                        hPhim->Fill(nTrackMomentum.Phi());
                        hRhoTP->Fill(ParentMomentum.Pt());
                        if (ParentMomentum.Pt()<0.060) {
                           hphi->Fill(phi);
                           hCos2phivsY->Fill(abs(ParentMomentum.Rapidity()), 2*cos(2*phi));
                           hCos2phi->Fill(2*cos(phi*2));
                        }
                        // Create 2-D histogram Pt (0, 1) cos2phi (0, 2) #bins 200
                        // Tprofile 
                        // Same thing for rapidity instead of 
                    }
                }
            } //end second (negative) track loop


        } // end first (positive) track loop

   } //end event loop

    cout << endl;
    cout << "Event loop done" << endl;
    cout << number_of_zeros << endl;
    cout << "Surviving Fraction of Events After Cuts" << endl;
    cout << SurvivingFraction/EventLimit << endl;
    hRhoMass->Scale(1.0 / hRhoMass->Integral());

    normalizeHistogramProbability(hMpiMass_0_100MeV);
    // CalculateDifferentialDistribution(hMpiMass_0_100MeV, hdNdM);

    normalizeHistogram(hphi);
    hphi->SetMarkerStyle(31);
    // Create the scatter plot
    // hphi->Scale(1.0 / hphi->Integral());

    // Fit function 1-Acos(2phi)
    // Create a TF1 object with your fit function
    TF1* fitFunc = new TF1("fitFunc", myFitFunction, 0, 10, 1); // One parameter (A)

    // Optionally, set some initial parameter values to guide the fit.
    fitFunc->SetParameter(0, 1.0); // Initial guess for A (amplitude)

    // Perform the fit
    hphi->Fit("fitFunc");

    // Create a TProfile to store the 1D projection of the mean value of Y for each bin in X
    TProfile* profileX = hCos2phivsPT->ProfileX();
    profileX->GetXaxis()->SetRangeUser(0, 0.23);

    // TProfile* profileX = hCos2phivsY->ProfileX();
    // profileX->GetXaxis()->SetRangeUser(0, 0.25);
    // Access the fit results
    double fitParameterA = fitFunc->GetParameter(0);
    double fitParameterErrorA = fitFunc->GetParError(0);

    // Print the fit resultss
    cout << "Fit parameter A: " << fitParameterA << " +/- " << fitParameterErrorA << endl;
    fitFunc->Write();
    OutputFile->Write();
    OutputFile->Close();

    // Now let's open the file and draw the histogram with the correct marker style
    TFile* InputFile = new TFile("/Users/daniel/Documents/Codes/C++/STAR_UPC/UpcOutput.root");
    TH1* hphi_drawn = dynamic_cast<TH1*>(InputFile->Get("hphi"));
    TF1* fitFunc_drawn = dynamic_cast<TF1*>(InputFile->Get("fitFunc"));

    if (!hphi_drawn || !fitFunc_drawn) {
        std::cout << "Error: Unable to retrieve histogram or fit function from the file!" << std::endl;
        InputFile->Close();
        return;
    }

    TCanvas* canvas = new TCanvas("canvas", "Histogram Canvas");
    gStyle->SetOptStat(0); // Remove statistics box from the canvas
    hphi_drawn->SetLineStyle(0);
    hphi_drawn->SetMarkerStyle(31); // Set the marker style when drawing
    hphi_drawn->Draw();
    // Optionally, draw the fit function as well
    fitFunc_drawn->SetLineColor(kRed); // Set the line color for the fit function
    fitFunc_drawn->Draw("same");

    canvas->Update();
    canvas->SaveAs("hphi.png");

    InputFile->Close();
}



//------------------------------------------------------
double GetMassSqr(double p, double beta){
  //function to calculate the mass squared from track 3 momentum magnitude and beta from the TOF
  if(beta == -999 or beta==0) return -99.; //"-99" just means no information

  float masssqr = p*p*(1.0/(beta*beta)-1);//mass square formula

  return masssqr ;
}
//--------------------------------------------------------
double GetDeltaTOF(double s, double p, double mass){
  /**
   * Summary: Computes the expected time of a particle moving along the path s in the laboratory
   frame.
   s (double): track lenght.
   p (double): momentum.
   mass (double): Invariant mass of the mass. 
   */
  if(s == -999 or s==0) return -99.; //"-99" just means no information

  float DeltaTOF = s/c*sqrt(1+mass*mass/(p*p));//mass square formula

  return DeltaTOF;
}
//--------------------------------------------------------------
// Function to fit: 1 - A*cos(2*phi)
Double_t myFitFunction(Double_t* x, Double_t* par) {
    return 1.0 - par[0] * cos(2.0 * x[0]);
}
//--------------------------------------------------------------
void normalizeHistogram(TH1* histogram) {
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
    }
}

void normalizeHistogramProbability(TH1F* originalHist) {
    // Calculate the total number of events in the histogram
    float totalEvents = originalHist->Integral();
    
    // Loop over each bin and normalize its content
    int numBins = originalHist->GetNbinsX();
    for (int bin = 1; bin <= numBins; ++bin) {
        float normalizedContent = originalHist->GetBinContent(bin) / totalEvents;
        float normalizedError = originalHist->GetBinError(bin) / totalEvents;

        originalHist->SetBinContent(bin, normalizedContent);
        originalHist->SetBinError(bin, normalizedError);
    }
}
