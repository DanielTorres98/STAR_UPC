#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770
#define RhoLifeTime        0.148
#define c                  29.9792

#include <complex>
#include <cmath>
#include <TMath.h>

#include "utils/utils.h"

void Analysis_Code() {
    TFile *myFile = TFile::Open("./Data/FemtoDst_Run10AuAu_wZDC.root");
    TTreeReader myReader("FemtoDst", myFile);

    // Output file where histograms go. Histograms defined below a TFile are automatically
    // associated with that TFile.
    TFile * OutputFile = new TFile( "UpcOutput_Rho.root", "RECREATE" );

    // 2D histograms. Defined by (name, title, #binsx, lower range x, upper range x, binsy, low y,
    // upper y)
    //
    TH2F * hdEdxVsP = new TH2F("hdEdxVsP", "Track energy loss (dE/dx) vs q*momentum (GeV/c)", 
                               1000, -3., 3., 1000, 0., 10.);
    TH2F * hMassVsPt = new TH2F("hMassVsPt", "Invariant Mass vs Transverse momentum", 
                               1000, 0, 0.1, 1000, 0.3, 1);
    hMassVsPt->GetXaxis()->SetTitle("P_{T}");
    hMassVsPt->GetYaxis()->SetTitle("M_{#rho}");
    TLine *line = new TLine(hMassVsPt->GetXaxis()->GetXmin(), 0.770, hMassVsPt->GetXaxis()->GetXmax(), 0.770);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
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

    TH1F * hRhoMass = new TH1F("hRhoMass", "#rho mass distribution", 0.8 / 0.0025, 0.5, 1.3);
    
    // Mass plots
    //
    TH1F * hRhoMass_0_20 = new TH1F("hRhoMass_0_20", "#pi^{+}#pi^{-} mass distribution", 100, 0.6, 0.9);
    TH1F * hRhoMass_20_40 = new TH1F("hRhoMass_20_40", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_40_60 = new TH1F("hRhoMass_40_60", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_60_80 = new TH1F("hRhoMass_60_80", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_80_100 = new TH1F("hRhoMass_80_100", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_100_150 = new TH1F("hRhoMass_100_150", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_150_200 = new TH1F("hRhoMass_150_200", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hRhoMass_200_1000 = new TH1F("hRhoMass_200_1000", "#pi^{+} #pi^{-} pairs", 100, 0.6, 0.9);
    TH1F * hMpiMass_0_100MeV = new TH1F("hMpiMass_0_100MeV", "#pi^{+}#pi^{-} mass distribution p_{T} < 100 MeV",
                                    static_cast<int>((1.3 - 0.5) / 0.0025), 0.5, 1.3);

    TH1F * hMpiMass_0_100MeV_v2 = new TH1F("hMpiMass_0_100MeV_v2", "#pi^{+}#pi^{-} mass distribution P_{T} < 100 MeV",
                                        100, 0.44, 1.1);
    TH1F * hMpiMass_0_100MeV_v3 = new TH1F("hMpiMass_0_100MeV_v3", "#pi^{+}#pi^{-} mass distribution P_{T} < 100 MeV",
                                        100, 0.44, 1.1);
    TH1F * hKaonMass = new TH1F("hKaonMass", "#pi^{+} #pi^{-} pairs", 100, 0.44, 0.56);
    TH1F * hRhoMass_NoCut = new TH1F("hRhoMass_NoCut", "#rho mass distribution", 100, 0.6, 0.9);
    TH1F * hRhoTP = new TH1F("hRhoTP", "#rho Transverse Momentum", 80, 0., 0.35);
    TH1F * hDDTOF = new TH1F("hDDTOF", "#Delta#Delta TOF", 100, 0, 100);
    TH1F * hDeltaPhi = new TH1F("hDeltaPhi", "#Delta#phi", 80, -3.14, 3.14);
    TH1F * hcosphi = new TH1F("hcosphi", "cos#phi", 80, -1, 1);
    TH1F * hPhipPPhim = new TH1F("hPhipPPhim", "#phi(#pi^{+} + #pi^{-})", 80, -3.14, 3.14);
    TH1F * hPhipMPhim = new TH1F("hPhipMPhim", "#phi(#pi^{+} - #pi^{-})", 80, -3.14, 3.14);
    TH1F * hPhip = new TH1F("hPhip", "#pi^{+}", 80, -3.14, 3.14);
    TH1F * hPhim = new TH1F("hPhim", "#pi^{-}", 80, -3.14, 3.14);
    TH1F * hphi = new TH1F("hphi", "#phi", 100, -TMath::Pi(), TMath::Pi());
    TH1F * hCos2phi = new TH1F("hCos2phi", "cos2#phi", 80, -1, 1);

    // Variables from Analysis Tree.
    //
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

    // ZDC
    TTreeReaderArray<Float_t> EventZDCEast(myReader, "mZDCEast");
    TTreeReaderArray<Float_t> EventZDCWest(myReader, "mZDCWest");

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
    TLorentzVector RhoMomentum;
    TLorentzVector Polarization;

    //"while" loop to go over all of the events. This could also be written as a "for" loop
    //
    int number_of_zeros = 0;
    float SurvivingFraction = 0.0;
    
    while (myReader.Next() && iEvent < EventLimit) {
        iEvent ++;
        TotTracks = trackPt.GetSize();

        //get out of loop if the event has anything other than 2 tracks. You will probably want to
        //uncomment this cut later, but if you want to look at more kinds of tracks it might be
        //good to comment out now. 
        // if (EventZDCEast[iEvent]<2 & EventZDCEast[iEvent]<2) continue;
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
        // I'll make nested loops here. So a loop of negative tracks inside of a loop of positive
        // tracks. In principle you only care about events with two tracks (I made that cut
        // earlier), but I want to keep the code general.
        // start first track loop, over + tracks.
        //
        bool found_zero = false;
        int count;
        count = 1;
        for (unsigned int iTp = 0; iTp < TotTracks; iTp++ ){
            //4 momentum is uniquely defined by pT, eta, phi, and mass. We assume the tracks are
            //pions and later evaluate the truth of this idea and correct accordingly.
            //
            TrackMomentum.SetPtEtaPhiM(trackPt[iTp], trackEta[iTp], trackPhi[iTp], PionPdgMass);
            float TrackMSqr = GetMassSqr(TrackMomentum.P(), trackBeta[iTp]);
            float t1 = GetDeltaTOF(trackLength[iTp], nTrackMomentum.P(), PionPdgMass);
            float NSigmaPion1 = trackNSigmaPion[iTp];
            float NSigmaKaon1 = trackNSigmaKaon[iTp];
            float NSigmaElectron1 = trackNSigmaElectron[iTp];
            float NSigmaProton1 = trackNSigmaProton[iTp];
            // cout << "GetMassSqr->" << TrackMSqr << " m^2->" << TrackMomentum.M2() <<endl;
            // cout << "P->" << TrackMomentum.P()<< " beta->" << trackBeta[iTp] <<endl;
            hTofVsP->Fill(TrackMomentum.P(), trackBeta[iTp]);

            if (!found_zero && trackBeta[iTp] == 0) {
                found_zero = true;
                number_of_zeros += 1;
            }

            // This loop is for positive pions, but before I make cuts you can fill general track 
            // information here fill 2d histogram with x value, yvalue
            // it's convention to separate negative and positive particles. There is no deep meaning.
            //
            pTrackMomentum.SetPtEtaPhiM(trackPt[iTp], trackEta[iTp], trackPhi[iTp], PionPdgMass);
            if(trackQ[iTp] > 0) hdEdxVsP->Fill(pTrackMomentum.P(), trackdEdx[iTp]);

            else hdEdxVsP->Fill(-1.*pTrackMomentum.P(), trackdEdx[iTp]);
            //get out of this loop if the charge is < 0, so beyond
            //this all tracks indexed by iTp should be positive
            if(trackQ[iTp] < 0) continue; 

            //call function defined above and detailed below
            float pTrackMSqr = GetMassSqr(pTrackMomentum.P(), trackBeta[iTp]); 
            // add cut on nSigma and (later) mass squared here.
            //negative track loop, nested
            for (unsigned int iTn = 0; iTn < TotTracks; iTn++ ){

                if(trackQ[iTn] > 0) continue; //negative track cut

                float NSigmaPion2 = trackNSigmaPion[iTn];
                float NSigmaKaon2 = trackNSigmaKaon[iTn];
                float NSigmaElectron2 = trackNSigmaElectron[iTn];
                float NSigmaProton2 = trackNSigmaProton[iTn];

                float NSigmaPion = NSigmaPion1*NSigmaPion1 + NSigmaPion2*NSigmaPion2;
                float NSigmaKaon = NSigmaKaon1*NSigmaKaon1 + NSigmaKaon2*NSigmaKaon2;
                float NSigmaElectron = NSigmaElectron1*NSigmaElectron1 +
                                       NSigmaElectron2*NSigmaElectron2;
                float NSigmaProton = NSigmaProton1*NSigmaProton1 + NSigmaProton2*NSigmaProton2;

                nTrackMomentum.SetPtEtaPhiM( trackPt[iTn], trackEta[iTn], trackPhi[iTn], PionPdgMass);

                float nTrackMSqr = GetMassSqr(nTrackMomentum.P(), trackBeta[iTn]);
                float t2 = GetDeltaTOF(trackLength[iTn], nTrackMomentum.P(), PionPdgMass);

                float DeltaTOF = abs(trackTOF[iTn] - trackTOF[iTp]);
                float DeltaTOF_expected = abs(t2-t1);

                float DeltaDeltaTOF = abs(DeltaTOF - DeltaTOF_expected);

                RhoMomentum = nTrackMomentum + pTrackMomentum;
                Polarization = nTrackMomentum - pTrackMomentum;
                float DeltaPhi = RhoMomentum.Phi() - Polarization.Phi();
                float RhoMass = RhoMomentum.M();
                float RhoMsqr = GetMassSqr(RhoMomentum.P(), trackBeta[iTn]);

                hXpionXkaon->Fill(NSigmaKaon, NSigmaPion);
                hXpionXee->Fill(NSigmaElectron, NSigmaPion);
                hXpionXpp->Fill(NSigmaProton, NSigmaPion);
                hRhoMass_NoCut->Fill(RhoMass);
                hMassVsPt->Fill(RhoMomentum.Pt(), RhoMass);
                if (RhoMomentum.Pt()<0.02){
                    hRhoMass_0_20->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.02 && RhoMomentum.Pt()<0.04){
                    hRhoMass_20_40->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.04 && RhoMomentum.Pt()<0.06){
                    hRhoMass_40_60->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.06 && RhoMomentum.Pt()<0.08){
                    hRhoMass_60_80->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.08 && RhoMomentum.Pt()<0.1){
                    hRhoMass_80_100->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.1 && RhoMomentum.Pt()<0.15){
                    hRhoMass_100_150->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.15 && RhoMomentum.Pt()<0.2){
                    hRhoMass_150_200->Fill(RhoMass);
                }
                if (RhoMomentum.Pt()>0.2 && RhoMomentum.Pt()<1){
                    hRhoMass_200_1000->Fill(RhoMass);
                }
                if (RhoMomentum.Pt() < 0.100 && std::abs(RhoMomentum.Rapidity())<1) {  
                    hRhoMass->Fill(RhoMass);
                }
                if ((DeltaDeltaTOF < 0.750 && DeltaDeltaTOF > 0) && (NSigmaPion < 8)) {
                    if (RhoMomentum.Pt()<0.100 && std::abs(RhoMomentum.Rapidity())<1){
                        count ++;
                        SurvivingFraction ++;

                        hMpiMass_0_100MeV->Fill(RhoMass);
                        hMpiMass_0_100MeV_v2->Fill(RhoMass);
                        hMpiMass_0_100MeV_v3->Fill(RhoMass);
                        hKaonMass->Fill(RhoMass);
                    }
                    if (RhoMass > 0.65 && RhoMass < 0.90) {
                        // Defining new reference frame variables.
                        //
                        float cosphi = (RhoMomentum.Px()*Polarization.Px() + RhoMomentum.Py()*
                               Polarization.Py())/(RhoMomentum.Pt()*Polarization.Pt());
                        float sinphi = (RhoMomentum.Px()*Polarization.Py() - RhoMomentum.Py()*
                               Polarization.Px())/(RhoMomentum.Pt()*Polarization.Pt());
                        float phi = TMath::ACos(cosphi);
                        if (RhoMomentum.Pt()*sinphi< 0){
                            phi = -1*phi;
                        }

                        // 2D Histogram
                        //
                        hPxPy->Fill(RhoMomentum.Pt()*cos(phi), RhoMomentum.Pt()*sin(phi));
                        hCos2phivsPT->Fill(RhoMomentum.Pt(), 2*cos(2*phi));

                        // 1D Histograms
                        //
                        hDeltaPhi->Fill(DeltaPhi);
                        hPhipPPhim->Fill(RhoMomentum.Phi());
                        hPhipMPhim->Fill(Polarization.Phi());
                        hcosphi->Fill(cosphi);
                        hPhip->Fill(pTrackMomentum.Phi());
                        hPhim->Fill(nTrackMomentum.Phi());
                        hRhoTP->Fill(RhoMomentum.Pt());
                        if (RhoMomentum.Pt()<0.060) {
                           hphi->Fill(phi);
                           hCos2phivsY->Fill(abs(RhoMomentum.Rapidity()), 2*cos(2*phi));
                           hCos2phi->Fill(2*cos(2*phi));
                        }
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
    hRhoMass_NoCut->Scale(1.0 / hRhoMass_NoCut->Integral());

    double nEntries = hRhoMass->GetEntries();
    hRhoMass->Scale(1/(1074600.0*(nEntries/EventLimit)));
    hRhoMass_0_20->Scale(1.0 / hRhoMass_0_20->Integral());
    hRhoMass_20_40->Scale(1.0 / hRhoMass_20_40->Integral());
    hRhoMass_40_60->Scale(1.0 / hRhoMass_40_60->Integral());
    hRhoMass_60_80->Scale(1.0 / hRhoMass_60_80->Integral());
    hRhoMass_80_100->Scale(1.0 / hRhoMass_80_100->Integral());
    hRhoMass_100_150->Scale(1.0 / hRhoMass_100_150->Integral());
    hRhoMass_150_200->Scale(1.0 / hRhoMass_150_200->Integral());
    hRhoMass_200_1000->Scale(1.0 / hRhoMass_200_1000->Integral());

    double nEntries_2 = hMpiMass_0_100MeV->GetEntries();
    hMpiMass_0_100MeV->Scale(1/(1074600.0*(nEntries_2/EventLimit)));

    // Perform the second fit and set its line and marker color.

    normalizeHistogram(hphi);

    OutputFile->Write();
    OutputFile->Close();
}
