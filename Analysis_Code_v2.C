#define PionPdgMass 0.139570
#define ProtonPdgMass 0.938272
#define KaonPdgMass 0.493677 // charged kaons, not K0s
#define RhoPdgMass 0.770
#define RhoLifeTime 0.148
#define c 29.9792

#include <TMath.h>

#include <TChain.h>

#include "utils/utils.h"

// TODO: Add Boolean function to check all requirements to be a pion
//

bool check_track(Float_t NSigma_Pion, Float_t NSigma_Kaon);

bool check_pair_of_tracks(float NSigma_Pion_1, float NSigma_Pion_2, float DeltaDeltaTOF, 
                          int Q1, int Q2, float rapidity);

bool check_event(int Total_Tracks);

float compute_phi(TLorentzVector TrackMomentum_1, TLorentzVector TrackMomentum_2, TLorentzVector TotalTrackMomentum);

void Analysis_Code_v2(const Char_t *outputFile = "/Users/daniel/Documents/Research/BNL/STAR_UPC/output_files/UPC_output.root") {

    TFile *myFile = TFile::Open("./Data/FemtoDst_Run10AuAu_wZDC.root");
    TTreeReader myReader("FemtoDst", myFile);

    // std::ofstream outFile("output.txt");
    // std::streambuf *coutbuf = std::cout.rdbuf(); // Save old buf
    // std::cout.rdbuf(outFile.rdbuf()); // Redirect std::cout to outFile

    /*** Track Histograms ***/ 
    
    // Cos(2phi) vs Transverse momentum vs Invariant Mass
    //
    TH3F* H3D2Cos2phivsPtvsMass = new TH3F("H3DCos2phivsPtvsMass", "2cos(2#phi) vs P_{T} vs M_{#pi#pi}",
                                            200, -2.0, 2.0,
                                            100, 0.0, 1.0,
                                            400, 0.4, 1.2); 

    TH3F * H3D_PhiVsPtVsMass = new TH3F("H3D_PhiVsPtVsMass", "#phi vs P_{T} vs M_{#pi#pi}",
                                        100, -TMath::Pi(), TMath::Pi(),
                                        500, 0., 0.5,                   
                                        static_cast<int>((1.3 - 0.4) / 0.0025), 0.4, 1.3);               

    // Mass vs Transverse Momentum signals
    //
    // Todo: 
    // Change the order of mass and Pt so that it matches the name of the histogram.
    TH2F *hMassVsPt = new TH2F("hMassVsPt", "#pi^{+}#pi^{-} invariant Mass vs Transverse momentum",
                                1000, 0, 1.0,   // Pt range
                                static_cast<int>((1.3 - 0.4) / 0.0025), 0.4, 1.3); // Mass range

    // Track Energy loss vs q momentum signals
    //
    TH2F * hdEdxVsP  = new TH2F("hdEdxVsP", "Track energy loss (dE/dx) vs q*momentum (GeV/c)", 
                                 1000, -2., 2.,  // Momentum Range
                                 1000, 2.0, 4.); // dEdx range 


    // Track nSigma Pion vs q momentum signals
    //
    TH2F * hNSigmaPionVsP  = new TH2F("hNSigmaPionVsP", "n#sigma_{#pi} vs q*p (GeV/c)", 
                                 1000, -3., 3.,  // Momentum range
                                 1000, -4., 5.); // NSigmaPion range
        
    // Charge signal
    //
    TH1F *hCharge = new TH1F("hCharge", "Charge Distribution", 11, -5.5, 5.5);

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

    unsigned int TotTracks;  // # of tracks in given event

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

    // TLorentzVector is a 4-vector class (https://root.cern.ch/doc/master/classTLorentzVector.html)
    //
    TLorentzVector TrackMomentum_1;
    TLorentzVector TrackMomentum_2;
    TLorentzVector TotalTrackMomentum;
    TLorentzVector Polarization;


    // Loop over all entries in the TChain
    //
    float SurvivingFraction = 0;
    unsigned int iEvent = 0; // event loop iterator index
    cout << EventsInFile << endl;
    while (myReader.Next() && iEvent < EventLimit) {
        iEvent ++;
        TotTracks = trackPt.GetSize();

        // Check event conditions to perform the analysis.
        //
        if (TotTracks != 2) continue;
        // if (check_event(TotTracks)) {
        //     continue;
        // }
        // cout << "Event # :" << iEvent << endl;
        // Little block of code to keep track of % done running code. Can be deleted.
        //
        fraction = iEvent * 100.0 / totalevents;
        if (fraction % 1 == 0 && fraction != prev_frac)
        {
            printf("\b\b\b%2d%%", fraction);
            std::cout << std::flush;
            prev_frac = fraction;
        }
        // Iterates over all tracks except the last one.
        //
        // for (size_t iTrack_1 = 0; iTrack_1 < (TotTracks - 1); iTrack_1++) {
        for (size_t iTrack_1 = 0; iTrack_1 < TotTracks; iTrack_1++) {
            if (trackQ[iTrack_1] < 0) continue;
            TrackMomentum_1.SetPtEtaPhiM(trackPt[iTrack_1], trackEta[iTrack_1], trackPhi[iTrack_1], PionPdgMass);
            hdEdxVsP->Fill(trackQ[iTrack_1]*TrackMomentum_1.P(), trackdEdx[iTrack_1]);
            hCharge->Fill(trackQ[iTrack_1]);
            hNSigmaPionVsP->Fill(trackQ[iTrack_1]*TrackMomentum_1.P(), trackNSigmaPion[iTrack_1]);
            // if (check_track(trackNSigmaPion[iTrack_1], trackNSigmaKaon[iTrack_1])) continue;
            float t1 = GetDeltaTOF(trackLength[iTrack_1], TrackMomentum_1.P(), PionPdgMass);
            // Iterates over all tracks starting from track_1 to avoid repeating two or more times the same pair of tracks.
            //
            // for (size_t iTrack_2 = iTrack_1 + 1; iTrack_2 < TotTracks; iTrack_2++) {
            for (size_t iTrack_2 = 0; iTrack_2 < TotTracks; iTrack_2++) {
                if (trackQ[iTrack_2] > 0) continue;
                int charge1 = static_cast<int>(trackQ[iTrack_1]);
                int charge2 = static_cast<int>(trackQ[iTrack_2]);
                TrackMomentum_2.SetPtEtaPhiM(trackPt[iTrack_2], trackEta[iTrack_2], trackPhi[iTrack_2], PionPdgMass);
                TotalTrackMomentum = TrackMomentum_1 + TrackMomentum_2;

                float t2 = GetDeltaTOF(trackLength[iTrack_2], TrackMomentum_2.P(), PionPdgMass);
                float DeltaTOF = abs(trackTOF[iTrack_2] - trackTOF[iTrack_1]);
                float DeltaTOF_expected = abs(t2-t1);
                float DeltaDeltaTOF = abs(DeltaTOF - DeltaTOF_expected);

                float rapidity = TotalTrackMomentum.Rapidity();

                if (check_pair_of_tracks(trackNSigmaPion[iTrack_1], trackNSigmaPion[iTrack_2], DeltaDeltaTOF,
                      trackQ[iTrack_1], trackQ[iTrack_2], rapidity)) continue;
                float phi = compute_phi(TrackMomentum_1, TrackMomentum_2, TotalTrackMomentum);
                SurvivingFraction ++;
                hMassVsPt->Fill(TotalTrackMomentum.Pt(), TotalTrackMomentum.M());
                H3D_PhiVsPtVsMass->Fill(phi, TotalTrackMomentum.Pt(), TotalTrackMomentum.M());
                H3D2Cos2phivsPtvsMass->Fill(2*cos(2*phi), TotalTrackMomentum.Pt(), TotalTrackMomentum.M());
            }
        } // end of track loop
    } // end event loop
    // cout << endl;
    // cout << "Event loop done" << endl;
    // cout << hMassVsPt->GetEntries() << endl;
    cout << "Surviving Fraction of Events After Cuts" << endl;
    cout << SurvivingFraction / EventLimit << endl;

    // std::cout.rdbuf(coutbuf); // Reset to standard output again
    // outFile.close();

    TFile * OutputFile = new TFile( outputFile, "RECREATE" );

    // Save Histograms
    //
    OutputFile->cd();

    hCharge       ->Write();
    hNSigmaPionVsP->Write();
    hdEdxVsP      ->Write();

    H3D2Cos2phivsPtvsMass->Write();

    hMassVsPt->Write();

    H3D_PhiVsPtVsMass->Write();

    OutputFile->Close();
}

/**
 * @brief Check if a single track has the requirements to go through the analysis.
 *
 * @param NSigma_Pion number of standard deviations from being a pion 
 * @param NSigma_Kaon number of standard deviations from being a kaon
 * @return true:  if the track does not satisfy the criteria to be part of the analysis
 * @return false: if the track satisfy the criteria to be part of the analysis
 */
bool check_track(Float_t NSigma_Pion, Float_t NSigma_Kaon) {
    // Check single track NSigma pion and NSigma Kaon
    //
    if (abs(NSigma_Pion) >= 3.0) return true;

    else return false;
}

/**
 * @brief Check if a pair of particles satisfy the criteria for the analysis.
 * 
 * @param NSigma_Pion_1 number of standard deviations from being a pion of particle 1
 * @param NSigma_Pion_2 number of standard deviations from being a pion of particle 2
 * @return true: if the pair does not satisfy the criteria to be part of the analysis.
 * @return false: if the pair satisfy the criteria to be part of the analysis.
 */
bool check_pair_of_tracks(float NSigma_Pion_1, float NSigma_Pion_2, float DeltaDeltaTOF, 
                          int Q1, int Q2, float rapidity) {
    float NSigmaPion = NSigma_Pion_1*NSigma_Pion_1 + NSigma_Pion_2*NSigma_Pion_2;
    if (NSigmaPion >= 8) return true;
    else if (abs(NSigma_Pion_1) >=3 || abs(NSigma_Pion_2) >=3) return true;
    else if (DeltaDeltaTOF >= 0.750 ||  DeltaDeltaTOF <= 0.) return true;
    else if (Q1*Q2 == 1) return true;
    else if (rapidity >= 1) return true;
    else return false;
}

/**
 * @brief Check if an event satisfy the criteria to go into the analysis.
 * 
 * @param Total_Tracks Total number of tracks of the event.
 * @return true: If the pair does not satisfy the criteria to be part of the analysis.
 * @return false: If the pair satisfy the criteria to be part of the analysis.
 */
bool check_event(int Total_Tracks) {
    // If there are no tracks just skip
    //
    if (Total_Tracks != 2) return true;

    else return false;
}

/**
 * @brief Compute the angle between the sum of the momentum of two tracks and the difference of the momentum of
 * both tracks.
 * 
 * @param TrackMomentum_1 Momentum track 1
 * @param TrackMomentum_2 Momentum track 2
 * @param TotalTrackMomentum TrackMomentum_1 + TrackMomentum_2
 * @return float phi
 */
float compute_phi(TLorentzVector TrackMomentum_1, TLorentzVector TrackMomentum_2, TLorentzVector TotalTrackMomentum) {
    TLorentzVector Polarization = TrackMomentum_1 - TrackMomentum_2;
    float cosphi = (TotalTrackMomentum.Px()*Polarization.Px() + TotalTrackMomentum.Py()*
                    Polarization.Py())/(TotalTrackMomentum.Pt()*Polarization.Pt());
    float sinphi = (TotalTrackMomentum.Px()*Polarization.Py() - TotalTrackMomentum.Py()*
                    Polarization.Px())/(TotalTrackMomentum.Pt()*Polarization.Pt());

    float phi = TMath::ATan2(sinphi, cosphi);
    return phi;
}