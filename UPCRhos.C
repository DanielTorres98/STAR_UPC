#define PionPdgMass        0.139570
#define ProtonPdgMass      0.938272
#define KaonPdgMass        0.493677  //charged kaons, not K0s
#define RhoPdgMass         0.770

//I want to define a function. I have to define it before the main function "FindRhos". I'll write
// what this function does *after* FindRhos. I would encourage you to define functions for any
// repeated code, e.g. track cuts
double GetMassSqr(double p, double beta);



void FindRhos() {
    // TChain tc = new TChain( "FemtoDst" );
    // tc->Add( "/Users/jdb/Dropbox/bnl/data/FemtoDst_Run10AuAu_tofMult_allTracks.root" );

    // TFile *myFile = TFile::Open("/Users/jdb/bnl/work/upc/data/FemtoDst_Run12CuAuPART1.root");
    // TFile *myFile = TFile::Open("./FemtoDst_Run10AuAu_tofMult_allTracks.root");
    TFile *myFile = TFile::Open("./Data/FemtoDst_Run10AuAu_wZDC.root");

    TTreeReader myReader("FemtoDst", myFile);

    // Output file where histograms go. Histograms defined below a TFile are automatically
    // associated with that TFile.
    TFile * OutputFile = new TFile( "UpcOutput.root", "RECREATE" );
    // 2D histogram. Defined by (name, title, #binsx, lower range x, upper range x, binsy, low y,
    // upper y)
    //
    TH2F * hdEdxVsP = new TH2F("hdEdxVsP", "Track energy loss (dE/dx) vs q*momentum (GeV/c)", 
                               1000, -3., 3., 1000, 0., 10.);

    TH2F * hTofVsP = new TH2F("hTofVsP", "TOF m^2 (GeV/c) vs momentum (GeV/c)", 1000, 0, 3.5,
                              1000, 0.0, 1.2);

    //1d histogram. Defined by (name, title, #binsx, lower range x, upper range x)
    //
    TH1F * hRhoMass = new TH1F("hRhoMass", "#rho mass distribution", 100, 0.6, 0.9);
    hRhoMass->GetXaxis()->SetTitle("m_{#rho}(GeV)");
    hRhoMass->GetYaxis()->SetTitle("Counts");

    // Transverse Momentum histogram
    //
    TH1F * hRhoTP = new TH1F("hRhoTP", "#rho Transverse Momentum", 10000, 0., 1.);

    // Deltaphi histogram: Δϕ = ϕ(π⁺ + π⁻) − ϕ(π⁺ − π⁻)
    //
    TH1F * hDeltaPhi = new TH1F("hDeltaPhi", "#Delta#phi", 80, -3.14, 3.14);
    hDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    hDeltaPhi->GetYaxis()->SetTitle("Counts");
    // PhipPPhim histogram ϕ(π⁺ + π⁻)
    TH1F * hPhipPPhim = new TH1F("hPhipPPhim", "#phi(#pi^{+} + #pi^{-})", 80, -3.14, 3.14);
    hPhipPPhim->GetXaxis()->SetTitle("#phi(#pi^{+} + #pi^{-})");
    hPhipPPhim->GetYaxis()->SetTitle("Counts");
    // PhipMPhim histogram ϕ(π⁺ - π⁻)
    TH1F * hPhipMPhim = new TH1F("hPhipMPhim", "#phi(#pi^{+} - #pi^{-})", 80, -3.14, 3.14);
    hPhipMPhim->GetXaxis()->SetTitle("#phi(#pi^{+} - #pi^{-})");
    hPhipMPhim->GetYaxis()->SetTitle("Counts");
    // Phip histogram
    TH1F * hPhip = new TH1F("hPhip", "#pi^{+}", 80, -3.14, 3.14);
    hPhip->GetXaxis()->SetTitle("#phi");
    hPhip->GetYaxis()->SetTitle("Counts");
    // Phim histogram
    TH1F * hPhim = new TH1F("hPhim", "#pi^{-}", 80, -3.14, 3.14);
    hPhim->GetXaxis()->SetTitle("#phi");
    hPhim->GetYaxis()->SetTitle("Counts");

    TTreeReaderArray<Float_t> trackPt(myReader, "Tracks.mPt"); //transverse momentum in GeV
    TTreeReaderArray<Float_t> trackEta(myReader, "Tracks.mEta"); //pseudorapidity
    TTreeReaderArray<Float_t> trackPhi(myReader, "Tracks.mPhi"); //azimuthal angle in STAR coordinate system
    TTreeReaderArray<Char_t> trackQ(myReader, "Tracks.mQ"); //track charge +/- 1 // This was written as a char, I guess > 0 means + and < 0 means -

    TTreeReaderArray<Float_t> trackDCA(myReader, "Tracks.mDCA"); //distance of closest approach to primary vertex
    TTreeReaderArray<Float_t> trackPi(myReader, "Tracks.mNSigmaPion"); //number of standard deviations from dE/dx fit for pions. A standard cut would be something like |nSigmaPion| < 2. The number can be negative for above/below the curve, but it's rare to care.
    TTreeReaderArray<Float_t> trackdEdx(myReader, "Tracks.mDedx"); //energy lost in TPC

    TTreeReaderArray<Float_t> trackBeta(myReader, "BTofPidTraits.mBTofBeta"); //velocity fraction of track as measured by TOF
    TTreeReaderArray<Float_t> trackTof(myReader, "BTofPidTraits.mBTof"); //TOF value measured

    vector<TLorentzVector> uusiglv; // uu pairs with daughters from k0s mass range
    vector<TLorentzVector> uunotsiglv; // uu pairs with daughters not from k0s mass range



    unsigned int nPos, nNeg;
    unsigned int iEvent = 0; // event loop iterator index
    unsigned int TotTracks; //# of tracks in given event

    int prev_frac = -1;
    int fraction;
    int EventsInFile = myReader.GetEntries();
    cout << "There are " << EventsInFile << " events in this file" << endl;

    //TLorentzVector is a 4-vector class (https://root.cern.ch/doc/master/classTLorentzVector.html)
    //I'll declare the + and - track vectors here and then reset the values for every track.
    //
    TLorentzVector pTrackMomentum;
    TLorentzVector nTrackMomentum;
    TLorentzVector TrackMomentum;
    TLorentzVector RhoMomentum;
    TLorentzVector Polarization;

    float pTrackMSqr, nTrackMSqr, TrackMSqr, RhoMsqr;

    // Set smaller if you want a small fraction of the data. If you want all just make this a big
    // number
    //
    int EventLimit = EventsInFile;
    cout << "The event limit is " << EventLimit << " events" << endl; 
    int totalevents = std::min(EventsInFile, EventLimit);

    //"while" loop to go over all of the events. This could also be written as a "for" loop
    //
    int number_of_zeros = 0;
    
    while (myReader.Next() && iEvent < EventLimit) {
        iEvent ++;
        TotTracks = trackPt.GetSize();
        std::cout << "The size of trackPt is: " << trackPt.GetSize() << std::endl;

        //get out of loop if the event has anything other than 2 tracks. You will probably want to
        //uncomment this cut later, but if you want to look at more kinds of tracks it might be
        //good to comment out now. 
        //if(TotTracks != 2) continue; 


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

        for (unsigned int iTp = 0; iTp < TotTracks; iTp++ ){
            //4 momentum is uniquely defined by pT, eta, phi, and mass. We assume the tracks are
            //pions and later evaluate the truth of this idea and correct accordingly.
            TrackMomentum.SetPtEtaPhiM(trackPt[iTp], trackEta[iTp], trackPhi[iTp], PionPdgMass);
            TrackMSqr = GetMassSqr(TrackMomentum.P(), trackBeta[iTp]);
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

                nTrackMomentum.SetPtEtaPhiM( trackPt[iTn], trackEta[iTn], trackPhi[iTn], PionPdgMass);

                nTrackMSqr = GetMassSqr(nTrackMomentum.P(), trackBeta[iTn]);
                // hTofVsP->Fill(nTrackMomentum.P(), nTrackMSqr);
                // add cut on nSigma and (later) mass squared here.

                RhoMomentum = nTrackMomentum + pTrackMomentum;
                Polarization = nTrackMomentum - pTrackMomentum;
                float DeltaPhi = RhoMomentum.Phi() - Polarization.Phi();
                float RhoMass = RhoMomentum.M();
                RhoMsqr = GetMassSqr(RhoMomentum.P(), trackBeta[iTn]);
                // cout << "GetMassSqr->" << RhoMsqr << " m^2->" << RhoMomentum.M2() <<endl;
                // cout << "P->" << RhoMomentum.P()<< " beta->" << trackBeta[iTn] <<endl;
                // hTofVsP->Fill(RhoMomentum.P(), RhoMomentum.M());
                if (abs(trackPi[iTp])>0.66) continue;
                // cout << trackPi[iTp] << endl;
                if (abs(RhoMass - 0.770)<0.1) {
                    hDeltaPhi->Fill(DeltaPhi);
                    hPhipPPhim->Fill(RhoMomentum.Phi());
                    hPhipMPhim->Fill(Polarization.Phi());
                    hPhip->Fill(pTrackMomentum.Phi());
                    hPhim->Fill(nTrackMomentum.Phi());
                    hRhoMass->Fill(RhoMass);
                    hRhoTP->Fill(RhoMomentum.Pt());
                }
                // here you can plot rho mass and try to see a peak

                // here you can add a cut around the rho mass that you want (it should be something like 0.6-0.95 or whatever)

                // here you would start doing other things with the rho, maybe plot rho pT

            } //end second (negative) track loop


        } // end first (positive) track loop

   } //end event loop
   cout << endl;
   cout << "Event loop done" << endl;
   cout << number_of_zeros << endl;
   OutputFile->Write();
   OutputFile->Close();
}



//------------------------------------------------------
double GetMassSqr(double p, double beta){
  //function to calculate the mass squared from track 3 momentum magnitude and beta from the TOF
  if(beta == -999 or beta==0) return -99.; //"-99" just means no information

  float masssqr = p*p*(1.0/(beta*beta)-1);//mass square formula

  return masssqr ;
}
//--------------------------------------------------------------