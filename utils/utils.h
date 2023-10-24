// Useful functions
#define c   29.9792

double GetMassSqr(double p, double beta){
    /**
    * @brief Function to calculate the mass squared from track 3 momentum magnitude and beta from
    * TOF.
    * 
    * @param p: 3 momentum magnitude.
    * @param beta: Beta value from TOF 
    * @return double: Mass squared value. 
    */

    if(beta == -999 or beta==0) return -99.; //"-99" just means no information

    float masssqr = p*p*(1.0/(beta*beta)-1);//mass square formula

    return masssqr;
}


double GetDeltaTOF(double s, double p, double mass){
    /**
    * @brief Computes the expected time of a particle moving along the path s in the laboratory
    * frame.
    * 
    * @param s: track lenght.
    * @param p: 3 momentum
    * @param mass 
    * @return double: Invariant mass of the mass. 
    */

    if(s == -999 or s==0) return -99.; //"-99" just means no information

    float DeltaTOF = s/c*sqrt(1+mass*mass/(p*p));//mass square formula

    return DeltaTOF;
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

void normalizeHistogramProbability(TH1F* histogram) {
    /**
     * @brief Normalize histogram so that the integral under the curve is equal to 1.
     * 
     * @param histogram: Histogram to normalize.
     */

    // Calculate the total number of events in the histogram
    float totalEvents = histogram->Integral();
    
    // Loop over each bin and normalize its content
    int numBins = histogram->GetNbinsX();
    for (int bin = 1; bin <= numBins; ++bin) {
        float normalizedContent = histogram->GetBinContent(bin) / totalEvents;
        float normalizedError = histogram->GetBinError(bin) / totalEvents;

        histogram->SetBinContent(bin, normalizedContent);
        histogram->SetBinError(bin, normalizedError);
    }
}


void AddEmptyLegendEntry(TLegend *legend, const char *text = "", int color = 0) {
    /**
    * @brief Funtion to add a legend with no indicator.
    * 
    * @param legend: TLegend object.
    * @param text: Legend text.
    * @param color: Color of text.
    */

    // Create a "dummy" histogram with no contents and assign it a color
    TH1F *dummy = new TH1F(Form("dummy%d", color), "", 0, 0, 0);
    dummy->SetLineColor(color);

    // Add the dummy entry to the legend with the specified text
    legend->AddEntry(dummy, text, "l");

    // Make sure the dummy histogram doesn't appear in the plot
    dummy->SetFillColor(0);
    dummy->SetLineColor(0);
    dummy->SetLineWidth(0);
    dummy->SetMarkerSize(0);
}

void SetCanvasStyle(TCanvas* canvas) {
    // Set ROOT style options
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLineWidth(2);

    // Set canvas margins and frame properties
    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.11);
    canvas->SetLeftMargin(0.11);
    canvas->SetRightMargin(0.05);
    canvas->SetFrameLineWidth(2);

    // Create and configure TLatex object
    TLatex tl;
    tl.SetTextSize(0.06);
    tl.SetNDC();
}