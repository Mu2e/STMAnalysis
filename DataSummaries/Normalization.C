// Generates the table of muon captures per incident signal photon and the normalization plot
// Usage example - $ root -l -q 'Normalization.C(true, {{2061, 411}, {0, 0}, {3280, 391}}, 2.69e16)'
// Original author - Pawel Plesniak

#include <limits>

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void plot(std::vector<std::vector<double>> capturedMuons, std::vector<double> expectedMuonCaptures, const bool highResolution, const std::string normalization_type) {
    /*
        Description
            Generates the normalization plot. If an entry is null, it is not included in th plot

        Arguments
            capturedMuons - as documented in CountMuCapPerMeasuredPhoton
            expectedMuonCaptures - expected number of captured muons and its uncertainty from the number of POTs
            highResolution - controls whether the plot is generated in low resolution or high resolution
            normalization_type - type of normalization used in the plot

        Variables
            lineWidth - thickness of the normalization estimate lines on the plots
            plotColors - colours of the normalization estimates and associated errors
            errorBoundOpacity - opacity of the normalization uncertainties in the plots
            plotTitle - plot title and axis labels
            plotFileName - name of the output file contianing the plot
            i - iterator
            px - number of x pixels for canvases
            py - number of y pixels for canvases
            c - TCanvas on which the plot is generated
            lowResAxisTextSize - updated text size for low resolution plots
            muonCaptureCount - number of caputred muons estimated from nSignalPhotons
            muonCaptureUncertainty - uncertainty of number of caputred muons estimated from nSignalPhotons
            normalizationSource - vector of strings containing the sources of the normalization estimates
            uncertaintyQuadrature - quadrature sum of uncertainties, used as a temporary variable before being assigned to muonCaptureUncertainty
            N - number of data points to plot for the estimate of the number of captured muons
            yMin - minimum value of y range
            yMax - maximum value of y range
            muonCaptureCountLines - vector of the TLines that define the expected number of captured muons
            errorX - x coordinates of the uncertainty error bound associated with the measurement
            errorY - y coordinates of the uncertainty error bound associated with the measurement
            errorBounds - vector of TGraphs that define the uncertainty bounds of the associated measurements
            muonCountMinusUncertainty - normalization count minus the associated uncertainty
            muonCountPlusUncertainty - normalization count plus the associated uncertainty
            xMin - smallest value to display on the normalization plot
            xMax - largest value to display on the normalization plot
            h - TGraph required to generate an axis range for the lines and error bounds
            lx1 - defines an x coordinate for the legend
            ly1 - defines a y coordinate for the legend
            lx2 - defines the other x coordinate for the legend
            ly2 - defines the other y coordinate for the legend
    */

    // Define the plot formatting
    const double lineWidth = 2;
    std::vector<int> plotColors = {kBlack, kRed, kBlue, kGreen + 2};
    const double errorBoundOpacity = 0.3;
    const std::string plotTitle = "Normalized muon capture count; Captured muon count;";
    const std::string plotFileName = std::string("MuonCaptures.") + normalization_type + "." + (highResolution ? "high" : "low") + ".png";

    // Set up iterator
    int i = 0;

    // Set up the canvas
    int px = (highResolution ? 1500 : 750), py = (highResolution ? 1000 : 500);
    TCanvas* c = new TCanvas("c", "c", px, py);

    // Format the canvas if the text is not high resolution
    double lowResAxisTextSize = 0.06;
    gStyle->SetTitleFontSize(lowResAxisTextSize);
    if (!highResolution) {
        gStyle->SetTitleAlign(33);
        gStyle->SetTitleX(0.75);
        c->SetBottomMargin(0.14);
        c->SetLeftMargin(0.175);
        c->Update();
    };

    // Determine which sources of noramlizations to plot, generate the averages for the plot
    std::vector<double> muonCaptureCount, muonCaptureUncertainty;
    std::vector<std::string> normalizationSource;
    double uncertaintyQuadrature = 0.0;

    if (expectedMuonCaptures[0] > std::numeric_limits<double>::epsilon()) {
        normalizationSource.push_back("POT");
        muonCaptureCount.push_back(expectedMuonCaptures[0]);
        muonCaptureUncertainty.push_back(expectedMuonCaptures[1]);
    };

    if (capturedMuons[0][0] > std::numeric_limits<double>::epsilon()) {
        normalizationSource.push_back("347 keV signal");
        muonCaptureCount.push_back(capturedMuons[0][0]);
        muonCaptureUncertainty.push_back(capturedMuons[0][1]);
    };
    if (capturedMuons[1][0] > std::numeric_limits<double>::epsilon()) {
        normalizationSource.push_back("844 keV signal");
        muonCaptureCount.push_back(capturedMuons[1][0]);
        muonCaptureUncertainty.push_back(capturedMuons[1][1]);
    };
    if (capturedMuons[2][0] > std::numeric_limits<double>::epsilon()) {
        normalizationSource.push_back("1809 keV signal");
        muonCaptureCount.push_back(capturedMuons[2][0]);
        muonCaptureUncertainty.push_back(capturedMuons[2][1]);
    };
    const int N = normalizationSource.size();
    if (N == 0)
        Fatal("plot", "There is no data to be plotted, exiting");
    const double yMin = 0, yMax = N;

    // Set up the average lines
    std::vector<TLine*> muonCaptureCountLines;
    for (i = 0; i < N; i++)
        muonCaptureCountLines.push_back(new TLine(muonCaptureCount[i], yMin, muonCaptureCount[i], yMax));

    // Set up the error bounds
    std::vector<double> errorX = {0, 0, 0, 0}, errorY = {yMin, yMin, yMax, yMax};
    std::vector<TGraph*> errorBounds;
    for (i = 0; i < N; i++) {
        errorX = {muonCaptureCount[i] - muonCaptureUncertainty[i], muonCaptureCount[i] + muonCaptureUncertainty[i], muonCaptureCount[i] + muonCaptureUncertainty[i], muonCaptureCount[i] - muonCaptureUncertainty[i]};
        errorBounds.emplace_back(new TGraph(4, errorX.data(), errorY.data()));
    };

    // Set up the parameters to draw a TH1F - this is required to get the axis
    std::vector<double> muonCountMinusUncertainty(N), muonCountPlusUncertainty(N);
    std::transform(muonCaptureCount.begin(), muonCaptureCount.end(), muonCaptureUncertainty.begin(), muonCountMinusUncertainty.begin(), std::minus<double>());
    std::transform(muonCaptureCount.begin(), muonCaptureCount.end(), muonCaptureUncertainty.begin(), muonCountPlusUncertainty.begin(),  std::plus<double>());
    const double xMin = *std::min_element(muonCountMinusUncertainty.begin(), muonCountMinusUncertainty.end()) * 0.5;
    const double xMax = *std::max_element(muonCountPlusUncertainty.begin(),  muonCountPlusUncertainty.end()) * 1.05;

    // Draw the TH1F
    TH1F *h = new TH1F("h", plotTitle.c_str(), 100, xMin, xMax);
    h->SetMinimum(yMin); // Y-axis range
    h->SetMaximum(yMax);
    h->Draw("HIST"); // Draw the histogram as an axis background
    h->SetStats(0);
    h->GetYaxis()->SetLabelSize(0);
    h->GetYaxis()->SetTickLength(0);


    // Set up the legend
    const double lx1 = (highResolution ? 0.7 : 0.6);
    const double ly1 = lx1;
    const double lx2 = 0.9;
    const double ly2 = 0.9;
    TLegend *legend = new TLegend(lx1, ly1, lx2, ly2, "Data");
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->SetHeader("Normalization source", "C");

    // Draw the muon capture normalizations and associated errors
    for (i = 0; i < N; i++) {
        muonCaptureCountLines[i]->Draw();
        muonCaptureCountLines[i]->SetLineColor(plotColors[i]);
        muonCaptureCountLines[i]->SetLineWidth(lineWidth);
        errorBounds[i]->Draw("F SAME");
        errorBounds[i]->SetFillColorAlpha(plotColors[i], errorBoundOpacity);
        legend->AddEntry(errorBounds[i], normalizationSource[i].c_str(), "F")->SetTextColor(plotColors[i]);
    };

    // Update and save the canvas
    legend->Draw();
    c->Update();
    c->SaveAs(plotFileName.c_str());

    // Cleanup
    h->Delete();
    c->Close();

    std::cout << "Finished" << std::endl;
    return;
};

std::string doubleToStringScientific(double value, int sigFigs) {
    /*
        Description
            Converts a double to a string with a defined number of significant figures

        Arguments
            value - value to convert to string
            sigFigs - number of significant figures to use

        Variables
            out - output string stream used to convert 'value' to a std::string
    */
    std::ostringstream out;
    out << std::scientific << std::setprecision(sigFigs - 1) << value;
    return out.str();
};

std::string doubleToStringFixed(double value, int sigFigs) {
    /*
        Description
            Converts a double to a string with a defined number of significant figures

        Arguments
            value - value to convert to string
            sigFigs - number of significant figures to use

        Variables
            out - output string stream used to convert 'value' to a std::string
    */
    std::ostringstream out;
    out << std::fixed << std::setprecision(sigFigs - 1) << value;
    return out.str();
};

void Normalization(const std::string normalization_type, std::vector<std::vector<double>> nSignalPhotons = {{0, 0}, {0, 0}, {0, 0}}, const double nPOTs = 2.69e16, bool makePlot = true) {
    /*
        Description
            Generates the table describing the number of measured photons per muon capture, and the table if
                - makePlot is true
                - at least one of the signal photon counters or POT count is non-zero
            File is named as
                MuonCaptures.<high/low>.png
            such that <high/low> defines high and low resolutions respectively
            The parameters used to generate the table are defined within this file. Update them with new values once the simulation work has been done

        Arguments
            makePlot - determines whether the plot is made or not
            nSignalPhotons - number of signal photons and their associated uncertainty, defined as {{n347, u347}, {n844, u844}, {n1809, u1809}}
            nPOTs - number of POTs

        Variables
            correctionNameColumnWidth - width of the column containing the correction name
            signal347ColumnWidth - width of the 347 signal column
            signal844ColumnWidth - width of the 844 signal column
            signal1809ColumnWidth - width of the 1809 signal column
            fullWidth - full table column width
            nSF - number of significant figures to display
            order - order in which the signal photons are presented in the table, and the order in which their correction parameters are defined
            nOrder - number of signal photons to analyze
            i - iterator variable
            j - iterator variable
            correctionFactorNames - names of the correction factors defined in the table, stored in presentation order
            nCorrectionFactors - the number of correction factors
            The following variables contain the correction factors and their associated uncertainties as {{c347, u347}, {c844, u844}, {c1809, u1809}}, with cX being the correction factor of signal X, with associated uncertainty uX
                pFinalState - probability of a final state given a muon stop
                absorberAcceptance - probability that the absorber does not shift the photon energy out of the measurement window
                detectorAcceptance - photopeak efficiency
                path attenuation - attenuation of beamline elements
                geantRateCorrection - correction from the simulated rate provided by geant to the experimentally measured rates
                signalInEnergyWindow - probability that the emitted signal has been measured within the selected measurement window
                clippingFactor - dead time correction due to beam intensity fluctuations
                geometricAcceptance - acceptance from the ST to the SSC aperture
                timeCutAcceptance - probability that the signal photon arrives in the selected time window
            measuredPhotonPerMuonCapture - probability and uncertainty of measuring a signal photon given a muon capture as {{m347, u347}, {m844, u844}, {m1809, u1809}}, with mX being the number of measured signal photons with energy X per muon capture, with associated uncertainty uX
            correctionFactors - vector of all the correction factors
            correctionFactor - iterator for correctionFactors
            muonCapturePerSignalPhoton - number and uncertainty of captured muons per signal photon, calculated as the inverse of measuredPhotonPerMuonCapture as {{n347, u347}, {n844, u844}, {n1809, u1809}}, with nX being the number of captured muon per signal photon X, with associated uncertainty uX
            capturedMuons - number of measured signal photons expected using nSignalPhotons
            signalColumnWidths - vector containing the signal column widths
            boolValues - vector of the available boolean values
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Sanity check
    if (nSignalPhotons.size() != 3)
        Fatal("CountMuCapPerMeasuredPhoton", "Incorrect format of nSignalPhotons, expected as {{n347, u347}, {n844, u844}, {n1809, u1809}}");
    if (nPOTs <= 0) {
        std::cerr << "Error: Number of POTs must be positive." << std::endl;
        return;
    };
    std::vector<std::string> normalization_types{"all", "VD90", "VD101", "DetectorMC", "DetectorReco"};
    if (std::find(normalization_types.begin(), normalization_types.end(), normalization_type) == normalization_types.end()) {
        std::cout << "Provided normalization type: " << normalization_type << std::endl;
        std::cerr << "Error: Normalization type must be one of 'all', 'VD90', 'VD101', 'DetectorMC', or 'DetectorReco'." << std::endl;
        return;
    };

    // Declare iterator variables
    int i = 0, j = 0;

    // Define the signal order
    std::vector<std::string> order = {"347 keV", "844 keV", "1809 keV"};
    const int nOrder = order.size();

    // Generate the expected number of muon captures and its uncertainty from the number of POTs
    // Number of events from POT.fcl
    const double nMDC2020_POTs = 2e8;
    const double nMDC2020_MuBeamCat = 869305.0;
    const double uMDC2020_MuBeamCat = std::sqrt(nMDC2020_MuBeamCat);
    // Number of events from MuBeamResampler.fcl
    const double nMDC2020_MuonBeamResampler = 4e9;
    const double nMDC2020_TargetStopsCat = 1432353.0 * 1000.0; // Factor 1000 is associated with the fact that MuBeamResampler has a prescale of 1000 applied, which is not documented :(
    const double uMDC2020_TargetStopsCat = 1000.0 * std::sqrt(1432353.0); // Propagate the error correctly
    // Resampling factors for MuBeamResampler.fcl
    const double R_MuBeamResampler = nMDC2020_MuonBeamResampler / nMDC2020_MuBeamCat;
    const double uR_MuBeamResampler = R_MuBeamResampler * uMDC2020_MuBeamCat / (nMDC2020_MuBeamCat * nMDC2020_MuBeamCat);
    // Number of POTs in sim.mu2e.TargetStopsCat.MDC2020p.art
    const double nMDC2020_TargetStopsCat_POTs = nMDC2020_POTs * R_MuBeamResampler;
    const double uMDC2020_TargetStopsCat_POTs = nMDC2020_TargetStopsCat_POTs * uR_MuBeamResampler / R_MuBeamResampler;
    // Probability of muon stop per POT
    const double pMDC2020_MuonStop_from_POT = nMDC2020_TargetStopsCat / nMDC2020_TargetStopsCat_POTs;
    const double uMDC2020_MuonStop_from_POT = pMDC2020_MuonStop_from_POT * std::sqrt(std::pow(uMDC2020_TargetStopsCat / nMDC2020_TargetStopsCat, 2) + std::pow(uMDC2020_TargetStopsCat_POTs / nMDC2020_TargetStopsCat_POTs, 2));
    // Probability of muon capture given a muon stop
    const double pMuonCapture = 0.61;
    const double uMuonCapture = 0.001;
    // Probability of muon capture per POT
    const double pMDC2020_MuonCapture_from_POT = pMDC2020_MuonStop_from_POT * pMuonCapture;
    const double uMDC2020_MuonCapture_from_POT = pMDC2020_MuonCapture_from_POT * std::sqrt(std::pow(uMDC2020_MuonStop_from_POT / pMDC2020_MuonStop_from_POT, 2) + std::pow(uMuonCapture / pMuonCapture, 2));
    // Expected number of captured muons in the simulation
    const double nExpectedMuonCaptures = nPOTs * pMDC2020_MuonCapture_from_POT;
    const double uExpectedMuonCaptures = nExpectedMuonCaptures * (uMDC2020_MuonCapture_from_POT / pMDC2020_MuonCapture_from_POT);
    // Construct the same format as the other parameters for easier plotting
    const std::vector<double> normalizationFromPOT = {nExpectedMuonCaptures, uExpectedMuonCaptures};

    // Define the parameters that contribute to the number of measured signal photons per muon capture, defined as {value, uncertainty} for each signal photon

    // Define the correction factors and their associated uncertianties
    //                                                          347 corr    uncert      844 corr    uncert      1809 corr   uncert
    std::vector<std::vector<double>> pFinalState            = {{1.31,       0.013},    {0.093,      0.007},     {0.51,      0.05}       };
    std::vector<std::vector<double>> geantRateCorrection    = {{1.0,        0},        {2.590e-1,   0},         {1.0,       0}          };
    std::vector<std::vector<double>> geometricAcceptance    = {{3.25e-9,    0},        {3.25e-9,    0},         {3.25e-9,   0}          };
    std::vector<std::vector<double>> pathAttenuation        = {{1,          0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> absorberAcceptance     = {{0.87,       0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> energyWindowAcceptance = {{0.67,       0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> timeCutAcceptance      = {{0.9976,     0},        {0.6298,     0},         {0.68,      0}          };
    std::vector<std::vector<double>> detectorAcceptance     = {{0.628,      1.528e-4}, {0.288,      1.432e-4},  {0.179,     1.212e-4}   };
    std::vector<std::vector<double>> clippingFactor         = {{0.85,       0},        {1,          0},         {0.85,      0}          };

    // Store all the correction factors in a single variable
    std::vector<std::vector<std::vector<double>>> correctionFactors = {pFinalState, geantRateCorrection, geometricAcceptance, pathAttenuation, absorberAcceptance, energyWindowAcceptance, timeCutAcceptance};

    // Define the names of the correction factors and store all the correction factors in a single variable for easier iteration when applying the correction factors and printing the table
    // This contains the names of all the relevant correction factors for all the normalization types, if others are needed they are appended
    std::vector<std::string> correctionFactorNames = {
        "Probability of final state",
        "GEANT rate correction",
        "Geometric acceptance",
        "Path attenuation",
        "Absorber acceptance",
        "Energy window acceptance",
        "Time cut acceptance"
    };

    // If the normalization is at the detector using MC truth, we need to include the detector acceptance correction factor
    // If the normalization is at the detector using reconstructed data, we need to include both the detector acceptance and clipping factor correction factors

    if (normalization_type != "VD90") {
        correctionFactors.push_back(detectorAcceptance);
        correctionFactorNames.push_back("Detector acceptance");
        if (normalization_type == "DetectorReco") {
            correctionFactors.push_back(clippingFactor);
            correctionFactorNames.push_back("Clipping factor");
        };
    }
    else
        geometricAcceptance    = {{7.74e-6,    0},        {7.74e-6,    0},         {7.74e-6,   0}          };
    const int nCorrectionFactors = correctionFactorNames.size();

    // Construct the variable to store the number of photons per muon capture and its uncertainty
    //                                                                  347 uncert  844     uncert  1809    uncert
    std::vector<std::vector<double>> measuredPhotonPerMuonCapture   = {{1,  0},     {1,     0},     {1,     0}};

    // Incorporate the effect of the correction factors into the number of measured photons per muon capture
    for (std::vector<std::vector<double>> correctionFactor : correctionFactors) {
        for (i = 0; i < nOrder; i++) {
            // Update the rate
            measuredPhotonPerMuonCapture[i][0] *= correctionFactor[i][0];

            // Update the quadratic sum of the uncertainty
            if (correctionFactor[i][0] < std::numeric_limits<double>::epsilon())
                Fatal("CountMuCapPerMeasuredPhoton", "Declared correction factor is zero, this would mean we don't get any signal. Did you mean to set this to unity? Fix this!");
            measuredPhotonPerMuonCapture[i][1] += std::pow(correctionFactor[i][1]/correctionFactor[i][0], 2);
        };
    };

    // Finalize the uncertainty
    for (i = 0; i < nOrder; i++)
        measuredPhotonPerMuonCapture[i][1] = measuredPhotonPerMuonCapture[i][0] * std::sqrt(measuredPhotonPerMuonCapture[i][1]);

    // Calculate the number of muon captures per measured signal photon and its errors
    std::vector<std::vector<double>> muonCapturePerSignalPhoton = {{0, 0}, {0, 0}, {0, 0}};
    for (i = 0; i < nOrder; i++) {
        muonCapturePerSignalPhoton[i][0] = 1.0/measuredPhotonPerMuonCapture[i][0];
        muonCapturePerSignalPhoton[i][1] = std::pow(muonCapturePerSignalPhoton[i][0], 2) * measuredPhotonPerMuonCapture[i][1];
    };

    // Determine the normalization
    std::vector<std::vector<double>> normalizationFromSignalPhotons = {{0, 0}, {0, 0}, {0, 0}};
    for (i = 0; i < nOrder; i++) {
        normalizationFromSignalPhotons[i][0] = nSignalPhotons[i][0] * muonCapturePerSignalPhoton[i][0];
        if (nSignalPhotons[i][0] < std::numeric_limits<double>::epsilon() || muonCapturePerSignalPhoton[i][0] < std::numeric_limits<double>::epsilon())
            normalizationFromSignalPhotons[i][1] = 0.0;
        else
            normalizationFromSignalPhotons[i][1] = normalizationFromSignalPhotons[i][0] * std::sqrt(std::pow(nSignalPhotons[i][1]/nSignalPhotons[i][0], 2) + std::pow(muonCapturePerSignalPhoton[i][1]/muonCapturePerSignalPhoton[i][0], 2));
    };

    // Determine the expected number of signal photons from the number of POTs and the expected number of captured muons from POTs
    std::vector<std::vector<double>> expectedPhotonCountFromPOT = {{0, 0}, {0, 0}, {0, 0}};
    for (i = 0; i < nOrder; i++) {
        expectedPhotonCountFromPOT[i][0] = nExpectedMuonCaptures * measuredPhotonPerMuonCapture[i][0];
        if (nExpectedMuonCaptures < std::numeric_limits<double>::epsilon() || measuredPhotonPerMuonCapture[i][0] < std::numeric_limits<double>::epsilon())
            expectedPhotonCountFromPOT[i][1] = 0.0;
        else
            expectedPhotonCountFromPOT[i][1] = expectedPhotonCountFromPOT[i][0] * std::sqrt(std::pow(nExpectedMuonCaptures/nExpectedMuonCaptures, 2) + std::pow(measuredPhotonPerMuonCapture[i][1]/measuredPhotonPerMuonCapture[i][0], 2));
    };

    // Determine the ratio of the expected number of signal photons from POTs to the determined number of signal photons from the measured signal photons
    std::vector<std::vector<double>> measuredToExpectedNormalizationRatio = {{0, 0}, {0, 0}, {0, 0}};
    for (i = 0; i < nOrder; i++) {
        measuredToExpectedNormalizationRatio[i][0] = nSignalPhotons[i][0] / expectedPhotonCountFromPOT[i][0];
        if (expectedPhotonCountFromPOT[i][0] < std::numeric_limits<double>::epsilon() || nSignalPhotons[i][0] < std::numeric_limits<double>::epsilon())
            measuredToExpectedNormalizationRatio[i][1] = 0.0;
        else 
            measuredToExpectedNormalizationRatio[i][1] = measuredToExpectedNormalizationRatio[i][0] * std::sqrt(std::pow(expectedPhotonCountFromPOT[i][1]/expectedPhotonCountFromPOT[i][0], 2) + std::pow(nSignalPhotons[i][1]/nSignalPhotons[i][0], 2));
    };

    // Define the table formatting
    const int correctionNameColumnWidth = 40, signal347ColumnWidth = 30, signal844ColumnWidth = 30, signal1809ColumnWidth = 30, fullWidth = correctionNameColumnWidth + signal347ColumnWidth + signal844ColumnWidth + signal1809ColumnWidth;
    const int nSF = 4;

    // Print the title line and rules
    const std::vector<int> signalColumnWidths = {signal347ColumnWidth, signal844ColumnWidth, signal1809ColumnWidth};
    std::string tableTitle = "Normalized muon capture count per measured signal photon";
    if (normalization_type == "VD90")
        tableTitle += " at VD90";
    else if (normalization_type == "DetectorMC")
        tableTitle += " at the detector using MC truth";
    else if (normalization_type == "DetectorReco")
        tableTitle += " at the detector using reconstructed data";

    std::cout << std::endl; // Buffer line
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::string((fullWidth - tableTitle.size())/2, ' ') << tableTitle << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Title line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Correction factor";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << order[i];
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the correction factors and bottom rule
    for (i = 0; i < nCorrectionFactors; i++) {
        std::cout << std::setw(correctionNameColumnWidth) << std::left << correctionFactorNames[i];
        for (j = 0; j < nOrder; j++)
            std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringScientific(correctionFactors[i][j][0], nSF) + " ± " + doubleToStringScientific(correctionFactors[i][j][1], nSF);
        std::cout << std::endl;
    };
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the normalized quantities and bottom rule
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "N signal photons per captured muon";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(measuredPhotonPerMuonCapture[i][0], nSF) + " ± " + doubleToStringScientific(measuredPhotonPerMuonCapture[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "N captured muons per signal photon";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(muonCapturePerSignalPhoton[i][0], nSF) + " ± " + doubleToStringScientific(muonCapturePerSignalPhoton[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // End line

    // Print the normalization estimate
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Measured signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(nSignalPhotons[i][0], nSF) + " ± " + doubleToStringFixed(nSignalPhotons[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(expectedPhotonCountFromPOT[i][0], nSF) + " ± " + doubleToStringFixed(expectedPhotonCountFromPOT[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Measured to expected ratio";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(measuredToExpectedNormalizationRatio[i][0], nSF) + " ± " + doubleToStringFixed(measuredToExpectedNormalizationRatio[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Normalized captured muons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(normalizationFromSignalPhotons[i][0], nSF) + " ± " + doubleToStringScientific(normalizationFromSignalPhotons[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // End line
    std::cout << std::endl; // Buffer line

    // Generate plot if relevant to do so
    if (makePlot) {
        if (std::abs(nSignalPhotons[0][0]) < std::numeric_limits<double>::epsilon() && std::abs(nSignalPhotons[1][0]) < std::numeric_limits<double>::epsilon() && std::abs(nSignalPhotons[2][0]) < std::numeric_limits<double>::epsilon() && nPOTs == 0)
            Fatal("CountMuCapPerMeasuredPhoton", "Plot has been requested, but the number of signal photons and POTs are all zero, skipping.");

        // Generate the plots
        std::vector<bool> boolValues = {true, false};
        for (bool highResolution : boolValues)
            plot(normalizationFromSignalPhotons, normalizationFromPOT, highResolution, normalization_type);
    };

    return;
};
