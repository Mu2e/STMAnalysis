// Calculates the expected number of signal photons assuming the MC truth normalizations are correct. The expected number of signal photons is calculated from the probability of final state for each signal photon energy, and the number of POTs simulated including resampling. 
// Usage example - $ root -l -q 'ExpectedPhotonCount.C({{1.679e10, 3.02e9}, {0.0, 0.0}, {1.583e11, 2.28e10}}, {{4.275e11, 8.56e10}, {0.0, 0.0}, {7.660e13, 1.180e13}})'
// Original author - Pawel Plesniak


void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
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

void ExpectedPhotonCount(double nPOTs = 2.69e16) {
    /*
        Description
            Generates the table describing the expected number of signal photons for each energy

        Arguments
            nPOTs - number of Protons on Target simulated including resampling

        Variables
            nSF - number of significant figures to use in the table
            correctionNameColumnWidth - width of the correction name column
            signal347ColumnWidth - width of the 347 keV signal column
            signal844ColumnWidth - width of the 844 keV signal column
            signal1809ColumnWidth - width of the 1809 keV signal column
            fullWidth - total width of the table
            i, j - iterators
            ratioMCToDetector - vector of vectors containing the ratio of MC truth normalization to detector response normalization and associated uncertainties for 347, 844, and 1809 keV
    */


    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Define the table formatting
    const int correctionNameColumnWidth = 40, signal347ColumnWidth = 30, signal844ColumnWidth = 30, signal1809ColumnWidth = 30, fullWidth = correctionNameColumnWidth + signal347ColumnWidth + signal844ColumnWidth + signal1809ColumnWidth;
    const int nSF = 4;

    // Define the signal order
    std::vector<std::string> order = {"347", "844", "1809"};
    const int nOrder = order.size();

    // Declare iterator variables
    int i = 0, j = 0;

    // Define the parameters that contribute to the number of measured signal photons per muon capture, defined as {value, uncertainty} for each signal photon
    std::vector<std::string> correctionFactorNames = {
        // "Probability of final state",
        "Absorber acceptance",
        "Detector acceptance",
        "Path attenuation",
        "GEANT rate correction",
        "Energy window acceptance",
        "Clipping factor",
        "Geometric acceptance",
        "Time cut acceptance"
    };

    // Define the correction factors and their associated uncertianties
    const int nCorrectionFactors = correctionFactorNames.size();
    //                                                          347 corr    uncert      844 corr    uncert      1809 corr   uncert
    std::vector<std::vector<double>> pFinalState            = {{1.31,       0.013},    {0.093,      0.007},     {0.51,      0.05}       };
    std::vector<std::vector<double>> absorberAcceptance     = {{0.87,       0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> detectorAcceptance     = {{0.628,      1.528e-4}, {0.288,      1.432e-4},  {0.179,     1.212e-4}   };
    std::vector<std::vector<double>> pathAttenuation        = {{1,          0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> geantRateCorrection    = {{1,          0},        {0.259,      0},         {1.0,       0}          };
    std::vector<std::vector<double>> signalInEnergyWindow   = {{0.67,       0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> clippingFactor         = {{0.85,       0},        {1,          0},         {0.85,      0}          };
    std::vector<std::vector<double>> geometricAcceptance    = {{3.25e-9,    0},        {3.25e-9,    0},         {3.25e-9,   0}          };
    std::vector<std::vector<double>> timeCutAcceptance      = {{0.9976,     0},        {0.6298,     0},         {0.68,      0}          };


    // Calculate the ratio of MC truth normalization to expected number of captured muons and its errors
    std::vector<std::vector<double>> expectation = {{0, 0}, {0, 0}, {0, 0}};
    const double pMuonStopMDC2020 = 1432535.0 / (2e8 * (4e8 / 869305)); // Based on the MDC2020 workflow
    const double uMuonStopMDC2020 = std::sqrt(1432535) / (2e8 * (4e8 / 869305)); // Based on the MDC2020 workflow
    const double pMuonCapture = 0.61;
    const double uMuonCapture = 0.001;
    const double nExpectedMuonCaptures = nPOTs * pMuonStopMDC2020 * pMuonCapture;
    const double uExpectedMuonCaptures = nExpectedMuonCaptures * std::sqrt(std::pow(uMuonStopMDC2020/pMuonStopMDC2020, 2) + std::pow(uMuonCapture/pMuonCapture, 2));

    for (int i = 0; i < 3; i++) {
        expectation[i][0] = pFinalState[i][0] * nExpectedMuonCaptures;
        expectation[i][1] = expectation[i][0] * std::sqrt(std::pow(pFinalState[i][1]/pFinalState[i][0], 2) + std::pow(uExpectedMuonCaptures/nExpectedMuonCaptures, 2));
    };

    // Store all the correction factors in a single variable
    std::vector<std::vector<std::vector<double>>> correctionFactors = {absorberAcceptance, detectorAcceptance, pathAttenuation, geantRateCorrection, signalInEnergyWindow, clippingFactor, geometricAcceptance, timeCutAcceptance};

    // Incorporate the effect of the correction factors into the number of measured signal photons per muon capture
    std::vector<std::vector<double>> signalPhotonCount(expectation.begin(), expectation.end()); // Initialize the number of measured signal photons per muon capture to the expected number of signal photons per muon capture before corrections
    for (i = 0; i < nOrder; i++) {
        signalPhotonCount[i][1] = std::pow(expectation[i][1]/expectation[i][0], 2); // Initialize the quadratic sum of the relative uncertainty to the relative uncertainty on the expected number of signal photons per muon capture
    };

    // Incorporate the effect of the correction factors into the number of measured photons per muon capture
    for (std::vector<std::vector<double>> correctionFactor : correctionFactors) {
        for (i = 0; i < nOrder; i++) {
            // Update the rate
            signalPhotonCount[i][0] *= correctionFactor[i][0];

            // Update the quadratic sum of the uncertainty
            if (correctionFactor[i][0] < std::numeric_limits<double>::epsilon())
                Fatal("CountMuCapPerMeasuredPhoton", "Declared correction factor is zero, this would mean we don't get any signal. Did you mean to set this to unity? Fix this!");
            signalPhotonCount[i][1] += std::pow(correctionFactor[i][1]/correctionFactor[i][0], 2);
        };
    };

    // Finalize the uncertainty
    for (i = 0; i < nOrder; i++)
        signalPhotonCount[i][1] = signalPhotonCount[i][0] * std::sqrt(signalPhotonCount[i][1]);


    // Print the title line and rules
    const std::vector<int> signalColumnWidths = {signal347ColumnWidth, signal844ColumnWidth, signal1809ColumnWidth};
    std::string tableTitle = "Expected signal photon count";
    std::cout << "Probability of a POT producing a captured muon: " << pMuonStopMDC2020 << " ± " << uMuonStopMDC2020 << std::endl;
    std::cout << std::endl; // Buffer line
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::string((fullWidth - tableTitle.size())/2, ' ') << tableTitle << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Title line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "nPOTs";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(nPOTs, nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "nCapturedMuons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(nExpectedMuonCaptures, nSF) + " ± " + doubleToStringScientific(uExpectedMuonCaptures, nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected generated signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(expectation[i][0], nSF) + " ± " + doubleToStringScientific(expectation[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the correction factors and bottom rule
    for (i = 0; i < nCorrectionFactors; i++) {
        std::cout << std::setw(correctionNameColumnWidth) << std::left << correctionFactorNames[i];
        if (correctionFactorNames[i] == "Geometric acceptance") {
            for (j = 0; j < nOrder; j++)
                std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringScientific(correctionFactors[i][j][0], nSF) + " ± " + doubleToStringScientific(correctionFactors[i][j][1], nSF);
        }
        else {
            for (j = 0; j < nOrder; j++)
                std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringFixed(correctionFactors[i][j][0], nSF) + " ± " + doubleToStringFixed(correctionFactors[i][j][1], nSF);
        }
        std::cout << std::endl;
    };
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the normalized quantities and bottom rule
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(signalPhotonCount[i][0], nSF) + " ± " + doubleToStringScientific(signalPhotonCount[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // End line
    std::cout << std::endl; // Buffer line
};