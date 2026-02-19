// Generates the table of ratio of any normalization type to the expected based on the POT count
// Usage example - $ root -l -q 'NormalizationComparisonToExpectedCount.C({{1.679e10, 3.02e9}, {0.0, 0.0}, {1.583e11, 2.28e10}}, 2.69e16)'
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

void NormalizationComparisonToExpectedCount(std::vector<std::vector<double>> normalizations, const double nPOTs = 2.69e16) {
    /*
        Description
            Generates the table describing the ratio of a normalization method to the expected number of captured muons

        Arguments
            normalizations - vector of vectors containing the normalizations and associated uncertainties for 347, 844, and 1809 keV
            nPOTs - number of Protons on Target simulated including resampling

        Variables
            nSF - number of significant figures to use in the table
            correctionNameColumnWidth - width of the correction name column
            signal347ColumnWidth - width of the 347 keV signal column
            signal844ColumnWidth - width of the 844 keV signal column
            signal1809ColumnWidth - width of the 1809 keV signal column
            fullWidth - total width of the table
            i, j - iterators
            expectation - vector of vectors containing the expected number of captured muons and associated uncertainties for 347, 844, and 1809 keV
            signalProbability - vector of vectors containing the probability of final state for 347, 844, and 1809 keV
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

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


    // Define the signal photon energy strings
    std::string e347  = "347 keV";
    std::string e844  = "844 keV";
    std::string e1809 = "1809 keV";

    // Define table formatting
    const int nSF = 4;
    const int correctionNameColumnWidth = 20, signalColumnWidth = 30;
    const int fullWidth = correctionNameColumnWidth + signalColumnWidth * 3;

    // Print the title line and rules
    std::string title = "Expected to normalization signal photon count comparison for " + doubleToStringScientific(nPOTs, nSF) + " POTs";
    std::cout << std::endl; // Buffer line
    std::cout << std::string(fullWidth, '=') << std::endl; // Title line
    std::cout << std::string((fullWidth - static_cast<int>(title.size())) / 2, ' ') << title << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // Title line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "" << std::setw(signalColumnWidth) << std::left << "347 keV" << std::setw(signalColumnWidth) << std::left << "844 keV" << std::setw(signalColumnWidth) << std::left << "1809 keV" << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected values";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringScientific(signalPhotonCount[j][0], nSF) + " ± " + doubleToStringScientific(signalPhotonCount[j][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Normalization";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringScientific(normalizations[j][0], nSF) + " ± " + doubleToStringScientific(normalizations[j][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Ratio";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringFixed(normalizations[j][0]/signalPhotonCount[j][0], nSF) + " ± " + doubleToStringFixed(std::sqrt(std::pow(normalizations[j][1]/normalizations[j][0], 2) + std::pow(signalPhotonCount[j][1]/signalPhotonCount[j][0], 2)), nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // End line
    std::cout << std::endl; // Buffer line
    return;
};