// Calculates the expected number of signal photons assuming the correction factors are correct. The expected number of signal photons is calculated from the probability of final state for each signal photon energy, and the number of POTs simulated including resampling. 
// Usage example - $ root -l -q 'ExpectedPhotonCount.C'
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

void tabulate(const double nPOTs, const std::string location, bool &printSummary) {
    /*
        Description
            Generates the table describing the expected number of signal photons for each energy

        Arguments
            nPOTs - number of Protons on Target simulated including resampling
            location - location of the normalization. One of "VD90", "DetectorMC", or "DetectorReco"
            printSummary - whether to print the summary information about the expected number of captured muons in the simulation

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

    // Declare iterator variables
    int i = 0, j = 0;

    // Define the signal order
    std::vector<std::string> order = {"347 keV", "844 keV", "1809 keV"};
    const int nOrder = order.size();

    // Define the parameters that contribute to the number of measured signal photons per muon capture, defined as {value, uncertainty} for each signal photon
    //                                                          347 corr    uncert      844 corr    uncert      1809 corr   uncert
    std::vector<std::vector<double>> pFinalState            = {{1.31,       0.013},    {0.093,      0.007},     {0.51,      0.05}       };

    // Define the correction factors and their associated uncertianties
    //                                                          347 corr    uncert      844 corr    uncert      1809 corr   uncert
    std::vector<std::vector<double>> geometricAcceptance    = {{3.25e-9,    0},        {3.25e-9,    0},         {3.25e-9,   0}          };
    std::vector<std::vector<double>> pathAttenuation        = {{0.774,      0.008},    {0.763,      0.018},     {0.777,     0.012}      };
    std::vector<std::vector<double>> absorberAcceptance     = {{0.634,      0.033},    {0.835,      0.025},     {0.952,     0.048}      };
    std::vector<std::vector<double>> energyWindowAcceptance = {{0.67,       0},        {1,          0},         {1,         0}          };
    std::vector<std::vector<double>> timeCutAcceptance      = {{0.9976,     0},        {0.6298,     0},         {0.68,      0}          };
    std::vector<std::vector<double>> detectorAcceptance     = {{0.628,      1.528e-4}, {0.288,      1.432e-4},  {0.179,     1.212e-4}   };
    std::vector<std::vector<double>> clippingFactor         = {{0.85,       0},        {1,          0},         {0.85,      0}          };

    // Correct the normalization if normalizing against VD101
    if (location == "VD101") {
        geometricAcceptance    = {{7.74e-6,    0},        {7.74e-6,    0},         {7.74e-6,   0}          };
    };

    // Store all the correction factors in a single variable
    std::vector<std::vector<std::vector<double>>> correctionFactors = {geometricAcceptance, pathAttenuation, absorberAcceptance, energyWindowAcceptance, timeCutAcceptance};

    // Define the names of the correction factors and store all the correction factors in a single variable for easier iteration when applying the correction factors and printing the table
    // This contains the names of all the relevant correction factors for all the normalization types, if others are needed they are appended
    std::vector<std::string> correctionFactorNames = {
        "Geometric acceptance",
        "Path attenuation",
        "Absorber acceptance",
        "Energy window acceptance",
        "Time cut acceptance"
    };

    // If the normalization is at the detector using MC truth, we need to include the detector acceptance correction factor
    // If the normalization is at the detector using reconstructed data, we need to include both the detector acceptance and clipping factor correction factors

    if ((location == "DetectorMC") || (location == "DetectorReco")) {
        correctionFactors.push_back(detectorAcceptance);
        correctionFactorNames.push_back("Detector acceptance");
        if (location == "DetectorReco") {
            correctionFactors.push_back(clippingFactor);
            correctionFactorNames.push_back("Clipping factor");
        };
    };
    const int nCorrectionFactors = correctionFactorNames.size();


    // Calculate the ratio of MC truth normalization to expected number of captured muons and its errors
    std::vector<std::vector<double>> expectation = {{0, 0}, {0, 0}, {0, 0}};
    // Number of events from POT.fcl
    const double nMDC2020_POTs = 2e8;
    const double nMDC2020_MuBeamCat = 869305.0;
    const double uMDC2020_MuBeamCat = std::sqrt(nMDC2020_MuBeamCat);
    // Number of events from MuBeamResampler.fcl
    const double nMDC2020_MuonBeamResampler = 4e9;
    const double nMDC2020_TargetStopsCat = 1432353.0 * 1000.0; // Factor 1000 is associated with the fact that MuBeamResampler has a prescale of 1000 applied, which is not documented :(
    const double uMDC2020_TargetStopsCat = 1000.0 * std::sqrt(1432353.0 ); // Propagate the error correctly
    // Resampling factors for MuBeamResampler.fcl
    const double R_MuBeamResampler = nMDC2020_MuonBeamResampler / nMDC2020_MuBeamCat;
    const double uR_MuBeamResampler = R_MuBeamResampler * (uMDC2020_MuBeamCat / nMDC2020_MuBeamCat);
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
    // Print summary information
    if (printSummary) {
        std::cout << "Probability of a POT producing a stopped muon: " << pMDC2020_MuonStop_from_POT << " +/- " << uMDC2020_MuonStop_from_POT << std::endl;
        std::cout << "Probability of a POT producing a captured muon: " << pMDC2020_MuonCapture_from_POT << " +/- " << uMDC2020_MuonCapture_from_POT << std::endl;
        printSummary = false; // Set to false so that the summary information is only printed once if all three normalizations are printed
    }

    // Calculate the expected number of signal photons per muon capture and its uncertainty before corrections
    for (i = 0; i < nOrder; i++) {
        expectation[i][0] = pFinalState[i][0] * nExpectedMuonCaptures;
        expectation[i][1] = expectation[i][0] * std::sqrt(std::pow(pFinalState[i][1]/pFinalState[i][0], 2) + std::pow(uExpectedMuonCaptures/nExpectedMuonCaptures, 2));
    };
    // Initialize the number of measured signal photons per muon capture to the expected number of signal photons per muon capture before corrections
    std::vector<std::vector<double>> signalPhotonCount(expectation.begin(), expectation.end()); 
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

    // Define the table formatting
    const int correctionNameColumnWidth = 40, signal347ColumnWidth = 30, signal844ColumnWidth = 30, signal1809ColumnWidth = 30, fullWidth = correctionNameColumnWidth + signal347ColumnWidth + signal844ColumnWidth + signal1809ColumnWidth;
    const std::vector<int> signalColumnWidths = {signal347ColumnWidth, signal844ColumnWidth, signal1809ColumnWidth};
    const int nSF = 4;

    // Define the table title based on the normalization location
    std::string tableTitle = "Expected signal photon count at ";
    if (location == "VD90")
        tableTitle += "VD90";
    else if (location == "DetectorMC")
        tableTitle += "detector using MC truth normalization";
    else if (location == "DetectorReco")
        tableTitle += "detector using reconstructed data normalization";
    tableTitle += " for " + doubleToStringScientific(nPOTs, nSF) + " POTs";

    // Print the title line and rules
    std::cout << std::endl; // Buffer line

    // Title and top rule
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::string((fullWidth - tableTitle.size())/2, ' ') << tableTitle << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Title line

    // Column titles and middle rule
    std::cout << std::setw(correctionNameColumnWidth) << std::left << " ";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << order[i];
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Title line

    // Print the expected number of signal photons per muon capture and the number of POTs simulated including resampling, and the number of expected muon captures in the simulation, and the middle rule
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "nPOTs";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(nPOTs, nSF); // Note the drift in the columns is due to C++ being wierd with +/- 
    std::cout << std::endl;

    std::cout << std::setw(correctionNameColumnWidth) << std::left << "nCapturedMuons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(nExpectedMuonCaptures, nSF) + " +/- " + doubleToStringScientific(uExpectedMuonCaptures, nSF);
    std::cout << std::endl;

    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Signal photon probability";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(pFinalState[i][0], nSF) + " +/- " + doubleToStringFixed(pFinalState[i][1], nSF);
    std::cout << std::endl;

    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected generated signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(expectation[i][0], nSF) + " +/- " + doubleToStringScientific(expectation[i][1], nSF);
    std::cout << std::endl;

    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the correction factors and bottom rule
    for (i = 0; i < nCorrectionFactors; i++) {
        std::cout << std::setw(correctionNameColumnWidth) << std::left << correctionFactorNames[i];
        if (correctionFactorNames[i] == "Geometric acceptance") {
            for (j = 0; j < nOrder; j++)
                std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringScientific(correctionFactors[i][j][0], 3) + " +/- " + doubleToStringScientific(correctionFactors[i][j][1], 3);
        }
        else {
            for (j = 0; j < nOrder; j++)
                std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringFixed(correctionFactors[i][j][0], nSF) + " +/- " + doubleToStringFixed(correctionFactors[i][j][1], nSF);
        }
        std::cout << std::endl;
        if (correctionFactorNames[i] == "Time cut acceptance")
            std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    };
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the normalized quantities and bottom rule
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Expected signal photons";
    for (i = 0; i < nOrder; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(signalPhotonCount[i][0], nSF) + " +/- " + doubleToStringFixed(signalPhotonCount[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // End line
    std::cout << std::endl; // Buffer line
};

void ExpectedPhotonCount(const double nPOTs = 2.69e16, const std::string location = "all", bool printSummary = true) {
    /*
        Description
            Generates the table describing the expected number of signal photons for each energy

        Arguments
            nPOTs - number of Protons on Target simulated including resampling
            location - location of the normalization. One of "all", "VD101", "VD90", "DetectorMC", or "DetectorReco"
            printSummary - whether to print the summary information about the expected number of captured muons in the simulation

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

    // Sanity check
    if (nPOTs <= 0) {
        std::cerr << "Error: Number of POTs must be positive." << std::endl;
        return;
    };
    std::vector<std::string> locationNames{"all", "VD101", "VD90", "DetectorMC", "DetectorReco"};
    if (std::find(locationNames.begin(), locationNames.end(), location) == locationNames.end()) {
        std::cout << "Provided location: " << location << std::endl;
        std::cerr << "Error: Location must be one of 'all', 'VD101', 'VD90', 'DetectorMC', or 'DetectorReco'." << std::endl;
        return;
    };

    bool summaryPrinted = printSummary; // Set the global variable to control whether the summary information is printed in the tabulate function

    // Generate the table
    if (location == "all") {
        std::vector<std::string> allLocationNames{"VD101", "VD90", "DetectorMC", "DetectorReco"};
        for (const std::string& locationName : allLocationNames) {
            std::cout << "Normalization location: " << locationName << std::endl;
            tabulate(nPOTs, locationName, summaryPrinted);
        };
    }
    else
        tabulate(nPOTs, location, summaryPrinted);
};