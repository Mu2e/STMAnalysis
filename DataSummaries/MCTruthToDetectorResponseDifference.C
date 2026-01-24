// Generates the table of ratio of MC truth normalization to detector response normalization differences
// Usage example - $ root -l -q 'MCTruthToDetectorResponseDifference.C({{1.679e10, 3.02e9}, {0.0, 0.0}, {1.583e11, 2.28e10}}, {{4.275e11, 8.56e10}, {0.0, 0.0}, {7.660e13, 1.180e13}})'
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

void MCTruthToDetectorResponseDifference(std::vector<std::vector<double>> mcTruthNormalizations, std::vector<std::vector<double>> detectorResponseNormalizations, double nPOTs = 2.69e16) {
    /*
        Description
            Generates the table describing the ratio of MC truth normalization to detector response normalization differences

        Arguments
            mcTruthNormalizations - vector of vectors containing the MC truth normalizations and associated uncertainties for 347, 844, and 1809 keV
            detectorResponseNormalizations - vector of vectors containing the detector response normalizations and associated uncertainties for 347, 844, and 1809 keV

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

    // Sanity check
    if (mcTruthNormalizations.size() != 3 || detectorResponseNormalizations.size() != 3)
        Fatal("MCTruthToDetectorResponseDifference", "Incorrect format of input normalizations, expected as {{value347, uncertainty347}, {value844, uncertainty844}, {value1809, uncertainty1809}}");
    if (mcTruthNormalizations == detectorResponseNormalizations)
        Fatal("MCTruthToDetectorResponseDifference", "Input normalizations are identical, no difference to calculate");
    std::vector<std::vector<double>> zeroVector = {{0, 0}, {0, 0}, {0, 0}};
    if (mcTruthNormalizations == zeroVector || detectorResponseNormalizations == zeroVector)
        Fatal("MCTruthToDetectorResponseDifference", "Input normalizations cannot be zero");

    // Calculate the ratio of MC truth normalization to detector response normalization and its errors
    std::vector<std::vector<double>> ratioMCToDetector = {{0, 0}, {0, 0}, {0, 0}};
    for (int i = 0; i < 3; i++) {
        ratioMCToDetector[i][0] = detectorResponseNormalizations[i][0] / mcTruthNormalizations[i][0];
        ratioMCToDetector[i][1] = ratioMCToDetector[i][0] * std::sqrt(std::pow(mcTruthNormalizations[i][1]/mcTruthNormalizations[i][0], 2) + std::pow(detectorResponseNormalizations[i][1]/detectorResponseNormalizations[i][0], 2));
    };


    // Define the signal photon energy strings
    std::string e347  = "347 keV";
    std::string e844  = "844 keV";
    std::string e1809 = "1809 keV";

    // Define table formatting
    const int nSF = 4;
    const int correctionNameColumnWidth = 20, signalColumnWidth = 30;
    const int fullWidth = correctionNameColumnWidth + signalColumnWidth * 3;

    // Print the title line and rules
    std::string title = "Normalization ratio (MC truth / Detector response)";
    std::cout << std::string(fullWidth, '=') << std::endl; // Title line
    std::cout << std::string((fullWidth - static_cast<int>(title.size())) / 2, ' ') << title << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // Title line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "" << std::setw(signalColumnWidth) << std::left << "347 keV" << std::setw(signalColumnWidth) << std::left << "844 keV" << std::setw(signalColumnWidth) << std::left << "1809 keV" << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "MC Truth";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringScientific(mcTruthNormalizations[j][0], nSF) + " ± " + doubleToStringScientific(mcTruthNormalizations[j][1], nSF);
    std::cout << std::endl;
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Detector Response";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringScientific(detectorResponseNormalizations[j][0], nSF) + " ± " + doubleToStringScientific(detectorResponseNormalizations[j][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Ratio";;
    for (int j = 0; j < 3; j++)
        std::cout << std::setw(signalColumnWidth) << std::left << doubleToStringFixed(ratioMCToDetector[j][0], nSF) + " ± " + doubleToStringFixed(ratioMCToDetector[j][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // End line
    std::cout << std::endl; // Buffer line
    return;
};