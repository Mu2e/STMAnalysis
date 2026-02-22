// Calculates the signal to total ratio from the MC truth information at VD90.
// Usage example - $ root -l -q 'SignalTotalFractionCalculation.C({22, 33, 116}, {44, 5, 109})'
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

std::string doubleToString(double value, int sigFigs) {
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

void SignalTotalFractionCalculation(std::vector<double> signal, std::vector<double> background, const double signalAcceptance = 0.1) {
    /*
        Description
            Generates the table describing the number of signal photons measured in the STM energy range

        Arguments
            signal - vector of signal photon counts for each energy bin
            background - vector of background photon counts for each energy bin

        Variables
    */

    // Sanity check the input vectors
    if (signal.size() != 3 || background.size() != 3) {
        std::cerr << "Error: Signal and background vectors must be of the same size." << std::endl;
        return;
    };
    for (size_t i = 0; i < signal.size(); ++i) {
        if (signal[i] < 0 || background[i] < 0) {
            std::cerr << "Error: Signal and background counts must be non-negative." << std::endl;
            return;
        };
    };


    // Construct the holding variables
    std::vector<std::vector<double>> signals = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}; // Vector of vectors to hold the signal counts and uncertainties for each energy
    std::vector<std::vector<double>> backgrounds = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}; // Vector of vectors to hold the background counts and uncertainties for each energy
    std::vector<std::vector<double>> totals = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}; // Vector of vectors to hold the total counts and uncertainties for each energy
    std::vector<std::vector<double>> signalToTotalRatios = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}; // Vector of vectors to hold the signal to total ratios and uncertainties for each energy

    for (size_t i = 0; i < signal.size(); ++i) {
        signals[i][0] = signal[i];
        backgrounds[i][0] = background[i];
        signals[i][1] = std::sqrt(signal[i]);
        backgrounds[i][1] = std::sqrt(background[i]);
        totals[i][0] = signal[i] + background[i];
        totals[i][1] = std::sqrt(signal[i] + background[i]);
        signalToTotalRatios[i][0] = signal[i]/totals[i][0];
        signalToTotalRatios[i][1] = signalToTotalRatios[i][0] * std::sqrt(std::pow(signals[i][1]/signals[i][0], 2) + std::pow(totals[i][1]/totals[i][0], 2));
    };

    // Define signal energy strings
    std::vector<std::string> energyStrings = {"347 keV", "844 keV", "1809 keV"};

    // Define the table title and column widths
    std::string tableTitle = "Signal to Total within " + doubleToString(signalAcceptance*100, 2) + "%";
    const int titleWidth = 40, signalWidth = 30, fullWidth = titleWidth + signalWidth * 3;

    // Print the table
    std::cout << std::endl; // Buffer line
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::string((fullWidth - tableTitle.size())/2, ' ') << tableTitle << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::setw(titleWidth) << std::left << " ";
    for (const std::string& energyString : energyStrings) {
        std::cout << std::setw(signalWidth) << std::left << energyString;
    };
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl;

    std::cout << std::setw(titleWidth) << std::left << "Signal photons";
    for (size_t i = 0; i < energyStrings.size(); ++i) {
        std::cout << std::setw(signalWidth) << std::left << doubleToString(signals[i][0], 1)+ " ± " + doubleToString(signals[i][1], 1);
    };
    std::cout << std::endl;

    std::cout << std::setw(titleWidth) << std::left << "Background photons";
    for (size_t i = 0; i < energyStrings.size(); ++i) {
        std::cout << std::setw(signalWidth) << std::left << doubleToString(backgrounds[i][0], 1)+ " ± " + doubleToString(backgrounds[i][1], 1);
    };
    std::cout << std::endl;

    std::cout<<std::string(fullWidth, '-')<<std::endl;

    std::cout << std::setw(titleWidth) << std::left << "Total photons";
    for (size_t i = 0; i < energyStrings.size(); ++i) {
        std::cout << std::setw(signalWidth) << std::left << doubleToString(totals[i][0], 1)+ " ± " + doubleToString(totals[i][1], 1);
    };
    std::cout << std::endl;

    std::cout<<std::string(fullWidth, '-')<<std::endl;

    std::cout << std::setw(titleWidth) << std::left << "Fraction";
    for (size_t i = 0; i < energyStrings.size(); ++i) {
        std::cout << std::setw(signalWidth) << std::left << doubleToString(signalToTotalRatios[i][0], 3)+ " ± " + doubleToString(signalToTotalRatios[i][1], 3);
    };
    std::cout << std::endl;

    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::endl; // Buffer line

    return;
};