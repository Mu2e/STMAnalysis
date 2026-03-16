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
void NormalizationComparisonToExpectedCount(std::vector<double> c, std::vector<double> cSim, std::vector<double> errSim) {
    // Define the parameters that contribute to the number of measured signal photons per muon capture, defined as {value, uncertainty} for each signal photon

    // Define the correction factors and their associated uncertianties
    //                                                          347 corr    uncert      844 corr    uncert      1809 corr   uncert
    std::vector<std::vector<double>> geantRateCorrection    = {{1.0,        0},        {2.590e-1,   0},         {1.0,       0}          };
    std::vector<std::vector<double>> resamplingCorrection   = {{1.0,        0},        {1.0,        0},         {0.01,      0}          };

    // Store all the correction factors in a single variable
    std::vector<std::vector<std::vector<double>>> correctionFactors = {geantRateCorrection, resamplingCorrection};
    std::vector<std::string> correctionFactorNames = {"GEANT correction", "Resampling factor correction"};
    const int nCorrectionFactors = correctionFactorNames.size();
    const int nSignals = geantRateCorrection.size();
    const std::vector<std::string> order = {"347 keV", "844 keV", "1809 keV"};
    int i = 0, j = 0;

    // Vectors to hold Data, Sim, and Ratio {value, absolute_uncertainty}
    std::vector<std::vector<double>> correctedSignalCounts = {{c[0], 0.0}, {c[1], 0.0}, {c[2], 0.0}};
    std::vector<std::vector<double>> correctedSimCounts    = {{cSim[0], 0.0}, {cSim[1], 0.0}, {cSim[2], 0.0}};
    std::vector<std::vector<double>> ratio                 = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    
// ====================================================================
    // Determine the counts & propagate errors (DATA ONLY CORRECTIONS)
    // ====================================================================
    for (int s = 0; s < nSignals; s++) {
        // 1. Initial fractional statistical uncertainty squared
        if (c[s] > 0)    correctedSignalCounts[s][1] += 1.0 / c[s];
        if (cSim[s] > 0) correctedSimCounts[s][1]    += std::pow(errSim[s] / cSim[s], 2);

        // 2. Loop through every correction factor for this specific signal
        for (int cf = 0; cf < nCorrectionFactors; cf++) {
            
            if (correctionFactors[cf][s][0] > 0) {
                // DIVIDE ONLY the DATA counts by the correction factor value
                // (This increases the Data count if the factor is < 1)
                correctedSignalCounts[s][0] /= correctionFactors[cf][s][0];
                
                // Add the fractional uncertainty squared to the DATA
                double frac_uncert_sq = std::pow(correctionFactors[cf][s][1] / correctionFactors[cf][s][0], 2);
                correctedSignalCounts[s][1] += frac_uncert_sq;
            } else {
                std::cerr << "Error: Correction factor for " << order[s] << " is zero! Cannot divide." << std::endl;
                exit(1); 
            }
        }
        
        // 3. Convert fractional uncertainty back to absolute uncertainty
        correctedSignalCounts[s][1] = correctedSignalCounts[s][0] * std::sqrt(correctedSignalCounts[s][1]);
        correctedSimCounts[s][1]    = correctedSimCounts[s][0] * std::sqrt(correctedSimCounts[s][1]);

        // 4. Compute the Ratio (Corrected Data / Raw Sim)
        if (correctedSimCounts[s][0] > 0) {
            ratio[s][0] = correctedSignalCounts[s][0] / correctedSimCounts[s][0]; 
            
            // sigma_R = R * sqrt( (sigma_DataCorr/DataCorr)^2 + (sigma_Sim/Sim)^2 )
            double data_frac_sq = (correctedSignalCounts[s][0] > 0) ? std::pow(correctedSignalCounts[s][1] / correctedSignalCounts[s][0], 2) : 0.0;
            double sim_frac_sq  = std::pow(correctedSimCounts[s][1] / correctedSimCounts[s][0], 2);
            
            ratio[s][1] = ratio[s][0] * std::sqrt(data_frac_sq + sim_frac_sq);
        }
    }
    // ====================================================================


    // Define the table formatting
    const int correctionNameColumnWidth = 45, signal347ColumnWidth = 30, signal844ColumnWidth = 30, signal1809ColumnWidth = 30;
    const int fullWidth = correctionNameColumnWidth + signal347ColumnWidth + signal844ColumnWidth + signal1809ColumnWidth;
    const int nSF = 4;

    // Print the title line and rules
    const std::vector<int> signalColumnWidths = {signal347ColumnWidth, signal844ColumnWidth, signal1809ColumnWidth};
    std::string tableTitle = "Data/Sim Comparison at VD90";

    std::cout << std::endl; // Buffer line
    std::cout << std::string(fullWidth, '=') << std::endl;
    std::cout << std::string((fullWidth - tableTitle.size())/2, ' ') << tableTitle << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Title line
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Parameter";
    for (i = 0; i < nSignals; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << order[i];
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the original signal counts (Data and Sim)
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Original count at VD90 (Sim)";
    for (i = 0; i < nSignals; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(c[i], nSF) + " ± " + doubleToStringScientific(std::sqrt(c[i]), nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the correction factors
    for (i = 0; i < nCorrectionFactors; i++) {
        std::cout << std::setw(correctionNameColumnWidth) << std::left << correctionFactorNames[i];
        for (j = 0; j < nSignals; j++)
            std::cout << std::setw(signalColumnWidths[j]) << std::left << doubleToStringScientific(correctionFactors[i][j][0], nSF) + " ± " + doubleToStringScientific(correctionFactors[i][j][1], nSF);
        std::cout << std::endl;
    };
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line

    // Print the corrected quantities
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Corrected count (Sim)";
    for (i = 0; i < nSignals; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(correctedSignalCounts[i][0], nSF) + " ± " + doubleToStringScientific(correctedSignalCounts[i][1], nSF);
    std::cout << std::endl;

    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Corrected count (Estimated)";
    for (i = 0; i < nSignals; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringScientific(correctedSimCounts[i][0], nSF) + " ± " + doubleToStringScientific(correctedSimCounts[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '-') << std::endl; // Section line
    
    // Print the Ratio
    std::cout << std::setw(correctionNameColumnWidth) << std::left << "Ratio (Sim / Estimated)";
    for (i = 0; i < nSignals; i++)
        std::cout << std::setw(signalColumnWidths[i]) << std::left << doubleToStringFixed(ratio[i][0], nSF) + " ± " + doubleToStringFixed(ratio[i][1], nSF);
    std::cout << std::endl;
    std::cout << std::string(fullWidth, '=') << std::endl; // Bottom rule

    return;
}