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

void NormalizationSourceComparison(std::vector<std::vector<double>> VD90Signals, std::vector<std::vector<double>> DetectorMCSignals, std::vector<std::vector<double>> DetectorRecoSignals, const double nPOTs = 2.69e16) {
    /*
        Description
            Generates the comparison of the normalization estimates from the different sources

        Arguments
            VD90Signals - signal photons measured in VD90
            DetectorMCSignals - signal photons measured in the detector Monte Carlo
            DetectorRecoSignals - signal photons measured in the detector reconstruction
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
    std::vector<std::string> normalization_types{"all", "VD90", "DetectorMC", "DetectorReco"};
    if (std::find(normalization_types.begin(), normalization_types.end(), normalization_type) == normalization_types.end()) {
        std::cout << "Provided normalization type: " << normalization_type << std::endl;
        std::cerr << "Error: Normalization type must be one of 'all', 'VD90', 'DetectorMC', or 'DetectorReco'." << std::endl;
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
    const double uMDC2020_TargetStopsCat = 1000.0 * std::sqrt(1432353.0); // Propagate tthe error correctly
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
    std::vector<std::vector<std::vector<double>>> correctionFactors = {pFinalState, geantRateCorrection, geometricAcceptance, pathAttenuation, absorberAcceptance, energyWindowAcceptance, timeCutAcceptance, detectorAcceptance, clippingFactor};


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
        measuredToExpectedNormalizationRatio[i][0] = (nSignalPhotons[i][0] / geantRateCorrection[i][0]) / expectedPhotonCountFromPOT[i][0];
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

    return;
};
