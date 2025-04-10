// Plots the stages of the ZS analysis
// Usage example - $ root -l -q 'plotZSAnalysis.C("FinalData/CatZSAnalysis.root", "FinalData/CatZSWaveforms.root", "ZSHPGe/ttree", 31858, -100, 27000000, 39000000)'
// Original author - Pawel Plesniak

#include <TError.h>
#include <iostream>
#include <limits>
#include <math.h>

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectAnalysisData(const std::string fileName, const std::string treeName, std::vector<int16_t> &ADCs, std::vector<double> &times, std::vector<int16_t> &gradients, std::vector<double> &averagedGradients, std::vector<unsigned int> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotZSAnalysis"
            treeName - as documented in function "plotZSAnalysis"
            ADCs - as documented in function "plotZSAnalysis"
            times - as documented in function "plotZSAnalysis"
            gradients - as documented in function "plotZSAnalysis"
            averagedGradients - as documented in function "plotZSAnalysis"
            eventIds - as documented in function "plotZSAnalysis"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            dataTime - time from file
            ADCMin - minimum ADC value allowed
            dataGradient - gradient from file
            dataAveragedGradient - averaged gradient from file
            dataEventId - event ID from file
            entries - number of entries in the TTree
            i - iterator variable
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    const int16_t ADCMin = (-1 * std::pow(2, 15) + 1), overlapThreshold = 100;
    int16_t dataADC, dataGradient;
    uint32_t dataTime;
    double dataAveragedGradient;
    unsigned int dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "gradient") != branchNames.end())
        tree->SetBranchAddress("gradient", &dataGradient);
    else
        Fatal("collectData", "Requested branch 'gradient' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "averagedGradient") != branchNames.end())
        tree->SetBranchAddress("averagedGradient", &dataAveragedGradient);
    else
        Fatal("collectData", "Requested branch 'averagedGradient' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "eventId") != branchNames.end())
        tree->SetBranchAddress("eventId", &dataEventId);
    else
        Fatal("collectData", "Requested branch 'eventId' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Convert the time from ADC clock ticks to [ns]
    const double tADC = 3.125;

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        ADCs.push_back(dataADC);
        times.push_back(dataTime * tADC);
        gradients.push_back(dataGradient > overlapThreshold ? ADCMin : dataGradient);
        averagedGradients.push_back(dataAveragedGradient > overlapThreshold ? ADCMin : dataAveragedGradient);
        eventIds.push_back(dataEventId);
    };

    // Check if data has been collected
    if (ADCs.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void collectResultData(std::string fileName, std::string treeName, std::vector<int16_t> &ADCs, std::vector<double> &times, std::vector<unsigned int> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotZSAnalysis"
            treeName - as documented in function "plotZSAnalysis"
            ADCs - as documented in function "plotZSAnalysis"
            times - as documented in function "plotZSAnalysis"
            eventIds - as documented in function "plotZSAnalysis"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            dataTime - time from file
            dataEventId - event ID from file
            entries - number of entries in the TTree
            tADC - ADC clock tick length [ns]
            i - iterator variable
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    int16_t dataADC;
    uint32_t dataTime;
    unsigned int dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "eventId") != branchNames.end())
        tree->SetBranchAddress("eventId", &dataEventId);
    else
        Fatal("collectData", "Requested branch 'eventId' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Convert the time from ADC clock ticks to [ns]
    const double tADC = 3.125;

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        ADCs.push_back(dataADC);
        times.push_back(dataTime * tADC);
        eventIds.push_back(dataEventId);
    };

    // Check if data has been collected
    if (eventIds.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void plot(std::vector<int16_t> &inputADCs, std::vector<double> &inputTimes, std::vector<int16_t> &gradients, std::vector<double> &averagedGradients, std::vector<int16_t> &outputADCs, std::vector<double> outputTimes, unsigned int eventId, const double threshold, const double tMin, const double tMax, const std::string type) {
    /*
        Description
            Plots the ZS waveforms, their analysis steps, and the results

        Arguments
            inputADCs - as documented in function "plotZSAnalysis"
            inputTimes - as documented in function "plotZSAnalysis"
            gradients - as documented in function "plotZSAnalysis"
            averagedGradients - as documented in function "plotZSAnalysis"
            outputADCs - as documented in function "plotZSAnalysis"
            outputTimes - as documented in function "plotZSAnalysis"
            eventId - as documented in function "plotZSAnalysis"
            tMin - as documented in function "plotZSAnalysis"
            tMax - as documented in function "plotZSAnalysis"
            type - as documented in function "plotZSAnalysis"

        Variables
            nInput - number of input ADC values
            nOutput - number of output ADC values
            time - data product used as x axis for plotting the full ADC range
            tADC - time interval of the ADC clock tick [ns]
            tMaxFull - time at the end of the input ADC. Note time is modelled from 0 for the input ADCs
            tMinCutIndex - index of the minimum time
            tMaxCutIndex - index of the maximum time
            tMaxCut - maximum time plotted
            time - vector of ADC times for input data
            inputADCsA - array of input ADCs
            outputADCsA - array of output ADCs
            gradientsA - array of analysis gradients
            averagedGradientsA - array of analysis averaged gradients
            outputTimesA - array of result output itmes
            nInput - number of accepted input data entries
            nOutput - numnber of accepted output data entries
            c - canvas used to generate the plot
            gInputADCs - graph of input ADCs
            gGradients - graph of analysis gradients
            gAveragedGradient - graph of analysis averaged gradients
            splitIndexes - indexes of discontinuities in the output ADC times
            splitOutputTime - vector of times used to generate the output waveforms
            splitOutputADCs - vector of ADCs used to generate the output waveofrms
            nOutputPlots - number of output waveform plots to generate. Need to be separate as otherwise the TGraph connects the discontinuities displaying data that is not there
            outputPlots - plots of output waveforms
            x1, x2, y1, y1 - coordinates of legend as a percentage of the canvas area
            legend - legend for the plot
            cFileName - filename the plot is saved as
    */

    // Count the number of entries used to generate this plot
    int nInput = inputADCs.size(), nOutput = outputADCs.size();

    // Construct the time data used for the ZS algoritm analysis stages
    const double tADC = 3.125, tEvent0 = 1696.88; // tEvent0 is needed as the first event is event 1 and not 0 and the start time of this event is tEvent0 ns
    const double tMinFull = std::min(*std::min_element(inputTimes.begin(), inputTimes.end()), *std::min_element(outputTimes.begin(), outputTimes.end())) - tEvent0;
    const double tMaxFull = std::max(*std::max_element(inputTimes.begin(), inputTimes.end()), *std::max_element(outputTimes.begin(), outputTimes.end())) - tEvent0;

    // Sanity checks
    if (tMin > tMaxFull)
        Fatal("plot", "Selected minimum cut is greater than the maximum, exiting.");

    // Convert the times to the number of clock ticks for selecting the range of entries to plot
    const int tMinCutIndex = (tMinFull > tMin || tMin < std::numeric_limits<double>::epsilon()) ? static_cast<int>(tMinFull / tADC) : static_cast<int>(tMin / tADC);
    const int tMaxCutIndex = (tMaxFull < tMax || tMax < std::numeric_limits<double>::epsilon()) ? static_cast<int>(tMaxFull / tADC) : static_cast<int>(tMax / tADC);
    const double tMinns = tMinCutIndex * tADC, tMaxns = tMaxCutIndex * tADC;

    // Sanity checks
    if (tMaxCutIndex - tMinCutIndex < 10)
        Fatal("plot", "Selected data range is such that there are very few entries to plot, exiting.");

    // Construct the variables used to store the data relevant for plotting
    std::vector<double> inputADCsA, outputADCsA, gradientsA, averagedGradientsA, inputTimesA, outputTimesA;

    // Select and convert the relevant data for plotting
    for (int i = tMinCutIndex; i < tMaxCutIndex; i++) {
        inputADCsA.push_back(static_cast<double>(inputADCs[i]));
        inputTimesA.push_back(inputTimes[i]);
        gradientsA.push_back(static_cast<double>(gradients[i]));
        averagedGradientsA.push_back(averagedGradients[i]);
    };
    nInput = inputADCsA.size();
    double t = 0.0;
    for (int i = 0; i < nOutput; i++) {
        t = outputTimes[i];
        if (tMinns < t && t < tMaxns) {
            outputADCsA.push_back(static_cast<double>(outputADCs[i]));
            outputTimesA.push_back(t);
        };
    };
    nOutput = outputADCsA.size();
    if (nOutput == 0)
        Fatal("plot", "Number of output points is zero, exiting");

    // Set up the canvas for the plot
    TCanvas* c = new TCanvas("c", "c", 1500, 1000);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.08);

    // Construct the legend
    const double x1 = 0.7, x2 = 0.9, y1 = 0.3, y2 = 0.6;
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->SetTextSize(0.04);

    if (type == "fit") {
        // Plot the gradients
        TGraph *gGradients = new TGraph(nInput, inputTimesA.data(), gradientsA.data());
        gGradients->GetXaxis()->SetLimits(tMinns, tMaxns);
        gGradients->SetLineColor(kRed);
        gGradients->SetLineWidth(2);
        gGradients->Draw("APL");
        gGradients->SetTitle(Form("Zero suppression %s;Time [ns];ADC [arb. unit]", type.c_str()));
        legend->AddEntry(gGradients, "Grad.", "l")->SetTextColor(kRed);

        // Plot the analysis gradients
        TGraph *gAveragedGradients = new TGraph(nInput, inputTimesA.data(), averagedGradientsA.data());
        gAveragedGradients->SetLineColor(kBlue);
        gAveragedGradients->SetLineWidth(2);
        gAveragedGradients->Draw("PL SAME");
        legend->AddEntry(gAveragedGradients, "Avg. Grad.", "l")->SetTextColor(kBlue);

        // Add the averaged gradient threshold
        TLine *averagedGradientThreshold = new TLine(tMinns, threshold, tMaxns, threshold);
        averagedGradientThreshold->SetLineColor(kGreen + 2);
        averagedGradientThreshold->SetLineWidth(2);
        averagedGradientThreshold->Draw("SAME");
        legend->AddEntry(averagedGradientThreshold, "Threshold", "l")->SetTextColor(kGreen + 2);
    }
    else {
        // Plot the zero-suppressed waveforms
        // Plot the input ADCs
        TGraph *gInputADCs = new TGraph(nInput, inputTimesA.data(), inputADCsA.data());
        gInputADCs->GetXaxis()->SetLimits(tMinns, tMaxns);
        gInputADCs->SetLineWidth(3);
        gInputADCs->Draw("APL");
        gInputADCs->SetTitle(Form("Zero suppression %s;Time [ns];ADC [arb. unit]", type.c_str()));
        legend->AddEntry(gInputADCs, "Input ADCs", "l");

        // Determine the start indexes of separate waveforms
        std::vector<uint> splitIndexes = {0};
        for (uint i = 1; i < nOutput; i++) {
            if ((outputTimesA[i] - outputTimesA[i - 1]) > tADC)
                splitIndexes.push_back(i);
        };
        splitIndexes.push_back(nOutput);
        const int nOutputPlots = splitIndexes.size();

        // Construct the plots of zero-suppressed waveforms
        std::vector<TGraph*> outputPlots;
        if (nOutputPlots > 2) { // Case where we have multiple peaks - needs to be split due to the ADC discontinuity
            // Construct variables to hold the plot data
            std::vector<double> splitOutputTimes, splitOutputADCs;

            // Generate the plots
            for (uint i = 1; i < nOutputPlots; i++) {
                splitOutputTimes.assign(outputTimesA.begin() + splitIndexes[i - 1], outputTimesA.begin() + splitIndexes[i]);
                splitOutputADCs.assign(outputADCsA.begin() + splitIndexes[i - 1], outputADCsA.begin() + splitIndexes[i]);
                outputPlots.emplace_back(new TGraph(splitIndexes[i] - splitIndexes[i - 1], splitOutputTimes.data(), splitOutputADCs.data()));
                outputPlots.back()->GetXaxis()->SetLimits(tMinns, tMaxns);
                outputPlots.back()->Draw("PL SAME");
                outputPlots.back()->SetLineColor(kGreen + 2);
                outputPlots.back()->SetLineWidth(10);
                splitOutputTimes.clear();
                splitOutputADCs.clear();
            };
        }
        else {
            outputPlots.emplace_back(new TGraph(nOutput, outputTimesA.data(), outputADCsA.data()));
            outputPlots.back()->GetXaxis()->SetLimits(tMinns, tMaxns);
            outputPlots.back()->Draw("PL SAME");
            outputPlots.back()->SetLineColor(kGreen + 2);
            outputPlots.back()->SetLineWidth(10);
        };

        // Add the relevant legend entry
        legend->AddEntry(outputPlots.back(), "Output ADCs", "l")->SetTextColor(kGreen + 2);
    };

    // Draw the legend
    legend->Draw();
    gPad->RedrawAxis();
    gPad->Update();
    c->Modified();
    c->Update();

    // Save the generated plot
    std::string cFileName = "ZS.Event" + std::to_string(eventId) + "." + type + ".png";
    c->SaveAs(cFileName.c_str());

    // Cleanup
    legend->Delete();
    c->Close();

    return;
};

std::vector<unsigned int> makeUniqueEventIds(std::vector<unsigned int> &inputEventIds, std::vector<unsigned int> &outputEventIds) {
    /*
        Description
            Generates a unique intersection vector of the event IDs retrieved from the input files

        Argument
            inputEventIds - as documented in function "plotZSAnalysis"
            outputEventIds - as documented in function "plotZSAnalysis"

        Variables
            inputEventIdSet - set of event IDs from ZS input data
            uniqueInputEventIds - vector of unique event IDs from ZS input data
            outputEventIdSet - set of event IDs from ZS output data
            uniqueOutputEventIds - vector of unique event IDs from ZS output data
            overlap - vector of unique event IDs present in both input and output data
    */
    std::cout << "\nMaking the event IDs unique" << std::endl;

    // Get a unqiue set of input event IDs
    std::set<unsigned int> inputEventIdSet(inputEventIds.begin(), inputEventIds.end());

    // Convert the set to a vector
    std::vector<unsigned int> uniqueInputEventIds(inputEventIdSet.begin(), inputEventIdSet.end());

    // Get a unique set of output event IDs
    std::set<unsigned int> outputEventIdSet(outputEventIds.begin(), outputEventIds.end());

    // Convert the set to a vector
    std::vector<unsigned int> uniqueOutputEventIds(outputEventIdSet.begin(), outputEventIdSet.end());

    // Find the intersection of the two unique event ID vectors
    std::vector<unsigned int> overlap;
    std::set_intersection(uniqueInputEventIds.begin(), uniqueInputEventIds.end(), uniqueOutputEventIds.begin(), uniqueOutputEventIds.end(), std::back_inserter(overlap));

    // Return the unqiue intersection
    return overlap;
}

void plotZSAnalysis(const std::string analysisFileName, const std::string resultFileName, const std::string treeName, const unsigned int eventID = 0, const double threshold = -100, const double tMin = 0.0, const double tMax = 0.0) {
    /*
        Description
            Generates a plot of the ZS waveform analysis with file name
                ZS.Event<EventID>.<type>.png
            <EventID> is allocated even if the parameter "eventID" is not used.
            <type> is either "fit" or "results"

        Arguments
            fileName - name of the root file generated with STMZeroSuppression_module.cc
            treeName - name of the ttree containing all the data in fileName
            eventID - controls what eventID to plot for. If 0, generates plots for all events
            threshold - ADC threshold for averaged gradient data
            tMin - minimum ttime to plot [ns]. If zero, does not apply a time cut [ns]
            tMax - maximum time to plot [ns]. If zero, does not apply a time cut [ns]

        Variables
            inputADCs - vector of ADC values used as input to ZS algorithm
            outputADCs - vector of ADC values stored by ZS algorithm
            gradients - vector of gradients calculated with ZS algorithm
            averagedGradients - vector of averaged gradients calculated with ZS algorithm
            inputTimes - vector of times associated with inputADCs
            outputTimes - vector of times associated with outputADCs
            inputEventIds - vector of eventIds used as input to ZS algorithm
            outputEventIds - vector of eventIds saved by ZS algorithm
            nInput - number of ADC entries imported from the ZS algorithm input
            nOutput - number of ADC entried imported from the ZS algorithm output
            plotInputADCs - selected input ADCs to use when generating the input plot
            plotOutputADCs - selected output ADCs to use when generating the input plot
            plotGradients - selected gradients to use when generating the input plot
            plotAveragedGradients - selected averaged gradients to use when generating the input plot
            plotOutputTimes - selected output times to use when generating the input plot
            plotEventIds - selected event IDss to use when generating the input plot
            plotEventId - iterator variable for plotEventIds
            plotType - vector of generated plot types
            type - iterator for plotType
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Sanity checks
    if (tMin < 0)
        Fatal("plotZSAnalysis", "tMin must be positive!");
    if (tMax < 0)
        Fatal("plotZSAnalysis", "tMax must be positive!");
    if (tMax < tMin)
        Fatal("plotZSAnalysis", "tMax must be greater than tMin!");

    // Construct the variables used to hold the input data
    std::vector<int16_t> inputADCs, outputADCs, gradients;
    std::vector<double> averagedGradients, inputTimes, outputTimes;
    std::vector<unsigned int> inputEventIds, outputEventIds;

    // Read the data from the input files
    collectAnalysisData(analysisFileName, treeName, inputADCs, inputTimes, gradients, averagedGradients, inputEventIds);
    collectResultData(resultFileName, treeName, outputADCs, outputTimes, outputEventIds);
    int nInput = inputADCs.size(), nOutput = outputADCs.size();

    // Construct the variables used to hold the data used for plotting
    std::vector<int16_t> plotInputADCs, plotOutputADCs, plotGradients;
    std::vector<double> plotAveragedGradients, plotInputTimes, plotOutputTimes;
    std::vector<unsigned int> plotEventIds = makeUniqueEventIds(inputEventIds, outputEventIds);

    // Validate the event ID in the case that it is not present in the input files
    if (eventID != 0) {
        if (std::find(plotEventIds.begin(), plotEventIds.end(), eventID) == plotEventIds.end())
            Fatal("plotZSAnalysis", "Requested eventID is not one of those from the file, exiting.");
        plotEventIds.clear();
        plotEventIds.push_back(eventID);
    };

    // Generate the plots
    std::vector<std::string> plotType = {"fit", "results"};
    for (unsigned int plotEventId : plotEventIds) {
        std::cout << "\nGenerating plot for event with ID " << plotEventId << std::endl;

        // Select the relevant data for plotting
        for (int i = 0; i < nInput; i++) {
            if (inputEventIds[i] == plotEventId) {
                plotInputADCs.push_back(inputADCs[i]);
                plotInputTimes.push_back(inputTimes[i]);
                plotGradients.push_back(gradients[i]);
                plotAveragedGradients.push_back(averagedGradients[i]);
            };
        };
        for (int i = 0; i < nOutput; i++) {
            if (outputEventIds[i] == plotEventId) {
                plotOutputADCs.push_back(outputADCs[i]);
                plotOutputTimes.push_back(outputTimes[i]);
            };
        }

        // Sanity checks
        if (plotInputADCs.empty())
            Fatal("plotZSAnalysis", "Empty input ADC vector, exiting!");
        if (plotInputTimes.empty())
            Fatal("plotZSAnalysis", "Empty input time vector, exiting!");
        if (plotGradients.empty())
            Fatal("plotZSAnalysis", "Empty gradient vector, exiting!");
        if (plotAveragedGradients.empty())
            Fatal("plotZSAnalysis", "Empty averaged gradient vector, exiting!");
        if (plotOutputADCs.empty())
            Fatal("plotZSAnalysis", "Empty output ADC vector, exiting!");
        if (plotOutputTimes.empty())
            Fatal("plotZSAnalysis", "Empty output time vector, exiting!");

        // Generate the plots
        for (std::string type : plotType)
            plot(plotInputADCs, plotInputTimes, plotGradients, plotAveragedGradients, plotOutputADCs, plotOutputTimes, plotEventId, threshold, tMin, tMax, type);

        // Prepare the buffer variables for the next event ID
        plotInputADCs.clear();
        plotInputTimes.clear();
        plotGradients.clear();
        plotAveragedGradients.clear();
        plotOutputADCs.clear();
        plotOutputTimes.clear();
    };
    return;
};
