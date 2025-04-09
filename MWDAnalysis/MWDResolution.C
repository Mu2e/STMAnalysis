// MWDResolution
// Extracts the resolution for the Cs peak from the ENERGY spectrum for each (M,L) in the LaBr detector

void MWDResolution(std::string filename)
{
  TFile* file = new TFile(filename.c_str(), "READ");
  TH1D* energySpectrum = (TH1D*) file->Get("plotSTMSpectrum/energySpectrum"); // "plotSTMDigisSpectrum/adcSpectrum"
  double resolution, res_error;
  double peak_mean_guess = 1.836; //0.662; // MeV, 2230.0 ADC;
  std::vector<double> adcPeaks;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TF1* fitGaus = new TF1("fitGaus", "[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",peak_mean_guess-0.15,peak_mean_guess+0.15);
  // pean_mean_guess +- 250 for ADC units
  energySpectrum->Rebin(2);
  //energySpectrum->Draw("HIST");

  //-------------------------------------------------------------------------------------------------------------------
  // Guessing the height for each TSpectrum peak
  int bin = energySpectrum->FindBin(peak_mean_guess);
  double guessHeight = energySpectrum->GetBinContent(bin);
  std::cout << "Guess for height of peak: " << guessHeight << std::endl;

  // Guess for FWHM
  double sgnlHeight = guessHeight;
  double halfMax = sgnlHeight/2.0;
  double x1, x2;
  bool threshold = false;
  //std::cout << "Half Max: " << halfMax << std::endl;
  for (int k = bin - 25; k < bin + 25; ++k) // Use +- 25 for LaBr and +-5 for HPGe
    {
      double k_content = energySpectrum->GetBinContent(k);
      //std::cout << "k_content: " << k_content << std::endl;
      if (k_content > halfMax && threshold == false)
        {
          x1 = energySpectrum->GetBinCenter(k);
          threshold = true;
          //std::cout << "x1: " << x1 << std::endl;
        }
     if (k_content < halfMax && threshold == true)
        {
          x2 = energySpectrum->GetBinCenter(k);
          //std::cout << "x2: " << x2 << std::endl;
          break;
        }
    }

  double FWHM = x2 - x1;
  //std::cout << "FWHM: " << FWHM << std::endl;
  double sigma = FWHM/2.35;
  std::cout << "Estimated sigma: " << sigma << std::endl;

  //------------------------------------------------------------------------------------------------------------------
  fitGaus->SetParameters(guessHeight,peak_mean_guess,sigma,-0.05,2);
  fitGaus->SetParLimits(0,0,3e2);
  fitGaus->SetParLimits(1,0,2.5);
  fitGaus->SetParLimits(2,0,50);

  TFitResultPtr fitresult = energySpectrum->Fit(fitGaus, "RS+");
  std::cout << fitGaus->Integral(1.1,1.2) << std::endl;
  resolution = fitresult->Parameter(2);
  res_error = fitresult->ParError(2);
  //fitGaus->Draw("LSAME");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Append to the list of (M,L) resolutions
  ofstream out;
  out.open("mwdResolution.log", ios::out | ios::app );
  out << filename << ", " << resolution << std::endl;
}
