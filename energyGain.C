void energyGain(std::string detector)
{
  std::string filename = "stmSpectrum.root";
  TFile* data = new TFile(filename.c_str(), "READ");

  // Without subset
  TH1D* energySpectrum = (TH1D*) data->Get("plotSTMDigisSpectrum/energySpectrum");
  energySpectrum->Rebin(10);

  // Initializing the energy peak and adc peak locations
  const double energyPeaks[4] = {0.898,1.173,1.333,1.836};
  std::vector<double> energyPeaks;

  // Finding the location of each peak
  const int n_peaks = 7;
  TSpectrum* spectrum = new TSpectrum(n_peaks);
  int n_found_peaks = spectrum->Search(energySpectrum);
  for (int i_peak = 0; i_peak < n_found_peaks; ++i_peak)
    {
      double peak_x_pos = *(spectrum->GetPositionX() + i_peak);
      std::cout << peak_x_pos << std::endl;
      energyPeaks.push_back(peak_x_pos);
    }
  // Initializing the output file
  ofstream out;
  out.open("energyPeaks.log", ios::out | ios::trunc);
  out << "peak_number, par_name, par_value, par_error" << std::endl;

  // Plotting the adc spectrum
  energySpectrum->Draw("HIST");

  // Estimating background using TSpectrum
  TH1* hb = spectrum->Background(energySpectrum, 20, "SAME");

  TF1* fline = new TF1("fline", "pol1", 0, 5000);
  energySpectrum->Fit(fline, "qn");

  //energySpectrum->Add(hb,-1);
  std::cout << detector << std::endl;
  //Fitting
  for (int j = 0; j < energyPeaks.size(); ++j)
    {
      TF1* fitGaus;
      TString fname(Form("fgaus_%d", j));
      if (detector == "H")
        fitGaus = new TF1(fname, "[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4] + [5]*TMath::Gaus(x,[6],[7])",energyPeaks[j]-50,energyPeaks[j]+50);
      else if (detector == "L")
        fitGaus = new TF1(fname, "[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",energyPeaks[j]-250,energyPeaks[j]+250);
      fitGaus->SetParName(0,"Amplitude");
      fitGaus->SetParName(1,"Mean");
      fitGaus->SetParName(2,"Sigma");
      fitGaus->SetParName(3,"Linear");
      fitGaus->SetParName(4,"Constant");

      // Guessing the height for each TSpectrum peak
      int bin = energySpectrum->FindBin(energyPeaks[j]);
      double guessHeight = energySpectrum->GetBinContent(bin);
      std::cout << "TSpectrum value for peak " << j  <<  ": " << int (energyPeaks[j]) << std::endl;
      std::cout << "Guess for height of peak " << j  <<  ": " << guessHeight       << std::endl;

      // Guess for FWHM
      double bkrndHeight = hb->GetBinContent(bin);
      double sgnlHeight = guessHeight - bkrndHeight;
      double halfMax = sgnlHeight/2.0;
      double x1;
      double x2;
      bool threshold = false;
      std::cout << "Half Max: " << halfMax << std::endl;
      for (int k = bin - 25; k < bin + 25; ++k)
        {
          double k_content = energySpectrum->GetBinContent(k) - hb->GetBinContent(k);
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
      std::cout << "Sigma: " << sigma << std::endl;

      fitGaus->SetParameters(guessHeight,energyPeaks[j],sigma,-0.05,2);
      fitGaus->SetParLimits(1,energyPeaks[j]-10,energyPeaks[j]+10);
      //fitGaus->SetParLimits(0,0,200);
      //fitGaus->SetParLimits(2,0.001,10);

      // If the found peak is less than 500 ADC sample then don't fit
      if(energyPeaks[j] > 500) {
        TFitResultPtr fitresult = energySpectrum->Fit(fitGaus, "RS+");
        int n_par = fitresult->NPar();
        //Draw fitted peaks
        fitGaus->Draw("LSAME");
        // Export fit parameters to a file
        for (int i_par = 0; i_par < n_par; ++i_par)
          {
            out   << j << ", "
                  << fitresult->GetParameterName(i_par)
                  << ", " << fitresult->Parameter(i_par)
                  << ", " << fitresult->ParError(i_par)
                  << std::endl;
          }
      }
    }
  out.close();
}
