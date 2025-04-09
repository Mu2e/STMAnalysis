void stmSpectrumRateOverlay()
{
  std::string low_rate = "stmSpectrumLowRate.root";
  std::string high_rate = "stmSpectrumHighRate.root";

  TFile* low = new TFile(low_rate.c_str(), "READ");
  TFile* high = new TFile(high_rate.c_str(), "READ");

  TH1D* low_spectrum = (TH1D*) low->Get("plotSTMSpectrum/energySpectrum");
  TH1D* high_spectrum = (TH1D*) high->Get("plotSTMSpectrum/energySpectrum");

  low_spectrum->Rebin(10);
  high_spectrum->Rebin(10);

  low_spectrum->Draw();
  high_spectrum->Draw("SAME");
  high_spectrum->SetLineColor(kRed);

}
