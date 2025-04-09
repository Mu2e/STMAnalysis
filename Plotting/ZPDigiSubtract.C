void ZPDigiSubtract() {
  std::string rawfile = "stmDigis.root";
  std::string ZPfile = "stmDigisZP.root";

  TFile* rawData = new TFile(rawfile.c_str(), "READ");
  TFile* ZPData = new TFile(ZPfile.c_str(), "READ");

  TH1D* rawWaveform = (TH1D*) rawData->Get("plotSTMDigis/digiSpectrum_1_0_Off");

  for(int i = 0; i < 189; ++i)
    {
      TH1D* ZPWaveform = (TH1D*) ZPData->Get(("plotSTMDigis/digiSpectrum_1_" + to_string(i) + "_Off").c_str());
      rawWaveform->Add(ZPWaveform,-1);
    }
  rawWaveform->Draw();
}
