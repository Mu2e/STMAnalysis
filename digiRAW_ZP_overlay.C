void digiRAW_ZP_overlay() {

  std::string rawfile = "stmDigis.root";
  std::string ZPfile = "stmDigisZP.root";

  TFile* rawData = new TFile(rawfile.c_str(), "READ");
  TFile* ZPData = new TFile(ZPfile.c_str(), "READ");

  TH1D* rawWaveform = (TH1D*) rawData->Get("plotSTMDigis/digiSpectrum_1_0_Off");
  TH1D* ZPWaveform = (TH1D*) ZPData->Get("plotSTMDigis/digiSpectrum_1_0_Off");

  rawWaveform->Draw();
  ZPWaveform->Draw("SAME");
  ZPWaveform->SetLineColor(kRed);
}
