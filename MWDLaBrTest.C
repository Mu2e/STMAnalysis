void MWDLaBrTest() {

  double tau = 35;//e3; // ns
  double T0 = 3.125; // clock tick, ns

  double min_t = 0;
  double max_t = 500;//e3;

  TF1* pulse = new TF1("pulse", "[0]*TMath::Exp(-x/[1])", min_t,max_t);
  TF1* pulse2 = new TF1("pulse2", "[0]*TMath::Exp(-x/[1])", min_t,max_t);

  int n_bins = (max_t - min_t) / T0;
  TH1F* hWaveform = new TH1F("hWaveform", "", n_bins,min_t,max_t);
  double pulse_start = 5;//e3; //ns
  double pulse2_start = 200;
  //  double pulse2_start = 50;//e3; //ns
  pulse->SetParameters(1000, tau);
  pulse2->SetParameter(1500, tau);
  double pedestal = 3800;
  for (int i_bin = 1; i_bin <= hWaveform->GetNbinsX(); ++i_bin) {
    double x = hWaveform->GetXaxis()->GetBinCenter(i_bin);
    double content = pedestal;

    if ((x > pulse_start) && (x <= pulse2_start)) { content += pulse->Eval(x); }
    //    if (x > pulse2_start) { content += pulse->Eval(x); }
    //    std::cout << "t = " << x << ", V = " << content << std::endl;
    if (x > pulse2_start) { content += pulse2->Eval(x); }
    hWaveform->SetBinContent(i_bin, content);
  }

  // Deconvolution
  TH1F* hDeconvolution = new TH1F("hDeconvolution", "", n_bins,min_t,max_t);
  hDeconvolution->SetBinContent(1, hWaveform->GetBinContent(1));
  for (int i_bin = 2; i_bin <= hDeconvolution->GetNbinsX(); ++i_bin) {
    double D = hWaveform->GetBinContent(i_bin)-pedestal - (1 - T0/tau)*(hWaveform->GetBinContent(i_bin-1)-pedestal) + hDeconvolution->GetBinContent(i_bin-1);
    hDeconvolution->SetBinContent(i_bin, D);
  }

  hDeconvolution->SetLineColor(kRed);

  hDeconvolution->Draw("HIST");
  hWaveform->Draw("HIST SAME");
}
