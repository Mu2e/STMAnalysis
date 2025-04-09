void CountMWDs(std::string filename) {
  TFile* file = new TFile(filename.c_str(), "READ");
  int count = 0;

  TDirectoryFile* data = (TDirectoryFile*) file->Get("plotSTMDigis");

  ofstream out;
  out.open("mwdContents.log", ios::out | ios::app );

  count = data->GetListOfKeys()->GetSize();

  std::cout << count << std::endl;
  out << filename << ", " << count << std::endl;
}
