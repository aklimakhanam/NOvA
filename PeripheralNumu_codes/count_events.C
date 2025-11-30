void count_events(const std::string& filename = "Cosmics_MuE_CosRej.root")
{
  TFile* file = TFile::Open(filename.c_str(), "READ");

  TH2* h2 = dynamic_cast<TH2*>(file->Get("hMuonE_HadEFrac_et_fc_cos_0p58"));

  //h2->Scale(1./13842.23);
  h2->Scale(1./1103.91);
 
  std::cout << "\n Event Counts Above EHadFrac Thresholds in MuonE Bins " << std::endl;

  // Define your RecoE bins and associated EHadFrac thresholds
  std::vector<std::pair<float, float>> recoEBins = {
    {0.8, 1.0},
    {1.0, 1.2},
    {1.2, 1.4},
    {1.4, 1.6},
    {1.6, 1.8},
    {1.8, 2.0},
  };
  std::vector<float> hadFracThresholds = {0.24, 0.215, 0.195, 0.175, 0.17, 0.17};

  for (size_t i = 0; i < recoEBins.size(); ++i) {
    float recoELo = recoEBins[i].first;
    float recoEHi = recoEBins[i].second;
    float hadFracMin = hadFracThresholds[i];

    int binX_lo = h2->GetXaxis()->FindBin(recoELo + 1e-4);
    int binX_hi = h2->GetXaxis()->FindBin(recoEHi - 1e-4);
    int binY_lo = h2->GetYaxis()->FindBin(hadFracMin + 1e-4);
    int binY_hi = h2->GetYaxis()->GetNbins();

    double integral = h2->Integral(binX_lo, binX_hi, binY_lo, binY_hi);

    
  }

  file->Close();
}
