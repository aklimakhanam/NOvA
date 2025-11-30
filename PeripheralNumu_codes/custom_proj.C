void custom_proj()
{
  // Open file and get histogram
  const std::string& filename_an = "AntiNumu_NumuE.root";
  const std::string& filename_cos = "Cosmics_NumuE.root";
  TFile* file_an = TFile::Open(filename_an.c_str(), "READ");
  TFile* file_cos = TFile::Open(filename_cos.c_str(), "READ");
  TH2* h1 = dynamic_cast<TH2*>(file_cos->Get("hMuonE_HadEFrac_fc_cos_0p56"));
  TH2* h2 = dynamic_cast<TH2*>(file_an->Get("hMuonE_HadEFrac_fc_0p56"));
  h1->Scale(1.0/1103.91);
  h2->Scale(1.0/13842.23);

  // RecoE bins and thresholds
  std::vector<std::pair<float,float>> recoEBins = {
    {0.8, 1.0},
    {1.0, 1.2},
    {1.2, 1.4},
    {1.4, 1.6},
    {1.6, 1.8},
    {1.8, 2.0},
  };
  std::vector<float> hadFracThresholds = {0.255, 0.230, 0.205, 0.185, 0.175, 0.165}; // w/o ET
  //std::vector<float> hadFracThresholds = {0.22, 0.2, 0.175, 0.16, 0.145, 0.135}; // w/ ET

  int nBins = recoEBins.size() + 2;
  // Make a histogram with nBins, bin edges matching recoE ranges
  double edges[9] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2};
  TH1D* hCounts_cos = new TH1D("hCounts_cos","Peripheral sample selection w/o ET ;Reco NumuE (GeV);Counts (POT normalized)", nBins, edges);
  TH1D* hCounts_an = new TH1D("hCounts_an","Peripheral sample selection w/o ET ;Reco NumuE (GeV);Counts (POT normalized)", nBins, edges);

  std::cout << "\n Event Counts (using Integral):\n";
  for (int i = 1; i < nBins - 1; ++i) {
    int binX_lo = h2->GetXaxis()->FindBin(recoEBins[i - 1].first + 1e-4);
    int binX_hi = h2->GetXaxis()->FindBin(recoEBins[i - 1].second - 1e-4);
    int binY_lo = h2->GetYaxis()->FindBin(hadFracThresholds[i - 1] + 1e-4);
    int binY_hi = h2->GetYaxis()->GetNbins();

    // Use ROOT Integral
    double count_cos = h1->Integral(binX_lo, binX_hi, binY_lo, binY_hi);
    double count_an = h2->Integral(binX_lo, binX_hi, binY_lo, binY_hi);

    hCounts_cos->SetBinContent(i+1, count_cos);
    hCounts_an->SetBinContent(i+1, count_an);

    std::cout << "Reco MuE [" << recoEBins[i - 1].first << ", " << recoEBins[i - 1].second << "], "
              << "EHadFrac > " << hadFracThresholds[i - 1]
              << " => Events (AntiNumu): " << count_an << std::endl;

    std::cout << "Reco MuE [" << recoEBins[i - 1].first << ", " << recoEBins[i - 1].second << "], "
              << "EHadFrac > " << hadFracThresholds[i - 1]
              << " => Events (Cosmics): " << count_cos << std::endl;
  }

  // Draw
  TCanvas* c = new TCanvas("c", "Event Counts", 800, 600);
  hCounts_cos->SetLineColor(kRed);
  hCounts_an->SetLineColor(kBlue);
  hCounts_cos->SetLineWidth(2);
  hCounts_an->SetLineWidth(2);
  hCounts_cos->SetStats(false);
  hCounts_an->SetStats(false);
  hCounts_an->Draw("HIST");
  hCounts_cos->Draw("HIST SAME");

  // Add a legend
  TLegend* legend = new TLegend(0.15, 0.7, 0.38, 0.88); // x1,y1,x2,y2 in NDC
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // transparent
  legend->AddEntry(hCounts_an, "#bar{#nu}_{#mu}", "l");
  legend->AddEntry(hCounts_cos, "Cosmics", "l");
  legend->Draw();

  c->SaveAs("event_counts_integral_numuE.png");

  //file->Close();
}
