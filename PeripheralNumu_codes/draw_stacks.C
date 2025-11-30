void draw_stacks() {
  // Open file
  TFile *f = TFile::Open("dumb_pred_test_output.root");
  TH1F *h_beam  = (TH1F*)f->Get("h_beam");
  TH1F *h_unosc = (TH1F*)f->Get("h_unosc");


  // Load cosmics
  TFile* fCosm = new TFile("spectra1D_nowgt_rhc_high_stat.root","READ");
  TH1* h_cosm = (TH1*)fCosm->Get("h_NuMuE");

  TH1* h_cosm_unosc = (TH1*)h_cosm->Clone("h_cosm_unosc");

  // Style
  // Remove fill
  h_cosm->SetFillStyle(0);
  h_cosm_unosc->SetFillStyle(0);
  h_beam->SetFillStyle(0);
  h_unosc->SetFillStyle(0);

  // Line colors
  h_cosm->SetLineColor(kGreen+2);
  h_cosm_unosc->SetLineColor(kGreen+2);
  h_beam->SetLineColor(kRed);
  h_unosc->SetLineColor(kBlack);

  // Line widths
  h_cosm->SetLineWidth(2); 
  h_cosm_unosc->SetLineWidth(2);
  h_beam->SetLineWidth(2);
  h_unosc->SetLineWidth(2);

  // Make stacks
  THStack* st_osc = new THStack("st_osc", "");
  st_osc->Add(h_cosm);
  st_osc->Add(h_beam);

  THStack* st_unosc = new THStack("st_unosc", "");
  st_unosc->Add(h_cosm_unosc);
  st_unosc->Add(h_unosc);

  // Draw
  TCanvas* c = new TCanvas("c", "Stacked comparison", 900, 700);
  st_unosc->Draw("HIST");
  st_unosc->GetXaxis()->SetTitle("Captured Reconstructed Neutrino Energy (GeV)");
  st_unosc->GetYaxis()->SetTitle("Events");

  st_osc->Draw("HIST SAME");

  // Legend
  TLegend* leg = new TLegend(0.62, 0.65, 0.88, 0.88);
  leg->AddEntry(h_cosm,  "Cosmics",            "l");
  leg->AddEntry(h_beam,  "Oscillated Beam",    "l");
  leg->AddEntry(h_unosc, "Unoscillated Beam",  "l");
  leg->Draw();

  c->SaveAs("stack_overlay.png");
  c->SaveAs("stack_overlay.pdf");

  std::cout << "Done." << std::endl;
}
