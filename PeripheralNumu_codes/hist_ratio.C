#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void hist_ratio() {
    // Open ROOT files
    TFile* file4 = TFile::Open("Cosmics_MuE_CosRej.root");
    TFile* file5 = TFile::Open("AntiNumu_MuE_CosRej.root");
    TFile* file6 = TFile::Open("AntiNumu_MuE_CosRej.root");
    TFile* file7 = TFile::Open("Cosmics_MuE_CosRej.root");

    // Get 2D histograms
    TH2F* hist4 = (TH2F*)file4->Get("hMuonE_HadEFrac_fc_cos_0p56");
    TH2F* hist5 = (TH2F*)file5->Get("hMuonE_HadEFrac_fc_0p56");
    TH2F* hist6_v1 = (TH2F*)file6->Get("hMuonE_HadEFrac_et_fc_0p56");
    TH2F* hist7_v1 = (TH2F*)file7->Get("hMuonE_HadEFrac_et_fc_cos_0p56");

    // Normalize histograms by POT/LiveTime
    hist4->Scale(1./1103.91);
    hist5->Scale(1./13842.23);
    hist6_v1->Scale(1./13842.23);
    hist7_v1->Scale(1./1103.91);

    // Create copies for ratio calculation
    TH2F* fc_cos = (TH2F*)hist4->Clone("fc_cos");
    TH2F* fc_numu = (TH2F*)hist5->Clone("fc_numu");
    TH2F* et_fc_numu = (TH2F*)hist6_v1->Clone("et_fc_numu");
    TH2F* et_fc_cos = (TH2F*)hist7_v1->Clone("et_fc_cos");

    fc_numu->Divide(fc_cos);
    et_fc_numu->Divide(et_fc_cos);

  /*
    hist4->GetYaxis()->SetRangeUser(0.3,1.0);  
    hist5->GetYaxis()->SetRangeUser(0.3,1.0);  
    hist6_v1->GetYaxis()->SetRangeUser(0.3,1.0); 
    hist7_v1->GetYaxis()->SetRangeUser(0.3,1.0); */

    // Project histograms onto x-axis
    TH1D* projX4 = hist4->ProjectionX("projX4");
    TH1D* projX5 = hist5->ProjectionX("projX5");
    TH1D* projX6_v1 = hist6_v1->ProjectionX("projX6_v1");
    TH1D* projX7_v1 = hist7_v1->ProjectionX("projX7_v1"); 

    // Set y-axis to start at zero
  /*
    double minY = 0;
    double maxY = std::max({projX4->GetMaximum(), projX5->GetMaximum()}) * 1.2; 

    // Create canvas

    TCanvas* c3 = new TCanvas("c3", "X Projections", 800, 600);
    c3->cd();
    projX4->GetYaxis()->SetRangeUser(1, maxY);
    projX4->SetMinimum(0);

    projX4->SetLineColor(kRed);
    projX4->SetLineWidth(2);
    projX4->Draw("HIST");

    projX5->SetLineColor(kBlue);
    projX5->SetLineWidth(2);
    projX5->Draw("HIST SAME");

    TLegend* legend1 = new TLegend(0.7, 0.6, 0.9, 0.9);
    legend1->SetTextSize(0.025);
    legend1->AddEntry(projX4, "Cosmics", "l");
    legend1->AddEntry(projX7_v1, "Cosmics_et_fc", "l");
    legend1->AddEntry(projX5, "AntiNumuCC", "l");
    legend1->AddEntry(projX6_v1, "AntiNumuCC_et_fc", "l");
    legend1->Draw();

    projX4->GetYaxis()->SetTitle("Events");  
    projX4->GetYaxis()->SetTitleOffset(1.6);  
    projX4->SetTitle("Failed containment, EHadFrac, PID > 0.48"); 
    //projX4->SetTitle("Failed containment variations"); 
   
    gStyle->SetStatX(0.25); 
    gStyle->SetStatY(0.9);  
    gStyle->SetStatW(0.15);  
    gStyle->SetStatH(0.12); 

    c3->SaveAs("antinumu_cosmics_MuE_fc_0p48.png");  */
    
    TCanvas *c4 = new TCanvas("c4", "S/B ratio", 1200, 600);
    c4->Divide(2,1);

    c4->cd(1);
    fc_numu->Draw("colz");
    fc_numu->SetTitle("Signal to Background ratio");
    fc_numu->GetYaxis()->SetTitle("E_{had} / E_{#nu}");

    c4->cd(2);
    et_fc_numu->Draw("colz");
    et_fc_numu->SetTitle("Signal to Background ratio");
    et_fc_numu->GetYaxis()->SetTitle("E_{had} / E_{#nu}");

    c4->SaveAs("S_B_ratio_0p58.png");

    // Save histograms to a ROOT file
    TFile* outFile = new TFile("S_B_ratio_histos_0p56.root", "RECREATE");
    fc_numu->Write("fc_numu_ratio");
    et_fc_numu->Write("et_fc_numu_ratio");
    outFile->Close();


  /*  TCanvas *c5 = new TCanvas("c5", "", 1200, 1200);
    c5->Divide(2, 2);

    c5->cd(1);
    hist4->Draw("colz");
    hist4->SetTitle("Uncontained cosmics w/o ET");

    c5->cd(3);
    hist5->Draw("colz");
    hist5->SetTitle("Uncontained AntiNumu w/o ET");

    c5->cd(4);
    hist6_v1->Draw("colz");
    hist6_v1->SetTitle("Uncontained AntiNumu w/ ET");

    c5->cd(2);
    hist7_v1->Draw("colz");
    hist7_v1->SetTitle("Uncontained cosmics w/ ET");

    // Save the canvas
    c5->SaveAs("spectra2D_MuE_0p4.png");  */

    TCanvas *c6 = new TCanvas("c6", "", 1200, 1200);
    c6->Divide(2, 1);

    c6->cd(1);
    hist4->Draw("colz");
    hist4->SetTitle("Uncontained cosmics");
    hist4->GetYaxis()->SetTitle("E_{had} / E_{#nu}"); 

    c6->cd(2);
    hist5->Draw("colz");
    hist5->SetTitle("Uncontained AntiNumu");
    hist5->GetYaxis()->SetTitle("E_{had} / E_{#nu}"); 

    // Save the canvas 
    c6->SaveAs("spectra2D_MuE.png");

    TFile* spec2D = new TFile("spectra2D_MuE.root", "RECREATE");
    hist4->Write("cosmics0p56");
    hist5->Write("antinumu0p56");
    spec2D->Close();


}
