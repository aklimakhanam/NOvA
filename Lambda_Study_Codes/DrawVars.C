#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/HistAxis.h"
#include "TLegend.h"
#include "TMath.h"

#include "Utilities/rootlogon.C"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"

#include <iostream>
using namespace ana;

double POT = 12.5511e20;

TCanvas* PlotThreeHists(
    TH1* h1,
    TH1* h2,
    TH1* h3,
    const std::string& title,
    const std::string& label1,
    const std::string& label2,
    const std::string& label3,
    const std::string& outputdir,
    const std::string& filename
)
{
    // Create canvas
    TCanvas* c = new TCanvas("c", title.c_str(), 1200, 800);
    c->SetLogy();

    // Style histograms
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);

    h2->SetLineColor(kRed);
    h2->SetLineWidth(2);

    if (h3) {
        h3->SetLineColor(kGreen+2);
        h3->SetLineWidth(2);
    }

    // Title
    h1->SetTitle(title.c_str());

    // // Area normalize
    // h1->Scale(1.0 / h1->Integral());
    // h2->Scale(1.0 / h2->Integral());
    // if (h3) {
    //     h3->Scale(1.0 / h3->Integral());
    // }

    double maxY = h1->GetMaximum();
    if (h2) maxY = std::max(maxY, h2->GetMaximum());
    if (h3) maxY = std::max(maxY, h3->GetMaximum());
    maxY *= 1.2; // add 20% margin
    h1->SetMaximum(maxY);

    // Draw
    h1->Draw("hist");
    h2->Draw("hist same");

    if (h3) {
        h3->Draw("hist same");
    }

    // Legend
    TLegend* leg = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg->AddEntry(h1, label1.c_str(), "l");
    leg->AddEntry(h2, label2.c_str(), "l");
    if (h3) {
        leg->AddEntry(h3, label3.c_str(), "l");
    }
    leg->Draw();

    // Save
    c->SaveAs((outputdir + "/Plots_log/" + filename + ".png").c_str());

    return c;
}


TCanvas* Plot2DHist(
    TH2* h2,
    const std::string& title,
    const std::string& xLabel,
    const std::string& yLabel,
    const std::string& outputdir,
    const std::string& filename
)
{
    // Create canvas
    TCanvas* c = new TCanvas("c2d", title.c_str(), 1200, 800);

    // Apply title and axis labels
    h2->SetTitle(title.c_str());
    h2->GetXaxis()->SetTitle(xLabel.c_str());
    h2->GetYaxis()->SetTitle(yLabel.c_str());
    //log colz
    c->SetLogz();

    // Draw the 2D histogram with colz palette
    h2->Draw("colz");

    // Save the canvas
    c->SaveAs((outputdir + "/Plots_log/" + filename + ".png").c_str());

    return c;
}

void DrawVars(std::string infile = "lambda_study1.root", std::string outputdir = ".")
{
    // Open ROOT file
    TFile* fin = TFile::Open(infile.c_str(), "READ");
    if(!fin || fin->IsZombie()){
        std::cerr << "Cannot open file: " << infile << std::endl;
        return;
    }

    // Load spectra
    std::unique_ptr<Spectrum> sRecoENu_sig     = Spectrum::LoadFrom(fin, "sRecoENu_sig");
    std::unique_ptr<Spectrum> sRecoENu_presel = Spectrum::LoadFrom(fin, "sRecoENu_presel");
    std::unique_ptr<Spectrum> sRecoENu_failedsig = Spectrum::LoadFrom(fin, "sRecoENu_failedsig");
    std::unique_ptr<Spectrum> sProtonKE_sig    = Spectrum::LoadFrom(fin, "sProtonKE_sig");
    std::unique_ptr<Spectrum> sPiMinusKE_sig   = Spectrum::LoadFrom(fin, "sPiMinusKE_sig");
    std::unique_ptr<Spectrum> sNPngs_sig       = Spectrum::LoadFrom(fin, "sNPngs_sig");
    std::unique_ptr<Spectrum> sProtonCalE_sig    = Spectrum::LoadFrom(fin, "sProtonCalE_sig");
    std::unique_ptr<Spectrum> sPiMinusCalE_sig   = Spectrum::LoadFrom(fin, "sPiMinusCalE_sig");
    std::unique_ptr<Spectrum> sNHits_Png2_sig       = Spectrum::LoadFrom(fin, "sNHits_Png2_sig");
    std::unique_ptr<Spectrum> sNHits_Png3_sig       = Spectrum::LoadFrom(fin, "sNHits_Png3_sig");
    std::unique_ptr<Spectrum> sMuonID0_sig     = Spectrum::LoadFrom(fin, "sMuonID0_sig");
    std::unique_ptr<Spectrum> sProtonID0_sig     = Spectrum::LoadFrom(fin, "sProtonID0_sig");
    std::unique_ptr<Spectrum> sPionID0_sig     = Spectrum::LoadFrom(fin, "sPionID0_sig");
    std::unique_ptr<Spectrum> sMuonID1_sig     = Spectrum::LoadFrom(fin, "sMuonID1_sig");
    std::unique_ptr<Spectrum> sProtonID1_sig     = Spectrum::LoadFrom(fin, "sProtonID1_sig");
    std::unique_ptr<Spectrum> sPionID1_sig     = Spectrum::LoadFrom(fin, "sPionID1_sig");
    std::unique_ptr<Spectrum> sMuonID2_sig     = Spectrum::LoadFrom(fin, "sMuonID2_sig");
    std::unique_ptr<Spectrum> sProtonID2_sig     = Spectrum::LoadFrom(fin, "sProtonID2_sig");
    std::unique_ptr<Spectrum> sPionID2_sig     = Spectrum::LoadFrom(fin, "sPionID2_sig");
    std::unique_ptr<Spectrum> sMuonID0_presel     = Spectrum::LoadFrom(fin, "sMuonID0_presel");
    std::unique_ptr<Spectrum> sProtonID0_presel     = Spectrum::LoadFrom(fin, "sProtonID0_presel");
    std::unique_ptr<Spectrum> sPionID0_presel     = Spectrum::LoadFrom(fin, "sPionID0_presel");
    std::unique_ptr<Spectrum> sMuonID1_presel     = Spectrum::LoadFrom(fin, "sMuonID1_presel");
    std::unique_ptr<Spectrum> sProtonID1_presel     = Spectrum::LoadFrom(fin, "sProtonID1_presel");
    std::unique_ptr<Spectrum> sPionID1_presel     = Spectrum::LoadFrom(fin, "sPionID1_presel");
    std::unique_ptr<Spectrum> sMuonID2_presel     = Spectrum::LoadFrom(fin, "sMuonID2_presel");
    std::unique_ptr<Spectrum> sProtonID2_presel     = Spectrum::LoadFrom(fin, "sProtonID2_presel");
    std::unique_ptr<Spectrum> sPionID2_presel     = Spectrum::LoadFrom(fin, "sPionID2_presel");
    std::unique_ptr<Spectrum> sMuonID0_failedsig     = Spectrum::LoadFrom(fin, "sMuonID0_failedsig");
    std::unique_ptr<Spectrum> sProtonID0_failedsig     = Spectrum::LoadFrom(fin, "sProtonID0_failedsig");
    std::unique_ptr<Spectrum> sPionID0_failedsig     = Spectrum::LoadFrom(fin, "sPionID0_failedsig");
    std::unique_ptr<Spectrum> sMuonID1_failedsig     = Spectrum::LoadFrom(fin, "sMuonID1_failedsig");
    std::unique_ptr<Spectrum> sProtonID1_failedsig     = Spectrum::LoadFrom(fin, "sProtonID1_failedsig");
    std::unique_ptr<Spectrum> sPionID1_failedsig     = Spectrum::LoadFrom(fin, "sPionID1_failedsig");
    std::unique_ptr<Spectrum> sMuonID2_failedsig     = Spectrum::LoadFrom(fin, "sMuonID2_failedsig");
    std::unique_ptr<Spectrum> sProtonID2_failedsig     = Spectrum::LoadFrom(fin, "sProtonID2_failedsig");
    std::unique_ptr<Spectrum> sPionID2_failedsig     = Spectrum::LoadFrom(fin, "sPionID2_failedsig");
    std::unique_ptr<Spectrum> sProtonID1_PionID1_sig = Spectrum::LoadFrom(fin, "sProtonID1_PionID1_sig");
    std::unique_ptr<Spectrum> sProtonID2_PionID2_sig = Spectrum::LoadFrom(fin, "sProtonID2_PionID2_sig");
    std::unique_ptr<Spectrum> sProtonID1_PionID1_failedsig = Spectrum::LoadFrom(fin, "sProtonID1_PionID1_failedsig");
    std::unique_ptr<Spectrum> sProtonID2_PionID2_failedsig = Spectrum::LoadFrom(fin, "sProtonID2_PionID2_failedsig");

    // Convert to TH1
    TH1* hRecoENu_sig   = sRecoENu_sig->ToTH1(POT);
    TH1* hRecoENu_presel = sRecoENu_presel->ToTH1(POT);
    TH1* hRecoENu_failedsig = sRecoENu_failedsig->ToTH1(POT);
    TH1* hProtonKE_sig  = sProtonKE_sig->ToTH1(POT);
    TH1* hPiMinusKE_sig = sPiMinusKE_sig->ToTH1(POT);
    TH1* hNPngs_sig     = sNPngs_sig->ToTH1(POT);
    TH1* hProtonCalE_sig  = sProtonCalE_sig->ToTH1(POT);
    TH1* hPiMinusCalE_sig = sPiMinusCalE_sig->ToTH1(POT);
    TH1* hNHits_Png2_sig       = sNHits_Png2_sig->ToTH1(POT);
    TH1* hNHits_Png3_sig       = sNHits_Png3_sig->ToTH1(POT);
    TH1* hMuonID0_sig    = sMuonID0_sig->ToTH1(POT);
    TH1* hProtonID0_sig  = sProtonID0_sig->ToTH1(POT);
    TH1* hPionID0_sig    = sPionID0_sig->ToTH1(POT);
    TH1* hMuonID1_sig    = sMuonID1_sig->ToTH1(POT);
    TH1* hProtonID1_sig  = sProtonID1_sig->ToTH1(POT);
    TH1* hPionID1_sig    = sPionID1_sig->ToTH1(POT);
    TH1* hMuonID2_sig    = sMuonID2_sig->ToTH1(POT);
    TH1* hProtonID2_sig  = sProtonID2_sig->ToTH1(POT);
    TH1* hPionID2_sig    = sPionID2_sig->ToTH1(POT);
    TH1* hMuonID0_presel    = sMuonID0_presel->ToTH1(POT);
    TH1* hProtonID0_presel  = sProtonID0_presel->ToTH1(POT);
    TH1* hPionID0_presel    = sPionID0_presel->ToTH1(POT);
    TH1* hMuonID1_presel    = sMuonID1_presel->ToTH1(POT);
    TH1* hProtonID1_presel  = sProtonID1_presel->ToTH1(POT);
    TH1* hPionID1_presel    = sPionID1_presel->ToTH1(POT);
    TH1* hMuonID2_presel    = sMuonID2_presel->ToTH1(POT);
    TH1* hProtonID2_presel  = sProtonID2_presel->ToTH1(POT);
    TH1* hPionID2_presel    = sPionID2_presel->ToTH1(POT); 
    TH1* hMuonID0_failedsig    = sMuonID0_failedsig->ToTH1(POT);
    TH1* hProtonID0_failedsig  = sProtonID0_failedsig->ToTH1(POT);
    TH1* hPionID0_failedsig    = sPionID0_failedsig->ToTH1(POT);
    TH1* hMuonID1_failedsig    = sMuonID1_failedsig->ToTH1(POT);
    TH1* hProtonID1_failedsig  = sProtonID1_failedsig->ToTH1(POT);
    TH1* hPionID1_failedsig    = sPionID1_failedsig->ToTH1(POT);
    TH1* hMuonID2_failedsig    = sMuonID2_failedsig->ToTH1(POT);
    TH1* hProtonID2_failedsig  = sProtonID2_failedsig->ToTH1(POT);
    TH1* hPionID2_failedsig    = sPionID2_failedsig->ToTH1(POT); 
    // TH2* hProtonID1_PionID1_sig = sProtonID1_PionID1_sig->ToTH2(sProtonID1_PionID1_sig->POT());
    // TH2* hProtonID2_PionID2_sig = sProtonID2_PionID2_sig->ToTH2(sProtonID2_PionID2_sig->POT());
    TH2* hProtonID1_PionID1_sig = sProtonID1_PionID1_sig->ToTH2(POT);
    TH2* hProtonID2_PionID2_sig = sProtonID2_PionID2_sig->ToTH2(POT);
    TH2* hProtonID1_PionID1_failedsig = sProtonID1_PionID1_failedsig->ToTH2(POT);
    TH2* hProtonID2_PionID2_failedsig = sProtonID2_PionID2_failedsig->ToTH2(POT);

    // Set general style
    gStyle->SetOptStat(1110);

    TCanvas* c1 = new TCanvas("c1", "Lambda Spectra", 1200, 800);

    auto save_hist = [&](TH1* h, const std::string& title, const std::string& filename){
        h->SetTitle(title.c_str());
        h->SetLineColor(kBlue+1);
        h->SetLineWidth(2);
        h->Draw("hist");
        c1->SaveAs((outputdir + "/Plots_log/" + filename + ".png").c_str());
        // h->Reset(); // Clear canvas for next histogram
    };

    save_hist(hRecoENu_sig,   "Reconstructed Neutrino Energy (Signal);E_{#nu}^{Reco} [GeV];Events", "RecoENu");
    save_hist(hRecoENu_presel,   "Reconstructed Neutrino Energy (Signal + Background);E_{#nu}^{Reco} [GeV];Events", "RecoENu_presel");
    save_hist(hRecoENu_failedsig,   "Reconstructed Neutrino Energy (Background);E_{#nu}^{Reco} [GeV];Events", "RecoENu_failedsig");
    save_hist(hProtonKE_sig,  "Proton True Kinetic Energy;T_{p} [GeV];Events", "ProtonKE");
    save_hist(hPiMinusKE_sig, "PiMinus True Kinetic Energy;T_{#pi^{-}} [GeV];Events", "PiMinusKE");
    save_hist(hNPngs_sig,     "Number of Prongs;N_{prongs};Events", "NPngs");
    save_hist(hProtonCalE_sig,   " ;E_{p}^{Cal} [GeV];Events", "ProtonCalE");
    save_hist(hPiMinusCalE_sig,  " ;E_{#pi^{-}}^{Cal} [GeV];Events", "PiMinusCalE");
    save_hist(hNHits_Png2_sig,   "Number of Hits in 2nd Prongs;N_{hits};Events", "NHits_Png2");
    save_hist(hNHits_Png3_sig,   "Number of Hits in 3rd Prongs;N_{hits};Events", "NHits_Png3");
    save_hist(hMuonID0_sig,   "Muon ID of Leading Prong;MuonID;Events", "MuonID0");
    save_hist(hProtonID0_sig, "Proton ID of Leading Prong;ProtonID;Events", "ProtonID0");
    save_hist(hPionID0_sig,   "Pion ID of Leading Prong;PionID;Events", "PionID0");
    save_hist(hMuonID1_sig,   "Muon ID of 2nd Prong;MuonID;Events", "MuonID1");
    save_hist(hProtonID1_sig, "Proton ID of 2nd Prong;ProtonID;Events", "ProtonID1");
    save_hist(hPionID1_sig,   "Pion ID of 2nd Prong;PionID;Events", "PionID1");
    save_hist(hMuonID2_sig,   "Muon ID of 3rd Prong;MuonID;Events", "MuonID2");
    save_hist(hProtonID2_sig, "Proton ID of 3rd Prong;ProtonID;Events", "ProtonID2");
    save_hist(hPionID2_sig,   "Pion ID of 3rd Prong;PionID;Events", "PionID2");
    save_hist(hMuonID0_presel,   "Muon ID of Leading Prong (Base Selection);MuonID;Events", "MuonID0_Base");
    save_hist(hProtonID0_presel, "Proton ID of Leading Prong (Base Selection);ProtonID;Events", "ProtonID0_Base");
    save_hist(hPionID0_presel,   "Pion ID of Leading Prong (Base Selection);PionID;Events", "PionID0_Base");
    save_hist(hMuonID1_presel,   "Muon ID of 2nd Prong (Base Selection);MuonID;Events", "MuonID1_Base");
    save_hist(hProtonID1_presel, "Proton ID of 2nd Prong (Base Selection);ProtonID;Events", "ProtonID1_Base");
    save_hist(hPionID1_presel,   "Pion ID of 2nd Prong (Base Selection);PionID;Events", "PionID1_Base");
    save_hist(hMuonID2_presel,   "Muon ID of 3rd Prong (Base Selection);MuonID;Events", "MuonID2_Base");
    save_hist(hProtonID2_presel, "Proton ID of 3rd Prong (Base Selection);ProtonID;Events", "ProtonID2_Base");
    save_hist(hPionID2_presel,   "Pion ID of 3rd Prong (Base Selection);PionID;Events", "PionID2_Base");  
    save_hist(hMuonID0_failedsig,   "Muon ID of Leading Prong (Failed Signal Selection);MuonID;Events", "MuonID0_FailedSig");
    save_hist(hProtonID0_failedsig, "Proton ID of Leading Prong (Failed Signal Selection);ProtonID;Events", "ProtonID0_FailedSig");
    save_hist(hPionID0_failedsig,   "Pion ID of Leading Prong (Failed Signal Selection);PionID;Events", "PionID0_FailedSig");
    save_hist(hMuonID1_failedsig,   "Muon ID of 2nd Prong (Failed Signal Selection);MuonID;Events", "MuonID1_FailedSig");
    save_hist(hProtonID1_failedsig, "Proton ID of 2nd Prong (Failed Signal Selection);ProtonID;Events", "ProtonID1_FailedSig");
    save_hist(hPionID1_failedsig,   "Pion ID of 2nd Prong (Failed Signal Selection);PionID;Events", "PionID1_FailedSig");
    save_hist(hMuonID2_failedsig,   "Muon ID of 3rd Prong (Failed Signal Selection);MuonID;Events", "MuonID2_FailedSig");
    save_hist(hProtonID2_failedsig, "Proton ID of 3rd Prong (Failed Signal Selection);ProtonID;Events", "ProtonID2_FailedSig");
    save_hist(hPionID2_failedsig,   "Pion ID of 3rd Prong (Failed Signal Selection);PionID;Events", "PionID2_FailedSig");


    PlotThreeHists(
        hRecoENu_presel,
        hRecoENu_sig,
        hRecoENu_failedsig,
        "",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "RecoENu_Comparison"
    );   
    
    PlotThreeHists(
        hMuonID0_presel,
        hMuonID0_sig,
        hMuonID0_failedsig,
        "Muon ID of Leading Prong Comparison",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "MuonID0_Comparison"
    );
    PlotThreeHists(
        hProtonID0_presel,
        hProtonID0_sig,
        hProtonID0_failedsig,
        "Proton ID of Leading Prong Comparison",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "ProtonID0_Comparison"
    );
    PlotThreeHists(
        hPionID0_presel,
        hPionID0_sig,
        hPionID0_failedsig,
        "Pion ID of Leading Prong Comparison",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "PionID0_Comparison"
    );
    PlotThreeHists(
        hMuonID1_presel,
        hMuonID1_sig,
        hMuonID1_failedsig,
        "Muon ID of 2nd Prong Comparison",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "MuonID1_Comparison"
    );
    PlotThreeHists(
        hProtonID1_presel,
        hProtonID1_sig,
        hProtonID1_failedsig,
        "Proton ID of 2nd Prong Comparison (Log)",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "ProtonID1_Comparison"
    );
    PlotThreeHists(
        hPionID1_presel,
        hPionID1_sig,
        hPionID1_failedsig,
        "Pion ID of 2nd Prong Comparison (Log)",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "PionID1_Comparison"
    );
    PlotThreeHists(
        hMuonID2_presel,
        hMuonID2_sig,
        hMuonID2_failedsig,
        "Muon ID of 3rd Prong Comparison (Log)",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "MuonID2_Comparison"
    );
    PlotThreeHists(
        hProtonID2_presel,
        hProtonID2_sig,
        hProtonID2_failedsig,
        "Proton ID of 3rd Prong Comparison (Log)",
        "sig+bkg",
        "sig",
        "failed sig",
        outputdir,
        "ProtonID2_Comparison"
    );
    PlotThreeHists(
        hPionID2_presel,
        hPionID2_sig,
        hPionID2_failedsig,
        "Pion ID of 3rd Prong Comparison (Log)",
        "sig+bkg",
        "sig",
        "bkg",
        outputdir,
        "PionID2_Comparison"
    ); 

    PlotThreeHists(
        hProtonID1_sig,
        hPionID1_sig,
        nullptr,
        "Proton ID and Pion ID of 2nd prong",
        "Proton ID",
        "Pion ID",
        "",
        outputdir,
        "Proton_PionID1_Comparison"
    );
    
    PlotThreeHists(
        hProtonID2_sig,
        hPionID2_sig,
        nullptr,
        "Proton ID and Pion ID of 3rd prong",
        "Proton ID",
        "Pion ID",
        "",
        outputdir,
        "Proton_PionID2_Comparison"
    );

    Plot2DHist(
        hProtonID1_PionID1_sig,
        "Proton ID vs Pion ID of 2nd Prong (Sig)",
        "Proton ID",
        "Pion ID",
        outputdir,
        "ProtonID1_vs_PionID1"
    );
    Plot2DHist(
        hProtonID2_PionID2_sig,
        "Proton ID vs Pion ID of 3rd Prong (Sig)",
        "Proton ID",
        "Pion ID",
        outputdir,
        "ProtonID2_vs_PionID2"
    );

    Plot2DHist(
        hProtonID1_PionID1_failedsig,
        "Proton ID vs Pion ID of 2nd Prong (Bkg)",
        "Proton ID",
        "Pion ID",
        outputdir,
        "ProtonID1_vs_PionID1_failedsig"
    );
    Plot2DHist(
        hProtonID2_PionID2_failedsig,
        "Proton ID vs Pion ID of 3rd Prong (Bkg)",
        "Proton ID",
        "Pion ID",
        outputdir,
        "ProtonID2_vs_PionID2_failedsig"
    );

    fin->Close();

    std::cout << "All spectra plotted and saved in: " << outputdir << std::endl;
}