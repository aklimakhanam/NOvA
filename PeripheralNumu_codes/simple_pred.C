// Use that prediction to make a contour for a really dumb numu analysis

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Fit/Fit.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Fit/FrequentistSurface.h"
#include "CAFAna/Analysis/Exposures.h"

#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/Vars/Binnings.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"
#include "3FlavorAna/Vars/NueVars.cxx"
#include "3FlavorAna/Vars/NumuVarsExtra.h"
//#include "3FlavorAna/Ana2024/joint_fit_2024_loader_tools.h"

#include "OscLib/OscCalcPMNSOpt.h"
#include "3FlavorAna/Ana2024/Predictions/plotting_header.h"
#include "3FlavorAna/Ana2024/joint_fit_2024_loader_tools.h"
#include "TRatioPlot.h"
#include <iomanip>
#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"


using namespace ana;

// inline std::vector <std::pair <Spectrum*, double> > GetNumuCosmics2024(const int nq = 4, std::string beam="rhc",
// 								  bool GetFromUPS=false, bool NERSC = false)
//   {

//     std::cout<<"\n============= You're loading " << beam << " Numu cosmics =============\n"<<std::endl;   

//     std::string dir = "/cvmfs/nova.osgstorage.org/analysis/3flavor/Ana2024/fit_inputs/v6_pothack/";
//     if(NERSC)  dir      = std::string(std::getenv("FCHELPERANA2020_LIB_PATH")) + "/Predictions/";
//     std::string filename = "cosm.nocnnhack.nova2024ana.v6.root";
    
//     TFile * fcosm = TFile::Open((dir+filename).c_str());

//     std::cout << "Loading " << beam << " Numu cosmics from the file at:" << std::endl;
//     std::cout << "\n" << dir + filename << "\n" << std::endl;

//     if(fcosm->IsZombie()) {std::cerr<< "bad cosmics\n"; exit(1);}
//     std::vector <std::pair <Spectrum*, double > > numu_cosmics;

//     double livetime;
//     if(beam=="rhc") livetime = RHCLivetime;
//     else            livetime = FHCLivetime;

//     for (int i = 1; i <=nq; ++i) {
//       Spectrum* cosm  = Spectrum::LoadFrom(fcosm, ("numu_"+beam+"_"+std::to_string(i)).c_str()).release();

//       double cosmFull  = cosm->Integral(cosm->Livetime(), 0, kLivetime);
//       double cosmScale = cosm->Integral(livetime        , 0, kLivetime);
    
//       std::cout << "Quantile " << i << ": A Total of " << cosmFull << " cosmics, which will scale to " << cosmScale << " in the spill window."
// 		<< "\n\t    The livetime is " << livetime <<", and the internal livetime is "<<cosm->Livetime()
// 		<< "\n\t    --> an uncertainty of " << 1/sqrt(cosmFull)
// 		<< "\n\t    ....Returning the spectrum.\n"
// 		<< std::endl;
    
    
//       numu_cosmics.push_back( { cosm, 1/sqrt(cosmFull) } );

//     }//end quantiles

//     for (size_t q=0; q<numu_cosmics.size(); ++q) {
//       std::cout << "Quantile " << q << ", Int = " << numu_cosmics[q].first->Integral(livetime        , 0, kLivetime) << std::endl;
//     }

//     return numu_cosmics;
//   }


void simple_pred()
{
  const std::string fname = "prod_caf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1";
  //const std::string fname = "prod_caf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";
  SpectrumLoader loader(fname);

  const Binning bins = Binning::Simple(50,0,5);

  const Cut kNumu2024CosRej_0p56(
    [](const caf::SRProxy* sr)
      {
        return (kNumuContPID(sr) > 0.56);
      });

  const Cut kThisCut_fc_0p56  = kNumuQuality && kNumu2024PID && !kNumuContainFD2024 
                                && kNumu2024CosRej_0p56 && k3flavor2024FDVeto 
                                && kNumuHadFracCut;


  // Beam prediction
  PredictionNoExtrap pred(loader, kNullLoader,
                          "Reconstructed Neutrino Energy (GeV)",
                          bins, kNumuE2024, kThisCut_fc_0p56);

  loader.Go();

  const double pot = kAna2024RHCPOT;
  const double livetime = kAna2024RHCLivetime;


  // Load cosmics (same beam mode as the file above: rhc)

  // auto numu_cosmics = GetNumuCosmics2024(4, "rhc");

  // // Merge into one spectrum
  // Spectrum cosmTot = *numu_cosmics[0].first;
  // for(size_t i=1; i<numu_cosmics.size(); ++i){
  //   cosmTot = cosmTot + *numu_cosmics[i].first;
  // }

  // Oscillation calculator from best-fit file

  std::string sHieOct   = "NOUO";
  std::string sFitDir   = "/exp/nova/data/groups/3flavor/Ana2024/Results/FrequentistFit/v8/NoFC/"; 
  std::string sFitFile  = sFitDir  + "bestfit_3Flavor_realdata_systs_numu_nue_nueLowE_rhcfhc_PtExtrap_Ana2024.root";

  osc::IOscCalcAdjustable* calc = nullptr;
  TFile *fBestFit = new TFile(sFitFile.c_str(), "READ");
  calc = dynamic_cast<osc::IOscCalcAdjustable*>(
            LoadFrom<osc::IOscCalc>(fBestFit, (sHieOct + "/calc").c_str()).release()
          );



  // Predictions: Beam only 
  Spectrum spred_beam  = pred.Predict(calc);


  // Plot beam + cosmics as separate histograms
  TCanvas* c1 = new TCanvas("c1","Beam + Cosmics",800,600);

  auto* h_beam  = spred_beam.ToTH1(pot); 
  // auto* h_cosm  = cosmTot.ToTH1(livetime, kLivetime);
  auto* h_unosc = pred.PredictUnoscillated().ToTH1(pot);

  h_beam->SetLineColor(kRed);
  // h_cosm->SetLineColor(kGreen+2);
  h_unosc->SetLineColor(kBlack);

  h_beam->Draw("hist");
  h_beam->GetYaxis()->SetRangeUser(0, 20);
  // h_cosm->Draw("hist same");
  h_unosc->Draw("hist same");

  TFile* fIn = new TFile("spectra1D_nowgt_rhc.root","READ");
  TH1* h_NuMuE = (TH1*)fIn->Get("h_NuMuE");
  if(h_NuMuE){
    h_NuMuE->SetLineColor(kGreen+2);
    h_NuMuE->SetLineWidth(2);
    h_NuMuE->Draw("hist same");
  }

  TLegend* leg = new TLegend(0.6,0.6,0.88,0.88);
  leg->AddEntry(h_beam,"Oscillated Beam","l");
  // leg->AddEntry(h_cosm,"Cosmic Background","l");
  leg->AddEntry(h_unosc,"Unoscillated Beam","l");
  leg->AddEntry(h_NuMuE, "Cosmic background", "l");
  leg->Draw("same");

  // Ratio: Oscillated / Unoscillated

  TH1* h_ratio = (TH1*)h_beam->Clone("h_ratio");
  // h_ratio->SetTitle("Oscillated / Unoscillated Beam Prediction");
  h_ratio->Divide(h_unosc);
  h_ratio->SetLineColor(kBlue);
  h_ratio->SetLineWidth(2);
  h_ratio->GetYaxis()->SetRangeUser(0, 2);
  h_ratio->GetYaxis()->SetTitle("Ratio to no Oscillation");
  h_ratio->GetXaxis()->SetTitle("Captured Reconstructed Neutrino Energy (GeV)");

  // Plot the ratio
  TCanvas* c3 = new TCanvas("c3","Osc to Unosc Ratio",800,600);
  h_ratio->Draw("hist");
  c3->SetGrid(0,0);
  c3->SaveAs("ratio_plot.png");
  c3->SaveAs("ratio_plot.pdf");

  // Save histograms to ROOT file
  TFile* fOut = new TFile("dumb_pred_test_output.root", "RECREATE");
  h_beam->Write("h_beam");
  h_unosc->Write("h_unosc");
  h_ratio->Write("h_ratio");
  fOut->Close();

  // Contour with beam only
  TCanvas* c2 = new TCanvas("c2","Contour",800,600);
  SingleSampleExperiment expt(&pred, spred_beam);
  FrequentistSurface surf(&expt, calc,
                          &kFitSinSqTheta23, 30, .35, .65,
                          &kFitDmSq32Scaled, 30, 2.2, 2.8);
  surf.Draw();
  surf.DrawContour(Gaussian68Percent2D(surf), 7, kRed); // dashed
  surf.DrawBestFit(kRed);
}
