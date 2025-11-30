#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/Cuts/QuantileCuts.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "3FlavorAna/Vars/Binnings.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Analysis/Exposures.h"
#include "3FlavorAna/Vars/NueVars.cxx"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"
#include "3FlavorAna/Vars/NumuVarsExtra.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "3FlavorAna/Ana2024/Predictions/syst_variations.h"

#include "TCanvas.h"                                                                                                                                                      
#include "TH2.h"
#include "TFile.h"
#include <iostream>

using namespace ana;

void EHadFrac_vs_RecoE_cosmics( bool isFHC = true, std::string period = "allperiods" ) {  
  std::string polarity = "";
  if ( isFHC) polarity = "fhc";
  if (!isFHC) polarity = "rhc";

  std::cout << "\n================================= \n"
	    << " make distribution for quantiles"
	    << "\n================================= \n"
	    << " polarity: " << polarity << ", period: " << period
	    << "\n================================= \n" 
	    << std::endl;
  
  // Which cut am I using?

  const Cut kNumu2024CosRej_0p56(
    [](const caf::SRProxy* sr)
      { 
        return (kNumuContPID(sr) > 0.56);
      });
 
  const Cut kThisCut_fc_0p56  = kNumuQuality && kNumu2024PID && !kNumuContainFD2024 
                                && kNumu2024CosRej_0p56 && k3flavor2024FDVeto 
                                && kNumuHadFracCut;
  // const Cut kThisCut_et_fc_0p56  = kNumuQuality && kNumu2024PID && kNumuProngsContainFD2024_et && !kNumuContainFD2024 && kNumu2024CosRej_0p56 && k3flavor2024FDVeto;

  //Declare a new variable to understand how many faces a track fails to be in the failed containment selection
  // const Var kNumuFaliedFaceCountFD = 
  //   Var([](const caf::SRProxy* sr) {
  //       int count = 0;
  //       if (kDistAllTop(sr)    <= 60) ++count;
  //       if (kDistAllBottom(sr) <= 12)  ++count;
  //       if (kDistAllEast(sr)   <= 16)  ++count;
  //       if (kDistAllWest(sr)   <= 12)  ++count;
  //       if (kDistAllFront(sr)  <= 18) ++count;
  //       if (kDistAllBack(sr)   <= 18)  ++count;
  //       return count;
  //   });

  // What XSec Weight am I using?
  const Weight weight = kXSecCVWgt2024;
  // Binning
  /*
  Binning kNearEdgeBins = Binning::Simple(100, -50., 50.);
  Binning kFullNearEdgeBins = Binning::Simple(1660, -60., 1600.);
  Binning kFFCBins = Binning::Simple(6,0,6); */

  Binning kMuEBins = Binning::Simple(50,0,5);

  /*
  std::string labelTop     = "Min. Distance to Top Face (cm), cosmics";
  std::string labelBottom  = "Minimum Distance to Bottom Detector Face (cm)";
  std::string labelEast    = "Minimum Distance to East Detector Face (cm)"; 
  std::string labelWest    = "Minimum Distance to West Detector Face (cm)";
  std::string labelFront   = "Minimum Distance to Front Detector Face (cm)";
  std::string labelBack    = "Minimum Distance to Back Detector Face (cm)";
  std::string labelfacecount    = "Failed face count (Cosmics)"; */

  // Figure out my definition, and then declare my spectrum loader.
  std::string MyDef = "prod_caf_R20-11-25-prod5.1reco.q.t_fd_cosmic_"+polarity+"_"+period+"_v1_goodruns_nostale"; 
  SpectrumLoader loader( MyDef );
  loader.SetSpillCut(kStandardSpillCuts);
  loader.SetCosmicHornCurrent(polarity); // for cosmics
  
  // Declare histaxis here
  /*
  const HistAxis kAxisTop    (labelTop,     kFullNearEdgeBins, kDistAllTop);
  const HistAxis kAxisBottom (labelBottom,  kNearEdgeBins, kDistAllBottom);
  const HistAxis kAxisEast   (labelEast,    kFullNearEdgeBins, kDistAllEast);
  const HistAxis kAxisWest   (labelWest,    kNearEdgeBins, kDistAllWest);
  const HistAxis kAxisFront  (labelFront,   kNearEdgeBins, kDistAllFront);
  const HistAxis kAxisBack   (labelBack,    kNearEdgeBins, kDistAllBack);
  const HistAxis kAxisFFC    (labelfacecount, kFFCBins, kNumuFaliedFaceCountFD); */

  const HistAxis kNumuEAxis("Reconstructed Neutrino Energy (GeV)", kMuEBins, kNumuE2024);
  const HistAxis kHadEFracAxis2024_0p56("E_{had.} / E_{#nu} (Cosmics)", Binning::Simple(200,0,1), kNumuHadEFrac2024);
  const HistAxis kHadEFracAxis2024_et_0p56("E_{had.} / E_{#nu}, PID > 0.56, w/ ET, Cosmics", Binning::Simple(200,0,1), kNumuHadEFrac2024);

  // Declare my spectra
/*
  Spectrum sTop_fc_ffc (loader, kAxisTop, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sBottom_fc_ffc (loader, kAxisBottom, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sEast_fc_ffc (loader, kAxisEast, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sWest_fc_ffc (loader, kAxisWest, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sFront_fc_ffc (loader, kAxisFront, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sBack_fc_ffc (loader, kAxisBack, kAxisFFC, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sEastTop_fc (loader, kAxisEast, kAxisTop, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight); 
  Spectrum sRecoE_HadEFrac_fc (loader, kNumuCCOptimisedAxis2024, kHadEFracAxis2024, kThisCut_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);
  Spectrum sRecoE_HadEFrac_et_fc (loader, kNumuCCOptimisedAxis2024, kHadEFracAxis2024, kThisCut_et_fc && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight); */

  Spectrum sMuonE_HadEFrac_fc_0p56 (loader, kNumuEAxis, kHadEFracAxis2024_0p56, kThisCut_fc_0p56 && kInCosmicTimingWindow , kNoShift, k3FAna2024Wgt);
  // Spectrum sMuonE_HadEFrac_et_fc_0p56 (loader, kNumuEAxis, kHadEFracAxis2024_et_0p56, kThisCut_et_fc_0p56 && kInCosmicTimingWindow , kNoShift, kPPFXFluxCVWgt*weight);

  // Set my loader going.
  loader.Go();

  // Convert spectra to TH2 histograms
/*
  TH2* hTop_fc_ffc    = sTop_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hBottom_fc_ffc = sBottom_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hEast_fc_ffc = sEast_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hWest_fc_ffc   = sWest_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hFront_fc_ffc = sFront_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hBack_fc_ffc = sBack_fc_ffc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hEastTop_fc = sEastTop_fc.ToTH2(kAna2024RHCLivetime, kLivetime);  
  TH2* hRecoE_HadEFrac_fc = sRecoE_HadEFrac_fc.ToTH2(kAna2024RHCLivetime, kLivetime);
  TH2* hRecoE_HadEFrac_et_fc = sRecoE_HadEFrac_et_fc.ToTH2(kAna2024RHCLivetime, kLivetime); */

  TH2* hMuonE_HadEFrac_fc_0p56 = sMuonE_HadEFrac_fc_0p56.ToTH2(sMuonE_HadEFrac_fc_0p56.Livetime(), kLivetime);
  // TH2* hMuonE_HadEFrac_et_fc_0p56 = sMuonE_HadEFrac_et_fc_0p56.ToTH2(sMuonE_HadEFrac_et_fc_0p56.Livetime(), kLivetime);

  //Make a single output file
  //std::string fileName = "quantiles_" + polarity + "_" + period + "_numu_spline_kHadEFracAxis2024_combined_nt_cos0.root";
  // std::string fileName = "spectra_cosmics_" + polarity + "_" + period + ".root";
  // TFile *OutFile = TFile::Open(fileName.c_str(), "RECREATE");

  const char* process_env = std::getenv("PROCESS");
  std::string process_str = process_env ? process_env : "0";
  std::string fileName = "spectra_cosmics_muE_cosRej_" + polarity + "_" + period + "_job" + process_str + ".root";
  TFile *OutFile = TFile::Open(fileName.c_str(), "RECREATE");

  // Write all histograms to the same file
/*
  hTop_fc_ffc->Write("Top_cos_fc_ffc");
  hBottom_fc_ffc->Write("Bottom_cos_fc_ffc");
  hEast_fc_ffc->Write("East_cos_fc_ffc");
  hWest_fc_ffc->Write("West_cos_fc_ffc");
  hFront_fc_ffc->Write("Front_cos_fc_ffc");
  hBack_fc_ffc->Write("Back_cos_fc_ffc");
  hEastTop_fc->Write("East_Top_cos_fc"); */
  sMuonE_HadEFrac_fc_0p56.SaveTo(OutFile, "sMuonE_HadEFrac_fc_cos_0p56");
  // sMuonE_HadEFrac_et_fc_0p56.SaveTo(OutFile, "sMuonE_HadEFrac_et_fc_cos_0p56");

  hMuonE_HadEFrac_fc_0p56->Write("hMuonE_HadEFrac_fc_cos_0p56");
  // hMuonE_HadEFrac_et_fc_0p56->Write("hMuonE_HadEFrac_et_fc_cos_0p56");

  OutFile->Close();


}
