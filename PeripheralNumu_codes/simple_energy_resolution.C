#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/Cuts/QuantileCuts.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "3FlavorAna/Vars/Binnings.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Analysis/Exposures.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"
#include "3FlavorAna/Vars/NueVars.cxx"
#include "3FlavorAna/Vars/NumuVarsExtra.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "3FlavorAna/Ana2024/Predictions/syst_variations.h"
#include "TCanvas.h"                                                                                                                                                      
#include "TH2.h"
#include "TFile.h"
#include <iostream>

using namespace ana;

const Var kTrueHadE = kTrueE - kTrueMuonE;

#define res(recoE, trueE) (recoE - (trueE))/(trueE)

void simple_energy_resolution( bool isFHC = true, std::string period = "full" ) {  
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
  
  const Cut kThisCut = kNumuQuality && kNumu2024PID && !kNumuContainFD2024 && kNumu2024CosRej_0p56 && k3flavor2024FDVeto 
                                && kNumuHadFracCut;

  
  // Binning
  const Binning Ebins = Binning::Simple(200, 0.0, 5.0);
  const Binning HadEbins = Binning::Simple(200, 0.0, 3.0);
  const Binning ResEbins = Binning::Simple(200, -1.0, 1.0);
  const Binning MuResEbins = Binning::Simple(200, -0.4, 0.4);
  const Binning NumuResEbins = Binning::Simple(200, -0.7, 0.7);
 
 // Figure out my definition, and then declare my spectrum loader.
  std::string MyDef = "prod_caf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_"+polarity+"_nova_v08_"+period+"_v1";
  SpectrumLoader loader( MyDef );
  loader.SetSpillCut(kStandardSpillCuts);
  
  // Declare histaxis here
  const HistAxis kTrueEAxis("TrueE (GeV)", Ebins, kTrueE);
  const HistAxis kRecoEAxis("RecoE (GeV)", Ebins, kNumuE2024);
  const HistAxis kMuEAxis("MuE (GeV)", Ebins, kNumuMuE2024);
  const HistAxis kHadEAxis("HadE (GeV)", HadEbins, kNumuHadE2024);
  const HistAxis kNumuEResolutionAxis("NumuE Resolution", NumuResEbins, res(kNumuE2024, kTrueE));
  const HistAxis kMuEResolutionAxis("MuE Resolution", MuResEbins, res(kNumuMuE2024, kTrueMuonE)); 
  const HistAxis kHadEResolutionAxis("HadE Resolution", ResEbins, res(kNumuHadE2024, kTrueHadE));

  // Declare my spectra
  Spectrum sNumuEResolution (loader, kNumuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);
  Spectrum sMuEResolution (loader, kMuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);  
  Spectrum sHadEResolution (loader, kHadEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt); 
  Spectrum sNumuEResolution_TrueE (loader, kTrueEAxis, kNumuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);
  Spectrum sNumuEResolution_RecoE (loader, kRecoEAxis, kNumuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);
  Spectrum sNumuEResolution_MuE (loader, kMuEAxis, kNumuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);
  Spectrum sNumuEResolution_HadE (loader, kHadEAxis, kNumuEResolutionAxis, kThisCut, kNoShift, k3FAna2024Wgt);
  Spectrum sTrueE_RecoE (loader, kTrueEAxis, kRecoEAxis, kThisCut, kNoShift, k3FAna2024Wgt);

  // Set my loader going.
  loader.Go();

  // Convert spectra to TH2 histograms
  TH1* hNumuEResolution  = sNumuEResolution .ToTH1(kAna2024RHCPOT);
  TH1* hMuEResolution  = sMuEResolution .ToTH1(kAna2024RHCPOT);
  TH1* hHadEResolution  = sHadEResolution .ToTH1(kAna2024RHCPOT);
  TH2* hNumuEResolution_TrueE  = sNumuEResolution_TrueE .ToTH2(kAna2024RHCPOT);
  TH2* hNumuEResolution_RecoE  = sNumuEResolution_RecoE .ToTH2(kAna2024RHCPOT); 
  TH2* hNumuEResolution_MuE  = sNumuEResolution_MuE .ToTH2(kAna2024RHCPOT);  
  TH2* hNumuEResolution_HadE  = sNumuEResolution_HadE .ToTH2(kAna2024RHCPOT);
  TH2* hTrueE_RecoE = sTrueE_RecoE .ToTH2(kAna2024RHCPOT);
  
  // Make a single output file
  std::string fileName = "energy_resolution_" + polarity + "_" + period + ".root";
  TFile *OutFile = TFile::Open(fileName.c_str(), "RECREATE");
  /*
  const char* process_env = std::getenv("PROCESS");
  std::string process_str = process_env ? process_env : "0";
  std::string fileName = "energy_resolution_" + polarity + "_" + period + "_job" + process_str + ".root";
  TFile *OutFile = TFile::Open(fileName.c_str(), "RECREATE");
  */

  // Write all histograms to the same file

  hNumuEResolution->Write("hNumuEResolution");
  hMuEResolution->Write("hMuEResolution");
  hHadEResolution->Write("hHadEResolution");
  hNumuEResolution_TrueE->Write("hNumuEResolution_TrueE");  
  hNumuEResolution_RecoE->Write("hNumuEResolution_RecoE");
  hNumuEResolution_MuE->Write("hNumuEResolution_MuE"); 
  hNumuEResolution_HadE->Write("hNumuEResolution_HadE");
  hTrueE_RecoE->Write("hTrueE_RecoE");

  OutFile->Close();


}
