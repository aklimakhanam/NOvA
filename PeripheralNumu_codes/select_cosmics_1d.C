#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Fit/Fit.h"
#include "CAFAna/Analysis/Plots.h"
#include "3FlavorAna/Plotting/NuePlotStyle.h"
#include "CAFAna/Fit/FrequentistSurface.h"
#include "CAFAna/Analysis/Style.h"
#include "CAFAna/Core/LabelsAndBins.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Experiment/Dmsq32Constraint.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/FC/FCSurface.h"
// #include "3FlavorAna/Prediction/PredictionSystJoint2018.h"
#include "3FlavorAna/Systs/3FlavorAna2024Systs.h"
#include "CAFAna/Prediction/PredictionCombinePeriods.h"
#include "3FlavorAna/Systs/NumuSysts.h"
#include "CAFAna/Vars/FitVars.h"
#include "OscLib/IOscCalc.h"
#include "3FlavorAna/Ana2024/joint_fit_2024_loader_tools.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "CAFAna/Vars/Vars.h"
#include "3FlavorAna/Vars/NueVars.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "CAFAna/Weights/GenieWeights.h"
#include "3FlavorAna/Ana2024/Predictions/syst_variations.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "3FlavorAna/Vars/Binnings.h"

#include "TCanvas.h"
#include "TBox.h"
#include "TColor.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGaxis.h"

#include <algorithm>
#include <vector>
#include <string>

using namespace ana;

void select_cosmics_1d (){

  const HistAxis axisNumu("Captured Reconstructed Neutrino Energy (GeV)", Binning::Simple(50,0,5), kNumuE2024);
  const HistAxis axisMuEnergy( "Captured Reconstructed Muon Energy (GeV)", Binning::Simple(50,0,5), kMuE);

  // Weight
  const Weight weights = k3FAna2024Wgt;

  // Definition
  std::string MyDef = "prod_caf_R20-11-25-prod5.1reco.q.t_fd_cosmic_rhc_allperiods_v1_goodruns_nostale";
  SpectrumLoader loader( MyDef );
  loader.SetSpillCut(kStandardSpillCuts);
  loader.SetCosmicHornCurrent("rhc");

  const Cut kNumu2024CosRej_0p56(
    [](const caf::SRProxy* sr)
      {
        return (kNumuContPID(sr) > 0.56);
      });

  // Cuts
  const Cut kThisCut_fc_0p56  = kNumuQuality && kNumu2024PID && !kNumuContainFD2024 && kNumu2024CosRej_0p56 && k3flavor2024FDVeto && kNumuHadFracCut;

  // Spectra
  Spectrum *s_NuMuE = new Spectrum(loader, axisNumu, kThisCut_fc_0p56 && kInCosmicTimingWindow, kNoShift, weights);
  Spectrum *s_MuE = new Spectrum(loader, axisMuEnergy, kThisCut_fc_0p56 && kInCosmicTimingWindow, kNoShift, weights);

  loader.Go();
  
  // Histograms
  TH1* h_MuE = s_MuE->ToTH1( kAna2024RHCLivetime , kLivetime);
  TH1* h_NuMuE = s_NuMuE->ToTH1( kAna2024RHCLivetime , kLivetime);
  
  // Save to file
  TFile* file = new TFile("spectra1D_nowgt_rhc_high_stat.root","recreate");
   s_MuE->SaveTo(file, "spectra_MuE_cosmics");
   s_NuMuE->SaveTo(file, "spectra_NuMuE_cosmics");
   h_MuE->Write("h_MuE");
   h_NuMuE->Write("h_NuMuE");

  file->Close();
}
