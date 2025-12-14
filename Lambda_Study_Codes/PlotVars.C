#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "NDAna/numubarcc_inc/NumubarCCIncVars.h"
#include "NDAna/numubarcc_inc/NumubarCCIncCuts.h"
#include "NDAna/numubarcc_inc/NumubarCCIncBins.h"
#include "NDAna/numubarcc_inc/NumubarCCIncSystDefs.h"

#include "TFile.h"

#include <string>

using namespace ana;

double POT = 12.5511e20;

const Cut kThisCut_sig = ana::xsec::numubarcc::kIsNumuCC_Lambda && 
                     ana::xsec::numubarcc::kFullSelectionCut &&
                     (ana::xsec::numubarcc::kRecoNPngs == 2 or ana::xsec::numubarcc::kRecoNPngs == 3)
                     ;

const Cut kThisCut_presel = ana::xsec::numubarcc::kFullSelectionCut &&
                          (ana::xsec::numubarcc::kRecoNPngs == 2 or ana::xsec::numubarcc::kRecoNPngs == 3)
                          ;

const Cut kThisCut_failedsig = !ana::xsec::numubarcc::kIsNumuCC_Lambda && 
                     ana::xsec::numubarcc::kFullSelectionCut &&
                     (ana::xsec::numubarcc::kRecoNPngs == 2 or ana::xsec::numubarcc::kRecoNPngs == 3)
                     ;

void PlotVars(string outfile = "lambda_study1.root"){
  
  SpectrumLoader * loader = new SpectrumLoader(ana::xsec::numubarcc::nominal_prod5_1_flatcaf.c_str());
  
  Spectrum * sRecoENu_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kRecoENuAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sRecoENu_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kRecoENuAxis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sRecoENu_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kRecoENuAxis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonKE_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonKEAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPiMinusKE_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPiMinusKEAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sNPngs_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kNPngsAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonCalE_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonCalEAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt); 
  Spectrum * sPiMinusCalE_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPiMinusCalEAxis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  
  Spectrum * sNHits_Png2_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kNHits_Png2Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sNHits_Png3_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kNHits_Png3Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID0_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID0Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID0_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID0Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID0_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID0Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);       

  Spectrum * sMuonID1_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID1Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID1_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID1Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID1_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID1Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);  

  Spectrum * sMuonID2_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID2Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID2_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID2Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID2_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID2Axis, 
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID0_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID0Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID0_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID0Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID0_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID0Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);  

  Spectrum * sMuonID1_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID1Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID1_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID1Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID1_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID1Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID2_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID2Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID2_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID2Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID2_presel = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID2Axis, 
    kThisCut_presel, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID0_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID0Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID0_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID0Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID0_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID0Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID1_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID1Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID1_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID1Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID1_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID1Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sMuonID2_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kMuon_ID2Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sProtonID2_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID2Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);
  Spectrum * sPionID2_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kPionID2Axis, 
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonID1_PionID1_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID1Axis, ana::xsec::numubarcc::kPionID1Axis,
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonID2_PionID2_sig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID2Axis, ana::xsec::numubarcc::kPionID2Axis,
    kThisCut_sig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonID1_PionID1_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID1Axis, ana::xsec::numubarcc::kPionID1Axis,
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  Spectrum * sProtonID2_PionID2_failedsig = new Spectrum(
    *loader, ana::xsec::numubarcc::kProtonID2Axis, ana::xsec::numubarcc::kPionID2Axis,
    kThisCut_failedsig, kNoShift, ana::xsec::numubarcc::std_wgt);

  loader->Go();

  TFile * fout = TFile::Open(outfile.c_str(), "RECREATE");
   sRecoENu_sig->SaveTo(fout, "sRecoENu_sig");
   sRecoENu_presel->SaveTo(fout, "sRecoENu_presel");
   sRecoENu_failedsig->SaveTo(fout, "sRecoENu_failedsig");
   sProtonKE_sig->SaveTo(fout, "sProtonKE_sig"); 
   sPiMinusKE_sig->SaveTo(fout, "sPiMinusKE_sig");
   sNPngs_sig->SaveTo(fout, "sNPngs_sig");
   sProtonCalE_sig->SaveTo(fout, "sProtonCalE_sig");
   sPiMinusCalE_sig->SaveTo(fout, "sPiMinusCalE_sig");
   sNHits_Png2_sig->SaveTo(fout, "sNHits_Png2_sig");
   sNHits_Png3_sig->SaveTo(fout, "sNHits_Png3_sig");
   sMuonID0_sig->SaveTo(fout, "sMuonID0_sig");
   sProtonID0_sig->SaveTo(fout, "sProtonID0_sig");
   sPionID0_sig->SaveTo(fout, "sPionID0_sig");
   sMuonID1_sig->SaveTo(fout, "sMuonID1_sig");
   sProtonID1_sig->SaveTo(fout, "sProtonID1_sig");
   sPionID1_sig->SaveTo(fout, "sPionID1_sig");
   sMuonID2_sig->SaveTo(fout, "sMuonID2_sig");
   sProtonID2_sig->SaveTo(fout, "sProtonID2_sig");
   sPionID2_sig->SaveTo(fout, "sPionID2_sig");
   sMuonID0_presel->SaveTo(fout, "sMuonID0_presel"); 
   sProtonID0_presel->SaveTo(fout, "sProtonID0_presel");
   sPionID0_presel->SaveTo(fout, "sPionID0_presel");
   sMuonID1_presel->SaveTo(fout, "sMuonID1_presel");
   sProtonID1_presel->SaveTo(fout, "sProtonID1_presel");
   sPionID1_presel->SaveTo(fout, "sPionID1_presel");
   sMuonID2_presel->SaveTo(fout, "sMuonID2_presel");
   sProtonID2_presel->SaveTo(fout, "sProtonID2_presel");
   sPionID2_presel->SaveTo(fout, "sPionID2_presel");
   sMuonID0_failedsig->SaveTo(fout, "sMuonID0_failedsig");
   sProtonID0_failedsig->SaveTo(fout, "sProtonID0_failedsig");
   sPionID0_failedsig->SaveTo(fout, "sPionID0_failedsig");
   sMuonID1_failedsig->SaveTo(fout, "sMuonID1_failedsig");
   sProtonID1_failedsig->SaveTo(fout, "sProtonID1_failedsig");
   sPionID1_failedsig->SaveTo(fout, "sPionID1_failedsig");
   sMuonID2_failedsig->SaveTo(fout, "sMuonID2_failedsig");
   sProtonID2_failedsig->SaveTo(fout, "sProtonID2_failedsig");
   sPionID2_failedsig->SaveTo(fout, "sPionID2_failedsig");
   sProtonID1_PionID1_sig->SaveTo(fout, "sProtonID1_PionID1_sig");
   sProtonID2_PionID2_sig->SaveTo(fout, "sProtonID2_PionID2_sig");
   sProtonID1_PionID1_failedsig->SaveTo(fout, "sProtonID1_PionID1_failedsig");
   sProtonID2_PionID2_failedsig->SaveTo(fout, "sProtonID2_PionID2_failedsig");

   fout->Close();


}
