#pragma once

#include "CAFAna/Core/EventList.h"

#include "3FlavorAna/Cuts/NueCuts2024.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/Cuts/QuantileCuts2024.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "3FlavorAna/Vars/NumuVars.h"
#include "CAFAna/Vars/Vars.h"

#include "StandardRecord/Proxy/SRProxy.h"
#include "3FlavorAna/Vars/NumuVars.h"

using namespace ana;

void GetEventList(){

// integer variables 
  std::vector<const Var*> intVars = {&kRun, &kSubrun, &kEvt, &kSlc};

  const Cut kNumu2024CosRej_0p56(
    [](const caf::SRProxy* sr)
      {
        return (kNumuContPID(sr) > 0.56);
      });

  std::vector<Cut> numuRHCcuts = {kIsNumuCC && kNumuQuality && kNumu2024PID && !kNumuContainFD2024 
                                && kNumu2024CosRej_0p56 && k3flavor2024FDVeto 
                                && kNumuHadFracCut};

  std::vector<std::string> numuRHClist = {"numu_rhc_events_list_optimized.txt"};

  MakeTextListFile("prod_caf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1",numuRHCcuts,numuRHClist,{},intVars,&kStandardSpillCuts);

}
