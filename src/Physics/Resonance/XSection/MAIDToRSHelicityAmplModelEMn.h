//____________________________________________________________________________
/*!
\class    genie::MAIDToRSHelicityAmplModelEMn

\brief    Neutron wrapper for MAID2007 A/S to RS-f Delta amplitudes.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_N_H_
#define _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_N_H_

#include "Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMp.h"

namespace genie {

class MAIDToRSHelicityAmplModelEMn : public MAIDToRSHelicityAmplModelEMp {

public:
  MAIDToRSHelicityAmplModelEMn();
  MAIDToRSHelicityAmplModelEMn(string config);
  virtual ~MAIDToRSHelicityAmplModelEMn();
};

}        // genie namespace
#endif   // _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_N_H_

