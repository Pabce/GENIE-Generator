//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmplModelEMn

\brief    MAID2007 Delta-only electromagnetic helicity amplitudes on neutrons.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_EM_N_H_
#define _MAID_HELICITY_AMPL_MODEL_EM_N_H_

#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h"

namespace genie {

class MAIDHelicityAmplModelEMn : public MAIDHelicityAmplModelEMp {

public:
  MAIDHelicityAmplModelEMn();
  MAIDHelicityAmplModelEMn(string config);
  virtual ~MAIDHelicityAmplModelEMn();
};

}        // genie namespace
#endif   // _MAID_HELICITY_AMPL_MODEL_EM_N_H_

