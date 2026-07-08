//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmplModelI

\brief    Pure abstract base class for MAID helicity-amplitude models.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_I_H_
#define _MAID_HELICITY_AMPL_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"

namespace genie {

class Interaction;

class MAIDHelicityAmplModelI : public Algorithm
{
public:
  virtual ~MAIDHelicityAmplModelI();

  virtual const MAIDHelicityAmpl & Compute(const Interaction * interaction) const = 0;

protected:
  MAIDHelicityAmplModelI();
  MAIDHelicityAmplModelI(string name);
  MAIDHelicityAmplModelI(string name, string config);
};

}        // namespace

#endif   // _MAID_HELICITY_AMPL_MODEL_I_H_

