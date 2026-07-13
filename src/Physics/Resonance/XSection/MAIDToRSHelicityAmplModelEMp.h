//____________________________________________________________________________
/*!
\class    genie::MAIDToRSHelicityAmplModelEMp

\brief    Converts Delta-only MAID2007 electromagnetic A/S amplitudes to
          Rein-Sehgal-style f amplitudes for native RS cross-section tests.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_P_H_
#define _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_P_H_

#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

namespace genie {

class MAIDHelicityAmplModelI;

class MAIDToRSHelicityAmplModelEMp : public RSHelicityAmplModelI {

public:
  MAIDToRSHelicityAmplModelEMp();
  MAIDToRSHelicityAmplModelEMp(string config);
  virtual ~MAIDToRSHelicityAmplModelEMp();

  void SetKinematics(double W, double Q2, double M) const;
  const RSHelicityAmpl & Compute(Resonance_t res, const FKR & fkr) const;

  void Configure(const Registry & config);
  void Configure(string config);

protected:
  MAIDToRSHelicityAmplModelEMp(string name, string config);

  void LoadConfig(void);
  double KCM(double W, double Q2, double M) const;

  const MAIDHelicityAmplModelI * fMAIDHelicityAmplModel;
  mutable RSHelicityAmpl fAmpl;
  mutable double fW;
  mutable double fQ2;
  mutable double fM;

  bool   fApplyLuisWidthScale;
  double fLuisWidthScale;
  double fDeltaMass;
  double fAM0;
  double fAE0;
  double fAS0;
  double fBetaM;
  double fBetaE;
  double fBetaS;
  double fGammaM;
  double fGammaE;
  double fGammaS;
  double fDeltaS;
  double fDipoleMass2;
};

}        // genie namespace
#endif   // _MAID_TO_RS_HELICITY_AMPL_MODEL_EM_P_H_
