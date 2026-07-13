//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmplModelEMp

\brief    MAID2007 electromagnetic helicity amplitudes on protons.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_MODEL_EM_P_H_
#define _MAID_HELICITY_AMPL_MODEL_EM_P_H_

#include <map>
#include <string>

#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelI.h"

namespace genie {

class FKR;
class RSHelicityAmplModelI;

typedef struct SMAIDAmpFit
{
  bool   has;
  double a0;
  double alpha;
  double beta;
}
MAIDAmpFit_t;

typedef struct SMAIDP11Fit
{
  double a0;
  double c1;
  double c4;
  double beta;
}
MAIDP11Fit_t;

typedef struct SMAIDTargetFit
{
  MAIDAmpFit_t a12;
  MAIDAmpFit_t a32;
  MAIDAmpFit_t s12;
}
MAIDTargetFit_t;

class MAIDHelicityAmplModelEMp : public MAIDHelicityAmplModelI {

public:
  MAIDHelicityAmplModelEMp();
  MAIDHelicityAmplModelEMp(string config);
  virtual ~MAIDHelicityAmplModelEMp();

  const MAIDHelicityAmpl & Compute(const Interaction * interaction) const;

  void Configure(const Registry & config);
  void Configure(string config);

protected:
  MAIDHelicityAmplModelEMp(string name, string config);

  void LoadConfig(void);
  bool ComputeMAID2007(Resonance_t res, double W, double Q2, double M) const;
  bool ComputeRSOriginal(Resonance_t res, double W, double Q2, double M) const;
  bool ComputeRSOriginalP33(Resonance_t res, double W, double Q2, double M) const;
  void ComputeEMFKR(Resonance_t res, double W, double Q2, double M,
                    FKR & fkr) const;
  double KCM(double W, double Q2, double M) const;
  double KGammaLab(double W, double M) const;
  void LoadMAID2007P11Fit(void);
  void LoadMAID2007P11AmpFit(const string & key, MAIDP11Fit_t & fit,
                             const MAIDP11Fit_t & def);
  void LoadMAID2007Fits(void);
  void LoadMAID2007Fit(Resonance_t res, const string & res_name,
                       const MAIDTargetFit_t & def);
  void LoadMAID2007AmpFit(const string & key, MAIDAmpFit_t & fit,
                          const MAIDAmpFit_t & def);

  mutable MAIDHelicityAmpl fAmpl;

  std::map<Resonance_t, MAIDTargetFit_t> fMAID2007Fits;
  MAIDP11Fit_t fP11A12;
  MAIDP11Fit_t fP11S12;
  const RSHelicityAmplModelI * fRSHelicityAmplModel;
  string fAmplitudeSource;
  string fRSQConvention;
  bool   fUseNeutronAmplitudes;
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
  double fRSVectorMass;
  double fRSOmega;
  double fRSZeta;
  double fRSOriginalPhase;
};

}        // genie namespace
#endif   // _MAID_HELICITY_AMPL_MODEL_EM_P_H_
