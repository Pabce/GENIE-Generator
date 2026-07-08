//____________________________________________________________________________
/*!
\class    genie::MAIDRESPXSec

\brief    MAID2007 electromagnetic resonance cross-section model.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_RES_PXSEC_H_
#define _MAID_RES_PXSEC_H_

#include <map>
#include <string>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResList.h"

namespace genie {

class MAIDHelicityAmpl;
class MAIDHelicityAmplModelI;
class XSecIntegratorI;

typedef struct SMAID2007ResParam
{
  double mass;
  double width;
  double beta_pi;
  double xr;
}
MAID2007ResParam_t;

class MAIDRESPXSec : public XSecAlgorithmI {

public:
  MAIDRESPXSec();
  MAIDRESPXSec(string config);
  virtual ~MAIDRESPXSec();

  double XSec(const Interaction * i, KinePhaseSpace_t k) const;
  double Integral(const Interaction * i) const;
  bool   ValidProcess(const Interaction * i) const;

  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);
  void SigmaTLDirectAS(const MAIDHelicityAmpl & ampl, double W, double Q2,
                       double M, Resonance_t res,
                       double & sigT, double & sigL) const;
  void SigmaTLRSFBridge(const MAIDHelicityAmpl & ampl, double W, double Q2,
                        double M, double & sigT, double & sigL) const;
  double KGammaLab(double W, double M) const;
  double QLab(double W, double Q2, double M) const;
  double QPiCM(double W, double M) const;
  double GammaRunning(Resonance_t res, double W, double M) const;
  double GammaPi(Resonance_t res, double W, double M) const;
  double GammaTotal(Resonance_t res, double W, double M) const;
  double BreitWigner(Resonance_t res, double W, double M) const;
  double BreitWignerMAIDW2(Resonance_t res, double W, double M) const;
  double BreitWignerRSOriginal(Resonance_t res, double W) const;
  double BreitWignerConventionCorrection(Resonance_t res, double W) const;
  bool   PassRSNativeWCut(Resonance_t res, double W) const;
  bool   SelectedResonance(Resonance_t res) const;
  double ResonanceMass(Resonance_t res) const;
  double ResonanceWidth(Resonance_t res) const;
  double ResonancePiBranching(Resonance_t res) const;
  double ResonanceXr(Resonance_t res) const;
  void   LoadMAID2007ResParams(void);
  void   LoadMAID2007ResParam(Resonance_t res, const string & res_name,
                              const MAID2007ResParam_t & def);

  const MAIDHelicityAmplModelI * fHAmplModelEMp;
  const MAIDHelicityAmplModelI * fHAmplModelEMn;
  const XSecAlgorithmI * fRSFallbackXSecModel;
  const XSecIntegratorI * fXSecIntegrator;

  BaryonResList fResList;
  std::map<Resonance_t, MAID2007ResParam_t> fMAID2007ResParams;
  string fXSecRoute;
  string fBreitWignerMode;
  string fBreitWignerConventionCorrection;
  string fMAIDRunningWidthMode;
  bool   fUseRSFallbackForMissingResonances;
  bool   fUseDRJoinScheme;
  bool   fApplyRSNativeWCut;
  bool   fApplyPiNBranchingToMAIDBW;
  double fWcut;
  double fXSecScaleEM;
  bool   fWghtBW;
  double fDeltaMass;
  double fDeltaWidth;
  double fDeltaXr;
  double fPionMass;
  double fN0ResMaxNWidths;
  double fN2ResMaxNWidths;
  double fGnResMaxNWidths;
};

}       // genie namespace

#endif  // _MAID_RES_PXSEC_H_
