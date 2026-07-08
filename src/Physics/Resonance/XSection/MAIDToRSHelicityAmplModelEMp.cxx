//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMp.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDToRSHelicityAmplModelEMp::MAIDToRSHelicityAmplModelEMp() :
RSHelicityAmplModelI("genie::MAIDToRSHelicityAmplModelEMp")
{
  fMAIDHelicityAmplModel = 0;
  fW = 0.0;
  fQ2 = 0.0;
  fM = kProtonMass;
}
//____________________________________________________________________________
MAIDToRSHelicityAmplModelEMp::MAIDToRSHelicityAmplModelEMp(string config) :
RSHelicityAmplModelI("genie::MAIDToRSHelicityAmplModelEMp", config)
{
  fMAIDHelicityAmplModel = 0;
  fW = 0.0;
  fQ2 = 0.0;
  fM = kProtonMass;
}
//____________________________________________________________________________
MAIDToRSHelicityAmplModelEMp::MAIDToRSHelicityAmplModelEMp(
  string name, string config) :
RSHelicityAmplModelI(name, config)
{
  fMAIDHelicityAmplModel = 0;
  fW = 0.0;
  fQ2 = 0.0;
  fM = kProtonMass;
}
//____________________________________________________________________________
MAIDToRSHelicityAmplModelEMp::~MAIDToRSHelicityAmplModelEMp()
{

}
//____________________________________________________________________________
void MAIDToRSHelicityAmplModelEMp::SetKinematics(
  double W, double Q2, double M) const
{
  fW = W;
  fQ2 = Q2;
  fM = M;
}
//____________________________________________________________________________
const RSHelicityAmpl &
MAIDToRSHelicityAmplModelEMp::Compute(Resonance_t res, const FKR & /*fkr*/) const
{
  fAmpl.Set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  if(!fMAIDHelicityAmplModel) return fAmpl;
  if(fW <= 0.0 || fQ2 < 0.0 || fM <= 0.0) return fAmpl;

  const int hit_nuc = (this->Id().Name() == "genie::MAIDToRSHelicityAmplModelEMn") ?
    kPdgNeutron : kPdgProton;
  const int target = (hit_nuc == kPdgNeutron) ? kPdgTgtFreeN : kPdgTgtFreeP;

  Interaction * interaction =
    Interaction::RESEM(target, hit_nuc, kPdgElectron, 1.0);
  interaction->ExclTagPtr()->SetResonance(res);
  interaction->KinePtr()->SetW(fW);
  interaction->KinePtr()->SetQ2(fQ2);

  const MAIDHelicityAmpl & as = fMAIDHelicityAmplModel->Compute(interaction);
  delete interaction;
  if(!as.IsValid()) return fAmpl;

  const double mr = (res == kP33_1232) ? fDeltaMass : utils::res::Mass(res);
  const double k  = this->KCM(mr, fQ2, fM);
  const double kr = this->KCM(mr, 0.0, fM);
  if(k <= 0.0 || kr <= 0.0) return fAmpl;

  const double norm =
    TMath::Sqrt((fM/mr) * kr / (2.0 * kPi * kAem));
  const double fm1 = -norm * as.AmplA12();
  const double fm3 = -norm * as.AmplA32();
  const double f0p = norm * (fQ2/(k*k)) * as.AmplS12();

  fAmpl.Set(fm1, -fm1, fm3, -fm3, -f0p, f0p);
  return fAmpl;
}
//____________________________________________________________________________
double MAIDToRSHelicityAmplModelEMp::KCM(double W, double Q2, double M) const
{
  if(W <= 0.0) return 0.0;
  const double omega = (W*W - M*M - Q2) / (2.0 * W);
  return TMath::Sqrt(TMath::Max(0.0, omega*omega + Q2));
}
//____________________________________________________________________________
void MAIDToRSHelicityAmplModelEMp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDToRSHelicityAmplModelEMp::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDToRSHelicityAmplModelEMp::LoadConfig(void)
{
  RgAlg as_alg("genie::MAIDHelicityAmplModelEMp", "Default");
  this->GetParamDef("MAIDHelicityAmplAlg", as_alg, as_alg);
  fMAIDHelicityAmplModel =
    dynamic_cast<const MAIDHelicityAmplModelI *>(
      AlgFactory::Instance()->GetAlgorithm(as_alg.name, as_alg.config));

  this->GetParamDef("ApplyLuisWidthScale", fApplyLuisWidthScale, false);
  this->GetParamDef("LuisWidthScale", fLuisWidthScale, TMath::Sqrt(130.0/115.0));
  this->GetParamDef("DeltaMass", fDeltaMass, 1.232);
  this->GetParamDef("AM@P33(1232)", fAM0, 300.0);
  this->GetParamDef("AE@P33(1232)", fAE0, -6.37);
  this->GetParamDef("AS@P33(1232)", fAS0, -12.40);
  this->GetParamDef("BetaM@P33(1232)", fBetaM, 0.01);
  this->GetParamDef("BetaE@P33(1232)", fBetaE, -0.021);
  this->GetParamDef("BetaS@P33(1232)", fBetaS, 0.12);
  this->GetParamDef("GammaM@P33(1232)", fGammaM, 0.23);
  this->GetParamDef("GammaE@P33(1232)", fGammaE, 0.16);
  this->GetParamDef("GammaS@P33(1232)", fGammaS, 0.23);
  this->GetParamDef("DeltaS@P33(1232)", fDeltaS, 4.9);
  this->GetParamDef("DipoleMass2", fDipoleMass2, 0.71);
}
//____________________________________________________________________________
