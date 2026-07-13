//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <algorithm>
#include <cctype>
#include <cmath>

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Physics/Resonance/XSection/FKR.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"

using namespace genie;
using namespace genie::constants;

namespace {
  MAIDAmpFit_t no_fit(void)
  {
    MAIDAmpFit_t fit = { false, 0.0, 0.0, 0.0 };
    return fit;
  }

  MAIDAmpFit_t fit(double a0, double alpha, double beta)
  {
    MAIDAmpFit_t fit = { true, a0, alpha, beta };
    return fit;
  }

  MAIDAmpFit_t real_photon_fit(double a0)
  {
    MAIDAmpFit_t fit = { true, a0, 0.0, 0.0 };
    return fit;
  }

  MAIDP11Fit_t p11_fit(double a0, double c1, double c4, double beta)
  {
    MAIDP11Fit_t fit = { a0, c1, c4, beta };
    return fit;
  }

  double eval_fit(const MAIDAmpFit_t & fit, double Q2)
  {
    return fit.a0 * (1.0 + fit.alpha * Q2) * TMath::Exp(-fit.beta * Q2);
  }

  double eval_p11_fit(const MAIDP11Fit_t & fit, double Q2)
  {
    const double q8 = TMath::Power(Q2, 4);
    return fit.a0 * (1.0 + fit.c1 * Q2 + fit.c4 * q8) *
      TMath::Exp(-fit.beta * Q2);
  }

  bool has_q2_fit(const MAIDAmpFit_t & fit)
  {
    return fit.has && (fit.alpha != 0.0 || fit.beta != 0.0);
  }

  bool target_fit_has_q2_dependence(const MAIDTargetFit_t & target)
  {
    return has_q2_fit(target.a12) || has_q2_fit(target.a32) ||
           has_q2_fit(target.s12);
  }

  const double kQ2ZeroTol = 1.0e-12;

  string uppercase_copy(string value)
  {
    std::transform(value.begin(), value.end(), value.begin(), ::toupper);
    return value;
  }
}

//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp() :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp")
{
  fRSHelicityAmplModel = 0;
  fUseNeutronAmplitudes = false;
}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp(string config) :
MAIDHelicityAmplModelI("genie::MAIDHelicityAmplModelEMp", config)
{
  fRSHelicityAmplModel = 0;
  fUseNeutronAmplitudes = false;
}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::MAIDHelicityAmplModelEMp(string name, string config) :
MAIDHelicityAmplModelI(name, config)
{
  fRSHelicityAmplModel = 0;
  fUseNeutronAmplitudes = false;
}
//____________________________________________________________________________
MAIDHelicityAmplModelEMp::~MAIDHelicityAmplModelEMp()
{

}
//____________________________________________________________________________
const MAIDHelicityAmpl &
MAIDHelicityAmplModelEMp::Compute(const Interaction * interaction) const
{
  fAmpl.Reset();

  if(!interaction) return fAmpl;
  const Resonance_t res = interaction->ExclTag().Resonance();

  const Kinematics & kine = interaction->Kine();
  const double W  = kine.W();
  const double Q2 = kine.Q2();
  const double M  = interaction->InitState().Tgt().HitNucMass();

  if(W <= 0.0 || Q2 < 0.0 || M <= 0.0) return fAmpl;

  const string source = uppercase_copy(fAmplitudeSource);
  bool ok = false;
  if(source == "MAID2007") {
    ok = this->ComputeMAID2007(res, W, Q2, M);
  }
  else if(source == "RSORIGINAL") {
    ok = this->ComputeRSOriginal(res, W, Q2, M);
  }
  else if(source == "RSORIGINALP33") {
    ok = this->ComputeRSOriginalP33(res, W, Q2, M);
  }
  else {
    LOG("MAIDHelicity", pWARN)
      << "Unknown AmplitudeSource = " << fAmplitudeSource;
  }

  if(!ok) fAmpl.SetValid(false);

  return fAmpl;
}
//____________________________________________________________________________
bool MAIDHelicityAmplModelEMp::ComputeMAID2007(
  Resonance_t res, double /*W*/, double Q2, double M) const
{
  const double scale = fApplyLuisWidthScale ? fLuisWidthScale : 1.0;

  if(res != kP33_1232) {
    if(res == kP11_1440) {
      const double A12 = eval_p11_fit(fP11A12, Q2);
      const double S12 = eval_p11_fit(fP11S12, Q2);
      fAmpl.SetAS(A12 * 1.0e-3 * scale, 0.0, S12 * 1.0e-3 * scale);
      return true;
    }

    std::map<Resonance_t, MAIDTargetFit_t>::const_iterator fit_iter =
      fMAID2007Fits.find(res);
    if(fit_iter == fMAID2007Fits.end()) return false;

    const MAIDTargetFit_t & target = fit_iter->second;
    if(Q2 > kQ2ZeroTol && !target_fit_has_q2_dependence(target)) return false;
    if(!target.a12.has && !target.a32.has && !target.s12.has) return false;

    const double A12 = target.a12.has ? eval_fit(target.a12, Q2) : 0.0;
    const double A32 = target.a32.has ? eval_fit(target.a32, Q2) : 0.0;
    const double S12 = target.s12.has ? eval_fit(target.s12, Q2) : 0.0;
    fAmpl.SetAS(A12 * 1.0e-3 * scale,
                A32 * 1.0e-3 * scale,
                S12 * 1.0e-3 * scale);
    return true;
  }

  const double k  = this->KCM(fDeltaMass, Q2, M);
  const double kr = this->KCM(fDeltaMass, 0.0, M);
  if(k <= 0.0 || kr <= 0.0) return false;

  const double gd  = 1.0 / TMath::Power(1.0 + Q2/fDipoleMass2, 2);
  const double tau = Q2 / (4.0 * M * M);

  const double AM =
    fAM0 * (1.0 + fBetaM * Q2) * TMath::Exp(-fGammaM * Q2) * gd * (k/kr);
  const double AE =
    fAE0 * (1.0 + fBetaE * Q2) * TMath::Exp(-fGammaE * Q2) * gd * (k/kr);
  const double AS =
    fAS0 * (1.0 + fBetaS * Q2) / (1.0 + fDeltaS * tau) *
    TMath::Exp(-fGammaS * Q2) * gd * (k*k/(kr*kr));

  const double A12 = -0.5 * (AM + 3.0*AE) * 1.0e-3 * scale;
  const double A32 =  0.5 * kSqrt3 * (AE - AM) * 1.0e-3 * scale;
  const double S12 = -kSqrt2 * AS * 1.0e-3 * scale;
  fAmpl.SetAS(A12, A32, S12);

  const double norm =
    TMath::Sqrt((M/fDeltaMass) * kr / (2.0 * kPi * kAem));
  const double qstar2 = k*k;

  const double fm1 = -norm * A12;
  const double fm3 = -norm * A32;
  double f0p = 0.0;
  if(qstar2 > 0.0) f0p = norm * (Q2/qstar2) * S12;

  fAmpl.SetRSF(fm1, -fm1, fm3, -fm3, -f0p, f0p);
  return true;
}
//____________________________________________________________________________
bool MAIDHelicityAmplModelEMp::ComputeRSOriginal(
  Resonance_t res, double W, double Q2, double M) const
{
  if(!fRSHelicityAmplModel) return false;
  if(W <= 0.0 || Q2 < 0.0 || M <= 0.0) return false;

  FKR fkr;
  this->ComputeEMFKR(res, W, Q2, M, fkr);
  fRSHelicityAmplModel->SetKinematics(W, Q2, M);
  const RSHelicityAmpl & rs = fRSHelicityAmplModel->Compute(res, fkr);

  const double mr = utils::res::Mass(res);
  const double kr = this->KCM(mr, 0.0, M);
  const double qstar = this->KCM(mr, Q2, M);
  if(mr <= 0.0 || kr <= 0.0 || qstar <= 0.0) return false;

  const double norm =
    TMath::Sqrt((M/mr) * kr / (2.0 * kPi * kAem));
  if(norm <= 0.0) return false;

  const double A12 = -rs.AmpMinus1() / norm;
  const double A32 = -rs.AmpMinus3() / norm;
  double S12 = 0.0;
  if(Q2 > 0.0) S12 = rs.Amp0Plus() / (norm * (Q2/(qstar*qstar)));

  fAmpl.SetAS(A12, A32, S12);
  fAmpl.SetRSF(rs.AmpMinus1(), rs.AmpPlus1(), rs.AmpMinus3(),
               rs.AmpPlus3(), rs.Amp0Minus(), rs.Amp0Plus());
  return true;
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::ComputeEMFKR(
  Resonance_t res, double W, double Q2, double M, FKR & fkr) const
{
  const int IR = utils::res::ResonanceIndex(res);
  const double q2 = -Q2;
  const double W2 = W*W;
  const double M2 = M*M;
  const double k = 0.5 * (W2 - M2) / M;
  const double v = k - 0.5 * q2 / M;
  const double q3sq = v*v - q2;
  const double q3 = TMath::Sqrt(TMath::Max(0.0, q3sq));

  const double go = TMath::Power(1.0 - 0.25*q2/M2, 0.5 - IR);
  const double gv = go * TMath::Power(1.0/(1.0 - q2/(fRSVectorMass*fRSVectorMass)), 2);

  const double d = TMath::Power(W + M, 2.0) - q2;
  const double sq2omg = TMath::Sqrt(2.0/fRSOmega);
  const double mq_w = M*q3/W;

  fkr.Reset();
  fkr.Lamda  = sq2omg * mq_w;
  fkr.Tv     = gv / (3.0*W*sq2omg);
  fkr.Rv     = kSqrt2 * mq_w*(W + M)*gv / d;
  if(q3sq > 0.0) {
    fkr.S = (-q2/q3sq) * (3.0*W*M + q2 - M2) * gv / (6.0*M2);
  }
  fkr.Ta     = 0.0;
  fkr.Ra     = 0.0;
  fkr.B      = 0.0;
  fkr.C      = 0.0;
  fkr.R      = fkr.Rv;
  fkr.T      = fkr.Tv;
  fkr.Tplus  = -fkr.Tv;
  fkr.Tminus = -fkr.Tv;
  fkr.Rplus  = -fkr.Rv;
  fkr.Rminus = -fkr.Rv;
}
//____________________________________________________________________________
bool MAIDHelicityAmplModelEMp::ComputeRSOriginalP33(
  Resonance_t res, double W, double Q2, double M) const
{
  if(res != kP33_1232) return false;

  double qrs = 0.0;
  const string qconv = uppercase_copy(fRSQConvention);
  if(qconv == "QSTAR") {
    qrs = this->KCM(W, Q2, M);
  }
  else {
    const double nu_lab = (W*W - M*M + Q2) / (2.0 * M);
    qrs = TMath::Sqrt(TMath::Max(0.0, Q2 + nu_lab*nu_lab));
  }

  const double gv =
    TMath::Sqrt(1.0 + Q2/(4.0*M*M)) /
    TMath::Power(1.0 + Q2/(fRSVectorMass*fRSVectorMass), 2);
  const double denom = TMath::Power(W + M, 2) + Q2;
  const double R =
    fRSOriginalPhase * kSqrt2 * (M/W) * ((W + M) * qrs / denom) * gv;

  const double fm1 = -kSqrt2 * R;
  const double fm3 = -kSqrt6 * R;
  fAmpl.SetRSF(fm1, -fm1, fm3, -fm3, 0.0, 0.0);

  const double kr = this->KCM(fDeltaMass, 0.0, M);
  if(kr <= 0.0) return true;
  const double norm =
    TMath::Sqrt((M/fDeltaMass) * kr / (2.0 * kPi * kAem));
  fAmpl.SetAS(-fm1/norm, -fm3/norm, 0.0);
  return true;
}
//____________________________________________________________________________
double MAIDHelicityAmplModelEMp::KCM(double W, double Q2, double M) const
{
  if(W <= 0.0) return 0.0;
  const double omega = (W*W - M*M - Q2) / (2.0 * W);
  return TMath::Sqrt(TMath::Max(0.0, omega*omega + Q2));
}
//____________________________________________________________________________
double MAIDHelicityAmplModelEMp::KGammaLab(double W, double M) const
{
  if(M <= 0.0) return 0.0;
  return (W*W - M*M) / (2.0 * M);
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadMAID2007P11Fit(void)
{
  const MAIDP11Fit_t a12_def =
    fUseNeutronAmplitudes ?
    p11_fit(54.1, 0.95, 0.0, 1.77) :
    p11_fit(-61.4, -1.22, -0.55, 1.51);
  const MAIDP11Fit_t s12_def =
    fUseNeutronAmplitudes ?
    p11_fit(-41.5, 2.98, 0.0, 1.55) :
    p11_fit(4.2, 40.0, 1.5, 1.75);

  this->LoadMAID2007P11AmpFit("A12@P11(1440)", fP11A12, a12_def);
  this->LoadMAID2007P11AmpFit("S12@P11(1440)", fP11S12, s12_def);
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadMAID2007P11AmpFit(
  const string & key, MAIDP11Fit_t & fit, const MAIDP11Fit_t & def)
{
  this->GetParamDef(key + "-A0",   fit.a0,   def.a0);
  this->GetParamDef(key + "-C1",   fit.c1,   def.c1);
  this->GetParamDef(key + "-C4",   fit.c4,   def.c4);
  this->GetParamDef(key + "-Beta", fit.beta, def.beta);
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadMAID2007Fits(void)
{
  fMAID2007Fits.clear();

  const MAIDTargetFit_t d13 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(-77.0, -0.53, 1.55), fit(-154.0, 0.58, 1.75),
                     fit(13.6, 15.7, 1.57) } :
    MAIDTargetFit_t{ fit(-27.0, 7.77, 1.09), fit(161.0, 0.69, 2.10),
                     fit(-63.6, 4.19, 3.40) };
  const MAIDTargetFit_t s11_1535 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(-51.0, 4.75, 1.69), no_fit(),
                     fit(28.5, 0.36, 1.55) } :
    MAIDTargetFit_t{ fit(66.0, 1.61, 0.70), no_fit(),
                     fit(-2.0, 23.9, 0.81) };
  const MAIDTargetFit_t s31 =
    { fit(66.0, 1.86, 2.50), no_fit(), fit(16.2, 2.83, 2.00) };
  const MAIDTargetFit_t s11_1650 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(9.0, 0.13, 1.55), no_fit(),
                     fit(10.1, -0.50, 1.55) } :
    MAIDTargetFit_t{ fit(33.0, 1.45, 0.62), no_fit(),
                     fit(-3.5, 2.88, 0.76) };
  const MAIDTargetFit_t d15 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(-62.0, 0.01, 2.00), fit(-84.0, 0.01, 2.00),
                     fit(0.0, 0.0, 0.0) } :
    MAIDTargetFit_t{ fit(15.0, 0.10, 2.00), fit(22.0, 0.10, 2.00),
                     fit(0.0, 0.0, 0.0) };
  const MAIDTargetFit_t f15 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(28.0, 0.0, 1.20), fit(-38.0, 4.09, 1.75),
                     fit(0.0, 0.0, 0.0) } :
    MAIDTargetFit_t{ fit(-25.0, 3.98, 1.20), fit(134.0, 1.00, 2.22),
                     fit(-44.0, 3.14, 1.68) };
  const MAIDTargetFit_t d33 =
    { fit(226.0, 1.91, 1.77), fit(210.0, 1.97, 2.20),
      fit(0.0, 0.0, 0.0) };
  const MAIDTargetFit_t p13 = fUseNeutronAmplitudes ?
    MAIDTargetFit_t{ fit(-3.0, 12.7, 1.55), fit(-31.0, 4.99, 1.55),
                     fit(0.0, 0.0, 0.0) } :
    MAIDTargetFit_t{ fit(73.0, 1.89, 1.55), fit(-11.0, 16.0, 1.55),
                     fit(-53.0, 2.46, 1.55) };
  const MAIDTargetFit_t f35 =
    { real_photon_fit(18.0), real_photon_fit(-28.0), fit(0.0, 0.0, 0.0) };
  const MAIDTargetFit_t f37 =
    { real_photon_fit(-94.0), real_photon_fit(-121.0), fit(0.0, 0.0, 0.0) };

  this->LoadMAID2007Fit(kD13_1520, "D13(1520)", d13);
  this->LoadMAID2007Fit(kS11_1535, "S11(1535)", s11_1535);
  this->LoadMAID2007Fit(kS31_1620, "S31(1620)", s31);
  this->LoadMAID2007Fit(kS11_1650, "S11(1650)", s11_1650);
  this->LoadMAID2007Fit(kD15_1675, "D15(1675)", d15);
  this->LoadMAID2007Fit(kF15_1680, "F15(1680)", f15);
  this->LoadMAID2007Fit(kD33_1700, "D33(1700)", d33);
  this->LoadMAID2007Fit(kP13_1720, "P13(1720)", p13);
  this->LoadMAID2007Fit(kF35_1905, "F35(1905)", f35);
  this->LoadMAID2007Fit(kF37_1950, "F37(1950)", f37);
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadMAID2007Fit(
  Resonance_t res, const string & res_name, const MAIDTargetFit_t & def)
{
  MAIDTargetFit_t target;
  this->LoadMAID2007AmpFit("A12@" + res_name, target.a12, def.a12);
  this->LoadMAID2007AmpFit("A32@" + res_name, target.a32, def.a32);
  this->LoadMAID2007AmpFit("S12@" + res_name, target.s12, def.s12);
  fMAID2007Fits[res] = target;
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadMAID2007AmpFit(
  const string & key, MAIDAmpFit_t & fit, const MAIDAmpFit_t & def)
{
  this->GetParamDef(key + "-Has",   fit.has,   def.has);
  this->GetParamDef(key + "-A0",    fit.a0,    def.a0);
  this->GetParamDef(key + "-Alpha", fit.alpha, def.alpha);
  this->GetParamDef(key + "-Beta",  fit.beta,  def.beta);
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDHelicityAmplModelEMp::LoadConfig(void)
{
  this->GetParamDef("AmplitudeSource", fAmplitudeSource, string("MAID2007"));
  this->GetParamDef("RSQConvention", fRSQConvention, string("RSLab"));
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
  this->GetParamDef("RSVectorMass", fRSVectorMass, 0.84);
  this->GetParamDef("RES-Omega", fRSOmega, 1.05);
  this->GetParamDef("RES-Zeta", fRSZeta, 0.76338);
  this->GetParamDef("RSOriginalPhase", fRSOriginalPhase, -1.0);

  RgAlg rs_alg(
    fUseNeutronAmplitudes ?
    "genie::RSHelicityAmplModelEMn" : "genie::RSHelicityAmplModelEMp",
    "Default");
  this->GetParamDef("RSHelicityAmplAlg", rs_alg, rs_alg);
  fRSHelicityAmplModel =
    dynamic_cast<const RSHelicityAmplModelI *>(
      AlgFactory::Instance()->GetAlgorithm(rs_alg.name, rs_alg.config));

  this->LoadMAID2007P11Fit();
  this->LoadMAID2007Fits();
}
//____________________________________________________________________________
