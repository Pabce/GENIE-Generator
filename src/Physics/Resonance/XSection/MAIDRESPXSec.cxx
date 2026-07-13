//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <algorithm>
#include <cctype>

#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/MAIDRESPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;
using namespace genie::constants;

namespace {
  const char * kMAID2007ResonanceNameList =
    "P33(1232),P11(1440),D13(1520),S11(1535),S31(1620),S11(1650),"
    "D15(1675),F15(1680),D33(1700),P13(1720),F35(1905),P31(1910),"
    "F37(1950)";

  MAID2007ResParam_t maid_res_param(
    double mass, double width, double beta_pi, double xr)
  {
    MAID2007ResParam_t param = { mass, width, beta_pi, xr };
    return param;
  }

  string uppercase_copy(string value)
  {
    std::transform(value.begin(), value.end(), value.begin(), ::toupper);
    return value;
  }
}

//____________________________________________________________________________
MAIDRESPXSec::MAIDRESPXSec() :
XSecAlgorithmI("genie::MAIDRESPXSec")
{

}
//____________________________________________________________________________
MAIDRESPXSec::MAIDRESPXSec(string config) :
XSecAlgorithmI("genie::MAIDRESPXSec", config)
{

}
//____________________________________________________________________________
MAIDRESPXSec::~MAIDRESPXSec()
{

}
//____________________________________________________________________________
double MAIDRESPXSec::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction->InitState();
  const Kinematics & kine = interaction->Kine();
  const Target & target = init_state.Tgt();

  const double W  = kine.W();
  const double Q2 = kine.Q2();
  const double W2 = W*W;
  const double M  = target.HitNucMass();
  const double M2 = M*M;
  const double E  = init_state.ProbeE(kRfHitNucRest);
  const Resonance_t res = interaction->ExclTag().Resonance();

  if(W <= 0.0 || Q2 <= 0.0 || M <= 0.0 || E <= 0.0) return 0.0;

  if(fUseDRJoinScheme && W >= fWcut) return 0.0;
  if(fApplyRSNativeWCut && !this->PassRSNativeWCut(res, W)) return 0.0;

  const bool is_p = pdg::IsProton(target.HitNucPdg());
  const MAIDHelicityAmplModelI * ampl_model =
    is_p ? fHAmplModelEMp : fHAmplModelEMn;
  if(!ampl_model) return 0.0;

  const MAIDHelicityAmpl & ampl = ampl_model->Compute(interaction);
  if(!ampl.IsValid()) {
    if(fUseRSFallbackForMissingResonances && fRSFallbackXSecModel) {
      return fRSFallbackXSecModel->XSec(interaction, kps);
    }
    return 0.0;
  }

  double sigT = 0.0;
  double sigL = 0.0;
  const string route = uppercase_copy(fXSecRoute);
  if(route == "DIRECTAS") {
    this->SigmaTLDirectAS(ampl, W, Q2, M, res, sigT, sigL);
  }
  else if(route == "RSFBRIDGE") {
    if(res != kP33_1232) {
      if(fUseRSFallbackForMissingResonances && fRSFallbackXSecModel) {
        return fRSFallbackXSecModel->XSec(interaction, kps);
      }
      return 0.0;
    }
    this->SigmaTLRSFBridge(ampl, W, Q2, M, sigT, sigL);
  }
  else {
    LOG("MAIDRESPXSec", pWARN) << "Unknown XSecRoute = " << fXSecRoute;
    return 0.0;
  }

  if(sigT <= 0.0 && sigL <= 0.0) return 0.0;

  const double nu = (W2 - M2 + Q2) / (2.0 * M);
  const double Eprime = E - nu;
  if(Eprime <= 0.0) return 0.0;

  const double sin2_theta_2 = Q2 / (4.0 * E * Eprime);
  if(sin2_theta_2 <= 0.0 || sin2_theta_2 >= 1.0) return 0.0;

  const double tan2_theta_2 = sin2_theta_2 / (1.0 - sin2_theta_2);
  const double q3sq = nu*nu + Q2;
  const double epsilon = 1.0 / (1.0 + 2.0 * (q3sq/Q2) * tan2_theta_2);
  if(epsilon >= 1.0) return 0.0;

  const double kgamma = this->KGammaLab(W, M);
  if(kgamma <= 0.0) return 0.0;

  const double gamma_flux =
    kAem / (2.0 * kPi2) * (Eprime/E) * kgamma / (Q2 * (1.0 - epsilon));
  double xsec = gamma_flux * (sigT + epsilon * sigL);

  const double jac = W * kPi / (E * Eprime * M);
  xsec *= jac;
  xsec = TMath::Max(0.0, xsec);

  xsec *= fXSecScaleEM;

  if(kps != kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction, kPSWQ2fE, kps);
    xsec *= J;
  }

  if(interaction->TestBit(kIAssumeFreeNucleon)) return xsec;

  const int nnucl = is_p ? target.Z() : target.N();
  xsec *= nnucl;

  return xsec;
}
//____________________________________________________________________________
void MAIDRESPXSec::SigmaTLDirectAS(
  const MAIDHelicityAmpl & ampl, double W, double Q2, double M,
  Resonance_t res,
  double & sigT, double & sigL) const
{
  sigT = 0.0;
  sigL = 0.0;

  const double mr = this->ResonanceMass(res);
  const double kgamma = this->KGammaLab(W, M);
  const double kr = this->KGammaLab(mr, M);
  const double qr = this->QLab(mr, Q2, M);
  if(kgamma <= 0.0 || kr <= 0.0 || qr <= 0.0) return;

  const double bw = fWghtBW ? this->BreitWigner(res, W, M) : 1.0;
  const double common = kPi * (kr/kgamma) * (2.0 * M) * bw;

  sigT = common * (ampl.Ampl2A12() + ampl.Ampl2A32());
  sigL =
    2.0 * common * (Q2/(qr*qr)) *
    ((mr*mr)/(M*M)) * ampl.Ampl2S12();
}
//____________________________________________________________________________
void MAIDRESPXSec::SigmaTLRSFBridge(
  const MAIDHelicityAmpl & ampl, double W, double Q2, double M,
  double & sigT, double & sigL) const
{
  sigT = 0.0;
  sigL = 0.0;

  const double kgamma = this->KGammaLab(W, M);
  const double mr = this->ResonanceMass(kP33_1232);
  const double qr = this->QLab(mr, Q2, M);
  if(kgamma <= 0.0 || qr <= 0.0 || Q2 <= 0.0) return;

  const double bw = fWghtBW ? this->BreitWigner(kP33_1232, W, M) : 1.0;
  const double pref = kPi * (mr/M) * bw / kgamma;

  const double sigma_left_rs =
    0.5 * pref * (ampl.Amp2Plus3() + ampl.Amp2Plus1());
  const double sigma_right_rs =
    0.5 * pref * (ampl.Amp2Minus3() + ampl.Amp2Minus1());
  const double sigma_scalar_rs =
    0.5 * pref * ((qr*qr)/Q2) * ((M*M)/(mr*mr)) *
    (ampl.Amp20Plus() + ampl.Amp20Minus());

  const double bridge = 4.0 * kPi * kAem * mr;
  sigT = bridge * (sigma_left_rs + sigma_right_rs);
  sigL = bridge * (2.0 * sigma_scalar_rs);
}
//____________________________________________________________________________
double MAIDRESPXSec::KGammaLab(double W, double M) const
{
  if(M <= 0.0) return 0.0;
  return (W*W - M*M) / (2.0 * M);
}
//____________________________________________________________________________
double MAIDRESPXSec::QLab(double W, double Q2, double M) const
{
  const double nu = (W*W - M*M + Q2) / (2.0 * M);
  return TMath::Sqrt(TMath::Max(0.0, nu*nu + Q2));
}
//____________________________________________________________________________
double MAIDRESPXSec::QPiCM(double W, double M) const
{
  if(W <= 0.0) return 0.0;
  const double W2 = W*W;
  const double mpi2 = fPionMass*fPionMass;
  const double m2 = M*M;
  const double kallen =
    W2*W2 + m2*m2 + mpi2*mpi2 - 2.0*(W2*m2 + W2*mpi2 + m2*mpi2);
  return TMath::Sqrt(TMath::Max(0.0, kallen)) / (2.0 * W);
}
//____________________________________________________________________________
double MAIDRESPXSec::GammaRunning(Resonance_t res, double W, double M) const
{
  const double mr = this->ResonanceMass(res);
  const double wr = this->ResonanceWidth(res);
  const double xr = this->ResonanceXr(res);
  const double qpi = this->QPiCM(W, M);
  const double qpir = this->QPiCM(mr, M);
  if(W <= 0.0 || mr <= 0.0 || wr <= 0.0 ||
     qpi <= 0.0 || qpir <= 0.0) {
    return 0.0;
  }

  const int L = utils::res::OrbitalAngularMom(res);
  const double power = 2.0 * L + 1.0;
  const double rs_width = wr * TMath::Power(qpi/qpir, power);

  const string width_mode = uppercase_copy(fMAIDRunningWidthMode);
  if(width_mode == "RSORIGINAL" || width_mode == "RS") {
    return rs_width;
  }
  if(width_mode == "MAID2007" || width_mode == "MAID") {
    double barrier = 1.0;
    if(L > 0 && xr > 0.0) {
      barrier =
        TMath::Power((qpir*qpir + xr*xr)/(qpi*qpi + xr*xr), L);
    }
    return rs_width * (mr/W) * barrier;
  }

  LOG("MAIDRESPXSec", pWARN)
    << "Unknown MAIDRunningWidthMode = " << fMAIDRunningWidthMode;
  return rs_width;
}
//____________________________________________________________________________
double MAIDRESPXSec::GammaPi(Resonance_t res, double W, double M) const
{
  const double beta_pi = this->ResonancePiBranching(res);
  const double gamma_running = this->GammaRunning(res, W, M);
  if(beta_pi <= 0.0 || gamma_running <= 0.0) return 0.0;
  return beta_pi * gamma_running;
}
//____________________________________________________________________________
double MAIDRESPXSec::GammaTotal(Resonance_t res, double W, double M) const
{
  return this->GammaRunning(res, W, M);
}
//____________________________________________________________________________
double MAIDRESPXSec::BreitWigner(Resonance_t res, double W, double M) const
{
  double bw = 0.0;
  const string bw_mode = uppercase_copy(fBreitWignerMode);
  if(bw_mode == "MAIDW2") {
    bw = this->BreitWignerMAIDW2(res, W, M);
  }
  else if(bw_mode == "RSORIGINAL" || bw_mode == "RSORIGINALLGAMMA") {
    bw = this->BreitWignerRSOriginal(res, W);
  }
  else {
    LOG("MAIDRESPXSec", pWARN)
      << "Unknown BreitWignerMode = " << fBreitWignerMode;
    return 0.0;
  }

  return bw * this->BreitWignerConventionCorrection(res, W);
}
//____________________________________________________________________________
double MAIDRESPXSec::BreitWignerMAIDW2(
  Resonance_t res, double W, double M) const
{
  const double gamma_tot = this->GammaTotal(res, W, M);
  const double gamma_num =
    fApplyPiNBranchingToMAIDBW ? this->GammaPi(res, W, M) : gamma_tot;
  if(gamma_num <= 0.0 || gamma_tot <= 0.0) return 0.0;
  const double mr = this->ResonanceMass(res);
  const double mr2 = mr*mr;
  const double w2 = W*W;
  return (1.0/kPi) * (mr * gamma_num) /
    (TMath::Power(w2 - mr2, 2) + mr2 * gamma_tot*gamma_tot);
}
//____________________________________________________________________________
double MAIDRESPXSec::BreitWignerRSOriginal(
  Resonance_t res, double W) const
{
  const int LR = utils::res::OrbitalAngularMom(res);
  const double MR = utils::res::Mass(res);
  const double WR = utils::res::Width(res);
  const double NR =
    utils::res::BWNorm(res, fN0ResMaxNWidths, fN2ResMaxNWidths, fGnResMaxNWidths);
  if(W <= 0.0 || MR <= 0.0 || WR <= 0.0 || NR <= 0.0) return 0.0;
  if(utils::res::IsDelta(res)) {
    return utils::bwfunc::BreitWignerLGamma(W, LR, MR, WR, NR);
  }
  return utils::bwfunc::BreitWignerL(W, LR, MR, WR, NR);
}
//____________________________________________________________________________
double MAIDRESPXSec::BreitWignerConventionCorrection(
  Resonance_t res, double W) const
{
  const double mr = this->ResonanceMass(res);
  const string correction = uppercase_copy(fBreitWignerConventionCorrection);
  if(correction == "NONE") return 1.0;
  if(correction == "W2TOLINEAR2MDELTA") return 2.0 * mr;
  if(correction == "LINEARTOW2ONEOVER2MDELTA") return 1.0/(2.0 * mr);
  if(correction == "W2TOLINEAR2W") return 2.0 * W;
  if(correction == "LINEARTOW2ONEOVER2W") return (W > 0.0) ? 1.0/(2.0 * W) : 0.0;
  if(correction == "LINEARTOW2RSNATIVEMATCH") {
    return (W > 0.0 && mr > 0.0) ? W/(2.0 * mr * mr) : 0.0;
  }

  LOG("MAIDRESPXSec", pWARN)
    << "Unknown BreitWignerConventionCorrection = "
    << fBreitWignerConventionCorrection;
  return 1.0;
}
//____________________________________________________________________________
bool MAIDRESPXSec::PassRSNativeWCut(Resonance_t res, double W) const
{
  const int IR = utils::res::ResonanceIndex(res);
  const double MR = utils::res::Mass(res);
  const double WR = utils::res::Width(res);
  if(W <= 0.0 || MR <= 0.0 || WR <= 0.0) return false;
  if(IR == 0) return W <= MR + fN0ResMaxNWidths * WR;
  if(IR == 2) return W <= MR + fN2ResMaxNWidths * WR;
  return W <= MR + fGnResMaxNWidths * WR;
}
//____________________________________________________________________________
bool MAIDRESPXSec::SelectedResonance(Resonance_t res) const
{
  return fResList.Find(res);
}
//____________________________________________________________________________
double MAIDRESPXSec::ResonanceMass(Resonance_t res) const
{
  std::map<Resonance_t, MAID2007ResParam_t>::const_iterator iter =
    fMAID2007ResParams.find(res);
  return (iter != fMAID2007ResParams.end()) ? iter->second.mass :
    utils::res::Mass(res);
}
//____________________________________________________________________________
double MAIDRESPXSec::ResonanceWidth(Resonance_t res) const
{
  std::map<Resonance_t, MAID2007ResParam_t>::const_iterator iter =
    fMAID2007ResParams.find(res);
  return (iter != fMAID2007ResParams.end()) ? iter->second.width :
    utils::res::Width(res);
}
//____________________________________________________________________________
double MAIDRESPXSec::ResonancePiBranching(Resonance_t res) const
{
  std::map<Resonance_t, MAID2007ResParam_t>::const_iterator iter =
    fMAID2007ResParams.find(res);
  return (iter != fMAID2007ResParams.end()) ? iter->second.beta_pi : 0.0;
}
//____________________________________________________________________________
double MAIDRESPXSec::ResonanceXr(Resonance_t res) const
{
  std::map<Resonance_t, MAID2007ResParam_t>::const_iterator iter =
    fMAID2007ResParams.find(res);
  return (iter != fMAID2007ResParams.end()) ? iter->second.xr : 0.5;
}
//____________________________________________________________________________
double MAIDRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this, interaction);
  return xsec;
}
//____________________________________________________________________________
bool MAIDRESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo & proc_info = interaction->ProcInfo();
  const XclsTag & xcls = interaction->ExclTag();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsEM()) return false;
  if(!xcls.KnownResonance()) return false;
  if(!this->SelectedResonance(xcls.Resonance())) return false;

  const int hitnuc = init_state.Tgt().HitNucPdg();
  if(!pdg::IsProton(hitnuc) && !pdg::IsNeutron(hitnuc)) return false;

  const int probe = init_state.ProbePdg();
  if(!pdg::IsChargedLepton(probe)) return false;

  return true;
}
//____________________________________________________________________________
void MAIDRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESPXSec::LoadMAID2007ResParams(void)
{
  fMAID2007ResParams.clear();

  this->LoadMAID2007ResParam(
    kP33_1232, "P33(1232)", maid_res_param(fDeltaMass, fDeltaWidth, 1.00, fDeltaXr));
  this->LoadMAID2007ResParam(
    kP11_1440, "P11(1440)", maid_res_param(1.440, 0.350, 0.70, 0.470));
  this->LoadMAID2007ResParam(
    kD13_1520, "D13(1520)", maid_res_param(1.530, 0.130, 0.60, 0.500));
  this->LoadMAID2007ResParam(
    kS11_1535, "S11(1535)", maid_res_param(1.535, 0.100, 0.40, 0.500));
  this->LoadMAID2007ResParam(
    kS31_1620, "S31(1620)", maid_res_param(1.620, 0.150, 0.25, 0.470));
  this->LoadMAID2007ResParam(
    kS11_1650, "S11(1650)", maid_res_param(1.690, 0.100, 0.85, 0.500));
  this->LoadMAID2007ResParam(
    kD15_1675, "D15(1675)", maid_res_param(1.675, 0.150, 0.45, 0.500));
  this->LoadMAID2007ResParam(
    kF15_1680, "F15(1680)", maid_res_param(1.680, 0.135, 0.70, 0.500));
  this->LoadMAID2007ResParam(
    kD33_1700, "D33(1700)", maid_res_param(1.740, 0.450, 0.15, 0.700));
  this->LoadMAID2007ResParam(
    kP13_1720, "P13(1720)", maid_res_param(1.740, 0.250, 0.20, 0.500));
  this->LoadMAID2007ResParam(
    kF35_1905, "F35(1905)", maid_res_param(1.905, 0.350, 0.10, 0.500));
  this->LoadMAID2007ResParam(
    kP31_1910, "P31(1910)", maid_res_param(1.910, 0.250, 0.25, 0.500));
  this->LoadMAID2007ResParam(
    kF37_1950, "F37(1950)", maid_res_param(1.945, 0.280, 0.40, 0.500));
}
//____________________________________________________________________________
void MAIDRESPXSec::LoadMAID2007ResParam(
  Resonance_t res, const string & res_name, const MAID2007ResParam_t & def)
{
  MAID2007ResParam_t param;
  this->GetParamDef("Mass@" + res_name,   param.mass,    def.mass);
  this->GetParamDef("Width@" + res_name,  param.width,   def.width);
  this->GetParamDef("PiBR@" + res_name,   param.beta_pi, def.beta_pi);
  this->GetParamDef("Xr@" + res_name,     param.xr,      def.xr);
  fMAID2007ResParams[res] = param;
}
//____________________________________________________________________________
void MAIDRESPXSec::LoadConfig(void)
{
  this->GetParamDef("XSecRoute", fXSecRoute, string("DirectAS"));
  this->GetParamDef("RES-EM-XSecScale", fXSecScaleEM, 1.0);
  this->GetParamDef("BreitWignerWeight", fWghtBW, true);
  this->GetParamDef("BreitWignerMode", fBreitWignerMode, string("MAIDW2"));
  this->GetParamDef(
    "BreitWignerConventionCorrection",
    fBreitWignerConventionCorrection,
    string("None"));
  this->GetParamDef(
    "MAIDRunningWidthMode",
    fMAIDRunningWidthMode,
    string("RSOriginal"));
  this->GetParamDef("UseDRJoinScheme", fUseDRJoinScheme, true);
  this->GetParamDef("Wcut", fWcut, 1.7);
  this->GetParamDef(
    "ApplyPiNBranchingToMAIDBW",
    fApplyPiNBranchingToMAIDBW,
    false);
  this->GetParamDef("DeltaMass", fDeltaMass, 1.232);
  this->GetParamDef("DeltaWidth", fDeltaWidth, 0.13);
  this->GetParamDef("DeltaXr", fDeltaXr, 0.57);
  this->GetParamDef("PionMass", fPionMass, kPionMass);
  this->LoadMAID2007ResParams();
  this->GetParamDef(
    "UseRSFallbackForMissingResonances",
    fUseRSFallbackForMissingResonances,
    false);
  this->GetParamDef("ApplyRSNativeWCut", fApplyRSNativeWCut, false);
  this->GetParamDef("MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0);
  this->GetParamDef("MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0);
  this->GetParamDef("MaxNWidthForGNRes", fGnResMaxNWidths, 4.0);

  fResList.Clear();
  string resonances;
  this->GetParamDef(
    "ResonanceNameList", resonances, string(kMAID2007ResonanceNameList));
  fResList.DecodeFromNameList(resonances);

  fHAmplModelEMp =
    dynamic_cast<const MAIDHelicityAmplModelI *>(this->SubAlg("HelicityAmplEMpAlg"));
  fHAmplModelEMn =
    dynamic_cast<const MAIDHelicityAmplModelI *>(this->SubAlg("HelicityAmplEMnAlg"));
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *>(this->SubAlg("XSec-Integrator"));
  fRSFallbackXSecModel = 0;
  if(fUseRSFallbackForMissingResonances) {
    fRSFallbackXSecModel =
      dynamic_cast<const XSecAlgorithmI *>(this->SubAlg("RSFallbackXSecAlg"));
  }

  assert(fHAmplModelEMp);
  assert(fHAmplModelEMn);
  assert(fXSecIntegrator);
  if(fUseRSFallbackForMissingResonances) assert(fRSFallbackXSecModel);
}
//____________________________________________________________________________
