//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/Registry.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMn.h"
#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h"

using namespace genie;
using namespace genie::constants;

namespace {

bool close(double got, double expected, double tol, const std::string & label)
{
  const double diff = std::fabs(got - expected);
  if(diff <= tol) return true;
  std::cerr << label << " got " << got << " expected " << expected
            << " diff " << diff << " tol " << tol << std::endl;
  return false;
}

Interaction * make_interaction(Resonance_t res, bool neutron, double Q2)
{
  const int target = neutron ? kPdgTgtFreeN : kPdgTgtFreeP;
  const int hitnuc = neutron ? kPdgNeutron : kPdgProton;
  Interaction * interaction =
    Interaction::RESEM(target, hitnuc, kPdgElectron, 4.0);
  interaction->ExclTagPtr()->SetResonance(res);
  interaction->KinePtr()->SetW(1.5);
  interaction->KinePtr()->SetQ2(Q2);
  return interaction;
}

bool check_delta_q2_zero(const MAIDHelicityAmplModelEMp & model)
{
  Interaction * interaction = make_interaction(kP33_1232, false, 0.0);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);

  const double AM = 300.0;
  const double AE = -6.37;
  const double AS = -12.40;
  const double a12 = -0.5 * (AM + 3.0 * AE) * 1.0e-3;
  const double a32 = 0.5 * kSqrt3 * (AE - AM) * 1.0e-3;
  const double s12 = -kSqrt2 * AS * 1.0e-3;

  bool ok = ampl.IsValid();
  ok = close(ampl.AmplA12(), a12, 1.0e-12, "P33 A12 Q2=0") && ok;
  ok = close(ampl.AmplA32(), a32, 1.0e-12, "P33 A32 Q2=0") && ok;
  ok = close(ampl.AmplS12(), s12, 1.0e-12, "P33 S12 Q2=0") && ok;
  delete interaction;
  return ok;
}

bool check_roper_q2(const MAIDHelicityAmplModelEMp & model)
{
  const double Q2 = 0.4;
  Interaction * interaction = make_interaction(kP11_1440, false, Q2);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);

  const double q8 = TMath::Power(Q2, 4);
  const double a12 =
    -61.4 * (1.0 - 1.22 * Q2 - 0.55 * q8) * TMath::Exp(-1.51 * Q2) * 1.0e-3;
  const double s12 =
    4.2 * (1.0 + 40.0 * Q2 + 1.5 * q8) * TMath::Exp(-1.75 * Q2) * 1.0e-3;

  bool ok = ampl.IsValid();
  ok = close(ampl.AmplA12(), a12, 1.0e-12, "P11 A12 Q2=0.4") && ok;
  ok = close(ampl.AmplA32(), 0.0, 1.0e-12, "P11 A32 Q2=0.4") && ok;
  ok = close(ampl.AmplS12(), s12, 1.0e-12, "P11 S12 Q2=0.4") && ok;
  delete interaction;
  return ok;
}

bool check_d13_q2(const MAIDHelicityAmplModelEMp & model)
{
  const double Q2 = 0.5;
  Interaction * interaction = make_interaction(kD13_1520, false, Q2);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);

  const double a12 = -27.0 * (1.0 + 7.77 * Q2) * TMath::Exp(-1.09 * Q2) * 1.0e-3;
  const double a32 = 161.0 * (1.0 + 0.69 * Q2) * TMath::Exp(-2.10 * Q2) * 1.0e-3;
  const double s12 = -63.6 * (1.0 + 4.19 * Q2) * TMath::Exp(-3.40 * Q2) * 1.0e-3;

  bool ok = ampl.IsValid();
  ok = close(ampl.AmplA12(), a12, 1.0e-12, "D13 A12 Q2=0.5") && ok;
  ok = close(ampl.AmplA32(), a32, 1.0e-12, "D13 A32 Q2=0.5") && ok;
  ok = close(ampl.AmplS12(), s12, 1.0e-12, "D13 S12 Q2=0.5") && ok;
  delete interaction;
  return ok;
}

bool check_s11_neutron_q2(const MAIDHelicityAmplModelEMn & model)
{
  const double Q2 = 0.6;
  Interaction * interaction = make_interaction(kS11_1535, true, Q2);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);

  const double a12 = -51.0 * (1.0 + 4.75 * Q2) * TMath::Exp(-1.69 * Q2) * 1.0e-3;
  const double s12 = 28.5 * (1.0 + 0.36 * Q2) * TMath::Exp(-1.55 * Q2) * 1.0e-3;

  bool ok = ampl.IsValid();
  ok = close(ampl.AmplA12(), a12, 1.0e-12, "n S11 A12 Q2=0.6") && ok;
  ok = close(ampl.AmplA32(), 0.0, 1.0e-12, "n S11 A32 Q2=0.6") && ok;
  ok = close(ampl.AmplS12(), s12, 1.0e-12, "n S11 S12 Q2=0.6") && ok;
  delete interaction;
  return ok;
}

bool check_f15_q2(const MAIDHelicityAmplModelEMp & model)
{
  const double Q2 = 0.3;
  Interaction * interaction = make_interaction(kF15_1680, false, Q2);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);

  const double a12 = -25.0 * (1.0 + 3.98 * Q2) * TMath::Exp(-1.20 * Q2) * 1.0e-3;
  const double a32 = 134.0 * (1.0 + 1.00 * Q2) * TMath::Exp(-2.22 * Q2) * 1.0e-3;
  const double s12 = -44.0 * (1.0 + 3.14 * Q2) * TMath::Exp(-1.68 * Q2) * 1.0e-3;

  bool ok = ampl.IsValid();
  ok = close(ampl.AmplA12(), a12, 1.0e-12, "F15 A12 Q2=0.3") && ok;
  ok = close(ampl.AmplA32(), a32, 1.0e-12, "F15 A32 Q2=0.3") && ok;
  ok = close(ampl.AmplS12(), s12, 1.0e-12, "F15 S12 Q2=0.3") && ok;
  delete interaction;
  return ok;
}

bool check_missing_q2_fit(const MAIDHelicityAmplModelEMp & model)
{
  Interaction * f35_q0 = make_interaction(kF35_1905, false, 0.0);
  const bool f35_real_photon_ok = model.Compute(f35_q0).IsValid();
  delete f35_q0;

  Interaction * f35_q2 = make_interaction(kF35_1905, false, 0.1);
  const bool f35_q2_missing = !model.Compute(f35_q2).IsValid();
  delete f35_q2;

  Interaction * p31 = make_interaction(kP31_1910, false, 0.1);
  const bool p31_missing = !model.Compute(p31).IsValid();
  delete p31;

  if(!f35_real_photon_ok) std::cerr << "F35 Q2=0 amplitude should be valid" << std::endl;
  if(!f35_q2_missing) std::cerr << "F35 Q2-dependent amplitude should be missing" << std::endl;
  if(!p31_missing) std::cerr << "P31 MAID amplitude should be missing" << std::endl;
  return f35_real_photon_ok && f35_q2_missing && p31_missing;
}

bool check_maid_integrator_conventions()
{
  const Algorithm * alg =
    AlgFactory::Instance()->GetAlgorithm("genie::MAIDRESPXSec", "Default");
  if(!alg) {
    std::cerr << "Could not load genie::MAIDRESPXSec/Default" << std::endl;
    return false;
  }

  const Registry & config = alg->GetConfig();
  bool ok = true;
  if(!config.Exists("BreitWignerNorm")) {
    std::cerr << "MAIDRESPXSec/Default must explicitly set BreitWignerNorm"
              << std::endl;
    ok = false;
  } else if(config.GetBool("BreitWignerNorm")) {
    std::cerr << "MAIDRESPXSec/Default must not inherit the RES integrator "
              << "BreitWignerNorm W truncation" << std::endl;
    ok = false;
  }

  if(!config.Exists("ApplyRSNativeWCut")) {
    std::cerr << "MAIDRESPXSec/Default must explicitly set ApplyRSNativeWCut"
              << std::endl;
    ok = false;
  } else if(config.GetBool("ApplyRSNativeWCut")) {
    std::cerr << "MAIDRESPXSec/Default should leave RS-native W cuts disabled"
              << std::endl;
    ok = false;
  }

  return ok;
}

} // namespace

int main(int /*argc*/, char ** /*argv*/)
{
  MAIDHelicityAmplModelEMp proton_model;
  MAIDHelicityAmplModelEMn neutron_model;
  Registry config("MAID2007HelicityTest", false);
  config.Set("AmplitudeSource", "MAID2007");
  config.Set("ApplyLuisWidthScale", false);
  proton_model.Configure(config);
  neutron_model.Configure(config);

  bool ok = true;
  ok = check_delta_q2_zero(proton_model) && ok;
  ok = check_roper_q2(proton_model) && ok;
  ok = check_d13_q2(proton_model) && ok;
  ok = check_s11_neutron_q2(neutron_model) && ok;
  ok = check_f15_q2(proton_model) && ok;
  ok = check_missing_q2_fit(proton_model) && ok;
  ok = check_maid_integrator_conventions() && ok;

  if(ok) {
    std::cout << "MAID2007 helicity checks passed" << std::endl;
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
  }
  std::fflush(stdout);
  std::fflush(stderr);
  std::_Exit(1);
}

//____________________________________________________________________________
