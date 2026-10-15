//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/RunOpt.h"
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

bool check_maid2008_proton_q2(const MAIDHelicityAmplModelEMp & model)
{
  bool ok = true;

  const double d13_q2 = 0.5;
  Interaction * d13 = make_interaction(kD13_1520, false, d13_q2);
  const MAIDHelicityAmpl & d13_ampl = model.Compute(d13);
  const double d13_q4 = d13_q2 * d13_q2;
  const double d13_q8 = d13_q4 * d13_q4;
  const double d13_a12 = -27.4 *
    (1.0 + 8.580*d13_q2 - 0.252*d13_q4 + 0.357*d13_q8) *
    TMath::Exp(-1.20*d13_q2) * 1.0e-3;
  const double d13_a32 = 160.6 *
    (1.0 - 0.820*d13_q2 + 0.541*d13_q4 - 0.016*d13_q8) *
    TMath::Exp(-1.06*d13_q2) * 1.0e-3;
  ok = d13_ampl.IsValid() && ok;
  ok = close(d13_ampl.AmplA12(), d13_a12, 1.0e-12,
             "MAID2008 D13 A12 Q2=0.5") && ok;
  ok = close(d13_ampl.AmplA32(), d13_a32, 1.0e-12,
             "MAID2008 D13 A32 Q2=0.5") && ok;
  delete d13;

  const double p13_q2 = 1.0;
  Interaction * p13 = make_interaction(kP13_1720, false, p13_q2);
  const MAIDHelicityAmpl & p13_ampl = model.Compute(p13);
  const double p13_a32 = -11.5 * (1.0 + 10.83*p13_q2 - 0.66*p13_q2*p13_q2) *
    TMath::Exp(-0.43*p13_q2) * 1.0e-3;
  ok = p13_ampl.IsValid() && ok;
  ok = close(p13_ampl.AmplA32(), p13_a32, 1.0e-12,
             "MAID2008 P13 A32 Q2=1.0") && ok;
  delete p13;

  const double p11_q2 = 0.4;
  Interaction * p11 = make_interaction(kP11_1440, false, p11_q2);
  const MAIDHelicityAmpl & p11_ampl = model.Compute(p11);
  const double p11_q4 = p11_q2 * p11_q2;
  const double p11_q8 = p11_q4 * p11_q4;
  const double p11_a12 = -61.4 *
    (1.0 + 0.871*p11_q2 - 3.516*p11_q4 - 0.158*p11_q8) *
    TMath::Exp(-1.36*p11_q2) * 1.0e-3;
  ok = p11_ampl.IsValid() && ok;
  ok = close(p11_ampl.AmplA12(), p11_a12, 1.0e-12,
             "MAID2008 P11 A12 Q2=0.4") && ok;
  delete p11;

  return ok;
}

bool check_p33_only_width_scale(const MAIDHelicityAmplModelEMp & model)
{
  Interaction * interaction = make_interaction(kP33_1232, false, 0.0);
  const MAIDHelicityAmpl & ampl = model.Compute(interaction);
  const double scale = TMath::Sqrt(130.0/115.0);
  const double a12 = -0.5 * (300.0 + 3.0*(-6.37)) * 1.0e-3 * scale;
  const bool ok = ampl.IsValid() &&
    close(ampl.AmplA12(), a12, 1.0e-12, "P33-only width scale");
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

bool check_maid2008_luis_comparison_config()
{
  const Algorithm * helicity = AlgFactory::Instance()->GetAlgorithm(
    "genie::MAIDHelicityAmplModelEMp", "MAID2008");
  const Algorithm * xsec = AlgFactory::Instance()->GetAlgorithm(
    "genie::MAIDRESPXSec", "DirectAS-MAID2008LuisComparison");
  if(!helicity || !xsec) {
    std::cerr << "Could not load a MAID2008 comparison configuration" << std::endl;
    return false;
  }

  const Registry & hconfig = helicity->GetConfig();
  const Registry & xconfig = xsec->GetConfig();
  bool ok = true;
  ok = hconfig.Exists("MAID2008-A32@P13(1720)-Beta") && ok;
  if(!hconfig.Exists("MAID2008-A32@P13(1720)-Beta") ||
     !close(hconfig.GetDouble("MAID2008-A32@P13(1720)-Beta"), 0.43, 1.0e-12,
            "MAID2008 P13 A32 beta")) ok = false;
  if(!xconfig.Exists("UseDRJoinScheme") || xconfig.GetBool("UseDRJoinScheme")) {
    std::cerr << "Luis comparison configuration must disable the RES/DIS join"
              << std::endl;
    ok = false;
  }
  if(!xconfig.Exists("MAIDTotalWidthMode") ||
     xconfig.GetString("MAIDTotalWidthMode") != "FullMAID") {
    std::cerr << "Luis comparison configuration must select FullMAID widths"
              << std::endl;
    ok = false;
  }
  if(!xconfig.Exists("HelicityAmplEMpAlg") ||
     xconfig.GetAlg("HelicityAmplEMpAlg").config != "MAID2008-P33WidthScaled") {
    std::cerr << "Luis comparison configuration must select MAID2008 P33 scaling"
              << std::endl;
    ok = false;
  }
  return ok;
}

// Exercise factory/XML loading, not just a bare Registry: Default used to
// overwrite the MAID2008 coefficients in the P33-scaled comparison variant.
bool compare_amplitudes(const MAIDHelicityAmplModelEMp & model,
                        const MAIDHelicityAmplModelEMp & reference,
                        bool neutron, const std::string & label)
{
  bool ok = true;
  for(auto res : {kP33_1232, kP11_1440, kD13_1520, kS11_1535,
                  kS31_1620, kS11_1650, kD15_1675, kF15_1680,
                  kD33_1700, kP13_1720, kF35_1905, kP31_1910, kF37_1950}) {
    for(double q2 : {0.0, 0.5, 2.0, 4.0}) {
      std::unique_ptr<Interaction> interaction(make_interaction(res, neutron, q2));
      const MAIDHelicityAmpl got = model.Compute(interaction.get());
      const MAIDHelicityAmpl expected = reference.Compute(interaction.get());
      const std::string point = label + " res=" + std::to_string(int(res)) +
        " Q2=" + std::to_string(q2);
      if(got.IsValid() != expected.IsValid()) {
        std::cerr << point << " validity mismatch" << std::endl;
        ok = false;
      }
      ok = close(got.AmplA12(), expected.AmplA12(), 1.0e-12, point + " A12") && ok;
      ok = close(got.AmplA32(), expected.AmplA32(), 1.0e-12, point + " A32") && ok;
      ok = close(got.AmplS12(), expected.AmplS12(), 1.0e-12, point + " S12") && ok;
    }
  }
  return ok;
}

bool check_factory_amplitudes()
{
  bool ok = true;
  for(bool neutron : {false, true}) {
    const std::string name = neutron ? "genie::MAIDHelicityAmplModelEMn" :
                                       "genie::MAIDHelicityAmplModelEMp";
    for(const std::string config : {"Default", "MAID2008", "MAID2008-P33WidthScaled"}) {
      const auto * model = dynamic_cast<const MAIDHelicityAmplModelEMp *>(
        AlgFactory::Instance()->GetAlgorithm(name, config));
      if(!model) return false;
      std::unique_ptr<MAIDHelicityAmplModelEMp> reference(neutron ?
        new MAIDHelicityAmplModelEMn : new MAIDHelicityAmplModelEMp);
      Registry settings("FactoryAmplitudeReference", false);
      settings.Set("AmplitudeSource", config == "Default" ? "MAID2007" : "MAID2008");
      settings.Set("ApplyLuisP33WidthScale", config == "MAID2008-P33WidthScaled");
      reference->Configure(settings);
      ok = compare_amplitudes(*model, *reference, neutron, name + "/" + config) && ok;
      if(!neutron && config != "Default") {
        ok = check_maid2008_proton_q2(*model) && ok;
      }
    }
    const auto * xsec = AlgFactory::Instance()->GetAlgorithm(
      "genie::MAIDRESPXSec", "DirectAS-MAID2008LuisComparison");
    const auto * actual = dynamic_cast<const MAIDHelicityAmplModelEMp *>(
      xsec->SubAlg(neutron ? "HelicityAmplEMnAlg" : "HelicityAmplEMpAlg"));
    if(!actual) return false;
    if(!neutron) {
      ok = check_maid2008_proton_q2(*actual) && ok;
      ok = check_p33_only_width_scale(*actual) && ok;
    }
  }

  // Selecting a source on an already configured object must replace all
  // source defaults, including polynomial coefficients absent in MAID2007.
  std::unique_ptr<Algorithm> owned(AlgFactory::Instance()->AdoptAlgorithm(
    "genie::MAIDHelicityAmplModelEMp", "Default"));
  auto * model = dynamic_cast<MAIDHelicityAmplModelEMp *>(owned.get());
  if(!model) return false;
  MAIDHelicityAmplModelEMp reference;
  for(const std::string source : {"MAID2008", "MAID2007", "MAID2008"}) {
    Registry settings("SwitchAmplitudeSource", false);
    settings.Set("AmplitudeSource", source);
    model->Configure(settings);
    reference.Configure(settings);
    ok = compare_amplitudes(*model, reference, false, "source switch " + source) && ok;
  }
  // The legacy unqualified keys remain explicit user overrides.
  Registry override("ExplicitCoefficientOverride", false);
  override.Set("A32@P13(1720)-Beta", 0.8);
  model->Configure(override);
  std::unique_ptr<Interaction> interaction(make_interaction(kP13_1720, false, 2.0));
  const double expected = -11.5 * (1.0 + 10.83*2.0 - 0.66*4.0) *
    std::exp(-0.8*2.0) * 1.0e-3;
  ok = close(model->Compute(interaction.get()).AmplA32(), expected, 1.0e-12,
             "explicit beta override on MAID2008") && ok;
  return ok;
}

} // namespace

int main(int /*argc*/, char ** /*argv*/)
{
  RunOpt::Instance()->SetTuneName("MAID07_00a_00_000");
  RunOpt::Instance()->BuildTune();

  MAIDHelicityAmplModelEMp proton_model;
  MAIDHelicityAmplModelEMn neutron_model;
  Registry config("MAID2007HelicityTest", false);
  config.Set("AmplitudeSource", "MAID2007");
  config.Set("ApplyLuisWidthScale", false);
  proton_model.Configure(config);
  neutron_model.Configure(config);

  MAIDHelicityAmplModelEMp maid2008_model;
  Registry maid2008_config("MAID2008HelicityTest", false);
  maid2008_config.Set("AmplitudeSource", "MAID2008");
  maid2008_config.Set("ApplyLuisWidthScale", false);
  maid2008_model.Configure(maid2008_config);

  MAIDHelicityAmplModelEMp p33_scaled_model;
  Registry p33_scaled_config("MAID2008P33ScaledHelicityTest", false);
  p33_scaled_config.Set("AmplitudeSource", "MAID2008");
  p33_scaled_config.Set("ApplyLuisWidthScale", false);
  p33_scaled_config.Set("ApplyLuisP33WidthScale", true);
  p33_scaled_model.Configure(p33_scaled_config);

  bool ok = true;
  ok = check_delta_q2_zero(proton_model) && ok;
  ok = check_roper_q2(proton_model) && ok;
  ok = check_d13_q2(proton_model) && ok;
  ok = check_s11_neutron_q2(neutron_model) && ok;
  ok = check_f15_q2(proton_model) && ok;
  ok = check_maid2008_proton_q2(maid2008_model) && ok;
  ok = check_p33_only_width_scale(p33_scaled_model) && ok;
  ok = check_missing_q2_fit(proton_model) && ok;
  ok = check_maid_integrator_conventions() && ok;
  ok = check_maid2008_luis_comparison_config() && ok;
  ok = check_factory_amplitudes() && ok;

  if(ok) {
    std::cout << "MAID2007/2008 helicity checks passed" << std::endl;
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
  }
  std::fflush(stdout);
  std::fflush(stderr);
  std::_Exit(1);
}

//____________________________________________________________________________
