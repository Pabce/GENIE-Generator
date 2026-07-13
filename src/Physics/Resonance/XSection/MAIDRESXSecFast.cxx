//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <vector>

#include <TMath.h>
#include <TLorentzVector.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/KPhaseSpace.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Physics/Resonance/XSection/MAIDRESXSecFast.h"

using std::ostringstream;
using std::vector;

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::units;

namespace {

class d2XSecMAIDRESFast_dWQ2_E :
  public ROOT::Math::IBaseFunctionMultiDim {

public:
  d2XSecMAIDRESFast_dWQ2_E(
    const XSecAlgorithmI * model, const Interaction * interaction) :
    ROOT::Math::IBaseFunctionMultiDim(),
    fModel(model),
    fInteraction(interaction),
    fWmin(0.0),
    fWmax(0.0),
    fWcutBelowWmin(false),
    fKPS(0)
  {
    fKPS = fInteraction->PhaseSpacePtr();
    Range1D_t Wl = fKPS->WLim();
    fWmin = Wl.min;
    fWmax = Wl.max;

    Registry config =
      (const_cast<XSecAlgorithmI *>(fModel))->GetConfig();

    const bool use_dr_join =
      config.GetBoolDef("UseDRJoinScheme", true);
    if(use_dr_join) {
      const double wcut = config.GetDoubleDef("Wcut", 1.7);
      fWmax = TMath::Min(wcut, fWmax);
      if(wcut < fWmin) fWcutBelowWmin = true;
    }

    const bool norm_bw =
      config.GetBoolDef("BreitWignerNorm", false);
    if(norm_bw) {
      const double n2_widths =
        config.GetDoubleDef("MaxNWidthForN2Res", 2.0);
      const double n0_widths =
        config.GetDoubleDef("MaxNWidthForN0Res", 6.0);
      const double gn_widths =
        config.GetDoubleDef("MaxNWidthForGNRes", 4.0);
      const Resonance_t resonance = fInteraction->ExclTag().Resonance();
      const int    ir = utils::res::ResonanceIndex(resonance);
      const double mr = utils::res::Mass(resonance);
      const double wr = utils::res::Width(resonance);
      double wcut = mr + gn_widths * wr;
      if(ir == 0) wcut = mr + n0_widths * wr;
      else if(ir == 2) wcut = mr + n2_widths * wr;
      fWmax = TMath::Min(wcut, fWmax);
      if(wcut < fWmin) fWcutBelowWmin = true;
    }

    if(fWmax < fWmin) fWcutBelowWmin = true;
  }

  ~d2XSecMAIDRESFast_dWQ2_E() {}

  unsigned int NDim(void) const { return 2; }

  double DoEval(const double * xin) const
  {
    if(fWcutBelowWmin) return 0.0;

    const double dW = fWmax - fWmin;
    if(dW <= 0.0) return 0.0;

    const double W = fWmin + dW * xin[0];
    fInteraction->KinePtr()->SetW(W);

    Range1D_t Q2l = fKPS->Q2Lim_W();
    if(Q2l.min < 0.0 || Q2l.max < 0.0 || Q2l.max < Q2l.min) {
      return 0.0;
    }

    const double dQ2 = Q2l.max - Q2l.min;
    const double Q2  = Q2l.min + dQ2 * xin[1];
    fInteraction->KinePtr()->SetQ2(Q2);

    const double xsec =
      fModel->XSec(fInteraction, kPSWQ2fE) * dW * dQ2;
    return xsec / (1E-38 * units::cm2);
  }

  ROOT::Math::IBaseFunctionMultiDim * Clone(void) const
  {
    return new d2XSecMAIDRESFast_dWQ2_E(fModel, fInteraction);
  }

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double fWmin;
  double fWmax;
  bool fWcutBelowWmin;
  KPhaseSpace * fKPS;
};

}

//____________________________________________________________________________
MAIDRESXSecFast::MAIDRESXSecFast() :
XSecIntegratorI("genie::MAIDRESXSecFast"),
fUsePauliBlocking(false),
fUseCache(false),
fCacheAllResonances(false),
fEMin(0.01),
fEMax(20.0),
fNKnots(50),
fSingleResXSecModel(0)
{

}
//____________________________________________________________________________
MAIDRESXSecFast::MAIDRESXSecFast(string config) :
XSecIntegratorI("genie::MAIDRESXSecFast", config),
fUsePauliBlocking(false),
fUseCache(false),
fCacheAllResonances(false),
fEMin(0.01),
fEMax(20.0),
fNKnots(50),
fSingleResXSecModel(0)
{

}
//____________________________________________________________________________
MAIDRESXSecFast::~MAIDRESXSecFast()
{

}
//____________________________________________________________________________
double MAIDRESXSecFast::Integrate(
  const XSecAlgorithmI * model, const Interaction * interaction) const
{
  if(!model->ValidProcess(interaction)) return 0.0;
  fSingleResXSecModel = model;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
    LOG("MAIDRESXSecFast", pDEBUG) << "*** Below energy threshold";
    return 0.0;
  }

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo & proc_info = interaction->ProcInfo();
  const Target & target = init_state.Tgt();

  const InteractionType_t it = proc_info.InteractionTypeId();
  const int nucleon_pdgc = target.HitNucPdg();
  const int probe_pdgc = init_state.ProbePdg();
  const double Ev = init_state.ProbeE(kRfHitNucRest);
  const Resonance_t res = interaction->ExclTag().Resonance();

  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty()) {
    Interaction * in = new Interaction(*interaction);
    if(pdg::IsProton(nucleon_pdgc)) {
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeP);
    } else {
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeN);
    }
    if(xsl->SplineExists(model, in)) {
      const Spline * spl = xsl->GetSpline(model, in);
      double xsec = spl->Evaluate(Ev);
      SLOG("MAIDRESXSecFast", pNOTICE)
        << "XSec[RES/" << utils::res::AsString(res) << "/free] (Ev = "
        << Ev << " GeV) = " << xsec/(1E-38 * cm2) << " x 1E-38 cm^2";
      if(!interaction->TestBit(kIAssumeFreeNucleon)) {
        const int nnucl =
          pdg::IsProton(nucleon_pdgc) ? target.Z() : target.N();
        xsec *= nnucl;
      }
      delete in;
      return xsec;
    }
    delete in;
  }

  const bool bare_xsec_pre_calc = RunOpt::Instance()->BareXSecPreCalc();
  if(fUseCache && bare_xsec_pre_calc && !fUsePauliBlocking) {
    Cache * cache = Cache::Instance();
    string key = this->CacheBranchName(res, it, probe_pdgc, nucleon_pdgc);
    LOG("MAIDRESXSecFast", pINFO)
      << "Finding cache branch with key: " << key;
    CacheBranchFx * cache_branch =
      dynamic_cast<CacheBranchFx *>(cache->FindCacheBranch(key));

    if(!cache_branch) {
      LOG("MAIDRESXSecFast", pNOTICE)
        << "Computing/caching MAID RES production xsec";
      if(fCacheAllResonances) {
        this->CacheSelectedResExcitationXSec(interaction);
      } else {
        this->CacheResExcitationXSec(interaction, res);
      }
      cache_branch =
        dynamic_cast<CacheBranchFx *>(cache->FindCacheBranch(key));
      assert(cache_branch);
    }

    const CacheBranchFx & cbranch = (*cache_branch);
    const double eps = 1E-12;
    const double eval_E =
      TMath::Max(fEMin, TMath::Min(Ev, fEMax * (1.0 - eps)));
    double rxsec = cbranch(eval_E);

    SLOG("MAIDRESXSecFast", pNOTICE)
      << "XSec[RES/" << utils::res::AsString(res) << "/free] (Ev = "
      << Ev << " GeV) = " << rxsec/(1E-38 * cm2) << " x 1E-38 cm^2";

    if(interaction->TestBit(kIAssumeFreeNucleon)) return rxsec;

    const int nnucl = pdg::IsProton(nucleon_pdgc) ? target.Z() : target.N();
    rxsec *= nnucl;
    return rxsec;
  }

  return this->IntegrateDirect(model, interaction);
}
//____________________________________________________________________________
double MAIDRESXSecFast::IntegrateDirect(
  const XSecAlgorithmI * model, const Interaction * interaction) const
{
  const Resonance_t res = interaction->ExclTag().Resonance();
  const double Ev = interaction->InitState().ProbeE(kRfHitNucRest);

  LOG("MAIDRESXSecFast", pINFO)
    << "*** Integrating d^2 XSec/dWdQ^2 for R: "
    << utils::res::AsString(res) << " at Ev = " << Ev;

  ROOT::Math::IBaseFunctionMultiDim * func =
    new d2XSecMAIDRESFast_dWQ2_E(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  ROOT::Math::IntegratorMultiDim ig(
    ig_type, 0, fGSLRelTol, fGSLMaxEval);
  ig.SetFunction(*func);

  double kine_min[2] = { 0.0, 0.0 };
  double kine_max[2] = { 1.0, 1.0 };
  double xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);

  delete func;
  return xsec;
}
//____________________________________________________________________________
void MAIDRESXSecFast::CacheSelectedResExcitationXSec(
  const Interaction * interaction) const
{
  const unsigned int nres = fResList.NResonances();
  for(unsigned int ires = 0; ires < nres; ires++) {
    this->CacheResExcitationXSec(interaction, fResList.ResonanceId(ires));
  }
}
//____________________________________________________________________________
void MAIDRESXSecFast::CacheResExcitationXSec(
  const Interaction * in, Resonance_t res) const
{
  assert(fSingleResXSecModel);

  Cache * cache = Cache::Instance();

  const int probe_code = in->InitState().ProbePdg();
  const int nuc_code = in->InitState().Tgt().HitNucPdg();
  const int tgt_code =
    pdg::IsProton(nuc_code) ? kPdgTgtFreeP : kPdgTgtFreeN;

  Interaction * interaction = new Interaction(*in);
  interaction->InitStatePtr()->SetPdgs(tgt_code, probe_code);
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(nuc_code);
  interaction->ExclTagPtr()->SetResonance(res);

  const InteractionType_t it =
    interaction->ProcInfo().InteractionTypeId();
  const string key = this->CacheBranchName(res, it, probe_code, nuc_code);

  CacheBranchFx * cache_branch =
    dynamic_cast<CacheBranchFx *>(cache->FindCacheBranch(key));
  if(cache_branch) {
    delete interaction;
    return;
  }

  LOG("MAIDRESXSecFast", pNOTICE)
    << "\n ** Creating cache branch - key = " << key;
  cache_branch = new CacheBranchFx("MAID RES Excitation XSec");
  cache->AddCacheBranch(key, cache_branch);
  assert(cache_branch);

  const KPhaseSpace & kps = interaction->PhaseSpace();
  const double Ethr = kps.Threshold();
  LOG("MAIDRESXSecFast", pNOTICE) << "E threshold = " << Ethr;

  const int nknots = TMath::Max(2, fNKnots);
  vector<double> E(nknots, fEMin);

  const bool threshold_in_range = (Ethr > fEMin && Ethr < fEMax);
  int nkb = (threshold_in_range && nknots > 7) ? 5 : 0;
  int nka = nknots - nkb;
  if(nka < 2) {
    nkb = 0;
    nka = nknots;
  }

  if(fEMax <= fEMin) {
    E[0] = fEMin;
    for(int i = 1; i < nknots; i++) E[i] = fEMin;
  }
  else if(!threshold_in_range && fEMax <= TMath::Max(Ethr, fEMin)) {
    const double dE = (fEMax - fEMin) / (nknots - 1);
    for(int i = 0; i < nknots; i++) E[i] = fEMin + i*dE;
  }
  else {
    if(nkb > 0) {
      const double dEb = (Ethr - fEMin) / nkb;
      for(int i = 0; i < nkb; i++) E[i] = fEMin + i*dEb;
    }

    const double E0 = TMath::Max(Ethr, fEMin);
    const double logE0 = TMath::Log10(E0);
    const double dEa = (TMath::Log10(fEMax) - logE0) / (nka - 1);
    for(int i = 0; i < nka; i++) {
      E[i+nkb] = TMath::Power(10.0, logE0 + i*dEa);
    }
  }

  TLorentzVector p4(0.0, 0.0, 0.0, 0.0);
  for(int ie = 0; ie < nknots; ie++) {
    double xsec = 0.0;
    const double Ev = E[ie];
    p4.SetPxPyPzE(0.0, 0.0, Ev, Ev);
    interaction->InitStatePtr()->SetProbeP4(p4);

    if(Ev > Ethr + kASmallNum) {
      xsec = this->IntegrateDirect(fSingleResXSecModel, interaction);
    } else {
      LOG("MAIDRESXSecFast", pINFO)
        << "** Below threshold E = " << Ev << " <= " << Ethr;
    }

    cache_branch->AddValues(Ev, xsec);
    SLOG("MAIDRESXSecFast", pNOTICE)
      << "RES XSec (R:" << utils::res::AsString(res)
      << ", E=" << Ev << ") = "
      << xsec/(1E-38 * units::cm2) << " x 1E-38 cm^2";
  }

  cache_branch->CreateSpline();
  delete interaction;
}
//____________________________________________________________________________
string MAIDRESXSecFast::CacheBranchName(
  Resonance_t res, InteractionType_t it, int probe_pdgc,
  int nucleon_pdgc) const
{
  Cache * cache = Cache::Instance();

  const string res_name = utils::res::AsString(res);
  const string it_name = InteractionType::AsString(it);
  const string nuc_name = pdg::IsProton(nucleon_pdgc) ? "p" : "n";

  Registry config =
    (const_cast<XSecAlgorithmI *>(fSingleResXSecModel))->GetConfig();
  const RgAlg hemp =
    config.GetAlgDef("HelicityAmplEMpAlg", RgAlg("", ""));
  const RgAlg hemn =
    config.GetAlgDef("HelicityAmplEMnAlg", RgAlg("", ""));

  ostringstream intk;
  intk << "MAIDResExcitationXSec/R:" << res_name
       << ";probe:" << probe_pdgc
       << ";int:" << it_name << nuc_name
       << ";route:" << config.GetStringDef("XSecRoute", "DirectAS")
       << ";hamp:" << hemp.name << "/" << hemp.config
       << ";hamn:" << hemn.name << "/" << hemn.config
       << ";bwght:" << config.GetBoolDef("BreitWignerWeight", true)
       << ";bwmode:" << config.GetStringDef("BreitWignerMode", "MAIDW2")
       << ";bwcorr:" << config.GetStringDef(
         "BreitWignerConventionCorrection", "None")
       << ";width:" << config.GetStringDef(
         "MAIDRunningWidthMode", "RSOriginal")
       << ";pin:" << config.GetBoolDef(
         "ApplyPiNBranchingToMAIDBW", false)
       << ";dr:" << config.GetBoolDef("UseDRJoinScheme", true)
       << ";wcut:" << config.GetDoubleDef("Wcut", 1.7)
       << ";rscut:" << config.GetBoolDef("ApplyRSNativeWCut", false)
       << ";bwnorm:" << config.GetBoolDef("BreitWignerNorm", false)
       << ";rsfb:" << config.GetBoolDef(
         "UseRSFallbackForMissingResonances", false)
       << ";mass:" << config.GetDoubleDef(
         "Mass@" + res_name, utils::res::Mass(res))
       << ";widthr:" << config.GetDoubleDef(
         "Width@" + res_name, utils::res::Width(res))
       << ";pibr:" << config.GetDoubleDef("PiBR@" + res_name, 0.0)
       << ";xr:" << config.GetDoubleDef("Xr@" + res_name, 0.5);

  const string algkey = fSingleResXSecModel->Id().Key();
  return cache->CacheBranchKey(algkey, intk.str());
}
//____________________________________________________________________________
void MAIDRESXSecFast::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESXSecFast::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESXSecFast::LoadConfig(void)
{
  GetParamDef("gsl-integration-type", fGSLIntgType, string("adaptive"));
  GetParamDef("gsl-relative-tolerance", fGSLRelTol, 0.01);
  GetParamDef("gsl-max-eval", fGSLMaxEval, 100000);
  GetParamDef("UsePauliBlockingForRES", fUsePauliBlocking, false);
  GetParamDef("UseCache", fUseCache, false);
  GetParamDef("CacheAllResonances", fCacheAllResonances, false);
  GetParamDef("CacheEMin", fEMin, 0.01);
  GetParamDef("ESplineMax", fEMax, 20.0);
  GetParamDef("CacheNKnots", fNKnots, 50);

  fEMin = TMath::Max(fEMin, 1E-6);
  fEMax = TMath::Max(fEMax, fEMin * (1.0 + 1E-6));
  fNKnots = TMath::Max(fNKnots, 2);

  fResList.Clear();
  string resonances;
  GetParamDef("ResonanceNameList", resonances, string(""));
  fResList.DecodeFromNameList(resonances);
}
//____________________________________________________________________________
