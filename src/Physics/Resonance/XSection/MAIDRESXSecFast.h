//____________________________________________________________________________
/*!

\class    genie::MAIDRESXSecFast

\brief    Fast W,Q2 cross-section integrator for MAID electromagnetic
          resonance production.

          This is the MAID analogue of the fast RES integrator: it transforms
          the W,Q2 integration area and can cache free-nucleon resonance
          excitation cross sections. The cache is intentionally MAID-specific
          so that MAID configuration choices, resonance lists and W-support
          conventions are controlled by the MAID model rather than by the
          Rein-Sehgal cache helper.

\author   GENIE Collaboration

\created  July 07, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _MAID_RES_XSEC_FAST_H_
#define _MAID_RES_XSEC_FAST_H_

#include "Framework/ParticleData/BaryonResList.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/Interaction/InteractionType.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class MAIDRESXSecFast : public XSecIntegratorI {

public:
  MAIDRESXSecFast();
  MAIDRESXSecFast(string param_set);
  virtual ~MAIDRESXSecFast();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig(void);
  double IntegrateDirect(const XSecAlgorithmI * model,
                         const Interaction * i) const;
  void   CacheSelectedResExcitationXSec(const Interaction * i) const;
  void   CacheResExcitationXSec(const Interaction * i, Resonance_t r) const;
  string CacheBranchName(Resonance_t r, InteractionType_t it,
                         int probe, int nuc) const;

  bool   fUsePauliBlocking;      ///< account for Pauli blocking?
  bool   fUseCache;              ///< use local free-nucleon xsec cache?
  bool   fCacheAllResonances;    ///< warm all selected resonances at once?
  double fEMin;                  ///< minimum cache-spline energy
  double fEMax;                  ///< maximum cache-spline energy
  int    fNKnots;                ///< number of cache-spline knots

  mutable const XSecAlgorithmI * fSingleResXSecModel;
  BaryonResList fResList;
};

}       // genie namespace

#endif  // _MAID_RES_XSEC_FAST_H_
