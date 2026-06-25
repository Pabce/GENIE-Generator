//____________________________________________________________________________
/*!

\program gmkspl

\brief   GENIE utility program building XML cross section splines that can
         be loaded into GENIE to speed-up event generation.
         The list of neutrino PDG codes is passed from the command line.
         The list of nuclear target PDG codes is either passed from the
         command line or extracted from the input ROOT/GEANT geometry.

         Syntax :

           gmkspl  -p nupdg
                  <-t tgtpdg, -f geomfile>
                  <-o | --output-cross-sections> xsec_xml_file_name
                  [-n nknots]
                  [-e max_energy]
                  [--no-copy]
                  [--seed seed_number]
                  [--input-cross-sections xml_file]
                  [--fast-spline-eval]
                  [--fast-spline-knots exact_knots]
                  [--fast-res-spline-knots exact_knots]
                  [--fast-dis-spline-knots exact_knots]
                  [--fast-spline-checks guard_knots]
                  [--fast-res-spline-checks guard_knots]
                  [--fast-dis-spline-checks guard_knots]
                  [--fast-spline-tol rel_tolerance]
                  [--fast-spline-area-tol rel_tolerance]
                  [--fast-coh-spline-tol rel_tolerance]
                  [--fast-res-spline-tol rel_tolerance]
                  [--fast-qel-spline-tol rel_tolerance]
                  [--fast-dis-spline-tol rel_tolerance]
                  [--fast-dis-spline-area-tol rel_tolerance]
                  [--fast-spline-refine-rounds rounds]
                  [--fast-res-spline-refine-rounds rounds]
                  [--fast-res-spline-refine-tol rel_tolerance]
                  [--fast-qel-spline-refine-rounds rounds]
                  [--fast-dis-spline-refine-rounds rounds]
                  [--fast-dis-spline-refine-tol rel_tolerance]
                  [--fast-spline-knot-jobs worker_processes]
                  [--fast-spline-quiet-knot-workers]
                  [--fast-qel-throws nucleon_throws]
                  [--fast-qel-low-e-throws nucleon_throws]
                  [--fast-qel-low-e-max energy]
                  [--fast-qel-rel-tol tolerance]
                  [--fast-qel-max-eval evals]
                  [--fast-qel-skip-kin-check]
                  [--fast-coh-t-steps t_steps]
                  [--fast-coh-analytic-t]
                  [--q2-min q2_value]
                  [--jobs worker_processes]

                  // command line args handled by RunOpt:
                  [--event-generator-list list_name] // default "Default"
                  [--tune tune_name]  // default "G18_02a_00_000"
                  [--xml-path path]
                  [--message-thresholds xml_file]


         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.

         Options :
           -p
               A comma separated list of nu PDG codes.
           -t
               A comma separated list of tgt PDG codes.
               PDG code format: 10LZZZAAAI
           -f
               A ROOT file containing a ROOT/GEANT geometry description.
           -o, --output-cross-sections
               Name of output XML file containing computed cross-section data.
               Default: `xsec_splines.xml'.
           -n
               Number of knots per spline.
               Default: 15 knots per decade of energy range with a minimum
               of 30 knots totally.
           -e
               Maximum energy in spline.
               Default: The max energy in the validity range of the spline
               generating thread.
           --no-copy
               Does not write out the input cross-sections in the output file
           --seed
              Random number seed.
           --input-cross-sections
              Name (incl. full path) of an XML file with pre-computed
              free-nucleon cross-section values. If loaded, it can speed-up
              cross-section calculation for nuclear targets.
           --fast-spline-eval
              Skip explicit integrations at spline knots below the interaction
              threshold and summarize per-spline timing instead of printing
              every knot evaluation.
           --fast-spline-knots
              Aggressive opt-in sparse mode: evaluate this many above-threshold
              knots exactly, interpolate the rest onto the normal output grid,
              and emit regular GENIE spline XML. Disabled by default.
           --fast-res-spline-knots
              RES-only exact-knot override for sparse mode. This can densify
              resonance turn-ons without changing the global sparse grid.
           --fast-dis-spline-knots
              DIS-only exact-knot override for sparse mode. This can let DIS
              use a denser or looser sparse grid without changing other rows.
           --fast-spline-checks
              In sparse mode, evaluate this many interpolated guard knots
              exactly. Passing guards are promoted into interpolation anchors.
              If the worst guard relative difference exceeds --fast-spline-tol,
              rebuild that spline exactly. Default: 0.
           --fast-res-spline-checks
              RES-only sparse guard-knot override. This lets resonance splines
              use a denser guard set without changing the global guard count.
           --fast-dis-spline-checks
              DIS-only sparse guard-knot override. This lets DIS rows use a
              denser guard set without changing the global guard count.
           --fast-spline-tol
              Sparse-mode guard relative-difference tolerance before exact
              fallback. Default: 0.25.
           --fast-spline-area-tol
              Sparse-mode guard area-difference tolerance before exact
              fallback. Disabled by default.
           --fast-coh-spline-tol
              COH-only sparse-mode guard tolerance before exact fallback.
              This leaves the global sparse fallback tolerance unchanged for
              all non-coherent splines.
           --fast-res-spline-tol
              RES-only sparse-mode guard tolerance before exact fallback.
              This leaves the global sparse fallback tolerance unchanged for
              all non-resonant splines.
           --fast-qel-spline-tol
              QEL-only sparse-mode guard tolerance before exact fallback.
              This lets expensive QEL rows use a looser sparse acceptance
              without changing DIS, RES, COH, or NuElectron policy.
           --fast-dis-spline-tol
              DIS-only sparse-mode guard tolerance before exact fallback.
           --fast-dis-spline-area-tol
              DIS-only sparse-mode guard area-difference tolerance before exact
              fallback. Useful for keeping aggressive DIS sparse rows honest
              in total area under the curve.
           --fast-spline-refine-rounds
              In sparse mode, promote failing guard knots and repeat guard
              selection for this many extra rounds before exact fallback.
           --fast-res-spline-refine-rounds
              RES-only sparse guard-promotion rounds. This can tighten
              resonance turn-ons without changing the global sparse policy.
           --fast-res-spline-refine-tol
              RES-only guard relative-difference tolerance that triggers
              refinement. Fallback still uses --fast-res-spline-tol.
              Default: 0.
           --fast-qel-spline-refine-rounds
              QEL-only sparse refinement override. This can avoid expensive
              QEL full-grid fallbacks without changing non-QEL sparse policy.
           --fast-dis-spline-refine-rounds
              DIS-only sparse guard-promotion rounds.
           --fast-dis-spline-refine-tol
              DIS-only guard relative-difference tolerance that triggers
              refinement. Fallback still uses --fast-dis-spline-tol.
              Default: 0.
           --fast-spline-knot-jobs
              Opt-in worker processes for exact knot integrations inside a
              single spline. This can accelerate one expensive interaction
              that --jobs cannot split. Default: 1.
           --fast-spline-quiet-knot-workers
              Redirect stdout/stderr to /dev/null inside exact-knot worker
              children. Parent process summaries and worker failure messages
              remain visible. Useful when forked worker notices interleave in
              production logs.
           --fast-qel-throws
              Override NewQELXSec/Default NumNucleonThrows for this process.
              Useful for aggressive gmkspl scans of nuclear QE splines.
           --fast-qel-low-e-throws
              Use this many NewQELXSec nucleon throws only for probe energies
              at or below --fast-qel-low-e-max. This can protect the threshold
              bend while keeping high-energy knots cheap.
           --fast-qel-low-e-max
              Probe-energy cutoff for --fast-qel-low-e-throws in GeV.
              Default inside NewQELXSec: 0.3.
           --fast-qel-rel-tol
              Override NewQELXSec/Default gsl-relative-tolerance for this
              process.
           --fast-qel-max-eval
              Override NewQELXSec/Default gsl-max-eval for this process.
           --fast-qel-skip-kin-check
              Ask NewQELXSec to skip the model-level ValidKinematics call.
              ComputeFullQELPXSec still rejects points outside the Q2 range
              before the model XSec call.
           --fast-coh-t-steps
              Override the BergerSehgalCOHPiPXSec2015 inner t quadrature step
              count for this process. Default physics configuration uses 100.
           --fast-coh-analytic-t
              Replace BergerSehgalCOHPiPXSec2015's inner t quadrature with the
              closed-form exponential integral for gmkspl speed scans.
           --q2-min
              Override the minimum Q^2 phase-space cut (in GeV^2). By default,
              the EM cut is configured in CommonPhaseSpaceCuts.xml [Default]
              as EM-Q2-min. For weak interactions, this explicitly opts the
              run into the Q2 cut.
           --jobs
              Build independent probe-target initial states in parallel worker
              processes. For a command with only one probe-target initial
              state, use the same worker count to build missing splines for
              that state in parallel. Default: 1, or GENIE_GMKSPL_JOBS if set.

           --event-generator-list
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "Default"].
           --xml-path
              A directory to load XML files from - overrides $GXMLPATH, and $GENIE/config
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.

        ***  See the User Manual for more details and examples. ***

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created September 27, 2005

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

#include <TSystem.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
//#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/CmdLnArgParser.h"

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#include "Tools/Geometry/ROOTGeomAnalyzer.h"
#endif

using std::string;
using std::vector;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::steady_clock;

using namespace genie;

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
using namespace genie::geometry;
#endif

// Prototypes:
void          GetCommandLineArgs (int argc, char ** argv);
void          PrintSyntax        (void);
PDGCodeList * GetNeutrinoCodes   (void);
PDGCodeList * GetTargetCodes     (void);
vector<InitialState> GetInitialStates(PDGCodeList * neutrinos,
                                      PDGCodeList * targets);
void          BuildSplinesForInitialState(const InitialState & init_state);
void          BuildSplinesSerial(const vector<InitialState> & init_states);
void          BuildSplinesParallel(const vector<InitialState> & init_states);
string        ParallelWorkerOutputName(size_t job_index);
double        SecondsSince(const steady_clock::time_point & start);
void          ApplyFastQELOverrides(void);
void          ApplyFastCOHOverrides(void);

// User-specified options:
string   gOptNuPdgCodeList  = "";
string   gOptTgtPdgCodeList = "";
string   gOptGeomFilename   = "";
int      gOptNKnots         = -1;
double   gOptMaxE           = -1.;
bool     gOptNoCopy         = false;
long int gOptRanSeed        = -1;   // random number seed
string   gOptInpXSecFile    = "";   // input cross-section file
string   gOptOutXSecFile    = "";   // output cross-section file
double   gOptQ2Min;                 // Q2 minimum override value
bool     gOptQ2MinSet       = false; // whether --q2-min was specified
bool     gOptFastSplineEval = false; // whether --fast-spline-eval was specified
int      gOptFastSplineKnots = 0;    // above-threshold exact knots in sparse mode
int      gOptFastRESSplineKnots = 0; // RES-only above-threshold exact knots
int      gOptFastDISSplineKnots = 0; // DIS-only above-threshold exact knots
int      gOptFastSplineChecks = 0;   // sparse-mode exact guard knots
int      gOptFastRESSplineChecks = 0; // RES-only sparse exact guard knots
int      gOptFastDISSplineChecks = 0; // DIS-only sparse exact guard knots
double   gOptFastSplineTol = 0.25;   // sparse-mode exact fallback tolerance
double   gOptFastSplineAreaTol = 0.; // sparse-mode area fallback tolerance
double   gOptFastCOHSplineTol = 0.;  // COH-only sparse fallback tolerance
double   gOptFastRESSplineTol = 0.;  // RES-only sparse fallback tolerance
double   gOptFastQELSplineTol = 0.;  // QEL-only sparse fallback tolerance
double   gOptFastDISSplineTol = 0.;  // DIS-only sparse fallback tolerance
double   gOptFastDISSplineAreaTol = 0.; // DIS-only area fallback tolerance
int      gOptFastSplineRefineRounds = 0; // extra guard-promotion rounds
int      gOptFastRESSplineRefineRounds = 0; // RES-only guard-promotion rounds
double   gOptFastRESSplineRefineTol = 0.; // RES-only refinement trigger
int      gOptFastQELSplineRefineRounds = 0; // QEL-only guard-promotion rounds
int      gOptFastDISSplineRefineRounds = 0; // DIS-only guard-promotion rounds
double   gOptFastDISSplineRefineTol = 0.; // DIS-only refinement trigger
int      gOptFastSplineKnotJobs = 1; // worker processes inside one spline
bool     gOptFastSplineQuietKnotWorkers = false; // silence exact-knot workers
int      gOptFastQELThrows = 0;      // optional NewQELXSec nucleon throws override
int      gOptFastQELLowEThrows = 0;  // optional low-energy NewQELXSec throws override
double   gOptFastQELLowEMax = 0.;    // optional low-energy throw cutoff
double   gOptFastQELRelTol = 0.;     // optional NewQELXSec relative tolerance override
int      gOptFastQELMaxEval = 0;     // optional NewQELXSec max eval override
bool     gOptFastQELSkipKinCheck = false; // skip redundant QEL kinematics check
int      gOptFastCOHTSteps = 0;      // optional COH t quadrature step override
bool     gOptFastCOHAnalyticT = false; // use analytic COH t integration
int      gOptJobs           = 1;     // number of gmkspl worker processes

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);
  ApplyFastQELOverrides();
  ApplyFastCOHOverrides();

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Apply Q2 minimum override before physics code queries phase-space cuts.
  if(gOptQ2MinSet) {
    KPhaseSpaceCuts::Instance()->SetQ2MinOverride(gOptQ2Min);
    LOG("gmkspl", pNOTICE)
      << "Overriding Q2 minimum phase-space cut from command line: "
      << gOptQ2Min << " GeV^2";
  }

  // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  // Init
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Get list of neutrinos and nuclear targets

  steady_clock::time_point list_start = steady_clock::now();
  LOG("gmkspl", pNOTICE) << "Building neutrino and target code lists";
  PDGCodeList * neutrinos = GetNeutrinoCodes();
  PDGCodeList * targets   = GetTargetCodes();
  LOG("gmkspl", pNOTICE)
    << "Built neutrino and target code lists in "
    << SecondsSince(list_start) << " s";

  if(!neutrinos || neutrinos->size() == 0 ) {
     LOG("gmkspl", pFATAL) << "Empty neutrino PDG code list";
     PrintSyntax();
     exit(2);
  }
  if(!targets || targets->size() == 0 ) {
     LOG("gmkspl", pFATAL) << "Empty target PDG code list";
     PrintSyntax();
     exit(3);
  }

  LOG("gmkspl", pINFO) << "Neutrinos: " << *neutrinos;
  LOG("gmkspl", pINFO) << "Targets: "   << *targets;

  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->SetCurrentTune(RunOpt::Instance()->Tune()->Name());
  xspl->SetFastSplineEval(gOptFastSplineEval || gOptFastSplineKnots > 0);
  xspl->SetFastSplineKnots(gOptFastSplineKnots);
  xspl->SetFastRESSplineKnots(gOptFastRESSplineKnots);
  xspl->SetFastDISSplineKnots(gOptFastDISSplineKnots);
  xspl->SetFastSplineChecks(gOptFastSplineChecks);
  xspl->SetFastRESSplineChecks(gOptFastRESSplineChecks);
  xspl->SetFastDISSplineChecks(gOptFastDISSplineChecks);
  xspl->SetFastSplineTolerance(gOptFastSplineTol);
  xspl->SetFastSplineAreaTolerance(gOptFastSplineAreaTol);
  xspl->SetFastCOHSplineTolerance(gOptFastCOHSplineTol);
  xspl->SetFastRESSplineTolerance(gOptFastRESSplineTol);
  xspl->SetFastQELSplineTolerance(gOptFastQELSplineTol);
  xspl->SetFastDISSplineTolerance(gOptFastDISSplineTol);
  xspl->SetFastDISSplineAreaTolerance(gOptFastDISSplineAreaTol);
  xspl->SetFastSplineRefineRounds(gOptFastSplineRefineRounds);
  xspl->SetFastRESSplineRefineRounds(gOptFastRESSplineRefineRounds);
  xspl->SetFastRESSplineRefineTolerance(gOptFastRESSplineRefineTol);
  xspl->SetFastQELSplineRefineRounds(gOptFastQELSplineRefineRounds);
  xspl->SetFastDISSplineRefineRounds(gOptFastDISSplineRefineRounds);
  xspl->SetFastDISSplineRefineTolerance(gOptFastDISSplineRefineTol);
  xspl->SetFastSplineKnotJobs(gOptFastSplineKnotJobs);
  xspl->SetFastSplineQuietKnotWorkers(gOptFastSplineQuietKnotWorkers);

  vector<InitialState> init_states = GetInitialStates(neutrinos, targets);
  int spline_jobs = 1;
  if(gOptJobs > 1 && init_states.size() == 1) {
    spline_jobs = gOptJobs;
    std::ostringstream value;
    value << spline_jobs;
    setenv("GENIE_GMKSPL_SPLINE_JOBS", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Using --jobs as per-spline worker count for the single "
      << "probe-target initial state: GENIE_GMKSPL_SPLINE_JOBS="
      << value.str();
  } else {
    setenv("GENIE_GMKSPL_SPLINE_JOBS", "1", 1);
  }
  if(gOptJobs > 1 && init_states.size() > 1) {
    BuildSplinesParallel(init_states);
  } else {
    BuildSplinesSerial(init_states);
  }

  // Save the splines at the requested XML file
  KPhaseSpaceCuts * cuts = KPhaseSpaceCuts::Instance();
  bool has_q2_metadata = false;
  double q2_metadata = -1.;
  string q2_metadata_source = "";
  PDGCodeList::const_iterator nuiter;
  for(nuiter = neutrinos->begin(); nuiter != neutrinos->end(); ++nuiter) {
    int nupdgc = *nuiter;
    if(!cuts->HasSplineQ2MinCutForProbe(nupdgc)) continue;

    double probe_q2_metadata = cuts->SplineQ2MinCutForProbe(nupdgc);
    string probe_q2_metadata_source =
      cuts->SplineQ2MinCutSourceForProbe(nupdgc);
    if(!has_q2_metadata) {
      has_q2_metadata = true;
      q2_metadata = probe_q2_metadata;
      q2_metadata_source = probe_q2_metadata_source;
    } else if(std::fabs(q2_metadata - probe_q2_metadata) > 1E-12) {
      LOG("gmkspl", pFATAL)
        << "Cannot write one spline Q2 metadata value for mixed probes "
        << "with different configured Q2 cuts.";
      exit(1);
    }
  }
  if(has_q2_metadata) {
    xspl->SetCurrentTuneQ2MinKinematics(q2_metadata, q2_metadata_source);
  }
  bool save_init = !gOptNoCopy;
  xspl->SaveAsXml(gOptOutXSecFile, save_init);

  delete neutrinos;
  delete targets;

  return 0;
}
//____________________________________________________________________________
double SecondsSince(const steady_clock::time_point & start)
{
  return duration_cast<duration<double>>(steady_clock::now() - start).count();
}
//____________________________________________________________________________
void ApplyFastQELOverrides(void)
{
  if(gOptFastQELThrows > 0) {
    std::ostringstream value;
    value << gOptFastQELThrows;
    setenv("GENIE_GMKSPL_FAST_QEL_THROWS", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec NumNucleonThrows via "
      << "GENIE_GMKSPL_FAST_QEL_THROWS=" << value.str();
  }
  if(gOptFastQELLowEMax > 0.) {
    std::ostringstream value;
    value << gOptFastQELLowEMax;
    setenv("GENIE_GMKSPL_FAST_QEL_LOW_E_MAX", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec low-energy throw cutoff via "
      << "GENIE_GMKSPL_FAST_QEL_LOW_E_MAX=" << value.str();
  }
  if(gOptFastQELLowEThrows > 0) {
    std::ostringstream value;
    value << gOptFastQELLowEThrows;
    setenv("GENIE_GMKSPL_FAST_QEL_LOW_E_THROWS", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec low-energy NumNucleonThrows via "
      << "GENIE_GMKSPL_FAST_QEL_LOW_E_THROWS=" << value.str();
  }
  if(gOptFastQELRelTol > 0.) {
    std::ostringstream value;
    value << gOptFastQELRelTol;
    setenv("GENIE_GMKSPL_FAST_QEL_REL_TOL", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec gsl-relative-tolerance via "
      << "GENIE_GMKSPL_FAST_QEL_REL_TOL=" << value.str();
  }
  if(gOptFastQELMaxEval > 0) {
    std::ostringstream value;
    value << gOptFastQELMaxEval;
    setenv("GENIE_GMKSPL_FAST_QEL_MAX_EVAL", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec gsl-max-eval via "
      << "GENIE_GMKSPL_FAST_QEL_MAX_EVAL=" << value.str();
  }
  if(gOptFastQELSkipKinCheck) {
    setenv("GENIE_GMKSPL_FAST_QEL_SKIP_KINEMATIC_CHECK", "1", 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding NewQELXSec kinematic checks via "
      << "GENIE_GMKSPL_FAST_QEL_SKIP_KINEMATIC_CHECK=1";
  }
}
//____________________________________________________________________________
void ApplyFastCOHOverrides(void)
{
  if(gOptFastCOHTSteps > 0) {
    std::ostringstream value;
    value << gOptFastCOHTSteps;
    setenv("GENIE_GMKSPL_FAST_COH_T_STEPS", value.str().c_str(), 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding BergerSehgalCOHPiPXSec2015 t quadrature steps via "
      << "GENIE_GMKSPL_FAST_COH_T_STEPS=" << value.str();
  }
  if(gOptFastCOHAnalyticT) {
    setenv("GENIE_GMKSPL_FAST_COH_ANALYTIC_T", "1", 1);
    LOG("gmkspl", pNOTICE)
      << "Overriding BergerSehgalCOHPiPXSec2015 t quadrature via "
      << "GENIE_GMKSPL_FAST_COH_ANALYTIC_T=1";
  }
}
//____________________________________________________________________________
vector<InitialState> GetInitialStates(
  PDGCodeList * neutrinos, PDGCodeList * targets)
{
  vector<InitialState> init_states;
  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;
  for(nuiter = neutrinos->begin(); nuiter != neutrinos->end(); ++nuiter) {
    for(tgtiter = targets->begin(); tgtiter != targets->end(); ++tgtiter) {
      int nupdgc  = *nuiter;
      int tgtpdgc = *tgtiter;
      init_states.push_back(InitialState(tgtpdgc, nupdgc));
    }
  }
  return init_states;
}
//____________________________________________________________________________
void BuildSplinesForInitialState(const InitialState & init_state)
{
  steady_clock::time_point total_start = steady_clock::now();
  LOG("gmkspl", pNOTICE)
    << "Configuring event-generation driver for initial state "
    << init_state.AsString();
  GEVGDriver driver;
  driver.SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  steady_clock::time_point configure_start = steady_clock::now();
  driver.Configure(init_state);
  LOG("gmkspl", pNOTICE)
    << "Configured event-generation driver for initial state "
    << init_state.AsString() << " in "
    << SecondsSince(configure_start) << " s";
  steady_clock::time_point spline_start = steady_clock::now();
  driver.CreateSplines(gOptNKnots, gOptMaxE);
  LOG("gmkspl", pNOTICE)
    << "Created splines for initial state " << init_state.AsString()
    << " in " << SecondsSince(spline_start)
    << " s; total initial-state time "
    << SecondsSince(total_start) << " s";
}
//____________________________________________________________________________
void BuildSplinesSerial(const vector<InitialState> & init_states)
{
  vector<InitialState>::const_iterator iter;
  for(iter = init_states.begin(); iter != init_states.end(); ++iter) {
    BuildSplinesForInitialState(*iter);
  }
}
//____________________________________________________________________________
string ParallelWorkerOutputName(size_t job_index)
{
  const char * tmpdir_env = std::getenv("TMPDIR");
  string tmpdir = (tmpdir_env && tmpdir_env[0] != '\0') ? tmpdir_env : ".";
  if(tmpdir[tmpdir.size()-1] != '/') tmpdir += "/";

  std::ostringstream filename;
  filename << tmpdir << "gmkspl-worker-" << getpid()
           << "-" << job_index << ".xml";
  return filename.str();
}
//____________________________________________________________________________
void BuildSplinesParallel(const vector<InitialState> & init_states)
{
  LOG("gmkspl", pNOTICE)
    << "Building splines with up to " << gOptJobs
    << " parallel gmkspl worker processes for "
    << init_states.size() << " initial states";

  struct WorkerInfo {
    WorkerInfo() : index(0), filename("") { }
    WorkerInfo(size_t i, const string & f) : index(i), filename(f) { }
    size_t index;
    string filename;
  };

  std::map<pid_t, WorkerInfo> workers;
  size_t next = 0;
  bool failed = false;

  while(next < init_states.size() || !workers.empty()) {
    while(!failed && next < init_states.size() &&
          (int)workers.size() < gOptJobs) {
      string worker_output = ParallelWorkerOutputName(next);
      pid_t pid = fork();

      if(pid < 0) {
        LOG("gmkspl", pERROR)
          << "Failed to fork gmkspl worker for initial state "
          << init_states[next].AsString();
        failed = true;
        break;
      }

      if(pid == 0) {
        LOG("gmkspl", pNOTICE)
          << "Worker building splines for initial state "
          << init_states[next].AsString();
        BuildSplinesForInitialState(init_states[next]);
        XSecSplineList::Instance()->SaveAsXml(worker_output, false);
        _exit(0);
      }

      workers.insert(std::map<pid_t, WorkerInfo>::value_type(
        pid, WorkerInfo(next, worker_output)));
      ++next;
    }

    if(workers.empty()) break;

    int status = 0;
    pid_t done = 0;
    do {
      done = waitpid(-1, &status, 0);
    } while(done < 0 && errno == EINTR);
    if(done < 0) {
      LOG("gmkspl", pERROR) << "Failed while waiting for gmkspl workers";
      failed = true;
      break;
    }

    std::map<pid_t, WorkerInfo>::iterator wit = workers.find(done);
    if(wit == workers.end()) continue;

    WorkerInfo info = wit->second;
    workers.erase(wit);

    bool exited_ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    if(!exited_ok) {
      LOG("gmkspl", pERROR)
        << "gmkspl worker failed for initial state index " << info.index;
      std::remove(info.filename.c_str());
      failed = true;
      continue;
    }

    XmlParserStatus_t load_status =
      XSecSplineList::Instance()->LoadFromXml(info.filename, true, false);
    std::remove(info.filename.c_str());
    if(load_status != kXmlOK) {
      LOG("gmkspl", pERROR)
        << "Could not import gmkspl worker output: " << info.filename;
      failed = true;
    }
  }

  if(failed) {
    std::map<pid_t, WorkerInfo>::const_iterator wit;
    for(wit = workers.begin(); wit != workers.end(); ++wit) {
      kill(wit->first, SIGTERM);
      std::remove(wit->second.filename.c_str());
    }
    gAbortingInErr = true;
    exit(1);
  }
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmkspl", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // output XML file name
  if( parser.OptionExists('o') ||
      parser.OptionExists("output-cross-sections") )
  {
    LOG("gmkspl", pINFO) << "Reading output filename";
    if( parser.OptionExists('o') ) {
      gOptOutXSecFile = parser.ArgAsString('o');
    }
    else {
      gOptOutXSecFile = parser.ArgAsString("output-cross-sections");
    }
  } else {
    LOG("gmkspl", pINFO) << "Unspecified filename - Using default";
    gOptOutXSecFile = "xsec_splines.xml";
  }

  // number of knots
  if( parser.OptionExists('n') ) {
    LOG("gmkspl", pINFO) << "Reading number of knots/spline";
    gOptNKnots = parser.ArgAsInt('n');
  } else {
    LOG("gmkspl", pINFO)
      << "Unspecified number of knots - Using default";
    gOptNKnots = -1;
  }

  // max spline energy (if < max of validity range)
  if( parser.OptionExists('e') ) {
    LOG("gmkspl", pINFO) << "Reading maximum spline energy";
    gOptMaxE = parser.ArgAsDouble('e');
  } else {
    LOG("gmkspl", pINFO)
       << "Unspecified maximum spline energy - Using default";
    gOptMaxE = -1;
  }

  // write out input splines?
  if( parser.OptionExists("no-copy") ) {
    LOG("gmkspl", pINFO) << "Not copying input splines to output";
    gOptNoCopy = true;
  }

  // comma-separated neutrino PDG code list
  if( parser.OptionExists('p') ) {
    LOG("gmkspl", pINFO) << "Reading neutrino PDG codes";
    gOptNuPdgCodeList = parser.ArgAsString('p');
  } else {
    LOG("gmkspl", pFATAL)
       << "Unspecified neutrino PDG code list - Exiting";
    PrintSyntax();
    exit(1);
  }

  // comma-separated target PDG code list or input geometry file
  bool tgt_cmd = true;
  if( parser.OptionExists('t') ) {
    LOG("gmkspl", pINFO) << "Reading target nuclei PDG codes";
    gOptTgtPdgCodeList = parser.ArgAsString('t');
  } else {
    LOG("gmkspl", pINFO) << "No code list specified from the command line";
    tgt_cmd = false;
  }

  bool tgt_geom = true;
  if( parser.OptionExists('f') ) {
    LOG("gmkspl", pINFO) << "Reading ROOT geometry filename";
    gOptGeomFilename = parser.ArgAsString('f');
  } else {
    LOG("gmkspl", pINFO) << "No geometry file was specified";
    tgt_cmd = false;
  }

  bool both =  tgt_geom &&  tgt_cmd;
  bool none = !tgt_geom && !tgt_cmd;
  if(none) {
    LOG("gmkspl", pFATAL)
          << "No geom file or cmd line target list was specified - Exiting";
    PrintSyntax();
    exit(1);
  }
  if(both) {
    LOG("gmkspl", pFATAL)
       << "You specified both a geom file and a cmd line target list "
         << "- Exiting confused";
    PrintSyntax();
    exit(1);
  }

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gmkspl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gmkspl", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // input cross-section file
  if( parser.OptionExists("input-cross-sections") ) {
    LOG("gmkspl", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("input-cross-sections");
  } else {
    LOG("gmkspl", pINFO) << "Unspecified input cross-section file";
    gOptInpXSecFile = "";
  }

  // low-risk fast spline evaluation mode
  gOptFastSplineEval = false;
  if(parser.OptionExists("fast-spline-eval")) {
    LOG("gmkspl", pINFO)
      << "Using fast spline evaluation mode";
    gOptFastSplineEval = true;
  }
  gOptFastSplineKnots = 0;
  if(parser.OptionExists("fast-spline-knots")) {
    gOptFastSplineKnots = parser.ArgAsInt("fast-spline-knots");
    if(gOptFastSplineKnots < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-spline-knots value: " << gOptFastSplineKnots
        << ". Sparse fast splines require at least 1 exact positive-energy knot.";
      PrintSyntax();
      exit(1);
    }
    LOG("gmkspl", pINFO)
      << "Using sparse fast spline mode with "
      << gOptFastSplineKnots << " exact above-threshold knots";
    gOptFastSplineEval = true;
  }
  gOptFastRESSplineKnots = 0;
  if(parser.OptionExists("fast-res-spline-knots")) {
    gOptFastRESSplineKnots = parser.ArgAsInt("fast-res-spline-knots");
    if(gOptFastRESSplineKnots < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-res-spline-knots value: "
        << gOptFastRESSplineKnots
        << ". RES sparse fast splines require at least 1 exact "
        << "positive-energy knot.";
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastDISSplineKnots = 0;
  if(parser.OptionExists("fast-dis-spline-knots")) {
    gOptFastDISSplineKnots = parser.ArgAsInt("fast-dis-spline-knots");
    if(gOptFastDISSplineKnots < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-knots value: "
        << gOptFastDISSplineKnots
        << ". DIS sparse fast splines require at least 1 "
        << "positive-energy knot.";
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastSplineChecks = 0;
  if(parser.OptionExists("fast-spline-checks")) {
    gOptFastSplineChecks = parser.ArgAsInt("fast-spline-checks");
    if(gOptFastSplineChecks < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-spline-checks value: " << gOptFastSplineChecks;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastRESSplineChecks = 0;
  if(parser.OptionExists("fast-res-spline-checks")) {
    gOptFastRESSplineChecks = parser.ArgAsInt("fast-res-spline-checks");
    if(gOptFastRESSplineChecks < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-res-spline-checks value: "
        << gOptFastRESSplineChecks;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastDISSplineChecks = 0;
  if(parser.OptionExists("fast-dis-spline-checks")) {
    gOptFastDISSplineChecks = parser.ArgAsInt("fast-dis-spline-checks");
    if(gOptFastDISSplineChecks < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-checks value: "
        << gOptFastDISSplineChecks;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastSplineTol = 0.25;
  if(parser.OptionExists("fast-spline-tol")) {
    gOptFastSplineTol = parser.ArgAsDouble("fast-spline-tol");
    if(gOptFastSplineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-spline-tol value: " << gOptFastSplineTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastSplineAreaTol = 0.;
  if(parser.OptionExists("fast-spline-area-tol")) {
    gOptFastSplineAreaTol = parser.ArgAsDouble("fast-spline-area-tol");
    if(gOptFastSplineAreaTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-spline-area-tol value: "
        << gOptFastSplineAreaTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastCOHSplineTol = 0.;
  if(parser.OptionExists("fast-coh-spline-tol")) {
    gOptFastCOHSplineTol = parser.ArgAsDouble("fast-coh-spline-tol");
    if(gOptFastCOHSplineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-coh-spline-tol value: "
        << gOptFastCOHSplineTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastRESSplineTol = 0.;
  if(parser.OptionExists("fast-res-spline-tol")) {
    gOptFastRESSplineTol = parser.ArgAsDouble("fast-res-spline-tol");
    if(gOptFastRESSplineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-res-spline-tol value: "
        << gOptFastRESSplineTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELSplineTol = 0.;
  if(parser.OptionExists("fast-qel-spline-tol")) {
    gOptFastQELSplineTol = parser.ArgAsDouble("fast-qel-spline-tol");
    if(gOptFastQELSplineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-spline-tol value: "
        << gOptFastQELSplineTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastDISSplineTol = 0.;
  if(parser.OptionExists("fast-dis-spline-tol")) {
    gOptFastDISSplineTol = parser.ArgAsDouble("fast-dis-spline-tol");
    if(gOptFastDISSplineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-tol value: "
        << gOptFastDISSplineTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastDISSplineAreaTol = 0.;
  if(parser.OptionExists("fast-dis-spline-area-tol")) {
    gOptFastDISSplineAreaTol =
      parser.ArgAsDouble("fast-dis-spline-area-tol");
    if(gOptFastDISSplineAreaTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-area-tol value: "
        << gOptFastDISSplineAreaTol;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastSplineChecks > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-spline-checks has no effect unless --fast-spline-knots "
      << "enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastRESSplineKnots > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-res-spline-knots has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineKnots > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-knots has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastRESSplineChecks > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-res-spline-checks has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineChecks > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-checks has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastSplineAreaTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-spline-area-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastCOHSplineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-coh-spline-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastRESSplineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-res-spline-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastQELSplineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-qel-spline-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineAreaTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-area-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastSplineRefineRounds = 0;
  if(parser.OptionExists("fast-spline-refine-rounds")) {
    gOptFastSplineRefineRounds =
      parser.ArgAsInt("fast-spline-refine-rounds");
    if(gOptFastSplineRefineRounds < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-spline-refine-rounds value: "
        << gOptFastSplineRefineRounds;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastSplineRefineRounds > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-spline-refine-rounds has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastRESSplineRefineRounds = 0;
  if(parser.OptionExists("fast-res-spline-refine-rounds")) {
    gOptFastRESSplineRefineRounds =
      parser.ArgAsInt("fast-res-spline-refine-rounds");
    if(gOptFastRESSplineRefineRounds < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-res-spline-refine-rounds value: "
        << gOptFastRESSplineRefineRounds;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastRESSplineRefineRounds > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-res-spline-refine-rounds has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastRESSplineRefineTol = 0.;
  if(parser.OptionExists("fast-res-spline-refine-tol")) {
    gOptFastRESSplineRefineTol =
      parser.ArgAsDouble("fast-res-spline-refine-tol");
    if(gOptFastRESSplineRefineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-res-spline-refine-tol value: "
        << gOptFastRESSplineRefineTol;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastRESSplineRefineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-res-spline-refine-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastQELSplineRefineRounds = 0;
  if(parser.OptionExists("fast-qel-spline-refine-rounds")) {
    gOptFastQELSplineRefineRounds =
      parser.ArgAsInt("fast-qel-spline-refine-rounds");
    if(gOptFastQELSplineRefineRounds < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-spline-refine-rounds value: "
        << gOptFastQELSplineRefineRounds;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastQELSplineRefineRounds > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-qel-spline-refine-rounds has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastDISSplineRefineRounds = 0;
  if(parser.OptionExists("fast-dis-spline-refine-rounds")) {
    gOptFastDISSplineRefineRounds =
      parser.ArgAsInt("fast-dis-spline-refine-rounds");
    if(gOptFastDISSplineRefineRounds < 0) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-refine-rounds value: "
        << gOptFastDISSplineRefineRounds;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineRefineRounds > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-refine-rounds has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastDISSplineRefineTol = 0.;
  if(parser.OptionExists("fast-dis-spline-refine-tol")) {
    gOptFastDISSplineRefineTol =
      parser.ArgAsDouble("fast-dis-spline-refine-tol");
    if(gOptFastDISSplineRefineTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-dis-spline-refine-tol value: "
        << gOptFastDISSplineRefineTol;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastSplineKnots == 0 && gOptFastDISSplineRefineTol > 0.) {
    LOG("gmkspl", pWARN)
      << "--fast-dis-spline-refine-tol has no effect unless "
      << "--fast-spline-knots enables sparse fast spline mode";
  }
  gOptFastSplineKnotJobs = 1;
  const char * knot_jobs_env = std::getenv("GENIE_GMKSPL_KNOT_JOBS");
  if(knot_jobs_env && knot_jobs_env[0] != '\0') {
    gOptFastSplineKnotJobs = atoi(knot_jobs_env);
  }
  if(parser.OptionExists("fast-spline-knot-jobs")) {
    gOptFastSplineKnotJobs =
      parser.ArgAsInt("fast-spline-knot-jobs");
  }
  if(gOptFastSplineKnotJobs < 1) {
    LOG("gmkspl", pFATAL)
      << "Invalid --fast-spline-knot-jobs value: "
      << gOptFastSplineKnotJobs
      << ". The exact-knot worker process count must be >= 1.";
    PrintSyntax();
    exit(1);
  }
  gOptFastSplineQuietKnotWorkers =
    parser.OptionExists("fast-spline-quiet-knot-workers");
  if(gOptFastSplineKnotJobs <= 1 && gOptFastSplineQuietKnotWorkers) {
    LOG("gmkspl", pWARN)
      << "--fast-spline-quiet-knot-workers has no effect unless "
      << "--fast-spline-knot-jobs uses worker processes";
  }

  gOptFastQELThrows = 0;
  if(parser.OptionExists("fast-qel-throws")) {
    gOptFastQELThrows = parser.ArgAsInt("fast-qel-throws");
    if(gOptFastQELThrows < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-throws value: " << gOptFastQELThrows;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELLowEThrows = 0;
  if(parser.OptionExists("fast-qel-low-e-throws")) {
    gOptFastQELLowEThrows = parser.ArgAsInt("fast-qel-low-e-throws");
    if(gOptFastQELLowEThrows < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-low-e-throws value: "
        << gOptFastQELLowEThrows;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELLowEMax = 0.;
  if(parser.OptionExists("fast-qel-low-e-max")) {
    gOptFastQELLowEMax = parser.ArgAsDouble("fast-qel-low-e-max");
    if(gOptFastQELLowEMax <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-low-e-max value: "
        << gOptFastQELLowEMax;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELRelTol = 0.;
  if(parser.OptionExists("fast-qel-rel-tol")) {
    gOptFastQELRelTol = parser.ArgAsDouble("fast-qel-rel-tol");
    if(gOptFastQELRelTol <= 0.) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-rel-tol value: " << gOptFastQELRelTol;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELMaxEval = 0;
  if(parser.OptionExists("fast-qel-max-eval")) {
    gOptFastQELMaxEval = parser.ArgAsInt("fast-qel-max-eval");
    if(gOptFastQELMaxEval < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-qel-max-eval value: " << gOptFastQELMaxEval;
      PrintSyntax();
      exit(1);
    }
  }
  gOptFastQELSkipKinCheck = parser.OptionExists("fast-qel-skip-kin-check");

  gOptFastCOHTSteps = 0;
  const char * coh_t_steps_env =
    std::getenv("GENIE_GMKSPL_FAST_COH_T_STEPS");
  if(coh_t_steps_env && coh_t_steps_env[0] != '\0') {
    gOptFastCOHTSteps = atoi(coh_t_steps_env);
  }
  if(parser.OptionExists("fast-coh-t-steps")) {
    gOptFastCOHTSteps = parser.ArgAsInt("fast-coh-t-steps");
    if(gOptFastCOHTSteps < 1) {
      LOG("gmkspl", pFATAL)
        << "Invalid --fast-coh-t-steps value: " << gOptFastCOHTSteps;
      PrintSyntax();
      exit(1);
    }
  }
  if(gOptFastCOHTSteps < 0) {
    LOG("gmkspl", pFATAL)
      << "Invalid --fast-coh-t-steps value: " << gOptFastCOHTSteps;
    PrintSyntax();
    exit(1);
  }

  gOptFastCOHAnalyticT = false;
  const char * coh_analytic_env =
    std::getenv("GENIE_GMKSPL_FAST_COH_ANALYTIC_T");
  if(coh_analytic_env && coh_analytic_env[0] != '\0' &&
     atoi(coh_analytic_env) != 0) {
    gOptFastCOHAnalyticT = true;
  }
  if(parser.OptionExists("fast-coh-analytic-t")) {
    gOptFastCOHAnalyticT = true;
  }
  if(gOptFastCOHAnalyticT && gOptFastCOHTSteps > 0) {
    LOG("gmkspl", pWARN)
      << "--fast-coh-t-steps is ignored by BergerSehgalCOHPiPXSec2015 "
      << "when --fast-coh-analytic-t is enabled.";
  }

  // Q2 minimum override
  gOptQ2MinSet = false;
  if(parser.OptionExists("q2-min")) {
    gOptQ2Min = parser.ArgAsDouble("q2-min");
    LOG("gmkspl", pINFO)
      << "Using command-line --q2-min override for the minimum Q2 "
      << "phase-space cut: " << gOptQ2Min << " GeV^2";
    gOptQ2MinSet = true;
  }

  const char * jobs_env = std::getenv("GENIE_GMKSPL_JOBS");
  if(jobs_env && jobs_env[0] != '\0') {
    gOptJobs = atoi(jobs_env);
  }
  if(parser.OptionExists("jobs")) {
    gOptJobs = parser.ArgAsInt("jobs");
  }
  if(gOptJobs < 1) {
    LOG("gmkspl", pFATAL)
      << "Invalid --jobs value: " << gOptJobs
      << ". The gmkspl worker process count must be >= 1.";
    PrintSyntax();
    exit(1);
  }

  //
  // print the command-line options
  //
  LOG("gmkspl", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gmkspl job configuration")
     << "\n Neutrino PDG codes : " << gOptNuPdgCodeList
     << "\n Target PDG codes : " << gOptTgtPdgCodeList
     << "\n Input ROOT geometry : " << gOptGeomFilename
     << "\n Output cross-section file : " << gOptOutXSecFile
     << "\n Input cross-section file : " << gOptInpXSecFile
     << "\n Random number seed : " << gOptRanSeed
     << "\n Fast spline evaluation : "
     << utils::print::BoolAsYNString(gOptFastSplineEval)
     << "\n Sparse fast spline exact knots : " << gOptFastSplineKnots
     << "\n Sparse fast spline RES exact knots : "
     << gOptFastRESSplineKnots
     << "\n Sparse fast spline DIS exact knots : "
     << gOptFastDISSplineKnots
     << "\n Sparse fast spline guard knots : " << gOptFastSplineChecks
     << "\n Sparse fast spline RES guard knots : "
     << gOptFastRESSplineChecks
     << "\n Sparse fast spline DIS guard knots : "
     << gOptFastDISSplineChecks
     << "\n Sparse fast spline guard tolerance : " << gOptFastSplineTol
     << "\n Sparse fast spline area tolerance : "
     << gOptFastSplineAreaTol
     << "\n Sparse fast spline COH guard tolerance : "
     << gOptFastCOHSplineTol
     << "\n Sparse fast spline RES guard tolerance : "
     << gOptFastRESSplineTol
     << "\n Sparse fast spline QEL guard tolerance : "
     << gOptFastQELSplineTol
     << "\n Sparse fast spline DIS guard tolerance : "
     << gOptFastDISSplineTol
     << "\n Sparse fast spline DIS area tolerance : "
     << gOptFastDISSplineAreaTol
     << "\n Sparse fast spline refine rounds : "
     << gOptFastSplineRefineRounds
     << "\n Sparse fast spline RES refine rounds : "
     << gOptFastRESSplineRefineRounds
     << "\n Sparse fast spline RES refine tolerance : "
     << gOptFastRESSplineRefineTol
     << "\n Sparse fast spline QEL refine rounds : "
     << gOptFastQELSplineRefineRounds
     << "\n Sparse fast spline DIS refine rounds : "
     << gOptFastDISSplineRefineRounds
     << "\n Sparse fast spline DIS refine tolerance : "
     << gOptFastDISSplineRefineTol
     << "\n Fast spline exact-knot worker processes : "
     << gOptFastSplineKnotJobs
     << "\n Fast spline quiet exact-knot workers : "
     << utils::print::BoolAsYNString(gOptFastSplineQuietKnotWorkers)
     << "\n Fast QEL nucleon throws override : " << gOptFastQELThrows
     << "\n Fast QEL low-energy throws override : " << gOptFastQELLowEThrows
     << "\n Fast QEL low-energy max override : " << gOptFastQELLowEMax
     << "\n Fast QEL relative tolerance override : " << gOptFastQELRelTol
     << "\n Fast QEL max eval override : " << gOptFastQELMaxEval
     << "\n Fast QEL skip kinematic check : "
     << utils::print::BoolAsYNString(gOptFastQELSkipKinCheck)
     << "\n Fast COH t quadrature steps override : " << gOptFastCOHTSteps
     << "\n Fast COH analytic t integration : "
     << utils::print::BoolAsYNString(gOptFastCOHAnalyticT)
     << "\n gmkspl worker processes : " << gOptJobs;
  if(gOptQ2MinSet) {
     LOG("gmkspl", pNOTICE)
       << "\n Q2 minimum override : " << gOptQ2Min << " GeV^2";
  }
  LOG("gmkspl", pNOTICE) << "\n";

  LOG("gmkspl", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmkspl", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gmkspl -p nupdg"
    << "\n    <-t tgtpdg, -f geomfile> "
    << "\n    <-o | --output-cross-sections> xsec_xml_file_name"
    << "\n    [-n nknots]"
    << "\n    [-e max_energy]"
    << "\n    [--no-copy]"
    << "\n    [--seed seed_number]"
    << "\n    [--input-cross-sections xml_file]"
    << "\n    [--fast-spline-eval]"
    << "\n    [--fast-spline-knots exact_knots]"
    << "\n    [--fast-res-spline-knots exact_knots]"
    << "\n    [--fast-dis-spline-knots exact_knots]"
    << "\n    [--fast-spline-checks guard_knots]"
    << "\n    [--fast-res-spline-checks guard_knots]"
    << "\n    [--fast-dis-spline-checks guard_knots]"
    << "\n    [--fast-spline-tol rel_tolerance]"
    << "\n    [--fast-spline-area-tol rel_tolerance]"
    << "\n    [--fast-coh-spline-tol rel_tolerance]"
    << "\n    [--fast-res-spline-tol rel_tolerance]"
    << "\n    [--fast-qel-spline-tol rel_tolerance]"
    << "\n    [--fast-dis-spline-tol rel_tolerance]"
    << "\n    [--fast-dis-spline-area-tol rel_tolerance]"
    << "\n    [--fast-spline-refine-rounds rounds]"
    << "\n    [--fast-res-spline-refine-rounds rounds]"
    << "\n    [--fast-res-spline-refine-tol rel_tolerance]"
    << "\n    [--fast-qel-spline-refine-rounds rounds]"
    << "\n    [--fast-dis-spline-refine-rounds rounds]"
    << "\n    [--fast-dis-spline-refine-tol rel_tolerance]"
    << "\n    [--fast-spline-knot-jobs worker_processes]"
    << "\n    [--fast-spline-quiet-knot-workers]"
    << "\n    [--fast-qel-throws nucleon_throws]"
    << "\n    [--fast-qel-low-e-throws nucleon_throws]"
    << "\n    [--fast-qel-low-e-max energy]"
    << "\n    [--fast-qel-rel-tol tolerance]"
    << "\n    [--fast-qel-max-eval evals]"
    << "\n    [--fast-qel-skip-kin-check]"
    << "\n    [--fast-coh-t-steps t_steps]"
    << "\n    [--fast-coh-analytic-t]"
    << "\n    [--q2-min q2_value]"
    << "\n    [--jobs worker_processes]"
    << RunOpt::RunOptSyntaxString(false)
    << "\n";

}
//____________________________________________________________________________
PDGCodeList * GetNeutrinoCodes(void)
{
  // split the comma separated list
  vector<string> nuvec = utils::str::Split(gOptNuPdgCodeList,  ",");

  // fill in the PDG code list
  PDGCodeList * list = new PDGCodeList;
  vector<string>::const_iterator iter;
  for(iter = nuvec.begin(); iter != nuvec.end(); ++iter) {
    list->push_back( atoi(iter->c_str()) );
  }
  return list;
}
//____________________________________________________________________________
PDGCodeList * GetTargetCodes(void)
{
  bool from_geom_file = ( gOptGeomFilename.size()   > 0 );
  bool from_cmd_line  = ( gOptTgtPdgCodeList.size() > 0 );

  if (from_cmd_line) {
     // split the comma separated list
     vector<string> tgtvec = utils::str::Split(gOptTgtPdgCodeList, ",");

     // fill in the PDG code list
     PDGCodeList * list = new PDGCodeList;
     vector<string>::const_iterator iter;
     for(iter = tgtvec.begin(); iter != tgtvec.end(); ++iter) {
        list->push_back( atoi(iter->c_str()) );
     }
     return list;
  }

  if (from_geom_file) {
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
     // create/configure a geometry driver
     LOG("gmkspl", pINFO) << "Creating/configuring a ROOT geom. driver";
     ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptGeomFilename);

     PDGCodeList * list = new PDGCodeList(geom->ListOfTargetNuclei());

     delete geom;
     return list;
#else
     LOG("gmkspl", pFATAL)
      << "To read-in a ROOT geometry you need to enable the geometry drivers!";
     gAbortingInErr = true;
     exit(1);
     return 0;
#endif

  }
  return 0;
}
//____________________________________________________________________________
