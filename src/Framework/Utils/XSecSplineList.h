//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 12, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _XSEC_SPLINE_LIST_H_
#define _XSEC_SPLINE_LIST_H_

#include <ostream>
#include <map>
#include <set>
#include <vector>
#include <string>

#include "Framework/Conventions/XmlParserStatus.h"

using std::map;
using std::set;
using std::pair;
using std::vector;
using std::string;
using std::ostream;

namespace genie {

class XSecAlgorithmI;
class Interaction;
class Spline;

class XSecSplineList;
ostream & operator << (ostream & stream, const XSecSplineList & xsl);

class XSecSplineList {

public:

  static XSecSplineList * Instance();

  // Save/load to/from XML file
  void               SaveAsXml   (const string & filename, bool save_init = true) const;
  XmlParserStatus_t  LoadFromXml (const string & filename, bool keep = false,
                                  bool mark_as_loaded = true);

  // Print available splines
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const XSecSplineList & xsl);

  // Set and query current tune.
  // An XSecSplineList can keep splines for numerous tunes and pick the appropriate
  // one for each process, as instructed.
  void   SetCurrentTune (const string & tune) { fCurrentTune = tune; }
  string CurrentTune    (void) const  { return fCurrentTune; }
  bool   HasSplineFromTune( const string & tune ) const { return fSplineMap.count(tune) > 0 ; }
  void   SetTuneQ2MinKinematics(const string & tune, double q2min, const string & source);
  bool   HasTuneQ2MinKinematics(const string & tune) const;
  bool   GetTuneQ2MinKinematics(const string & tune, double & q2min,
                                string * unit = 0, string * source = 0) const;
  void   SetCurrentTuneQ2MinKinematics(double q2min, const string & source)
         { this->SetTuneQ2MinKinematics(fCurrentTune, q2min, source); }
  bool   HasCurrentTuneQ2MinKinematics(void) const
         { return this->HasTuneQ2MinKinematics(fCurrentTune); }
  bool   GetCurrentTuneQ2MinKinematics(double & q2min,
                                string * unit = 0, string * source = 0) const
         { return this->GetTuneQ2MinKinematics(fCurrentTune, q2min, unit, source); }

  // Query the existence, access or create a spline
  // The results of the following methods depend on the current tune setting
  bool           SplineExists (const XSecAlgorithmI * alg, const Interaction * i) const;
  bool           SplineExists (string spline_key) const;
  const Spline * GetSpline    (const XSecAlgorithmI * alg, const Interaction * i) const;
  const Spline * GetSpline    (string spline_key) const;
  void           CreateSpline (const XSecAlgorithmI * alg, const Interaction * i,
                               int nknots = -1, double e_min = -1, double e_max = -1);
  int  NSplines (void) const;
  bool IsEmpty  (void) const;

  // Methods for building / getting keys
  // The results of the following methods depend on the current tune setting
  string BuildSplineKey(const XSecAlgorithmI * alg, const Interaction * i) const;
  const vector<string> * GetSplineKeys(void) const;


  // XSecSplineList options
  void   SetLogE   (bool   on); ///< set opt to build splines as f(E) or as f(logE)
  void   SetNKnots (int    nk); ///< set default number of knots for building the spline
  void   SetMinE   (double Ev); ///< set default minimum energy for xsec splines
  void   SetMaxE   (double Ev); ///< set default maximum energy for xsec splines
  void   SetFastSplineEval (bool on); ///< skip below-threshold integrations and summarize per-knot timing
  void   SetFastSplineKnots (int nk); ///< evaluate only this many above-threshold knots and interpolate the rest
  void   SetFastRESSplineKnots (int nk); ///< RES-only exact-knot override
  void   SetFastDISSplineKnots (int nk); ///< DIS-only exact-knot override
  void   SetFastSplineChecks (int nk); ///< exact guard knots for sparse fast splines
  void   SetFastRESSplineChecks (int nk); ///< RES-only exact guard knot override
  void   SetFastDISSplineChecks (int nk); ///< DIS-only exact guard knot override
  void   SetFastSplineTolerance (double tol); ///< max guard relative difference before exact fallback
  void   SetFastSplineAreaTolerance (double tol); ///< max guard area-difference before exact fallback
  void   SetFastCOHSplineTolerance (double tol); ///< COH-only sparse guard tolerance override
  void   SetFastRESSplineTolerance (double tol); ///< RES-only sparse guard tolerance override
  void   SetFastQELSplineTolerance (double tol); ///< QEL-only sparse guard tolerance override
  void   SetFastDISSplineTolerance (double tol); ///< DIS-only sparse guard tolerance override
  void   SetFastDISSplineAreaTolerance (double tol); ///< DIS-only sparse area tolerance override
  void   SetFastSplineRefineRounds (int nr); ///< extra sparse guard-promotion rounds before exact fallback
  void   SetFastRESSplineRefineRounds (int nr); ///< RES-only sparse guard-promotion rounds
  void   SetFastRESSplineRefineTolerance (double tol); ///< RES-only guard tolerance that triggers refinement
  void   SetFastQELSplineRefineRounds (int nr); ///< QEL-only sparse guard-promotion rounds
  void   SetFastDISSplineRefineRounds (int nr); ///< DIS-only sparse guard-promotion rounds
  void   SetFastDISSplineRefineTolerance (double tol); ///< DIS-only guard tolerance that triggers refinement
  void   SetFastSplineKnotJobs (int nj); ///< worker processes for exact knot integrations inside one spline
  void   SetFastSplineQuietKnotWorkers (bool on); ///< silence stdout/stderr in exact-knot worker children
  bool   UseLogE   (void) const { return fUseLogE;  }
  bool   FastSplineEval (void) const { return fFastSplineEval; }
  int    FastSplineKnots (void) const { return fFastSplineKnots; }
  int    FastRESSplineKnots (void) const { return fFastRESSplineKnots; }
  int    FastDISSplineKnots (void) const { return fFastDISSplineKnots; }
  int    FastSplineChecks (void) const { return fFastSplineChecks; }
  int    FastRESSplineChecks (void) const { return fFastRESSplineChecks; }
  int    FastDISSplineChecks (void) const { return fFastDISSplineChecks; }
  int    FastSplineRefineRounds (void) const { return fFastSplineRefineRounds; }
  int    FastSplineKnotJobs (void) const { return fFastSplineKnotJobs; }
  double FastSplineTolerance (void) const { return fFastSplineTolerance; }
  double FastSplineAreaTolerance (void) const { return fFastSplineAreaTolerance; }
  double FastCOHSplineTolerance (void) const { return fFastCOHSplineTolerance; }
  double FastRESSplineTolerance (void) const { return fFastRESSplineTolerance; }
  double FastQELSplineTolerance (void) const { return fFastQELSplineTolerance; }
  double FastDISSplineTolerance (void) const { return fFastDISSplineTolerance; }
  double FastDISSplineAreaTolerance (void) const { return fFastDISSplineAreaTolerance; }
  int    FastRESSplineRefineRounds (void) const { return fFastRESSplineRefineRounds; }
  double FastRESSplineRefineTolerance (void) const { return fFastRESSplineRefineTolerance; }
  int    FastQELSplineRefineRounds (void) const { return fFastQELSplineRefineRounds; }
  int    FastDISSplineRefineRounds (void) const { return fFastDISSplineRefineRounds; }
  double FastDISSplineRefineTolerance (void) const { return fFastDISSplineRefineTolerance; }
  bool   FastSplineQuietKnotWorkers (void) const { return fFastSplineQuietKnotWorkers; }
  int    NKnots    (void) const { return fNKnots;   }
  double Emin      (void) const { return fEmin;     }
  double Emax      (void) const { return fEmax;     }

private:

  XSecSplineList();
  XSecSplineList(const XSecSplineList & spline_list);
  virtual ~XSecSplineList();

  static XSecSplineList * fInstance;

  bool   fUseLogE;
  bool   fFastSplineEval;
  int    fFastSplineKnots;
  int    fFastRESSplineKnots;
  int    fFastDISSplineKnots;
  int    fFastSplineChecks;
  int    fFastRESSplineChecks;
  int    fFastDISSplineChecks;
  int    fFastSplineRefineRounds;
  int    fFastSplineKnotJobs;
  double fFastSplineTolerance;
  double fFastSplineAreaTolerance;
  double fFastCOHSplineTolerance;
  double fFastRESSplineTolerance;
  double fFastQELSplineTolerance;
  double fFastDISSplineTolerance;
  double fFastDISSplineAreaTolerance;
  int    fFastRESSplineRefineRounds;
  double fFastRESSplineRefineTolerance;
  int    fFastQELSplineRefineRounds;
  int    fFastDISSplineRefineRounds;
  double fFastDISSplineRefineTolerance;
  bool   fFastSplineQuietKnotWorkers;
  int    fNKnots;
  double fEmin;
  double fEmax;

  string fCurrentTune; ///< The `active' tune, out the many that can co-exist

  struct KinematicsMetadata {
    KinematicsMetadata() : has_q2_min(false), q2_min(-1.), q2_unit("GeV^2"), q2_source("") { }
    bool   has_q2_min;
    double q2_min;
    string q2_unit;
    string q2_source;
  };

  map<string, map<string, Spline *> > fSplineMap;       ///< tune -> { xsec_alg/xsec_config/interaction -> Spline }
  map<string, set<string>           > fLoadedSplineSet; ///< tune -> { set of initialy loaded splines             }
  map<string, KinematicsMetadata     > fKinematicsMetadata; ///< tune -> optional kinematic cut metadata           }

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (XSecSplineList::fInstance !=0) {
            delete XSecSplineList::fInstance;
            XSecSplineList::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _XSEC_SPLINE_LIST_H_
