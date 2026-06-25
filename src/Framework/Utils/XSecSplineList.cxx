//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <fenv.h>  //provides: int feenableexcept(int excepts);
#include <cmath>   //provides: std::isnan()

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <map>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include <signal.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::ofstream;
using std::endl;

namespace genie {

  using namespace std::chrono ;

namespace {

std::vector<int> SelectSparseIndices(int first, int last, int count)
{
  std::vector<int> indices;
  if(first > last || count <= 0) return indices;

  const int n = last - first + 1;
  if(count >= n) {
    for(int i = first; i <= last; ++i) indices.push_back(i);
    return indices;
  }

  std::vector<char> selected(n, false);
  if(count == 1) {
    selected[n - 1] = true;
  } else {
    // Sparse spline accuracy is most fragile just above threshold.
    // Bias exact knots there, while still keeping the high-energy endpoint.
    const double low_energy_bias = 3.0;
    for(int k = 0; k < count; ++k) {
      double frac = double(k) / double(count - 1);
      double pos = TMath::Power(frac, low_energy_bias) * double(n - 1);
      int offset = int(pos + 0.5);
      if(offset < 0) offset = 0;
      if(offset >= n) offset = n - 1;
      selected[offset] = true;
    }
  }

  for(int offset = 0; offset < n && int(indices.size()) < count; ++offset) {
    if(selected[offset]) indices.push_back(first + offset);
  }
  for(int offset = 0; offset < n && int(indices.size()) < count; ++offset) {
    if(!selected[offset]) indices.push_back(first + offset);
  }
  std::sort(indices.begin(), indices.end());
  return indices;
}

std::vector<int> SelectGuardIndices(
  int first, int last, const std::vector<int> & exact, int count)
{
  std::vector<int> guards;
  if(first > last || count <= 0) return guards;

  std::vector<int> candidates;
  for(int i = first; i <= last; ++i) {
    if(!std::binary_search(exact.begin(), exact.end(), i)) {
      candidates.push_back(i);
    }
  }
  if(candidates.empty()) return guards;

  if(count >= int(candidates.size())) return candidates;
  if(count == 1) {
    guards.push_back(candidates[candidates.size()/2]);
    return guards;
  }

  std::vector<char> selected(candidates.size(), false);
  for(int k = 0; k < count; ++k) {
    double pos = double(k) * double(candidates.size() - 1) / double(count - 1);
    int offset = int(pos + 0.5);
    if(offset < 0) offset = 0;
    if(offset >= int(candidates.size())) offset = int(candidates.size()) - 1;
    selected[offset] = true;
  }
  for(size_t i = 0; i < candidates.size() && int(guards.size()) < count; ++i) {
    if(selected[i]) guards.push_back(candidates[i]);
  }
  for(size_t i = 0; i < candidates.size() && int(guards.size()) < count; ++i) {
    if(!selected[i]) guards.push_back(candidates[i]);
  }
  std::sort(guards.begin(), guards.end());
  return guards;
}

double InterpCoord(double energy, bool use_log_e)
{
  return use_log_e ? TMath::Log10(energy) : energy;
}

double InterpolateSparseXSec(
  int index, const std::vector<double> & energy,
  const std::vector<double> & xsec, const std::vector<int> & exact,
  bool use_log_e)
{
  std::vector<int>::const_iterator upper =
    std::upper_bound(exact.begin(), exact.end(), index);

  if(upper == exact.begin()) return xsec[*upper];
  if(upper == exact.end()) return xsec[exact.back()];

  int hi = *upper;
  int lo = *(upper - 1);
  double x0 = InterpCoord(energy[lo], use_log_e);
  double x1 = InterpCoord(energy[hi], use_log_e);
  if(x1 == x0) return xsec[lo];

  double x = InterpCoord(energy[index], use_log_e);
  double t = (x - x0) / (x1 - x0);
  double y = xsec[lo] + t * (xsec[hi] - xsec[lo]);
  return TMath::Max(0., y);
}

double RelativeDifference(double a, double b)
{
  double denom = TMath::Max(TMath::Max(std::fabs(a), std::fabs(b)), 1.0e-300);
  return std::fabs(a - b) / denom;
}

double SparseArea(
  const std::vector<int> & indices,
  const std::vector<double> & energy,
  const std::vector<double> & xsec,
  const std::map<int, double> * overrides)
{
  if(indices.size() < 2) return 0.;

  double area = 0.;
  int last = indices.front();
  double last_y = xsec[last];
  if(overrides) {
    std::map<int, double>::const_iterator it = overrides->find(last);
    if(it != overrides->end()) last_y = it->second;
  }
  for(size_t i = 1; i < indices.size(); ++i) {
    int index = indices[i];
    double y = xsec[index];
    if(overrides) {
      std::map<int, double>::const_iterator it = overrides->find(index);
      if(it != overrides->end()) y = it->second;
    }
    area += 0.5 * (last_y + y) * (energy[index] - energy[last]);
    last = index;
    last_y = y;
  }
  return area;
}

struct KnotWorkerResult {
  int index;
  double xsec;
  double eval_time;
  int had_nan;
};

bool WriteAll(int fd, const void * data, size_t size)
{
  const char * ptr = static_cast<const char *>(data);
  size_t done = 0;
  while(done < size) {
    ssize_t nwrite = write(fd, ptr + done, size - done);
    if(nwrite < 0 && errno == EINTR) continue;
    if(nwrite <= 0) return false;
    done += size_t(nwrite);
  }
  return true;
}

bool ReadAll(int fd, void * data, size_t size)
{
  char * ptr = static_cast<char *>(data);
  size_t done = 0;
  while(done < size) {
    ssize_t nread = read(fd, ptr + done, size - done);
    if(nread < 0 && errno == EINTR) continue;
    if(nread <= 0) return false;
    done += size_t(nread);
  }
  return true;
}

} // namespace

//____________________________________________________________________________
ostream & operator << (ostream & stream, const XSecSplineList & list)
{
  list.Print(stream);
  return stream;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::fInstance = 0;
//____________________________________________________________________________
XSecSplineList::XSecSplineList()
{
  fInstance    =  0;
  fCurrentTune = "";
  fUseLogE     = true;
  fFastSplineEval = false;
  fFastSplineKnots = 0;
  fFastRESSplineKnots = 0;
  fFastDISSplineKnots = 0;
  fFastSplineChecks = 0;
  fFastRESSplineChecks = 0;
  fFastDISSplineChecks = 0;
  fFastSplineRefineRounds = 0;
  fFastSplineKnotJobs = 1;
  fFastSplineTolerance = 0.25;
  fFastSplineAreaTolerance = 0.;
  fFastCOHSplineTolerance = 0.;
  fFastRESSplineTolerance = 0.;
  fFastQELSplineTolerance = 0.;
  fFastDISSplineTolerance = 0.;
  fFastDISSplineAreaTolerance = 0.;
  fFastRESSplineRefineRounds = 0;
  fFastRESSplineRefineTolerance = 0.;
  fFastQELSplineRefineRounds = 0;
  fFastDISSplineRefineRounds = 0;
  fFastDISSplineRefineTolerance = 0.;
  fFastSplineQuietKnotWorkers = false;
  fNKnots      = 100;
  fEmin        =   0.01; // GeV
  fEmax        = 100.00; // GeV
}
//____________________________________________________________________________
XSecSplineList::~XSecSplineList()
{
// Clean up.

  map<string,  map<string, Spline *> >::iterator mm_iter = fSplineMap.begin();
  for( ; mm_iter != fSplineMap.end(); ++mm_iter) {
    // loop over splines for given tune
    map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::iterator m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      Spline * spline = m_iter->second;
      delete spline;
      spline = 0;
    }
    spl_map_curr_tune.clear();
  }
  fSplineMap.clear();
  fInstance = 0;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::Instance()
{
  if(fInstance == 0) {
    static XSecSplineList::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new XSecSplineList;
  }
  return fInstance;
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);
  return this->SplineExists(key);
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(string key) const
{

  if ( fCurrentTune.size() == 0 ) {
    SLOG("XSecSplLst", pERROR) << "Spline requested while CurrentTune not set" ;
    return false ;
  }

  SLOG("XSecSplLst", pDEBUG)
    << "Checking for spline: " << key << " in tune: " << fCurrentTune;

  map<string,  map<string, Spline *> >::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return false;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  bool exists = (spl_map_curr_tune.count(key) == 1);
  SLOG("XSecSplLst", pDEBUG)
    << "Spline found?...." << utils::print::BoolAsYNString(exists);
  return exists;
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);
  return this->GetSpline(key);
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(string key) const
{

  if ( fCurrentTune.size() == 0 ) {
    SLOG("XSecSplLst", pFATAL) << "Spline requested while CurrentTune not set" ;
    exit(0) ;
  }

  SLOG("XSecSplLst", pDEBUG)
    << "Getting spline: " << key << " in tune: " << fCurrentTune;

  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  map<string, Spline *>::const_iterator //\/
  m_iter = spl_map_curr_tune.find(key);
  if(m_iter == spl_map_curr_tune.end()) {
    SLOG("XSecSplLst", pWARN)
      << "Couldn't find spline: " << key << " in tune: " << fCurrentTune;
    return 0;
  }
  return m_iter->second;
}
//____________________________________________________________________________
void XSecSplineList::CreateSpline(const XSecAlgorithmI * alg,
        const Interaction * interaction, int nknots, double e_min, double e_max)
{
// Build a cross section spline for the input interaction using the input
// cross section algorithm and store in the list.
// For building this specific entry of the spline list, the user is allowed
// to override the list-wide nknots,e_min,e_max

  // FE_ALL_EXCEPT FE_INEXACT FE_UNDERFLOW
  // FE_DIVBYZERO  FE_INVALID FE_OVERFLO
  // rwh -- uncomment to catch NaN
  // feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);


  SLOG("XSecSplLst", pNOTICE)
     << "Creating cross section spline using the algorithm: " << *alg;

  string key = this->BuildSplineKey(alg,interaction);

  // If any of the nknots,e_min,e_max was not set or its value is not acceptable
  // use the list values
  //
  if (e_min   < 0.) e_min = this->Emin();
  if (e_max   < 0.) e_max = this->Emax();
  if (nknots <= 2) nknots = this->NKnots();
  assert( e_min < e_max );

  std::vector<double> xsec( nknots, 0. ) ;
  std::vector<double> E( nknots, 0. ) ;

  // Distribute the knots in the energy range (e_min,e_max) :
  // - Will use 5 knots linearly spaced below the energy thresholds so that the
  //   spline behaves correctly in (e_min,Ethr)
  // - Place 1 knot exactly on the input interaction threshold
  // - Place the remaining n-6 knots spaced either linearly or logarithmically
  //   above the input interaction threshold
  // The above scheme schanges appropriately if Ethr<e_min (i.e. no knots
  // are computed below threshold)
  //
  double Ethr = interaction->PhaseSpace().Threshold();
  SLOG("XSecSplLst", pNOTICE)
    << "Energy threshold for current interaction = " << Ethr << " GeV";

  if (Ethr>e_max) {
    SLOG("XSecSplLst", pFATAL) << "Energy threshold higher than maximum.";
    SLOG("XSecSplLst", pFATAL) << "Energy threshold = " << Ethr << " GeV";
    SLOG("XSecSplLst", pFATAL) << "Energy maximum = " << e_max << " GeV";
    return;
  }

  int nkb = (Ethr>e_min) ? 5 : 0; // number of knots <  threshold
  int nka = nknots-nkb;           // number of knots >= threshold

  // knots < energy threshold
  double dEb =  (Ethr>e_min) ? (Ethr - e_min) / nkb : 0;
  for(int i=0; i<nkb; i++) {
     E[i] = e_min + i*dEb;
  }
  // knots >= energy threshold
  double E0  = TMath::Max(Ethr,e_min);
  bool use_log_e = this->UseLogE();
  double dEa = 0;
  double log10_E0 = 0;
  if(use_log_e) {
    log10_E0 = TMath::Log10(E0);
    dEa = (TMath::Log10(e_max) - log10_E0) /(nka-1);
  }
  else
    dEa = (e_max-E0) /(nka-1);

  for(int i=0; i<nka; i++) {
     if(use_log_e)
       E[i+nkb] = TMath::Power(10., log10_E0 + i * dEa);
     else
       E[i+nkb] = E0 + i * dEa;
  }
  // force last point to avoid floating point cumulative slew
  E[nknots-1] = e_max;

  // Compute cross sections for the input interaction at the selected
  // set of energies
  //
  double pr_mass = interaction->InitStatePtr()->Probe()->Mass();
  double pr_mass2 = pr_mass * pr_mass;
  bool fast_eval = this->FastSplineEval();
  bool sparse_eval = (this->FastSplineKnots() > 0);
  int sparse_exact_knots = this->FastSplineKnots();
  double sparse_fallback_tolerance = this->FastSplineTolerance();
  double sparse_refine_tolerance = sparse_fallback_tolerance;
  double sparse_area_tolerance = this->FastSplineAreaTolerance();
  int sparse_guard_checks = this->FastSplineChecks();
  int sparse_refine_rounds = this->FastSplineRefineRounds();
  if(this->FastCOHSplineTolerance() > 0. &&
     interaction->ProcInfo().IsCoherentProduction()) {
    sparse_fallback_tolerance = this->FastCOHSplineTolerance();
    sparse_refine_tolerance = sparse_fallback_tolerance;
  } else if(this->FastRESSplineTolerance() > 0. &&
            interaction->ProcInfo().IsResonant()) {
    sparse_fallback_tolerance = this->FastRESSplineTolerance();
    sparse_refine_tolerance = sparse_fallback_tolerance;
  } else if(this->FastQELSplineTolerance() > 0. &&
            interaction->ProcInfo().IsQuasiElastic()) {
    sparse_fallback_tolerance = this->FastQELSplineTolerance();
    sparse_refine_tolerance = sparse_fallback_tolerance;
  } else if(this->FastDISSplineTolerance() > 0. &&
            interaction->ProcInfo().IsDeepInelastic()) {
    sparse_fallback_tolerance = this->FastDISSplineTolerance();
    sparse_refine_tolerance = sparse_fallback_tolerance;
  }
  if(this->FastDISSplineAreaTolerance() > 0. &&
     interaction->ProcInfo().IsDeepInelastic()) {
    sparse_area_tolerance = this->FastDISSplineAreaTolerance();
  }
  if(this->FastRESSplineChecks() > 0 &&
     interaction->ProcInfo().IsResonant()) {
    sparse_guard_checks = this->FastRESSplineChecks();
  }
  if(this->FastDISSplineChecks() > 0 &&
     interaction->ProcInfo().IsDeepInelastic()) {
    sparse_guard_checks = this->FastDISSplineChecks();
  }
  if(this->FastRESSplineKnots() > 0 &&
     interaction->ProcInfo().IsResonant()) {
    sparse_exact_knots = this->FastRESSplineKnots();
  }
  if(this->FastDISSplineKnots() > 0 &&
     interaction->ProcInfo().IsDeepInelastic()) {
    sparse_exact_knots = this->FastDISSplineKnots();
  }
  if(interaction->ProcInfo().IsResonant()) {
    if(this->FastRESSplineRefineRounds() > 0) {
      sparse_refine_rounds = this->FastRESSplineRefineRounds();
    }
    if(this->FastRESSplineRefineTolerance() > 0.) {
      sparse_refine_tolerance = this->FastRESSplineRefineTolerance();
    }
  }
  if(this->FastQELSplineRefineRounds() > 0 &&
     interaction->ProcInfo().IsQuasiElastic()) {
    sparse_refine_rounds = this->FastQELSplineRefineRounds();
  }
  if(interaction->ProcInfo().IsDeepInelastic()) {
    if(this->FastDISSplineRefineRounds() > 0) {
      sparse_refine_rounds = this->FastDISSplineRefineRounds();
    }
    if(this->FastDISSplineRefineTolerance() > 0.) {
      sparse_refine_tolerance = this->FastDISSplineRefineTolerance();
    }
  }
  if(sparse_refine_tolerance > sparse_fallback_tolerance) {
    sparse_refine_tolerance = sparse_fallback_tolerance;
  }
  int knot_jobs = this->FastSplineKnotJobs();
  int n_integrated = 0;
  int n_skipped_below_threshold = 0;
  int n_interpolated = 0;
  double total_integral_time = 0.;
  steady_clock::time_point spline_start = steady_clock::now();

  std::vector<char> exact_knot(nknots, false);
  std::vector<int> sparse_indices;
  std::vector<int> interpolation_indices;
  std::vector<int> guard_indices;
  int threshold_index = -1;
  int first_positive_index = -1;
  bool have_threshold_anchor = false;
  if(sparse_eval) {
    for(int i = 0; i < nknots; ++i) {
      if(E[i] <= Ethr) {
        threshold_index = i;
      } else {
        break;
      }
    }
    have_threshold_anchor = (threshold_index >= 0);
    first_positive_index = have_threshold_anchor ?
      TMath::Min(threshold_index + 1, nknots - 1) : 0;
    sparse_indices = SelectSparseIndices(first_positive_index, nknots - 1,
      TMath::Max(1, sparse_exact_knots));
    if(have_threshold_anchor) {
      interpolation_indices.push_back(threshold_index);
    }
    interpolation_indices.insert(
      interpolation_indices.end(), sparse_indices.begin(), sparse_indices.end());
    std::sort(interpolation_indices.begin(), interpolation_indices.end());
    interpolation_indices.erase(
      std::unique(interpolation_indices.begin(), interpolation_indices.end()),
      interpolation_indices.end());
    guard_indices = SelectGuardIndices(
      first_positive_index, nknots - 1, sparse_indices,
      sparse_guard_checks);
    if(have_threshold_anchor) {
      exact_knot[threshold_index] = true;
    }
    for(size_t i = 0; i < sparse_indices.size(); ++i) {
      exact_knot[sparse_indices[i]] = true;
    }
  }

  auto set_probe_p4 = [&](int i) {
    TLorentzVector p4(0,0,E[i],E[i]);
    if (pr_mass > 0.) {
      double pz = TMath::Max(0.,E[i]*E[i] - pr_mass2);
      pz = TMath::Sqrt(pz);
      p4.SetPz(pz);
    }
    interaction->InitStatePtr()->SetProbeP4(p4);
  };

  auto record_integral = [&](int i, double eval_time, bool had_nan) {
    total_integral_time += eval_time;
    ++n_integrated;

    if (!fast_eval && !sparse_eval) {
      SLOG("XSecSplLst", pNOTICE)
                         << "xsec(E = " << E[i] << ") =  "
                         << (1E+38/units::cm2)*xsec[i] << " x 1E-38 cm^2, evaluated in " << eval_time << " s";
    }
    if (sparse_eval) {
      SLOG("XSecSplLst", pNOTICE)
                         << "sparse exact xsec(knot = " << i
                         << ", E = " << E[i] << ") = "
                         << (1E+38/units::cm2)*xsec[i]
                         << " x 1E-38 cm^2, evaluated in "
                         << eval_time << " s";
    }
    if (had_nan) {
      SLOG("XSecSplLst", pWARN)
                       << "xsec(E = " << E[i] << ") is NaN"
                       << " : converting NaN to 0.0";
    }
  };

  auto evaluate_integral = [&](int i) -> double {
    set_probe_p4(i);

    if (fast_eval && E[i] <= Ethr) {
      xsec[i] = 0.0;
      ++n_skipped_below_threshold;
      return xsec[i];
    }

    steady_clock::time_point start = steady_clock::now();
    xsec[i] = alg->Integral(interaction);
    steady_clock::time_point end = steady_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    double eval_time = time_span.count();
    bool had_nan = std::isnan(xsec[i]);
    if ( had_nan ) {
      xsec[i] = 0.0;
    }
    record_integral(i, eval_time, had_nan);

    return xsec[i];
  };

  auto evaluate_integrals_parallel =
    [&](const std::vector<int> & indices) -> bool {
    if(knot_jobs <= 1 || indices.size() <= 1) return false;

    struct KnotWorkerInfo {
      KnotWorkerInfo() : index(-1), read_fd(-1) { }
      KnotWorkerInfo(int i, int fd) : index(i), read_fd(fd) { }
      int index;
      int read_fd;
    };

    std::map<pid_t, KnotWorkerInfo> workers;
    std::map<int, KnotWorkerResult> results;
    size_t next = 0;
    bool failed = false;

    while(next < indices.size() || !workers.empty()) {
      while(!failed && next < indices.size() &&
            int(workers.size()) < knot_jobs) {
        int index = indices[next];
        int fds[2] = {-1, -1};
        if(pipe(fds) != 0) {
          SLOG("XSecSplLst", pERROR)
            << "Failed to create pipe for spline-knot worker";
          failed = true;
          break;
        }

        pid_t pid = fork();
        if(pid < 0) {
          close(fds[0]);
          close(fds[1]);
          SLOG("XSecSplLst", pERROR)
            << "Failed to fork spline-knot worker for knot " << index;
          failed = true;
          break;
        }

        if(pid == 0) {
          close(fds[0]);
          if(this->FastSplineQuietKnotWorkers()) {
            int null_fd = open("/dev/null", O_WRONLY);
            if(null_fd >= 0) {
              dup2(null_fd, STDOUT_FILENO);
              dup2(null_fd, STDERR_FILENO);
              if(null_fd > STDERR_FILENO) close(null_fd);
            }
          }
          set_probe_p4(index);

          KnotWorkerResult result;
          result.index = index;
          result.xsec = 0.;
          result.eval_time = 0.;
          result.had_nan = 0;

          steady_clock::time_point start = steady_clock::now();
          result.xsec = alg->Integral(interaction);
          steady_clock::time_point end = steady_clock::now();
          duration<double> time_span =
            duration_cast<duration<double> >(end - start);
          result.eval_time = time_span.count();
          if(std::isnan(result.xsec)) {
            result.xsec = 0.;
            result.had_nan = 1;
          }

          bool wrote = WriteAll(fds[1], &result, sizeof(result));
          close(fds[1]);
          _exit(wrote ? 0 : 2);
        }

        close(fds[1]);
        workers.insert(std::map<pid_t, KnotWorkerInfo>::value_type(
          pid, KnotWorkerInfo(index, fds[0])));
        ++next;
      }

      if(workers.empty()) break;

      int status = 0;
      pid_t done = 0;
      do {
        done = waitpid(-1, &status, 0);
      } while(done < 0 && errno == EINTR);
      if(done < 0) {
        SLOG("XSecSplLst", pERROR)
          << "Failed while waiting for spline-knot workers";
        failed = true;
        break;
      }

      std::map<pid_t, KnotWorkerInfo>::iterator wit = workers.find(done);
      if(wit == workers.end()) continue;

      KnotWorkerInfo info = wit->second;
      workers.erase(wit);

      KnotWorkerResult result;
      bool read_ok = ReadAll(info.read_fd, &result, sizeof(result));
      close(info.read_fd);

      bool exited_ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;
      if(!read_ok || !exited_ok || result.index != info.index) {
        SLOG("XSecSplLst", pERROR)
          << "Spline-knot worker failed for knot " << info.index;
        failed = true;
        continue;
      }

      results.insert(std::map<int, KnotWorkerResult>::value_type(
        result.index, result));
    }

    if(failed) {
      std::map<pid_t, KnotWorkerInfo>::const_iterator wit;
      for(wit = workers.begin(); wit != workers.end(); ++wit) {
        kill(wit->first, SIGTERM);
        if(wit->second.read_fd >= 0) close(wit->second.read_fd);
      }
      return false;
    }

    for(size_t i = 0; i < indices.size(); ++i) {
      int index = indices[i];
      std::map<int, KnotWorkerResult>::const_iterator rit =
        results.find(index);
      if(rit == results.end()) return false;
      const KnotWorkerResult & result = rit->second;
      xsec[index] = result.xsec;
      record_integral(index, result.eval_time, result.had_nan != 0);
    }
    return true;
  };

  auto evaluate_integrals = [&](const std::vector<int> & indices) {
    if(indices.empty()) return;
    if(evaluate_integrals_parallel(indices)) return;
    if(knot_jobs > 1 && indices.size() > 1) {
      SLOG("XSecSplLst", pWARN)
        << "Falling back to serial exact knot integrations";
    }
    for(size_t i = 0; i < indices.size(); ++i) {
      evaluate_integral(indices[i]);
    }
  };

  if(sparse_eval) {
    for (int i = 0; i < nknots; i++) {
      if (E[i] <= Ethr) {
        xsec[i] = 0.0;
        exact_knot[i] = true;
        ++n_skipped_below_threshold;
      }
    }

    evaluate_integrals(sparse_indices);

    auto promote_guards = [&]() {
      interpolation_indices.insert(
        interpolation_indices.end(), guard_indices.begin(), guard_indices.end());
      std::sort(interpolation_indices.begin(), interpolation_indices.end());
      interpolation_indices.erase(
        std::unique(interpolation_indices.begin(), interpolation_indices.end()),
        interpolation_indices.end());
    };

    auto refresh_interpolation = [&]() {
      n_interpolated = 0;
      for(int i = nkb; i < nknots; ++i) {
        if(exact_knot[i]) continue;
        xsec[i] = InterpolateSparseXSec(
          i, E, xsec, interpolation_indices, use_log_e);
        ++n_interpolated;
      }
    };

    auto evaluate_guards =
      [&](double & worst_rel, int & worst_index, double & area_rel) {
      worst_rel = 0.;
      worst_index = -1;
      area_rel = 0.;
      std::vector<double> approx(guard_indices.size(), 0.);
      std::map<int, double> approx_map;
      for(size_t i = 0; i < guard_indices.size(); ++i) {
        int index = guard_indices[i];
        approx[i] = xsec[index];
        approx_map[index] = approx[i];
      }
      evaluate_integrals(guard_indices);
      for(size_t i = 0; i < guard_indices.size(); ++i) {
        int index = guard_indices[i];
        double exact = xsec[index];
        exact_knot[index] = true;
        double rel = RelativeDifference(exact, approx[i]);
        if(rel > worst_rel) {
          worst_rel = rel;
          worst_index = index;
        }
      }
      if(sparse_area_tolerance > 0. && !guard_indices.empty()) {
        std::vector<int> area_indices = interpolation_indices;
        area_indices.insert(area_indices.end(),
          guard_indices.begin(), guard_indices.end());
        std::sort(area_indices.begin(), area_indices.end());
        area_indices.erase(std::unique(area_indices.begin(), area_indices.end()),
          area_indices.end());
        double approx_area = SparseArea(area_indices, E, xsec, &approx_map);
        double exact_area = SparseArea(area_indices, E, xsec, 0);
        area_rel = RelativeDifference(exact_area, approx_area);
      }
    };

    refresh_interpolation();

    double worst_guard_rel = 0.;
    int worst_guard_index = -1;
    double guard_area_rel = 0.;
    int guard_anchor_count = 0;
    int refine_rounds_done = 0;
    evaluate_guards(worst_guard_rel, worst_guard_index, guard_area_rel);
    guard_anchor_count += guard_indices.size();

    while(!guard_indices.empty() &&
          (worst_guard_rel > sparse_refine_tolerance ||
           (sparse_area_tolerance > 0. &&
            guard_area_rel > sparse_area_tolerance)) &&
          refine_rounds_done < sparse_refine_rounds) {
      promote_guards();
      refresh_interpolation();
      ++refine_rounds_done;
      guard_indices = SelectGuardIndices(
        first_positive_index, nknots - 1, interpolation_indices,
        sparse_guard_checks);
      evaluate_guards(worst_guard_rel, worst_guard_index, guard_area_rel);
      guard_anchor_count += guard_indices.size();
    }

    bool fell_back = false;
    if(!guard_indices.empty() &&
       (worst_guard_rel > sparse_fallback_tolerance ||
        (sparse_area_tolerance > 0. &&
         guard_area_rel > sparse_area_tolerance))) {
      fell_back = true;
      n_interpolated = 0;
      std::vector<int> fallback_indices;
      for(int i = nkb; i < nknots; ++i) {
        if(exact_knot[i]) continue;
        fallback_indices.push_back(i);
      }
      evaluate_integrals(fallback_indices);
      for(size_t i = 0; i < fallback_indices.size(); ++i) {
        int index = fallback_indices[i];
        exact_knot[index] = true;
      }
    }
    if(!fell_back) {
      promote_guards();
      refresh_interpolation();
    }

    steady_clock::time_point spline_end = steady_clock::now();
    duration<double> spline_time =
      duration_cast<duration<double>>(spline_end - spline_start);
    SLOG("XSecSplLst", pNOTICE)
      << "Sparse fast spline " << key << ": evaluated " << n_integrated
      << "/" << nknots << " knots in " << total_integral_time
      << " s integral time (" << spline_time.count() << " s wall time), "
      << "interpolated " << n_interpolated << ", skipped "
      << n_skipped_below_threshold << " at/below-threshold knots, guard worst_rel="
      << worst_guard_rel << " at knot " << worst_guard_index
      << ", guard tolerance=" << sparse_fallback_tolerance
      << ", guard area_rel=" << guard_area_rel
      << ", area tolerance=" << sparse_area_tolerance
      << ", refine tolerance=" << sparse_refine_tolerance
      << ", exact sparse knots=" << sparse_exact_knots
      << ", guard checks=" << sparse_guard_checks
      << ", guard anchors=" << guard_anchor_count
      << ", refine rounds=" << refine_rounds_done
      << "/" << sparse_refine_rounds
      << ", fallback=" << utils::print::BoolAsYNString(fell_back);
  } else {
    std::vector<int> exact_indices;
    for (int i = 0; i < nknots; i++) {
      if (fast_eval && E[i] <= Ethr) {
        xsec[i] = 0.0;
        ++n_skipped_below_threshold;
        continue;
      }
      exact_indices.push_back(i);
    }
    evaluate_integrals(exact_indices);
  }
  if (fast_eval && !sparse_eval) {
    steady_clock::time_point spline_end = steady_clock::now();
    duration<double> spline_time = duration_cast<duration<double>>(spline_end - spline_start);
    SLOG("XSecSplLst", pNOTICE)
      << "Spline " << key << ": evaluated " << n_integrated
      << "/" << nknots << " knots in " << total_integral_time
      << " s integral time (" << spline_time.count() << " s wall time), skipped "
      << n_skipped_below_threshold << " at/below-threshold knots";
  }

  // Warn about odd case of decreasing cross section
  //    but allow for small variation due to integration errors
  const double eps_xsec = 1.0e-5;
  const double xsec_scale = (1.0-eps_xsec);
  if ( xsec[nknots-1] < xsec[nknots-2]*xsec_scale ) {
    SLOG("XSecSplLst", pWARN)
      << "Last point oddity: " << key <<  " has "
      << " xsec[nknots-1] " << xsec[nknots-1] << " < "
      << " xsec[nknots-2] " << xsec[nknots-2];
  }

  // Build
  //
  Spline * spline = new Spline(nknots, E.data(), xsec.data());

  // Save
  //
  map<string,  map<string, Spline *> >::iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    map<string, Spline *> spl_map_curr_tune;
    fSplineMap.insert( map<string, map<string, Spline *> >::value_type(
      fCurrentTune, spl_map_curr_tune) );
    mm_iter = fSplineMap.find(fCurrentTune);
  }
  map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  spl_map_curr_tune.insert( map<string, Spline *>::value_type(key, spline) );
}
//____________________________________________________________________________
int XSecSplineList::NSplines(void) const
{
  map<string,  map<string, Spline *> >::const_iterator //
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  return (int) spl_map_curr_tune.size();
}
//____________________________________________________________________________
bool XSecSplineList::IsEmpty(void) const
{
  int n = this->NSplines();
  return (n == 0);
}
//____________________________________________________________________________
void XSecSplineList::SetLogE(bool on)
{
  fUseLogE = on;
}
//____________________________________________________________________________
void XSecSplineList::SetNKnots(int nk)
{
  fNKnots = nk;
  if(fNKnots<10) fNKnots = 10; // minimum acceptable number of knots
}
//____________________________________________________________________________
void XSecSplineList::SetMinE(double Ev)
{
  if(Ev>0) fEmin = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SetMaxE(double Ev)
{
  if(Ev>0) fEmax = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineEval(bool on)
{
  fFastSplineEval = on;
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineKnots(int nk)
{
  fFastSplineKnots = nk;
  if(fFastSplineKnots < 0) fFastSplineKnots = 0;
}
//____________________________________________________________________________
void XSecSplineList::SetFastRESSplineKnots(int nk)
{
  fFastRESSplineKnots = nk;
  if(fFastRESSplineKnots < 0) fFastRESSplineKnots = 0;
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineKnots(int nk)
{
  fFastDISSplineKnots = nk;
  if(fFastDISSplineKnots < 0) fFastDISSplineKnots = 0;
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineChecks(int nk)
{
  fFastSplineChecks = TMath::Max(0, nk);
}
//____________________________________________________________________________
void XSecSplineList::SetFastRESSplineChecks(int nk)
{
  fFastRESSplineChecks = TMath::Max(0, nk);
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineChecks(int nk)
{
  fFastDISSplineChecks = TMath::Max(0, nk);
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineRefineRounds(int nr)
{
  fFastSplineRefineRounds = TMath::Max(0, nr);
}
//____________________________________________________________________________
void XSecSplineList::SetFastRESSplineRefineRounds(int nr)
{
  fFastRESSplineRefineRounds = TMath::Max(0, nr);
}
//____________________________________________________________________________
void XSecSplineList::SetFastQELSplineRefineRounds(int nr)
{
  fFastQELSplineRefineRounds = TMath::Max(0, nr);
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineRefineRounds(int nr)
{
  fFastDISSplineRefineRounds = TMath::Max(0, nr);
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineKnotJobs(int nj)
{
  fFastSplineKnotJobs = TMath::Max(1, nj);
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineQuietKnotWorkers(bool on)
{
  fFastSplineQuietKnotWorkers = on;
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineTolerance(double tol)
{
  if(tol > 0.) fFastSplineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastSplineAreaTolerance(double tol)
{
  if(tol > 0.) fFastSplineAreaTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastCOHSplineTolerance(double tol)
{
  if(tol > 0.) fFastCOHSplineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastRESSplineTolerance(double tol)
{
  if(tol > 0.) fFastRESSplineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastQELSplineTolerance(double tol)
{
  if(tol > 0.) fFastQELSplineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineTolerance(double tol)
{
  if(tol > 0.) fFastDISSplineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineAreaTolerance(double tol)
{
  if(tol > 0.) fFastDISSplineAreaTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastRESSplineRefineTolerance(double tol)
{
  if(tol > 0.) fFastRESSplineRefineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetFastDISSplineRefineTolerance(double tol)
{
  if(tol > 0.) fFastDISSplineRefineTolerance = tol;
}
//____________________________________________________________________________
void XSecSplineList::SetTuneQ2MinKinematics(
  const string & tune, double q2min, const string & source)
{
  if(tune.empty()) return;

  KinematicsMetadata & meta = fKinematicsMetadata[tune];
  meta.has_q2_min = true;
  meta.q2_min     = q2min;
  meta.q2_unit    = "GeV^2";
  meta.q2_source  = source;
}
//____________________________________________________________________________
bool XSecSplineList::HasTuneQ2MinKinematics(const string & tune) const
{
  map<string, KinematicsMetadata>::const_iterator it =
    fKinematicsMetadata.find(tune);
  return it != fKinematicsMetadata.end() && it->second.has_q2_min;
}
//____________________________________________________________________________
bool XSecSplineList::GetTuneQ2MinKinematics(
  const string & tune, double & q2min, string * unit, string * source) const
{
  map<string, KinematicsMetadata>::const_iterator it =
    fKinematicsMetadata.find(tune);
  if(it == fKinematicsMetadata.end() || !it->second.has_q2_min) return false;

  q2min = it->second.q2_min;
  if(unit)   *unit   = it->second.q2_unit;
  if(source) *source = it->second.q2_source;
  return true;
}
//____________________________________________________________________________
void XSecSplineList::SaveAsXml(const string & filename, bool save_init) const
{
//! Save XSecSplineList to XML file

  SLOG("XSecSplLst", pNOTICE)
       << "Saving XSecSplineList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    SLOG("XSecSplLst", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by genie::XSecSplineList::SaveSplineList() -->";
  outxml << endl << endl;

  int uselog = (fUseLogE ? 1 : 0);
  outxml << "<genie_xsec_spline_list "
                << "version=\"3.00\" uselog=\"" << uselog << "\">";
  outxml << endl << endl;

  // loop over tunes
  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.begin();
  for( ; mm_iter != fSplineMap.end(); ++mm_iter) {

    string tune_name = mm_iter->first;
    outxml << "  <genie_tune name=\"" << tune_name << "\">";
    outxml << endl << endl;

    map<string, KinematicsMetadata>::const_iterator kiter =
      fKinematicsMetadata.find(tune_name);
    if(kiter != fKinematicsMetadata.end() && kiter->second.has_q2_min) {
      outxml << "    <kinematics>" << endl;
      outxml << "      <q2 min=\"" << std::setprecision(16)
             << kiter->second.q2_min
             << "\" unit=\"" << kiter->second.q2_unit
             << "\" source=\"" << kiter->second.q2_source
             << "\"/>" << endl;
      outxml << "    </kinematics>" << endl << endl;
    }

    // loop over splines for given tune
    const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::const_iterator //\/
    m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      string key = m_iter->first;

      // If current spline is from the initial loaded set,
      // look-up input option to decide whether to write out in
      // new output file or not
      bool from_init_set = false;
      map<string, set<string> >::const_iterator //\/
      it = fLoadedSplineSet.find(tune_name);
      if(it != fLoadedSplineSet.end()) {
         const set<string> & init_set_curr_tune = it->second;
         from_init_set = (init_set_curr_tune.count(key) == 1);
      }
      if(from_init_set && !save_init) continue;

      // Add current spline to output file
      Spline * spline = m_iter->second;
      spline->SaveAsXml(outxml,"E","xsec", key);
    }//spline loop

    outxml << "  </genie_tune>" << endl;
  }//tune loop

  outxml << "</genie_xsec_spline_list>" << endl;

  outxml.close();
}
//____________________________________________________________________________
XmlParserStatus_t XSecSplineList::LoadFromXml(
  const string & filename, bool keep, bool mark_as_loaded)
{
//! Load XSecSplineList from ROOT file. If keep = true, then the loaded splines
//! are added to the existing list. If false, then the existing list is reset
//! before loading the splines.

  SLOG("XSecSplLst", pNOTICE)
    << "Loading splines from: " << filename;
  SLOG("XSecSplLst", pINFO)
    << "Option to keep pre-existing splines is switched "
    << ( (keep) ? "ON" : "OFF" );

  if(!keep) {
    fSplineMap.clear();
    fLoadedSplineSet.clear();
    fKinematicsMetadata.clear();
  }

  const int kNodeTypeStartElement = 1;
  const int kNodeTypeEndElement   = 15;
  const int kKnotX                = 0;
  const int kKnotY                = 1;

  xmlTextReaderPtr reader;

  int ret = 0, val_type = -1, iknot = 0, nknots = 0;
  double * E = 0, * xsec = 0;
  string spline_name = "";
  string temp_tune ;
  bool in_kinematics = false;

  reader = xmlNewTextReaderFilename(filename.c_str());
  if (reader != NULL) {
        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            xmlChar * name  = xmlTextReaderName     (reader);
            xmlChar * value = xmlTextReaderValue    (reader);
            int       type  = xmlTextReaderNodeType (reader);
            int       depth = xmlTextReaderDepth    (reader);

            if(depth==0 && type==kNodeTypeStartElement) {
               LOG("XSecSplLst", pDEBUG) << "Root element = " << name;
               if(xmlStrcmp(name, (const xmlChar *) "genie_xsec_spline_list")) {
                   LOG("XSecSplLst", pERROR)
                     << "\nXML doc. has invalid root element! [filename: " << filename << "]";
                   return kXmlInvalidRoot;
               }

               xmlChar * xvrs   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"version");
               xmlChar * xinlog = xmlTextReaderGetAttribute(reader,(const xmlChar*)"uselog");
               string svrs      = utils::str::TrimSpaces((const char *)xvrs);
               string sinlog    = utils::str::TrimSpaces((const char *)xinlog);

               LOG("XSecSplLst", pNOTICE)
                   << "Input x-section spline XML file format version: " << svrs;

               if (atoi(sinlog.c_str()) == 1) this->SetLogE(true);
               else this->SetLogE(false);

               xmlFree(xvrs);
               xmlFree(xinlog);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "genie_tune")) && type==kNodeTypeStartElement) {
               xmlChar * xtune = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
               temp_tune    = utils::str::TrimSpaces((const char *)xtune);
               SLOG("XSecSplLst", pNOTICE) << "Loading x-section splines for GENIE tune: " << temp_tune;
               xmlFree(xtune);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "kinematics")) && type==kNodeTypeStartElement) {
               in_kinematics = true;
            }

            if( in_kinematics && (!xmlStrcmp(name, (const xmlChar *) "q2")) &&
                type==kNodeTypeStartElement) {
               xmlChar * xmin    = xmlTextReaderGetAttribute(reader,(const xmlChar*)"min");
               xmlChar * xunit   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"unit");
               xmlChar * xsource = xmlTextReaderGetAttribute(reader,(const xmlChar*)"source");

               if(xmin) {
                 string smin = utils::str::TrimSpaces((const char *)xmin);
                 double q2min = atof(smin.c_str());
                 if(std::isfinite(q2min) && q2min > 0.) {
                   KinematicsMetadata & meta = fKinematicsMetadata[temp_tune];
                   meta.has_q2_min = true;
                   meta.q2_min = q2min;
                   meta.q2_unit = xunit ?
                     utils::str::TrimSpaces((const char *)xunit) : "GeV^2";
                   meta.q2_source = xsource ?
                     utils::str::TrimSpaces((const char *)xsource) : "";
                   SLOG("XSecSplLst", pNOTICE)
                     << "Loaded Q2 minimum metadata for tune " << temp_tune
                     << ": " << meta.q2_min << " " << meta.q2_unit;
                 } else {
                   SLOG("XSecSplLst", pWARN)
                     << "Ignoring invalid Q2 minimum metadata in " << filename
                     << ": " << smin;
                 }
               }

               if(xmin)    xmlFree(xmin);
               if(xunit)   xmlFree(xunit);
               if(xsource) xmlFree(xsource);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeStartElement) {
               xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
               xmlChar * xnkn  = xmlTextReaderGetAttribute(reader,(const xmlChar*)"nknots");
               string sname    = utils::str::TrimSpaces((const char *)xname);
               string snkn     = utils::str::TrimSpaces((const char *)xnkn);

               spline_name = sname;
               SLOG("XSecSplLst", pNOTICE) << "Loading spline: " << spline_name;

               nknots = atoi( snkn.c_str() );
               iknot=0;
               E     = new double[nknots];
               xsec  = new double[nknots];

               xmlFree(xname);
               xmlFree(xnkn);
            }

            if( (!xmlStrcmp(name, (const xmlChar *) "E"))    && type==kNodeTypeStartElement) { val_type = kKnotX; }
            if( (!xmlStrcmp(name, (const xmlChar *) "xsec")) && type==kNodeTypeStartElement) { val_type = kKnotY; }

            if( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==5) {
                if      (val_type==kKnotX) E   [iknot] = atof((const char *)value);
                else if (val_type==kKnotY) xsec[iknot] = atof((const char *)value);
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "knot")) && type==kNodeTypeEndElement) {
               iknot++;
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeEndElement) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
               LOG("XSecSplLst", pINFO) << "Done with current spline";
               for(int i=0; i<nknots; i++) {
                  LOG("XSecSplLst", pINFO) << "xsec[E = " << E[i] << "] = " << xsec[i];
               }
#endif
               // done looping over knots - build the spline
               Spline * spline = new Spline(nknots, E, xsec);
               delete [] E;
               delete [] xsec;

               // insert the spline to the map
               map<string,  map<string, Spline *> >::iterator //\/
               mm_iter = fSplineMap.find( temp_tune );
               if(mm_iter == fSplineMap.end()) {
                 map<string, Spline *> spl_map_curr_tune;
                 fSplineMap.insert( map<string, map<string, Spline *> >::value_type(
                    temp_tune, spl_map_curr_tune) );
                 mm_iter = fSplineMap.find( temp_tune );
               }
               map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
               spl_map_curr_tune.insert(
                  map<string, Spline *>::value_type(spline_name, spline) );
               if(mark_as_loaded) fLoadedSplineSet[temp_tune].insert(spline_name);
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "kinematics")) && type==kNodeTypeEndElement) {
               in_kinematics = false;
            }
            xmlFree(name);
            xmlFree(value);
            ret = xmlTextReaderRead(reader);
        }
        xmlFreeTextReader(reader);
        if (ret != 0) {
          LOG("XSecSplLst", pERROR)
            << "\nXML file could not be parsed! [filename: " << filename << "]";
          return kXmlNotParsed;
        }
  } else {
    LOG("XSecSplLst", pERROR)
          << "\nXML file could not be found! [filename: " << filename << "]";
  }

  return kXmlOK;
}
//____________________________________________________________________________
string XSecSplineList::BuildSplineKey(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  if(!alg) {
    LOG("XSecSplLst", pWARN)
            << "Null XSecAlgorithmI - Returning empty spline key";
    return "";
  }

  if(!interaction) {
    LOG("XSecSplLst", pWARN)
            << "Null Interaction - Returning empty spline key";
    return "";
  }

  string alg_name  = alg->Id().Name();
  string param_set = alg->Id().Config();
  string intkey    = interaction->AsString();

  string key = alg_name + "/" + param_set + "/" + intkey;

  return key;
}
//____________________________________________________________________________
const vector<string> * XSecSplineList::GetSplineKeys(void) const
{
  map<string,  map<string, Spline *> >::const_iterator //\/
  mm_iter = fSplineMap.find(fCurrentTune);
  if(mm_iter == fSplineMap.end()) {
    SLOG("XSecSplLst", pWARN)
       << "No splines for tune " << fCurrentTune << " were found!";
    return 0;
  }
  const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
  vector<string> * keyv = new vector<string>(spl_map_curr_tune.size());
  unsigned int i=0;
  map<string, Spline *>::const_iterator m_iter = spl_map_curr_tune.begin();
  for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
    string key = m_iter->first;
    (*keyv)[i++]=key;
  }
  return keyv;
}
//____________________________________________________________________________
void XSecSplineList::Print(ostream & stream) const
{
  stream << "\n ******************* XSecSplineList *************************";
  stream << "\n [-] Options:";
  stream << "\n  |";
  stream << "\n  |-----o  UseLogE..................." << fUseLogE;
  stream << "\n  |-----o  Spline Emin..............." << fEmin;
  stream << "\n  |-----o  Spline Emax..............." << fEmax;
  stream << "\n  |-----o  Spline NKnots............." << fNKnots;
  stream << "\n  |";

  map<string, map<string, Spline *> >::const_iterator mm_iter;
  for(mm_iter = fSplineMap.begin(); mm_iter != fSplineMap.end(); ++mm_iter) {

    string curr_tune = mm_iter->first;
    stream << "\n [-] Available x-section splines for tune: " << curr_tune ;
    stream << "\n  |";

    const map<string, Spline *> & spl_map_curr_tune = mm_iter->second;
    map<string, Spline *>::const_iterator m_iter = spl_map_curr_tune.begin();
    for( ; m_iter != spl_map_curr_tune.end(); ++m_iter) {
      string key = m_iter->first;
      stream << "\n  |-----o  " << key;
    }
    stream << "\n";
  }
}
//___________________________________________________________________________

} // genie namespace
