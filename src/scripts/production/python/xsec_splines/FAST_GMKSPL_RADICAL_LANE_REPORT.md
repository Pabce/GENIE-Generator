# Fast gmkspl Radical Lane Report

## Status

The shared worktree already contains uncommitted edits in:

- `src/Apps/gMakeSplines.cxx`
- `src/Framework/Utils/XSecSplineList.cxx`
- `src/Framework/Utils/XSecSplineList.h`

After the coordinator checkpoint, I did not continue editing those files. I had
briefly added sparse-mode declarations to `XSecSplineList.h`, then removed only
that partial change so the shared C++ edits were left as they were.

New visible files added by this lane:

- `src/scripts/production/python/xsec_splines/compare_spline_xml_fastlane.py`
- `src/scripts/production/python/xsec_splines/FAST_GMKSPL_RADICAL_LANE_REPORT.md`

## Radical Acceleration Route

The highest-upside route is a combined opt-in mode:

1. Fan out independent probe-target initial states into worker processes.
2. Sparse-evaluate each spline with a small number of exact
   `alg->Integral(interaction)` knots.
3. Interpolate onto the usual dense energy grid and emit normal GENIE XML.
4. Evaluate guard knots exactly; if a spline exceeds tolerance, fall back to
   exact construction for that spline.

Suggested opt-in CLI:

```bash
gmkspl ... --q2-min 0.02 --jobs 8 \
  --fast-splines --fast-spline-knots 10 \
  --fast-spline-checks 5 --fast-spline-tol 0.15
```

Expected speedup:

- Sparse mode alone: roughly `n_default_knots / (fast_knots + checks)`, so
  around `2x` for 30-knot splines and higher for denser splines.
- Worker fan-out: near-linear over independent probe-target states until memory
  or disk contention dominates.
- Combined aggressive mode can plausibly reach `10x-30x` wall-clock reduction
  for multi-target neutrino campaigns if guard fallback rates stay low.

## Build Risks In Current Shared C++ State

The current shared C++ edits were updated by the shared lane while this report
was being written. On the latest symbol check, the worker functions and option
plumbing are present:

- `gMakeSplines.cxx` declares `GetInitialStates`, `BuildSplinesForInitialState`,
  `BuildSplinesSerial`, `BuildSplinesParallel`, and `ParallelWorkerOutputName`,
  and now appears to define them.
- `--jobs` and `GENIE_GMKSPL_JOBS` appear wired into parsing and logging.
- `--fast-spline-eval` appears wired to `XSecSplineList::SetFastSplineEval`.
- `XSecSplineList` has a new `LoadFromXml(..., mark_as_loaded)` signature for
  worker XML import.

Remaining build/runtime risks:

- The shared C++ changes were not compile-checked by this lane after the
  checkpoint.
- `--fast-spline-eval` is currently a low-risk exact-mode optimization: it skips
  below-threshold integrations and summarizes timing. It is not the sparse
  surrogate mode proposed above.
- Forking after GENIE/ROOT initialization may expose macOS-specific singleton or
  buffered-logging issues.
- Worker output merge needs an integration test with `--input-cross-sections`
  and `--no-copy`, because loaded/generated spline bookkeeping changed.

## Validation Plan

Generate an exact reference and candidate fast spline with the same probe,
target, tune, q2 cut, energy max, and nominal knot count:

```bash
gmkspl -p 14 -t 1000180400 -n 30 -e 20 \
  --tune G18_10a_02_11a --event-generator-list Default \
  --q2-min 0.02 -o exact.xml

gmkspl -p 14 -t 1000180400 -n 30 -e 20 \
  --tune G18_10a_02_11a --event-generator-list Default \
  --q2-min 0.02 -o fast.xml \
  --jobs 4 --fast-splines --fast-spline-knots 10 \
  --fast-spline-checks 5 --fast-spline-tol 0.15
```

Compare with the new dependency-free harness:

```bash
python3 src/scripts/production/python/xsec_splines/compare_spline_xml_fastlane.py \
  exact.xml fast.xml --point-rel-warn 0.25 --area-rel-warn 0.05 \
  --csv fast_vs_exact.csv
```

Initial acceptance gates:

- No missing candidate splines.
- Trapezoid-area relative difference below `5%` for every spline.
- Pointwise relative difference below `25%`, with manual review for isolated
  near-zero threshold points.
- Any physics-important channel failing guard checks should automatically fall
  back to exact construction.

## Failure Modes

- Sparse interpolation can miss narrow threshold structure, resonant turn-ons,
  or channel-specific non-monotonic features.
- Cubic interpolation over sparse exact knots can overshoot near zeros; dense
  emission should clamp negative approximations to zero before XML writing.
- Forked workers inherit GENIE and ROOT singleton state; this should remain
  opt-in and needs lxplus/macOS smoke tests separately.
- Worker XML merges must avoid marking generated worker splines as preloaded,
  otherwise `--no-copy` can accidentally drop new splines.
- Guard knots may not catch area drift between knots; always run the XML
  comparator before trusting production fast splines.
