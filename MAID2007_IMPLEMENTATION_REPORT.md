# MAID2007 Resonance Implementation Report

Date: 2026-07-09

Checkout: `/Users/pbarham/opt/genie-dev`

Implementation inventory source: commit `bfb82da02` (`Add MAID resonance configs and resonance-list override`) on branch `maid-dev`.

## Executive Summary

This implementation adds an opt-in MAID2007 electromagnetic resonance model to GENIE. The production configuration name is:

```text
genie::MAIDRESPXSec/Default
```

The model uses MAID2007 proton and neutron helicity amplitudes for the 13-state MAID set:

```text
P33(1232),P11(1440),D13(1520),S11(1535),S31(1620),S11(1650),
D15(1675),F15(1680),D33(1700),P13(1720),F35(1905),P31(1910),
F37(1950)
```

GENIE-only resonances remain in the standard `CommonParam[Resonances]` list, but are not part of `CommonParam[MAID2007Resonances]`. The model has an explicit RS fallback switch for missing MAID amplitudes, but the fallback is disabled by default.

The default MAID cross-section path is:

- `XSecRoute = DirectAS`
- MAID A/S helicity amplitudes
- MAID W2 Breit-Wigner convention
- no pi-N branching factor in the MAID Breit-Wigner numerator
- `UseDRJoinScheme = true`
- `Wcut = 1.7`
- `ApplyRSNativeWCut = false`
- `BreitWignerNorm = false`

The dedicated integrator is:

```text
genie::MAIDRESXSecFast/Default
```

It was added so MAID spline generation and event generation do not accidentally inherit the Rein-Sehgal cache warmup, resonance support truncation, or hard ESpline assumptions.

Two tune directories were added:

- `MAID07_00a`: a minimal RES-EM-only MAID tune for direct model checks.
- `G21_11m`: a G21_11a-style electromagnetic tune variant that swaps only the RES-EM cross-section model to MAID while leaving the event-generator module chain and resonance decay/final-state plumbing aligned with the normal RES path.

No new resonance decay code was needed. The implementation changes the RES-EM cross-section model, not the downstream RES event-generation and decay modules.

## Physics And Runtime Behavior

The MAID helicity implementation is table-driven from XML:

- `P33(1232)` uses the special MAID2007 M/E/S reduced-amplitude form.
- `P11(1440)` uses the Roper-specific MAID2007 A/S equations.
- Higher resonances use the MAID2007 `A0`, `alpha`, `beta`, and `S0` style coefficient tables where coefficients are published.

The `P31(1910)` resonance is included in the MAID2007 resonance selection and in the MAID resonance mass/width table, but no MAID A/S coefficient block is loaded for it. With `UseRSFallbackForMissingResonances=false`, it contributes zero in the MAID path. If `UseRSFallbackForMissingResonances=true`, it can be routed through the configured RS EM fallback.

For comparisons, the implementation also provides a native RS route that injects MAID A/S amplitudes through an A/S-to-RS-f bridge:

```text
genie::ReinSehgalRESPXSec/EM-NoPauliBlock-MAID2007F
```

This is for controlled RS/MAID comparisons only. The MAID default route is `genie::MAIDRESPXSec/Default`.

## File Inventory

### Core MAID Cross Section And Integrator

| File | Status | What it achieves |
| --- | --- | --- |
| `src/Physics/Resonance/XSection/MAIDRESPXSec.h` | Added | Declares the MAID2007 RES-EM single-resonance differential cross-section model, including MAID resonance parameters, Breit-Wigner controls, W-limit controls, selected-resonance filtering, and optional RS fallback state. |
| `src/Physics/Resonance/XSection/MAIDRESPXSec.cxx` | Added | Implements the MAID RES-EM cross section. It selects proton/neutron MAID helicity amplitudes, evaluates `DirectAS` and comparison `RSFBridge` routes, applies the MAID/RS Breit-Wigner variants, loads the 13-state MAID resonance mass/width/branching/vertex table, enforces the configurable `Wcut`, and delegates total integration to the configured integrator. |
| `src/Physics/Resonance/XSection/MAIDRESXSecFast.h` | Added | Declares the dedicated MAID fast RES integrator. It mirrors the useful W,Q2 integration transformation from the fast RES path while keeping MAID-specific cache and support decisions local to MAID. |
| `src/Physics/Resonance/XSection/MAIDRESXSecFast.cxx` | Added | Implements MAID total cross-section integration and optional local free-nucleon excitation caches. The cache key includes the MAID resonance list, helicity algorithms, Breit-Wigner mode, Wcut, fallback flag, and related model controls so cached splines track the actual MAID configuration. |

### MAID Helicity Amplitudes

| File | Status | What it achieves |
| --- | --- | --- |
| `src/Physics/Resonance/XSection/MAIDHelicityAmpl.h` | Added | Declares a compact container for MAID A/S amplitudes and the corresponding RS-style f amplitudes used in bridge comparisons. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmpl.cxx` | Added | Implements reset, assignment, validity, and printing for the MAID helicity-amplitude container. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelI.h` | Added | Defines the abstract MAID helicity-amplitude model interface used by `MAIDRESPXSec`. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelI.cxx` | Added | Implements constructors/destructor for the MAID helicity-amplitude model interface. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.h` | Added | Declares the proton MAID2007 helicity-amplitude model, coefficient containers, P33/P11/higher-resonance evaluation helpers, RS-original comparison modes, and configuration loading. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelEMp.cxx` | Added | Implements the proton MAID2007 amplitude evaluator. It supports the special P33 equations, the Roper equations, table-driven higher-resonance fits, optional RS-original comparison amplitudes, and conversion to RS-style f amplitudes where needed. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelEMn.h` | Added | Declares the neutron MAID2007 helicity-amplitude model as the neutron-specialized subclass of the proton implementation. |
| `src/Physics/Resonance/XSection/MAIDHelicityAmplModelEMn.cxx` | Added | Enables neutron coefficient loading/evaluation by setting the neutron-amplitude mode on the shared implementation. |

### MAID To Native RS Bridge

| File | Status | What it achieves |
| --- | --- | --- |
| `src/Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMp.h` | Added | Declares the proton bridge model that exposes MAID A/S amplitudes through the native `RSHelicityAmplModelI` interface. |
| `src/Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMp.cxx` | Added | Implements the A/S-to-f mapping for native RS comparisons. It receives external W, Q2, and nucleon mass through `SetKinematics`, asks the MAID model for A/S, and fills an `RSHelicityAmpl` object. |
| `src/Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMn.h` | Added | Declares the neutron bridge model. |
| `src/Physics/Resonance/XSection/MAIDToRSHelicityAmplModelEMn.cxx` | Added | Provides the neutron bridge subclass and routes it through the shared proton bridge implementation with the neutron MAID model configured. |
| `src/Physics/Resonance/XSection/RSHelicityAmpl.h` | Modified | Adds a public `Set(...)` method so a bridge model can directly populate RS f amplitudes without pretending they came from the original FKR calculation. |
| `src/Physics/Resonance/XSection/RSHelicityAmpl.cxx` | Modified | Implements the new `RSHelicityAmpl::Set(...)` assignment helper. |
| `src/Physics/Resonance/XSection/RSHelicityAmplModelI.h` | Modified | Adds an optional `SetKinematics(W,Q2,M)` hook. Default RS models ignore it, while the MAID bridge uses it because EM Q2 and W are not contained in the older FKR-only interface. |
| `src/Physics/Resonance/XSection/ReinSehgalRESPXSec.cxx` | Modified | Calls `SetKinematics` before computing helicity amplitudes and loads CC/NC/EM helicity-amplitude sub-algorithms from XML. This makes the native RS model configurable enough to run with MAID-derived EM f amplitudes for comparisons. |

### Source Registration

| File | Status | What it achieves |
| --- | --- | --- |
| `src/Physics/Resonance/XSection/LinkDef.h` | Modified | Registers MAID cross-section, integrator, helicity-amplitude, and bridge classes with the ROOT dictionary. |
| `config/master_config.xml` | Modified | Registers the new MAID and MAID-to-RS XML configuration files with GENIE's algorithm configuration system. |

### MAID XML Configuration

| File | Status | What it achieves |
| --- | --- | --- |
| `config/MAIDRESPXSec.xml` | Added | Defines `genie::MAIDRESPXSec` configurations. `Default` is the production MAID2007 route; additional parameter sets provide comparison modes such as RS-original amplitudes, P33-only bridge checks, MAID2007 running-width variants, RS-original Breit-Wigner variants, and explicit RS fallback. |
| `config/MAIDRESXSecFast.xml` | Added | Configures the dedicated MAID integrator. The default uses the MAID resonance list, no Pauli blocking, no local cache by default, 50 optional cache knots, and adaptive GSL integration settings. |
| `config/MAIDHelicityAmplModelEMp.xml` | Added | Stores proton MAID2007 helicity-amplitude coefficients and comparison-mode controls. It contains P33 reduced-amplitude parameters, P11 Roper parameters, and higher-resonance A/S fit coefficients. |
| `config/MAIDHelicityAmplModelEMn.xml` | Added | Stores the neutron MAID2007 helicity-amplitude coefficients and comparison-mode controls. |
| `config/MAIDToRSHelicityAmplModelEMp.xml` | Added | Configures the proton MAID A/S-to-RS-f bridge and its MAID source algorithm. |
| `config/MAIDToRSHelicityAmplModelEMn.xml` | Added | Configures the neutron MAID A/S-to-RS-f bridge and its MAID source algorithm. |
| `config/CommonParam.xml` | Modified | Adds the shared `CommonParam[MAID2007Resonances]` list containing the 13-state MAID set, while leaving GENIE's larger `CommonParam[Resonances]` list intact. |
| `config/RSPPResonanceSelector.xml` | Modified | Adds a `MAID2007` selector parameter set that reads `CommonParam[MAID2007Resonances]`. |
| `config/ReinSehgalRESPXSec.xml` | Modified | Adds configurable helicity-amplitude algorithm keys and the `EM-NoPauliBlock-MAID2007F` native-RS comparison configuration. |
| `config/ReinSehgalRESXSec.xml` | Modified | Adds `NoPauliBlock-MAID2007`, an RS integrator configuration restricted to `CommonParam[MAID2007Resonances]` for bridge/comparison studies. |
| `config/ReinSehgalRESXSecFast.xml` | Modified | Adds the matching fast-integrator `NoPauliBlock-MAID2007` configuration for native-RS comparisons on the MAID subset. |
| `config/ReinSehgalSPPPXSec.xml` | Modified | Adds a `MAID2007` parameter set so single-pion RES comparisons can also restrict to the MAID resonance subset. |

### Tune Directories

| File | Status | What it achieves |
| --- | --- | --- |
| `config/MAID07_00a/ModelConfiguration.xml` | Added | Defines the minimal MAID tune's physics model override: `XSecModel@genie::EventGenerator/RES-EM = genie::MAIDRESPXSec/Default`. |
| `config/MAID07_00a/TuneGeneratorList.xml` | Added | Makes the minimal MAID tune run only `genie::EventGenerator/RES-EM`. |
| `config/MAID07_00a/CommonPhaseSpaceCuts.xml` | Added | Sets the EM phase-space cut used by the minimal MAID tune, including `EM-Q2-min = 0.02`. |
| `config/G21_11m/ModelConfiguration.xml` | Added | Creates the G21 EM variant and swaps only the `RES-EM` cross-section model to `genie::MAIDRESPXSec/Default`, leaving the rest of the G21 EM model choices in place. |
| `config/G21_11m/TuneGeneratorList.xml` | Added | Defines the G21_11m default EM generator list as QEL-EM, MEC-EM, RES-EM, and DIS-EM. |
| `config/G21_11m/EventGenerator.xml` | Added | Provides the event-generator module definitions for the G21_11m overlay. This keeps resonance generation, decays, transport, and final-state handling aligned with the standard GENIE RES thread rather than introducing MAID-specific decay code. |
| `config/G21_11m/CommonParam.xml` | Added | Provides the tune-local common-parameter file for G21_11m. Because tune-local `CommonParam.xml` shadows the base file, this includes the MAID2007 resonance list needed by MAID configs in this tune. |
| `config/G21_11m/CommonPhaseSpaceCuts.xml` | Added | Sets the G21_11m EM phase-space cut, including `EM-Q2-min = 0.02`. |
| `config/G21_11m/MECInteractionListGenerator.xml` | Added | Carries the G21 EM MEC interaction-list settings needed by the tune overlay, including EM MEC mode and dinucleon-code choices. |

### Command-Line And Framework Plumbing

| File | Status | What it achieves |
| --- | --- | --- |
| `src/Framework/Utils/AppInit.h` | Modified | Declares the shared resonance-list override helper used by GENIE apps. |
| `src/Framework/Utils/AppInit.cxx` | Modified | Implements `utils::app_init::ResonanceNameList(...)`. It parses comma-separated resonance names, validates them with GENIE's resonance-name conversion, canonicalizes names, and overrides both `CommonParam[Resonances]` and `CommonParam[MAID2007Resonances]` for the current run. |
| `src/Apps/gEvGen.cxx` | Modified | Adds the `--resonances` command-line option to `gevgen` and applies the override after tune construction. This enables one-resonance or selected-resonance event generation. |
| `src/Apps/gMakeSplines.cxx` | Modified | Adds the same `--resonances` option to `gmkspl`, so spline generation can be restricted to one or more selected resonances consistently with `gevgen`. |
| `src/Framework/EventGen/PhysInteractionSelector.cxx` | Modified | Clamps negative cross-section values by assigning `xsec = TMath::Max(0., xsec)`. This protects interaction selection from small negative interpolation/integration artifacts. |
| `src/Framework/Numerical/Spline.cxx` | Modified | Hardens closest-knot lookup at spline boundaries and invalid knot indices. This prevents boundary pathologies when spline values are queried near or outside the tabulated support. |

### Tests

| File | Status | What it achieves |
| --- | --- | --- |
| `src/contrib/test/Makefile` | Modified | Adds the `gtestMAID2007Helicity` test target to the contrib test build, install, and clean lists. |
| `src/contrib/test/gtestMAID2007Helicity.cxx` | Added | Adds focused MAID2007 helicity/config checks. The test covers representative MAID amplitude behavior and asserts important default configuration choices, including `BreitWignerNorm=false` and `ApplyRSNativeWCut=false` for `MAIDRESPXSec/Default`. |

## Configuration Recipes

To configure an existing tune to use MAID for electromagnetic resonance production, modify the tune's `ModelConfiguration.xml` rather than the RES event-generator module chain:

```xml
<param type="alg" name="XSecModel@genie::EventGenerator/RES-EM">
  genie::MAIDRESPXSec/Default
</param>
```

Keep the normal `genie::EventGenerator/RES-EM` thread in `TuneGeneratorList.xml` or the selected `--event-generator-list`. That is what preserves the usual RES final-state and decay wiring.

To run only selected resonances:

```bash
gevgen ... --resonances 'P33(1232),P11(1440)'
gmkspl ... --resonances 'P33(1232),P11(1440)'
```

The override updates both the standard GENIE resonance list and the MAID2007 resonance list for that process.

## Validation Notes

Development validation for this implementation included:

- XML registration and config-resolution checks for `MAIDRESPXSec/Default`, `MAID07_00a`, and `G21_11m`.
- Local resonance cross-section comparisons between MAID DirectAS, MAID/RS Breit-Wigner variants, native RS, and native RS with MAID-derived f amplitudes.
- Local spline generation up to 10 GeV with the dedicated MAID integrator.
- Event-generation checks against direct cross-section curves after correcting the histogram normalization by the W factor.
- Focused build/test checks for the resonance xsec code and `gtestMAID2007Helicity`.

Generated artifacts from those studies, such as ROOT event files, spline XML outputs, and plots, are not part of this implementation inventory.

## Documentation Added By This Request

| File | Status | What it achieves |
| --- | --- | --- |
| `docs/handoffs/maid2007-implementation-report.md` | Added | This report. It records the implementation shape and file-by-file purpose so the MAID changes can be reviewed without relying on the IDE diff sidebar or prior chat context. |
