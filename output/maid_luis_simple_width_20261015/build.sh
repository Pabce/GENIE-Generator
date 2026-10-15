#!/usr/bin/env bash
set -euo pipefail
export GENIE=/Users/pbarham/opt/genie-dev
export GXMLPATH="$GENIE/config"
export PATH="$GENIE/bin:$PATH"
export DYLD_LIBRARY_PATH="$GENIE/lib:/Users/pbarham/opt/ROOTEGPythia6/lib:$(root-config --libdir):${DYLD_LIBRARY_PATH:-}"
ROOT_FLAGS="$(root-config --cflags --glibs) -lMinuit -lGeom -lEG -lGenVector -lMathMore"
GENIE_LIBS="$("$GENIE/bin/genie-config" --libs | sed 's/-lGPhHadTransp//g')"
EXTRA_FLAGS="$(lhapdf-config --cflags --libs) $(xml2-config --cflags --libs) $(gsl-config --cflags --libs) $(pkg-config --cflags --libs log4cpp)"
clang++ -O2 -I"$GENIE/src" $ROOT_FLAGS $EXTRA_FLAGS \
  "$GENIE/validation/maid_luis_comparison/scan_before_after.cxx" \
  -o "$GENIE/output/maid_luis_simple_width_20261015/scan_before_after" \
  $GENIE_LIBS -L/Users/pbarham/opt/ROOTEGPythia6/lib -lEGPythia6 \
  -Wl,-rpath,"$GENIE/lib" -Wl,-rpath,/Users/pbarham/opt/ROOTEGPythia6/lib \
  -Wl,-rpath,"$(root-config --libdir)"
