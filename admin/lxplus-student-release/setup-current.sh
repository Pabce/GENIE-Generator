#!/usr/bin/env bash
# Public entry point. Must be sourced.

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "ERROR: source this file instead of executing it:" >&2
  echo "  source /afs/cern.ch/work/p/pbarhama/public/genie/setup.sh" >&2
  exit 1
fi

_genie_student_root="/afs/cern.ch/work/p/pbarhama/public/genie"
_genie_student_release="$(readlink -f "${_genie_student_root}/current")"

if [[ -z "${_genie_student_release}" || ! -r "${_genie_student_release}/setup.sh" ]]; then
  echo "ERROR: no readable approved GENIE student release is published" >&2
  unset _genie_student_root _genie_student_release
  return 1
fi

# Resolve current once. The sourced versioned setup exports the physical path,
# so a later current-symlink update cannot alter this shell.
source "${_genie_student_release}/setup.sh"
_genie_student_setup_status=$?
unset _genie_student_root _genie_student_release
return "${_genie_student_setup_status}"

