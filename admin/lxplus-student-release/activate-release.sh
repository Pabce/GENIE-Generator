#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 vYYYY.N" >&2
  exit 2
fi

version="$1"
root="/afs/cern.ch/work/p/pbarhama/public/genie"
release="${root}/releases/${version}"

[[ -d "${release}" ]]
[[ -r "${release}/RELEASE.json" ]]
"${root}/source/GENIE-Student/admin/lxplus-student-release/validate-release.sh" "${release}"

temporary="${root}/.current.${version}.$$"
ln -s "releases/${version}" "${temporary}"
mv -Tf "${temporary}" "${root}/current"

echo "current=${release}"

