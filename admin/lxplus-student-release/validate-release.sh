#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 /absolute/release/prefix" >&2
  exit 2
fi

release="$(readlink -f "$1")"
[[ -r "${release}/setup.sh" ]]

source "${release}/setup.sh" >/dev/null

[[ "${GENIE}" == "${release}" ]]
[[ "$(command -v gevgen)" == "${release}/bin/gevgen" ]]
[[ "$(command -v gmkspl)" == "${release}/bin/gmkspl" ]]
[[ "$(command -v genie-config)" == "${release}/bin/genie-config" ]]
[[ -r "${release}/config/Messenger.xml" ]]
[[ -d "${release}/data" ]]

genie-config --features

for executable in gevgen gmkspl gntpc gconfigdump; do
  unresolved="$(ldd "${release}/bin/${executable}" 2>&1 | grep -E 'not found|/src/GENIE-Generator|/sw/genie_p6' || true)"
  if [[ -n "${unresolved}" ]]; then
    echo "ERROR: non-release dependency found for ${executable}:" >&2
    echo "${unresolved}" >&2
    exit 1
  fi
done

gconfigdump -h >/dev/null 2>&1 || true
echo "validated=${release}"
