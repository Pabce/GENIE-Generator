#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 vYYYY.N git-ref" >&2
  exit 2
fi

version="$1"
source_ref="$2"
work="/afs/cern.ch/work/p/pbarhama"
public_root="${work}/public/genie"
source_repo="${public_root}/source/GENIE-Student"
release="${public_root}/releases/${version}"
build_dir="${work}/build/genie-student-${version}"
tooling="${source_repo}/admin/lxplus-student-release"
lcg_view="/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh"
private_p6="${work}/sw/genie_p6/pythia6"
private_tp6="${work}/sw/genie_p6/TPythia6"
jobs="${JOBS:-8}"

if [[ -e "${release}" || -e "${build_dir}" ]]; then
  echo "ERROR: refusing to overwrite an existing release or build directory" >&2
  echo "release=${release}" >&2
  echo "build=${build_dir}" >&2
  exit 1
fi

git -C "${source_repo}" rev-parse --verify "${source_ref}^{commit}" >/dev/null
mkdir -p "${public_root}/releases" "${work}/build" "${release}"
fs setacl -dir "${release}" -acl system:anyuser none

mkdir -p \
  "${release}/external/pythia6/lib" \
  "${release}/external/tpythia6/include" \
  "${release}/external/tpythia6/lib"

install -m 0755 "${private_p6}/lib/libPythia6.so" "${release}/external/pythia6/lib/"
ln -s libPythia6.so "${release}/external/pythia6/lib/libpythia6.so"
install -m 0755 "${private_tp6}/build/libEGPythia6.so" "${release}/external/tpythia6/lib/"
install -m 0644 \
  "${private_tp6}/build/libEGPythia6.rootmap" \
  "${private_tp6}/build/libEGPythia6_rdict.pcm" \
  "${private_tp6}/build/module.modulemap" \
  "${release}/external/tpythia6/lib/"
cp -a "${private_tp6}/inc/." "${release}/external/tpythia6/include/"

git clone --no-local --no-checkout "${source_repo}" "${build_dir}"
git -C "${build_dir}" checkout --detach "${source_ref}"

unset ROOTSYS GXMLPATH
set +u
source "${lcg_view}"
set -u

export GENIE="${build_dir}"
export PYTHIA6="${release}/external/pythia6"
export TP6="${release}/external/tpythia6"
export CPATH="${TP6}/include${CPATH:+:${CPATH}}"
export ROOT_INCLUDE_PATH="${TP6}/include${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}"
export LIBRARY_PATH="${TP6}/lib:${PYTHIA6}/lib${LIBRARY_PATH:+:${LIBRARY_PATH}}"
export LD_LIBRARY_PATH="${TP6}/lib:${PYTHIA6}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

pythia8_inc="$(pythia8-config --includedir)"
pythia8_lib="$(pythia8-config --libdir)"
lcg_prefix="$(dirname "${lcg_view}")"
lhapdf6_inc="${lcg_prefix}/include"
lhapdf6_lib="${lcg_prefix}/lib"
log4cpp_inc="${lcg_prefix}/include"
log4cpp_lib="${lcg_prefix}/lib"

configure_args=(
  "--prefix=${release}"
  --enable-pythia6
  --enable-pythia8
  --disable-lhapdf5
  --enable-lhapdf6
  --enable-flux-drivers
  --enable-geom-drivers
  --enable-dylibversion
  "--with-pythia6-lib=${PYTHIA6}/lib"
  "--with-pythia8-inc=${pythia8_inc}"
  "--with-pythia8-lib=${pythia8_lib}"
  "--with-lhapdf6-inc=${lhapdf6_inc}"
  "--with-lhapdf6-lib=${lhapdf6_lib}"
  "--with-log4cpp-inc=${log4cpp_inc}"
  "--with-log4cpp-lib=${log4cpp_lib}"
  --with-libxml2-inc=/usr/include/libxml2
  --with-libxml2-lib=/usr/lib64
)

(
  cd "${build_dir}"
  ./configure "${configure_args[@]}"
  make -j "${jobs}"
  make install
  mkdir -p "${release}/src/make"
  install -m 0644 \
    src/make/Make.config \
    src/make/Make.config_no_paths \
    "${release}/src/make/"
) >"${release}/build.log" 2>&1

sed \
  -e "s|@RELEASE_VERSION@|${version}|g" \
  -e "s|@RELEASE_PREFIX@|${release}|g" \
  -e "s|@LCG_VIEW@|${lcg_view}|g" \
  "${tooling}/setup-release.sh.in" >"${release}/setup.sh"
chmod 0755 "${release}/setup.sh"
install -m 0755 "${tooling}/validate-release.sh" "${release}/validate-release.sh"

commit="$(git -C "${build_dir}" rev-parse HEAD)"
build_time="$(date --utc +%Y-%m-%dT%H:%M:%SZ)"
root_version="$(root-config --version)"
compiler="$(g++ --version | head -n 1)"
features="$("${release}/bin/genie-config" --features 2>/dev/null || true)"
configure_command="./configure ${configure_args[*]}"

python3 - "${release}/RELEASE.json" <<PY
import json
import sys

path = sys.argv[1]
manifest = {
    "schema": "genie-student-release-v1",
    "version": ${version@Q},
    "git_ref": ${source_ref@Q},
    "git_commit": ${commit@Q},
    "source_repository": ${source_repo@Q},
    "installation_prefix": ${release@Q},
    "built_at_utc": ${build_time@Q},
    "lcg_view": ${lcg_view@Q},
    "root_version": ${root_version@Q},
    "compiler": ${compiler@Q},
    "features": ${features@Q},
    "configure_command": ${configure_command@Q},
}
with open(path, "w", encoding="utf-8") as stream:
    json.dump(manifest, stream, indent=2, sort_keys=True)
    stream.write("\n")
PY

"${tooling}/validate-release.sh" "${release}" | tee "${release}/validation.log"

(
  cd "${release}"
  find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum >SHA256SUMS
)

chmod -R u+rwX,go+rX,go-w "${release}"
while IFS= read -r -d '' directory; do
  fs setacl -dir "${directory}" -acl system:anyuser rl
done < <(find "${release}" -type d -print0)

echo "built=${release}"
echo "source_ref=${source_ref}"
echo "commit=${commit}"
