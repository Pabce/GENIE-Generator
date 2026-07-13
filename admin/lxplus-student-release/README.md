# LXPLUS student release

This directory contains the version-controlled tooling used to build and
publish the read-only GENIE installation for students on LXPLUS.

## Public paths

- Source repository: `/afs/cern.ch/work/p/pbarhama/public/genie/source/GENIE-Student`
- Versioned installs: `/afs/cern.ch/work/p/pbarhama/public/genie/releases/<version>`
- Approved release: `/afs/cern.ch/work/p/pbarhama/public/genie/current`
- Student entry point: `/afs/cern.ch/work/p/pbarhama/public/genie/setup.sh`

Students should use a versioned setup path when exact reproducibility matters:

```bash
source /afs/cern.ch/work/p/pbarhama/public/genie/releases/v2026.1/setup.sh
```

To modify GENIE, students should clone the public source into their own work
area rather than editing the shared read-only checkout:

```bash
git clone /afs/cern.ch/work/p/pbarhama/public/genie/source/GENIE-Student \
  "$HOME/GENIE-Student"
```

## Release policy

Each release is an immutable directory built from an annotated
`student-vYYYY.N` tag on the `student/stable` branch. Never rebuild or edit an
existing release directory. Build and validate a new version, then change the
`current` symlink. Rollback is the inverse symlink change.

`build-release.sh <version> <git-ref>` builds a new candidate and refuses to
overwrite an existing release. `activate-release.sh <version>` changes the
approved release only after the versioned installation validates.
