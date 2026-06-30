#!/usr/bin/env bash
# 30k Zig-vs-C++ partition parity gate (issue #375).
#
# Runs an unbounded paired-loci partition with the Zig and C++ backends on the
# same 30k input, then compares cluster membership and per-event annotation
# fields. Exits 0 if identical, 1 otherwise. The 5k partis-test.py rung doesn't
# exercise enough cache round-trips to surface tiebreak-density-dependent
# divergences (#375 was latent until 30k).
#
# Usage:
#   bin/partis-30k-parity-test.sh [INFILE [OUTBASE [N]]]
#
# Defaults:
#   INFILE = /tmp/paired-paper-all-seqs.fa
#   OUTBASE = /tmp/parity-30k-test
#   N = 30000
#
# Wall budget on pax (32-core Zen 5):
#   ~50 min Zig + ~100 min C++ + a few min compare = ~2.5 h total when sequential.
#   Backends run sequentially to keep the host responsive; bump --n-procs in
#   this script to parallelize within a backend.

set -euo pipefail

INFILE="${1:-/tmp/paired-paper-all-seqs.fa}"
OUTBASE="${2:-/tmp/parity-30k-test}"
N="${3:-30000}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Match bin/zig-perf-counters-run.sh: per project CLAUDE.md, activate the venv
# before invoking partis so an unactivated shell doesn't fall back to system
# python and fail with an import error.
# shellcheck disable=SC1091
source "${REPO_ROOT}/.venv/bin/activate"

COMPARE="${REPO_ROOT}/bin/zig-compare.py"
if [[ ! -x "${COMPARE}" ]]; then
  echo "ERROR: cannot find zig-compare.py at ${COMPARE}" >&2
  exit 2
fi

if [[ ! -f "${INFILE}" ]]; then
  echo "ERROR: input file ${INFILE} does not exist" >&2
  exit 2
fi

CPP_DIR="${OUTBASE}-cpp"
ZIG_DIR="${OUTBASE}-zig"

run_backend() {
  local backend="$1" outdir="$2"
  shift 2
  rm -rf "${outdir}"
  echo ">> running ${backend} backend → ${outdir}"
  PYTHONHASHSEED=0 python3 "${REPO_ROOT}/bin/partis" partition --paired-loci \
    --paired-outdir "${outdir}" \
    --infname "${INFILE}" --n-max-queries "${N}" \
    --random-seed 1 --n-procs 10 --no-time-based-n-proc-reduction \
    --dont-write-git-info "$@"
}

# C++ first (slower), then Zig — sequential keeps the host responsive at the
# usual 10 procs. Either order works; the partition is deterministic.
run_backend "C++"  "${CPP_DIR}"
run_backend "Zig"  "${ZIG_DIR}"  --zig

echo ">> comparing"
"${COMPARE}" "${CPP_DIR}" "${ZIG_DIR}"
