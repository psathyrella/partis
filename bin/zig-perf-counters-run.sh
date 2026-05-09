#!/usr/bin/env bash
# Build partis-zig-core with -Dperf-counters=true and run a partition workload
# with PARTIS_ZIG_PERF_LOG set, so each query emits one PERFCOUNTER line plus
# one PERFCOUNTER_CACHE line (chunk-cache prefix-list-length histogram) into
# the log. See packages/zig-core/src/ham/perf_counters.zig for the metric set
# and packages/zig-core/PERF-NEXT.md Phase 0.1/0.2 for what these are for.
#
# Pair with bin/zig-perf-counters-aggregate.py to compute median/p95/max
# distributions and chunk-cache hit rate from the resulting log.
#
# Usage:
#   bin/zig-perf-counters-run.sh [INFILE [OUTBASE [N [N_PROCS]]]]
#
# Defaults:
#   INFILE  = /tmp/paired-paper-all-seqs.fa
#   OUTBASE = /tmp/perf-counters-test
#   N       = 1000
#   N_PROCS = 1   (PERF-NEXT.md 0.1 specifies n-procs 1 for clean wall-time
#                 numbers; for distribution stats only, n_procs > 1 is fine)
#
# After the run, the perf log path is:
#   $OUTBASE/zig-perf.log
# and the aggregated summary path is:
#   $OUTBASE/aggregate.txt
#
# The script restores the un-instrumented ReleaseFast build at the end so the
# working tree's binary is back to bit-equality-baseline state.

set -euo pipefail

INFILE="${1:-/tmp/paired-paper-all-seqs.fa}"
OUTBASE="${2:-/tmp/perf-counters-test}"
N="${3:-1000}"
N_PROCS="${4:-1}"

PARTIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ZIG_CORE_DIR="$PARTIS_DIR/packages/zig-core"
PERF_LOG="$OUTBASE/zig-perf.log"
AGG_OUT="$OUTBASE/aggregate.txt"

mkdir -p "$OUTBASE"
rm -f "$PERF_LOG" "$AGG_OUT"

# ── 1. Build with perf counters enabled ─────────────────────────────────────
echo "==> building partis-zig-core with -Dperf-counters=true"
(cd "$ZIG_CORE_DIR" && zig build -Doptimize=ReleaseFast -Dperf-counters=true)

# ── 2. Run partition with the perf log set ──────────────────────────────────
RUN_OUTBASE="$OUTBASE/run-output"
rm -rf "$RUN_OUTBASE"
mkdir -p "$RUN_OUTBASE"

echo "==> running partition: N=$N, n_procs=$N_PROCS, infile=$INFILE"
echo "==> perf log: $PERF_LOG"

# shellcheck disable=SC1091
source "$PARTIS_DIR/.venv/bin/activate"

PARTIS_ZIG_PERF_LOG="$PERF_LOG" PYTHONHASHSEED=0 python3 "$PARTIS_DIR/bin/partis" \
    partition --paired-loci \
    --paired-outdir "$RUN_OUTBASE" \
    --infname "$INFILE" \
    --n-max-queries "$N" \
    --random-seed 1 \
    --n-procs "$N_PROCS" \
    --no-time-based-n-proc-reduction \
    --dont-write-git-info \
    --zig

# ── 3. Aggregate ────────────────────────────────────────────────────────────
echo "==> aggregating $PERF_LOG → $AGG_OUT"
python3 "$PARTIS_DIR/bin/zig-perf-counters-aggregate.py" "$PERF_LOG" | tee "$AGG_OUT"

# ── 4. Restore un-instrumented build ────────────────────────────────────────
echo "==> restoring un-instrumented partis-zig-core (ReleaseFast, no perf counters)"
(cd "$ZIG_CORE_DIR" && zig build -Doptimize=ReleaseFast)

echo "==> done. raw log: $PERF_LOG    summary: $AGG_OUT"
