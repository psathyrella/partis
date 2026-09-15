#!/usr/bin/env bash
# Nightly Debug-mode GPA leak check for packages/zig-core (issue #405).
#
# Invoked by a systemd user timer (see deploy/nightly-ci/). Fetches main into a
# dedicated clone, builds packages/zig-core in Debug mode -- the only build
# mode where main.zig routes through the leak-tracking GeneralPurposeAllocator
# (GPA) instead of the C allocator; see issue #405 for why ReleaseSafe and
# ReleaseFast don't catch leaks in this codebase -- runs `zig build test` plus
# the full non-quick, non-paired partis-test.py tier against Debug binaries,
# and posts to ntfy if either step fails or a leak is detected.
#
# Environment overrides (all optional):
#   PARTIS_ZIG_NIGHTLY_CLONE   clone directory        (default ~/re/pz-zig-nightly-ci)
#   PARTIS_ZIG_NIGHTLY_LOGDIR  log directory          (default ~/re/partis-zig-nightly-logs)
#   PARTIS_ZIG_NIGHTLY_TOPIC   ntfy topic             (default bip-matsen-partis-zig-nightly)
#   PARTIS_ZIG_NIGHTLY_RETAIN  log retention days     (default 14)
#   PARTIS_ZIG_NIGHTLY_REF     git ref to reset to    (default origin/main)
#
# Exit code reflects the outcome: 0 on success (no failures, no leaks
# detected), non-zero on any build/test failure or on a detected leak.
# Systemd records this via ExecStart.
#
# See deploy/nightly-ci/README.md for one-time clone setup.

set -uo pipefail

CLONE="${PARTIS_ZIG_NIGHTLY_CLONE:-$HOME/re/pz-zig-nightly-ci}"
LOGDIR="${PARTIS_ZIG_NIGHTLY_LOGDIR:-$HOME/re/partis-zig-nightly-logs}"
TOPIC="${PARTIS_ZIG_NIGHTLY_TOPIC:-bip-matsen-partis-zig-nightly}"
RETAIN_DAYS="${PARTIS_ZIG_NIGHTLY_RETAIN:-14}"
REF="${PARTIS_ZIG_NIGHTLY_REF:-origin/main}"

if [[ ! -d "$CLONE/.git" ]]; then
    echo "FATAL: $CLONE is not a git clone. See deploy/nightly-ci/README.md for setup." >&2
    exit 2
fi

if ! command -v zig &>/dev/null; then
    echo "FATAL: zig not found on PATH." >&2
    exit 2
fi

cd "$CLONE" || { echo "FATAL: cd $CLONE failed" >&2; exit 2; }
git fetch origin --quiet || { echo "FATAL: git fetch failed" >&2; exit 2; }
git reset --hard "$REF" --quiet || { echo "FATAL: git reset --hard $REF failed" >&2; exit 2; }
git submodule update --init --quiet packages/zig-core || { echo "FATAL: submodule update failed" >&2; exit 2; }

mkdir -p "$LOGDIR" || { echo "FATAL: could not create $LOGDIR" >&2; exit 2; }
find "$LOGDIR" -maxdepth 1 -name '*.log' -type f -mtime +"$RETAIN_DAYS" -delete 2>/dev/null || true

SHA="$(git rev-parse --short HEAD)"
FULL_SHA="$(git rev-parse HEAD)"
DATE="$(date -u +%Y-%m-%d)"
LOG="$LOGDIR/$DATE-$SHA.log"

# Truncate (not append) at the start of this run: the filename already
# changes day-by-day (date + SHA), so truncating just keeps a same-SHA
# manual rerun from doubling up its own output.
run() {
    local start_ts end_ts status
    start_ts="$(date +%s)"

    echo "=== partis zig-core nightly Debug leak check ==="
    echo "commit:   $FULL_SHA"
    echo "short:    $SHA"
    echo "ref:      $REF"
    echo "clone:    $CLONE"
    echo "zig:      $(command -v zig) ($(zig version 2>/dev/null))"
    echo "host:     $(hostname)"
    echo "started:  $(date -Iseconds)"
    echo

    echo "--- zig build test (packages/zig-core, Debug -- default optimize mode) ---"
    (cd packages/zig-core && zig build test)
    status=$?

    if (( status == 0 )); then
        echo
        echo "--- zig build -Doptimize=Debug (packages/zig-core) ---"
        (cd packages/zig-core && zig build -Doptimize=Debug)
        status=$?
    fi

    if (( status == 0 )); then
        echo
        echo "--- partis-test.py --no-tree-gen --no-per-base-mutation (Debug binaries) ---"
        # shellcheck source=/dev/null
        source .venv/bin/activate
        partis-test.py --no-tree-gen --no-per-base-mutation \
            --bcrham-binary "$CLONE/packages/zig-core/zig-out/bin/partis-zig-core" \
            --ig-sw-binary "$CLONE/packages/zig-core/zig-out/bin/partis-zig-igsw"
        status=$?
    fi

    end_ts="$(date +%s)"
    echo
    echo "finished: $(date -Iseconds)"
    echo "elapsed:  $(( end_ts - start_ts ))s"
    echo "exit:     $status"
    return $status
}

run >"$LOG" 2>&1
STATUS=$?

# GPA prints "[gpa] (err): memory address 0x... leaked" (and an
# "error: '<test>' leaked: ..." wrapper in `zig build test` output) --
# grep case-insensitively for "leaked" rather than anchoring on either
# exact shape, since both contain the same word.
LEAK_LINES="$(grep -ic 'leaked' "$LOG" || true)"

if (( STATUS != 0 )) || (( LEAK_LINES > 0 )); then
    msg="partis zig-core nightly FAILED on ${HOSTNAME:-$(hostname)} @ $SHA"$'\n'"exit:       $STATUS"$'\n'"leak lines: $LEAK_LINES"$'\n'"log:        $LOG"$'\n'$'\n'"$(tail -40 "$LOG")"
    curl --silent --show-error --max-time 20 \
        -H "Title: partis zig-core nightly FAILED @ $SHA" \
        -H "Priority: high" \
        -H "Tags: warning,partis,zig,nightly" \
        -d "$msg" \
        "https://ntfy.sh/$TOPIC" \
        >>"$LOG" 2>&1 || echo "WARN: ntfy post failed" | tee -a "$LOG"
fi

exit $STATUS
