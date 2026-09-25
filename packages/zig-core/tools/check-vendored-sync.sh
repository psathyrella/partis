#!/bin/sh
# Verify that a C source vendored into packages/zig-core is byte-identical to
# the upstream copy in the sibling submodule that the C/C++ backends compile
# (issue #403). Drift means the two backends silently run different code.
#
# Usage: check-vendored-sync.sh <vendored-file> <upstream-file>
#
# Exits 0 with a note when the upstream file is missing: an uninitialized
# submodule is an ordinary state for a fresh clone (docs/install.md has people
# init them explicitly), and `zig build test` should not read as a broken build
# just because the check has nothing to compare against.
set -eu

vendored=$1
upstream=$2

if [ ! -f "$upstream" ]; then
    echo "note: skipping sync check for $vendored" >&2
    echo "      ($upstream not found -- run 'git submodule update --init' in the partis root to enable it)" >&2
    exit 0
fi

if ! diff -u "$vendored" "$upstream"; then
    echo "error: $vendored has drifted from $upstream" >&2
    echo "       these must stay byte-identical (issue #403); reconcile them and re-run." >&2
    exit 1
fi
