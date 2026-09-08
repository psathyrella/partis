#!/usr/bin/env bash
# Build the in-tree Zig backend (drop-in replacements for bcrham and ig-sw).
#
# Usage:
#   bin/zig-build.sh
#
# After building, pass --zig to partis:
#   partis --zig annotate ...
set -euo pipefail

ZIG_VERSION="0.15.2"
PARTIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ZIG_CORE_DIR="$PARTIS_DIR/packages/zig-core"
ZIG_EXE_FINAL="$ZIG_CORE_DIR/zig-out/bin/partis-zig-core"
ZIG_IGSW_EXE_FINAL="$ZIG_CORE_DIR/zig-out/bin/partis-zig-igsw"

echo "==> Building Zig backend in: $ZIG_CORE_DIR"

# ── 1. Ensure Zig is available ──────────────────────────────────────────────

find_zig() {
    # Prefer a zig in PATH that matches the required version
    if command -v zig &>/dev/null; then
        local v
        v=$(zig version 2>/dev/null || echo "")
        if [[ "$v" == "$ZIG_VERSION" ]]; then
            echo "zig"
            return
        fi
    fi
    # Check if we already downloaded it
    local downloaded="$ZIG_CORE_DIR/.zig-install/zig-$ZIG_VERSION/zig"
    if [[ -x "$downloaded" ]]; then
        echo "$downloaded"
        return
    fi
    echo ""
}

ZIG=$(find_zig)
if [[ -z "$ZIG" ]]; then
    echo "==> Downloading Zig $ZIG_VERSION..."
    OS=$(uname -s | tr '[:upper:]' '[:lower:]')
    ARCH=$(uname -m)
    case "$ARCH" in
        x86_64)  ARCH="x86_64" ;;
        aarch64|arm64) ARCH="aarch64" ;;
        *) echo "ERROR: Unsupported architecture: $ARCH" >&2; exit 1 ;;
    esac
    case "$OS" in
        linux)  OS_TAG="linux" ;;
        darwin) OS_TAG="macos" ;;
        *) echo "ERROR: Unsupported OS: $OS" >&2; exit 1 ;;
    esac
    ZIG_TARBALL="zig-${ARCH}-${OS_TAG}-${ZIG_VERSION}.tar.xz"
    ZIG_URL="https://ziglang.org/download/${ZIG_VERSION}/${ZIG_TARBALL}"
    ZIG_DIR="$ZIG_CORE_DIR/.zig-install/zig-${ZIG_VERSION}"
    mkdir -p "$ZIG_DIR"
    curl -fsSL "$ZIG_URL" | tar -xJ -C "$ZIG_DIR" --strip-components=1
    ZIG="$ZIG_DIR/zig"
    echo "    Zig $ZIG_VERSION installed at: $ZIG"
else
    echo "==> Using Zig: $ZIG ($($ZIG version))"
fi

# ── 2. Build ─────────────────────────────────────────────────────────────────

export PARTIS_DIR
echo "==> Building partis-zig-core (this takes ~30s)..."
(cd "$ZIG_CORE_DIR" && "$ZIG" build -Doptimize=ReleaseFast)

if [[ ! -x "$ZIG_EXE_FINAL" ]]; then
    echo "ERROR: build succeeded but executable not found at $ZIG_EXE_FINAL" >&2
    exit 1
fi

echo ""
echo "Build complete:"
echo "  $ZIG_EXE_FINAL"
echo "  $ZIG_IGSW_EXE_FINAL"
echo ""
echo "To use the Zig backend, pass --zig to partis:"
echo "  partis --zig annotate ..."
