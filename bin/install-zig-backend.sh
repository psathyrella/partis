#!/usr/bin/env bash
# Install the partis-zig-core Zig backend (drop-in replacements for bcrham and ig-sw).
#
# Usage:
#   bin/install-zig-backend.sh [--install-dir DIR]
#
# After running, pass these flags to partis:
#   --bcrham-binary <install-dir>/partis-zig-core/zig-out/bin/partis-zig-core
#   --ig-sw-binary  <install-dir>/partis-zig-core/zig-out/bin/partis-zig-igsw
#
# The Zig backend requires no external libraries (GSL, OpenBLAS, etc.) and
# produces bit-for-bit identical annotations to the C++ bcrham and C ig-sw binaries.
set -euo pipefail

ZIG_VERSION="0.15.2"
ZIG_REPO="https://github.com/matsengrp/partis-zig-core.git"
DEFAULT_INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/zig-backend"

INSTALL_DIR="$DEFAULT_INSTALL_DIR"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-dir) INSTALL_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

ZIG_CORE_DIR="$INSTALL_DIR/partis-zig-core"
ZIG_EXE_FINAL="$ZIG_CORE_DIR/zig-out/bin/partis-zig-core"

echo "==> Installing partis Zig backend to: $INSTALL_DIR"
mkdir -p "$INSTALL_DIR"

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
    local downloaded="$INSTALL_DIR/zig-$ZIG_VERSION/zig"
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
    # Normalize arch names
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
    ZIG_TARBALL="zig-${OS_TAG}-${ARCH}-${ZIG_VERSION}.tar.xz"
    ZIG_URL="https://ziglang.org/download/${ZIG_VERSION}/${ZIG_TARBALL}"
    ZIG_DIR="$INSTALL_DIR/zig-${ZIG_VERSION}"
    mkdir -p "$ZIG_DIR"
    curl -fsSL "$ZIG_URL" | tar -xJ -C "$ZIG_DIR" --strip-components=1
    ZIG="$ZIG_DIR/zig"
    echo "    Zig $ZIG_VERSION installed at: $ZIG"
else
    echo "==> Using Zig: $ZIG ($(zig version 2>/dev/null || $ZIG version))"
fi

# ── 2. Clone or update partis-zig-core ──────────────────────────────────────

if [[ -d "$ZIG_CORE_DIR/.git" ]]; then
    echo "==> Updating partis-zig-core..."
    git -C "$ZIG_CORE_DIR" pull --ff-only
else
    echo "==> Cloning partis-zig-core..."
    git clone "$ZIG_REPO" "$ZIG_CORE_DIR"
fi

# ── 3. Build ─────────────────────────────────────────────────────────────────

# build.zig references ig-sw sources inside the partis repo.
# Set PARTIS_DIR so build.zig can find them.
PARTIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PARTIS_DIR

echo "==> Building partis-zig-core (this takes ~30s)..."
(cd "$ZIG_CORE_DIR" && "$ZIG" build -Doptimize=ReleaseFast)

if [[ ! -x "$ZIG_EXE_FINAL" ]]; then
    echo "ERROR: build succeeded but executable not found at $ZIG_EXE_FINAL" >&2
    exit 1
fi

ZIG_IGSW_EXE_FINAL="$ZIG_CORE_DIR/zig-out/bin/partis-zig-igsw"

echo ""
echo "Build complete:"
echo "  $ZIG_EXE_FINAL"
echo "  $ZIG_IGSW_EXE_FINAL"
echo ""
echo "To use the Zig backend, pass these flags to partis:"
echo "  --bcrham-binary $ZIG_EXE_FINAL"
echo "  --ig-sw-binary $ZIG_IGSW_EXE_FINAL"
echo ""
echo "Omit these flags to use the default C++ bcrham and C ig-sw binaries."
