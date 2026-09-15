# partis zig-core nightly / weekly Debug leak checks

Setup for the systemd user timers described in issue #405. Both jobs build `packages/zig-core` with `-Doptimize=Debug`, the only build mode where `main.zig` routes allocations through `std.heap.GeneralPurposeAllocator` (GPA) instead of `std.heap.c_allocator` — see the issue for why `ReleaseSafe`/`ReleaseFast` don't catch leaks in this codebase. Installing and enabling these units is a manual, one-time step; it is *not* done automatically by merging the PR that adds these files.

## Nightly: `partis-zig-nightly-test.{service,timer}`

Runs `zig build test` plus the full non-quick, non-paired `partis-test.py` tier against Debug binaries. Takes well under a minute on pax (measured 21s for the `partis-test.py` tier alone, issue #405). Posts to ntfy only on failure or a detected leak.

### One-time setup

```bash
# 1. Dedicated clone (kept separate from any interactive checkout).
#    packages/zig-core is plain tracked content, not a git submodule
#    (brought in-tree by fa6f199dc) -- a normal clone is all that's needed.
git clone git@github.com:psathyrella/partis.git ~/re/pz-zig-nightly-ci
cd ~/re/pz-zig-nightly-ci

# 2. Python venv the nightly script activates before calling partis-test.py.
#    Backend choice doesn't matter here -- the nightly script always passes
#    --bcrham-binary/--ig-sw-binary explicitly, overriding whichever backend
#    got built at install time -- but PARTIS_BACKEND=zig avoids needing
#    gcc/scons/gsl/yaml-cpp on the box at all (see docs/install.md).
python -m venv .venv
source .venv/bin/activate
PARTIS_BACKEND=zig pip install -e .
deactivate

# 3. Install the systemd units
mkdir -p ~/.config/systemd/user
cp deploy/nightly-ci/partis-zig-nightly-test.service ~/.config/systemd/user/
cp deploy/nightly-ci/partis-zig-nightly-test.timer   ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now partis-zig-nightly-test.timer

# 4. Verify
systemctl --user list-timers partis-zig-nightly-test.timer
systemctl --user start partis-zig-nightly-test.service   # run once by hand
journalctl --user -u partis-zig-nightly-test.service -n 50
```

### Configuration

All overridable via environment variables read by `scripts/nightly_zig_debug_test.sh` (set them in the `.service` file's `Environment=` lines, or a drop-in):

| Variable | Default | Meaning |
|---|---|---|
| `PARTIS_ZIG_NIGHTLY_CLONE` | `~/re/pz-zig-nightly-ci` | the dedicated clone from step 1 |
| `PARTIS_ZIG_NIGHTLY_LOGDIR` | `~/re/partis-zig-nightly-logs` | per-night log files, `<date>-<sha>.log`, retained 14 days |
| `PARTIS_ZIG_NIGHTLY_TOPIC` | `bip-matsen-partis-zig-nightly` | ntfy.sh topic to post failures/leaks to |
| `PARTIS_ZIG_NIGHTLY_RETAIN` | `14` | log retention, days |
| `PARTIS_ZIG_NIGHTLY_REF` | `origin/main` | git ref to reset the clone to before each run |

Subscribe to the ntfy topic (e.g. the ntfy Android/iOS app, or `curl -s ntfy.sh/bip-matsen-partis-zig-nightly/json`) to receive failure/leak notifications.

## Weekly: `partis-zig-weekly-leak-test.{service,timer}`

Not yet implemented — see issue #405 Phase 2. Will run rungs 3-4 of issue #342's validation ladder (500-seq / 5000-seq paired-loci partition) built `-Doptimize=Debug`, at a cadence that can absorb a multi-hour run.
