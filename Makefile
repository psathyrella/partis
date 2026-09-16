export PATH := $(if $(wildcard .venv/bin),$(CURDIR)/.venv/bin:$(PATH),$(PATH))

.PHONY: all zig zig-fast zig-safe zig-debug test-safe test-fast test-quick-zig test-zig-safe test-zig-fast

all: zig

zig:
	./bin/zig-build.sh

zig-fast:
	./bin/zig-build.sh --release-fast

zig-safe:
	./bin/zig-build.sh --release-safe

zig-debug:
	./bin/zig-build.sh --debug

test-safe: zig-safe
	./bin/partis-test.py --quick --zig

test-fast: zig-fast
	./bin/partis-test.py --quick --zig

test-quick-zig: test-fast

test-zig-safe: test-safe

test-zig-fast: test-fast
