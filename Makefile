.PHONY: all zig zig-fast zig-safe zig-debug test-safe test-quick-zig

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

test-quick-zig: zig-fast
	./bin/partis-test.py --quick --zig
