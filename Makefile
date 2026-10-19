export PATH := $(if $(wildcard .venv/bin),$(CURDIR)/.venv/bin:$(PATH),$(PATH))
.DEFAULT_GOAL := all

PIPELINE_TESTS ?= $(wildcard test/test_*_pipeline.sh)

.PHONY: test unit-tests pipeline-tests

test: unit-tests pipeline-tests

unit-tests:
	python -m pytest test/

pipeline-tests:
	@if [ -z "$(strip $(PIPELINE_TESTS))" ]; then echo "no pipeline tests"; fi
	@set -e; for script in $(PIPELINE_TESTS); do echo "running $$script"; bash $$script; done

.PHONY: all zig zig-fast zig-safe zig-debug zig-test test-safe test-fast test-quick-zig test-zig-safe test-zig-fast

all: zig

zig:
	./bin/zig-build.sh

zig-fast:
	./bin/zig-build.sh --release-fast

zig-safe:
	./bin/zig-build.sh --release-safe

zig-debug:
	./bin/zig-build.sh --debug

zig-test:  # zig unit tests + vendored-source sync checks (needs packages/{ham,ig-sw} checked out to compare anything)
	./bin/zig-build.sh --test

test-safe: zig-safe
	./bin/partis-test.py --quick --zig

test-fast: zig-fast
	./bin/partis-test.py --quick --zig

test-quick-zig: test-fast

test-zig-safe: test-safe

test-zig-fast: test-fast
