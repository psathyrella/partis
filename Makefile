PIPELINE_TESTS ?= $(wildcard test/test_*_pipeline.sh)

.PHONY: test unit-tests pipeline-tests

test: unit-tests pipeline-tests

unit-tests:
	python -m pytest test/

pipeline-tests:
	@if [ -z "$(strip $(PIPELINE_TESTS))" ]; then echo "no pipeline tests"; fi
	@set -e; for script in $(PIPELINE_TESTS); do echo "running $$script"; bash $$script; done
