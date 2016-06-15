# ighutil -- Create executable script

CLJ := $(shell find src -name "*.clj" -or -name "*.java")
IGHUTIL ?= bin/ighutil

all: $(IGHUTIL)

$(IGHUTIL): $(CLJ) project.clj
	rm -f target/*.jar
	dev/lein uberjar
	cat bin/ighutil.template target/*-standalone.jar > bin/ighutil
	chmod +x bin/ighutil

test_reset.bam: testdata/test_all.bam $(IGHUTIL)
	$(IGHUTIL) reset-primary --no-filter-alleles -i testdata/test_all.bam -o test_reset.bam

test_with_prob.bam: test_reset.bam
	 $(IGHUTIL) calculate-match-probability -r testdata/ighvdj.fasta -i test_reset.bam --no-compress -o /dev/stdout | \
		samtools view -bF260 -o test_with_prob.bam -

test_with_prob.csv.bz2: test_with_prob.bam
	$(IGHUTIL) annotate-vdj -i test_with_prob.bam -o test_with_prob.csv.bz2

run-test: test_reset.bam test_with_prob.bam test_with_prob.csv.bz2

.PHONY: all run-test
