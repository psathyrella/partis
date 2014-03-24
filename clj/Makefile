# ighutil -- Create executable script


all:
	rm -f target/*.jar
	lein uberjar
	cat bin/ighutil.template target/*-standalone.jar > bin/ighutil
	chmod +x bin/ighutil

run-test: all
	bin/ighutil reset-primary -i testdata/test_all.bam -o test_reset.bam
	bin/ighutil calculate-match-probability -r testdata/ighvdj.fasta -i test_reset.bam -o test_with_prob.bam

.PHONY: all run-test
