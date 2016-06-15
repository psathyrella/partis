#!/bin/sh

set -e
set -u

lein uberjar
JAR=target/ighutil-0.1.0-SNAPSHOT-standalone.jar
java -jar $JAR reset-primary -i testdata/test_all.bam -o test_reset.bam
java -jar $JAR calculate-match-probability -r testdata/ighvdj.fasta -i test_reset.bam -o test_with_prob.bam
