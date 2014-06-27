#!/usr/bin/env python

import sys
import csv
from subprocess import check_output
from opener import opener

ref_seq = 'GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGACAGGCTGCAGGGAAGGGGCTGAAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATCCTACGGAGACTCCGTGAAAGGCCGGTTCACCATCTCAAGAGAGAATTCCGAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAACCGAGGACACGGCCGTATATTACTGTGCGAAAGTACTACGGTGACTGCCCCACTACTTTGACTGCTGGGGCCCGGGA'
with opener('r')('/home/dralph/Dropbox/work/recombinator/output/C/M/simu.csv') as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        with opener('w')('bcell/seq.fa') as seqfile:
            seqfile.write('>seq_1 NUKES\n')
            seqfile.write(ref_seq + '\n')
            seqfile.write('>seq_2 NUKES\n')
            seqfile.write(line['seq'] + '\n')
            seqfile.close()
            output = check_output('./stochhmm -viterbi -hmmtype pair', shell=True)
            print output
        sys.exit()
