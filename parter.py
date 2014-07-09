#!/usr/bin/env python
import sys
from subprocess import check_output,check_call
import utils

# "IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F"
# v: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
# d: GTATTACTATGGTTCGGGGAGTTATTATAAC
# j: ACTACTTTGACTACTGGGGCCAGGGA
try:
    check_call('./stochhmm -viterbi -hmmtype single -debug -k_v_guess 100 -k_d_guess 32 -v_fuzz 1 -d_fuzz 1 -only_genes \'IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F\'', shell=True)
except:
    print 'hrg'
sys.exit()
output = check_output('./stochhmm -seq bcell/seq.fa -model bcell/d.hmm -viterbi -label', shell=True)
for line in output.splitlines():
    if line.find('>>') == 0:  # header line
        score_str = line[line.find('Score:') + 6 : ]
        score = float(score_str)
    elif line.find('IGH') == 0:  # list of states
        previous_gene_name = ''
        for state in line.split():  # list of states we passed through
            if state == 'i':  # insert state
                continue
            gene_name = state[:state.rfind('_')]  # strip off the position label from the state name
            if previous_gene_name == '':
                previous_gene_name = gene_name
            else:
                assert gene_name == previous_gene_name
