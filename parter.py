#!/usr/bin/env python

from subprocess import check_output
import utils

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
