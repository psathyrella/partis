#!/usr/bin/env python
import sys
import os
import thread
from subprocess import check_output,check_call
import utils
from hmmwriter import HmmWriter

# "IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F"
# v: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
# d: GTATTACTATGGTTCGGGGAGTTATTATAAC
# j: ACTACTTTGACTACTGGGGCCAGGGA

human = 'A'
naivety = 'M'
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')
gene_set = ['IGHV1-18*01','IGHD3-10*01','IGHJ4*02_F']
tmp_hmmdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())
if not os.path.exists(tmp_hmmdir):  # use a tmp dir specific to this process for the hmm file fifos
    os.makedirs(tmp_hmmdir)
for gene in gene_set:
    writer = HmmWriter('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety, 'bcell', gene, naivety, germline_seqs[utils.get_region(gene)][gene])
    outfname = tmp_hmmdir + '/' + utils.sanitize_name(gene) + '.hmm'
    if os.path.exists(outfname):
        os.remove(outfname)
    os.mkfifo(outfname)
    thread.start_new_thread(writer.write, (outfname,))  # start a new thread so the fifo doesn't block

try:
    check_call('./stochhmm -viterbi -hmmtype single -k_v_guess 100 -k_d_guess 32 -v_fuzz 3 -d_fuzz 3 -only_genes \'' + ':'.join(gene_set) + '\' -hmmdir ' + tmp_hmmdir + '', shell=True)
except:
    print 'hrg'
for fname in os.listdir(tmp_hmmdir):
    if fname.endswith(".hmm"):
        os.remove(tmp_hmmdir + "/" + fname)
os.rmdir(tmp_hmmdir)
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
