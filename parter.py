#!/usr/bin/env python
import sys
import os
import thread
import csv
import itertools
import operator
import pysam
import contextlib
from subprocess import check_output,check_call
import utils
from opener import opener
from hmmwriter import HmmWriter

# "IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F"
# v: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
# d: GTATTACTATGGTTCGGGGAGTTATTATAAC
# j: ACTACTTTGACTACTGGGGCCAGGGA

# ----------------------------------------------------------------------------------------
human = 'A'
naivety = 'M'
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')

seqfname = 'bcell/seq.fa'
bamfile = seqfname.replace('.fa', '.bam')

# ----------------------------------------------------------------------------------------
def run_gene_set(gene_set, k_v, k_d, v_right_length, best_j):
    tmp_hmmdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())
    if not os.path.exists(tmp_hmmdir):  # use a tmp dir specific to this process for the hmm file fifos
        os.makedirs(tmp_hmmdir)
    # write hmm fifos
    for gene in gene_set:
        writer = HmmWriter('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety,
                           'bcell', gene, naivety, germline_seqs[utils.get_region(gene)][gene], v_right_length=v_right_length)
        outfname = tmp_hmmdir + '/' + utils.sanitize_name(gene) + '.hmm'
        if os.path.exists(outfname):
            os.remove(outfname)
        # os.mkfifo(outfname)
        # thread.start_new_thread(writer.write, (outfname,))  # start a new thread so the fifo doesn't block
        writer.write(outfname)


    d_fuzz = 5
    if 'IGHJ4*' in best_j and germline_seqs['d'][gene_set[1]][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
        d_fuzz = 10
    check_call('./stochhmm -viterbi -hmmtype single -k_v_guess ' + str(k_v) + ' -k_d_guess ' + str(k_d) + ' -v_fuzz 3 -d_fuzz ' + str(d_fuzz) + ' -only_genes \'' + ':'.join(gene_set) + '\' -hmmdir ' + tmp_hmmdir + '', shell=True)

    # remove hmm fifos
    for fname in os.listdir(tmp_hmmdir):
        if fname.endswith(".hmm"):
            os.remove(tmp_hmmdir + "/" + fname)
    os.rmdir(tmp_hmmdir)

# ----------------------------------------------------------------------------------------
def run():
    check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 5 --min-score 5 ' + seqfname + ' ' + bamfile + ' 2>/dev/null', shell=True)
    match_names = {}
    query_bounds, germline_bounds = {}, {}
    with contextlib.closing(pysam.Samfile(bamfile)) as bam:
        grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
        for _, reads in grouped:
            reads = list(reads)
            primary = next((r for r in reads if not r.is_secondary), None)
            query_seq = primary.seq
            best_j = ''
            for read in reads:
                read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
                gene = bam.references[read.tid]
                region = utils.get_region(gene)
                if region == 'j' and best_j == '':  # keep track of the primary j, since it is excised before aligning the ds, which complicates things a bit when choosing d_fuzz
                    assert 'J1P' not in gene  # I think we can remove this version (never see it), but I'm putting a check in here just in case
                    best_j = gene
                if 'J1P' in gene or 'J3P' in gene:
                    continue
                # TODO I think connor's code excises the *best* v match, so this isn't right
                if region not in match_names:
                    match_names[region] = []
                match_names[region].append(gene)
                query_bounds[gene] = (read.qstart, read.qend)
                germline_bounds[gene] = (read.pos, read.aend)
                # score = read.tags[0][1]
                # print gene,score
                # print '  %4d%4d  %s' % (read.qstart, read.qend, germline_seqs[region][gene][read.pos:read.aend])
                # print '  %4d%4d  %s' % (read.pos, read.aend, query_seq[read.qstart:read.qend])

    icombo = 0
    for gene_set in itertools.product(match_names['v'], match_names['d'], match_names['j']):
        if icombo > 5:  # TODO why stop at five? maybe more?
            break
        k_v = query_bounds[gene_set[0]][1]  # end of v match  combo[0][2]
        k_d = query_bounds[gene_set[1]][1] - query_bounds[gene_set[0]][1]  # end of d minus end of v  combo[1][2] - combo[0][2]
        # print '  ',gene_set,k_v,k_d
        v_right_length = len(germline_seqs['v'][gene_set[0]]) - germline_bounds[gene_set[0]][0]  # germline v length minus (germline) start of v match
        run_gene_set(gene_set, k_v, k_d, v_right_length, best_j)
        icombo += 1

if len(sys.argv) > 1 and sys.argv[1] == 'single':
    run()
    sys.exit()

with opener('r')('/home/dralph/Dropbox/work/recombinator/output/' + human + '/' + naivety + '/simu.csv') as simfile:
    reader = csv.DictReader(simfile)
    last_reco_id = -1
    for line in reader:
        if line['reco_id'] == last_reco_id:
            continue
        last_reco_id = line['reco_id']
        print '\n\n'
        print 'true:'
        utils.print_reco_event(germline_seqs, line, 0, 0)
        print 'inferred:'
        with opener('w')(seqfname) as seqfile:
            seqfile.write('>seq_1 NUKES\n')
            seqfile.write(line['seq'])
        run()

    # output = check_output('./stochhmm -seq ' + seqfile + ' -model bcell/d.hmm -viterbi -label', shell=True)
    # for line in output.splitlines():
    #     if line.find('>>') == 0:  # header line
    #         score_str = line[line.find('Score:') + 6 : ]
    #         score = float(score_str)
    #     elif line.find('IGH') == 0:  # list of states
    #         previous_gene_name = ''
    #         for state in line.split():  # list of states we passed through
    #             if state == 'i':  # insert state
    #                 continue
    #             gene_name = state[:state.rfind('_')]  # strip off the position label from the state name
    #             if previous_gene_name == '':
    #                 previous_gene_name = gene_name
    #             else:
    #                 assert gene_name == previous_gene_name
