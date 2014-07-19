#!/usr/bin/env python
import time
import_start = time.time()
import sys
import os
import csv
import itertools
import argparse
import StringIO
import operator
import pysam
import contextlib
from subprocess import Popen, check_call, PIPE
import utils
from opener import opener
from hmmwriter import HmmWriter

print 'import time: %.3f' % (time.time() - import_start)

# "IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F"
# v: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
# d: GTATTACTATGGTTCGGGGAGTTATTATAAC
# j: ACTACTTTGACTACTGGGGCCAGGGA

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simfile')  #, default='/home/dralph/Dropbox/work/recombinator/output/' + human + '/' + naivety + '/simu.csv')
parser.add_argument('--debug', type=int, default=0)
# parser.add_argument('--seq')
args = parser.parse_args()

human = 'A'
naivety = 'M'
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')
datadir = '/home/dralph/Dropbox/work/recombinator/data'

# ----------------------------------------------------------------------------------------
def run_hmm(hmm_seqfname, tmp_hmmdir, match_names, best_gene, k_v, k_d, v_right_length):
    start = time.time()

    # write hmm files
    all_genes = match_names['v'] + match_names['d'] + match_names['j']
    for gene in all_genes:  # TODO note that this takes more time than actually running the hmm
        # TODO note that if you're using all_genes here, you need to add an assertion that the different v (or dj) versions don't need drastically different k_v (or k_d). On the other hand, since the s-w only excises the *best*, I guess I'm screwed anyway if they're drastically different
        writer = HmmWriter(datadir + '/human-beings/' + human + '/' + naivety,
                           'bcell', gene, naivety, germline_seqs[utils.get_region(gene)][gene], v_right_length=v_right_length)
        outfname = tmp_hmmdir + '/' + utils.sanitize_name(gene) + '.hmm'
        if os.path.exists(outfname):
            os.remove(outfname)
        # os.mkfifo(outfname)  # what the hell? This is actually almost *twice* as *slow* with fifos instead of files
        # thread.start_new_thread(writer.write, (outfname,))  # start a new thread so the fifo doesn't block
        writer.write(outfname)

    write_stop = time.time()
    print 'hmm write time: %.3f' % (write_stop - start)

    assert 'J1P' not in best_gene['j'] and 'J2P' not in best_gene['j']  # I think we can remove this version (we never see it), but I'm putting a check in here just in case
    d_fuzz = 5
    if 'IGHJ4*' in best_gene['j'] and germline_seqs['d'][best_gene['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
        d_fuzz = 10
    # TODO since the s-w stuff excises the *best* v and *best* j, these k_v and k_d can be quite far off. Nonetheless seems ok a.t.m.
    # NOTE this call costs about 0.01-0.02 seconds if StochHMM.cpp does more or less nothing
    hmm_proc = Popen('./stochhmm -viterbi -hmmtype single -debug ' + str(args.debug) + ' -seq ' + hmm_seqfname + ' -k_v_guess ' + str(k_v) + ' -k_d_guess ' + str(k_d) + ' -v_fuzz 4 -d_fuzz ' + str(d_fuzz) \
                     + ' -only_genes \'' + ':'.join(all_genes) + '\' -hmmdir ' + tmp_hmmdir, shell=True, stdout=PIPE, stderr=PIPE)
    hmm_proc.wait()

    run_time = time.time()
    print 'hmm run time: %.3f' % (run_time - write_stop)

    print '\ninferred:'
    hmm_out, hmm_err = hmm_proc.communicate()
    if True:  #args.debug:
        print 'OUT\n',hmm_out

    hmmreader = csv.DictReader(StringIO.StringIO(hmm_err))
    for hmmline in hmmreader:
        utils.print_reco_event(germline_seqs, hmmline, 0, 0)

    # remove hmm fifos
    for fname in os.listdir(tmp_hmmdir):
        if fname.endswith(".hmm"):
            os.remove(tmp_hmmdir + "/" + fname)

# ----------------------------------------------------------------------------------------
def run(reco_info_local, seqfname):
    """
    Run smith-waterman alignment on the seqs in <seqfname>, and toss all the top matches into <bamfname>.
    Then run through <bamfname> to get the top hits and their locations to pass to the hmm.
    Then run the hmm on each gene set.
    """
    start = time.time()
    bamfname = seqfname.replace('.fa', '.bam')
    # large gap-opening penalty: we want *no* gaps in the middle of the alignments
    # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
    check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 3 --mismatch 1 --gap-open 100 ' + seqfname + ' ' + bamfname + ' 2>/dev/null', shell=True)
    print 's-w time: %.3f' % (time.time()-start)
    gene_choice_probs = utils.read_overall_gene_prob(datadir + '/human-beings/' + human + '/' + naivety)
    with contextlib.closing(pysam.Samfile(bamfname)) as bam:
        grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
        for _, reads in grouped:  # loop over each query sequence
            reads = list(reads)
            primary = next((r for r in reads if not r.is_secondary), None)
            query_seq = primary.seq
            query_name = primary.qname
            raw_best = {}
            all_match_names = {}
            for region in utils.regions:
                all_match_names[region] = []
            all_query_bounds, all_germline_bounds = {}, {}
            for read in reads:  # loop over the matches found for each query sequence
                read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
                gene = bam.references[read.tid]
                region = utils.get_region(gene)

                if region not in raw_best:  # best v, d, and j before multiplying by gene choice probs. needed 'cause *these* are the v and j that get excised
                    raw_best[region] = gene

                if 'J1P' in gene or 'J3P' in gene:
                    continue

                raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs
                choice_prob = 0.0  # set to zero if we didn't see it in data. kinda hacky, I suppose
                if gene in gene_choice_probs[region]:
                    choice_prob = gene_choice_probs[region][gene]
                score = choice_prob * raw_score  # multiply by the probability to choose this gene
                all_match_names[region].append((score,gene))
                all_query_bounds[gene] = (read.qstart, read.qend)
                all_germline_bounds[gene] = (read.pos, read.aend)
                    
                # if True:  #args.debug:
                #     buff_str = (17 - len(gene)) * ' '
                #     print '%8s%s%s%6.1f' % (' ', utils.color_gene(gene), buff_str, score),
                #     print ' %4d%4d   %s' % (read.pos, read.aend, germline_seqs[region][gene][read.pos:read.aend])
                #     print '%31s  %4d%4d   %s' % ('', read.qstart, read.qend, utils.color_mutants(germline_seqs[region][gene][read.pos:read.aend], query_seq[read.qstart:read.qend]))

            best = {}
            match_names = {}
            n_matches = {'v':0, 'd':0, 'j':0}
            n_used = {'v':0, 'd':0, 'j':0}
            v_fuzz = 4
            for region in utils.regions:
                all_match_names[region] = sorted(all_match_names[region], reverse=True)
                match_names[region] = []
            for region in utils.regions:
                for score,gene in all_match_names[region]:
                    n_matches[region] += 1
                    if n_matches[region] > 5:  # only take the top few from each region
                        # TODO should use *lots* of d matches, but fewer vs and js
                        # assert False  # TODO also should loop over *way* more k_d than k_v
                        # TODO also need to assert that subleading matches have k_v k_d consistent with leading k_v k_d and fuzzes
                        continue
                    n_used[region] += 1

                    match_names[region].append(gene)
                    # check consistency with best match (since the best match is excised in s-w code, and because stochhmm is run with *one* k_v k_d set)
                    if region not in best:
                        best[region] = gene
                    else:
                        if region == 'v':
                            necessary_fuzz = abs(all_query_bounds[gene][1] - all_query_bounds[raw_best[region]][1])
                            if necessary_fuzz > v_fuzz-1:  # if this v match ends at a wildly different position than did the best v match, expand v_fuzz accordingly
                                v_fuzz = necessary_fuzz + 1
                    if True:  #args.debug:
                        buff_str = (17 - len(gene)) * ' '
                        print '%8s%s%s%6.1f' % (' ', utils.color_gene(gene), buff_str, score),
                        glbounds = all_germline_bounds[gene]
                        qrbounds = all_query_bounds[gene]
                        print ' %4d%4d   %s' % (glbounds[0], glbounds[1], germline_seqs[region][gene][glbounds[0]:glbounds[1]])
                        print '%31s  %4d%4d   %s' % ('', qrbounds[0], qrbounds[1], utils.color_mutants(germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]))
                            
            print '  used',
            for region in utils.regions:
                if region != 'v':
                    print '      ',
                print ' %d / %d in %s' % (n_used[region], n_matches[region], region)

            print 'true:'
            utils.print_reco_event(germline_seqs, reco_info_local[query_name], 0, 0)

            # write seq to file for stochhmm to read
            tmp_hmmdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())
            if not os.path.exists(tmp_hmmdir):  # use a tmp dir specific to this process for the hmm file fifos
                os.makedirs(tmp_hmmdir)
            hmm_seqfname = tmp_hmmdir + '/seq.fa'
            with opener('w')(hmm_seqfname) as hmm_seqfile:
                hmm_seqfile.write('>' + reco_info_local[query_name]['unique_id'] + ' NUKES\n')
                hmm_seqfile.write(reco_info_local[query_name]['seq'] + '\n')

            k_v = all_query_bounds[best['v']][1]  # end of v match
            k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v
            v_right_length = len(germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match
            # TODO add checks that none of the matches have k_v k_d v_right_length very different from the first one
            # for gene_set in itertools.product(match_names['v'], match_names['d'], match_names['j']):  #get_gene_set(match_names['v'], match_names['d'], match_names['j'])):
            start = time.time()
            run_hmm(hmm_seqfname, tmp_hmmdir, match_names, best, k_v, k_d, v_right_length)
            print 'hmm time: %.3f' % (time.time()-start)

            os.remove(hmm_seqfname)
            os.rmdir(tmp_hmmdir)

if args.simfile != None:
    print 'read infile'
    with opener('r')(args.simfile) as simfile:
        reader = csv.DictReader(simfile)
        reco_info = {}
        last_reco_id = -1  # only run on the first seq in each reco event. they're all pretty similar
        tmpseqfname = 'bcell/seq.fa'  # file for input to s-w step
        with opener('w')(tmpseqfname) as seqfile:
            for line in reader:
                if line['reco_id'] == last_reco_id:
                    continue
                reco_info[line['unique_id']] = line
                last_reco_id = line['reco_id']
                seqfile.write('>' + line['unique_id'] + ' NUKES\n')
                seqfile.write(line['seq'] + '\n')

        print '  done'
        run(reco_info, tmpseqfname)
