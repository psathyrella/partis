import time
# import_start = time.time()
import sys
import os
import csv
import itertools
import StringIO
import operator
import pysam
import contextlib
from subprocess import Popen, check_call, PIPE
import utils
from opener import opener
from hmmwriter import HmmWriter

# print 'import time: %.3f' % (time.time() - import_start)

# "IGHV1-18*01:IGHD3-10*01:IGHJ4*02_F"
# v: CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA
# d: GTATTACTATGGTTCGGGGAGTTATTATAAC
# j: ACTACTTTGACTACTGGGGCCAGGGA

# ----------------------------------------------------------------------------------------

class PartitionDriver(object):
    def __init__(self, datadir, args):
        self.datadir = datadir
        self.args = args
        self.germline_seqs = utils.read_germlines()
        self.workdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())  # use a tmp dir specific to this process for the hmm input file
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.reco_info = self.read_sim_file()  # read simulation info and write sw input file
        self.pair_scores = {}
        self.default_v_fuzz = 4
        self.default_d_fuzz = 5
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        swinfname = self.workdir + '/seq.fa'  # file for input to s-w step
        swoutfname = swinfname.replace('.fa', '.bam')
        print 's-w'
        self.run_smith_waterman(swinfname, swoutfname)
        sw_csv_fname = self.workdir + '/sw_out.csv'  # csv file with info from sw step
        self.read_smith_waterman(swoutfname, sw_csv_fname)  # read sw bam output, collate it, and write to csv for hmm input
        hmm_seqfname = self.workdir + '/seq.fa'
        for query_name in self.reco_info:  # loop over query seqs, running hmm for each
            if not self.args.pair:
                print 'true:'
                utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], 0, 0)
                if self.args.algorithm == 'viterbi':
                    self.write_hmm_seqfile(hmm_seqfname, query_name, '')
                    self.run_hmm(hmm_seqfname, sw_csv_fname, query_name, '')
            else:
                for second_query_name in self.reco_info:
                    if second_query_name == query_name:
                        continue
                    if int(query_name) + int(second_query_name) in self.pair_scores:  # already did this pair
                        continue
                    print '%20s %20s   %d' % (query_name, second_query_name, self.reco_info[query_name]['reco_id'] == self.reco_info[second_query_name]['reco_id'])
                    self.write_hmm_seqfile(hmm_seqfname, query_name, second_query_name)
                    self.run_hmm(hmm_seqfname, sw_csv_fname, query_name, second_query_name)

        os.remove(swoutfname)
        os.remove(sw_csv_fname)
        os.remove(hmm_seqfname)
        os.rmdir(self.workdir)
                
    
    # ----------------------------------------------------------------------------------------
    def read_sim_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        reco_info = {}  # generator truth information
        with opener('r')(self.args.simfile) as simfile:
            last_reco_id = -1  # only run on the first seq in each reco event. they're all pretty similar
            reader = csv.DictReader(simfile)
            for line in reader:
                # if line['reco_id'] == last_reco_id:
                #     continue
                reco_info[line['unique_id']] = line
                last_reco_id = line['reco_id']
        return reco_info

    # ----------------------------------------------------------------------------------------
    def run_stochhmm(self, cmd_str, query_name, second_query_name):
        matchlist = []
    
        if self.args.debug == 2:  # not sure why, but popen thing hangs with debug 2
            check_call(cmd_str, shell=True)
            return matchlist
    
        hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        hmm_proc.wait()
        hmm_out, hmm_err = hmm_proc.communicate()
        if self.args.debug:
            print 'OUT\n',hmm_out
        # run_time = time.time()
        # print 'hmm run time: %.3f' % (run_time - write_stop)
    
        try:  # most likely reason for failure is something got kicked into stderr besides the csv info
            if self.args.algorithm == 'viterbi':
                hmmreader = csv.DictReader(StringIO.StringIO(hmm_err))
                for hmmline in hmmreader:
                    matchlist.append(hmmline)
            elif self.args.algorithm == 'forward':
                total_score = float(hmm_err)
                print '    total score: ',total_score
                assert total_score < 0.0 and total_score > -9999999.9
                self.pair_scores[int(query_name) + int(second_query_name)] = total_score  # is that a hacky way to hash it so I can reverse the order and get the same entry? seems ok, I guess
            return matchlist
        except:
            print 'ERR\n',hmm_err
            sys.exit()
    
    # ----------------------------------------------------------------------------------------
    def run_hmm(self, seqfname, sw_csv_fname, query_name, second_query_name):
        # start = time.time()
        # read info from s-w
        all_gene_str = ''
        best_genes = {}
        k_v, k_d, v_fuzz, d_fuzz, v_right_length = 0, 0, 0, 0, 0
        with opener('r')(sw_csv_fname) as sw_outfile:
            reader = csv.DictReader(sw_outfile)
            for line in reader:
                if line['unique_id'] != query_name:  # TODO note that I'm only using the s-w info for the *first* sequence. I *need* to add in the info from the second
                    continue
                k_v = int(line['k_v'])
                k_d = int(line['k_d'])
                v_fuzz = int(line['v_fuzz'])
                d_fuzz = int(line['d_fuzz'])
                v_right_length = int(line['v_right_length'])
                assert k_v > 0 and k_d > 0 and v_fuzz > 0 and d_fuzz > 0 and v_right_length > 0
                for region in utils.regions:
                    best_genes[region] = line['best_' + region]
                all_gene_str = line['all']
                
        # write hmm files
        for gene in all_gene_str.split(':'):  # TODO note that this takes more time than actually running the hmm
            writer = HmmWriter(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety,
                               'bcell', gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=v_right_length)
            outfname = self.workdir + '/' + utils.sanitize_name(gene) + '.hmm'
            if os.path.exists(outfname):
                os.remove(outfname)
            # os.mkfifo(outfname)  # what the hell? This is actually almost *twice* as *slow* with fifos instead of files
            # thread.start_new_thread(writer.write, (outfname,))  # start a new thread so the fifo doesn't block
            writer.write(outfname)
    
        # write_stop = time.time()
        # print 'hmm write time: %.3f' % (write_stop - start)
        if self.args.algorithm == 'viterbi':
            print '\ninferred:'
    
        assert 'J1P' not in best_genes['j'] and 'J2P' not in best_genes['j']  # I think we can remove this version (we never see it), but I'm putting a check in here just in case
        if 'IGHJ4*' in best_genes['j'] and self.germline_seqs['d'][best_genes['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
            d_fuzz = 10
        # TODO since the s-w stuff excises the *best* v and *best* j, these k_v and k_d can be quite far off. Nonetheless seems ok a.t.m.
        # NOTE this call costs about 0.01-0.02 seconds if StochHMM.cpp does more or less nothing
        cmd_str = './stochhmm'
        cmd_str += ' -' + self.args.algorithm
        if self.args.pair:
            cmd_str += ' -hmmtype pair'
        else:
            cmd_str += ' -hmmtype single'
        cmd_str += ' -debug ' + str(self.args.debug)
        cmd_str += ' -seq ' + seqfname
        cmd_str += ' -k_v_guess ' + str(k_v) + ' -k_d_guess ' + str(k_d)
        cmd_str += ' -v_fuzz ' + str(v_fuzz) + ' -d_fuzz ' + str(d_fuzz)
        cmd_str += ' -only_genes \'' + all_gene_str + '\''  # hm, wait, do I need these escaped quotes?
        cmd_str += ' -hmmdir ' + self.workdir

        matchlist = self.run_stochhmm(cmd_str, query_name, second_query_name)

        for matchline in matchlist:
            utils.print_reco_event(self.germline_seqs, matchline, 0, 0)
            if self.args.pair:
                matchline['seq'] = 'CACCATTTCCAGAGACAACTCCATGAGCTCCCTGTCTCTTCAAATGAACAGTCTGAGAGCCGAGGTCACGTCTGTGTATTACTGTGCGTTAGACAGTGGCCGGTTCTGCAGAGGCTCCTGGTTCCAGGGA'
                utils.print_reco_event(self.germline_seqs, matchline, 0, 0, True)
                
    
        # remove hmm model files
        for fname in os.listdir(self.workdir):
            if fname.endswith(".hmm"):
                os.remove(self.workdir + "/" + fname)
    
    # ----------------------------------------------------------------------------------------
    def write_hmm_seqfile(self, seqfname, query_name, second_query_name):
        # write seq to file for stochhmm to read
        with opener('w')(seqfname) as seqfile:
            seqfile.write('>' + query_name + ' NUKES\n')
            seqfile.write(self.reco_info[query_name]['seq'] + '\n')
            if second_query_name != '':
                seqfile.write('>' + second_query_name + ' NUKES\n')
                seqfile.write(self.reco_info[second_query_name]['seq'] + '\n')
    
    # ----------------------------------------------------------------------------------------
    def run_smith_waterman(self, infname, outfname):
        """
        Run smith-waterman alignment on the seqs in <infname>, and toss all the top matches into <outfname>.
        Then run through <outfname> to get the top hits and their locations to pass to the hmm.
        Then run the hmm on each gene set.
        """
        with opener('w')(infname) as swinfile:  # write *all* the input seqs to file, i.e. run s-w on all of 'em at once
            for query_name,line in self.reco_info.iteritems():
                swinfile.write('>' + query_name + ' NUKES\n')
                swinfile.write(line['seq'] + '\n')
        # start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 3 --mismatch 1 --gap-open 100 ' + infname + ' ' + outfname + ' 2>/dev/null', shell=True)
        # print 's-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_smith_waterman(self, infname, outfname):
        gene_choice_probs = utils.read_overall_gene_prob(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety)
        with contextlib.closing(pysam.Samfile(infname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))

            with opener('w')(outfname) as outfile:
                columns = ('unique_id', 'k_v', 'k_d', 'v_fuzz', 'd_fuzz', 'v_right_length', 'best_v', 'best_d', 'best_j', 'all')
                writer = csv.DictWriter(outfile, columns)
                writer.writeheader()

                for _, reads in grouped:  # loop over query sequences
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
        
                        raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs. TODO oh wait shit this isn't right. raw_score isn't a prob. what the hell is it, anyway?
                        choice_prob = 0.0  # set to zero if we didn't see it in data. kinda hacky, I suppose
                        if gene in gene_choice_probs[region]:
                            choice_prob = gene_choice_probs[region][gene]
                        score = choice_prob * raw_score  # multiply by the probability to choose this gene
                        all_match_names[region].append((score,gene))
                        all_query_bounds[gene] = (read.qstart, read.qend)
                        all_germline_bounds[gene] = (read.pos, read.aend)
        
                    # append the best matches to match_names and work out how much v_fuzz we need
                    best = {}
                    match_names = {}
                    n_matches = {'v':0, 'd':0, 'j':0}
                    n_used = {'v':0, 'd':0, 'j':0}
                    v_fuzz = self.default_v_fuzz  # TODO I have the defaults adjusted in two places at the moment
                    d_fuzz = self.default_d_fuzz
                    for region in utils.regions:
                        all_match_names[region] = sorted(all_match_names[region], reverse=True)
                        match_names[region] = []
                    for region in utils.regions:
                        for score,gene in all_match_names[region]:
                            n_matches[region] += 1
                            if n_matches[region] > self.args.n_max_per_region:  # only take the top few from each region
                                # TODO should use *lots* of d matches, but fewer vs and js
                                # assert False  # TODO also should loop over *way* more k_d than k_v
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
        
                            if self.args.debug:
                                buff_str = (17 - len(gene)) * ' '
                                print '%8s%s%s%6.1f' % (' ', utils.color_gene(gene), buff_str, score),
                                glbounds = all_germline_bounds[gene]
                                qrbounds = all_query_bounds[gene]
                                print ' %4d%4d   %s' % (glbounds[0], glbounds[1], self.germline_seqs[region][gene][glbounds[0]:glbounds[1]])
                                print '%31s  %4d%4d   %s' % ('', qrbounds[0], qrbounds[1], utils.color_mutants(self.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]))
                                    
                    # print how many of the available matches we used
                    if self.args.debug:
                        print '  used',
                        for region in utils.regions:
                            if region != 'v':
                                print '      ',
                            print ' %d / %d in %s' % (n_used[region], n_matches[region], region)
        
                    # write sw info to file for hmm
                    k_v = all_query_bounds[best['v']][1]  # end of v match
                    k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v
                    if k_d < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment. 
                        k_d = max(5, k_d)
                        d_fuzz = max(d_fuzz, 10)
                        
                    v_right_length = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match


                    row = {}
                    row['unique_id'] = query_name
                    row['k_v'] = k_v
                    row['k_d'] = k_d
                    row['v_fuzz'] = v_fuzz
                    row['d_fuzz'] = d_fuzz
                    row['v_right_length'] = v_right_length
                    for region in utils.regions:
                        row['best_' + region] = best[region]
                    row['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])
                    writer.writerow(row)
