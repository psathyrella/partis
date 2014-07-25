import time
import_start = time.time()
import sys
import os
import csv
import re
import itertools
import StringIO
import operator
import pysam
import contextlib
from subprocess import Popen, check_call, PIPE
import utils
from opener import opener
from hmmwriter import HmmWriter

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, datadir, args, default_v_right_length=90):
        self.datadir = datadir
        self.args = args
        self.germline_seqs = utils.read_germlines()
        self.workdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())  # use a tmp dir specific to this process for the hmm input file
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.reco_info = self.read_sim_file()  # read simulation info and write sw input file
        self.pair_scores = {}
        self.default_v_fuzz = 4
        self.default_d_fuzz = 10  # note that the k_d guess can be almost worthless, since the s-w step tends to expand the j match over a lot of the d
        self.pairscorefname = 'pairwise-scores.csv'
        self.default_hmm_dir = 'bcell/hmms/' + self.args.human + '/' + self.args.naivety
        self.default_v_right_length = default_v_right_length
        self.sw_info = {}

        # self.fracs_in_common = {}
        # self.fracs_in_common['v'] = 0.0
        # self.fracs_in_common['d'] = 0.0
        # self.fracs_in_common['j'] = 0.0
        # self.delta_k_v = 0
        # self.delta_k_d = 0
        # self.n_pairs = 0
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        with opener('w')(self.pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id_1,unique_id_2,score\n')
        swinfname = self.workdir + '/seq.fa'  # file for input to s-w step
        swoutfname = swinfname.replace('.fa', '.bam')
        if self.args.run_sw:
            print 's-w...',
            sys.stdout.flush()
            self.run_smith_waterman(swinfname, swoutfname)
            print 'done'
        self.read_smith_waterman(swoutfname)  # read sw bam output, collate it, and write to csv for hmm input
        hmm_seqfname = self.workdir + '/seq.fa'
        for query_name in self.reco_info:  # loop over query seqs, running hmm for each
            if not self.args.pair:
                if self.args.algorithm == 'viterbi':
                    print 'true:'
                    utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], 0, 0)
                    self.write_hmm_seqfile(hmm_seqfname, query_name, '')
                    self.run_hmm(hmm_seqfname, query_name, '')
            else:
                for second_query_name in self.reco_info:  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
                    if second_query_name == query_name:
                        continue
                    if self.get_score_index(query_name, second_query_name) in self.pair_scores:  # already did this pair
                        continue
                    from_same_event = self.reco_info[query_name]['reco_id'] == self.reco_info[second_query_name]['reco_id']
                    # if self.have_different_gene_matches(query_name, second_query_name):
                    #     continue
                    print '%20s %20s   %d' % (query_name, second_query_name, from_same_event)
                    if self.args.algorithm == 'viterbi':
                        print 'true:'
                        utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], 0, 0)
                        print '  and maybe'
                        utils.print_reco_event(self.germline_seqs, self.reco_info[second_query_name], 0, 0, True)
                    self.write_hmm_seqfile(hmm_seqfname, query_name, second_query_name)
                    self.run_hmm(hmm_seqfname, query_name, second_query_name)

        os.remove(swoutfname)
        os.remove(hmm_seqfname)
        os.rmdir(self.workdir)
        # print self.fracs_in_common['v'] / self.n_pairs
        # print self.fracs_in_common['d'] / self.n_pairs
        # print self.fracs_in_common['j'] / self.n_pairs
        # print float(self.delta_k_v) / self.n_pairs
        # print float(self.delta_k_d) / self.n_pairs
                
    # # ----------------------------------------------------------------------------------------
    # def have_different_gene_matches(self, query_name, second_query_name):
    #     first_info = self.sw_info[query_name]
    #     second_info = self.sw_info[second_query_name]
    #     for region in utils.regions:
    #         matches = list((gene for gene in first_info['all'].split(':') if ('IGH' + region.upper()) in gene))
    #         second_matches = list((gene for gene in second_info['all'].split(':') if ('IGH' + region.upper()) in gene))
    #         print region,
    #         for im in range(min(len(matches), len(second_matches))):
    #             if matches[im] == second_matches[im]:
    #                 print '1',
    #             else:
    #                 print '0',
    #         denominator = min(len(matches), len(second_matches))
    #         n_in_common = len(set(matches).intersection(set(second_matches)))
    #         frac_in_common = float(n_in_common) / denominator
    #         self.fracs_in_common[region] += frac_in_common
    #         # print '   %s %d/%d = %.2f' % (region, n_in_common, denominator, frac_in_common),
    #         # if frac_in_common < 0.5:  # if, for any single region, the two queries share fewer than half of the matches, don't bother to run the hmm
    #         #     return True
    #     print '   %4d%4d' % (first_info['k_v'],second_info['k_v']),
    #     print '   %4d%4d' % (first_info['k_d'],second_info['k_d']),
    #     self.delta_k_v += abs(first_info['k_v'] - second_info['k_v'])
    #     self.delta_k_d += abs(first_info['k_d'] - second_info['k_d'])
    #     self.n_pairs += 1
    #     print '   ',
    #     return False
            
    # ----------------------------------------------------------------------------------------
    def get_score_index(self, query_name, second_query_name):
        """
        Return a hashable combination of the two query names that's the same if we reverse their order.
        At the moment, just add 'em up.
        """
        if second_query_name == '':
            second_query_name = '0'
        return int(query_name) + int(second_query_name)

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
            sys.exit()  # um, not sure which I want here, but it doesn't really matter. TODO kinda
            return matchlist

        hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        hmm_proc.wait()
        hmm_out, hmm_err = hmm_proc.communicate()
        if self.args.debug:
            print 'OUT\n',hmm_out
        if hmm_proc.returncode != 0:
            print 'aaarrrrrrrgggggh\n',hmm_err
            sys.exit()
    
        if self.args.algorithm == 'viterbi':
            # try:  # most likely reason for failure is something got kicked into stderr besides the csv info
            hmmreader = csv.DictReader(StringIO.StringIO(hmm_err))
            for hmmline in hmmreader:
                matchlist.append(hmmline)
            self.pair_scores[self.get_score_index(query_name, second_query_name)] = 0  # just to keep track of the fact that we did it already
            # except:
            #     print 'ERR\n',hmm_err
            #     sys.exit()
        elif self.args.algorithm == 'forward':
            total_score = float(hmm_err)  # TODO oh, right, I need to divide by the total prob of each individual sequence
            print '    total score: ',total_score
            assert total_score < 0.0 and total_score > -9999999.9
            self.pair_scores[self.get_score_index(query_name, second_query_name)] = total_score  # is that a hacky way to hash it so I can reverse the order and get the same entry? seems ok, I guess
            with opener('a')(self.pairscorefname) as pairscorefile:
                pairscorefile.write('%s,%s,%f\n' % (query_name, second_query_name, total_score))

        return matchlist
    
    # ----------------------------------------------------------------------------------------
    def write_specified_hmms(self, hmmdir, gene_list, v_right_length):
        for gene in gene_list:
            if len(re.findall('J[123]P', gene)) > 0:  # pretty sure these are crap
                print '  poof'
                continue
            writer = HmmWriter(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety,
                               hmmdir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=v_right_length)
            writer.write()

    # ----------------------------------------------------------------------------------------
    def write_all_hmms(self, v_right_length):
        if not os.path.exists(self.default_hmm_dir):
            os.makedirs(self.default_hmm_dir)
        for region in utils.regions:
            self.write_specified_hmms(self.default_hmm_dir, self.germline_seqs[region], v_right_length)

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, seqfname, query_name, second_query_name):
        start = time.time()
        # TODO dammit I'm still only using info from the first query
        all_gene_str = self.sw_info[query_name]['all']
        best_genes = self.sw_info[query_name]['best']
        k_v = self.sw_info[query_name]['k_v']
        k_d = self.sw_info[query_name]['k_d']
        v_fuzz = self.sw_info[query_name]['v_fuzz']
        d_fuzz = self.sw_info[query_name]['d_fuzz']
        v_right_length = self.sw_info[query_name]['v_right_length']

        assert abs(v_right_length - self.default_v_right_length) < 5  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this

        hmmdir = self.default_hmm_dir
        # write hmm files
        if self.args.write_hmms_on_fly:
            hmmdir = self.workdir
            self.write_specified_hmms(hmmdir, all_gene_str.split(':'), v_right_length)  # TODO note that this takes more time than actually running the hmm
    
        if self.args.algorithm == 'viterbi':
            print '\ninferred:'
    
        assert len(re.findall('J[123]P', best_genes['j'])) == 0  # I think we can remove this version (we never see it), but I'm putting a check in here just in case

        hmm_csv_infname = self.workdir + '/hmm_input.csv'

        # TODO since the s-w stuff excises the *best* v and *best* j, these k_v and k_d can be quite far off. Nonetheless seems ok a.t.m.
        # NOTE this call costs about 0.01-0.02 seconds if StochHMM.cpp does more or less nothing
        cmd_str = './stochhmm'
        cmd_str += ' -' + self.args.algorithm
        if self.args.pair:
            cmd_str += ' -hmmtype pair'
        else:
            cmd_str += ' -hmmtype single'
        cmd_str += ' -debug ' + str(self.args.debug)
        # cmd_str += ' -seq ' + seqfname
        # cmd_str += ' -k_v_guess ' + str(k_v) + ' -k_d_guess ' + str(k_d)
        # cmd_str += ' -v_fuzz ' + str(v_fuzz) + ' -d_fuzz ' + str(d_fuzz)
        # cmd_str += ' -only_genes \'' + all_gene_str + '\''  # hm, wait, do I need these escaped quotes?
        cmd_str += ' -hmmdir ' + hmmdir
        cmd_str += ' -infile ' + hmm_csv_infname

        with opener('w')(hmm_csv_infname) as hmm_csv_infile:
            header = ['k_v_guess', 'k_d_guess', 'v_fuzz', 'd_fuzz', 'only_genes', 'name', 'seq', 'second_name', 'second_seq']
            # writer = csv.DictWriter(hmm_csv_infile, header)
            # writer.writeheader()
            # row = {}  # TODO really wet here. see above
            # row['k_v_guess'] = k_v
            # row['k_d_guess'] = k_d
            # row['v_fuzz'] = v_fuzz
            # row['d_fuzz'] = d_fuzz
            # row['only_genes'] = all_gene_str
            # row['seq'] = self.reco_info[query_name]['seq']
            # row['name'] = query_name
            hmm_csv_infile.write(' '.join(header) + '\n')
            hmm_csv_infile.write('%d %d %d %d %s %s %s' % (k_v, k_d, v_fuzz, d_fuzz, all_gene_str, query_name, self.reco_info[query_name]['seq']))
            if second_query_name != '':
                # row['seq'] = self.reco_info[second_query_name]['second_seq']
                # row['second_name'] = second_query_name
                hmm_csv_infile.write(' %s %s\n' % (second_query_name, self.reco_info[second_query_name]['seq']))
            else:
                hmm_csv_infile.write(' %s %s\n' % ('x', 'x'))
                # row['seq'] = ''
                # row['second_name'] = ''
            # writer.writerow(row)

        matchlist = self.run_stochhmm(cmd_str, query_name, second_query_name)

        for matchline in matchlist:  # will be empty list if not viterbi
            utils.print_reco_event(self.germline_seqs, matchline, 0, 0)
            if self.args.pair:
                tmpseq = matchline['seq']  # TODO oh, man, that's a cludge
                matchline['seq'] = matchline['second_seq']
                utils.print_reco_event(self.germline_seqs, matchline, 0, 0, True)
                matchline['seq'] = tmpseq
                
    
        # remove hmm model files
        os.remove(hmm_csv_infname)
        for fname in os.listdir(self.workdir):
            if fname.endswith(".hmm"):
                os.remove(self.workdir + "/" + fname)
        print 'hmm run time: %.3f' % (time.time() - start)
    
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
        start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 3 --mismatch 1 --gap-open 100 ' + infname + ' ' + outfname + ' 2>/dev/null', shell=True)
        print 's-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_smith_waterman(self, infname):
        """
        Read bamfile output by s-w step, and write the info (after a bit of collation) to a csv.
        Note that the only reason to bother writing the csv is so you can avoid rerunning the s-w step every time you run the hmm.
        """
        gene_choice_probs = utils.read_overall_gene_prob(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety)
        with contextlib.closing(pysam.Samfile(infname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
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
                # if 'IGHJ4*' in best_genes['j'] and self.germline_seqs['d'][best_genes['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
                #     d_fuzz = 10

                assert query_name not in self.sw_info
                self.sw_info[query_name] = {}
                self.sw_info[query_name]['k_v'] = k_v
                self.sw_info[query_name]['k_d'] = k_d
                self.sw_info[query_name]['v_fuzz'] = v_fuzz
                self.sw_info[query_name]['d_fuzz'] = d_fuzz
                self.sw_info[query_name]['v_right_length'] = v_right_length
                self.sw_info[query_name]['best'] = best
                self.sw_info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

                assert k_v > 0 and k_d > 0 and v_fuzz > 0 and d_fuzz > 0 and v_right_length > 0
