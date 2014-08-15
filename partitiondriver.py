import time
import_start = time.time()
import sys
import os
import math
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
from clusterer import Clusterer

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, datadir, args, default_v_right_length=90):
        self.datadir = datadir
        self.args = args
        self.germline_seqs = utils.read_germlines()
        self.workdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())  # use a tmp dir specific to this process for the hmm input file
        self.dbgfname = self.workdir + '/dbg.txt'
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.reco_info = self.read_sim_file()  # read simulation info and write sw input file
        self.default_v_fuzz = 2  # TODO play around with these default fuzzes
        self.default_d_fuzz = 2
        self.safety_buffer_that_I_should_not_need = 35  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this
        self.default_hmm_dir = 'bcell/hmms/' + self.args.human + '/' + self.args.naivety
        self.default_v_right_length = default_v_right_length
        self.sw_info = {}
        self.precluster_info = {}
        self.match_names_for_all_queries = set()
            
    # ----------------------------------------------------------------------------------------
    def from_same_event(self, query_name, second_query_name):
        if self.args.pair:
            return self.reco_info[query_name]['reco_id'] == self.reco_info[second_query_name]['reco_id']
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def read_sim_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        reco_info = {}  # generator truth information
        with opener('r')(self.args.simfile) as simfile:
            last_reco_id = -1  # only run on the first seq in each reco event. they're all pretty similar
            reader = csv.DictReader(simfile)
            n_queries = 0
            for line in reader:
                if self.args.queries != None and line['unique_id'] not in self.args.queries:
                    continue
                if self.args.reco_ids != None and line['reco_id'] not in self.args.reco_ids:
                    continue
                reco_info[line['unique_id']] = line
                last_reco_id = line['reco_id']
                n_queries += 1
                if self.args.n_total_queries > 0 and n_queries >= self.args.n_total_queries:
                    break
        assert len(reco_info) > 0
        return reco_info
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        # run smith-waterman
        swoutfname = self.workdir + '/query-seqs.bam'
        if not self.args.skip_sw:
            self.run_smith_waterman(swoutfname)
        else:
            assert False  # doesn't actually work swoutfname = swoutfname.replace(os.path.dirname(swoutfname), '.')
        self.read_smith_waterman(swoutfname)  # read sw bam output, collate it, and write to csv for hmm input

        # hamming preclustering
        with opener('w')(self.dbgfname) as dbgfile:
            dbgfile.write('hamming\n')
        hamming_clusters = self.hamming_precluster()

        stripped_clusters = None
        if self.args.algorithm == 'forward':
            # run a stripped-down hmm to precluster (one gene version per query, no v or d fuzz). The non-related seqs creep up toward scores above zero, but pairs that *are* related are super good about staying at large positive score when you strip down the HMM (because their scores are dominated by the *correc* choice)
            stripped_hmm_csv_infname = self.workdir + '/stripped_hmm_input.csv'
            stripped_hmm_csv_outfname = self.workdir + '/stripped_hmm_output.csv'
            stripped_pairscorefname = self.workdir + '/stripped-hmm-pairscores.csv'
            print 'stripped hmm'
            with opener('a')(self.dbgfname) as dbgfile:
                dbgfile.write('stripped hmm\n')
            self.write_hmm_input(stripped_hmm_csv_infname, preclusters=hamming_clusters, stripped=True)
            self.run_stochhmm(stripped_hmm_csv_infname, stripped_hmm_csv_outfname)
            self.read_stochhmm_output(stripped_hmm_csv_outfname, stripped_pairscorefname);
            stripped_clusters = Clusterer(0, greater_than=True)  # TODO better way to set threshhold?
            stripped_clusters.cluster(stripped_pairscorefname, debug=False)
            for query_name in self.sw_info:  # check for singletons that got split out in the preclustering step
                if query_name not in stripped_clusters.query_clusters:
                    print '    singleton ',query_name
            if self.args.no_clean:
                os.remove(stripped_hmm_csv_infname)
                os.remove(stripped_hmm_csv_outfname)
                os.remove(stripped_pairscorefname)

        hmm_csv_infname = self.workdir + '/hmm_input.csv'
        hmm_csv_outfname = self.workdir + '/hmm_output.csv'
        pairscorefname = self.workdir + '/hmm-pairscores.csv'
        print 'hmm'
        with opener('a')(self.dbgfname) as dbgfile:
            dbgfile.write('hmm\n')
        self.write_hmm_input(hmm_csv_infname, preclusters=stripped_clusters)
        self.run_stochhmm(hmm_csv_infname, hmm_csv_outfname)
        self.read_stochhmm_output(hmm_csv_outfname, pairscorefname);

        clusters = Clusterer(0, greater_than=True)  # TODO better way to set threshhold?
        clusters.cluster(pairscorefname, debug=False)
        for query_name in self.sw_info:  # check for singletons that got split out in the preclustering step
            if query_name not in clusters.query_clusters:
                print '    singleton ',query_name

        if self.args.no_clean:
            os.remove(hmm_csv_infname)
            os.remove(hmm_csv_outfname)
            os.remove(pairscorefname)
            os.rename(self.dbgfname, '/tmp/dralph/debug-' + self.args.human + '.txt')
            os.rmdir(self.workdir)

    # ----------------------------------------------------------------------------------------
    def hamming_precluster(self):
        if not self.args.algorithm == 'forward':
            return
        start = time.time()
        print 'hamming clustering'
        hammingfname = self.workdir + '/fractional-hamming-scores.csv'
        pair_scores = {}
        with opener('w')(hammingfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'score'))
            writer.writeheader()
            n_lines = 0
            for query_name in self.reco_info:
                for second_query_name in self.reco_info:
                    if second_query_name == query_name:
                        continue
                    if utils.get_key(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                        continue
                    pair_scores[utils.get_key(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                    n_lines += 1
                    query_seq = self.reco_info[query_name]['seq']
                    second_query_seq = self.reco_info[second_query_name]['seq']
                    mutation_frac = utils.hamming(query_seq, second_query_seq) / float(len(query_seq))
                    writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'score':mutation_frac})
                    if self.args.debug:
                        print '    %20s %20s %8.2f' % (query_name, second_query_name, mutation_frac)
                    with opener('a')(self.dbgfname) as dbgfile:
                        dbgfile.write('%20s %20s   %d  %f\n' % (query_name, second_query_name,  self.from_same_event(query_name, second_query_name), mutation_frac))

        clust = Clusterer(0.5, greater_than=False)  # TODO this 0.5 number isn't gonna be the same if the query sequences change length
        clust.cluster(hammingfname, debug=False)
        os.remove(hammingfname)
        print '    %d lines' % n_lines
        print '    hamming time: %.3f' % (time.time()-start)
        return clust

    # ----------------------------------------------------------------------------------------
    def run_smith_waterman(self, outfname):
        """
        Run smith-waterman alignment on the seqs in <infname>, and toss all the top matches into <outfname>.
        Then run through <outfname> to get the top hits and their locations to pass to the hmm.
        Then run the hmm on each gene set.
        """
        infname = self.workdir + '/query-seqs.fa'
        with opener('w')(infname) as swinfile:  # write *all* the input seqs to file, i.e. run s-w on all of 'em at once
            for query_name,line in self.reco_info.iteritems():
                swinfile.write('>' + query_name + ' NUKES\n')
                swinfile.write(line['seq'] + '\n')
        start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 5 --mismatch 3 --gap-open 100 ' + infname + ' ' + outfname + ' 2>/dev/null', shell=True)
        os.remove(infname)
        print '    s-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_smith_waterman(self, swinfoname):
        """
        Read bamfile output by s-w step, and write the info (after a bit of collation) to a csv.
        Note that the only reason to bother writing the csv is so you can avoid rerunning the s-w step every time you run the hmm.
        """
        start = time.time()
        gene_choice_probs = utils.read_overall_gene_prob(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety)
        with contextlib.closing(pysam.Samfile(swinfoname)) as bam:
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
                    # TODO the s-w allows the j right edge to be chopped off -- I should skip the matches where different amounts are chopped off in the query and germline
                    all_germline_bounds[gene] = (read.pos, read.aend)
    
                best = {}
                match_names = {}
                n_matches = {'v':0, 'd':0, 'j':0}
                n_used = {'v':0, 'd':0, 'j':0}
                k_d_min = 999
                k_d_max = 0
                k_v_min = 999
                k_v_max = 0
                for region in utils.regions:
                    all_match_names[region] = sorted(all_match_names[region], reverse=True)
                    match_names[region] = []
                for region in utils.regions:
                    for score,gene in all_match_names[region]:
                        n_matches[region] += 1
                        if n_matches[region] > self.args.n_max_per_region:  # only take the top few from each region
                            # TODO should use *lots* of d matches, but fewer vs and js
                            continue
                        n_used[region] += 1
                        match_names[region].append(gene)

                        if region == 'v':
                            this_k_v = all_query_bounds[gene][1]
                            k_v_min = min(this_k_v, k_v_min)
                            k_v_max = max(this_k_v, k_v_max)
                        if region == 'd':
                            this_k_d = all_query_bounds[gene][1] - all_query_bounds[raw_best['v']][1]  # end of d minus end of v
                            k_d_min = min(this_k_d, k_d_min)
                            k_d_max = max(this_k_d, k_d_max)
    
                        # check consistency with best match (since the best match is excised in s-w code, and because stochhmm is run with *one* k_v k_d set)
                        if region not in best:
                            best[region] = gene
    
                        if self.args.debug:
                            buff_str = (17 - len(gene)) * ' '
                            print '%8s%s%s%6.1e * %2.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, gene_choice_probs[region][gene], score / gene_choice_probs[region][gene], score),
                            glbounds = all_germline_bounds[gene]
                            qrbounds = all_query_bounds[gene]
                            print ' %4d%4d   %s' % (glbounds[0], glbounds[1], self.germline_seqs[region][gene][glbounds[0]:glbounds[1]])
                            print '%46s  %4d%4d   %s' % ('', qrbounds[0], qrbounds[1], utils.color_mutants(self.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]))
                                
                # print how many of the available matches we used
                if self.args.debug:
                    print '  used',
                    for region in utils.regions:
                        if region != 'v':
                            print '      ',
                        print ' %d / %d in %s' % (n_used[region], n_matches[region], region)

                for region in utils.regions:
                    if region not in best:
                        print ' no',region,'match found for',query_name
                        print '    true:'
                        utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], 0, 0, extra_str='    ')
                # write sw info to file for hmm
                # best k_v, k_d:
                k_v = all_query_bounds[best['v']][1]  # end of v match
                k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

                if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
                    # if self.args.debug:
                    print '  expanding k_d'
                    k_d_max = max(8, k_d_max)
                    
                v_right_length = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match
                if 'IGHJ4*' in best['j'] and self.germline_seqs['d'][best['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
                    # print '  doubly expanding k_d'
                    if k_d_max-k_d_min < 8:
                        k_d_min -= 5
                        k_d_max += 2

                k_v_min = max(0, k_v_min - self.default_v_fuzz)  # ok, so I don't *actually* want it to be zero... oh, well
                k_v_max += self.default_v_fuzz
                k_d_min = max(1, k_d_min - self.default_d_fuzz)
                k_d_max += self.default_d_fuzz
                assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0 and v_right_length > 0

                # print '  k_v min %d max %d (best %d)' % (k_v_min, k_v_max, k_v)
                # print '  k_d min %d max %d (best %d)' % (k_d_min, k_d_max, k_d)
                self.match_names_for_all_queries |= set(match_names['v'] + match_names['d'] + match_names['j'])

                assert query_name not in self.sw_info
                self.sw_info[query_name] = {}
                self.sw_info[query_name]['k_v'] = k_v
                self.sw_info[query_name]['k_d'] = k_d
                self.sw_info[query_name]['k_v_min'] = k_v_min
                self.sw_info[query_name]['k_v_max'] = k_v_max
                self.sw_info[query_name]['k_d_min'] = k_d_min
                self.sw_info[query_name]['k_d_max'] = k_d_max
                self.sw_info[query_name]['v_right_length'] = v_right_length
                self.sw_info[query_name]['best'] = best
                self.sw_info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])
        print '    sw read time: %.3f' % (time.time() - start)
        os.remove(swinfoname)
    
    # ----------------------------------------------------------------------------------------
    def write_specified_hmms(self, hmmdir, gene_list, v_right_length):
        for gene in gene_list:
            if len(re.findall('J[123]P', gene)) > 0:  # pretty sure these versions are crap
                print '  poof'
                continue
            if 'OR16' in gene or 'OR21' in gene or 'V3/OR15' in gene:
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
    def write_hmm_input(self, csv_fname, preclusters=None, stripped=False):  # TODO use different input files for the two hmm steps
        with opener('w')(csv_fname) as csvfile:
            # write header
            header = ['name', 'second_name', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seq', 'second_seq']  # I wish I had a good c++ csv reader 
            csvfile.write(' '.join(header) + '\n')

            # then write a line for each query sequence (or pair of them)
            pair_scores = {}
            n_lines, n_preclustered, n_previously_preclustered, n_removable, n_singletons = 0, 0, 0, 0, 0
            for query_name in self.reco_info:
                info = self.sw_info[query_name]

                # make sure we don't need to rewrite the hmm model files
                if abs(info['v_right_length'] - self.default_v_right_length) >= self.safety_buffer_that_I_should_not_need:
                    print 'WARNING VRIGHT ', info['v_right_length'], self.default_v_right_length, (info['v_right_length']-self.default_v_right_length)
                # assert abs(info['v_right_length'] - self.default_v_right_length) < self.safety_buffer_that_I_should_not_need

                # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                assert len(re.findall('J[123]P', info['best']['j'])) == 0

                only_genes = info['all']
                k_v_min = info['k_v_min']
                k_v_max = info['k_v_max']
                k_d_min = info['k_d_min']
                k_d_max = info['k_d_max']
                if self.args.pair:
                    for second_query_name in self.reco_info:  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
                        if second_query_name == query_name:
                            continue
                        if utils.get_key(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                            continue
                        pair_scores[utils.get_key(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                        if preclusters != None:  # if we've already run preclustering, skip the pairs that we know aren't matches
                            if query_name not in preclusters.query_clusters or second_query_name not in preclusters.query_clusters:  # singletons!
                                n_singletons += 1
                                continue
                            if utils.get_key(query_name, second_query_name) not in preclusters.pairscores:
                                n_previously_preclustered += 1
                                continue
                            if preclusters.query_clusters[query_name] != preclusters.query_clusters[second_query_name]:  # not in same cluster
                                n_preclustered += 1
                                continue
                            if preclusters.is_removable(preclusters.pairscores[utils.get_key(query_name, second_query_name)]):  # in same cluster, but score (link) is long. i.e. *this* pair is far apart, but other seqs to which they are linked are close to each other
                                n_removable += 1
                                continue

                        second_info = self.sw_info[second_query_name]
                        second_only_genes = second_info['all']
                        # TODO the current method of expanding the k space bounds to encompass all the gene matches, and both sequences, is *correct*, but I think it's unnecessarily slow
                        k_v_min = min(k_v_min, second_info['k_v_min'])
                        k_v_max = max(k_v_max, second_info['k_v_max'])
                        k_d_min = min(k_d_min, second_info['k_d_min'])
                        k_d_max = max(k_d_max, second_info['k_d_max'])
                        if stripped:  # strip down the hmm -- only use the single best gene for each sequence, and don't fuzz at all
                            only_genes = info['best']['v'] + ':' + info['best']['d'] + ':' + info['best']['j']
                            second_only_genes = second_info['best']['v'] + ':' + second_info['best']['d'] + ':' + second_info['best']['j']
                            k_v_min = info['k_v']
                            k_v_max = k_v_min + 1
                            k_d_min = info['k_d']
                            k_d_max = k_d_min + 1

                        only_genes = ':'.join(set(only_genes.split(':')) | set(second_only_genes.split(':')))  # NOTE using both sets of genes (from both query seqs) like this *really* helps, 
                        n_lines += 1
                        csvfile.write('%s %s %d %d %d %d %s %s %s\n' %
                                      (query_name, second_query_name, k_v_min, k_v_max, k_d_min, k_d_max, only_genes,
                                       self.reco_info[query_name]['seq'], self.reco_info[second_query_name]['seq']))
                else:
                    # assert self.args.algorithm == 'viterbi'  # TODO hm, was that really all I had to do to allow non-pair forward?
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, k_v_min, k_v_max, k_d_min, k_d_max, only_genes, self.reco_info[query_name]['seq']))
        print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)

    # ----------------------------------------------------------------------------------------
    def run_stochhmm(self, csv_infname, csv_outfname):
        start = time.time()

        # build the command line
        cmd_str = './stochhmm'
        cmd_str += ' --algorithm ' + self.args.algorithm
        if self.args.pair:
            cmd_str += ' --pair 1'
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + self.default_hmm_dir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname

        if self.args.debug == 2:  # not sure why, but popen thing hangs with debug 2
            check_call(cmd_str, shell=True)  # um, not sure which I want here, but it doesn't really matter. TODO kinda
            sys.exit()

        # run hmm
        hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        hmm_proc.wait()
        hmm_out, hmm_err = hmm_proc.communicate()
        print hmm_out,
        print hmm_err,
        if hmm_proc.returncode != 0:
            print 'aaarrrrrrrgggggh\n', hmm_out, hmm_err
            sys.exit()

        print '    hmm run time: %.3f' % (time.time() - start)

    # ----------------------------------------------------------------------------------------
    def read_stochhmm_output(self, hmm_csv_outfname, pairscorefname):
        """ Read hmm output file """
        # write header for pairwise score file
        with opener('w')(pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id,second_unique_id,score\n')
        with opener('r')(hmm_csv_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            last_id = None
            n_boundary_errors = 0
            for line in reader:
                boundary_error = False
                if 'boundary' in line['errors'].split(':'):
                    n_boundary_errors += 1
                    boundary_error = True
                else:
                    assert len(line['errors']) == 0
                if self.args.algorithm == 'viterbi':
                    # if this is the first line for this query (or query pair), print the true event
                    if last_id != utils.get_key(line['unique_id'], line['second_unique_id']):
                        print '%20s %20s   %d' % (line['unique_id'], line['second_unique_id'], self.from_same_event(line['unique_id'], line['second_unique_id']))
                        print '    true:'
                        utils.print_reco_event(self.germline_seqs, self.reco_info[line['unique_id']], 0, 0, extra_str='    ')
                        if self.args.pair:
                            # print '      and maybe'
                            utils.print_reco_event(self.germline_seqs, self.reco_info[line['second_unique_id']], 0, 0, self.from_same_event(line['unique_id'], line['second_unique_id']), '    ')
                        print '    inferred:'
                    utils.print_reco_event(self.germline_seqs, line, 0, 0, extra_str='    ')
                    if self.args.pair:
                        tmpseq = line['seq']  # TODO oh, man, that's a cludge
                        line['seq'] = line['second_seq']
                        utils.print_reco_event(self.germline_seqs, line, 0, 0, True, extra_str='    ')
                        line['seq'] = tmpseq
                else:
                    with opener('a')(self.dbgfname) as dbgfile:
                        dbgfile.write('%20s %20s   %d  %7.2f  %d\n' %
                                      (line['unique_id'], line['second_unique_id'], self.from_same_event(line['unique_id'], line['second_unique_id']), float(line['score']), boundary_error))
                    with opener('a')(pairscorefname) as pairscorefile:
                        pairscorefile.write('%s,%s,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))
                last_id = utils.get_key(line['unique_id'], line['second_unique_id'])
        if n_boundary_errors > 0:
            print '    %d boundary errors' % n_boundary_errors
