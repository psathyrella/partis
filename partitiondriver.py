import time
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

import_start = time.time()
from utils import utils
from utils.opener import opener
from hmmwriter import HmmWriter
from clusterer import Clusterer
from waterer import Waterer

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, datadir, args, default_v_right_length=90, stochhmm_dir=''):
        self.datadir = datadir
        self.stochhmm_dir = stochhmm_dir
        self.args = args
        self.germline_seqs = utils.read_germlines(self.datadir)
        self.workdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())  # use a tmp dir specific to this process for the hmm input file
        self.dbgfname = self.workdir + '/dbg.txt'
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.default_v_fuzz = 2  # TODO play around with these default fuzzes
        self.default_d_fuzz = 2
        self.safety_buffer_that_I_should_not_need = 35  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this
        self.default_hmm_dir = self.stochhmm_dir + '/bcell/hmms/' + self.args.human + '/' + self.args.naivety
        self.default_v_right_length = default_v_right_length
        self.precluster_info = {}
        self.match_names_for_all_queries = set()

        self.is_data = False  # data or simulation
        self.input_info = {}
        self.reco_info = {}  # generator truth information
        self.read_input_file()  # read simulation info and write sw input file
        self.sw_info = {}

    # ----------------------------------------------------------------------------------------
    def from_same_event(self, query_name, second_query_name):
        if self.args.pair:
            return self.reco_info[query_name]['reco_id'] == self.reco_info[second_query_name]['reco_id']
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def read_input_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        with opener('r')(self.args.simfile) as simfile:
            last_reco_id = -1  # only run on the first seq in each reco event. they're all pretty similar
            reader = csv.DictReader(simfile)
            n_queries = 0
            for line in reader:
                # if command line specified query or reco ids, skip other ones
                if self.args.queries != None and line['unique_id'] not in self.args.queries:
                    continue
                if self.args.reco_ids != None and line['reco_id'] not in self.args.reco_ids:
                    continue

                self.input_info[line['unique_id']] = {'unique_id':line['unique_id'], 'seq':line['seq']}
                if not self.is_data:
                    self.reco_info[line['unique_id']] = line
                last_reco_id = line['reco_id']
                n_queries += 1
                if self.args.n_total_queries > 0 and n_queries >= self.args.n_total_queries:
                    break
        assert len(self.input_info) > 0
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        with opener('w')(self.dbgfname) as dbgfile:
            dbgfile.write('hamming\n')

        # run smith-waterman
        waterer = Waterer(self)
        self.sw_info = waterer.sw_info
        # if not self.args.skip_sw:
        #     self.run_smith_waterman(swoutfname)
        # else:
        #     assert False  # doesn't actually work swoutfname = swoutfname.replace(os.path.dirname(swoutfname), '.')
        # self.read_smith_waterman(swoutfname)  # read sw bam output, collate it, and write to csv for hmm input

        # # cdr3 length partitioning
        # cdr3_length_clusters = self.cdr3_length_precluster()

        # hamming preclustering
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
            self.read_hmm_output(stripped_hmm_csv_outfname, stripped_pairscorefname);
            stripped_clusters = Clusterer(0, greater_than=True)  # TODO better way to set threshhold?
            stripped_clusters.cluster(stripped_pairscorefname, debug=False, reco_info=self.reco_info)
            for query_name in self.sw_info:  # check for singletons that got split out in the preclustering step
                if query_name not in stripped_clusters.query_clusters:
                    print '    singleton ',query_name
            if not self.args.no_clean:
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
        self.read_hmm_output(hmm_csv_outfname, pairscorefname);

        clusters = Clusterer(0, greater_than=True)  # TODO better way to set threshhold?
        clusters.cluster(pairscorefname, debug=False, reco_info=self.reco_info)
        for query_name in self.sw_info:  # check for singletons that got split out in the preclustering step
            if query_name not in clusters.query_clusters:
                print '    singleton ',query_name

        if not self.args.no_clean:
            os.remove(hmm_csv_infname)
            os.remove(hmm_csv_outfname)
            os.remove(pairscorefname)
            os.rename(self.dbgfname, '/tmp/dralph/debug-' + self.args.human + '.txt')
            os.rmdir(self.workdir)

    # # ----------------------------------------------------------------------------------------
    # def cdr3_length_precluster(self):
    #     cdr3_length_infname = self.workdir + '/cdr3_length_precluster_input.csv'
    #     with opener('w')(cdr3_length_infname) as infile:
    #         columns = ('unique_id', 'seq')
    #         # ARGGGG don't use reco_info for this... segregate out the simulation info
    #         csv.DictWriter()
        
    #     connordir = '/shared/silo_researcher/Matsen_F/MatsenGrp/working/cmccoy/20140414-bcell-validation'
    #     cmd_str = connordir + '/venv/bin/linsim dists-by-cdr3'
    #     cmd_str += ' --name-column unique_id'
    #     cmd_str += ' --sequence-column seq'
    #     cmd_str += cdr3_length_infname
    #     cmd_str += '/tmp/out'
    #     hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
    #     hmm_proc.wait()
    #     hmm_out, hmm_err = hmm_proc.communicate()

    #     if not self.no_clean:
    #         os.remove(cdr3_length_infname)

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
            for query_name in self.input_info:
                for second_query_name in self.input_info:
                    if second_query_name == query_name:
                        continue
                    if utils.get_key(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                        continue
                    pair_scores[utils.get_key(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                    n_lines += 1
                    query_seq = self.input_info[query_name]['seq']
                    second_query_seq = self.input_info[second_query_name]['seq']
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
            for query_name in self.input_info:
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
                    for second_query_name in self.input_info:  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
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
                                       self.input_info[query_name]['seq'], self.input_info[second_query_name]['seq']))
                else:
                    # assert self.args.algorithm == 'viterbi'  # TODO hm, was that really all I had to do to allow non-pair forward?
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, k_v_min, k_v_max, k_d_min, k_d_max, only_genes, self.input_info[query_name]['seq']))
        print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)

    # ----------------------------------------------------------------------------------------
    def run_stochhmm(self, csv_infname, csv_outfname):
        start = time.time()

        # build the command line
        cmd_str = self.stochhmm_dir + '/stochhmm'
        cmd_str += ' --algorithm ' + self.args.algorithm
        if self.args.pair:
            cmd_str += ' --pair 1'
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + self.default_hmm_dir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        # cmd_str += ' >/dev/null'

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
    def read_hmm_output(self, hmm_csv_outfname, pairscorefname):
        # TODO the input and output files for this function are almost identical at this point

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
