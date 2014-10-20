import time
import sys
import json
import os
import glob
import csv
import re
from subprocess import Popen, check_call, PIPE

import_start = time.time()
import utils
from opener import opener
from hmmwriter import HmmWriter
from clusterer import Clusterer
from waterer import Waterer
from parametercounter import ParameterCounter
import plotting

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
def from_same_event(is_data, pair_hmm, reco_info, query_name, second_query_name):
    if is_data:
        return False
    if pair_hmm:
        return reco_info[query_name]['reco_id'] == reco_info[second_query_name]['reco_id']
    else:
        return False

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, args):
        self.args = args
        self.germline_seqs = utils.read_germlines(self.args.datadir)

        self.safety_buffer_that_I_should_not_need = 35  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this

        with opener('r')(self.args.datadir + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with opener('r')(self.args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

        self.precluster_info = {}
        self.input_info = {}
        self.reco_info = None
        if not self.args.is_data:
            self.reco_info = {}  # generator truth information

        self.read_input_file()  # read simulation info and write sw input file

    # ----------------------------------------------------------------------------------------
    def clean(self, waterer):
        if self.args.no_clean:
            return
        if not self.args.read_cached_parameters:  # ha ha I don't have this arg any more
            waterer.clean()  # remove all the parameter files that the waterer's ParameterCounter made
            for fname in glob.glob(waterer.parameter_dir + '/hmms/*.yaml'):  # remove all the hmm files from the same directory
                os.remove(fname)
            os.rmdir(waterer.parameter_dir + '/hmms')
            os.rmdir(waterer.parameter_dir)

        # check to see if we missed anything
        leftovers = [ fname for fname in os.listdir(self.args.workdir) ]
        for over in leftovers:
            print self.args.workdir + '/' + over
        os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def read_input_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        with opener('r')(self.args.seqfile) as seqfile:
            delimiter = ','
            name_column = 'unique_id'
            seq_column = 'seq'
            # n_nukes_skip = 0  # skip all the TTTTTTTTs in data. TODO er that seems hackey doesn't it?
            if '.tsv' in self.args.seqfile:
                delimiter = '\t'
                name_column = 'name'
                seq_column = 'nucleotide'
                # n_nukes_skip = 9  # TODO ehhhhh I don't like this
            reader = csv.DictReader(seqfile, delimiter=delimiter)
            n_queries = 0
            for line in reader:
                # if command line specified query or reco ids, skip other ones
                if self.args.queries != None and line[name_column] not in self.args.queries:
                    continue
                if self.args.reco_ids != None and line['reco_id'] not in self.args.reco_ids:
                    continue

                self.input_info[line[name_column]] = {'unique_id':line[name_column], 'seq':line[seq_column][self.args.n_bases_skip:]}
                if not self.args.is_data:
                    self.reco_info[line['unique_id']] = line
                n_queries += 1
                if self.args.n_max_queries > 0 and n_queries >= self.args.n_max_queries:
                    break
        assert len(self.input_info) > 0
    
    # ----------------------------------------------------------------------------------------
    def cache_parameters(self, parameter_dir=''):

        sw_plotdir, hmm_plotdir = '', ''
        if os.getenv('www') != None:
            sw_plotdir = os.getenv('www') + '/partis/sw_parameters'
            hmm_plotdir = os.getenv('www') + '/partis/hmm_parameters'

        if parameter_dir == '':
            parameter_dir = self.args.parameter_dir

        sw_parameter_dir = parameter_dir + '/sw_parameters'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, from_scratch=True, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=sw_plotdir)
        waterer.run()

        print 'writing hmms with sw info'
        self.write_hmms(sw_parameter_dir, waterer.info['all_best_matches'])

        hmm_parameter_dir = parameter_dir + '/hmm_parameters'
        assert not self.args.pair  # er, could probably do it for forward pair, but I'm not really sure what it would mean
        self.run_hmm('viterbi', waterer.info, parameter_in_dir=sw_parameter_dir, parameter_out_dir=hmm_parameter_dir, count_parameters=True, plotdir=hmm_plotdir)

        print 'writing hmms with hmm info'
        self.write_hmms(hmm_parameter_dir, waterer.info['all_best_matches'])

        # self.clean(waterer)  # TODO get this working again. *man* it's a total bitch keeping track of all the files I'm writing

    # ----------------------------------------------------------------------------------------
    def partition(self):

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, from_scratch=False, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        # cdr3 length partitioning
        cdr3_cluster = False  # don't precluster on cdr3 length for the moment -- I cannot accurately infer cdr3 length in some sequences, so I need a way to pass query seqs to the clusterer with several possible cdr3 lengths. TODO fix that
        cdr3_length_clusters = None
        if cdr3_cluster:
            cdr3_length_clusters = self.cdr3_length_precluster(waterer)

        assert self.args.pair
        hamming_clusters = self.hamming_precluster(cdr3_length_clusters)
        stripped_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hamming_clusters, stripped=True)
        final_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=stripped_clusters, stripped=False)

        # self.clean(waterer)

    # ----------------------------------------------------------------------------------------
    def point_estimate(self):
        """ get a point estimate for each query sequence (or pair). i.e. run viterbi """
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, from_scratch=False, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.run_hmm('viterbi', waterer.info, parameter_in_dir=self.args.parameter_dir, plot_performance=self.args.plot_performance)
        # self.clean(waterer)  # TODO get this working again. *man* it's a total bitch keeping track of all the files I'm writing

        if self.args.plot_performance:
            plotting.compare_directories('/var/www/sharing/dralph/partis/performance/plots/', 'hmm', '/var/www/sharing/dralph/partis/sw_performance/plots/', 'smith-water', xtitle='inferred - true', stats='rms')

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, sw_info, parameter_in_dir, parameter_out_dir='', preclusters=None, stripped=False, prefix='', count_parameters=False, plotdir='', plot_performance=False):
        pcounter = None
        if count_parameters:
            pcounter = ParameterCounter(self.germline_seqs, parameter_out_dir, plotdir=plotdir)
        if plotdir != '':
            utils.prep_dir(plotdir + '/plots', '*.svg')

        perfplotter = None
        if plot_performance:
            assert not self.args.is_data
            assert not self.args.pair  # well, you *could* check the performance of pair viterbi, but I haven't implemented it yet
            assert algorithm == 'viterbi'  # same deal as previous assert
            from performanceplotter import PerformancePlotter
            perfplotter = PerformancePlotter(self.germline_seqs, os.getenv('www') + '/partis/performance', 'hmm')

        if prefix == '' and stripped:
            prefix = 'stripped'
        print prefix,'hmm'
        csv_infname = self.args.workdir + '/' + prefix + '_hmm_input.csv'
        csv_outfname = self.args.workdir + '/' + prefix + '_hmm_output.csv'
        pairscorefname = self.args.workdir + '/' + prefix + '_hmm_pairscores.csv'
        self.write_hmm_input(csv_infname, sw_info, preclusters=preclusters, stripped=stripped, parameter_dir=parameter_in_dir)
        self.run_hmm_binary(algorithm, csv_infname, csv_outfname, parameter_dir=parameter_in_dir)
        self.read_hmm_output(algorithm, csv_outfname, pairscorefname, pcounter, perfplotter)

        if count_parameters:
            pcounter.write_counts()

        clusters = None
        if self.args.pair and algorithm == 'forward':
            clusters = Clusterer(0, greater_than=True)
            clusters.cluster(pairscorefname, debug=False, reco_info=self.reco_info)
            if preclusters != None:
                for query_name in sw_info:  # check for singletons that got split out in the preclustering step
                    if query_name not in clusters.query_clusters and query_name != 'all_best_matches':
                        print '    singleton ', query_name

        if not self.args.no_clean:
            os.remove(csv_infname)
            os.remove(csv_outfname)
            os.remove(pairscorefname)

        return clusters

    # ----------------------------------------------------------------------------------------
    def cdr3_length_precluster(self, waterer, preclusters=None):
        cdr3lengthfname = self.args.workdir + '/cdr3lengths.csv'
        with opener('w')(cdr3lengthfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'cdr3_length', 'second_cdr3_length', 'score'))
            writer.writeheader()
            for query_name, second_query_name in self.get_pairs(preclusters):
                cdr3_length = waterer.info[query_name]['cdr3_length']
                second_cdr3_length = waterer.info[second_query_name]['cdr3_length']
                same_length = cdr3_length == second_cdr3_length
                if not self.args.is_data:
                    assert cdr3_length == int(self.reco_info[query_name]['cdr3_length'])
                    if second_cdr3_length != int(self.reco_info[second_query_name]['cdr3_length']):
                        print 'WARNING did not infer correct cdr3 length'
                        assert False
                writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'cdr3_length':cdr3_length, 'second_cdr3_length':second_cdr3_length, 'score':int(same_length)})

        clust = Clusterer(0.5, greater_than=True)  # i.e. cluster together if same_length == True
        clust.cluster(cdr3lengthfname, debug=False)
        os.remove(cdr3lengthfname)
        return clust

    # ----------------------------------------------------------------------------------------
    def get_pairs(self, preclusters=None):
        """ Get all unique the pairs of sequences in input_info, skipping where preclustered out """
        pair_scores = {}
        pairs = []
        n_lines, n_preclustered, n_previously_preclustered, n_removable, n_singletons = 0, 0, 0, 0, 0
        for query_name in self.input_info:
            for second_query_name in self.input_info:
                if second_query_name == query_name:
                    continue
                if utils.get_key(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                    continue
                pair_scores[utils.get_key(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                if preclusters != None:  # if we've already run preclustering, skip the pairs that we know aren't matches
                    if query_name not in preclusters.query_clusters or second_query_name not in preclusters.query_clusters:  # singletons (i.e. they were already preclustered into their own group)
                        n_singletons += 1
                        continue
                    if utils.get_key(query_name, second_query_name) not in preclusters.pairscores:  # preclustered out in a previous preclustering step
                        n_previously_preclustered += 1
                        continue
                    if preclusters.query_clusters[query_name] != preclusters.query_clusters[second_query_name]:  # not in same cluster
                        n_preclustered += 1
                        continue
                    if preclusters.is_removable(preclusters.pairscores[utils.get_key(query_name, second_query_name)]):  # in same cluster, but score (link) is long. i.e. *this* pair is far apart, but other seqs to which they are linked are close to each other
                        n_removable += 1
                        continue
                pairs.append((query_name, second_query_name))
                n_lines += 1

        print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)
        return pairs

    # ----------------------------------------------------------------------------------------
    def hamming_precluster(self, preclusters=None):
        start = time.time()
        print 'hamming clustering'
        hammingfname = self.args.workdir + '/fractional-hamming-scores.csv'
        with opener('w')(hammingfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'score'))
            writer.writeheader()
            for query_name, second_query_name in self.get_pairs(preclusters):
                query_seq = self.input_info[query_name]['seq']
                second_query_seq = self.input_info[second_query_name]['seq']
                mutation_frac = utils.hamming(query_seq, second_query_seq) / float(len(query_seq))
                writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'score':mutation_frac})
                if self.args.debug:
                    print '    %20s %20s %8.2f' % (query_name, second_query_name, mutation_frac)

        clust = Clusterer(0.5, greater_than=False)  # TODO this 0.5 number isn't gonna be the same if the query sequences change length
        clust.cluster(hammingfname, debug=False)
        os.remove(hammingfname)
        print '    hamming time: %.3f' % (time.time()-start)
        return clust
    
    # ----------------------------------------------------------------------------------------
    def write_hmms(self, parameter_dir, sw_matches):
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')

        gene_list = self.args.only_genes
        if gene_list == None:  # if specific genes weren't specified, do the ones for which we have matches
            gene_list = []
            for region in utils.regions:
                for gene in self.germline_seqs[region]:
                    if sw_matches == None or gene in sw_matches:  # shouldn't be None really, but I'm testing something
                        gene_list.append(gene)

        for gene in gene_list:
            # if self.args.debug:
            print '   ', utils.color_gene(gene),
            sys.stdout.flush()
            # if len(re.findall('J[123]P', gene)) > 0 or 'OR16' in gene or 'OR21' in gene or 'OR15' in gene or 'IGHV1-46*0' in gene or 'IGHV1-68' in gene or 'IGHV1-NL1' in gene or 'IGHV1-c' in gene or 'IGHV1-f' in gene:
            #     try:
            #         writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=self.args.v_right_length)
            #         writer.write()
            #     except AssertionError:
            #         print '  giving up on',gene
            #         continue
            # else:
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=self.args.v_right_length)
            writer.write()

    # ----------------------------------------------------------------------------------------
    def check_v_right(self, sw_v_right, query_name):
        """ make sure the v_right_length from smith waterman isn't too different to the one we assume when writing the hmm files """
        if abs(sw_v_right - self.args.v_right_length) >= self.safety_buffer_that_I_should_not_need:
            print 'WARNING %s sw v right: %d  (expecting about %d)' % (query_name, sw_v_right, self.args.v_right_length)

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, parameter_dir):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> and print a warning """
        for gene in gene_list:
            hmmfname = parameter_dir + '/hmms/' + utils.sanitize_name(gene) + '.yaml'
            if not os.path.exists(hmmfname):
                print 'WARNING removing match %s from gene list (d.n.e.: %s)' % (gene, hmmfname)
                gene_list.remove(gene)

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, csv_fname, sw_info, parameter_dir, preclusters=None, stripped=False):  # TODO use different input files for the two hmm steps
        with opener('w')(csv_fname) as csvfile:
            # write header
            header = ['name', 'second_name', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seq', 'second_seq']  # I wish I had a good c++ csv reader 
            csvfile.write(' '.join(header) + '\n')

            skipped_gene_matches = set()
            # then write a line for each query sequence (or pair of them)
            if self.args.pair:
                for query_name, second_query_name in self.get_pairs(preclusters):  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
                    info = sw_info[query_name]
                    second_info = sw_info[second_query_name]
    
                    self.check_v_right(info['v_right_length'], query_name)
    
                    # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                    assert len(re.findall('J[123]P', info['j_gene'])) == 0
    
                    # TODO the current method of expanding the k space bounds to encompass all the gene matches, and both sequences, is *correct*, but I think it's unnecessarily slow
                    k_v, k_d = {}, {}
                    k_v['min'] = min(info['k_v']['min'], second_info['k_v']['min'])
                    k_v['max'] = max(info['k_v']['max'], second_info['k_v']['max'])
                    k_d['min'] = min(info['k_d']['min'], second_info['k_d']['min'])
                    k_d['max'] = max(info['k_d']['max'], second_info['k_d']['max'])

                    only_genes = info['all']
                    second_only_genes = second_info['all']
                    if stripped:  # strip down the hmm -- only use the single best gene for each sequence, and don't fuzz at all
                        only_genes = info['v_gene'] + ':' + info['d_gene'] + ':' + info['j_gene']
                        second_only_genes = second_info['v_gene'] + ':' + second_info['d_gene'] + ':' + second_info['j_gene']
                        k_v['min'] = info['k_v']['best']
                        k_v['max'] = k_v['min'] + 1
                        k_d['min'] = info['k_d']['best']
                        k_d['max'] = k_d['min'] + 1

                    only_genes = []
                    tmp_set = set(info['all'].split(':')) | set(second_info['all'].split(':'))  # NOTE using both sets of genes (from both query seqs) like this *really* helps,
                    for gene in tmp_set:  # we only write hmm model files for genes that showed up at least once as a best match, skip genes for which this didn't happen
                        if gene not in sw_info['all_best_matches']:
                            skipped_gene_matches.add(gene)
                            continue
                        only_genes.append(gene)
                    self.check_hmm_existence(only_genes, parameter_dir)
                    csvfile.write('%s %s %d %d %d %d %s %s %s\n' %
                                  (query_name, second_query_name, k_v['min'], k_v['max'], k_d['min'], k_d['max'], ':'.join(only_genes),
                                   self.input_info[query_name]['seq'], self.input_info[second_query_name]['seq']))
            else:
                for query_name in self.input_info:
                    if query_name not in sw_info:
                        print '    %s not found in sw info' % query_name
                        continue
                    info = sw_info[query_name]
    
                    self.check_v_right(info['v_right_length'], query_name)
    
                    # # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                    # assert len(re.findall('J[123]P', info['j_gene'])) == 0
    
                    only_genes = []
                    tmp_set = set(info['all'].split(':'))
                    for gene in tmp_set:  # we only write hmm model files for genes that showed up at least once as a best match, so skip genes for which this didn't happen
                        if self.args.cache_parameters and gene not in sw_info['all_best_matches']:  # except if we're *not* caching parameters, then we hopefully *do* have the gene's hmm written... *sigh* this is getting complicated
                            skipped_gene_matches.add(gene)
                            continue
                        only_genes.append(gene)
                    self.check_hmm_existence(only_genes, parameter_dir)
                    # assert self.args.algorithm == 'viterbi'  # TODO hm, was that really all I had to do to allow non-pair forward?
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, info['k_v']['min'], info['k_v']['max'], info['k_d']['min'], info['k_d']['max'], ':'.join(only_genes), self.input_info[query_name]['seq']))
            if len(skipped_gene_matches) > 0:
                print '    not in best sw matches, skipping: %s' % ','.join([utils.color_gene(gene) for gene in skipped_gene_matches])

    # ----------------------------------------------------------------------------------------
    def run_hmm_binary(self, algorithm, csv_infname, csv_outfname, parameter_dir):
        start = time.time()

        # build the command line
        cmd_str = './ham/ham'
        cmd_str += ' --algorithm ' + algorithm
        if self.args.pair:
            cmd_str += ' --pair '
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + parameter_dir + '/hmms'
        cmd_str += ' --datadir ' + self.args.datadir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        # cmd_str += ' >/dev/null'
        # print cmd_str

        if self.args.debug == 2:  # not sure why, but popen thing hangs with debug 2
            check_call(cmd_str, shell=True)  # um, not sure which I want here, but it doesn't really matter. TODO kinda
            sys.exit()

        # run hmm
        check_call(cmd_str, shell=True)
        # hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        # hmm_proc.wait()
        # hmm_out, hmm_err = hmm_proc.communicate()
        # print hmm_out,
        # print hmm_err,
        # if hmm_proc.returncode != 0:
        #     print 'aaarrrrrrrgggggh\n', hmm_out, hmm_err
        #     sys.exit()

        print '    hmm run time: %.3f' % (time.time() - start)

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, algorithm, hmm_csv_outfname, pairscorefname, pcounter, perfplotter):
        # TODO the input and output files for this function are almost identical at this point

        # write header for pairwise score file
        with opener('w')(pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id,second_unique_id,score\n')
        with opener('r')(hmm_csv_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            last_id = None
            n_boundary_errors = 0
            for line in reader:
                # check for errors
                boundary_error = False
                if line['errors'] != None and 'boundary' in line['errors'].split(':'):
                    n_boundary_errors += 1
                    boundary_error = True
                else:
                    assert len(line['errors']) == 0

                if algorithm == 'viterbi':
                    this_id = utils.get_key(line['unique_id'], line['second_unique_id'])
                    if last_id != this_id:  # if this is the first line (match) for this query (or query pair), print the true event
                        print '%-20s %20s   %d' % (line['unique_id'], line['second_unique_id'], from_same_event(self.args.is_data, self.args.pair, self.reco_info, line['unique_id'], line['second_unique_id']))
                        if pcounter != None:  # increment counters
                            utils.add_match_seqs(self.germline_seqs, line, self.cyst_positions, self.tryp_positions)
                            pcounter.increment(line)
                        if perfplotter != None:
                            perfplotter.evaluate(self.reco_info[line['unique_id']], line, line['unique_id'])

                            # cpos = -1
                            #     cpos = self.cyst_positions[self.reco_info[line['unique_id']]['v_gene']]['cysteine-position']
                            #     # TODO fix this, i.e. figure out the number to add to it so it's correct: tpos = int(self.tryp_positions[self.reco_info[line['unique_id']]['j_gene']]) +
                    if self.args.debug:
                        self.print_hmm_output(line, print_true=(last_id != this_id))
                else:  # for forward, write the pair scores to file to be read by the clusterer
                    with opener('a')(pairscorefname) as pairscorefile:
                        pairscorefile.write('%s,%s,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))

                last_id = utils.get_key(line['unique_id'], line['second_unique_id'])
        if n_boundary_errors > 0:
            print '    %d boundary errors' % n_boundary_errors
        if perfplotter != None:
            perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def print_hmm_output(self, line, print_true=False):
        if print_true and not self.args.is_data:
            print '    true:'
            utils.print_reco_event(self.germline_seqs, self.reco_info[line['unique_id']], extra_str='    ')
            if self.args.pair:
                utils.print_reco_event(self.germline_seqs, self.reco_info[line['second_unique_id']], 0, 0, from_same_event(self.args.is_data, self.args.pair, self.reco_info, line['unique_id'], line['second_unique_id']), '    ')
            print '    inferred:'

        utils.print_reco_event(self.germline_seqs, line, 0, 0, extra_str='    ')
        if self.args.pair:
            tmpseq = line['seq']  # temporarily set 'seq' to the second query's seq. TODO oh, man, that's a cludge
            line['seq'] = line['second_seq']
            utils.print_reco_event(self.germline_seqs, line, 0, 0, True, extra_str='    ')
            line['seq'] = tmpseq

