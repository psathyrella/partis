import time
import sys
import json
from multiprocessing import Process, active_children
import math
import os
import glob
import csv
import re
from subprocess import Popen, check_call, PIPE

import_start = time.time()
import utils
from opener import opener
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
        assert self.args.seqfile != None
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
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=sw_plotdir)
        waterer.run()

        print 'writing hmms with sw info'
        self.write_hmms(sw_parameter_dir, waterer.info['all_best_matches'])

        hmm_parameter_dir = parameter_dir + '/hmm_parameters'
        assert not self.args.pair  # er, could probably do it for forward pair, but I'm not really sure what it would mean
        self.run_hmm('viterbi', waterer.info, parameter_in_dir=sw_parameter_dir, parameter_out_dir=hmm_parameter_dir, count_parameters=True, plotdir=hmm_plotdir)

        print 'writing hmms with hmm info'
        self.write_hmms(hmm_parameter_dir, waterer.info['all_best_matches'])

        # self.clean(waterer)  # TODO get this working again. *man* it's a total bitch keeping track of all the files I'm writing
        os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        assert os.path.exists(self.args.parameter_dir)

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        # cdr3 length partitioning
        cdr3_cluster = False  # don't precluster on cdr3 length for the moment -- I cannot accurately infer cdr3 length in some sequences, so I need a way to pass query seqs to the clusterer with several possible cdr3 lengths. TODO fix that
        cdr3_length_clusters = None
        if cdr3_cluster:
            cdr3_length_clusters = self.cdr3_length_precluster(waterer)

        assert self.args.pair
        hamming_clusters = self.hamming_precluster(cdr3_length_clusters)
        # stripped_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hamming_clusters, stripped=True)
        final_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hamming_clusters, stripped=False)

        # self.clean(waterer)
        os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def point_estimate(self):
        assert os.path.exists(self.args.parameter_dir)
        """ get a point estimate for each query sequence (or pair). i.e. run viterbi """
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.run_hmm('viterbi', waterer.info, parameter_in_dir=self.args.parameter_dir, plot_performance=self.args.plot_performance)
        # self.clean(waterer)  # TODO get this working again. *man* it's a total bitch keeping track of all the files I'm writing

        if self.args.plot_performance:
            plotting.compare_directories('/var/www/sharing/dralph/partis/performance/plots/', 'hmm', '/var/www/sharing/dralph/partis/sw_performance/plots/', 'smith-water', xtitle='inferred - true', stats='rms')
        os.rmdir(self.args.workdir)

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
        print '\nhmm %s' % prefix
        csv_infname = self.args.workdir + '/' + prefix + '_hmm_input.csv'
        csv_outfname = self.args.workdir + '/' + prefix + '_hmm_output.csv'
        pairscorefname = self.args.workdir + '/' + prefix + '_hmm_pairscores.csv'
        self.write_hmm_input(csv_infname, sw_info, preclusters=preclusters, stripped=stripped, parameter_dir=parameter_in_dir)
        if self.args.n_procs > 1:
            self.split_hmm_input(csv_infname)
            for iproc in range(self.args.n_procs):
                p = Process(target=self.run_hmm_binary, args=(algorithm, csv_infname, csv_outfname), kwargs={'parameter_dir':parameter_in_dir, 'iproc':iproc})
                p.start()
            while len(active_children()) > 0:
                print ' wait %s' % len(active_children()),
                sys.stdout.flush()
                time.sleep(1)
            self.merge_hmm_outputs(csv_outfname)
        else:
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
            if os.path.exists(csv_infname):  # if only one proc, this will already be deleted
                os.remove(csv_infname)
            os.remove(csv_outfname)
            os.remove(pairscorefname)

        return clusters

    # ----------------------------------------------------------------------------------------
    def split_hmm_input(self, infname):
        """ split the hmm input in <infname> into <self.args.n_procs> input files in subdirectories within <self.args.workdir> """
        info = []
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                info.append(line)
        queries_per_proc = float(len(info)) / self.args.n_procs
        n_queries_per_proc = int(math.ceil(queries_per_proc))
        for iproc in range(self.args.n_procs):
            workdir = self.args.workdir + '/hmm-' + str(iproc)
            utils.prep_dir(workdir)
            with opener('w')(workdir + '/' + os.path.basename(infname)) as sub_outfile:
                writer = csv.DictWriter(sub_outfile, reader.fieldnames)
                writer.writeheader()
                for iquery in range(iproc*n_queries_per_proc, (iproc + 1)*n_queries_per_proc):
                    if iquery >= len(info):
                        break
                    writer.writerow(info[iquery])

    # ----------------------------------------------------------------------------------------
    def merge_hmm_outputs(self, outfname):
        header = None
        outfo = []
        for iproc in range(self.args.n_procs):
            workdir = self.args.workdir + '/hmm-' + str(iproc)
            with opener('r')(workdir + '/' + os.path.basename(outfname)) as sub_outfile:
                reader = csv.DictReader(sub_outfile)
                header = reader.fieldnames
                for line in reader:
                    outfo.append(line)
            if not self.args.no_clean:
                os.remove(workdir + '/' + os.path.basename(outfname))
                os.rmdir(workdir)

        with opener('w')(outfname) as outfile:
            writer = csv.DictWriter(outfile, header)
            writer.writeheader()
            for line in outfo:
                writer.writerow(line)

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
        from hmmwriter import HmmWriter
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
            print '   ', utils.color_gene(gene),
            sys.stdout.flush()
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene])
            writer.write()

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, skipped_gene_matches, parameter_dir, query_name, second_query_name=None):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> and print a warning """
        # first get the list of genes for which we don't have hmm files
        genes_to_remove = []
        for gene in gene_list:
            hmmfname = parameter_dir + '/hmms/' + utils.sanitize_name(gene) + '.yaml'
            if not os.path.exists(hmmfname):
                if self.args.debug:
                    print '    WARNING %s removed from match list for %s %s (not in %s)' % (utils.color_gene(gene), query_name, '' if second_query_name==None else second_query_name, os.path.dirname(hmmfname))
                skipped_gene_matches.add(gene)
                genes_to_remove.append(gene)

        # then remove 'em from <gene_list>
        for gene in genes_to_remove:
            gene_list.remove(gene)

        # and finally, make sure we're left with at least one gene in each region
        for region in utils.regions:
            if 'IGH' + region.upper() not in ':'.join(gene_list):
                print 'ERROR no %s genes in %s for %s %s' % (region, ':'.join(gene_list), query_name, '' if second_query_name==None else second_query_name)
                print '    skipped %s' % (':'.join(skipped_gene_matches))
                assert False

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
    
                    # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                    assert len(re.findall('J[123]P', info['j_gene'])) == 0
    
                    # TODO the current method of expanding the k space bounds to encompass all the gene matches, and both sequences, is *correct*, but I think it's unnecessarily slow
                    k_v, k_d = {}, {}
                    k_v['min'] = min(info['k_v']['min'], second_info['k_v']['min'])
                    k_v['max'] = max(info['k_v']['max'], second_info['k_v']['max'])
                    k_d['min'] = min(info['k_d']['min'], second_info['k_d']['min'])
                    k_d['max'] = max(info['k_d']['max'], second_info['k_d']['max'])

                    only_genes = info['all'].split(':')
                    second_only_genes = second_info['all'].split(':')
                    self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir, query_name)
                    self.check_hmm_existence(second_only_genes, skipped_gene_matches, parameter_dir, second_query_name)
                    if stripped:  # strip down the hmm -- only use the single best gene for each sequence, and don't fuzz at all
                        assert False  # need to check some things here
                        only_genes = [ info[region + '_gene'] for region in utils.regions ]
                        second_only_genes = [ second_info[region + '_gene'] for region in utils.regions ]
                        k_v['min'] = info['k_v']['best']
                        k_v['max'] = k_v['min'] + 1
                        k_d['min'] = info['k_d']['best']
                        k_d['max'] = k_d['min'] + 1

                    final_only_genes = []
                    tmp_set = set(only_genes) | set(second_only_genes)  # NOTE using both sets of genes (from both query seqs) like this *really* helps,
                    # for gene in tmp_set:  # we only write hmm model files for genes that showed up at least once as a best match in at least one query. Skip genes for which this didn't happen
                    #     if gene not in sw_info['all_best_matches']:
                    #         # NOTE this shouldn't happen much any more since now I'm skipping matches for which we don't have hmm files
                    #         if self.args.debug:
                    #             print '     WARNING removing %s (not in best sw matches)' % utils.color_gene(gene)
                    #         skipped_gene_matches.add(gene)
                    #         continue
                    #     final_only_genes.append(gene)
                    # print only_genes, second_only_genes
                    final_only_genes = list(tmp_set)
                    self.check_hmm_existence(final_only_genes, skipped_gene_matches, parameter_dir, query_name, second_query_name)
                    csvfile.write('%s %s %d %d %d %d %s %s %s\n' %  # ha ha, cool! I just realized csv.DictWriter can handle space-separated files TODO switch over
                                  (query_name, second_query_name, k_v['min'], k_v['max'], k_d['min'], k_d['max'], ':'.join(final_only_genes),
                                   self.input_info[query_name]['seq'], self.input_info[second_query_name]['seq']))
            else:
                for query_name in self.input_info:
                    if query_name not in sw_info:
                        if self.args.debug:
                            print '    %s not found in sw info' % query_name
                        continue
                    info = sw_info[query_name]
    
                    only_genes = []
                    tmp_set = set(info['all'].split(':'))
                    # for gene in tmp_set:  # we only write hmm model files for genes that showed up at least once as a best match, so skip genes for which this didn't happen
                    #     if self.args.cache_parameters and gene not in sw_info['all_best_matches']:  # except if we're *not* caching parameters, then we hopefully *do* have the gene's hmm written... *sigh* this is getting complicated
                    #         # NOTE this shouldn't happen much any more since now I'm skipping matches for which we don't have hmm files
                    #         skipped_gene_matches.add(gene)
                    #         continue
                    #     only_genes.append(gene)
                    only_genes = list(tmp_set)
                    self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir, query_name)
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, info['k_v']['min'], info['k_v']['max'], info['k_d']['min'], info['k_d']['max'], ':'.join(only_genes), self.input_info[query_name]['seq']))
            if len(skipped_gene_matches) > 0:
                print '    not found in %s, i.e. were never the best sw match for any query, so removing from consideration for hmm: %s' % (parameter_dir, ','.join([utils.color_gene(gene) for gene in skipped_gene_matches]))

    # ----------------------------------------------------------------------------------------
    def run_hmm_binary(self, algorithm, csv_infname, csv_outfname, parameter_dir, iproc=-1):
        start = time.time()

        # build the command line
        cmd_str = './packages/ham/ham'
        cmd_str += ' --algorithm ' + algorithm
        if self.args.pair:
            cmd_str += ' --pair '
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + parameter_dir + '/hmms'
        cmd_str += ' --datadir ' + self.args.datadir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname

        workdir = self.args.workdir
        if iproc >= 0:
            workdir += '/hmm-' + str(iproc)
            cmd_str = cmd_str.replace(self.args.workdir, workdir)

        check_call(cmd_str, shell=True)

        if not self.args.no_clean:
            os.remove(csv_infname.replace(self.args.workdir, workdir))

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
                        if self.args.debug:
                            print '%-20s %20s' % (line['unique_id'], line['second_unique_id']),
                            if self.args.pair:
                                print '   %d' % from_same_event(self.args.is_data, self.args.pair, self.reco_info, line['unique_id'], line['second_unique_id']),
                            print ''
                        utils.add_match_info(self.germline_seqs, line, self.cyst_positions, self.tryp_positions, self.args.skip_unproductive, debug=(self.args.debug > 0))
                        if not (self.args.skip_unproductive and line['cdr3_length'] == -1):
                            if pcounter != None:  # increment counters (but only for the best [first] match)
                                pcounter.increment(line)
                            if perfplotter != None:
                                perfplotter.evaluate(self.reco_info[line['unique_id']], line, line['unique_id'])
                    if self.args.debug:
                        if last_id == this_id:
                            utils.add_match_info(self.germline_seqs, line, self.cyst_positions, self.tryp_positions, self.args.skip_unproductive, debug=False)
                        if line['cdr3_length'] == -1:
                            print '      ERROR %s failed to add match info' % line['unique_id']
                        self.print_hmm_output(line, print_true=(last_id != this_id))
                    last_id = utils.get_key(line['unique_id'], line['second_unique_id'])
                else:  # for forward, write the pair scores to file to be read by the clusterer
                    with opener('a')(pairscorefname) as pairscorefile:
                        pairscorefile.write('%s,%s,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))

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
                same_event = from_same_event(self.args.is_data, self.args.pair, self.reco_info, line['unique_id'], line['second_unique_id'])
                utils.print_reco_event(self.germline_seqs, self.reco_info[line['second_unique_id']], one_line=same_event, extra_str='    ')
            print '    inferred:'

        utils.print_reco_event(self.germline_seqs, line, extra_str='    ')
        if self.args.pair:
            tmpseq = line['seq']  # temporarily set 'seq' to the second query's seq. TODO oh, man, that's a cludge
            line['seq'] = line['second_seq']
            utils.print_reco_event(self.germline_seqs, line, one_line=True, extra_str='    ')
            line['seq'] = tmpseq

