import time
import sys
import json
import itertools
from Bio import SeqIO
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
def from_same_event(is_data, pair_hmm, reco_info, query_names):
    if is_data:
        return False
    if len(query_names) > 1:
        reco_id = reco_info[query_names[0]]['reco_id']  # the first one's reco id
        for iq in range(1, len(query_names)):  # then loop through the rest of 'em to see if they're all the same
            if reco_id != reco_info[query_names[iq]]['reco_id']:
                return False
        return True
    else:
        return False

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, args):
        self.args = args
        self.germline_seqs = utils.read_germlines(self.args.datadir)  #, add_fp=True)

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

        if self.args.seqfile != None:
            self.read_input_file()  # read simulation info and write sw input file

        self.outfile = None
        if self.args.outfname != None:
            if os.path.exists(self.args.outfname):
                os.remove(self.args.outfname)
            self.outfile = open(self.args.outfname, 'a')

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
        if '.csv' in self.args.seqfile:
            delimiter = ','
            name_column = 'unique_id'
            seq_column = 'seq'
            seqfile = opener('r')(self.args.seqfile)
            reader = csv.DictReader(seqfile, delimiter=delimiter)
        elif '.tsv' in self.args.seqfile:
            delimiter = '\t'
            name_column = 'name'
            seq_column = 'nucleotide'
            seqfile = opener('r')(self.args.seqfile)
            reader = csv.DictReader(seqfile, delimiter=delimiter)
        elif '.fasta' in self.args.seqfile:
            name_column = 'unique_id'
            seq_column = 'seq'
            reader = []
            n_fasta_queries = 0
            for seq_record in SeqIO.parse(self.args.seqfile, 'fasta'):
                reader.append({})
                reader[-1][name_column] = seq_record.name
                reader[-1][seq_column] = str(seq_record.seq).upper()
                n_fasta_queries += 1
                if self.args.n_max_queries > 0 and n_fasta_queries >= self.args.n_max_queries:
                    break
        n_queries = 0
        for line in reader:
            utils.intify(line)
            # if command line specified query or reco ids, skip other ones
            if self.args.queries != None and line[name_column] not in self.args.queries:
                continue
            if self.args.reco_ids != None and line['reco_id'] not in self.args.reco_ids:
                continue

            self.input_info[line[name_column]] = {'unique_id':line[name_column], 'seq':line[seq_column]}
            if not self.args.is_data:
                if 'fv_insertion' not in line:  # NOTE should be able to remove these lines now
                    line['fv_insertion'] = ''
                if 'jf_insertion' not in line:
                    line['jf_insertion'] = ''
                self.reco_info[line['unique_id']] = line
                utils.add_match_info(self.germline_seqs, line, self.cyst_positions, self.tryp_positions)
            n_queries += 1
            if self.args.n_max_queries > 0 and n_queries >= self.args.n_max_queries:
                break
        assert len(self.input_info) > 0
    
    # ----------------------------------------------------------------------------------------
    def cache_parameters(self, parameter_dir=''):
        assert self.args.plotdir != None
        sw_plotdir, hmm_plotdir = '', ''
        if os.getenv('www') != None:
            sw_plotdir = self.args.plotdir + '/sw_parameters'
            hmm_plotdir = self.args.plotdir + '/hmm_parameters'

        if parameter_dir == '':
            parameter_dir = self.args.parameter_dir

        sw_parameter_dir = parameter_dir + '/sw_parameters'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=sw_plotdir)
        waterer.run()

        print 'writing hmms with sw info'
        self.write_hmms(sw_parameter_dir, waterer.info['all_best_matches'])

        hmm_parameter_dir = parameter_dir + '/hmm_parameters'
        assert self.args.n_sets == 1  # er, could probably do it for forward pair, but I'm not really sure what it would mean
        self.run_hmm('viterbi', waterer.info, parameter_in_dir=sw_parameter_dir, parameter_out_dir=hmm_parameter_dir, hmm_type='k=1', count_parameters=True, plotdir=hmm_plotdir)

        print 'writing hmms with hmm info'
        self.write_hmms(hmm_parameter_dir, waterer.info['all_best_matches'])

        if not self.args.no_clean:
            os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        assert os.path.exists(self.args.parameter_dir)

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        # cdr3 length partitioning
        cdr3_cluster = False  # don't precluster on cdr3 length for the moment -- I cannot accurately infer cdr3 length in some sequences, so I need a way to pass query seqs to the clusterer with several possible cdr3 lengths (and I don't know how to do that!)
        cdr3_length_clusters = None
        if cdr3_cluster:
            cdr3_length_clusters = self.cdr3_length_precluster(waterer)

        hamming_clusters = self.hamming_precluster(cdr3_length_clusters)
        # stripped_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hamming_clusters, stripped=True)
        hmm_clusters = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hamming_clusters, hmm_type='k=2', make_clusters=True)

        self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=hmm_clusters, hmm_type='k=preclusters', prefix='k-')

        # self.clean(waterer)
        if not self.args.no_clean:
            os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self):
        assert os.path.exists(self.args.parameter_dir)
        assert self.args.run_algorithm != None
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.run_hmm(self.args.run_algorithm, waterer.info, parameter_in_dir=self.args.parameter_dir, hmm_type='k=nsets', plot_performance=self.args.plot_performance)
        # self.clean(waterer)
        if not self.args.no_clean:
            os.rmdir(self.args.workdir)

        # if self.args.plot_performance:
        #     plotting.compare_directories(self.args.plotdir + '/hmm-vs-sw', self.args.plotdir + '/hmm/plots', 'hmm', self.args.plotdir + '/sw/plots', 'smith-water', xtitle='inferred - true', stats='rms')

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, sw_info, parameter_in_dir, parameter_out_dir='', preclusters=None, hmm_type='', stripped=False, prefix='', count_parameters=False, plotdir='', plot_performance=False, make_clusters=False):
        pcounter = None
        if count_parameters:
            pcounter = ParameterCounter(self.germline_seqs, parameter_out_dir, plotdir=plotdir)
        true_pcounter = None
        if count_parameters and not self.args.is_data:
            true_pcounter = ParameterCounter(self.germline_seqs, parameter_out_dir, plotdir=plotdir + '/true')
        if plotdir != '':
            utils.prep_dir(plotdir + '/plots', '*.svg')
            if count_parameters and not self.args.is_data:
                utils.prep_dir(plotdir + '/true/plots', '*.svg')

        perfplotter = None
        if plot_performance:
            assert self.args.plotdir != None
            assert not self.args.is_data
            assert self.args.n_sets == 1  # well, you *could* check the performance on k-sets, but I haven't implemented it yet
            assert algorithm == 'viterbi'  # same deal as previous assert
            from performanceplotter import PerformancePlotter
            perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir + '/hmm_performance', 'hmm')

        if prefix == '' and stripped:
            prefix = 'stripped'
        print '\n%shmm' % prefix
        csv_infname = self.args.workdir + '/' + prefix + '_hmm_input.csv'
        csv_outfname = self.args.workdir + '/' + prefix + '_hmm_output.csv'
        pairscorefname = self.args.workdir + '/' + prefix + '_hmm_pairscores.csv'
        self.write_hmm_input(csv_infname, sw_info, preclusters=preclusters, hmm_type=hmm_type, stripped=stripped, parameter_dir=parameter_in_dir)
        if self.args.n_procs > 1:
            self.split_hmm_input(csv_infname)
            for iproc in range(self.args.n_procs):
                proc = Process(target=self.run_hmm_binary, args=(algorithm, csv_infname, csv_outfname), kwargs={'parameter_dir':parameter_in_dir, 'iproc':iproc})
                proc.start()
            while len(active_children()) > 0:
                # print ' wait %s' % len(active_children()),
                sys.stdout.flush()
                time.sleep(1)
            self.merge_hmm_outputs(csv_outfname)
        else:
            self.run_hmm_binary(algorithm, csv_infname, csv_outfname, parameter_dir=parameter_in_dir)
        # if self.outfile != None and algorithm == 'forward':
        #     self.outfile.write('hmm pairscores\n')
        self.read_hmm_output(algorithm, csv_outfname, pairscorefname, pcounter, perfplotter, true_pcounter)

        if count_parameters:
            pcounter.write_counts()
        if count_parameters and not self.args.is_data:
            true_pcounter.write_counts()

        clusters = None
        if make_clusters:  # a.t.m. don't cluster on the k-hmm output -- just use it as a check if we should split any of the clusters
            clusters = Clusterer(self.args.pair_hmm_cluster_cutoff, greater_than=True)
            if self.outfile != None:
                self.outfile.write('hmm clusters\n')
            else:
                print'hmm clusters'
            clusters.cluster(pairscorefname, debug=self.args.debug, reco_info=self.reco_info, outfile=self.outfile, plotdir=self.args.plotdir+'/pairscores')
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
    def run_hmm_binary(self, algorithm, csv_infname, csv_outfname, parameter_dir, iproc=-1):
        # start = time.time()

        # build the command line
        cmd_str = os.getenv('PWD') + '/packages/ham/bcrham'
        if self.args.slurm:
            cmd_str = 'srun ' + cmd_str
        cmd_str += ' --algorithm ' + algorithm
        cmd_str += ' --chunk-cache '
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

        print cmd_str
        check_call(cmd_str, shell=True)

        if not self.args.no_clean:
            os.remove(csv_infname.replace(self.args.workdir, workdir))

        # print '    hmm run time: %.3f' % (time.time() - start)

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
        pair_keys = []  # keep track of the pairs we've done already
        pairs = []
        n_lines, n_preclustered, n_previously_preclustered, n_removable, n_singletons = 0, 0, 0, 0, 0
        for a_name in self.input_info:
            for b_name in self.input_info:
                if b_name == a_name:
                    continue
                key = utils.get_key((a_name, b_name))
                if key in pair_keys:  # already did this pair
                    continue
                pair_keys.append(key)
                if preclusters != None:  # if we've already run preclustering, skip the pairs that we know aren't matches
                    if a_name not in preclusters.query_clusters or b_name not in preclusters.query_clusters:  # singletons (i.e. they were already preclustered into their own group)
                        n_singletons += 1
                        continue
                    if key not in preclusters.pairscores:  # preclustered out in a previous preclustering step
                        n_previously_preclustered += 1
                        continue
                    if preclusters.query_clusters[a_name] != preclusters.query_clusters[b_name]:  # not in same cluster
                        n_preclustered += 1
                        continue
                    if preclusters.is_removable(preclusters.pairscores[key]):  # in same cluster, but score (link) is long. i.e. *this* pair is far apart, but other seqs to which they are linked are close to each other
                        n_removable += 1
                        continue
                pairs.append((a_name, b_name))
                n_lines += 1

        print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)
        return pairs

    # ----------------------------------------------------------------------------------------
    def hamming_precluster(self, preclusters=None):
        start = time.time()
        print 'hamming clustering'
        hammingfname = self.args.workdir + '/fractional-hamming-scores.csv'
        chopped_off_left_sides = False
        with opener('w')(hammingfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'score'))
            writer.writeheader()
            for query_name, second_query_name in self.get_pairs(preclusters):
                query_seq = self.input_info[query_name]['seq']
                second_query_seq = self.input_info[second_query_name]['seq']
                if self.args.truncate_pairs:  # chop off the left side of the longer one if they're not the same length
                    min_length = min(len(query_seq), len(second_query_seq))
                    query_seq = query_seq[-min_length : ]
                    second_query_seq = second_query_seq[-min_length : ]
                    chopped_off_left_sides = True
                mutation_frac = utils.hamming(query_seq, second_query_seq) / float(len(query_seq))
                writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'score':mutation_frac})
                # utils.color_mutants(query_seq, second_query_seq, ref_label=str(query_name) + ' ' + str(second_query_name) + '   ', print_result=True)
                # if self.args.debug:
                #     print '    %20s %20s %8.2f' % (query_name, second_query_name, mutation_frac)

        clust = Clusterer(self.args.hamming_cluster_cutoff, greater_than=False)  # NOTE this 0.5 is reasonable but totally arbitrary
        if self.outfile != None:
            self.outfile.write('hamming clusters\n')
        clust.cluster(hammingfname, debug=self.args.debug, outfile=self.outfile)
        os.remove(hammingfname)
        if chopped_off_left_sides:
            print 'WARNING encountered unequal-length sequences, so chopped off the left-hand sides of each'
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
            print '  %s' % utils.color_gene(gene)
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety,
                               self.germline_seqs[utils.get_region(gene)][gene],
                               self.args)
            writer.write()

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, skipped_gene_matches, parameter_dir, query_name, second_query_name=None):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> and print a warning """
        # first get the list of genes for which we don't have hmm files
        if len(glob.glob(parameter_dir + '/hmms/*.yaml')) == 0:
            print 'ERROR no yamels in %s' % parameter_dir
            sys.exit()

        genes_to_remove = []
        for gene in gene_list:
            hmmfname = parameter_dir + '/hmms/' + utils.sanitize_name(gene) + '.yaml'
            if not os.path.exists(hmmfname):
                # if self.args.debug:
                #     print '    WARNING %s removed from match list for %s %s (not in %s)' % (utils.color_gene(gene), query_name, '' if second_query_name==None else second_query_name, os.path.dirname(hmmfname))
                skipped_gene_matches.add(gene)
                genes_to_remove.append(gene)

        # then remove 'em from <gene_list>
        for gene in genes_to_remove:
            gene_list.remove(gene)

    # ----------------------------------------------------------------------------------------
    def all_regions_present(self, gene_list, skipped_gene_matches, query_name, second_query_name=None):
        # and finally, make sure we're left with at least one gene in each region
        for region in utils.regions:
            if 'IGH' + region.upper() not in ':'.join(gene_list):
                print '       no %s genes in %s for %s %s' % (region, ':'.join(gene_list), query_name, '' if second_query_name==None else second_query_name)
                print '          skipped %s' % (':'.join(skipped_gene_matches))
                print 'ERROR giving up on query'
                return False

        return True

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, sw_info, query_names, parameter_dir, stripped=False, skipped_gene_matches=None):
        combo = {
            'k_v':{'min':99999, 'max':-1},
            'k_d':{'min':99999, 'max':-1},
            'only_genes':[],
            'seqs':[]
        }
        min_length = -1
        for name in query_names:  # first find the min length, so we know how much we'll have to chop off of each sequence
            if min_length == -1 or len(sw_info[name]['seq']) < min_length:
                min_length = len(sw_info[name]['seq'])
        for name in query_names:
            info = sw_info[name]
            query_seq = self.input_info[name]['seq']
            chop = 0
            if self.args.truncate_pairs:  # chop off the left side of the sequence if it's longer than min_length
                chop = max(0, len(query_seq) - min_length)
                query_seq = query_seq[ : min_length]
            combo['seqs'].append(query_seq)
    
            combo['k_v']['min'] = min(info['k_v']['min'] - chop, combo['k_v']['min'])
            combo['k_v']['max'] = max(info['k_v']['max'] - chop, combo['k_v']['max'])
            combo['k_d']['min'] = min(info['k_d']['min'], combo['k_d']['min'])
            combo['k_d']['max'] = max(info['k_d']['max'], combo['k_d']['max'])
    
            only_genes = info['all'].split(':')
            self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir, name)
            if stripped:  # strip down the hmm -- only use the single best gene for each sequence, and don't fuzz at all
                assert False  # need to check some things here
                # only_genes = [ a_info[region + '_gene'] for region in utils.regions ]
                # b_only_genes = [ b_info[region + '_gene'] for region in utils.regions ]
                # k_v['min'] = a_info['k_v']['best'] - a_chop  # NOTE since this only uses <a_info>, you shouldn't expect it to be very accurate
                # k_v['max'] = k_v['min'] + 1
                # k_d['min'] = a_info['k_d']['best']
                # k_d['max'] = k_d['min'] + 1
    
            combo['only_genes'] = list(set(only_genes) | set(combo['only_genes']))  # NOTE using both sets of genes (from both query seqs) like this *really* helps,

        # self.check_hmm_existence(combo['only_genes'], skipped_gene_matches, parameter_dir, name)  # this should be superfluous now
        if not self.all_regions_present(combo['only_genes'], skipped_gene_matches, query_names):
            return {}

        return combo

    # ----------------------------------------------------------------------------------------
    def remove_sw_failures(self, query_names, sw_info):
        non_failed_names = []
        for name in query_names:
            if name in sw_info:
                non_failed_names.append(name)
            else:
                print '    %s not found in sw info' % query_name

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, csv_fname, sw_info, parameter_dir, preclusters=None, hmm_type='', pair_hmm=False, stripped=False):
        csvfile = opener('w')(csv_fname)

        # write header
        header = ['names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seqs']  # I wish I had a good c++ csv reader 
        csvfile.write(' '.join(header) + '\n')

        skipped_gene_matches = set()
        assert hmm_type != ''
        if hmm_type == 'k=1':  # single vanilla hmm
            nsets = [[qn] for qn in self.input_info.keys()]
        elif hmm_type == 'k=2':  # pair hmm
            nsets = self.get_pairs(preclusters)
        elif hmm_type == 'k=preclusters':  # run the k-hmm on each cluster in <preclusters>
            assert preclusters != None
            nsets = [ val for key, val in preclusters.id_clusters.items() ]  # <nsets> is a list of sets (well, lists) of query names
        elif hmm_type == 'k=nsets':  # run on *every* combination of queries which has length <self.args.n_sets>
            nsets = itertools.combinations(self.input_info.keys(), self.args.n_sets)  # I <3 python
        else:
            assert False

        for query_names in nsets:
            self.remove_sw_failures(query_names, sw_info)
            combined_query = self.combine_queries(sw_info, query_names, parameter_dir, stripped=stripped, skipped_gene_matches=skipped_gene_matches)
            if len(combined_query) == 0:  # didn't find all regions
                continue
            csvfile.write('%s %d %d %d %d %s %s\n' %  # NOTE csv.DictWriter can handle tsvs, so this should really be switched to use that
                          (':'.join([str(qn) for qn in query_names]),
                           combined_query['k_v']['min'], combined_query['k_v']['max'],
                           combined_query['k_d']['min'], combined_query['k_d']['max'],
                           ':'.join(combined_query['only_genes']),
                           ':'.join(combined_query['seqs'])))

        if len(skipped_gene_matches) > 0:
            print '    not found in %s, i.e. were never the best sw match for any query, so removing from consideration for hmm:' % (parameter_dir)
            for region in utils.regions:
                print '      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in skipped_gene_matches if utils.get_region(gene) == region]))

        csvfile.close()

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, algorithm, hmm_csv_outfname, pairscorefname, pcounter, perfplotter, true_pcounter):
        # NOTE the input and output files for this function are almost identical at this point

        n_processed = 0
        # write header for pairwise score file
        with opener('w')(pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id,second_unique_id,score\n')
        with opener('r')(hmm_csv_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            last_key = None
            n_boundary_errors = 0
            for line in reader:
                utils.intify(line, splitargs=('unique_ids', 'seqs'))

                # check for errors
                boundary_error = False
                if line['errors'] != None and 'boundary' in line['errors'].split(':'):
                    n_boundary_errors += 1
                    boundary_error = True
                else:
                    assert len(line['errors']) == 0

                ids = line['unique_ids']
                dbg_str = '%s   %d' % (''.join(['%20s ' % i for i in ids]), from_same_event(self.args.is_data, True, self.reco_info, ids))
                if algorithm == 'viterbi':
                    line['seq'] = line['seqs'][0]  # add info for the best match as 'seq'
                    line['unique_id'] = ids[0]
                    utils.add_match_info(self.germline_seqs, line, self.cyst_positions, self.tryp_positions, debug=(self.args.debug > 0))

                    this_key = utils.get_key(ids)
                    if last_key != this_key:  # if this is the first line (i.e. the best viterbi path) for this query (or query pair), print the true event
                        n_processed += 1
                        if self.args.debug:
                            print dbg_str
                        if not (self.args.skip_unproductive and line['cdr3_length'] == -1):
                            if pcounter != None:  # increment counters (but only for the best [first] match)
                                pcounter.increment(line)
                            if true_pcounter != None:  # increment true counters
                                true_pcounter.increment(self.reco_info[ids[0]])
                            if perfplotter != None:
                                perfplotter.evaluate(self.reco_info[ids[0]], line)
                    if self.args.debug:
                        self.print_hmm_output(line, print_true=(last_key != this_key), perfplotter=perfplotter)
                    last_key = utils.get_key(ids)
                else:  # for forward, write the pair scores to file to be read by the clusterer
                    if self.args.debug and not self.args.partition:
                        print dbg_str + '   %10.3f' % float(line['score'])
                    if line['score'] == '-nan':
                        print '    WARNING encountered -nan, setting to -999999.0'
                        score = -999999.0
                    else:
                        score = float(line['score'])
                    with opener('a')(pairscorefname) as pairscorefile:
                        pairscorefile.write('%d,%d,%f\n' % (line['unique_ids'][0], line['unique_ids'][1], score))
                    # if self.args.outfname != None:
                    #     self.outfile.write('%d,%d,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))
                    n_processed += 1

        print '  processed %d queries' % n_processed
        if n_boundary_errors > 0:
            print '    %d boundary errors' % n_boundary_errors
        if perfplotter != None:
            perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def get_true_clusters(self, ids):
        clusters = {}
        for uid in ids:
            rid = self.reco_info[uid]['reco_id']
            found = False
            for clid in clusters:
                if rid == clid:
                    clusters[clid].append(uid)
                    found = True
                    break
            if not found:
                clusters[rid] = [uid]
        return clusters

    # ----------------------------------------------------------------------------------------
    def print_hmm_output(self, line, print_true=False, perfplotter=None):
        out_str_list = []
        ilabel = ''
        if print_true and not self.args.is_data:  # first print true event (if this is simulation)
            for reco_id, uids in self.get_true_clusters(line['unique_ids']).items():
                for iid in range(len(uids)):
                    out_str_list.append(utils.print_reco_event(self.germline_seqs, self.reco_info[uids[iid]], extra_str='    ', return_string=True, label='true:', one_line=(iid!=0)))
            # if self.args.n_sets > 1:
            #     same_event = from_same_event(self.args.is_data, self.args.n_sets > 1, self.reco_info, line['unique_ids'])
            #     out_str_list.append(utils.print_reco_event(self.germline_seqs, self.reco_info[line['unique_ids'][1]], one_line=same_event, extra_str='    ', return_string=True))
            ilabel = 'inferred:'

        out_str_list.append(utils.print_reco_event(self.germline_seqs, line, extra_str='    ', return_string=True, label=ilabel))
        for iextra in range(1, len(line['unique_ids'])):
            line['seq'] = line['seqs'][iextra]
            out_str_list.append(utils.print_reco_event(self.germline_seqs, line, extra_str='    ', return_string=True, one_line=True))
        if self.args.n_sets > 1:
            # assert False  # need to update this
            # tmpseq = line['seqs'].split(':')[0]  # temporarily set 'seq' to the second query's seq. NOTE oh, man, that's a cludge
            line['seq'] = line['seqs'][1]
            out_str_list.append(utils.print_reco_event(self.germline_seqs, line, one_line=True, extra_str='    ', return_string=True))
        line['seq'] = None

        # if not self.args.is_data:
        #     self.print_performance_info(line, perfplotter=perfplotter)

        if self.args.outfname == None:
            print ''.join(out_str_list),
        else:
            self.outfile.write(''.join(out_str_list))

    # ----------------------------------------------------------------------------------------
    def print_performance_info(self, line, perfplotter=None):
        true_line = self.reco_info[line['unique_id']]
        genes_ok = [ 'ok'  if line[region+'_gene']==true_line[region+'_gene'] else 'no' for region in utils.regions]
        print '         v  d  j   hamming      erosions      insertions'
        print '        %3s%3s%3s' % tuple(genes_ok),
        print '  %3d' % (perfplotter.hamming_distance_to_true_naive(true_line, line, line['unique_id']) if perfplotter != None else -1),
        print '   %4d%4d%4d%4d' % tuple([ int(line[ero+'_del']) - int(true_line[ero+'_del']) for ero in utils.real_erosions]),
        print '   %4d%4d' % tuple([len(line[bound+'_insertion']) - len(true_line[bound+'_insertion']) for bound in utils.boundaries])
