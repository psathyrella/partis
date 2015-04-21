import time
import sys
import json
import itertools
import shutil
import math
import os
import glob
import csv
import random
from collections import OrderedDict
from subprocess import Popen, check_call, PIPE

import utils
from opener import opener
from seqfileopener import get_seqfile_info
from clusterer import Clusterer
from glomerator import Glomerator
from waterer import Waterer
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    """ Class to parse input files, start bcrham jobs, and parse/interpret bcrham output for annotation and partitioning """
    def __init__(self, args):
        self.args = args
        self.germline_seqs = utils.read_germlines(self.args.datadir)
        with opener('r')(self.args.datadir + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with opener('r')(self.args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

        randomize_order = self.args.action == 'partition' and not self.args.force_dont_randomize_input_order
        if self.args.seqfile is not None:
            self.input_info, self.reco_info = get_seqfile_info(self.args.seqfile, self.args.is_data, self.germline_seqs, self.cyst_positions, self.tryp_positions,
                                                               self.args.n_max_queries, self.args.queries, self.args.reco_ids, randomize_order=randomize_order)  # if we're partitioning, we need to randomize input order (at least for simulation)

        self.cached_results = None

        self.sw_info = None

        utils.prep_dir(self.args.workdir)
        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'

    # ----------------------------------------------------------------------------------------
    def __del__(self):
        if not self.args.no_clean:
            try:
                os.rmdir(self.args.workdir)
            except OSError:
                raise Exception('workdir (%s) not empty: %s' % (self.args.workdir, ' '.join(os.listdir(self.args.workdir))))  # hm... you get weird recursive exceptions if you get here. Oh, well, it still works

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        """ Infer full parameter sets and write hmm files for sequences from <self.input_info>, first with Smith-Waterman, then using the SW output as seed for the HMM """
        # assert self.args.n_sets == 1  # er, could do it for n > 1, but I'd want to think through a few things first
        assert self.args.plotdir is not None

        sw_parameter_dir = self.args.parameter_dir + '/sw'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=self.args.plotdir + '/sw')
        waterer.run()
        self.sw_info = waterer.info
        self.write_hmms(sw_parameter_dir)

        parameter_out_dir = self.args.parameter_dir + '/hmm'
        self.run_hmm('viterbi', parameter_in_dir=sw_parameter_dir, parameter_out_dir=parameter_out_dir, hmm_type='k=nsets', count_parameters=True, plotdir=self.args.plotdir + '/hmm')
        self.write_hmms(parameter_out_dir)

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self, algorithm):
        """ Just run <algorithm> (either 'forward' or 'viterbi') on sequences in <self.input_info> and exit. You've got to already have parameters cached in <self.args.parameter_dir> """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir (' + self.args.parameter_dir + ') d.n.e')
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()
        self.sw_info = waterer.info
        self.run_hmm(algorithm, parameter_in_dir=self.args.parameter_dir, hmm_type='k=nsets', \
                     count_parameters=self.args.plot_parameters, plotdir=self.args.plotdir)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir %s d.n.e.' % self.args.parameter_dir)

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()
        self.sw_info = waterer.info

        self.glomclusters = []
        self.list_of_preclusters = []
        n_procs = self.args.n_procs
        n_proc_list = []  # list of the number of procs we used for each run
        # self.partitions = []
        while n_procs > 0:
            if len(self.glomclusters) == 0:
                list_of_preclusters = None
                tmpnprocs = len(self.input_info)
                hmm_type = 'k=1'
            else:
                # list_of_preclusters = self.glomclusters[-1].best_minus_ten_partitions
                list_of_preclusters = self.list_of_preclusters[-1]
                tmpnprocs = len(list_of_preclusters[0])
                hmm_type = 'k=preclusters'
            # print '--> %d clusters with %d procs' % (len(self.input_info) if len(self.glomclusters)==0 else len(self.glomclusters[-1].best_minus_ten_partitions[0]), n_procs)  # write_hmm_input uses the best-minus-ten partition
            print '--> %d clusters with %d procs' % (tmpnprocs, n_procs)  # write_hmm_input uses the best-minus-ten partition
            # hmm_type = 'k=1' if len(self.glomclusters)==0 else 'k=preclusters'
            self.run_hmm('forward', self.args.parameter_dir, preclusters=list_of_preclusters, hmm_type=hmm_type, n_procs=n_procs, shuffle_input_order=True)
            n_proc_list.append(n_procs)
            if n_procs == 1:
                break

            # if we already ran with this number of procs, or if we wouldn't be running with too many clusters per process, then reduce <n_procs> for the next run
            if len(n_proc_list) > 1 and n_proc_list[-1] == n_proc_list[-2] or \
               len(self.glomclusters[-1].best_partitions[0]) / n_procs < self.args.max_clusters_per_proc:  # just use the first path/particle, they should all be pretty similar
                if n_procs > 20:
                    n_procs = n_procs / 2
                elif n_procs > 6:
                    n_procs = int(n_procs / 1.5)
                else:
                    n_procs -= 1

        final_partitions = [[] for _ in range(self.args.smc_particles)]
        for glom in self.glomclusters:
            for ipath in range(self.args.smc_particles):
                final_partitions[ipath] += glom.partitions[ipath]
        if self.args.outfname is not None:
            self.write_partitions(final_partitions, self.args.outfname)

        final_glom = self.glomclusters[-1]
        for ipath in range(self.args.smc_particles):  # print the final partitions
            final_glom.print_partition(final_glom.best_partitions[ipath], final_glom.max_log_probs[ipath], 'final')
        final_glom.print_true_partition()

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, partitions, outfname):
        with opener('w')(outfname) as outfile:
            writer = csv.DictWriter(outfile, ('path_index', 'score', 'adj_mi', 'clusters'))  #'normalized_score'
            writer.writeheader()
            for ipath in range(len(partitions)):
                for part in partitions[ipath]:
                    cluster_str = ''
                    for ic in range(len(part['clusters'])):
                        if ic > 0:
                            cluster_str += ';'
                        cluster_str += ':'.join(part['clusters'][ic])
                    writer.writerow({'path_index' : ipath,
                                     'score' : part['score'],
                                     # 'normalized_score' : part['score'] / self.max_log_probs[ipath],
                                     'adj_mi' : part['adj_mi'],
                                     'clusters' : cluster_str})

    # ----------------------------------------------------------------------------------------
    def get_hmm_cmd_str(self, algorithm, csv_infname, csv_outfname, parameter_dir):
        """ Return the appropriate bcrham command string """
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
        cmd_str += ' --hamming-fraction-cutoff ' + str(self.args.hamming_cluster_cutoff)
        if self.args.smc_particles > 1:
            os.environ['GSL_RNG_TYPE'] = 'ranlux'
            os.environ['GSL_RNG_SEED'] = str(random.randint(0, 99999))
            cmd_str += ' --smc-particles ' + str(self.args.smc_particles)
        if self.args.rescale_emissions:
            cmd_str += ' --rescale-emissions'
        if self.args.action == 'partition':
            cmd_str += ' --partition'
            cmd_str += ' --cachefile ' + self.hmm_cachefname

        # print cmd_str
        # sys.exit()
        return cmd_str

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, parameter_in_dir, parameter_out_dir='', preclusters=None, hmm_type='', prefix='', \
                count_parameters=False, plotdir=None, n_procs=None,  # NOTE the local <n_procs>, which overrides the one inside <self.args>
                shuffle_input_order=False):  # @parameterfetishist
        """ Run bcrham, possibly with many processes, and parse and interpret the output """
        print 'hmm'
        if n_procs is None:
            n_procs = self.args.n_procs

        self.write_hmm_input(preclusters=preclusters, hmm_type=hmm_type, parameter_dir=parameter_in_dir, shuffle_input_order=shuffle_input_order)

        print '    running'
        sys.stdout.flush()
        # start = time.time()
        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir)
        if n_procs == 1:
            check_call(cmd_str.split())
            if self.args.action == 'partition':
                # TODO this is messy, you should combine this with the multicore stuff
                glomerer = Glomerator(self.reco_info)
                glomerer.read_cached_agglomeration([self.hmm_outfname,], self.args.smc_particles, debug=False, clean_up=(not self.args.no_clean))  #, outfname=self.hmm_outfname)
                self.glomclusters.append(glomerer)
                self.list_of_preclusters.append(glomerer.combined_conservative_best_minus_ten_partitions)
        else:
            self.split_input(n_procs, infname=self.hmm_infname, prefix='hmm')
            procs = []
            for iproc in range(n_procs):
                workdir = self.args.workdir + '/hmm-' + str(iproc)
                # print cmd_str.replace(self.args.workdir, workdir)
                # continue
                proc = Popen(cmd_str.replace(self.args.workdir, workdir).split(), stdout=PIPE, stderr=PIPE)
                procs.append(proc)
                time.sleep(0.1)
            for iproc in range(len(procs)):
                out, err = procs[iproc].communicate()
                if out != '':  # or err != '':
                    print '  proc %d' % iproc
                    print out
                    # print err
            self.merge_hmm_outputs(n_procs)

        sys.stdout.flush()
        # print '      hmm run time: %.3f' % (time.time()-start)

        if self.args.action == 'partition':
            self.read_cachefile()
            hmminfo = None
            # glom = Glomerator(self.reco_info)
            # glom.read_cached_agglomeration([self.hmm_outfname,], self.args.smc_particles, clean_up=(not self.args.no_clean), debug=True, outfname=self.args.outfname)
            # hmminfo = glom
        else:
            hmminfo = self.read_annotation_output(algorithm, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, plotdir=plotdir)

        if not self.args.no_clean:
            os.remove(self.hmm_infname)

        if self.args.vollmers_clustering:
            vollmers_clusterer = Clusterer()
            vollmers_clusterer.vollmers_cluster(hmminfo, reco_info=self.reco_info)

        return hmminfo

    # ----------------------------------------------------------------------------------------
    def divvy_up_queries(self, n_procs, info):
        # print 'divvy'
        naive_seqs = {}
        for line in info:
            query = line['names']
            # print ' ', query,
            if self.cached_results is not None and query in self.cached_results:
                # print '   cached'  #, self.cached_results[query]
                naive_seqs[query] = self.cached_results[query]['naive-seq']
            elif query in self.sw_info:
                # print '   sw'  #, self.sw_info[query]
                naive_seqs[query] = utils.get_full_naive_seq(self.germline_seqs, self.sw_info[query])
            else:
                print '   NOOOOOope', query
                assert False
            if naive_seqs[query] == '':
                print 'ack', query
                sys.exit()

        clust = Glomerator()
        divvied_queries = clust.naive_seq_glomerate(naive_seqs, n_clusters=n_procs)
        print 'divvy lengths'
        for dq in divvied_queries:
            print len(dq),
        print ''

        if len(divvied_queries) != n_procs:
            raise Exception('Wrong number of clusters')

        return divvied_queries

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname=None, info=None, prefix='sub'):
        """
        If <infname> is specified split the csv info from it into <n_procs> input files in subdirectories labelled with '<prefix>-' within <self.args.workdir>
        If <info> is specified, instead split the list <info> into pieces and return a list of the resulting lists
        """
        if info is None:
            assert infname is not None  # make sure only *one* of 'em is specified
            info = []
            with opener('r')(infname) as infile:
                reader = csv.DictReader(infile, delimiter=' ')
                for line in reader:
                    info.append(line)
        else:
            assert infname is None
            outlists = []
        queries_per_proc = float(len(info)) / n_procs
        if self.args.action == 'partition':
            divvied_queries = self.divvy_up_queries(n_procs, info)
        for iproc in range(n_procs):
            if infname is None:
                outlists.append([])
            else:
                subworkdir = self.args.workdir + '/' + prefix + '-' + str(iproc)
                utils.prep_dir(subworkdir)
                sub_outfile = opener('w')(subworkdir + '/' + os.path.basename(infname))
                writer = csv.DictWriter(sub_outfile, reader.fieldnames, delimiter=' ')
                writer.writeheader()
            for iquery in range(len(info)):
                if self.args.action == 'partition':
                    if info[iquery]['names'] not in divvied_queries[iproc]:  # NOTE I think the reason this doesn't seem to be speeding things up is that our hierarhical agglomeration time is dominated by the distance calculation, and that distance calculation time is roughly proportional to the number of sequences in the cluster (i.e. larger clusters take longer)
                        continue
                else:
                    if iquery % n_procs != iproc:
                        continue
                if infname is None:
                    outlists[-1].append(info[iquery])
                else:
                    writer.writerow(info[iquery])
            if infname is not None:
                sub_outfile.close()

            if os.path.exists(self.hmm_cachefname):
                check_call(['cp', self.hmm_cachefname, subworkdir + '/'])  # NOTE this is kind of wasteful to write it to each subdirectory (it could be large) but it's cleaner this way, 'cause then the subdirs are independent

        if infname is None:
            return outlists
        # sys.exit()

    # # ----------------------------------------------------------------------------------------
    # def merge_partition_files(self, fname, n_procs):
    #     """ 
    #     Merge the output partition files from several independent bcrham hierarchical agglomeration processes.
    #     The method is to 'rewind' each individual partition back by 10 units of log probability, then merge these rewound partitions.
    #     This is conservative, and doesn't seem to kick up any problems, but in the end is of course dependent on how accurate our model is.
    #     """

    #     # merged_log_probs = [0.0 for _ in range(self.args.smc_particles)]
    #     # merged_partitions = [[] for _ in range(self.args.smc_particles)]
    #     # for iproc in range(n_procs):
    #     #     workdir = self.args.workdir + '/hmm-' + str(iproc)
    #     #     infnames.append(workdir + '/' + os.path.basename(fname))
    #     #     # glomerer = Glomerator()
    #     #     # glomerer.read_cached_agglomeration(infname=workdir + '/' + os.path.basename(fname), partitions=self.partitions, debug=False, clean_up=(not self.args.no_clean))
    #     #     # for ipath in range(self.args.smc_particles):
    #     #     #     if math.isnan(glomerer.max_minus_ten_log_probs[ipath]):  # NOTE this should really have a way of handling -INFINITY
    #     #     #         raise Exception('nan while merging outputs ' + str(glomerer.max_minus_ten_log_probs[ipath]))
    #     #     #     merged_log_probs[ipath] += glomerer.max_minus_ten_log_probs[ipath]
    #     #     #     for cluster in glomerer.best_minus_ten_partitions[ipath]:
    #     #     #         merged_partitions[ipath].append(cluster)
    #     #     # # if not self.args.no_clean:
    #     #     # #     os.remove(workdir + '/' + os.path.basename(fname))
    #     infnames = [self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(fname) for iproc in range(n_procs)]
    #     glomerer = Glomerator(self.reco_info)
    #     glomerer.read_cached_agglomeration(infnames, self.args.smc_particles, debug=False, clean_up=(not self.args.no_clean))

    #     with opener('w')(fname) as partitionfile:
    #         writer = csv.DictWriter(partitionfile, ('path_index', 'partition', 'score'))
    #         writer.writeheader()
    #         for ipath in range(self.args.smc_particles):
    #             writer.writerow({'path_index' : ipath,
    #                              'partition' : ';'.join([':'.join([uid for uid in uids]) for uids in merged_partitions[ipath]]),
    #                              'score' : merged_log_probs[ipath]})

    # ----------------------------------------------------------------------------------------
    def merge_csv_files(self, fname, n_procs):
        """ Merge the output csv files from subsidiary bcrham processes, remaining agnostic about the csv content """
        header = None
        outfo = []
        for iproc in range(n_procs):
            workdir = self.args.workdir + '/hmm-' + str(iproc)
            with opener('r')(workdir + '/' + os.path.basename(fname)) as sub_outfile:
                reader = csv.DictReader(sub_outfile)
                header = reader.fieldnames
                for line in reader:
                    outfo.append(line)
            if not self.args.no_clean:
                os.remove(workdir + '/' + os.path.basename(fname))

        with opener('w')(fname) as outfile:
            writer = csv.DictWriter(outfile, header)
            writer.writeheader()
            for line in outfo:
                writer.writerow(line)

    # ----------------------------------------------------------------------------------------
    def merge_hmm_outputs(self, n_procs):
        """ Merge any/all output files from subsidiary bcrham processes """
        if self.args.action == 'partition':  # merge partitions from several files
            infnames = [self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(self.hmm_outfname) for iproc in range(n_procs)]
            glomerer = Glomerator(self.reco_info)
            glomerer.read_cached_agglomeration(infnames, self.args.smc_particles, debug=False, clean_up=(not self.args.no_clean))  #, outfname=self.hmm_outfname)
            self.glomclusters.append(glomerer)
            self.list_of_preclusters.append(glomerer.combined_conservative_best_minus_ten_partitions)
            # self.merge_partition_files(self.hmm_outfname, n_procs)
            self.merge_csv_files(self.hmm_cachefname, n_procs)
        else:
            self.merge_csv_files(self.hmm_outfname, n_procs)

        if not self.args.no_clean:
            for iproc in range(n_procs):
                workdir = self.args.workdir + '/hmm-' + str(iproc)
                os.remove(workdir + '/' + os.path.basename(self.hmm_infname))
                os.rmdir(workdir)

    # ----------------------------------------------------------------------------------------
    def get_pairs(self, preclusters=None):
        """ Get all unique the pairs of sequences in input_info, skipping where preclustered out """
        all_pairs = itertools.combinations(self.input_info.keys(), 2)
        if preclusters == None:
            print '    ?? lines (no preclustering)'  # % len(list(all_pairs)) NOTE I'm all paranoid the list conversion will be slow (although it doesn't seem to be a.t.m.)
            return all_pairs
        else:  # if we've already run preclustering, skip the pairs that we know aren't matches
            preclustered_pairs = []
            n_lines, n_preclustered, n_previously_preclustered, n_removable, n_singletons = 0, 0, 0, 0, 0
            for a_name, b_name in all_pairs:
                key = utils.get_key((a_name, b_name))
                # NOTE shouldn't need this any more:
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
                preclustered_pairs.append((a_name, b_name))
                n_lines += 1
            print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)
            return preclustered_pairs

    # ----------------------------------------------------------------------------------------
    def write_hmms(self, parameter_dir):
        """ Write hmm model files to <parameter_dir>/hmms, using information from <parameter_dir> """
        print 'writing hmms with info from %s' % parameter_dir
        # start = time.time()
        from hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')

        gene_list = self.args.only_genes
        if gene_list == None and self.sw_info is not None:  # if specific genes weren't specified, do the ones for which we have sw matches
            gene_list = []
            for region in utils.regions:
                for gene in self.germline_seqs[region]:
                    if gene in self.sw_info['all_best_matches']:
                        gene_list.append(gene)

        assert gene_list is not None
        for gene in gene_list:
            if self.args.debug:
                print '  %s' % utils.color_gene(gene)
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety,
                               self.germline_seqs[utils.get_region(gene)][gene],
                               self.args)
            writer.write()

        # print '    time to write hmms: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, skipped_gene_matches, parameter_dir):  #, query_name, second_query_name=None):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> and print a warning """
        # first get the list of genes for which we don't have hmm files
        if len(glob.glob(parameter_dir + '/hmms/*.yaml')) == 0:
            raise Exception('no yamels in %s' % parameter_dir)

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
        """ Check that we have at least one gene for each region """
        for region in utils.regions:
            if 'IGH' + region.upper() not in ':'.join(gene_list):
                print '       no %s genes in %s for %s %s' % (region, ':'.join(gene_list), query_name, '' if (second_query_name == None) else second_query_name)
                print '          skipped %s' % (':'.join(skipped_gene_matches))
                print 'giving up on query'
                return False

        return True

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, query_names, parameter_dir, skipped_gene_matches=None):
        """ Return the 'logical OR' of the queries in <query_names>, i.e. the maximal extent in k_v/k_d space and OR of only_gene sets """
        combo = {
            'k_v':{'min':99999, 'max':-1},
            'k_d':{'min':99999, 'max':-1},
            'only_genes':[],
            'seqs':[],
            'mute-freqs':[]
        }
        min_length = -1
        for name in query_names:  # first find the min length, so we know how much we'll have to chop off of each sequence
            if min_length == -1 or len(self.sw_info[name]['seq']) < min_length:
                min_length = len(self.sw_info[name]['seq'])
        for name in query_names:
            info = self.sw_info[name]
            query_seq = self.input_info[name]['seq']
            chop = 0
            if self.args.truncate_pairs:  # chop off the left side of the sequence if it's longer than min_length
                chop = max(0, len(query_seq) - min_length)
                query_seq = query_seq[ : min_length]
            combo['seqs'].append(query_seq)
            # for region in utils.regions:
            #     print '  ', region, name, utils.get_mutation_rate(self.germline_seqs, self.sw_info[name], restrict_to_region=region)
            combo['mute-freqs'].append(utils.get_mutation_rate(self.germline_seqs, self.sw_info[name]))  # TODO this just always uses the SW mutation rate, but I should really update it with the (multi-)hmm-derived ones (same goes for k space boundaries)

            combo['k_v']['min'] = min(info['k_v']['min'] - chop, combo['k_v']['min'])
            combo['k_v']['max'] = max(info['k_v']['max'] - chop, combo['k_v']['max'])
            combo['k_d']['min'] = min(info['k_d']['min'], combo['k_d']['min'])
            combo['k_d']['max'] = max(info['k_d']['max'], combo['k_d']['max'])

            only_genes = info['all'].split(':')  # sw matches for this query
            self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir)  #, name)
            used_only_genes = []
            for region in utils.regions:
                reg_genes = [g for g in only_genes if utils.get_region(g) == region]
                n_genes = min(len(reg_genes), int(self.args.n_max_per_region[utils.regions.index(region)]))  # minimum of [the number of gene matches for this region] and [the number we want for this region]
                for ig in range(n_genes):
                    used_only_genes.append(reg_genes[ig])

            combo['only_genes'] = list(set(used_only_genes) | set(combo['only_genes']))  # NOTE using the OR of all sets of genes (from all query seqs) like this *really* helps,

        # self.check_hmm_existence(combo['only_genes'], skipped_gene_matches, parameter_dir, name)  # this should be superfluous now
        if not self.all_regions_present(combo['only_genes'], skipped_gene_matches, query_names):
            return {}

        return combo

    # ----------------------------------------------------------------------------------------
    def remove_sw_failures(self, query_names):
        """ If any of the queries in <query_names> was unproductive, return an empty list (which will be skipped entirely), otherwise return the original name list """
        unproductive = False
        for qrn in query_names:
            if qrn in self.sw_info['skipped_unproductive_queries']:
                # print '    skipping unproductive %s along with %s' % (query_names[0], ' '.join(query_names[1:]))
                unproductive = True
        if unproductive:
            return []

        # otherwise they should be in self.sw_info, but doesn't hurt to check
        return_names = []
        for name in query_names:
            if name in self.sw_info:
                return_names.append(name)
            else:
                print '    %s not found in sw info' % ' '.join([qn for qn in query_names])
        return return_names

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, parameter_dir, preclusters=None, hmm_type='', shuffle_input_order=False):
        """ Write input file for bcrham """
        print '    writing input'
        if self.cached_results is not None:
            with opener('w')(self.hmm_cachefname) as cachefile:
                writer = csv.DictWriter(cachefile, ('unique_ids', 'score', 'naive-seq'))
                writer.writeheader()
                for uids, cachefo in self.cached_results.items():
                    writer.writerow({'unique_ids':uids, 'score':cachefo['logprob'], 'naive-seq':cachefo['naive-seq']})

        csvfile = opener('w')(self.hmm_infname)
        header = ['path_index', 'names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seqs', 'mute_freqs']
        writer = csv.DictWriter(csvfile, header, delimiter=' ')  # NOTE should eventually rewrite arg parser in ham to handle csvs (like in glomerator cache reader)
        writer.writeheader()
        # start = time.time()

        skipped_gene_matches = set()
        list_of_nsets = []  # will be length one unless we're partitioning with smc
        if hmm_type == 'k=1':  # single vanilla hmm
            list_of_nsets.append([[qn] for qn in self.input_info.keys()])
        elif hmm_type == 'k=preclusters':  # run the k-hmm on each cluster in <preclusters>
            # list_of_nsets = preclusters.best_minus_ten_partitions
            list_of_nsets = preclusters
        elif hmm_type == 'k=nsets':
            if self.args.all_combinations:  # run on *every* combination of queries which has length <self.args.n_sets>
                list_of_nsets.append(itertools.combinations(self.input_info.keys(), self.args.n_sets))
            else:  # put the first n together, and the second group of n (note that self.input_info is an OrderedDict)
                nsets = []
                keylist = self.input_info.keys()
                this_set = []
                for iquery in range(len(keylist)):
                    if iquery % self.args.n_sets == 0:  # every nth query, start a new group
                        if len(this_set) > 0:
                            nsets.append(this_set)
                        this_set = []
                    this_set.append(keylist[iquery])
                if len(this_set) > 0:
                    nsets.append(this_set)
                list_of_nsets.append(nsets)
        else:
            assert False

        # randomize the order of the query list in <nsets>. Note that the list gets split into chunks for parallelization later
        if shuffle_input_order:  # This can be *really* important if you're partitioning in parallel
            for ins in range(len(list_of_nsets)):  # this will be the trivial loop unless we're partitioning with smc
                random_nsets = []
                while len(list_of_nsets[ins]) > 0:
                    irand = random.randint(0, len(list_of_nsets[ins]) - 1)  # NOTE interval is inclusive
                    random_nsets.append(list_of_nsets[ins][irand])
                    list_of_nsets[ins].remove(list_of_nsets[ins][irand])
                list_of_nsets[ins] = random_nsets

        ins = 0
        for nsets in list_of_nsets:
            for query_names in nsets:
                non_failed_names = self.remove_sw_failures(query_names)
                if len(non_failed_names) == 0:
                    continue
                combined_query = self.combine_queries(non_failed_names, parameter_dir, skipped_gene_matches=skipped_gene_matches)
                if len(combined_query) == 0:  # didn't find all regions
                    continue
                writer.writerow({
                    'path_index' : ins,
                    'names' : ':'.join([qn for qn in non_failed_names]),
                    'k_v_min' : combined_query['k_v']['min'],
                    'k_v_max' : combined_query['k_v']['max'],
                    'k_d_min' : combined_query['k_d']['min'],
                    'k_d_max' : combined_query['k_d']['max'],
                    'only_genes' : ':'.join(combined_query['only_genes']),
                    'seqs' : ':'.join(combined_query['seqs']),
                    'mute_freqs' : ':'.join([str(f) for f in combined_query['mute-freqs']])
                })
            ins += 1

        if len(skipped_gene_matches) > 0:
            print '    not found in %s, i.e. were never the best sw match for any query, so removing from consideration for hmm:' % (parameter_dir)
            for region in utils.regions:
                print '      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in sorted(skipped_gene_matches) if utils.get_region(gene) == region]))

        csvfile.close()
        # print '        input write time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self):
        """ Read cached bcrham partition info """
        if self.cached_results is None:
            self.cached_results = {}

        with opener('r')(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                if line['errors'] != '':
                    raise Exception('in bcrham output for %s: %s ' % (line['unique_ids'], line['errors']))
                if line['unique_ids'] not in self.cached_results:
                    self.cached_results[line['unique_ids']] = {'logprob':float(line['score']), 'naive-seq':line['naive-seq']}
                    if line['naive-seq'] == '':  # I forget why this was happening, but it shouldn't any more (note that I don't actually *remove* the check, though...)
                        raise Exception(line['unique_ids'])

        if not self.args.no_clean:
            os.remove(self.hmm_cachefname)

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, algorithm, count_parameters=False, parameter_out_dir=None, plotdir=None):
        """ Read bcrham annotation output """
        print '    read output'

        if count_parameters:
            assert parameter_out_dir is not None
            assert plotdir is not None
        pcounter = ParameterCounter(self.germline_seqs) if count_parameters else None
        true_pcounter = ParameterCounter(self.germline_seqs) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter(self.germline_seqs, plotdir + '/hmm/performance', 'hmm') if self.args.plot_performance else None

        n_seqs_processed, n_events_processed = 0, 0
        # hmminfo = OrderedDict()  # no longer used, I *think*
        with opener('r')(self.hmm_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            boundary_error_queries = []
            for line in reader:
                utils.process_input_line(line,
                                         splitargs=('unique_ids', 'seqs'),
                                         int_columns=('nth_best', 'v_5p_del', 'd_5p_del', 'cdr3_length', 'j_5p_del', 'j_3p_del', 'd_3p_del', 'v_3p_del'),
                                         float_columns=('score'))
                ids = line['unique_ids']
                same_event = utils.from_same_event(self.args.is_data, self.reco_info, ids)
                if same_event is None:
                    same_event = -1
                id_str = ''.join(['%20s ' % i for i in ids])

                # check for errors
                if line['nth_best'] == 0:  # if this is the first line for this set of ids (i.e. the best viterbi path or only forward score)
                    if line['errors'] is not None and 'boundary' in line['errors'].split(':'):
                        boundary_error_queries.append(':'.join([uid for uid in ids]))
                    else:
                        assert len(line['errors']) == 0

                utils.add_cdr3_info(self.germline_seqs, self.cyst_positions, self.tryp_positions, line)
                if self.args.debug:
                    if line['nth_best'] == 0:  # if this is the first line (i.e. the best viterbi path) for this query (or query pair), print the true event
                        print '%s   %d' % (id_str, same_event)
                    self.print_hmm_output(line, print_true=(line['nth_best']==0))  #, perfplotter=perfplotter)
                if line['nth_best'] == 0 and (line['cdr3_length'] != -1 or not self.args.skip_unproductive):  # if it's productive, or if we're not skipping unproductive rearrangements
                    if pcounter is not None:
                        pcounter.increment_reco_params(line)
                    if true_pcounter is not None:
                        true_pcounter.increment_reco_params(self.reco_info[ids[0]])  # NOTE doesn't matter which id you pass it, since they all have the same reco parameters
                    n_events_processed += 1
                    for iseq in range(len(ids)):
                        tmp_line = dict(line)  # make a copy of the info, into which we'll insert the sequence-specific stuff
                        tmp_line['seq'] = line['seqs'][iseq]
                        tmp_line['unique_id'] = ids[iseq]
                        utils.add_match_info(self.germline_seqs, tmp_line, self.cyst_positions, self.tryp_positions, debug=(self.args.debug > 0))
                        if pcounter is not None:
                            pcounter.increment_mutation_params(tmp_line)
                        if true_pcounter is not None:
                            true_pcounter.increment_mutation_params(self.reco_info[ids[iseq]])  # NOTE doesn't matter which id you pass it, since they all have the same reco parameters
                        if perfplotter is not None:
                            perfplotter.evaluate(self.reco_info[ids[iseq]], tmp_line)
                        n_seqs_processed += 1

        if pcounter is not None:
            pcounter.write(parameter_out_dir)
            if not self.args.no_plot:
                pcounter.plot(plotdir, subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)
        if true_pcounter is not None:
            true_pcounter.write(parameter_out_dir + '/true')
            if not self.args.no_plot:
                true_pcounter.plot(plotdir + '/true', subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)
        if perfplotter is not None:
            perfplotter.plot()

        print '  processed %d sequences (%d events)' % (n_seqs_processed, n_events_processed)
        if len(boundary_error_queries) > 0:
            print '    %d boundary errors (%s)' % (len(boundary_error_queries), ', '.join(boundary_error_queries))

        if self.args.outfname is not None:
            outpath = self.args.outfname
            if self.args.outfname[0] != '/':  # if full output path wasn't specified on the command line
                outpath = os.getcwd() + '/' + outpath
            shutil.copyfile(self.hmm_outfname, outpath)

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)

        # return hmminfo

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
    def print_hmm_output(self, line, print_true=False):  #, perfplotter=None):
        out_str_list = []
        ilabel = ''
        if print_true and not self.args.is_data:  # first print true event (if this is simulation)
            for uids in self.get_true_clusters(line['unique_ids']).values():
                for iid in range(len(uids)):
                    true_event_str = utils.print_reco_event(self.germline_seqs, self.reco_info[uids[iid]], extra_str='    ', return_string=True, label='true:', one_line=(iid != 0))
                    out_str_list.append(true_event_str)
            ilabel = 'inferred:'

        for iseq in range(0, len(line['unique_ids'])):
            tmpline = dict(line)
            tmpline['seq'] = line['seqs'][iseq]
            label = ilabel if iseq==0 else ''
            event_str = utils.print_reco_event(self.germline_seqs, tmpline, extra_str='    ', return_string=True, label=label, one_line=(iseq>0))
            out_str_list.append(event_str)

        # if not self.args.is_data:
        #     self.print_performance_info(line, perfplotter=perfplotter)

        print ''.join(out_str_list),

    # ----------------------------------------------------------------------------------------
    def print_performance_info(self, line, perfplotter=None):
        true_line = self.reco_info[line['unique_id']]
        genes_ok = ['ok'  if (line[region+'_gene'] == true_line[region+'_gene']) else 'no' for region in utils.regions]
        print '         v  d  j   hamming      erosions      insertions'
        print '        %3s%3s%3s' % tuple(genes_ok),
        print '  %3d' % (perfplotter.hamming_distance_to_true_naive(true_line, line, line['unique_id']) if perfplotter != None else -1),
        print '   %4d%4d%4d%4d' % tuple([int(line[ero+'_del']) - int(true_line[ero+'_del']) for ero in utils.real_erosions]),
        print '   %4d%4d' % tuple([len(line[bound+'_insertion']) - len(true_line[bound+'_insertion']) for bound in utils.boundaries])
