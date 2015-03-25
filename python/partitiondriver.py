import time
import sys
import json
import itertools
import math
import os
import glob
import csv
import random
from collections import OrderedDict
from subprocess import Popen, check_call

import utils
from opener import opener
from seqfileopener import get_seqfile_info
from clusterer import Clusterer
from waterer import Waterer
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
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
                raise Exception('ERROR workdir (%s) not empty: %s' % (self.args.workdir, ' '.join(os.listdir(self.args.workdir))))  # hm... you get weird recursive exceptions if you get here. Oh, well, it still works

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        assert self.args.n_sets == 1  # er, could do it for n > 1, but I'd want to think through a few things first
        assert self.args.plotdir is not None

        sw_parameter_dir = self.args.parameter_dir + '/sw'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=self.args.plotdir + '/sw')
        waterer.run()
        self.write_hmms(sw_parameter_dir, waterer.info['all_best_matches'])

        parameter_out_dir = self.args.parameter_dir + '/hmm'
        self.run_hmm('viterbi', waterer.info, parameter_in_dir=sw_parameter_dir, parameter_out_dir=parameter_out_dir, hmm_type='k=1', count_parameters=True, plotdir=self.args.plotdir + '/hmm')
        self.write_hmms(parameter_out_dir, waterer.info['all_best_matches'])

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self, algorithm):
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('ERROR parameter dir (' + self.args.parameter_dir + ') d.n.e')
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.run_hmm(algorithm, waterer.info, parameter_in_dir=self.args.parameter_dir, hmm_type='k=nsets', \
                     count_parameters=self.args.plot_parameters, plotdir=self.args.plotdir)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('ERROR parameter dir %s d.n.e.' % self.args.parameter_dir)

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        partition, glomclusters = None, None
        n_procs = self.args.n_procs
        n_proc_list = []  # list of the number of procs we used for each run
        while n_procs > 0:
            print '--> %d clusters with %d procs' % (len(self.input_info) if glomclusters is None else len(glomclusters.best_minus_ten_partition), n_procs)  # write_hmm_input uses the best-minus-ten partition
            hmm_type = 'k=1' if glomclusters is None else 'k=preclusters'
            partition = self.run_hmm('forward', waterer.info, self.args.parameter_dir, preclusters=glomclusters, hmm_type=hmm_type,
                                     n_procs=n_procs, shuffle_input_order=True)
            n_proc_list.append(n_procs)
            glomclusters = Clusterer()
            glomclusters.hierarch_agglom(log_probs=self.cached_results, partitions=partition, reco_info=self.reco_info, workdir=self.args.workdir, debug=True)
            if n_procs == 1:
                break

            # if we already ran with this number of procs, or if we wouldn't be running with too many clusters per process, then reduce <n_procs> for the next run
            # TODO I think there's really no way around the fact that eventually I'm going to have to delve into the cached log probs to decide how many procs
            if len(n_proc_list) > 1 and n_proc_list[-1] == n_proc_list[-2] or \
               len(glomclusters.best_partition) / n_procs < self.args.max_clusters_per_proc:
                if n_procs > 20:
                    n_procs = n_procs / 2
                elif n_procs > 6:
                    n_procs = int(n_procs / 1.5)
                else:
                    n_procs -= 1

    # ----------------------------------------------------------------------------------------
    def get_hmm_cmd_str(self, algorithm, csv_infname, csv_outfname, parameter_dir):
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
        cmd_str += ' --naive-preclustering'  # TODO remove the option to not do this from ham code
        if self.args.action == 'partition':
            cmd_str += ' --partition'
            cmd_str += ' --cachefile ' + self.hmm_cachefname

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, sw_info, parameter_in_dir, parameter_out_dir='', preclusters=None, hmm_type='', prefix='', \
                count_parameters=False, plotdir=None, n_procs=None,  # NOTE the local <n_procs>, which overrides the one inside <self.args>
                shuffle_input_order=False):  # @parameterfetishist
        print 'hmm'
        if n_procs is None:
            n_procs = self.args.n_procs

        self.write_hmm_input(sw_info, preclusters=preclusters, hmm_type=hmm_type, parameter_dir=parameter_in_dir, shuffle_input_order=shuffle_input_order)

        print '    running'
        sys.stdout.flush()
        # start = time.time()
        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir)
        if n_procs == 1:
            check_call(cmd_str.split())
        else:
            self.split_input(n_procs, infname=self.hmm_infname, prefix='hmm')
            procs = []
            for iproc in range(n_procs):
                workdir = self.args.workdir + '/hmm-' + str(iproc)
                procs.append(Popen(cmd_str.replace(self.args.workdir, workdir).split()))
                time.sleep(0.1)
            for proc in procs:
                proc.wait()
            self.merge_hmm_outputs(n_procs)

        sys.stdout.flush()
        # print '      hmm run time: %.3f' % (time.time()-start)

        if self.args.action == 'partition':
            hmminfo = self.read_partition_output()
        else:
            hmminfo = self.read_annotation_output(algorithm, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, plotdir=plotdir)

        if not self.args.no_clean:
            os.remove(self.hmm_infname)

        if self.args.vollmers_clustering:
            vollmers_clusterer = Clusterer()
            vollmers_clusterer.vollmers_cluster(hmminfo, reco_info=self.reco_info, workdir=self.args.workdir)

        return hmminfo

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname=None, info=None, prefix='sub'):
        """ 
        If <infname> is specified split the csv info from it into <n_procs> input files in subdirectories labelled with '<prefix>-' within <self.args.workdir>
        If <info> is specified, instead split the list <info> into pieces and return a list of the resulting lists
        """
        if info is None:
            assert infname is not None
            info = []
            with opener('r')(infname) as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    info.append(line)
        else:
            assert infname is None  # make sure only *one* of 'em is specified
            outlists = []
        queries_per_proc = float(len(info)) / n_procs
        # n_queries_per_proc = int(math.ceil(queries_per_proc))
        for iproc in range(n_procs):
            if infname is None:
                outlists.append([])
            else:
                subworkdir = self.args.workdir + '/' + prefix + '-' + str(iproc)
                utils.prep_dir(subworkdir)
                sub_outfile = opener('w')(subworkdir + '/' + os.path.basename(infname))
                writer = csv.DictWriter(sub_outfile, reader.fieldnames)
                writer.writeheader()
            # for iquery in range(iproc*n_queries_per_proc, (iproc + 1)*n_queries_per_proc):  # NOTE this is the old version
            for iquery in range(len(info)):  # TODO this is probably pretty wasteful
                # if iquery >= len(info):  # NOTE this is the old version
                #     break
                if iquery % n_procs != iproc:
                    continue
                if infname is None:
                    outlists[-1].append(info[iquery])
                else:
                    writer.writerow(info[iquery])

            if os.path.exists(self.hmm_cachefname):
                check_call(['cp', self.hmm_cachefname, subworkdir + '/'])  # NOTE this is kind of wasteful to write it to each subdirectory (it could be large) but it's cleaner this way, 'cause then the subdirs are independent

        if infname is None:
            return outlists

    # ----------------------------------------------------------------------------------------
    def merge_partition_files(self, fname, n_procs):
        merged_log_prob = 0
        merged_partition = []
        for iproc in range(n_procs):
            workdir = self.args.workdir + '/hmm-' + str(iproc)
            glomerer = Clusterer()
            partition_info = self.read_partition_info(workdir + '/' + os.path.basename(fname))
            glomerer.hierarch_agglom(partitions=partition_info, debug=False)
            if math.isnan(glomerer.max_minus_ten_log_prob):  # NOTE this should really have a way of handling -INFINITY
                raise Exception('ERROR nan while merging outputs ' + str(glomerer.max_minus_ten_log_prob))
            merged_log_prob += glomerer.max_minus_ten_log_prob
            for cluster in glomerer.best_minus_ten_partition:
                merged_partition.append(cluster)

            if not self.args.no_clean:
                os.remove(workdir + '/' + os.path.basename(fname))

        with opener('w')(fname) as partitionfile:
            writer = csv.DictWriter(partitionfile, ('partition', 'score'))
            writer.writeheader()
            writer.writerow({'partition' : ';'.join([':'.join([str(uid) for uid in uids]) for uids in merged_partition]),
                             'score' : merged_log_prob})

    # ----------------------------------------------------------------------------------------
    def merge_csv_files(self, fname, n_procs):
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
        if self.args.action == 'partition':  # merge partitions from several files
            self.merge_partition_files(self.hmm_outfname, n_procs)
            self.merge_csv_files(self.hmm_cachefname, n_procs)
        else:
            self.merge_csv_files(self.hmm_outfname, n_procs)

        if not self.args.no_clean:
            for iproc in range(n_procs):
                workdir = self.args.workdir + '/hmm-' + str(iproc)
                os.remove(workdir + '/' + os.path.basename(self.hmm_infname))
                os.rmdir(workdir)

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
        clust.single_link(cdr3lengthfname, debug=False)
        os.remove(cdr3lengthfname)
        return clust

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
    def get_hamming_distances(self, pairs):  #, return_info):
        # NOTE duplicates a function in utils
        return_info = []
        for query_a, query_b in pairs:
            seq_a = self.input_info[query_a]['seq']
            seq_b = self.input_info[query_b]['seq']
            if self.args.truncate_pairs:  # chop off the left side of the longer one if they're not the same length
                min_length = min(len(seq_a), len(seq_b))
                seq_a = seq_a[-min_length : ]
                seq_b = seq_b[-min_length : ]
            mutation_frac = utils.hamming(seq_a, seq_b) / float(len(seq_a))
            return_info.append({'id_a':query_a, 'id_b':query_b, 'score':mutation_frac})

        return return_info

    # ----------------------------------------------------------------------------------------
    def write_hmms(self, parameter_dir, sw_matches):
        print 'writing hmms with info from %s' % parameter_dir
        start = time.time()
        from hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')

        gene_list = self.args.only_genes
        if gene_list == None:  # if specific genes weren't specified, do the ones for which we have sw matches
            gene_list = []
            for region in utils.regions:
                for gene in self.germline_seqs[region]:
                    if gene in sw_matches:
                        gene_list.append(gene)

        for gene in gene_list:
            if self.args.debug:
                print '  %s' % utils.color_gene(gene)
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety,
                               self.germline_seqs[utils.get_region(gene)][gene],
                               self.args)
            writer.write()

        print '    time to write hmms: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, skipped_gene_matches, parameter_dir):  #, query_name, second_query_name=None):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> and print a warning """
        # first get the list of genes for which we don't have hmm files
        if len(glob.glob(parameter_dir + '/hmms/*.yaml')) == 0:
            raise Exception('ERROR no yamels in %s' % parameter_dir)

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
                print '       no %s genes in %s for %s %s' % (region, ':'.join(gene_list), query_name, '' if (second_query_name  == None) else second_query_name)
                print '          skipped %s' % (':'.join(skipped_gene_matches))
                print 'ERROR giving up on query'
                return False

        return True

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, sw_info, query_names, parameter_dir, skipped_gene_matches=None):
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
            self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir)  #, name)

            combo['only_genes'] = list(set(only_genes) | set(combo['only_genes']))  # NOTE using both sets of genes (from both query seqs) like this *really* helps,

        # self.check_hmm_existence(combo['only_genes'], skipped_gene_matches, parameter_dir, name)  # this should be superfluous now
        if not self.all_regions_present(combo['only_genes'], skipped_gene_matches, query_names):
            return {}

        return combo

    # ----------------------------------------------------------------------------------------
    def remove_sw_failures(self, query_names, sw_info):
        # if any of the queries in <query_names> was unproductive, skip the whole kitnkaboodle
        unproductive = False
        for qrn in query_names:
            if qrn in sw_info['skipped_unproductive_queries']:
                # print '    skipping unproductive %s along with %s' % (query_names[0], ' '.join(query_names[1:]))
                unproductive = True
        if unproductive:
            return []

        # otherwise they should be in sw_info, but doesn't hurt to check
        return_names = []
        for name in query_names:
            if name in sw_info:
                return_names.append(name)
            else:
                print '    %s not found in sw info' % ' '.join([str(qn) for qn in query_names])
        return return_names

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, sw_info, parameter_dir, preclusters=None, hmm_type='', shuffle_input_order=False):
        print '    writing input'
        if self.cached_results is not None:
            with opener('w')(self.hmm_cachefname) as cachefile:
                writer = csv.DictWriter(cachefile, ('unique_ids', 'score', 'naive-seq'))
                writer.writeheader()
                for uids, cachefo in self.cached_results.items():
                    writer.writerow({'unique_ids':uids, 'score':cachefo['logprob'], 'naive-seq':cachefo['naive-seq']})

        csvfile = opener('w')(self.hmm_infname)
        # start = time.time()

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
            assert preclusters is not None
            if self.args.action == 'partition':
                nsets = preclusters.best_minus_ten_partition
            else:
                nsets = [val for val in preclusters.id_clusters.values() if len(val) > 1]  # <nsets> is a list of sets (well, lists) of query names
        elif hmm_type == 'k=nsets':
            if self.args.all_combinations:  # run on *every* combination of queries which has length <self.args.n_sets>
                nsets = itertools.combinations(self.input_info.keys(), self.args.n_sets)
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
        else:
            assert False

        if shuffle_input_order:  # This is *really* important if you're partitioning in parallel
            random_nsets = []
            while len(nsets) > 0:
                irand = random.randint(0, len(nsets) - 1)  # NOTE interval is inclusive
                random_nsets.append(nsets[irand])
                nsets.remove(nsets[irand])
            nsets = random_nsets

        for query_names in nsets:
            non_failed_names = self.remove_sw_failures(query_names, sw_info)
            if len(non_failed_names) == 0:
                continue
            combined_query = self.combine_queries(sw_info, non_failed_names, parameter_dir, skipped_gene_matches=skipped_gene_matches)
            if len(combined_query) == 0:  # didn't find all regions
                continue
            csvfile.write('%s %d %d %d %d %s %s\n' %  # NOTE csv.DictWriter can handle tsvs, so this should really be switched to use that
                          (':'.join([str(qn) for qn in non_failed_names]),
                           combined_query['k_v']['min'], combined_query['k_v']['max'],
                           combined_query['k_d']['min'], combined_query['k_d']['max'],
                           ':'.join(combined_query['only_genes']),
                           ':'.join(combined_query['seqs'])))

        if len(skipped_gene_matches) > 0:
            print '    not found in %s, i.e. were never the best sw match for any query, so removing from consideration for hmm:' % (parameter_dir)
            for region in utils.regions:
                print '      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in skipped_gene_matches if utils.get_region(gene) == region]))

        csvfile.close()
        # print '        input write time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def read_partition_info(self, fname):
        partition_info = []
        with opener('r')(fname) as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                if line['partition'] == '':
                    raise Exception('ERROR null partition (one of the processes probably got passed zero sequences')  # shouldn't happen any more
                uids = []
                for cluster in line['partition'].split(';'):
                    try:
                        uids.append([int(unique_id) for unique_id in cluster.split(':')])  # TODO remove this try/except bullshit
                    except ValueError:
                        uids.append([unique_id for unique_id in cluster.split(':')])
                partition_info.append({'clusters':uids, 'score':float(line['score'])})
        return partition_info

    # ----------------------------------------------------------------------------------------
    def read_partition_output(self):
        partition_info = self.read_partition_info(self.hmm_outfname)

        if self.cached_results is None:
            self.cached_results = {}
        with opener('r')(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                if line['unique_ids'] not in self.cached_results:
                    self.cached_results[line['unique_ids']] = {'logprob':float(line['score']), 'naive-seq':line['naive-seq']}

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)
            os.remove(self.hmm_cachefname)

        return partition_info

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, algorithm, count_parameters=False, parameter_out_dir=None, plotdir=None):
        print '    read output'

        if count_parameters:
            assert parameter_out_dir is not None
            assert plotdir is not None
        pcounter = ParameterCounter(self.germline_seqs) if count_parameters else None
        true_pcounter = ParameterCounter(self.germline_seqs) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter(self.germline_seqs, plotdir + '/hmm/performance', 'hmm') if self.args.plot_performance else None

        n_processed = 0
        hmminfo = OrderedDict()
        with opener('r')(self.hmm_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            last_key = None
            boundary_error_queries = []
            for line in reader:
                utils.intify(line, splitargs=('unique_ids', 'seqs'))
                ids = line['unique_ids']
                this_key = utils.get_key(ids)
                same_event = utils.from_same_event(self.args.is_data, self.reco_info, ids)
                if same_event is None:
                    same_event = -1
                id_str = ''.join(['%20s ' % i for i in ids])

                # check for errors
                if last_key != this_key:  # if this is the first line for this set of ids (i.e. the best viterbi path or only forward score)
                    if line['errors'] != None and 'boundary' in line['errors'].split(':'):
                        boundary_error_queries.append(':'.join([str(uid) for uid in ids]))
                    else:
                        assert len(line['errors']) == 0

                if algorithm == 'viterbi':
                    line['seq'] = line['seqs'][0]  # add info for the best match as 'seq'
                    line['unique_id'] = ids[0]
                    utils.add_match_info(self.germline_seqs, line, self.cyst_positions, self.tryp_positions, debug=(self.args.debug > 0))

                    if last_key != this_key or self.args.plot_all_best_events:  # if this is the first line (i.e. the best viterbi path) for this query (or query pair), print the true event
                        n_processed += 1
                        if self.args.debug:
                            print '%s   %d' % (id_str, same_event)
                        if line['cdr3_length'] != -1 or not self.args.skip_unproductive:  # if it's productive, or if we're not skipping unproductive rearrangements
                            hmminfo[':'.join([str(uid) for uid in line['unique_ids']])] = line
                            if pcounter is not None:  # increment counters (but only for the best [first] match)
                                pcounter.increment(line)
                            if true_pcounter is not None:  # increment true counters
                                true_pcounter.increment(self.reco_info[ids[0]])
                            if perfplotter is not None:
                                perfplotter.evaluate(self.reco_info[ids[0]], line)

                    if self.args.debug:
                        self.print_hmm_output(line, print_true=(last_key != this_key))  #, perfplotter=perfplotter)
                    line['seq'] = None
                    line['unique_id'] = None

                else:  # for forward, write the pair scores to file to be read by the clusterer
                    print '%3d %10.3f    %s' % (same_event, float(line['score']), id_str)
                    if line['score'] == '-nan':
                        print '    WARNING encountered -nan, setting to -999999.0'
                        score = -999999.0
                    else:
                        score = float(line['score'])
                    if len(ids) == 2:
                        raise Exception('DEPRECATED')
                        # hmminfo.append({'id_a':line['unique_ids'][0], 'id_b':line['unique_ids'][1], 'score':score})
                    n_processed += 1

                last_key = utils.get_key(ids)

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

        print '  processed %d queries' % n_processed
        if len(boundary_error_queries) > 0:
            print '    %d boundary errors (%s)' % (len(boundary_error_queries), ', '.join(boundary_error_queries))

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)

        return hmminfo

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
                    out_str_list.append(utils.print_reco_event(self.germline_seqs, self.reco_info[uids[iid]], extra_str='    ', return_string=True, label='true:', one_line=(iid != 0)))
            ilabel = 'inferred:'

        out_str_list.append(utils.print_reco_event(self.germline_seqs, line, extra_str='    ', return_string=True, label=ilabel))
        for iextra in range(1, len(line['unique_ids'])):
            line['seq'] = line['seqs'][iextra]
            out_str_list.append(utils.print_reco_event(self.germline_seqs, line, extra_str='    ', return_string=True, one_line=True))

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
