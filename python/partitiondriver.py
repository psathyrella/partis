import time
import sys
import json
import itertools
import shutil
import math
import os
import glob
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import random
from collections import OrderedDict
from subprocess import Popen, check_call, PIPE, check_output

import utils
from opener import opener
from seqfileopener import get_seqfile_info
import annotationclustering
from glomerator import Glomerator
from clusterpath import ClusterPath
from waterer import Waterer
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter
from hist import Hist

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    """ Class to parse input files, start bcrham jobs, and parse/interpret bcrham output for annotation and partitioning """
    def __init__(self, args):
        self.args = args
        self.germline_seqs = utils.read_germlines(self.args.datadir)
        self.cyst_positions = utils.read_cyst_positions(self.args.datadir)
        with opener('r')(self.args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

        if self.args.seqfile is not None:
            self.input_info, self.reco_info = get_seqfile_info(self.args.seqfile, self.args.is_data, self.germline_seqs, self.cyst_positions, self.tryp_positions,
                                                               self.args.n_max_queries, self.args.queries, self.args.reco_ids)

        self.cached_results = None
        self.bcrham_divvied_queries = None
        self.n_max_divvy = 100  # if input info is longer than this, divvy with bcrham
        self.n_likelihoods_calculated = None
        self.n_max_calc_per_process = 250

        self.sw_info = None

        utils.prep_dir(self.args.workdir)
        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'

        if self.args.outfname is not None:
            outdir = os.path.dirname(self.args.outfname)
            if outdir != '' and not os.path.exists(outdir):
                os.makedirs(outdir)

        if self.args.persistent_cachefname is not None and os.path.exists(self.args.persistent_cachefname):
            check_call(['cp', '-v', self.args.persistent_cachefname, self.hmm_cachefname])

    # ----------------------------------------------------------------------------------------
    def clean(self):
        # merge persistent and current cache files into the persistent cache file
        if self.args.persistent_cachefname is not None:
            lockfname = self.args.persistent_cachefname + '.lock'
            while os.path.exists(lockfname):
                print '  waiting for lock on %s' % lockfname
                time.sleep(0.5)
            lockfile = open(lockfname, 'w')
            if not os.path.exists(self.args.persistent_cachefname):
                open(self.args.persistent_cachefname, 'w').close()
            tmp_cachefname = self.args.persistent_cachefname + '.tmp'
            check_call(['cp', self.args.persistent_cachefname, tmp_cachefname])
            self.merge_files(infnames=[tmp_cachefname, self.hmm_cachefname], outfname=self.args.persistent_cachefname)
            lockfile.close()
            os.remove(lockfname)
        if not self.args.no_clean and os.path.exists(self.hmm_cachefname):
            os.remove(self.hmm_cachefname)

        if not self.args.no_clean:
            try:
                os.rmdir(self.args.workdir)
            except OSError:
                raise Exception('workdir (%s) not empty: %s' % (self.args.workdir, ' '.join(os.listdir(self.args.workdir))))  # hm... you get weird recursive exceptions if you get here. Oh, well, it still works

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        """ Infer full parameter sets and write hmm files for sequences from <self.input_info>, first with Smith-Waterman, then using the SW output as seed for the HMM """
        sw_parameter_dir = self.args.parameter_dir + '/sw'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True)
        waterer.run()
        self.sw_info = waterer.info
        self.write_hmms(sw_parameter_dir)
        parameter_out_dir = self.args.parameter_dir + '/hmm'
        self.run_hmm('viterbi', parameter_in_dir=sw_parameter_dir, parameter_out_dir=parameter_out_dir, count_parameters=True)
        self.write_hmms(parameter_out_dir)

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self, algorithm):
        """ Just run <algorithm> (either 'forward' or 'viterbi') on sequences in <self.input_info> and exit. You've got to already have parameters cached in <self.args.parameter_dir> """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir (' + self.args.parameter_dir + ') d.n.e')
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.sw_info = waterer.info
        self.run_hmm(algorithm, parameter_in_dir=self.args.parameter_dir)

    # ----------------------------------------------------------------------------------------
    # get number of clusters based on sum of last paths in <self.smc_info>
    def get_n_clusters(self):
        if self.args.smc_particles == 1:
            return len(self.paths[-1].partitions[self.paths[-1].i_best_minus_x])

        nclusters = 0
        for iproc in range(len(self.smc_info[-1])):  # number of processes
            path = self.smc_info[-1][iproc][0]  # uses the first smc particle, but the others will be similar
            nclusters += len(path.partitions[path.i_best_minus_x])
        return nclusters

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir %s d.n.e.' % self.args.parameter_dir)

        if self.args.print_partitions:
            cp = ClusterPath(-1)
            cp.readfile(self.args.outfname)
            cp.print_partitions(reco_info=self.reco_info, one_line=True, abbreviate=True, n_to_print=100)
            return

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()
        self.sw_info = waterer.info
        if not self.args.dont_pad_sequences:  # have to do this before we divvy... sigh TODO clean this up and only call it in one place
            # TODO oh, duh, just do it in waterer
            self.pad_seqs_to_same_length()  # adds padded info to sw_info

        n_procs = self.args.n_procs
        n_proc_list = []  # list of the number of procs we used for each run

        # add initial lists of paths
        if self.args.smc_particles == 1:
            cp = ClusterPath(-1)
            cp.add_partition([[cl, ] for cl in self.input_info.keys()], logprob=0., n_procs=n_procs)
            self.paths = [cp, ]
        else:
            initial_divvied_queries = self.divvy_up_queries(n_procs, [line for line in self.sw_info.values() if 'unique_id' in line], 'unique_id', 'seq')
            self.smc_info = [[], ]
            for clusters in initial_divvied_queries:  # one set of <clusters> for each process
                self.smc_info[-1].append([])
                for iptl in range(self.args.smc_particles):
                    cp = ClusterPath(-1)
                    cp.add_partition([[cl, ] for cl in clusters], logprob=0., n_procs=n_procs)
                    self.smc_info[-1][-1].append(cp)

        # cache hmm naive seqs for each single query
        if len(self.input_info) > 50 or self.args.naive_vsearch or self.args.naive_swarm:
            self.run_hmm('viterbi', self.args.parameter_dir, n_procs=int(math.ceil(float(len(self.input_info)) / 50)), cache_naive_seqs=True)

        if self.args.naive_vsearch or self.args.naive_swarm:
            self.cluster_with_naive_vsearch_or_swarm(self.args.parameter_dir)
            return

        # run that shiznit
        while n_procs > 0:
            start = time.time()
            nclusters = self.get_n_clusters()
            print '--> %d clusters with %d procs' % (nclusters, n_procs)  # write_hmm_input uses the best-minus-ten partition
            self.run_hmm('forward', self.args.parameter_dir, n_procs=n_procs, divvy_with_bcrham=(self.get_n_clusters() > self.n_max_divvy and not self.args.random_divvy))
            n_proc_list.append(n_procs)

            print '      partition step time: %.3f' % (time.time()-start)
            if n_procs == 1 or len(n_proc_list) >= self.args.n_partition_steps:
                break

            if self.args.smc_particles == 1:  # for smc, we merge pairs of processes; otherwise, we do some heuristics to come up with a good number of clusters for the next iteration
                n_calcd_per_process = self.get_n_calculated_per_process()
                if n_calcd_per_process < self.n_max_calc_per_process:
                    # if n_procs > 4 or n_proc_list[-1] == n_proc_list[-2]:
                    n_procs = int(n_procs / 1.2)
            else:
                n_procs = len(self.smc_info[-1])  # if we're doing smc, the number of particles is determined by the file merging process

        # deal with final partition
        if self.args.smc_particles == 1:
            for path in self.paths:
                self.check_path(path)
            if self.args.debug:
                print 'final'
                for ipath in range(len(self.paths)):  # one path for each glomeration step
                    self.paths[ipath].print_partitions(self.reco_info, one_line=True, print_header=(ipath==0))
                    print ''
            if self.args.outfname is not None:
                self.write_partitions(self.args.outfname, [self.paths[-1], ])  # [last agglomeration step]
        else:
            # self.merge_pairs_of_procs(1)  # DAMMIT why did I have this here? I swear there was a reason but I can't figure it out, and it seems to work without it
            final_paths = self.smc_info[-1][0]  # [last agglomeration step][first (and only) process in the last step]
            for path in final_paths:
                self.check_path(path)
            if self.args.debug:
                for ipath in range(self.args.smc_particles):
                    path = final_paths[ipath]
                    path.print_partition(path.i_best, self.reco_info, extrastr=str(ipath) + ' final')
            if self.args.outfname is not None:
                self.write_partitions(self.args.outfname, final_paths)

        if self.args.debug and not self.args.is_data:
            tmpglom = Glomerator(self.reco_info)
            tmpglom.print_true_partition()

    # ----------------------------------------------------------------------------------------
    def check_path(self, path):
        missing_ids = set()
        def check_partition(partition):
            for uid in self.input_info:
                found = False
                for cluster in partition:
                    if uid in cluster:
                        found = True
                        break
                if not found and uid not in self.sw_info['skipped_unproductive_queries']:
                    missing_ids.add(uid)
            for cluster in partition:
                for uid in cluster:
                    if uid not in self.input_info:
                        missing_ids.add(uid)

        check_partition(path.partitions[path.i_best])
        # for ipart in range(len(path.partitions)):
        #     check_partition(path.partitions[ipart])

        if len(missing_ids) > 0:
            print 'WARNING not found in merged partitions: ' + ' '.join(missing_ids)

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, outfname, paths):
        with opener('w')(outfname) as outfile:
            headers = ['logprob', 'n_clusters', 'n_procs', 'clusters']
            if self.args.smc_particles > 1:
                headers += ['path_index', 'logweight']
            if not self.args.is_data:
                headers += ['adj_mi', 'n_true_clusters', 'bad_clusters']
            writer = csv.DictWriter(outfile, headers)
            writer.writeheader()
            true_partition = None if self.args.is_data else utils.get_true_partition(self.reco_info)
            for ipath in range(len(paths)):
                paths[ipath].write_partitions(writer, self.args.is_data, self.reco_info, true_partition, self.args.smc_particles, path_index=self.args.seed + ipath, n_to_write=self.args.n_partitions_to_write, calc_adj_mi='best')

    # ----------------------------------------------------------------------------------------
    def cluster_with_naive_vsearch_or_swarm(self, parameter_dir):  # TODO change name of function if you switch to just swarm
        start = time.time()
        # read cached naive seqs
        naive_seqs = {}
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                unique_ids = line['unique_ids'].split(':')
                assert len(unique_ids) == 1
                unique_id = unique_ids[0]
                naive_seqs[unique_id] = line['naive_seq']

        # make a fasta file
        fastafname = self.args.workdir + '/simu.fasta'
                
        # if not os.path.exists(fastafname):
        if self.args.naive_swarm:
            print '    NOTE: replacing N with A for input to swarm'
        with open(fastafname, 'w') as fastafile:
            for query, naive_seq in naive_seqs.items():
                if self.args.naive_swarm:
                    query += '_1'
                    naive_seq = utils.remove_ambiguous_ends(naive_seq, '', '')
                    naive_seq = naive_seq.replace('N', 'A')
                fastafile.write('>' + query + '\n' + naive_seq + '\n')

        if self.args.naive_vsearch:
            bound = self.get_naive_hamming_auto_bounds(parameter_dir) /  2.  # yay for heuristics! (I did actually optimize this...)
            id_fraction = 1. - bound
            clusterfname = self.args.workdir + '/vsearch-clusters.txt'
            cmd = './bin/vsearch-1.1.3-linux-x86_64 --uc ' + clusterfname + ' --cluster_fast ' + fastafname + ' --id ' + str(id_fraction) + ' --maxaccept 0 --maxreject 0'
            # print cmd
            check_call(cmd.split())
    
        elif self.args.naive_swarm:
            clusterfname = self.args.workdir + '/vsearch-clusters.txt'
            cmd = './bin/swarm-2.1.1-linux-x86_64 ' + fastafname
            cmd += ' -t 5'  # five threads TODO set this more intelligently
            cmd += ' -f'
            cmd += ' --match-reward ' + str(self.args.match_mismatch[0])
            cmd += ' --mismatch-penalty ' + str(self.args.match_mismatch[1])
            cmd += ' --gap-opening-penalty ' + str(self.args.gap_open_penalty)
            # cmd += ' --gap-extension-penalty'
            cmd += ' --uclust-file ' + clusterfname
            print cmd
            check_call(cmd.split())
        else:
            assert False

        # read output
        id_clusters = {}
        with open(clusterfname) as clusterfile:
            reader = csv.DictReader(clusterfile, fieldnames=['type', 'cluster_id', '3', '4', '5', '6', '7', 'crap', 'query', 'morecrap'], delimiter='\t')
            for line in reader:
                if line['type'] == 'C':  # batshit output format: some lines are a cluster, and some are a query sequence. Skip the cluster ones.
                    continue
                cluster_id = int(line['cluster_id'])
                if cluster_id not in id_clusters:
                    id_clusters[cluster_id] = []
                uid = line['query']
                if self.args.naive_swarm and uid[-2:] == '_1':  # remove (dummy) abundance information
                    uid = uid[:-2]
                id_clusters[cluster_id].append(uid)
        partition = id_clusters.values()
        adj_mi = -1
        if not self.args.is_data:
            adj_mi = utils.mutual_information(partition, self.reco_info, debug=True)
        cp = ClusterPath(-1)
        cp.add_partition(partition, logprob=0.0, n_procs=1, adj_mi=adj_mi)
        if self.args.outfname is not None:
            self.write_partitions(self.args.outfname, [cp, ])

        if not self.args.no_clean:
            os.remove(fastafname)
            os.remove(clusterfname)

        print '      vsearch/swarm time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def get_naive_hamming_auto_bounds(self, parameter_dir):
        mutehist = Hist(fname=parameter_dir + '/all-mean-mute-freqs.csv')
        mute_freq = mutehist.get_mean()
        # just use a line based on two points (mute_freq, threshold)
        x1, x2 = 0.1, 0.25
        y1, y2 = 0.04, 0.08
        m = (y2 - y1) / (x2 - x1);
        b = 0.5 * (y1 + y2 - m*(x1 + x2));
        # for x in [x1, x2]:
        #     print '%f x + %f = %f' % (m, b, m*x + b)
        bound = m * mute_freq + b
        return bound

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
        cmd_str += ' --datadir ' + os.getcwd() + '/' + self.args.datadir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        cmd_str += ' --max-logprob-drop ' + str(self.args.max_logprob_drop)

        if self.args.auto_hamming_fraction_bounds:
            bound = self.get_naive_hamming_auto_bounds(parameter_dir)
            print '    auto naive hamming bound %f' % bound
            self.args.hamming_fraction_bounds = [bound, bound]
            cmd_str += ' --no-fwd'  # assume that auto hamming bounds means we're naive hamming clustering (which is a good assumption, since we set the lower and upper bounds to the same thing)

        cmd_str += ' --hamming-fraction-bound-lo ' + str(self.args.hamming_fraction_bounds[0])
        cmd_str += ' --hamming-fraction-bound-hi ' + str(self.args.hamming_fraction_bounds[1])
        if self.args.smc_particles > 1:
            os.environ['GSL_RNG_TYPE'] = 'ranlux'
            os.environ['GSL_RNG_SEED'] = str(random.randint(0, 99999))
            cmd_str += ' --smc-particles ' + str(self.args.smc_particles)
        if self.args.rescale_emissions:
            cmd_str += ' --rescale-emissions'
        if self.args.action == 'partition':
            cmd_str += ' --partition'
            cmd_str += ' --cachefile ' + self.hmm_cachefname
        if self.args.truncate_n_sets:
            cmd_str += ' --truncate-seqs'
        if not self.args.dont_allow_unphysical_insertions:
            cmd_str += ' --unphysical-insertions'
        assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
        cmd_str += ' --ambig-base ' + utils.ambiguous_bases[0]

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def execute_iproc(self, cmd_str):
        proc = Popen(cmd_str.split(), stdout=PIPE, stderr=PIPE)
        return proc

    # ----------------------------------------------------------------------------------------
    def get_n_calculated_per_process(self):
        if self.n_likelihoods_calculated is None:
            return
        print 'n calcd:'
        total = 0
        for procinfo in self.n_likelihoods_calculated:
            print '    %d' % procinfo['fwd']
            total += procinfo['fwd']
        print '  total: %d (%.1f per proc)' % (total, float(total) / len(self.n_likelihoods_calculated))
        return float(total) / len(self.n_likelihoods_calculated)

    # ----------------------------------------------------------------------------------------
    def execute(self, cmd_str, n_procs, total_naive_hamming_cluster_procs=None):
        print '    running'
        start = time.time()
        if n_procs == 1:
            if total_naive_hamming_cluster_procs is not None:
                cmd_str = cmd_str.replace('XXX', str(total_naive_hamming_cluster_procs))
            # print cmd_str
            # sys.exit()
            check_call(cmd_str.split())
        else:
            if total_naive_hamming_cluster_procs is not None:
                n_leftover = total_naive_hamming_cluster_procs - (total_naive_hamming_cluster_procs / n_procs) * n_procs
            cmd_strs, workdirs = [], []
            for iproc in range(n_procs):
                workdirs.append(self.args.workdir + '/hmm-' + str(iproc))
                cmd_strs.append(cmd_str.replace(self.args.workdir, workdirs[-1]))
                if total_naive_hamming_cluster_procs is not None:
                    clusters_this_proc = total_naive_hamming_cluster_procs / n_procs
                    if n_leftover > 0:
                        clusters_this_proc += 1
                        n_leftover -= 1
                    cmd_strs[-1] = cmd_strs[-1].replace('XXX', str(clusters_this_proc))
                # print cmd_strs[-1]
                # sys.exit()
            procs, n_tries = [], []
            self.n_likelihoods_calculated = []
            for iproc in range(n_procs):
                procs.append(self.execute_iproc(cmd_strs[iproc]))
                n_tries.append(1)
                self.n_likelihoods_calculated.append({})
            while procs.count(None) != len(procs):
                for iproc in range(n_procs):
                    if procs[iproc] is None:
                        continue
                    out, err = procs[iproc].communicate()
                    utils.process_out_err(out, err, extra_str=str(iproc), info=self.n_likelihoods_calculated[iproc])
                    outfname = self.hmm_outfname.replace(self.args.workdir, workdirs[iproc])
                    if os.path.exists(outfname):  # TODO also check cachefile, if necessary
                        procs[iproc] = None  # job succeeded
                    elif n_tries[iproc] > 5:
                        raise Exception('exceeded max number of tries for command\n    %s\nlook for output in %s' % (cmd_strs[iproc], workdirs[iproc]))
                    else:
                        print '    rerunning proc %d' % iproc
                        procs[iproc] = self.execute_iproc(cmd_strs[iproc])
                        n_tries[iproc] += 1
                        time.sleep(0.1)

        sys.stdout.flush()
        print '      hmm run time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, parameter_in_dir, parameter_out_dir='', count_parameters=False, n_procs=None, cache_naive_seqs=False, divvy_with_bcrham=False):
        """ 
        Run bcrham, possibly with many processes, and parse and interpret the output.
        NOTE the local <n_procs>, which overrides the one from <self.args>
        """
        print 'hmm'
        if n_procs is None:
            n_procs = self.args.n_procs
        if n_procs < 1 or n_procs > 9999:
            raise Exception('bad n_procs %s' % n_procs)

        # if not naive_hamming_cluster:  # should already be there
        self.write_hmm_input(parameter_dir=parameter_in_dir)  # TODO don't keep rewriting it

        # sys.stdout.flush()
        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir)
        if cache_naive_seqs:
            print '      caching all naive sequences'
            assert '--partition' in cmd_str
            cmd_str = cmd_str.replace('--partition', '--cache-naive-seqs')

        if n_procs > 1 and self.args.smc_particles == 1:  # if we're doing smc (i.e. if > 1), we have to split things up more complicatedly elsewhere
            if divvy_with_bcrham:
                print '      bcrham naive hamming clustering'
                assert '--partition' in cmd_str and algorithm == 'forward'
                n_divvy_procs = max(1, self.get_n_clusters() / 500)  # number of bcrham procs used to divvy up queries with naive hamming clustering
                self.split_input(n_divvy_procs, self.hmm_infname, 'hmm', algorithm, cache_naive_seqs, bcrham_naive_hamming_cluster=True)
                self.execute(cmd_str.replace('--partition', '--naive-hamming-cluster XXX'), n_procs=n_divvy_procs, total_naive_hamming_cluster_procs=n_procs)
                self.read_naive_hamming_clusters(n_procs=n_divvy_procs)
            self.split_input(n_procs, self.hmm_infname, 'hmm', algorithm, cache_naive_seqs, bcrham_naive_hamming_cluster=False)

        self.execute(cmd_str, n_procs)

        self.read_hmm_output(algorithm, n_procs, count_parameters, parameter_out_dir, cache_naive_seqs)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self):
        """ a.t.m. just want to know which values we have """
        cachefo = {}
        if not os.path.exists(self.hmm_cachefname):
            return cachefo
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                cachefo[line['unique_ids']] = {}
        return cachefo

    # ----------------------------------------------------------------------------------------
    def get_expected_number_of_forward_calculations(self, info, namekey, seqkey):
        start = time.time()
        def join_names(name1, name2):  # mimics function in glomeraor.cc
            sortedlist = sorted([name1, name2])
            return ':'.join(sortedlist)

        naive_seqs = self.get_naive_seqs(info, namekey, seqkey)
        cachefo = self.read_cachefile()
        n_total, n_cached = 0, 0
        for id_a, id_b in itertools.combinations(naive_seqs.keys(), 2):
            seq_a, seq_b = naive_seqs[id_a], naive_seqs[id_b]
            hfrac = utils.hamming_fraction(seq_a, seq_b)
            if hfrac >= self.args.hamming_fraction_bounds[0] and hfrac <= self.args.hamming_fraction_bounds[1]:  # NOTE not sure the equals match up exactly with what's in ham, but it's an estimate, so it doesn't matter
                n_total += 1
                if join_names(id_a, id_b) in cachefo:
                    n_cached += 1
                    assert ':'.join(sorted([id_a, id_b], reverse=True)) not in cachefo
                    assert id_a in cachefo
                    assert id_b in cachefo

        print 'expected total: %d  (cached: %d) --> %d' % (n_total, n_cached, n_total - n_cached)
        print '      expected calc time: %.3f' % (time.time()-start)
        return n_total - n_cached

    # ----------------------------------------------------------------------------------------
    def get_naive_seqs(self, info, namekey, seqkey):
        def get_query_from_sw(qry):
            assert qry in self.sw_info
            naive_seq = utils.get_full_naive_seq(self.germline_seqs, self.sw_info[qry])
            if not self.args.dont_pad_sequences:
                padleft = self.sw_info[qry]['padded']['padleft']  # we're padding the *naive* seq corresponding to qry now, but it'll be the same length as the qry seq
                padright = self.sw_info[qry]['padded']['padright']
                assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
                naive_seq = padleft * utils.ambiguous_bases[0] + naive_seq + padright * utils.ambiguous_bases[0]
            return naive_seq

        naive_seqs = {}
        for line in info:
            query = line[namekey]
            seqstr = line['padded'][seqkey] if 'padded' in line else line[seqkey]
            # NOTE cached naive seqs should all be the same length
            if self.cached_results is not None and seqstr in self.cached_results and self.cached_results[seqstr]['naive_seq'] is not None:  # first try to used cached hmm results
                naive_seqs[query] = self.cached_results[seqstr]['naive_seq']
            elif len(query.split(':')) == 1:  # ...but if we don't have them, use smith-waterman (should only be for single queries)
               naive_seqs[query] = get_query_from_sw(query)
            elif len(query.split(':')) > 1:
                naive_seqs[query] = get_query_from_sw(query.split(':')[0])  # just arbitrarily use the naive seq from the first one. This is ok partly because if we cache the logprob but not the naive seq, that's because we thought about merging two clusters but did not -- so they're naive seqs should be similar. Also, this is just for divvying queries.
            else:
                raise Exception('no naive sequence found for ' + str(query))
            if naive_seqs[query] == '':
                raise Exception('zero-length naive sequence found for ' + str(query))
        return naive_seqs

    # ----------------------------------------------------------------------------------------
    def divvy_up_queries(self, n_procs, info, namekey, seqkey, debug=True):
        if self.bcrham_divvied_queries is not None:
            print 'using bcrham_divvied_queries'
            if len(self.bcrham_divvied_queries) != n_procs:
                raise Exception('Wrong number of clusters %d for %d procs' % (len(self.bcrham_divvied_queries), n_procs))
            return self.bcrham_divvied_queries

        print 'no bcrham divvies, divvying with python glomerator'
        naive_seqs = self.get_naive_seqs(info, namekey, seqkey)

        if self.args.truncate_n_sets:
            assert False  # deprecated and broken
            # print '  truncate in divvy'
            # self.truncate_seqs(naive_seqs, kvinfo=None, cyst_positions=cyst_positions)

        clust = Glomerator()
        divvied_queries = clust.naive_seq_glomerate(naive_seqs, n_clusters=n_procs)
        if debug:
            print '  divvy lengths'
            for dq in divvied_queries:
                print '  ', len(dq),
            print ''

        if len(divvied_queries) != n_procs:
            raise Exception('Wrong number of clusters %d for %d procs' % (len(divvied_queries), n_procs))

        return divvied_queries

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname, prefix, algorithm, cache_naive_seqs, bcrham_naive_hamming_cluster):
        """ Do stuff. Probably correctly. """
        assert self.args.smc_particles == 1

        # read single input file
        info = []
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile, delimiter=' ')
            for line in reader:
                info.append(line)

        # initialize output files
        sub_outfiles, writers = [], []
        for iproc in range(n_procs):
            subworkdir = self.args.workdir + '/' + prefix + '-' + str(iproc)
            utils.prep_dir(subworkdir)
            # prep each suboutput file
            sub_outfiles.append(opener('w')(subworkdir + '/' + os.path.basename(infname)))
            writers.append(csv.DictWriter(sub_outfiles[-1], reader.fieldnames, delimiter=' '))
            writers[-1].writeheader()
            # copy cachefile to this subdir
            if os.path.exists(self.hmm_cachefname):
                check_call(['cp', self.hmm_cachefname, subworkdir + '/'])

        cluster_divvy = False
        if self.args.action == 'partition' and algorithm == 'forward':
            if not self.args.random_divvy and not cache_naive_seqs and not bcrham_naive_hamming_cluster:
                cluster_divvy = True
        # self.get_expected_number_of_forward_calculations(info, 'names', 'seqs')
        if cluster_divvy:  # cluster similar sequences together (otherwise just do it in order)
            print 'cluster divvy in split_input'
            divvied_queries = self.divvy_up_queries(n_procs, info, 'names', 'seqs')
        else:
            print 'modulo divvy'
        for iproc in range(n_procs):
            for iquery in range(len(info)):
                if cluster_divvy:
                    if info[iquery]['names'] not in divvied_queries[iproc]:  # NOTE I think the reason this doesn't seem to be speeding things up is that our hierarhical agglomeration time is dominated by the distance calculation, and that distance calculation time is roughly proportional to the number of sequences in the cluster (i.e. larger clusters take longer)
                        continue
                else:
                    if iquery % n_procs != iproc:
                        continue
                writers[iproc].writerow(info[iquery])

        for iproc in range(n_procs):
            sub_outfiles[iproc].close()

        if self.bcrham_divvied_queries is not None:
            self.bcrham_divvied_queries = None

    # ----------------------------------------------------------------------------------------
    def merge_subprocess_files(self, fname, n_procs, csv_header=True):
        subfnames = []
        for iproc in range(n_procs):
            subfnames.append(self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(fname))
        self.merge_files(subfnames, fname, csv_header)

    # ----------------------------------------------------------------------------------------
    def merge_files(self, infnames, outfname, csv_header=True):
        """ 
        Merge <infnames> into <outfname>.
        NOTE that <outfname> is overwritten with the zero-length file if it exists, otherwise it is created.
        Some of <infnames> may not exist.
        """
        assert outfname not in infnames
        start = time.time()

        header = ''
        with open(outfname, 'w') as outfile:
            if csv_header:  # if not <csv_header>, we'll just end up with a zero-length file
                for fname in infnames:
                    if not os.path.exists(fname) or os.stat(fname).st_size == 0:
                        continue
                    with open(fname) as headfile:
                        reader = csv.DictReader(headfile)
                        writer = csv.DictWriter(outfile, reader.fieldnames)
                        writer.writeheader()
                        header = ','.join(reader.fieldnames)
                    break

        cmd = 'cat ' + ' '.join(infnames) + ' | grep -v \'' + header + '\' | sort | uniq >>' + outfname
        check_call(cmd, shell=True)

        if not self.args.no_clean:
            for infname in infnames:
                os.remove(infname)

        print '    time to merge csv files: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def merge_all_hmm_outputs(self, n_procs, cache_naive_seqs):
        """ Merge any/all output files from subsidiary bcrham processes (used when *not* doing smc) """
        assert self.args.smc_particles == 1  # have to do things more complicatedly for smc
        if self.args.action == 'partition':  # merge partitions from several files
            if n_procs > 1:
                self.merge_subprocess_files(self.hmm_cachefname, n_procs)

            if not cache_naive_seqs:
                if n_procs == 1:
                    infnames = [self.hmm_outfname, ]
                else:
                    infnames = [self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(self.hmm_outfname) for iproc in range(n_procs)]
                previous_info = None
                if len(self.paths) > 1:
                    previous_info = self.paths[-1]
                glomerer = Glomerator(self.reco_info)
                glomerer.read_cached_agglomeration(infnames, smc_particles=1, previous_info=previous_info, calc_adj_mi=self.args.debug, debug=self.args.debug)  #, outfname=self.hmm_outfname)
                assert len(glomerer.paths) == 1
                # self.check_path(glomerer.paths[0])
                self.paths.append(glomerer.paths[0])
        else:
            self.merge_subprocess_files(self.hmm_outfname, n_procs)

        if not self.args.no_clean:
            if n_procs == 1:
                # print 'removing ', self.hmm_outfname
                os.remove(self.hmm_outfname)
            else:
                for iproc in range(n_procs):
                    subworkdir = self.args.workdir + '/hmm-' + str(iproc)
                    os.remove(subworkdir + '/' + os.path.basename(self.hmm_infname))
                    if os.path.exists(subworkdir + '/' + os.path.basename(self.hmm_outfname)):
                        # print 'removing ', subworkdir + '/' + os.path.basename(self.hmm_outfname)
                        os.remove(subworkdir + '/' + os.path.basename(self.hmm_outfname))
                    os.rmdir(subworkdir)

    # ----------------------------------------------------------------------------------------
    def merge_pairs_of_procs(self, n_procs):
        """ Merge the output from pairs of processes (used when doing smc)"""
        assert self.args.action == 'partition'
        assert self.args.smc_particles > 1
        if n_procs > 1:
            groups_to_merge = [[i, i+1] for i in range(0, n_procs-1, 2)]  # e.g. for n_procs = 5, we merge the groups [0, 1], [2, 3, 4]
        else:
            groups_to_merge = [[], ]
        if n_procs % 2 != 0:  # if it's odd, add the last proc to the last group
            groups_to_merge[-1].append(n_procs-1)
        self.smc_info.append([])
        for group in groups_to_merge:
            if n_procs == 1:
                infnames = [self.hmm_outfname, ]
            else:
                infnames = [self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(self.hmm_outfname) for iproc in group]
            assert len(self.smc_info[-2]) == n_procs
            previous_info = None
            if len(self.smc_info) > 2:
                previous_info = [self.smc_info[-2][iproc] for iproc in group]
            glomerer = Glomerator(self.reco_info)
            paths = glomerer.read_cached_agglomeration(infnames, self.args.smc_particles, previous_info=previous_info, calc_adj_mi=self.args.debug, debug=self.args.debug)  #, outfname=self.hmm_outfname)
            self.smc_info[-1].append(paths)

            # ack? self.glomclusters.append(glomerer)
            # boof? self.list_of_preclusters.append(glomerer.combined_conservative_best_minus_ten_partitions)

        if n_procs > 1:  # TODO I don't think this is right any more...
            self.merge_subprocess_files(self.hmm_cachefname, n_procs)
            
        if not self.args.no_clean:
            if n_procs == 1:
                os.remove(self.hmm_outfname)
            else:
                for iproc in range(n_procs):
                    subworkdir = self.args.workdir + '/hmm-' + str(iproc)
                    os.remove(subworkdir + '/' + os.path.basename(self.hmm_infname))
                    os.remove(subworkdir + '/' + os.path.basename(self.hmm_outfname))
                    os.rmdir(subworkdir)

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
        print '  writing hmms with info from %s' % parameter_dir
        # start = time.time()
        from hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')

        # gene_list = self.args.only_genes
        # if gene_list is None and self.sw_info is not None:  # if specific genes weren't specified, do the ones for which we have sw matches
        #     print 'only-gene s argument not specified, writing hmms using sw matches'
        #     gene_list = []
        #     for region in utils.regions:
        #         for gene in self.germline_seqs[region]:
        #             if gene in self.sw_info['all_best_matches']:
        #                 gene_list.append(gene)

        # if gene_list is None:  # ack, just do 'em all
        #     print 'just do them all'
        #     gene_list = []
        #     for region in utils.regions:
        #         gene_list += list(self.germline_seqs[region].keys())

        if self.args.only_genes is None:  # make a list of all the genes for which we have counts in <parameter_dir> (a.tm., this is all the genes that appeared as a best match at least once)
            gene_list = []
            for region in utils.regions:
                with opener('r')(parameter_dir + '/' + region + '_gene-probs.csv') as pfile:
                    reader = csv.DictReader(pfile)
                    for line in reader:
                        gene_list.append(line[region + '_gene'])
        else:
            gene_list = self.args.only_genes

        for gene in gene_list:
            if self.args.debug:
                print '  %s' % utils.color_gene(gene)
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety,
                               self.germline_seqs, self.args,
                               self.cyst_positions, self.tryp_positions)
            writer.write()

        # print '    time to write hmms: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def check_hmm_existence(self, gene_list, skipped_gene_matches, parameter_dir):
        """ Check if hmm model file exists, and if not remove gene from <gene_list> """
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
    def get_truncation_parameters(self, seqinfo, cyst_positions, extra_info=None, debug=False):
        # kinda obscurely written 'cause I generalized it to work with padding, then switched over to another function
        # find min length [over sequences in <seqinfo>] to left and right of the cysteine positions
        typelist = ['left', 'right']
        if extra_info is not None:
            typelist += ['v_5p_del', 'j_3p_del']
        lengths = {'min' : {t : None for t in typelist},
                   'max' : {t : None for t in typelist}}
        for query, seq in seqinfo.items():
            cpos = cyst_positions[query]
            if cpos < 0 or cpos >= len(seq):
                raise Exception('cpos %d invalid for %s (%s)' % (cpos, query, seq))
            thislen = {'left' : cpos, 'right' : len(seq) - cpos}  # NOTE right-hand one includes <cpos>, i.e. dleft + dright = len(seq)
            if extra_info is not None:
                thislen['v_5p_del'] = extra_info[query]['v_5p_del']
                thislen['j_3p_del'] = extra_info[query]['j_3p_del']
            for tp in typelist:
                if lengths['min'][tp] is None or delta[tp] < lengths['min'][tp]:
                    lengths['min'][tp] = delta[tp]
                if lengths['max'][tp] is None or delta[tp] > lengths['max'][tp]:
                    lengths['max'][tp] = delta[tp]
            if debug:
                print '        %d %d  (%d, %d - %d)   %s' % (dleft, dright, cpos, len(seq), cpos, query)

        return lengths

    # ----------------------------------------------------------------------------------------
    def get_padding_parameters(self, debug=False):
        all_v_matches, all_j_matches = set(), set()
        for query in self.sw_info['queries']:
            swfo = self.sw_info[query]
            for match in swfo['all'].split(':'):
                if utils.get_region(match) == 'v':
                    all_v_matches.add(match)
                elif utils.get_region(match) == 'j':
                    all_j_matches.add(match)

        maxima = {'gl_cpos' : None, 'gl_cpos_to_j_end' : None}  #, 'fv_insertion_len' : None, 'jf_insertion_len' : None}
        for query in self.sw_info['queries']:
            fvstuff = max(0, len(swfo['fv_insertion']) - swfo['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
            jfstuff = max(0, len(swfo['jf_insertion']) - swfo['j_3p_del'])

            for v_match in all_v_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
                gl_cpos = self.cyst_positions[v_match]['cysteine-position'] + fvstuff
                if maxima['gl_cpos'] is None or gl_cpos > maxima['gl_cpos']:
                    maxima['gl_cpos'] = gl_cpos

            swfo = self.sw_info[query]
            seq = swfo['seq']
            cpos = swfo['cyst_position']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
            for j_match in all_j_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
                # TODO this is totally wrong -- I'm only storing j_3p_del for the best match... but hopefully it'll give enough padding for the moment
                gl_cpos_to_j_end = len(seq) - cpos + swfo['j_3p_del'] + jfstuff
                if maxima['gl_cpos_to_j_end'] is None or gl_cpos_to_j_end > maxima['gl_cpos_to_j_end']:
                    maxima['gl_cpos_to_j_end'] = gl_cpos_to_j_end

            # if maxima['fv_insertion_len'] is None or len(swfo['fv_insertion']) > maxima['fv_insertion_len']:
            #     maxima['fv_insertion_len'] = len(swfo['fv_insertion'])
            # if maxima['jf_insertion_len'] is None or len(swfo['jf_insertion']) > maxima['jf_insertion_len']:
            #     maxima['jf_insertion_len'] = len(swfo['jf_insertion'])

        if debug:
            print '    maxima:',
            for k, v in maxima.items():
                print '%s %d    ' % (k, v),
            print ''
        return maxima

    # ----------------------------------------------------------------------------------------
    def pad_seqs_to_same_length(self, debug=False):
        """
        Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
        Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
        """

        maxima = self.get_padding_parameters(debug)

        for query in self.sw_info['queries']:
            swfo = self.sw_info[query]
            if 'padded' in swfo:  # already added padded information (we're probably partitioning, and this is not the first step)
                return
            seq = swfo['seq']
            cpos = swfo['cyst_position']
            if cpos < 0 or cpos >= len(seq):
                print 'hm now what do I want to do here?'
            k_v = swfo['k_v']

            # padleft = maxima['fv_insertion_len'] + maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
            # padright = maxima['gl_cpos_to_j_end'] + maxima['jf_insertion_len'] - (len(seq) - cpos)
            padleft = maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
            padright = maxima['gl_cpos_to_j_end'] - (len(seq) - cpos)
            if padleft < 0 or padright < 0:
                raise Exception('bad padding %d %d for %s' % (padleft, padright, query))

            swfo['padded'] = {}
            padfo = swfo['padded']  # shorthand
            assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
            padfo['seq'] = padleft * utils.ambiguous_bases[0] + seq + padright * utils.ambiguous_bases[0]
            if query in self.sw_info['indels']:
                print '    also padding reversed sequence'
                self.sw_info['indels'][query]['reversed_seq'] = padleft * utils.ambiguous_bases[0] + self.sw_info['indels'][query]['reversed_seq'] + padright * utils.ambiguous_bases[0]
            padfo['k_v'] = {'min' : k_v['min'] + padleft, 'max' : k_v['max'] + padleft}
            padfo['cyst_position'] = swfo['cyst_position'] + padleft
            padfo['padleft'] = padleft
            padfo['padright'] = padright
            if debug:
                print '      pad %d %d   %s' % (padleft, padright, query)
                print '     %d --> %d (%d-%d --> %d-%d)' % (len(seq), len(padfo['seq']),
                                                            k_v['min'], k_v['max'],
                                                            padfo['k_v']['min'], padfo['k_v']['max'])

        if debug:
            for query in self.sw_info['queries']:
                print '%20s %s' % (query, self.sw_info[query]['padded']['seq'])

    # ----------------------------------------------------------------------------------------
    def truncate_seqs(self, seqinfo, kvinfo, cyst_positions, debug=False):
        """ 
        Truncate <seqinfo> to have the same length to the left and right of the conserved cysteine.
        """

        lengths = self.get_truncation_parameters(seqinfo, cyst_positions, debug)

        # truncate all the sequences to these lengths
        for query, seq in seqinfo.items():
            cpos = cyst_positions[query]
            istart = cpos - lengths['min']['left']
            istop = cpos + lengths['min']['right']
            chopleft = istart
            chopright = len(seq) - istop
            if debug:
                print '      chop %d %d   %s' % (chopleft, chopright, query)
                print '     %d --> %d (%d-%d --> %d-%d)      %s' % (len(seq), len(seq[istart : istop]),
                                                                    -1 if kvinfo is None else kvinfo[query]['min'], -1 if kvinfo is None else kvinfo[query]['max'],
                                                                    -1 if kvinfo is None else (kvinfo[query]['min'] - chopleft), -1 if kvinfo is None else (kvinfo[query]['max'] - chopleft),
                                                                    query)
            seqinfo[query] = seq[istart : istop]
            if kvinfo is not None:
                kvinfo[query]['min'] -= chopleft
                kvinfo[query]['max'] -= chopleft

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, query_names, parameter_dir, skipped_gene_matches=None):
        """ 
        Return the 'logical OR' of the queries in <query_names>, i.e. the maximal extent in k_v/k_d space and OR of only_gene sets.
        Also truncates sequences.
        """

        combo = {
            'k_v':{'min':99999, 'max':-1},
            'k_d':{'min':99999, 'max':-1},
            'only_genes':[],
            'seqs':[],
            'mute-freqs':[],
            'cyst_positions':[]
        }

        # TODO this whole thing probably ought to use cached hmm info if it's available
        # TODO this just always uses the SW mutation rate, but I should really update it with the (multi-)hmm-derived ones (same goes for k space boundaries)
        # TODO/NOTE this is the mutation rate *before* truncation

        for name in query_names:
            swfo = self.sw_info[name]
            if 'padded' in swfo:
                assert not self.args.dont_pad_sequences
                k_v = swfo['padded']['k_v']
                seq = swfo['padded']['seq']
                cpos = swfo['padded']['cyst_position']
            else:
                k_v = swfo['k_v']
                seq = swfo['seq']
                cpos = swfo['cyst_position']
            k_d = swfo['k_d']  # don't need to adjust k_d for padding
            combo['seqs'].append(seq)
            combo['mute-freqs'].append(utils.get_mutation_rate(self.germline_seqs, swfo))
            combo['cyst_positions'].append(cpos)  # TODO use cached hmm values instead of SW
            combo['k_v']['min'] = min(k_v['min'], combo['k_v']['min'])
            combo['k_v']['max'] = max(k_v['max'], combo['k_v']['max'])
            combo['k_d']['min'] = min(k_d['min'], combo['k_d']['min'])
            combo['k_d']['max'] = max(k_d['max'], combo['k_d']['max'])

            # work out which genes to tell the hmm to use
            only_genes = swfo['all'].split(':')  # start with all the sw matches for this query
            self.check_hmm_existence(only_genes, skipped_gene_matches, parameter_dir)  # remove the ones for which we don't have hmm files (we only write hmms for genes that appeared as the best sw match for at least one query, but swfo['all'] in general includes genes that were never the *best* match for any one query)
            genes_to_use = []
            for region in utils.regions:  # take the best <self.args.n_max_per_region> from each region
                reg_genes = [g for g in only_genes if utils.get_region(g) == region]
                n_genes = min(len(reg_genes), self.args.n_max_per_region[utils.regions.index(region)])  # minimum of [the number of gene matches for this region] and [the number we want for this region]
                for ig in range(n_genes):  # take the first <n_genes> matches (they're ordered by sw match score)
                    genes_to_use.append(reg_genes[ig])

            # and finally OR this query's genes into the ones from previous queries
            combo['only_genes'] = list(set(genes_to_use) | set(combo['only_genes']))  # NOTE using the OR of all sets of genes (from all query seqs) like this *really* helps,

        if not self.all_regions_present(combo['only_genes'], skipped_gene_matches, query_names):
            return {}

        return combo

    # ----------------------------------------------------------------------------------------
    def remove_sw_failures(self, query_names):
        """ If any of the queries in <query_names> was unproductive, return an empty list (which will be skipped entirely), otherwise return the original name list """
        unproductive, unknown = False, False
        for qrn in query_names:
            if qrn in self.sw_info['skipped_unproductive_queries']:
                unproductive = True
            if qrn in self.sw_info['skipped_unknown_queries']:
                unknown = True
        if unproductive or unknown:
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
    def write_to_single_input_file(self, fname, mode, nsets, parameter_dir, skipped_gene_matches, path_index=0, logweight=0.):
        csvfile = opener(mode)(fname)
        header = ['path_index', 'logweight', 'names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seqs', 'mute_freqs', 'cyst_positions']  # NOTE logweight is for the whole partition
        writer = csv.DictWriter(csvfile, header, delimiter=' ')  # NOTE should eventually rewrite arg parser in ham to handle csvs (like in glomerator cache reader)
        if mode == 'w':
            writer.writeheader()
        # start = time.time()

        # for k in nsets:
        #     print k
        if self.args.random_divvy:  #randomize_input_order:  # NOTE nsets is a list of *lists* of ids
            print 'randomizing input order'
            random_nsets = []
            while len(nsets) > 0:
                irand = random.randint(0, len(nsets) - 1)  # NOTE interval is inclusive
                random_nsets.append(nsets[irand])
                nsets.remove(nsets[irand])
            nsets = random_nsets
        # print '---'
        # for k in nsets:
        #     print k
        
        for query_names in nsets:
            non_failed_names = self.remove_sw_failures(query_names)
            if len(non_failed_names) == 0:
                continue
            combined_query = self.combine_queries(non_failed_names, parameter_dir, skipped_gene_matches=skipped_gene_matches)
            if len(combined_query) == 0:  # didn't find all regions
                continue
            writer.writerow({
                'path_index' : path_index,
                'logweight' : logweight,  # NOTE same for all lines with the same <path_index> (since they're all from the same partition)
                'names' : ':'.join([qn for qn in non_failed_names]),
                'k_v_min' : combined_query['k_v']['min'],
                'k_v_max' : combined_query['k_v']['max'],
                'k_d_min' : combined_query['k_d']['min'],
                'k_d_max' : combined_query['k_d']['max'],
                'only_genes' : ':'.join(combined_query['only_genes']),
                'seqs' : ':'.join(combined_query['seqs']),  # may be truncated, and thus not the same as those in <input_info> or <sw_info>
                'mute_freqs' : ':'.join([str(f) for f in combined_query['mute-freqs']]),  # a.t.m., not corrected for truncation
                'cyst_positions' : ':'.join([str(cpos) for cpos in combined_query['cyst_positions']]),  # TODO should really use the hmm cpos if it's available
                # 'cyst_positions' : ':'.join([str(self.sw_info[qn]['cyst_position']) for qn in non_failed_names])  # TODO should really use the hmm cpos if it's available
            })

        csvfile.close()
        # print '        input write time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, parameter_dir):
        """ Write input file for bcrham """
        print '    writing input'

        # if (self.args.action == 'partition' or self.args.n_sets > 1) and not self.args.dont_pad_sequences:
        if not self.args.dont_pad_sequences:
            self.pad_seqs_to_same_length()  # adds padded info to sw_info (returns if stuff has already been padded)

        skipped_gene_matches = set()

        if self.args.smc_particles > 1:
            assert self.args.action == 'partition'
            n_procs = len(self.smc_info[-1])
            for iproc in range(n_procs):
                if n_procs == 1:
                    fname = self.hmm_infname
                else:
                    subworkdir = self.args.workdir + '/hmm-' + str(iproc)
                    utils.prep_dir(subworkdir)
                    if os.path.exists(self.hmm_cachefname):  # copy cachefile to this subdir
                        check_call(['cp', self.hmm_cachefname, subworkdir + '/'])  # NOTE this is kind of wasteful to write it to each subdirectory (it could be large) but it's cleaner this way, 'cause then the subdirs are independent
                    fname = subworkdir + '/' + os.path.basename(self.hmm_infname)
                procinfo = self.smc_info[-1][iproc]  # list of ClusterPaths, one for each smc particle
                for iptl in range(len(procinfo)):
                    path = procinfo[iptl]
                    self.write_to_single_input_file(fname, 'w' if iptl==0 else 'a', list(path.partitions[path.i_best_minus_x]), parameter_dir,  #  list() is important since we may modify <nsets>
                                                    skipped_gene_matches, path_index=iptl, logweight=path.logweights[path.i_best_minus_x])
        else:
            if self.args.action == 'partition':
                nsets = list(self.paths[-1].partitions[self.paths[-1].i_best_minus_x])  #  list() is important since we modify <nsets>
            else:
                if self.args.n_sets == 1:  # single vanilla hmm (does the same thing as the below for n=1, but is more transparent)
                    nsets = [[qn] for qn in self.input_info.keys()]
                else:
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

            self.write_to_single_input_file(self.hmm_infname, 'w', nsets, parameter_dir, skipped_gene_matches)

        if self.args.debug and len(skipped_gene_matches) > 0:
            print '    not found in %s, so removing from consideration for hmm (i.e. were only the nth best, but never the best sw match for any query):' % (parameter_dir),
            for region in utils.regions:
                # print '  %s: %d' % (region, len([gene for gene in skipped_gene_matches if utils.get_region(gene) == region])),
                print '\n      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in sorted(skipped_gene_matches) if utils.get_region(gene) == region]))
            print ''

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, algorithm, n_procs, count_parameters, parameter_out_dir, cache_naive_seqs):
        if self.args.smc_particles == 1:
            if self.args.action == 'partition' or n_procs > 1:
                self.merge_all_hmm_outputs(n_procs, cache_naive_seqs)
        else:
            self.merge_pairs_of_procs(n_procs)

        # if os.path.exists(self.hmm_cachefname):
        #     self.read_cachefile(self.hmm_cachefname)

        if self.args.action != 'partition':
            self.read_annotation_output(algorithm, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir)

        if not self.args.no_clean and os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)

    # ----------------------------------------------------------------------------------------
    def get_bcrham_truncations(self, line):
        chops = []
        for iq in range(len(line['unique_ids'])):
            # print line['unique_ids'][iq]
            original_seq = self.input_info[line['unique_ids'][iq]]['seq']
            bcrham_seq = line['seqs'][iq]
            chopleft = original_seq.find(bcrham_seq)
            assert chopleft >= 0 and chopleft <= len(original_seq)
            chopright = len(original_seq) - len(bcrham_seq) - chopleft
            assert chopright >= 0 and chopright <= len(original_seq)

            chops.append({'left' : chopleft, 'right' : chopright})

            # print '    ', line['v_5p_del'], line['j_3p_del']
            # TODO this isn't right yet

            truncated_seq = original_seq
            if chopleft > 0:
                truncated_seq = truncated_seq[chopleft : ]
                # line['chopleft'] = chopleft
                # if line['v_5p_del'] < chopleft:
                #     print 'ERROR v_5p_del %d smaller than chopleft %d, i.e. bcrham probably couldn\'t delete some things it wanted to' % (line['v_5p_del'], chopleft)
                # else:
                #     line['v_5p_del'] -= chopleft
                # line['seqs'][iq] = original_seq
            if chopright > 0:
                truncated_seq = truncated_seq[ : -chopright]
                # line['chopright'] = chopright
                # if line['j_3p_del'] < chopright:
                #     print 'ERROR j_3p_del %d smaller than chopright %d, i.e. bcrham probably couldn\'t delete some things it wanted to' % (line['j_3p_del'], chopright)
                # else:
                #     line['j_3p_del'] -= chopright
                # line['seqs'][iq] = original_seq
            # print '    ', line['v_5p_del'], line['j_3p_del']

            # print '    ', original_seq
            # print '    ', bcrham_seq
            # print '    ', truncated_seq
            assert bcrham_seq == truncated_seq

        return chops

    # ----------------------------------------------------------------------------------------
    def read_naive_hamming_clusters(self, n_procs):
        if n_procs > 1:
            self.merge_subprocess_files(self.hmm_outfname, n_procs)
        self.bcrham_divvied_queries = []
        with opener('r')(self.hmm_outfname) as hmm_csv_outfile:
            lines = [l.strip() for l in hmm_csv_outfile.readlines()]
            if len(lines) != n_procs + 1:
                raise Exception('expected %d lines (%d procs), but got %d!' % (n_procs + 1, n_procs, len(lines)))
            if lines[0] != 'partition':
                raise Exception('bad bcrham naive clustering output header: %s' % lines[0])
            for line in lines[1:]:
                clusters = line.split('|')
                self.bcrham_divvied_queries += [cl.split(';') for cl in clusters]

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)

        print '  read divvies from bcrham:'
        print '      ', ' '.join([str(len(cl)) for cl in self.bcrham_divvied_queries])

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, algorithm, count_parameters=False, parameter_out_dir=None):
        """ Read bcrham annotation output """
        print '    read output'

        if count_parameters:
            assert parameter_out_dir is not None
        pcounter = ParameterCounter(self.germline_seqs) if count_parameters else None
        true_pcounter = ParameterCounter(self.germline_seqs) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter(self.germline_seqs, 'hmm') if self.args.plot_performance else None

        n_seqs_processed, n_events_processed = 0, 0
        hmminfo = {}
        with opener('r')(self.hmm_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            boundary_error_queries = []
            for line in reader:
                utils.process_input_line(line,
                                         splitargs=('unique_ids', 'seqs'),
                                         int_columns=('nth_best', 'v_5p_del', 'd_5p_del', 'cdr3_length', 'j_5p_del', 'j_3p_del', 'd_3p_del', 'v_3p_del'),
                                         float_columns=('logprob'))
                ids = line['unique_ids']
                same_event = utils.from_same_event(self.args.is_data, self.reco_info, ids)
                if same_event is None:
                    same_event = -1

                # check for errors
                if line['nth_best'] == 0:  # if this is the first line for this set of ids (i.e. the best viterbi path or only forward score)
                    if line['errors'] is not None and 'boundary' in line['errors'].split(':'):
                        boundary_error_queries.append(':'.join([uid for uid in ids]))
                    else:
                        assert len(line['errors']) == 0

                utils.add_cdr3_info(self.germline_seqs, self.cyst_positions, self.tryp_positions, line)
                if self.args.debug:
                    if line['nth_best'] == 0:  # if this is the first line (i.e. the best viterbi path) for this query (or query pair), print the true event
                        print '      %s' % ':'.join(ids),
                        if not self.args.is_data:
                            print '   %d' % same_event,
                        print ''
                    self.print_hmm_output(line, print_true=(line['nth_best']==0))  #, perfplotter=perfplotter)
                if line['nth_best'] == 0 and (line['cdr3_length'] != -1 or not self.args.skip_unproductive):  # if it's productive, or if we're not skipping unproductive rearrangements
                    if pcounter is not None:
                        pcounter.increment_reco_params(line)
                    if true_pcounter is not None:
                        true_pcounter.increment_reco_params(self.reco_info[ids[0]])  # NOTE doesn't matter which id you pass it, since they all have the same reco parameters
                    n_events_processed += 1
                    for iseq in range(len(ids)):
                        uid = ids[iseq]
                        hmminfo[uid] = dict(line)  # make a copy of the info, into which we'll insert the sequence-specific stuff
                        hmminfo[uid]['seq'] = line['seqs'][iseq]
                        hmminfo[uid]['unique_id'] = uid
                        utils.add_match_info(self.germline_seqs, hmminfo[uid], self.cyst_positions, self.tryp_positions, debug=(self.args.debug > 0))
                        if pcounter is not None:
                            pcounter.increment_mutation_params(hmminfo[uid])
                        if true_pcounter is not None:
                            true_pcounter.increment_mutation_params(self.reco_info[uid])  # NOTE doesn't matter which id you pass it, since they all have the same reco parameters
                        if perfplotter is not None:
                            perfplotter.evaluate(self.reco_info[uid], hmminfo[uid], None if self.args.dont_pad_sequences else self.sw_info[uid]['padded'])
                        n_seqs_processed += 1

        if pcounter is not None:
            pcounter.write(parameter_out_dir)
            if self.args.plotdir is not None:
                pcounter.plot(self.args.plotdir + '/hmm', subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)
        if true_pcounter is not None:
            true_pcounter.write(parameter_out_dir + '/true')
            if self.args.plotdir is not None:
                true_pcounter.plot(self.args.plotdir + '/hmm/true', subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)
        if perfplotter is not None:
            assert self.args.plotdir is not None
            perfplotter.plot(self.args.plotdir + '/hmm/performance')

        print '    processed %d sequences (%d events)' % (n_seqs_processed, n_events_processed)
        if len(boundary_error_queries) > 0:
            print '      %d boundary errors (%s)' % (len(boundary_error_queries), ', '.join(boundary_error_queries))

        if self.args.outfname is not None:
            outpath = self.args.outfname
            if self.args.outfname[0] != '/':  # if full output path wasn't specified on the command line
                outpath = os.getcwd() + '/' + outpath
            shutil.copyfile(self.hmm_outfname, outpath)
            with open(outpath) as outfile:
                reader = csv.DictReader(outfile)
                outfo = []
                for line in reader:
                    outfo.append(line)
                    outfo[-1]['naive_seq'] = utils.get_full_naive_seq(self.germline_seqs, line)
            with open(outpath, 'w') as outfile:
                writer = csv.DictWriter(outfile, outfo[0].keys())
                writer.writeheader()
                for line in outfo:
                    writer.writerow(line)

        if self.args.annotation_clustering == 'vollmers':
            if self.args.outfname is not None:
                outfile = open(self.args.outfname, 'w')  # NOTE overwrites annotation info that's already been written to <self.args.outfname>
                headers = ['n_clusters', 'threshold', 'clusters']  #, 'true_clusters']
                if not self.args.is_data:
                    headers += ['adj_mi', ]  #, 'n_true_clusters']
                writer = csv.DictWriter(outfile, headers)
                writer.writeheader()

            for thresh in self.args.annotation_clustering_thresholds:
                adj_mi, partition = annotationclustering.vollmers(hmminfo, threshold=thresh, reco_info=self.reco_info)
                n_clusters = len(partition)
                if self.args.outfname is not None:
                    row = {'n_clusters' : n_clusters, 'threshold' : thresh, 'clusters' : utils.get_str_from_partition(partition)}
                    if not self.args.is_data:
                        row['adj_mi'] = adj_mi
                        # row['n_true_clusters'] = len(utils.get_true_partition(self.reco_info))
                        # true_partition = [cl for cl in utils.get_true_partition(self.reco_info).values()]
                        # row['true_clusters'] = utils.get_str_from_partition(true_partition)
                    writer.writerow(row)
            if self.args.outfname is not None:
                outfile.close()

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)

    # ----------------------------------------------------------------------------------------
    def print_hmm_output(self, line, print_true=False):
        out_str_list = []
        if print_true and not self.args.is_data:  # first print true event (if this is simulation)
            for uids in utils.get_true_clusters(line['unique_ids'], self.reco_info).values():
                synthetic_true_line = dict(self.reco_info[uids[0]])
                synthetic_true_line['unique_ids'] = uids
                synthetic_true_line['seqs'] = [self.reco_info[iid]['seq'] for iid in uids]
                del synthetic_true_line['unique_id']
                del synthetic_true_line['seq']
                # del synthetic_true_line['indels']
                indelfos = [self.reco_info[iid]['indels'] for iid in uids]
                event_str = utils.print_reco_event(self.germline_seqs, synthetic_true_line, extra_str='    ', return_string=True, label='true:', indelfos=indelfos)
                out_str_list.append(event_str)

        event_str = utils.print_reco_event(self.germline_seqs, line, extra_str='    ', return_string=True, label='inferred:', indelfos=[self.sw_info['indels'].get(uid, None) for uid in line['unique_ids']])
        out_str_list.append(event_str)

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
