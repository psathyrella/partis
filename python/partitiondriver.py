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
from clusterpath import ClusterPath
from waterer import Waterer
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter

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

        randomize_order = self.args.action == 'partition' and not self.args.force_dont_randomize_input_order
        if self.args.seqfile is not None:
            self.input_info, self.reco_info = get_seqfile_info(self.args.seqfile, self.args.is_data, self.germline_seqs, self.cyst_positions, self.tryp_positions,
                                                               self.args.n_max_queries, self.args.queries, self.args.reco_ids, randomize_order=randomize_order, replace_N_with=None)  # if we're partitioning, we need to randomize input order (at least for simulation)

        self.cached_results = None

        self.sw_info = None

        utils.prep_dir(self.args.workdir)
        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'

        if self.args.outfname is not None:
            outdir = os.path.dirname(self.args.outfname)
            if outdir != '' and not os.path.exists(outdir):
                os.makedirs(outdir)

    # ----------------------------------------------------------------------------------------
    def merge_cachefiles(self, infnames, outfname):
        for fname in infnames:
            self.read_cachefile(fname)
        self.write_cachefile(outfname)

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.args.initial_cachefname is not None:
            self.merge_cachefiles(infnames=[self.hmm_cachefname, self.args.initial_cachefname], outfname=self.args.initial_cachefname)
            # check_call(['cp', '-v', self.hmm_cachefname, self.args.initial_cachefname])
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
        # assert self.args.n_sets == 1  # er, could do it for n > 1, but I'd want to think through a few things first
        assert self.args.plotdir is not None

        sw_parameter_dir = self.args.parameter_dir + '/sw'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=self.args.plotdir + '/sw')
        waterer.run()
        self.sw_info = waterer.info
        self.write_hmms(sw_parameter_dir)
        parameter_out_dir = self.args.parameter_dir + '/hmm'
        self.run_hmm('viterbi', parameter_in_dir=sw_parameter_dir, parameter_out_dir=parameter_out_dir, count_parameters=True, plotdir=self.args.plotdir + '/hmm')
        self.write_hmms(parameter_out_dir)

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self, algorithm):
        """ Just run <algorithm> (either 'forward' or 'viterbi') on sequences in <self.input_info> and exit. You've got to already have parameters cached in <self.args.parameter_dir> """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir (' + self.args.parameter_dir + ') d.n.e')
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()

        self.sw_info = waterer.info
        self.run_hmm(algorithm, parameter_in_dir=self.args.parameter_dir, count_parameters=self.args.plot_parameters, plotdir=self.args.plotdir)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        if not os.path.exists(self.args.parameter_dir):
            raise Exception('parameter dir %s d.n.e.' % self.args.parameter_dir)

        # run smith-waterman
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs, parameter_dir=self.args.parameter_dir, write_parameters=False)
        waterer.run()
        self.sw_info = waterer.info
        if self.args.pad_sequences:  # have to do this before we divvy... sigh TODO clean this up and only call it in one place
            # TODO oh, duh, just do it in waterer
            self.pad_seqs_to_same_length()  # adds padded info to sw_info

        n_procs = self.args.n_procs
        n_proc_list = []  # list of the number of procs we used for each run

        # add initial lists of paths
        if self.args.smc_particles == 1:
            cp = ClusterPath(-1)
            cp.add_partition([[cl, ] for cl in self.input_info.keys()], 0., 0., -1.)
            self.paths = [cp, ]
        else:
            initial_divvied_queries = self.divvy_up_queries(n_procs, self.sw_info)
            self.smc_info = [[], ]
            for clusters in initial_divvied_queries:  # one set of <clusters> for each process
                self.smc_info[-1].append([])
                for iptl in range(self.args.smc_particles):
                    cp = ClusterPath(-1)
                    cp.add_partition([[cl, ] for cl in clusters], 0., 0., -1.)
                    self.smc_info[-1][-1].append(cp)

        # get number of clusters based on sum of last paths in <self.smc_info>
        def get_n_clusters():
            if self.args.smc_particles == 1:
                return len(self.paths[-1].partitions[self.paths[-1].i_best_minus_x])

            nclusters = 0
            for iproc in range(len(self.smc_info[-1])):  # number of processes
                path = self.smc_info[-1][iproc][0]  # uses the first smc particle, but the others will be similar
                nclusters += len(path.partitions[path.i_best_minus_x])
            return nclusters

        # run that shiznit
        while n_procs > 0:
            nclusters = get_n_clusters()
            print '--> %d clusters with %d procs' % (nclusters, n_procs)  # write_hmm_input uses the best-minus-ten partition
            self.run_hmm('forward', self.args.parameter_dir, n_procs=n_procs, shuffle_input_order=(self.args.smc_particles == 1))  # don't shuffle sequences if we have multiple paths, 'cause, you know, if you do ALL HELL BREAKS LOOSE
            n_proc_list.append(n_procs)

            if n_procs == 1:
                break

            if self.args.smc_particles == 1:  # for smc, we merge pairs of processes; otherwise, we do some heuristics to come up with a good number of clusters for the next iteration
                n_procs = n_procs / 2
            else:
                n_procs = len(self.smc_info[-1])  # if we're doing smc, the number of particles is determined by the file merging process

        # deal with final partition
        tmpglom = Glomerator(self.reco_info)
        if self.args.smc_particles == 1:
            for path in self.paths:
                self.check_path(path)
            if self.args.debug:
                print 'final'
                for ipath in range(len(self.paths)):  # one path for each glomeration step
                    self.paths[ipath].print_partitions(self.reco_info, one_line=True, header=(ipath==0))
                    print ''
            if self.args.outfname is not None:
                self.write_partitions(self.args.outfname, [self.paths[-1], ], tmpglom)  # [last agglomeration step]
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
                self.write_partitions(self.args.outfname, final_paths, tmpglom)

        if self.args.debug and not self.args.is_data:
            tmpglom.print_true_partition()

    # ----------------------------------------------------------------------------------------
    def check_path(self, path):
        def check_partition(partition, ipart):
            for uid in self.input_info:
                found = False
                for cluster in partition:
                    if uid in cluster:
                        found = True
                        break
                if not found and uid not in self.sw_info['skipped_unproductive_queries']:
                    path.print_partition(ipart, self.reco_info, one_line=False, abbreviate=False)
                    raise Exception('%s not found in merged partition' % uid)
            for cluster in partition:
                for uid in cluster:
                    if uid not in self.input_info:
                        raise Exception('%s not found in merged partition' % uid)

        for ipart in range(len(path.partitions)):
            check_partition(path.partitions[ipart], ipart)

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, outfname, paths, tmpglom):
        with opener('w')(outfname) as outfile:
            writer = csv.DictWriter(outfile, ('path_index', 'score', 'logweight', 'adj_mi', 'n_clusters', 'n_procs', 'bad_clusters'))  #'normalized_score'
            writer.writeheader()
            true_partition = tmpglom.get_true_partition()
            for ipath in range(len(paths)):
                for ipart in range(len(paths[ipath].partitions)):
                    # print 'ipart', ipart
                    part = paths[ipath].partitions[ipart]
                    cluster_str = ''
                    bad_clusters = []  # inferred clusters that aren't really all from the same event
                    for ic in range(len(part)):
                        # print '  %d:  %s' % (ic, ':'.join(part[ic]))
                        if ic > 0:
                            cluster_str += ';'
                        cluster_str += ':'.join(part[ic])

                        same_event = utils.from_same_event(self.args.is_data, self.reco_info, part[ic])  # are all the sequences from the same event?
                        entire_cluster = True  # ... and if so, are they the entire true cluster?
                        if same_event:
                            reco_id = self.reco_info[part[ic][0]]['reco_id']  # they've all got the same reco_id then, so pick an aribtrary one
                            true_cluster = true_partition[reco_id]
                            # print '    true %s:  %s' % (reco_id, ':'.join(true_cluster))
                            for uid in true_cluster:
                                if uid not in part[ic]:
                                    # print '        missing %s (and maybe others)' % uid
                                    entire_cluster = False
                                    break
                        else:
                            entire_cluster = False
                        if not same_event or not entire_cluster:
                            bad_clusters.append(':'.join(part[ic]))

                    if len(bad_clusters) > 25:
                        bad_clusters = ['too', 'long']
                    writer.writerow({'path_index' : os.getpid() + ipath,
                                     'score' : paths[ipath].logprobs[ipart],
                                     'logweight' : paths[ipath].logweights[ipart],
                                     # 'normalized_score' : part['score'] / self.max_log_probs[ipath],
                                     'adj_mi' : paths[ipath].adj_mis[ipart],
                                     'n_clusters' : len(part),
                                     'n_procs' : paths[ipath].n_procs[ipart],
                                     'bad_clusters' : ';'.join(bad_clusters)
                                     # 'clusters' : cluster_str
                                 })

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
        if self.args.truncate_n_sets:
            cmd_str += ' --truncate-seqs'
        if self.args.allow_unphysical_insertions:
            cmd_str += ' --unphysical-insertions'
        assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
        cmd_str += ' --ambig-base ' + utils.ambiguous_bases[0]

        # print cmd_str
        # sys.exit()
        return cmd_str

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, parameter_in_dir, parameter_out_dir='', count_parameters=False, plotdir=None, n_procs=None, shuffle_input_order=False):
        """ 
        Run bcrham, possibly with many processes, and parse and interpret the output.
        NOTE the local <n_procs>, which overrides the one from <self.args>
        """
        print 'hmm'
        if n_procs is None:
            n_procs = self.args.n_procs

        self.write_hmm_input(parameter_dir=parameter_in_dir, shuffle_input_order=shuffle_input_order)

        print '    running'
        sys.stdout.flush()
        # start = time.time()
        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir)
        if n_procs == 1:
            check_call(cmd_str.split())
        else:
            if self.args.smc_particles == 1:  # if we're doing smc (i.e. if > 1), we have to split things up more complicatedly elsewhere
                self.split_input(n_procs, infname=self.hmm_infname, prefix='hmm')

            procs = []
            for iproc in range(n_procs):
                workdir = self.args.workdir + '/hmm-' + str(iproc)
                proc = Popen(cmd_str.replace(self.args.workdir, workdir).split(), stdout=PIPE, stderr=PIPE)
                procs.append(proc)
                time.sleep(0.1)
            for iproc in range(len(procs)):
                out, err = procs[iproc].communicate()
                utils.process_out_err(out, err, extra_str=str(iproc))

        sys.stdout.flush()
        # print '      hmm run time: %.3f' % (time.time()-start)

        self.read_hmm_output(algorithm, n_procs, count_parameters, parameter_out_dir, plotdir)

        if self.args.vollmers_clustering:
            vollmers_clusterer = Clusterer()
            vollmers_clusterer.vollmers_cluster(hmminfo, reco_info=self.reco_info)

    # ----------------------------------------------------------------------------------------
    def divvy_up_queries(self, n_procs, info, debug=True):
        naive_seqs = {}
        # , cyst_positions = {}, {}
        for line in info:
            query = line['names'] if 'names' in line else line['unique_id']  # the first time through, we just pass in <self.sw_info>, so we need 'unique_id' instead of 'names'
            seqstr = line['padded']['seqs'] if 'padded' in line else line['seqs']
            # NOTE cached naive seqs should all be the same length
            if self.cached_results is not None and seqstr in self.cached_results:  # first try to used cached hmm results
                naive_seqs[query] = self.cached_results[seqstr]['naive_seq']
                # cyst_positions[query] = self.cached_results[seqstr]['cyst_position']
            elif query in self.sw_info:  # ...but if we don't have them, use smith-waterman (should only be for single queries)
                naive_seqs[query] = utils.get_full_naive_seq(self.germline_seqs, self.sw_info[query])
                if self.args.pad_sequences:
                    padleft = self.sw_info[query]['padded']['padleft']  # we're padding the *naive* seq corresponding to query now, but it'll be the same length as the query seq
                    padright = self.sw_info[query]['padded']['padright']
                    assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
                    naive_seqs[query] = padleft * utils.ambiguous_bases[0] + naive_seqs[query] + padright * utils.ambiguous_bases[0]
                # cyst_positions[query] = self.sw_info[query]['cyst_position']
            elif len(query.split(':')) > 1:  # hmm... multiple queries without cached hmm naive sequences... that shouldn't happen
                print 'ERROR no naive sequence for %s' % query
                for qry in query.split(':'):
                    print '    %s %s' % (qry, utils.get_full_naive_seq(self.germline_seqs, self.sw_info[qry]))
            else:
                raise Exception('no naive sequence found for ' + str(query))
            if naive_seqs[query] == '':
                raise Exception('zero-length naive sequence found for ' + str(query))

        if self.args.truncate_n_sets:
            assert False  # deprecated and broken
            # print '  truncate in divvy'
            # self.truncate_seqs(naive_seqs, kvinfo=None, cyst_positions=cyst_positions)

        # print 'naive'
        # for qr, s in naive_seqs.items():
        #     print qr, s
        # sys.exit()

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
    def split_input(self, n_procs, infname, prefix):
        """ Do stuff. Probably correctly. """
        # read single input file
        assert self.args.smc_particles == 1
        info = []
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile, delimiter=' ')
            for line in reader:
                info.append(line)

        # initialize
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

        if self.args.action == 'partition':
            divvied_queries = self.divvy_up_queries(n_procs, info)
        for iproc in range(n_procs):
            for iquery in range(len(info)):
                if self.args.action == 'partition':
                    if info[iquery]['names'] not in divvied_queries[iproc]:  # NOTE I think the reason this doesn't seem to be speeding things up is that our hierarhical agglomeration time is dominated by the distance calculation, and that distance calculation time is roughly proportional to the number of sequences in the cluster (i.e. larger clusters take longer)
                        continue
                else:
                    if iquery % n_procs != iproc:
                        continue
                writers[iproc].writerow(info[iquery])

        for iproc in range(n_procs):
            sub_outfiles[iproc].close()

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
    def merge_all_hmm_outputs(self, n_procs):
        """ Merge any/all output files from subsidiary bcrham processes (used when *not* doing smc) """
        assert self.args.smc_particles == 1  # have to do things more complicatedly for smc
        if self.args.action == 'partition':  # merge partitions from several files
            if n_procs == 1:
                infnames = [self.hmm_outfname, ]
            else:
                infnames = [self.args.workdir + '/hmm-' + str(iproc) + '/' + os.path.basename(self.hmm_outfname) for iproc in range(n_procs)]
            previous_info = None
            if len(self.paths) > 1:
                previous_info = self.paths[-1]
            glomerer = Glomerator(self.reco_info)
            glomerer.read_cached_agglomeration(infnames, smc_particles=1, previous_info=previous_info, debug=self.args.debug)  #, outfname=self.hmm_outfname)
            assert len(glomerer.paths) == 1
            self.check_path(glomerer.paths[0])
            self.paths.append(glomerer.paths[0])

            if n_procs > 1:
                self.merge_csv_files(self.hmm_cachefname, n_procs)
        else:
            self.merge_csv_files(self.hmm_outfname, n_procs)

        if not self.args.no_clean:
            if n_procs == 1:
                os.remove(self.hmm_outfname)
            else:
                for iproc in range(n_procs):
                    subworkdir = self.args.workdir + '/hmm-' + str(iproc)
                    os.remove(subworkdir + '/' + os.path.basename(self.hmm_infname))
                    if os.path.exists(subworkdir + '/' + os.path.basename(self.hmm_outfname)):
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
            paths = glomerer.read_cached_agglomeration(infnames, self.args.smc_particles, previous_info=previous_info, debug=self.args.debug)  #, outfname=self.hmm_outfname)
            self.smc_info[-1].append(paths)

            # ack? self.glomclusters.append(glomerer)
            # boof? self.list_of_preclusters.append(glomerer.combined_conservative_best_minus_ten_partitions)

        if n_procs > 1:
            self.merge_csv_files(self.hmm_cachefname, n_procs)
            
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

        maxima = {'gl_cpos' : None, 'gl_cpos_to_j_end' : None}
        for query in self.sw_info['queries']:
            swfo = self.sw_info[query]
            seq = swfo['seq']
            cpos = swfo['cyst_position']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
            if cpos < 0 or cpos >= len(seq):
                raise Exception('cpos %d invalid for %s (%s)' % (cpos, query, seq))
            for v_match in all_v_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
                gl_cpos = self.cyst_positions[v_match]['cysteine-position']
                if maxima['gl_cpos'] is None or gl_cpos > maxima['gl_cpos']:
                    maxima['gl_cpos'] = gl_cpos
            for j_match in all_j_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
                gl_cpos_to_j_end = len(seq) - cpos + swfo['j_3p_del']  # TODO this is totally wrong -- I'm only storing j_3p_del for the best match... but hopefully it'll give enough padding for the moment
                if maxima['gl_cpos_to_j_end'] is None or gl_cpos_to_j_end > maxima['gl_cpos_to_j_end']:
                    maxima['gl_cpos_to_j_end'] = gl_cpos_to_j_end
            # if debug:
            #     print '        %d %d  (%d, %d - %d)   %s' % (dleft, dright, cpos, len(seq), cpos, query)

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

            padleft = maxima['gl_cpos'] - cpos
            padright = maxima['gl_cpos_to_j_end'] - (len(seq) - cpos)

            swfo['padded'] = {}
            padfo = swfo['padded']  # shorthand
            assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
            padfo['seq'] = padleft * utils.ambiguous_bases[0]+ seq + padright * utils.ambiguous_bases[0]
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
                print '%12s %s' % (query, self.sw_info[query]['padded']['seq'])

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
            if self.args.pad_sequences:
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
        unproductive, indel = False, False
        for qrn in query_names:
            if qrn in self.sw_info['skipped_unproductive_queries']:
                unproductive = True
            if qrn in self.sw_info['skipped_indel_queries']:
                indel = True
        if unproductive or indel:
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
    def write_hmm_input(self, parameter_dir, shuffle_input_order=False):
        """ Write input file for bcrham """
        print '    writing input'
        if self.cached_results is None:
            if self.args.initial_cachefname is not None:
                check_call(['cp', '-v', self.args.initial_cachefname, self.args.workdir + '/'])
        else:
            self.write_cachefile(self.hmm_cachefname)

        if self.args.pad_sequences:
            self.pad_seqs_to_same_length()  # adds padded info to sw_info

        def shuffle_nset_order(tmp_nsets):
            # randomize the order of the query list in <tmp_nsets>. Note that the list gets split into chunks for parallelization later
            assert self.args.smc_particles == 1
            random_nsets = []
            while len(tmp_nsets) > 0:
                irand = random.randint(0, len(tmp_nsets) - 1)  # NOTE interval is inclusive
                random_nsets.append(tmp_nsets[irand])
                tmp_nsets.remove(tmp_nsets[irand])
            return random_nsets

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

            if shuffle_input_order:  # TODO make sure this is ok, and doesn't overwrite anything untoward (<self.paths>)
                nsets = shuffle_nset_order(nsets)

            self.write_to_single_input_file(self.hmm_infname, 'w', nsets, parameter_dir, skipped_gene_matches)

        if len(skipped_gene_matches) > 0:
            print '    not found in %s, so removing from consideration for hmm (i.e. were only the nth best, but never the best sw match for any query):' % (parameter_dir),
            for region in utils.regions:
                # print '  %s: %d' % (region, len([gene for gene in skipped_gene_matches if utils.get_region(gene) == region])),
                print '\n      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in sorted(skipped_gene_matches) if utils.get_region(gene) == region]))
            print ''

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, algorithm, n_procs, count_parameters, parameter_out_dir, plotdir):
        if self.args.smc_particles == 1:
            if self.args.action == 'partition' or n_procs > 1:
                self.merge_all_hmm_outputs(n_procs)
        else:
            self.merge_pairs_of_procs(n_procs)

        if os.path.exists(self.hmm_cachefname):
            self.read_cachefile(self.hmm_cachefname)

        if self.args.action != 'partition':
            self.read_annotation_output(algorithm, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, plotdir=plotdir)

        if not self.args.no_clean and os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self, fname):
        """ Read cached bcrham partition info """
        if self.cached_results is None:
            self.cached_results = {}

        with opener('r')(fname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                # query = line['unique_ids']
                seqstr = line['query_seqs']  # colon-separated list of the query sequences corresponding to this cache
                if 'errors' in line and line['errors'] != '':  # not sure why this needed to be an exception
                    print 'error in bcrham output %s for sequences:' % (line['errors'])
                    print '        %s' % seqstr.replace(':', '\n')

                score, naive_seq, cpos = None, None, None
                if line['score'] != '':
                    score = float(line['score'])
                if line['naive_seq'] != '':
                    naive_seq = line['naive_seq']
                    cpos = int(line['cyst_position'])

                if seqstr in self.cached_results:  # make sure we don't get contradicting info
                    if abs(score - self.cached_results[seqstr]['logprob']) > 1e-4:
                        print 'WARNING unequal logprobs for %s: %f %f' % (seqstr, score, self.cached_results[seqstr]['logprob'])
                    if naive_seq != self.cached_results[seqstr]['naive_seq']:
                        raise Exception('different naive seqs for %s: %s %s' % (seqstr, naive_seq, self.cached_results[seqstr]['naive_seq']))
                    if cpos != self.cached_results[seqstr]['cyst_position']:
                        raise Exception('unequal cyst positions for %s: %d %d' % (seqstr, cpos, self.cached_results[seqstr]['cyst_position']))
                else:
                    self.cached_results[seqstr] = {'logprob' : score, 'naive_seq' : naive_seq, 'cyst_position' : cpos}

    # ----------------------------------------------------------------------------------------
    def write_cachefile(self, fname):
        # write everything we've cached so far to file for bcrham (or a future process) to read
        lockfname = fname + '.lock'
        while os.path.exists(lockfname):
            print '  waiting for lock to clear'
            time.sleep(0.5)
        lockfile = open(lockfname, 'w')
        with opener('w')(fname) as cachefile:
            writer = csv.DictWriter(cachefile, ('query_seqs', 'score', 'naive_seq', 'cyst_position'))
            writer.writeheader()
            for seqstr, cachefo in self.cached_results.items():
                writer.writerow({'query_seqs':seqstr, 'score':cachefo['logprob'], 'naive_seq':cachefo['naive_seq'], 'cyst_position':cachefo['cyst_position']})

        lockfile.close()
        os.remove(lockfname)

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

        print '    processed %d sequences (%d events)' % (n_seqs_processed, n_events_processed)
        if len(boundary_error_queries) > 0:
            print '      %d boundary errors (%s)' % (len(boundary_error_queries), ', '.join(boundary_error_queries))

        if self.args.outfname is not None:
            outpath = self.args.outfname
            if self.args.outfname[0] != '/':  # if full output path wasn't specified on the command line
                outpath = os.getcwd() + '/' + outpath
            shutil.copyfile(self.hmm_outfname, outpath)

        if not self.args.no_clean:
            os.remove(self.hmm_outfname)

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

        if self.args.truncate_n_sets:
            chops = self.get_bcrham_truncations(line)
        for iseq in range(0, len(line['unique_ids'])):
            tmpline = dict(line)
            tmpline['seq'] = line['seqs'][iseq]
            if self.args.truncate_n_sets:
                tmpline['chops'] = chops[iseq]
            label = ilabel if iseq==0 else ''
            event_str = utils.print_reco_event(self.germline_seqs, tmpline, extra_str='    ', return_string=True, label=label, one_line=(iseq>0))
            out_str_list.append(event_str)

            # if iseq == 0:
            #     true_naive = utils.get_full_naive_seq(self.germline_seqs, self.reco_info[tmpline['unique_ids'][iseq]])
            #     inf_naive = utils.get_full_naive_seq(self.germline_seqs, tmpline)
            #     utils.color_mutants(true_naive, inf_naive, print_result=True, extra_str='       mistaken: ')

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
