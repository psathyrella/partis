import itertools
import os
import sys
import math
import csv
from sklearn.metrics.cluster import adjusted_mutual_info_score

import utils
from opener import opener
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class Glomerator(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, reco_info=None):
        self.reco_info = reco_info
        # self.max_log_probs, self.best_partitions = [], []
        # self.max_minus_ten_log_probs, self.best_minus_ten_partitions = [], []
        self.paths = None

    # ----------------------------------------------------------------------------------------
    def naive_seq_glomerate(self, naive_seqs, n_clusters, debug=False):
        """ Perform hierarchical agglomeration (with naive hamming distance as the distance), stopping at <n_clusters> """
        clusters = [[names,] for names in naive_seqs.keys()]
        # for seq_a, seq_b in itertools.combinations(naive_seqs.values(), 2):
        #     if len(seq_a) != len(seq_b):
        #         print '  different lengths'
        #         continue
        #     print seq_a, seq_b, utils.hamming(seq_a, seq_b)

        distances = {}
        def glomerate():
            smallest_min_distance = None
            clusters_to_merge = None
            for clust_a, clust_b in itertools.combinations(clusters, 2):
                min_distance = None  # find the minimal hamming distance between any two sequences in the two clusters
                for query_a in clust_a:
                    for query_b in clust_b:
                        joint_key = ';'.join(sorted([query_a, query_b]))
                        if joint_key not in distances:
                            distances[joint_key] = utils.hamming(naive_seqs[query_a], naive_seqs[query_b])
                        if debug:
                            print '    %25s %25s   %4d   (%s)' % (query_a, query_b, distances[joint_key], joint_key)
                        if min_distance is None or distances[joint_key] < min_distance:
                            min_distance = distances[joint_key]
                if smallest_min_distance is None or min_distance < smallest_min_distance:
                    smallest_min_distance = min_distance
                    clusters_to_merge = (clust_a, clust_b)
            if debug:
                print 'merging', clusters_to_merge
            clusters.append(clusters_to_merge[0] + clusters_to_merge[1])
            clusters.remove(clusters_to_merge[0])
            clusters.remove(clusters_to_merge[1])

        while len(clusters) > n_clusters:
            glomerate()

        # roughly equalize the cluster sizes
        if len(clusters) > 1:
            # mean_length = sum([len(cl) for cl in clusters]) / float(len(clusters))
            clusters.sort(key=len)

            def homogenize():
                # print 'before'
                # for cl in clusters:
                #     print len(cl),
                # print ''
                if len(clusters) > 2:
                    clusters[0] = clusters[0] + clusters[1]
                    clusters[1] = clusters[-1][ : len(clusters[-1])/2]
                    clusters[-1] = clusters[-1][len(clusters[-1])/2 : ]
                else:  # only two clusters
                    together = clusters[0] + clusters[1]
                    clusters[0] = together[ : len(together)/2]
                    clusters[1] = together[len(together)/2 : ]
                # print 'after'
                # for cl in clusters:
                #     print len(cl),
                # print ''
                clusters.sort(key=len)

            itries = 0
            while len(clusters[0]) < len(clusters[-1]) / 1.75:
                if debug:
                    print '  homogenizing', len(clusters[0]), len(clusters[-1])  #, len(clusters[-1]) / 3.
                homogenize()
                itries += 1
                if itries > 100:
                    break

        return clusters

    # ----------------------------------------------------------------------------------------
    def print_partitions(self, path, extrastr='', one_line=False):
        for ip in range(len(path.partitions)):
            self.print_partition(path.partitions[ip], path.logprobs[ip], path.adj_mis[ip], extrastr=extrastr, one_line=one_line)

    # ----------------------------------------------------------------------------------------
    def print_partition(self, partition, logprob, adj_mi=-1., extrastr='', one_line=False):
        if one_line:
            print '    %s %-15.2f   %-8.3f   %3d clusters (%d ways):' % (extrastr, logprob, adj_mi, len(partition), -1),  #1./math.exp(self.logweights[ip])),
        else:
            print '  %s partition   %-15.2f    %-8.2f' % (extrastr, logprob, adj_mi)
            print '   clonal?   ids'
        for cluster in partition:
            same_event = utils.from_same_event(self.reco_info is None, self.reco_info, cluster)
            if same_event is None:
                same_event = -1
            cluster_str = ':'.join([str(uid) for uid in cluster])
            if one_line:
                print '   %s' % cluster_str,
            else:
                print '     %d    %s' % (int(same_event), cluster_str)
        if one_line:
            print ''

    # ----------------------------------------------------------------------------------------
    def print_true_partition(self):
        print '  true partition'
        print '   clonal?   ids'
        true_partition = {}
        for key in self.reco_info:
            reco_id = self.reco_info[key]['reco_id']
            if reco_id not in true_partition:
                true_partition[reco_id] = []
            true_partition[reco_id].append(self.reco_info[key]['unique_id'])
        for cluster in true_partition.values():
            same_event = utils.from_same_event(self.reco_info is None, self.reco_info, cluster)
            if same_event is None:
                same_event = -1
            print '     %d    %s' % (int(same_event), ':'.join([str(uid) for uid in cluster]))

    # ----------------------------------------------------------------------------------------
    def mutual_information(self, partition, debug=False):
        if self.reco_info is None:
            return -1.0
        true_cluster_list, inferred_cluster_list = [], []
        for iclust in range(len(partition)):
            for uid in partition[iclust]:
                if uid not in self.reco_info:
                    raise Exception('ERROR %s' % str(uid))
                true_cluster_list.append(self.reco_info[uid]['reco_id'])
                inferred_cluster_list.append(iclust)
        adj_mi = adjusted_mutual_info_score(true_cluster_list, inferred_cluster_list)
        if debug:
            print '       true clusters %d' % len(set(true_cluster_list))
            print '   inferred clusters %d' % len(set(inferred_cluster_list))
            print '         adjusted mi %.2f' % adj_mi
        return adj_mi

    # ----------------------------------------------------------------------------------------
    def read_file_info(self, infname, n_paths, clean_up):
        paths = [None for _ in range(n_paths)]
        with opener('r')(infname) as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                if line['partition'] == '':
                    raise Exception('ERROR null partition (one of the processes probably got passed zero sequences')  # shouldn't happen any more FLW
                uids = []
                for cluster in line['partition'].split(';'):
                    uids.append([unique_id for unique_id in cluster.split(':')])
                path_index = int(line['path_index'])
                if paths[path_index] is None:
                    paths[path_index] = ClusterPath(int(line['initial_path_index']))
                else:
                    assert paths[path_index].initial_path_index == int(line['initial_path_index'])
                paths[path_index].add_partition(uids, float(line['score']), float(line['logweight']), self.mutual_information(uids, debug=False))

        for cp in paths:
            for ptn in cp.partitions:
                if len(ptn) == 0:
                    raise Exception('zero length partition read from %s' % infname)

        if clean_up:
            os.remove(infname)

        return paths

    # ----------------------------------------------------------------------------------------
    def merge_partitions(self, partitions, logprobs, logweights, adj_mis):
        pass

    # ----------------------------------------------------------------------------------------
    def merge_fileinfos(self, fileinfos, smc_particles, previous_info=None, debug=False):
        self.paths = [ClusterPath(None) for _ in range(smc_particles)]  # each path's initial_path_index is None since we're merging paths that, in general, have different initial path indices

        if previous_info is not None:  # DEAR FUTURE SELF this won't make any sense until you find that picture you took of the white board
            assert len(previous_info) == len(fileinfos)  # both are the number of processes we're merging into one
            for ifile in range(len(fileinfos)):
                if debug:
                    print 'ifile', ifile
                for ipath in range(smc_particles):
                    if debug:
                        print '  ipath', ipath
                        print '    before'
                        self.print_partitions(fileinfos[ifile][ipath], one_line=True)
                    initial_path_index = fileinfos[ifile][ipath].initial_path_index  # which previous path are we hooking up to?
                    previous_path = previous_info[ifile][initial_path_index]
                    current_path = fileinfos[ifile][ipath]
                    first_new_logprob = current_path.logprobs[0]
                    extended_path = ClusterPath(None)
                    for ip in range(len(previous_path.partitions)):
                        # if previous_path.logprobs[ip] >= first_new_logprob:  # skip the merges past which we rewound
                        #     continue
                        extended_path.add_partition(previous_path.partitions[ip], previous_path.logprobs[ip], previous_path.logweights[ip], previous_path.adj_mis[ip])
                    for ip in range(len(current_path.partitions)):
                        extended_path.add_partition(current_path.partitions[ip], current_path.logprobs[ip], current_path.logweights[ip], current_path.adj_mis[ip])
                    fileinfos[ifile][ipath] = extended_path
                    if debug:
                        print '    after'
                        self.print_partitions(fileinfos[ifile][ipath], one_line=True)

        if debug:
            print '   find combined conservative'
        for ipath in range(smc_particles):
            combined_conservative_max_minus_ten_logprob = 0.
            combined_conservative_best_minus_ten_partition = []
            total_n_ways = 0
            # find the combination of the best-minus-ten for *each* file, which is more conservative than combining the files and *then* rewinding by 10.
            for ifile in range(len(fileinfos)):
                if debug:
                    # self.print_partitions(fileinfos[ifile][ipath], extrastr=('%d' % (ifile)), one_line=True)
                    fileinfos[ifile][ipath].print_partitions(self.reco_info, extrastr=('%d' % (ifile)), one_line=True)
                for cluster in fileinfos[ifile][ipath].best_minus_ten_partition:
                    # first make sure we didn't already add any of the uids in <cluster>
                    for uid in cluster:
                        for existing_cluster in combined_conservative_best_minus_ten_partition:
                            if uid in existing_cluster:
                                raise Exception('%s already in cluster %s' % (uid, ':'.join([qn for qn in cluster])))
                    # then append
                    combined_conservative_best_minus_ten_partition.append(cluster)
                combined_conservative_max_minus_ten_logprob += fileinfos[ifile][ipath].max_minus_ten_logprob
                # total_n_ways += 1. / math.exp(fileinfos[ifile][ipath].max_minus_ten_logweight)

            # combined_conservative_max_minus_ten_logweight = math.log(1. / total_n_ways)

            # then merge all the steps in each path
            if debug:
                print 'merge path %d with %d processes' % (ipath, len(fileinfos))
            def last_one():
                last = True
                for ifile in range(len(fileinfos)):  # we're finished when all the files are out of glomeration steps (i.e. they all only have one [the last] line left)
                    last &= len(fileinfos[ifile][ipath].partitions) == 1
                return last

            def remove_one_of_the_first_partitions():
                maxdelta, ibestfile = None, None
                for ifile in range(len(fileinfos)):
                    if len(fileinfos[ifile][ipath].partitions) == 1:  # if this is the last line (i.e. there aren't any more glomeration steps in this file), leave it alone
                        continue
                    thisdelta = fileinfos[ifile][ipath].logprobs[1] - fileinfos[ifile][ipath].logprobs[0]  # logprob difference between the next partition and this one
                    if maxdelta is None or thisdelta > maxdelta:
                        maxdelta = thisdelta
                        ibestfile = ifile
                fileinfos[ibestfile][ipath].remove_first_partition()

            def add_next_global_partition():
                global_partition = []
                global_logprob = 0.
                # total_n_ways = 0
                for ifile in range(len(fileinfos)):  # combine the first line in each file to make a global partition
                    for cluster in fileinfos[ifile][ipath].partitions[0]:
                        global_partition.append(cluster)
                    global_logprob += fileinfos[ifile][ipath].logprobs[0]
                    # total_n_ways += 1. / math.exp(fileinfos[ifile][ipath].logweights[0])
                # global_logweight = math.log(1. / total_n_ways)
                global_adj_mi = self.mutual_information(global_partition, debug=False)
                self.paths[ipath].add_partition(global_partition, global_logprob, 0., global_adj_mi)  # don't know the logweight yet (or maybe at all!)

                # if debug:
                #     self.print_partition(global_partition, global_logprob, global_adj_mi, one_line=True)

            while not last_one():
                add_next_global_partition()
                remove_one_of_the_first_partitions()
            add_next_global_partition()

            self.paths[ipath].set_synthetic_logweight_history(self.reco_info)
            if debug:
                self.paths[ipath].print_partitions(self.reco_info, one_line=True)

            # replace the default one with the more conservative one
            # NOTE the combined conservative partition doesn't necessarily occur in the merged string
            self.paths[ipath].best_minus_ten_partition = combined_conservative_best_minus_ten_partition
            self.paths[ipath].max_minus_ten_logprob = combined_conservative_max_minus_ten_logprob
            # self.paths[ipath].max_minus_ten_logweight = combined_conservative_max_minus_ten_logweight
            if debug:
                print '  combined conservative'
                self.print_partition(self.paths[ipath].best_minus_ten_partition, self.paths[ipath].max_minus_ten_logprob, one_line=True)

    # ----------------------------------------------------------------------------------------
    def read_cached_agglomeration(self, infnames, smc_particles, previous_info=None, debug=False, clean_up=True):
        """ Read the partitions output by bcrham. If <all_partitions> is specified, add the info to it """
        fileinfos = []
        for fname in infnames:
            fileinfos.append(self.read_file_info(fname, smc_particles, clean_up))
        self.merge_fileinfos(fileinfos, smc_particles, previous_info=previous_info, debug=True)

        return self.paths

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, outfname, mode, paths):
        with open(outfname, mode) as outfile:
            writer = csv.DictWriter(outfile, ('path_index', 'score', 'logweight', 'adj_mi'))
            writer.writeheader()
            for ipath in range(len(paths)):
                for ipart in range(len(paths[ipath].partitions)):
                    writer.writerow({'path_index' : ipath,
                                     'score' : paths[ipath].logprobs[ipart],
                                     'logweight' : paths[ipath].logweights[ipart],
                                     # 'normalized_score' : part['score'] / self.max_log_probs[ipath],
                                     'adj_mi' : paths[ipath].adj_mis[ipart]})
