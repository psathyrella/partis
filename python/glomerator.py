import itertools
import os
import sys
import math
import csv
from sklearn.metrics.cluster import adjusted_mutual_info_score

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
class ClusterPath(object):
    def __init__(self, initial_path_index):
        self.partitions = []
        self.logprobs = []
        self.logweights = []
        self.adj_mis = []
        self.initial_path_index = initial_path_index
        self.max_logprob, self.max_minus_ten_logprob, self.max_minus_ten_logweight = None, None, None
        self.best_adj_mi = None
        self.best_partition, self.best_minus_ten_partition = None, None

    def update_best_minus_ten_partition(self):
        for ip in range(len(self.partitions)):  # they should be in order of increasing logprob
            if self.logprobs[ip] > self.max_logprob - 10.:
                self.max_minus_ten_logprob = self.logprobs[ip]
                self.max_minus_ten_logweight = self.logweights[ip]
                self.best_minus_ten_partition = self.partitions[ip]
                break

    def add_partition(self, partition, logprob, logweight, adj_mi):
        self.partitions.append(partition)
        self.logprobs.append(logprob)
        if len(self.logprobs) > 1:
            assert self.logprobs[-1] >= self.logprobs[-2]
        self.logweights.append(logweight)
        self.adj_mis.append(adj_mi)
        if self.max_logprob is None or logprob > self.max_logprob:
            self.max_logprob = logprob
            self.best_partition = partition
            self.best_adj_mi = adj_mi
            self.update_best_minus_ten_partition()

    def remove_first_partition(self):
        self.partitions.pop(0)
        self.logprobs.pop(0)
        self.logweights.pop(0)
        self.adj_mis.pop(0)

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
    def print_partition(self, partition, logprob, adj_mi=-1., extrastr=''):
        print '  %s partition   %-15.2f    %-8.2f' % (extrastr, logprob, adj_mi)
        print '   clonal?   ids'
        for cluster in partition:
            same_event = utils.from_same_event(self.reco_info is None, self.reco_info, cluster)
            if same_event is None:
                same_event = -1
            print '     %d    %s' % (int(same_event), ':'.join([str(uid) for uid in cluster]))

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

    # # ----------------------------------------------------------------------------------------
    # def find_best_partition(self, partitions, debug=False):
    #     if len(partitions) == 0:
    #         raise Exception('zero length partition')  # not sure *why*...

    #     max_log_prob, best_partition = None, None
    #     for part in partitions:  # NOTE these are sorted in order of agglomeration, with the initial partition first
    #         # print '%d  %10.2f  %.2f' % (ipath, part['score'], part['adj_mi'])
    #         if max_log_prob is None or part['score'] > max_log_prob:
    #             max_log_prob = part['score']
    #             best_partition = part['clusters']

    #     # if debug:
    #     #     self.print_partition(best_partition, max_log_prob, 'best')
    #     #     self.mutual_information(best_partition, debug=True)
    #     if max_log_prob is None:
    #         print partitions
    #         raise Exception('couldn\'t find max log prob')

    #     return max_log_prob, best_partition

    # # ----------------------------------------------------------------------------------------
    # def find_best_minus_ten_partition(self, path):
    #     # find the best-minus-ten, i.e.: reel back glomeration by ten units of log prob to be conservative before we pass to the multiple-process merge
    #     max_minus_ten_log_prob, best_minus_ten_partition = None, None
    #     for part in partitions:  # they're 
    #         if part['score'] > max_log_prob - 10.0:
    #             max_minus_ten_log_prob = part['score']
    #             best_minus_ten_partition = part['clusters']
    #             break

    #     if best_minus_ten_partition is None or len(best_minus_ten_partition) == 0:
    #         print max_log_prob
    #         print partitions
    #         raise Exception('couldn\'t find max minus ten log prob')

    #     return max_minus_ten_log_prob, best_minus_ten_partition

    # ----------------------------------------------------------------------------------------
    def read_file_info(self, infname, n_paths, clean_up):
        paths = [None for _ in range(n_paths)]
        # partitions = [[] for _ in range(n_paths)]
        # initial_path_indices = {}
        with opener('r')(infname) as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                if line['partition'] == '':
                    raise Exception('ERROR null partition (one of the processes probably got passed zero sequences')  # shouldn't happen any more FLW
                uids = []
                for cluster in line['partition'].split(';'):
                    uids.append([unique_id for unique_id in cluster.split(':')])
                path_index = int(line['path_index'])
                # partitions[path_index].append({'clusters' : uids,
                #                                'score' : float(line['score']),
                #                                'logweight' : float(line['logweight']),
                #                                'adj_mi' : self.mutual_information(uids, debug=False)})
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

    # # ----------------------------------------------------------------------------------------
    # def combine_log_weights(self, logw1, logw2):
    #     n_ways1 = 1. / math.exp(logw1) if logw1 is not None else 0
    #     n_ways2 = 1. / math.exp(logw2)
    #     return math.log(1. / (n_ways1 + n_ways2))

    # ----------------------------------------------------------------------------------------
    def merge_fileinfos(self, fileinfos, smc_particles, debug=False):
        self.paths = [ClusterPath(None) for _ in range(smc_particles)]  # each path's initial_path_index is None since we're merging paths that, in general, have different initial path indices
        for ipath in range(smc_particles):

            combined_conservative_max_minus_ten_logprob = 0.
            # combined_conservative_max_minus_ten_logweight = None
            combined_conservative_best_minus_ten_partition = []
            total_n_ways = 0
            # find the combination of the best-minus-ten for *each* file, which is more conservative than combining the files and *then* rewinding by 10.
            for ifile in range(len(fileinfos)):
                for cluster in fileinfos[ifile][ipath].best_minus_ten_partition:
                    # first make sure we didn't already add any of the uids in <cluster>
                    for uid in cluster:
                        for existing_cluster in combined_conservative_best_minus_ten_partition:
                            if uid in existing_cluster:
                                raise Exception('%s already in cluster %s' % (uid, ':'.join([qn for qn in cluster])))
                    # then append
                    combined_conservative_best_minus_ten_partition.append(cluster)
                combined_conservative_max_minus_ten_logprob += fileinfos[ifile][ipath].max_minus_ten_logprob
                # print 'combine %f %f --> %f %f' % (combined_conservative_max_minus_ten_logweight, fileinfos[ifile][ipath].max_minus_ten_logweight,
                #                                    1./math.exp(combined_conservative_max_minus_ten_logweight), 1./math.exp(fileinfos[ifile][ipath].max_minus_ten_logweight))
                total_n_ways += 1. / math.exp(fileinfos[ifile][ipath].max_minus_ten_logweight)
                # combined_conservative_max_minus_ten_logweight = self.combine_log_weights(combined_conservative_max_minus_ten_logweight, fileinfos[ifile][ipath].max_minus_ten_logweight)
            combined_conservative_max_minus_ten_logweight = math.log(1. / total_n_ways)

            # then merge all the steps in each path
            if debug:
                print 'IPATH', ipath
            def last_one():
                last = True
                for ifile in range(len(fileinfos)):  # we're finished when all the files are out of glomeration steps (i.e. they all only have one [the last] line left)
                    last &= len(fileinfos[ifile][ipath].partitions) == 1
                return last

            def add_next_global_partition():
                # NOTE this adds a merge for *each* *file*, i.e. each agglomerations step has n_files merges
                if debug:
                    print '  ADD'
                global_partition = []
                global_logprob, global_logweight = 0., 0.
                for ifile in range(len(fileinfos)):  # combine the first line in each file to make a global partition
                    if debug:
                        print '    ifile', ifile
                    for cluster in fileinfos[ifile][ipath].partitions[0]:
                        global_partition.append(cluster)
                    global_logprob += fileinfos[ifile][ipath].logprobs[0]
                    global_logweight += fileinfos[ifile][ipath].logweights[0]
                    if len(fileinfos[ifile][ipath].partitions) > 1:  # if this isn't the last line (i.e. if there's more glomeration steps in this file), move on to the next line
                        fileinfos[ifile][ipath].remove_first_partition()
                self.paths[ipath].add_partition(global_partition, global_logprob, global_logweight,
                                                self.mutual_information(global_partition, debug=False))

                if debug:
                    for ptn in self.paths[ipath].partitions:
                        pl = []
                        for cl in ptn:
                            pl.append(':'.join(cl))
                        print '    ', ptn['score'], '  '.join(pl)

            while not last_one():
                add_next_global_partition()
            add_next_global_partition()
            self.paths[ipath].best_minus_ten_partition = combined_conservative_best_minus_ten_partition  # replace the default one with the more conservative one
            self.paths[ipath].max_minus_ten_logprob = combined_conservative_max_minus_ten_logprob
            self.paths[ipath].max_minus_ten_logweight = combined_conservative_max_minus_ten_logweight

    # ----------------------------------------------------------------------------------------
    def read_cached_agglomeration(self, infnames, smc_particles, debug=False, clean_up=True):
        """ Read the partitions output by bcrham. If <all_partitions> is specified, add the info to it """
        fileinfos = []
        for fname in infnames:
            fileinfos.append(self.read_file_info(fname, smc_particles, clean_up))
        self.merge_fileinfos(fileinfos, smc_particles, debug=False)

        return self.paths

        # for ipath in range(len(self.partitions)):
        #     logprob, ptn = self.find_best_partition(self.partitions[ipath], debug=True)
        #     self.max_log_probs.append(logprob)
        #     self.best_partitions.append(ptn)

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, outfname, mode):
        with open(outfname, mode) as outfile:
            writer = csv.DictWriter(outfile, ('path_index', 'score', 'normalized_score', 'adj_mi'))
            writer.writeheader()
            for ipath in range(len(self.paths)):
                for ipart in range(len(self.paths[ipath].partitions)):
                    writer.writerow({'path_index' : ipath,
                                     'score' : self.paths[ipath].logprobs[ipart],
                                     'logweight' : self.paths[ipath].logweights[ipart],
                                     # 'normalized_score' : part['score'] / self.max_log_probs[ipath],
                                     'adj_mi' : self.paths[ipath].adj_mis[ipart]})
