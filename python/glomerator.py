import itertools
import os
import sys
import csv
from sklearn.metrics.cluster import adjusted_mutual_info_score

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
class Glomerator(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, reco_info=None):
        self.reco_info = reco_info
        self.max_log_probs, self.best_partitions = [], []
        self.max_minus_ten_log_probs, self.best_minus_ten_partitions = [], []
        self.partitions = None

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
            while len(clusters[0]) < len(clusters[-1]) / 3.:
                print 'homogenizing', len(clusters[0]), len(clusters[-1]), len(clusters[-1]) / 3.
                homogenize()
                itries += 1
                if itries > 10:
                    break

        return clusters

    # ----------------------------------------------------------------------------------------
    def print_partition(self, partition, logprob, extrastr=''):
        print '  %s partition %.2f' % (extrastr, logprob)
        print '   clonal?   ids'
        for cluster in partition:
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
    def find_best_partition(self, partitions, debug=False):
        max_log_prob, best_partition = None, None
        for part in partitions:  # NOTE these are sorted in order of agglomeration, with the initial partition first
            # print '%d  %10.2f  %.2f' % (ipath, part['score'], part['adj_mi'])
            if max_log_prob is None or part['score'] > max_log_prob:
                max_log_prob = part['score']
                best_partition = part['clusters']

        if debug:
            self.print_partition(best_partition, max_log_prob, 'best')
            self.mutual_information(best_partition, debug=True)

        return max_log_prob, best_partition

    # ----------------------------------------------------------------------------------------
    def find_best_minus_ten_partition(self, max_log_prob, partitions):
        # find the best-minus-ten, i.e.: reel back glomeration by ten units of log prob to be conservative before we pass to the multiple-process merge
        max_minus_ten_log_probs, best_minus_ten_partitions = None, None
        for part in partitions:
            if part['score'] > max_log_prob - 10.0:
                max_minus_ten_log_prob = part['score']
                best_minus_ten_partition = part['clusters']
                break
        assert len(best_minus_ten_partition) > 0
        return max_minus_ten_log_prob, best_minus_ten_partition

    # ----------------------------------------------------------------------------------------
    def read_file_info(self, infname, n_paths, clean_up):
        partitions = [[] for _ in range(n_paths)]
        with opener('r')(infname) as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                if line['partition'] == '':
                    raise Exception('ERROR null partition (one of the processes probably got passed zero sequences')  # shouldn't happen any more FLW
                # if len(partitions) < int(line['path_index']) + 1:  # should only happen the first time through
                #     partitions.append([])
                #     if debug:
                #         print 'ADD', len(all_partitions), int(line['path_index'])
                uids = []
                for cluster in line['partition'].split(';'):
                    uids.append([unique_id for unique_id in cluster.split(':')])
                partitions[int(line['path_index'])].append({'clusters' : uids,
                                                            'score' : float(line['score']),
                                                            'adj_mi' : self.mutual_information(uids, debug=False)})
                # if int(line['path_index']) == 0 and debug:
                #     print 'appended', float(line['score']), len(all_partitions), len(all_partitions[int(line['path_index'])])
        if clean_up:
            os.remove(infname)
        return partitions

    # ----------------------------------------------------------------------------------------
    def merge_fileinfos(self, fileinfos, smc_particles, debug=False):
        self.partitions = [[] for _ in range(smc_particles)]
        self.combined_conservative_max_minus_ten_logprobs = [0.0 for _ in range(smc_particles)]
        self.combined_conservative_best_minus_ten_partitions = [[] for _ in range(smc_particles)]
        for ipath in range(smc_particles):
            # find the combination of the best-minus-ten for *each* file, which is more conservative
            for ifile in range(len(fileinfos)):
                max_logprob, best_partition = self.find_best_partition(fileinfos[ifile][ipath])
                max_minus_ten_logprob, best_minus_ten_partition = self.find_best_minus_ten_partition(max_logprob, fileinfos[ifile][ipath])
                for cluster in best_minus_ten_partition:
                    self.combined_conservative_best_minus_ten_partitions[ipath].append(cluster)
                self.combined_conservative_max_minus_ten_logprobs[ipath] += max_minus_ten_logprob
            self.print_partition(self.combined_conservative_best_minus_ten_partitions[ipath], self.combined_conservative_max_minus_ten_logprobs[ipath], 'combined conservative')

            if debug:
                print 'IPATH', ipath
            def last_one():
                last = True
                for ifile in range(len(fileinfos)):  # we're finished when all the files are out of glomeration steps (i.e. they all only have one [the last] line left)
                    last &= len(fileinfos[ifile][ipath]) == 1
                return last

            def add_next_global_partition():
                # NOTE this adds a merge for *each* *file*, i.e. each agglomerations step has n_files merges
                if debug:
                    print '  ADD'
                global_partition = []
                global_logprob = 0.0
                for ifile in range(len(fileinfos)):  # combine the first line in each file to make a global partition
                    if debug:
                        print '    ifile', ifile
                    for cluster in fileinfos[ifile][ipath][0]['clusters']:
                        global_partition.append(cluster)
                    global_logprob += fileinfos[ifile][ipath][0]['score']
                    if len(fileinfos[ifile][ipath]) > 1:  # if this isn't the last line (i.e. if there's more glomeration steps in this file), move on to the next line
                        fileinfos[ifile][ipath].pop(0)
                self.partitions[ipath].append({'clusters' : global_partition,
                                          'score' : global_logprob,
                                          'adj_mi' : self.mutual_information(global_partition, debug=False)})

                if debug:
                    for ptn in self.partitions[ipath]:
                        # pl = ':'.join(cl for cl in ptn['clusters'])
                        pl = []
                        for cl in ptn['clusters']:
                            pl.append(':'.join(cl))
                        print '    ', ptn['score'], '  '.join(pl)

            while not last_one():
                add_next_global_partition()
            add_next_global_partition()

    # ----------------------------------------------------------------------------------------
    def read_cached_agglomeration(self, infnames, smc_particles, debug=False, clean_up=True):
        """ Read the partitions output by bcrham. If <all_partitions> is specified, add the info to it """
        fileinfos = []
        for fname in infnames:
            fileinfos.append(self.read_file_info(fname, smc_particles, clean_up))
        self.merge_fileinfos(fileinfos, smc_particles, debug=False)

        for ipath in range(len(self.partitions)):
            logprob, ptn = self.find_best_partition(self.partitions[ipath], debug=True)
            self.max_log_probs.append(logprob)
            self.best_partitions.append(ptn)
            # NOTE you don't want to use this, it's safer to use the combined_conservative stuff
            # logprob, ptn = self.find_best_minus_ten_partition(self.max_log_probs[ipath], self.partitions[ipath])
            # self.max_minus_ten_log_probs.append(logprob)
            # self.best_minus_ten_partitions.append(ptn)

    # ----------------------------------------------------------------------------------------
    def write_partitions(self, outfname, mode):
        with open(outfname, mode) as outfile:
            writer = csv.DictWriter(outfile, ('path_index', 'score', 'normalized_score', 'adj_mi'))
            writer.writeheader()
            for ipath in range(len(self.partitions)):
                for part in self.partitions[ipath]:
                    writer.writerow({'path_index' : ipath,
                                     'score' : part['score'],
                                     'normalized_score' : part['score'] / self.max_log_probs[ipath],
                                     'adj_mi' : part['adj_mi']})
