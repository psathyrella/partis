import sys
import csv
import math
import itertools
from operator import itemgetter
from subprocess import check_call
from sklearn.metrics.cluster import adjusted_mutual_info_score

import utils
import plotting
from opener import opener

class Glomerator(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        pass

    # ----------------------------------------------------------------------------------------
    def naive_seq_glomerate(self, naive_seqs, n_clusters):
        """ Perform hierarchical agglomeration (with naive hamming distance as the distance), stopping at <n_clusters> """
        clusters = [ [names, ] for names in naive_seqs.keys() ]
        # for seq_a, seq_b in itertools.combinations(naive_seqs.values(), 2):
        #     if len(seq_a) != len(seq_b):
        #         print '  different lengths'
        #         continue
        #     print seq_a, seq_b, utils.hamming(seq_a, seq_b)

        distances = {}
        def glomerate(debug=False):
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
            glomerate(debug=False)

        return clusters

    # ----------------------------------------------------------------------------------------
    def read_cached_agglomeration(self, log_probs=None, partitions=None, infname=None, debug=False, reco_info=None, outfile=None, plotdir='', workdir=None):
        self.max_log_prob, self.best_partition = None, None
        for part in partitions:  # NOTE these are sorted in order of agglomeration, with the initial partition first
            if self.max_log_prob is None or part['score'] > self.max_log_prob:
                self.max_log_prob = part['score']
                self.best_partition = part['clusters']

        if debug:
            print '  best partition ', self.max_log_prob
            print '   clonal?   ids'
            for cluster in self.best_partition:
                same_event = utils.from_same_event(reco_info is None, reco_info, cluster)
                if same_event is None:
                    same_event = -1
                print '     %d    %s' % (int(same_event), ':'.join([ str(uid) for uid in cluster ]))

        self.max_minus_ten_log_prob, self.best_minus_ten_partition = None, None  # reel back glomeration by ten units of log prob to be conservative before we pass to the multiple-process merge
        for part in partitions:
            if part['score'] > self.max_log_prob - 10.0:
                self.max_minus_ten_log_prob = part['score']
                self.best_minus_ten_partition = part['clusters']
                break

        if debug:
            # print '        best minus ten ', self.max_minus_ten_log_prob
            # for cluster in self.best_minus_ten_partition:
            #     print '           ', ':'.join([ str(uid) for uid in cluster ])

            if reco_info is not None:
                true_cluster_list, inferred_cluster_list = [], []
                for iclust in range(len(self.best_partition)):
                    for uid in self.best_partition[iclust]:
                        if uid not in reco_info:
                            raise Exception('ERROR %s' % str(uid))
                        true_cluster_list.append(reco_info[uid]['reco_id'])
                        inferred_cluster_list.append(iclust)
                print '       true clusters %d' % len(set(true_cluster_list))
                print '   inferred clusters %d' % len(set(inferred_cluster_list))
                print '         adjusted mi %.2f' % adjusted_mutual_info_score(true_cluster_list, inferred_cluster_list)

