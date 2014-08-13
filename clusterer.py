import sys
import csv
import math
from opener import opener
# ./venv/bin/linsim compare-clustering --true-name-column unique_id --inferred-name-column unique_id  --true-group-column reco_id --inferred-group-column reco_id /tmp/dralph/true.csv /tmp/dralph/inf.csv 

class Clusterer(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, threshold, greater_than=True):  # put in same cluster if greater than threshold, or less than equal to?
        self.threshold = threshold
        self.debug = False
        self.greater_than = greater_than
        self.max_id = -1  # maximum previously used id
        self.cluster_ids = []
        self.query_clusters = {}  # map from query name to cluster id
        self.id_clusters = {}  # map from cluster id to query name list

    # ----------------------------------------------------------------------------------------
    def cluster(self, infname, debug=False):
        self.debug = debug
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                if self.debug:
                    print '%22s %22s   %.3f' % (line['unique_id'], line['second_unique_id'], float(line['score'])),
                self.incorporate_into_clusters(line['unique_id'], line['second_unique_id'], float(line['score']))
                if self.debug:
                    print ''

        for query, cluster_id in self.query_clusters.iteritems():
            if cluster_id not in self.id_clusters:
                self.id_clusters[cluster_id] = []
            self.id_clusters[cluster_id].append(query)
        
        if True:  #self.debug:
            for cluster_id in self.id_clusters:
                print self.id_clusters[cluster_id]
            # print 'unique_id,reco_id'
            # for name,cluster_id in self.query_clusters.iteritems():
            #     print '%s,%d' % (name, cluster_id)

    # ----------------------------------------------------------------------------------------
    def add_new_cluster(self, query_name):
        if self.debug:
            print '    new cluster ',query_name,
        assert query_name not in self.query_clusters
        self.max_id += 1
        self.query_clusters[query_name] = self.max_id
        self.cluster_ids.append(self.max_id)

    # ----------------------------------------------------------------------------------------
    def merge_clusters(self, query_name, second_query_name):
        """ move all queries with same id as <second_query_name> to <query_name>'s cluster """
        if self.debug:
            print '     merging ',self.query_clusters[query_name], ' and ',self.query_clusters[second_query_name],
        first_cluster_id = self.query_clusters[query_name]
        second_cluster_id = self.query_clusters[second_query_name]

        if first_cluster_id == second_cluster_id:  # already in the same cluster
            return
        for name,cluster_id in self.query_clusters.iteritems():
            if cluster_id == second_cluster_id:
                self.query_clusters[name] = first_cluster_id

        if second_cluster_id in self.cluster_ids:
            self.cluster_ids.remove(second_cluster_id)
        else:
            print 'oh, man, something\'s wrong'
            print 'uniqe_id,reco_id'
            for name,cluster_id in self.query_clusters.iteritems():
                print '%s,%d' % (name, cluster_id)
            sys.exit()

    # ----------------------------------------------------------------------------------------
    def add_to_cluster(self, cluster_id, query_name):
        if self.debug:
            print '    adding ',query_name,' to ',cluster_id,
        self.query_clusters[query_name] = cluster_id

    # ----------------------------------------------------------------------------------------
    def is_removable(self, score):
        if score == float('nan') or score == float('-nan'):
            assert False
            return True
        if self.greater_than:
            return score <= self.threshold
        else:
            return score > self.threshold

    # ----------------------------------------------------------------------------------------
    def incorporate_into_clusters(self, query_name, second_query_name, score):
        if query_name in self.query_clusters and second_query_name in self.query_clusters:  # if both seqs are already in clusters
            if self.is_removable(score):
                if self.debug:
                    print '    removing link',
            else:
                if self.query_clusters[query_name] != self.query_clusters[second_query_name]:
                    self.merge_clusters(query_name, second_query_name)
                else:
                    if self.debug:
                        print '     already together',
        elif query_name in self.query_clusters:
            if self.is_removable(score):
                self.add_new_cluster(second_query_name)
            else:
                self.add_to_cluster(self.query_clusters[query_name], second_query_name)
        elif second_query_name in self.query_clusters:
            if self.is_removable(score):
                self.add_new_cluster(query_name)
            else:
                self.add_to_cluster(self.query_clusters[second_query_name], query_name)
        else:
            if self.is_removable(score):
                self.add_new_cluster(query_name)
                self.add_new_cluster(second_query_name)
            else:
                self.add_new_cluster(query_name)
                self.add_to_cluster(self.query_clusters[query_name], second_query_name)
