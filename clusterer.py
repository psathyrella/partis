import sys
import csv
from opener import opener

class Clusterer(object):
    def __init__(self):
        self.debug = False
        self.max_id = -1  # maximum previously used id
        self.cluster_ids = []
        self.query_clusters = {}

    # ----------------------------------------------------------------------------------------
    def cluster(self, infname):
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                if self.debug:
                    print line['unique_id_1'], line['unique_id_2'], float(line['score'])
                self.incorporate_into_clusters(line['unique_id_1'], line['unique_id_2'], float(line['score']))

        print 'unique_id,reco_id'
        for name,cluster_id in self.query_clusters.iteritems():
            print '%s,%d' % (name, cluster_id)

    # ----------------------------------------------------------------------------------------
    def add_new_cluster(self, query_name):
        if self.debug:
            print '    new cluster ',query_name
        assert query_name not in self.query_clusters
        self.max_id += 1
        self.query_clusters[query_name] = self.max_id
        self.cluster_ids.append(self.max_id)

    # ----------------------------------------------------------------------------------------
    def merge_clusters(self, query_name, second_query_name):
        """ move all queries with same id as <second_query_name> to <query_name>'s cluster """
        if self.debug:
            print '    merging ',self.query_clusters[query_name], ' and ',self.query_clusters[second_query_name]
        first_cluster_id = self.query_clusters[query_name]
        second_cluster_id = self.query_clusters[second_query_name]

        if first_cluster_id == second_cluster_id:  # already in the same cluster
            return
        for name,cluster_id in self.query_clusters.iteritems():
            if cluster_id == second_cluster_id:
                self.query_clusters[name] = first_cluster_id

        if self.debug:
            print '  removing ',second_cluster_id, ' from ',self.cluster_ids,
        if second_cluster_id in self.cluster_ids:
            self.cluster_ids.remove(second_cluster_id)
            if self.debug:
                print '  --> ',self.cluster_ids, '   ', self.query_clusters
        else:
            print 'oh, man, something\'s wrong'
            print 'uniqe_id,reco_id'
            for name,cluster_id in self.query_clusters.iteritems():
                print '%s,%d' % (name, cluster_id)
            sys.exit()

    # ----------------------------------------------------------------------------------------
    def add_to_cluster(self, cluster_id, query_name):
        if self.debug:
            print '    adding ',query_name,' to ',cluster_id
        self.query_clusters[query_name] = cluster_id

    # ----------------------------------------------------------------------------------------
    def incorporate_into_clusters(self, query_name, second_query_name, total_score):
        together = total_score > -300.0
        if query_name in self.query_clusters and second_query_name in self.query_clusters:
            if together:
                self.merge_clusters(query_name, second_query_name)
            # else:
            #     assert False  # um, I think so
        elif query_name in self.query_clusters:
            if together:  # add second query to first query's existing cluster
                self.add_to_cluster(self.query_clusters[query_name], second_query_name)
            else:
                self.add_new_cluster(second_query_name)
        elif second_query_name in self.query_clusters:
            if together:
                self.add_to_cluster(self.query_clusters[second_query_name], query_name)
            else:
                self.add_new_cluster(query_name)
        else:
            self.add_new_cluster(query_name)
            if together:
                self.add_to_cluster(self.query_clusters[query_name], second_query_name)
            else:
                self.add_new_cluster(second_query_name)

