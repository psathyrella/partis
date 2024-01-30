from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import csv
import math
from subprocess import check_call
import time

from . import utils
from io import open
import six

# ----------------------------------------------------------------------------------------
def vollmers(info, threshold, debug=False):
    """
    From Vollmers paper:
        Lineage Clustering. IGH sequences were clustered into IGH lineages according
        to similarity in their junctional region. Lineages were created according to the
        following steps. A lineage is formed and populated with one IGH sequence (seed). Then, all
        IGH sequences in the lineages (initially only the seed) are compared with all
        other IGH sequences of the same length using the same V and J segments. If
        their junctional regions (untemplated nucleotides and D segments) are at
        least 90% identical, the IGH sequence is added to the lineage. This process is
        repeated until the lineage does not grow.

    NOTE I'm interpreting this to mean
      - if *any* sequence already in the cluster is 90% to the prospective sequence that it's added to the cluster
      - 'sequences the same length' means cdr3 the same length (entire sequence the same length only made sense for their primers
      - since the 90% is on d + insertions, also have to not merge if the d + insertions aren't the same length
    """

    id_clusters = {}  # map from cluster id to list of query names

    def get_d_plus_insertions(uid):
        assert len(info[uid]['d_qr_seqs']) == 1
        return info[uid]['vd_insertion'] + info[uid]['d_qr_seqs'][0] + info[uid]['dj_insertion']

    def get_cdr3_seq(uid):
        cpos = info[uid]['codon_positions']['v']
        tpos = info[uid]['codon_positions']['j']
        assert len(info[uid]['seqs']) == 1
        seq = info[uid]['seqs'][0]
        cdr3_seq = seq[cpos : tpos+3]
        if len(cdr3_seq) != info[uid]['cdr3_length']:
            raise Exception('ERROR bad cdr3 sequence %s %d' % (cdr3_seq, info[uid]['cdr3_length']))
        return cdr3_seq

    def from_same_lineage(cluster_id, uid):
        for clid in id_clusters[cluster_id]:  # loop over seqs already in the cluster (it only has to match one of 'em)
            is_match = True
            for key in ('cdr3_length', 'v_gene', 'j_gene'):  # same cdr3 length, v gene, and j gene
                if info[clid][key] != info[uid][key]:
                    is_match = False
                    break
            if not is_match:
                continue
            cl_seq = get_d_plus_insertions(clid)
            u_seq = get_d_plus_insertions(uid)
            if len(cl_seq) != len(u_seq):
                continue
            hamming_frac = utils.hamming_fraction(cl_seq, u_seq)
            if hamming_frac > 1. - threshold:
                continue

            return True  # if we get to here, it's a match

        return False

    def check_unclustered_seqs():
        """ loop through all unclustered sequences, adding them to the most recently created cluster """
        uids_to_remove = []
        for unique_id in unclustered_seqs:
            assert unique_id not in id_clusters[last_cluster_id]  # not sure why I had the below, but I swear it's impossible. Remove this assertion when it fails to get triggered for a while
            # if unique_id in id_clusters[last_cluster_id]:  # sequence is already in this cluster
            #     continue
            if from_same_lineage(last_cluster_id, unique_id):
                if debug:
                    print('     adding', unique_id)
                id_clusters[last_cluster_id].append(unique_id)
                uids_to_remove.append(unique_id)
        for uid in uids_to_remove:
            unclustered_seqs.remove(uid)

    def add_cluster(clid):
        if debug:
            print('  starting cluster %d' % clid)
        id_clusters[clid] = [unclustered_seqs[0],]
        unclustered_seqs.remove(unclustered_seqs[0])
        while True:
            last_size = len(id_clusters[clid])
            check_unclustered_seqs()
            if last_size == len(id_clusters[clid]):  # stop when cluster stops growing
                break
            if debug:
                print('    running again (%d --> %d)' % (last_size, len(id_clusters[clid])))

    # ----------------------------------------------------------------------------------------
    # the business
    start = time.time()
    if any(len(l['unique_ids']) > 1 for l in info.values()):
        raise Exception('all initial annotations in annotation clustering need to have length 1, but got: %s' % [l['unique_ids'] for l in info.values() if len(l['unique_ids'])>1])
    print('  vj-cdr3 clustering %d sequences with threshold %.3f' % (len(info), threshold))
    unclustered_seqs = list(info.keys())
    last_cluster_id = 0
    while len(unclustered_seqs) > 0:
        add_cluster(last_cluster_id)
        last_cluster_id += 1

    partition = [uids for uids in id_clusters.values()]  # convert to list of lists (no clid info)
    print('  vjcdr3 time: %.1f' % (time.time() - start))
    return partition

# ----------------------------------------------------------------------------------------
def bashford_rogers(info, reco_info=None, debug=False):
    """
    Fig 2a, methods, and Fig S1 here: http://genome.cshlp.org/content/23/11/1874

    to cluster:
      - precluster sequences with more than 95% sequence identity
      - then within each of these clusters, each seq is a node
      - connect all nodes that differ by one base (substitution or indel)
      - clusters are groups of connected nodes
    """

# ----------------------------------------------------------------------------------------
class SingleLinkClusterer(object):
    def __init__(self, threshold=0.0, greater_than=True, singletons=None):  # put in same cluster if greater than threshold, or less than equal to?
        self.threshold = threshold
        self.debug = False
        self.greater_than = greater_than
        self.max_id = -1  # maximum previously used id
        self.cluster_ids = []
        self.query_clusters = {}  # map from query name to cluster id
        self.id_clusters = {}  # map from cluster id to list of query names
        if singletons is None:
            singletons = []
        else:
            for st in singletons:
                self.add_new_cluster(st, dbg_str_list=[])
        self.singletons = singletons
        self.pairscores = {}  # used by external code to see if we saw a given pair
        self.plotscores = {'all':[], 'same':[], 'diff':[]}  # keep track of scores for plotting

    # ----------------------------------------------------------------------------------------
    def add_new_cluster(self, query_name, dbg_str_list):
        dbg_str_list.append('    new cluster ' + str(query_name))
        assert query_name not in self.query_clusters
        self.max_id += 1
        self.query_clusters[query_name] = self.max_id
        self.cluster_ids.append(self.max_id)

    # ----------------------------------------------------------------------------------------
    def merge_clusters(self, query_name, second_query_name, dbg_str_list):
        """ move all queries with same id as <second_query_name> to <query_name>'s cluster """
        if self.query_clusters[query_name] == self.query_clusters[second_query_name]:
            dbg_str_list.append('     already together')
            return
        dbg_str_list.append('     merging ' + str(self.query_clusters[query_name]) + ' and ' + str(self.query_clusters[second_query_name]))
        first_cluster_id = self.query_clusters[query_name]
        second_cluster_id = self.query_clusters[second_query_name]

        if first_cluster_id == second_cluster_id:  # already in the same cluster
            return
        for name, cluster_id in six.iteritems(self.query_clusters):
            if cluster_id == second_cluster_id:
                self.query_clusters[name] = first_cluster_id

        if second_cluster_id in self.cluster_ids:
            self.cluster_ids.remove(second_cluster_id)
        else:
            print('oh, man, something\'s wrong')
            print('uniqe_id,reco_id')
            for name, cluster_id in six.iteritems(self.query_clusters):
                print('%s,%d' % (name, cluster_id))
            sys.exit()

    # ----------------------------------------------------------------------------------------
    def add_to_cluster(self, cluster_id, query_name, dbg_str_list):
        dbg_str_list.append('    adding ' + str(query_name) + ' to ' + str(cluster_id))
        self.query_clusters[query_name] = cluster_id

    # ----------------------------------------------------------------------------------------
    def is_removable(self, score):
        if math.isnan(score):
            assert False
        if self.greater_than:
            return score <= self.threshold
        else:
            return score > self.threshold

    # ----------------------------------------------------------------------------------------
    def incorporate_into_clusters(self, query_name, second_query_name, score, dbg_str_list):
        """ figure out how to add query pair into clusters using single-link"""
        if math.isnan(score):
            print('ERROR nan passed for %d %d (dbg %s)' %(query_name, second_query_name, dbg_str_list))
            sys.exit()
        if self.is_removable(score):
            dbg_str_list.append('    removing link')
            if query_name not in self.query_clusters:
                self.add_new_cluster(query_name, dbg_str_list)
            if second_query_name not in self.query_clusters:
                self.add_new_cluster(second_query_name, dbg_str_list)
            return
        if query_name in self.query_clusters and second_query_name in self.query_clusters:  # if both seqs are already in clusters
            self.merge_clusters(query_name, second_query_name, dbg_str_list)
        elif query_name in self.query_clusters:
            self.add_to_cluster(self.query_clusters[query_name], second_query_name, dbg_str_list)
        elif second_query_name in self.query_clusters:
            self.add_to_cluster(self.query_clusters[second_query_name], query_name, dbg_str_list)
        else:
            self.add_new_cluster(query_name, dbg_str_list)
            self.add_to_cluster(self.query_clusters[query_name], second_query_name, dbg_str_list)

    # ----------------------------------------------------------------------------------------
    def single_link(self, input_scores=None, infname=None, debug=False, reco_info=None, outfile=None):
        if infname is None:
            assert input_scores is not None
        else:
            assert input_scores is None  # should only specify <input_scores> *or* <infname>
            input_scores = []
            with open(infname, 'r') as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    input_scores.append(line)
        sorted_lines = sorted(input_scores, key=lambda k: float(k['logprob']))
        for line in sorted_lines:
            a_name = line['id_a']
            b_name = line['id_b']
            score = float(line['logprob'])
            dbg_str_list = ['%22s %22s   %8.3f' % (a_name, b_name, score), ]
            if reco_info is None:
                dbg_str_list[-1] += '   %s' % ('-')
            else:
                from_same_event = utils.from_same_event(reco_info, [a_name, b_name])
                dbg_str_list[-1] += '   %d' % (from_same_event)
            self.incorporate_into_clusters(a_name, b_name, score, dbg_str_list)
            self.pairscores[(utils.get_key((a_name, b_name)))] = score
            self.plotscores['all'].append(score)
            if reco_info is not None:
                if from_same_event:
                    self.plotscores['same'].append(score)
                else:
                    self.plotscores['diff'].append(score)
            if debug:
                outstr = ''.join(dbg_str_list)
                if outfile == None:
                    print(outstr)
                else:
                    outfile.write(outstr + '\n')

        for query, cluster_id in six.iteritems(self.query_clusters):
            if cluster_id not in self.id_clusters:
                self.id_clusters[cluster_id] = []
            self.id_clusters[cluster_id].append(query)
        for cluster_id, queries in self.id_clusters.items():
            if len(queries) == 1:
                self.singletons.append(queries[0])

        # print 'nearest',self.nearest_true_mate
        out_str_list = ['  %d clusters:\n'%len(self.id_clusters), ]
        for cluster_id in self.id_clusters:
            out_str_list.append('   ' + ' '.join([str(x) for x in self.id_clusters[cluster_id]]) + '\n')
        if outfile == None:
            print(''.join(out_str_list))
        else:
            outfile.write(''.join(out_str_list))
