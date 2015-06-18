import sys
import csv
import math
from subprocess import check_call
from sklearn.metrics.cluster import adjusted_mutual_info_score

import utils
import plotting
from opener import opener

class Clusterer(object):
    # ----------------------------------------------------------------------------------------
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

        # self.nearest_true_mate = {}  #

    # ----------------------------------------------------------------------------------------
    def single_link(self, input_scores=None, infname=None, debug=False, reco_info=None, outfile=None, plotdir=''):
        if infname is None:
            assert input_scores is not None
        else:
            assert input_scores is None  # should only specify <input_scores> *or* <infname>
            input_scores = []
            with opener('r')(infname) as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    input_scores.append(line)
        sorted_lines = sorted(input_scores, key=lambda k: float(k['logprob']))
        for line in sorted_lines:
            a_name = line['id_a']
            b_name = line['id_b']
            score = float(line['logprob'])
            from_same_event = -1 if (reco_info == None or a_name not in reco_info or b_name not in reco_info) else reco_info[a_name]['reco_id'] == reco_info[b_name]['reco_id']
            dbg_str_list = ['%22s %22s   %8.3f   %d' % (a_name, b_name, score, from_same_event), ]
            self.incorporate_into_clusters(a_name, b_name, score, dbg_str_list)
            self.pairscores[(utils.get_key((a_name, b_name)))] = score
            self.plotscores['all'].append(score)
            if reco_info != None:
                if from_same_event:
                    self.plotscores['same'].append(score)
                else:
                    self.plotscores['diff'].append(score)
            # if reco_info != None and reco_info[a_name]['reco_id'] == reco_info[b_name]['reco_id']:
            #     for query,score in {a_name:score, b_name:score}.iteritems():
            #         if query not in self.nearest_true_mate:
            #             self.nearest_true_mate[query] = score
            #         elif self.greater_than and score > self.nearest_true_mate[query]:
            #             self.nearest_true_mate[query] = score
            #         elif not self.greater_than and score < self.nearest_true_mate[query]:
            #             self.nearest_true_mate[query] = score
            if debug:
                outstr = ''.join(dbg_str_list)
                if outfile == None:
                    print outstr
                else:
                    outfile.write(outstr + '\n')

        if plotdir != '':
            utils.prep_dir(plotdir + '/plots', '*.svg')
            hists = {}
            for htype in ['all', 'same', 'diff']:
                raise Exception('deprecated')
                # # ----------------------------------------------------------------------------------------
                # def make_hist_from_list(values, hist_label, n_bins=30):
                #     """ Fill a histogram with float values in a list """
                #     if len(values) == 0:
                #         print 'WARNING no values for %s in make_hist' % hist_label
                #         return TH1D(hist_label, '', 1, 0, 1)
                #     xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
                #     set_bins(values, n_bins, is_log_x=False, xbins=xbins, var_type='float')
                #     hist = TH1D(hist_label, '', n_bins, xbins)
                #     for val in values:
                #         hist.Fill(val)
                #     return hist
                # ----------------------------------------------------------------------------------------
                # hists[htype] = plotting.make_hist_from_list(self.plotscores[htype], htype + '_pairscores')
                # hists[htype].SetTitle(htype)
            plotting.draw(hists['all'], 'float', plotdir=plotdir, plotname='pairscores', more_hists=[hists['same'], hists['diff']])
            check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
            check_call(['./bin/permissify-www', plotdir])

        for query, cluster_id in self.query_clusters.iteritems():
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
            print ''.join(out_str_list)
        else:
            outfile.write(''.join(out_str_list))

    # ----------------------------------------------------------------------------------------
    def vollmers_cluster(self, info, reco_info=None):
        # ./bin/partis.py --action run-viterbi --vollmers-clustering --seqfile test/regression/parameters/simu.csv --parameter-dir test/regression/parameters/simu/hmm --n-max-queries -1 --n-procs 10  --debug 0 --truncate-pairs
        """
        Cluster together sequences with similar rearrangement parameters

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

        def get_d_plus_insertions(uid):
            return info[uid]['vd_insertion'] + info[uid]['d_qr_seq'] + info[uid]['dj_insertion']

        def get_cdr3_seq(uid):
            cpos = info[uid]['cyst_position']
            tpos = info[uid]['tryp_position']
            assert len(info[uid]['seqs']) == 1
            seq = info[uid]['seqs'][0]
            cdr3_seq = seq[cpos : tpos+3]
            if len(cdr3_seq) != info[uid]['cdr3_length']:
                raise Exception('ERROR bad cdr3 sequence %s %d' % (cdr3_seq, info[uid]['cdr3_length']))
            return cdr3_seq

        def from_same_lineage(cluster_id, uid):
            for clid in self.id_clusters[cluster_id]:  # loop over seqs already in the cluster (it only has to match one of 'em)
                is_match = True
                for key in ('cdr3_length', 'v_gene', 'j_gene'):  # same cdr3 length, v gene, and d gene
                    if info[clid][key] != info[uid][key]:
                        # print '      %s doesn\'t match' % key
                        is_match = False
                if not is_match:
                    continue
                # cl_cdr3_seq = get_cdr3_seq(clid)
                # u_cdr3_seq = get_cdr3_seq(uid)
                # print utils.color_mutants(cl_cdr3_seq, u_cdr3_seq, print_result=False)
                cl_seq = get_d_plus_insertions(clid)
                u_seq = get_d_plus_insertions(uid)
                if len(cl_seq) != len(u_seq):
                    continue
                hamming_frac = utils.hamming_fraction(cl_seq, u_seq)
                if hamming_frac > 0.1:  # if cdr3 is more than 10 percent different we got no match
                    # print '      hamming too large', hamming_frac
                    continue

                return True  # if we get to here, it's a match

            return False

        def check_unclustered_seqs():
            uids_to_remove = []
            for unique_id in unclustered_seqs:
                if unique_id in self.id_clusters[last_cluster_id]:  # sequence is already in this cluster
                    continue
                if from_same_lineage(last_cluster_id, unique_id):
                    print '     adding', unique_id
                    self.id_clusters[last_cluster_id].append(unique_id)
                    uids_to_remove.append(unique_id)
            for uid in uids_to_remove:
                unclustered_seqs.remove(uid)

        def add_cluster(clid):
            print '  starting cluster %d' % clid
            self.id_clusters[clid] = [unclustered_seqs[0],]
            unclustered_seqs.remove(unclustered_seqs[0])
            while True:
                last_size = len(self.id_clusters[clid])
                check_unclustered_seqs()
                if last_size == len(self.id_clusters[clid]):  # stop when cluster stops growing
                    break
                print '    running again (%d --> %d)' % (last_size, len(self.id_clusters[clid]))

        unclustered_seqs = info.keys()
        last_cluster_id = 0
        while len(unclustered_seqs) > 0:
            add_cluster(last_cluster_id)
            last_cluster_id += 1

        for val in self.id_clusters.values():
            print ':'.join(val)

        if reco_info is not None:
            true_cluster_list, inferred_cluster_list = [], []
            for clid, uids in self.id_clusters.items():
                for uid in uids:
                    try:
                        true_cluster_list.append(reco_info[uid]['reco_id'])
                    except KeyError:
                        true_cluster_list.append(reco_info[int(uid)]['reco_id'])
                    inferred_cluster_list.append(clid)
            print '       true clusters %d' % len(set(true_cluster_list))
            print '   inferred clusters %d' % len(set(inferred_cluster_list))
            print '         adjusted mi %.2f' % adjusted_mutual_info_score(true_cluster_list, inferred_cluster_list)

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
        for name, cluster_id in self.query_clusters.iteritems():
            if cluster_id == second_cluster_id:
                self.query_clusters[name] = first_cluster_id

        if second_cluster_id in self.cluster_ids:
            self.cluster_ids.remove(second_cluster_id)
        else:
            print 'oh, man, something\'s wrong'
            print 'uniqe_id,reco_id'
            for name, cluster_id in self.query_clusters.iteritems():
                print '%s,%d' % (name, cluster_id)
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
            print 'ERROR nan passed for %d %d (dbg %s)' %(query_name, second_query_name, dbg_str_list)
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
