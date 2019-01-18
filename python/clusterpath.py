import dendropy
import os
import sys
import math
import csv

import utils
import treeutils

# ----------------------------------------------------------------------------------------
class ClusterPath(object):
    def __init__(self, initial_path_index=0, seed_unique_id=None, partition=None, fname=None, partition_lines=None):  # <partition> is a fully-formed partition, while <partition_lines> is straight from reading a file (perhaps could combine them, but I don't want to think through it now)
        # could probably remove path index since there's very little chance of doing smc in the future, but the path-merging code in glomerator was _very_ difficult to write, so I'm reluctant to nuke it
        self.initial_path_index = initial_path_index  # NOTE this is set to None if it's nonsensical, e.g. if we're merging several paths with different indices

        # NOTE make *damn* sure if you add another list here that you also take care of it in remove_first_partition()
        self.partitions = []  # it would of course be damn nice to glomph these into a class at some point
        self.logprobs = []
        self.n_procs = []
        self.ccfs = []  # pair of floats (not just a float) for each partition
        self.logweights = []
        self.n_lists = 5  # just to make sure you don't forget

        self.best_minus = 30.  # rewind by this many units of log likelihood when merging separate processes (note that this should really depend on the number of sequences)
        self.i_best, self.i_best_minus_x = None, None
        self.we_have_a_ccf = False  # did we read in at least one adj mi value from a file?
        self.trees = None  # list of trees corresponding to clusters in most likely partition

        self.seed_unique_id = seed_unique_id

        if partition is not None:
            self.add_partition(partition, logprob=0., n_procs=1)
        elif fname is not None:
            self.readfile(fname)
        elif partition_lines is not None:
            self.readlines(partition_lines)

    # ----------------------------------------------------------------------------------------
    def n_seqs(self, ip=0):  # number of sequences in each partition (shouldn't depend on which partition, but I'm not sure that I absolutely forbid the number to change)
        if len(self.partitions) == 0:
            return 0
        return len([u for c in self.partitions[ip] for u in c])

    # ----------------------------------------------------------------------------------------
    def get_headers(self, is_data):
        headers = ['logprob', 'n_clusters', 'n_procs', 'partition']
        if not is_data:
            headers += ['n_true_clusters', 'ccf_under', 'ccf_over']
            # headers += ['bad_clusters']  # uncomment to also write the clusters that aren't perfect
        if self.seed_unique_id is not None:
            headers += ['seed_unique_id', ]
        return headers

    # ----------------------------------------------------------------------------------------
    def update_best_minus_x_partition(self):
        if self.i_best is None:  # only happens if we just removed the only partition
            self.i_best_minus_x = None
            return
        if math.isinf(self.logprobs[self.i_best]):  # if logprob is infinite, set best and best minus x to the latest one
            self.i_best_minus_x = self.i_best
            return
        for ip in range(len(self.partitions)):  # they should be in order of increasing logprob (at least within a give number of procs)
            if self.n_procs[ip] != self.n_procs[self.i_best]:  # only consider partitions with the same number of procs (e.g. if best partition is for 1 proc, we want the best-minus-x to also be for 1 proc)
                continue
            if self.logprobs[ip] > self.logprobs[self.i_best] - self.best_minus:  # pick the first one that is above threshold
                self.i_best_minus_x = ip
                break

    # ----------------------------------------------------------------------------------------
    def add_partition(self, partition, logprob, n_procs, logweight=None, ccfs=None):
        if ccfs is None:
            ccfs = [None, None]
        # NOTE you typically want to allow duplicate (in terms of log prob) partitions, since they can have different n procs
        self.partitions.append(partition)  # NOTE not deep copied
        self.logprobs.append(logprob)
        self.n_procs.append(n_procs)
        self.logweights.append(logweight)
        if len(ccfs) != 2:
            raise Exception('tried to add partition with ccfs of length %d (%s)' % (len(ccfs), ccfs))
        self.ccfs.append(ccfs)
        if ccfs.count(None) != len(ccfs):
            self.we_have_a_ccf = True
        # set this as the best partition if 1) we haven't set i_best yet, 2) this partition is more likely than i_best, or 3) i_best is set for a partition with a larger number of procs or 4) logprob is infinite (i.e. we didn't calculate the full partitions logprob))
        if self.i_best is None or logprob > self.logprobs[self.i_best] or n_procs < self.n_procs[self.i_best] or math.isinf(logprob):
            self.i_best = len(self.partitions) - 1
        self.update_best_minus_x_partition()
        self.trees = None  # they'll be out of date if the best partition changed (so I guess in principle we could only do this if self.i_best changes, but this seems tidier)

    # ----------------------------------------------------------------------------------------
    def remove_partition(self, ip_to_remove):  # NOTE doesn't update self.we_have_a_ccf, but it probably won't change, right?
        self.partitions.pop(ip_to_remove)
        self.logprobs.pop(ip_to_remove)
        self.n_procs.pop(ip_to_remove)
        self.ccfs.pop(ip_to_remove)
        self.logweights.pop(ip_to_remove)
        assert self.n_lists == 5  # make sure we didn't add another list and forget to put it in here

        self.i_best = None  # indices are all off now, anyway, since we just removed a partition, so may as well start from scratch
        for ip in range(len(self.partitions)):  # update self.i_best
            if self.i_best is None or self.logprobs[ip] > self.logprobs[self.i_best] or self.n_procs[ip] < self.n_procs[self.i_best] or math.isinf(self.logprobs[ip]):  # NOTE duplicates code in add_partition()
                self.i_best = ip
        self.update_best_minus_x_partition()
        self.trees = None  # they'll be out of date if the best partition changed (so I guess in principle we could only do this if self.i_best changes, but this seems tidier)

    # ----------------------------------------------------------------------------------------
    def readfile(self, fname):
        if fname is None:
            raise Exception('can\'t read NoneType partition file')
        if os.stat(fname).st_size == 0:
            raise Exception('partition file %s has size zero' % fname)

        if utils.getsuffix(fname) == '.csv':
            with open(fname, 'r') as infile:
                reader = csv.DictReader(infile)
                if 'partition' not in reader.fieldnames:
                    raise Exception('\'partition\' not among headers in %s, maybe this isn\'t a partition file?' % fname)
                lines = [line for line in reader]  # not sure that I really need this step
            self.readlines(lines, process_csv=True)
        elif utils.getsuffix(fname) == '.yaml':
            utils.read_yaml_output(fname, cpath=self)
        else:
            raise Exception('unhandled annotation file suffix %s' % outfname)

    # ----------------------------------------------------------------------------------------
    def readlines(self, lines, process_csv=False):
        for line in lines:
            if 'path_index' in line and int(line['path_index']) != self.initial_path_index:  # if <lines> contains more than one path_index, that means they represent more than one path, so you need to use glomerator, not just one ClusterPath
                raise Exception('path index in lines %d doesn\'t match my initial path index %d' % (int(line['path_index']), self.initial_path_index))
            if 'seed_unique_id' in line and line['seed_unique_id'] != '':
                if self.seed_unique_id is None:
                    self.seed_unique_id = line['seed_unique_id']
                if line['seed_unique_id'] != self.seed_unique_id:
                    print '%s seed uids for each line not all the same %s %s' % (utils.color('yellow', 'warning'), line['seed_unique_id'], self.seed_unique_id)

            if process_csv:
                line['partition'] = [cluster_str.split(':') for cluster_str in line['partition'].split(';')]

            ccfs = [None, None]
            if 'ccf_under' in line and 'ccf_over' in line:  # I don't know what I want to do if there's one but not the other, but it shouldn't be possible
                if line['ccf_under'] != '' and line['ccf_over'] != '':
                    ccfs = [float(line['ccf_under']), float(line['ccf_over'])]
                self.we_have_a_ccf = True

            self.add_partition(line['partition'], float(line['logprob']), int(line.get('n_procs', 1)), logweight=float(line.get('logweight', 0)), ccfs=ccfs)

    # ----------------------------------------------------------------------------------------
    def calculate_missing_values(self, reco_info, only_ip=None):
        for ip in range(len(self.partitions)):
            if only_ip is not None and ip != only_ip:
                continue

            if self.ccfs[ip][0] is not None and self.ccfs[ip][1] is not None:  # already have them
                continue

            true_partition = utils.get_true_partition(reco_info, ids=[uid for cluster in self.partitions[ip] for uid in cluster])  # can't use the true partition we might already have in the calling function, since we need this to only have uids from this partition (i.e. this isn't really the true partition)
            new_vals = utils.new_ccfs_that_need_better_names(self.partitions[ip], true_partition, reco_info, seed_unique_id=self.seed_unique_id)
            if None not in new_vals:  # if the function finds messed up partitions, it returns None, None (at this point, this seems to just happens when a uid was found in multiple clusters, which happens for earlier partitions (n_procs > 1) when seed_unique_id is set, since we pass seed_unique_id to all the subprocs)
                self.ccfs[ip] = new_vals
                self.we_have_a_ccf = True

    # ----------------------------------------------------------------------------------------
    def get_ccf_str(self, ip):
        ccf_str = ''
        if self.we_have_a_ccf:
            ccf_str_list = [('%5s' % '-') if ccf is None else ('%5.2f' % ccf) for ccf in self.ccfs[ip]]
            ccf_str = ' %s ' % ' '.join(ccf_str_list)
            # if self.ccfs[ip][0] is None and self.ccfs[ip][1] is None:
            #     ccf_str = '   -  -    '
            # else:
            #     ccf_str = ' %5.2f %5.2f    ' % tuple(self.ccfs[ip])
        else:
            ccf_str = '   -  -    '

        return ccf_str

    # ----------------------------------------------------------------------------------------
    def print_partition(self, ip, reco_info=None, extrastr='', abbreviate=True, highlight_cluster_indices=None, right_extrastr=''):
        #  NOTE it's nicer to *not* sort by cluster size here, since preserving the order tends to frequently make it obvious which clusters are merging as your eye scans downwards through the output
        if ip > 0:  # delta between this logprob and the previous one
            delta_str = '%.1f' % (self.logprobs[ip] - self.logprobs[ip-1])
        else:
            delta_str = ''
        print '      %s  %-12.2f%-7s   %-5d  %4d' % (extrastr, self.logprobs[ip], delta_str, len(self.partitions[ip]), self.n_procs[ip]),

        print '    ' + self.get_ccf_str(ip),

        # clusters
        sorted_clusters = sorted(self.partitions[ip], key=lambda c: len(c), reverse=True)
        for iclust in range(len(sorted_clusters)):
            cluster = sorted_clusters[iclust]
            if abbreviate:
                cluster_str = ':'.join(['o' if len(uid) > 3 else uid for uid in cluster])
            else:
                # cluster_str = ':'.join(sorted([str(uid) for uid in cluster]))
                cluster_str = ':'.join([str(uid) for uid in cluster])

            if reco_info is not None and not utils.from_same_event(reco_info, cluster):
                cluster_str = utils.color('red', cluster_str)

            if self.seed_unique_id is not None and self.seed_unique_id in cluster:
                cluster_str = utils.color('reverse_video', cluster_str)

            if highlight_cluster_indices is not None and iclust in highlight_cluster_indices:
                cluster_str = utils.color('red', cluster_str)
            
            if abbreviate:
                print ' %s' % cluster_str,
            else:
                print '   %s' % cluster_str,
        print '%s' % right_extrastr,
        print ''

    # ----------------------------------------------------------------------------------------
    def print_partitions(self, reco_info=None, extrastr='', abbreviate=True, print_header=True, n_to_print=None, calc_missing_values='none', highlight_cluster_indices=None):
        assert calc_missing_values in ['none', 'all', 'best']
        if reco_info is not None and calc_missing_values == 'all':
            self.calculate_missing_values(reco_info)

        if print_header:
            print '    %7s %10s   %-7s %5s  %4s' % ('', 'logprob', 'delta', 'clusters', 'n_procs'),
            if reco_info is not None or self.we_have_a_ccf:
                print ' %5s %5s' % ('purity', 'completeness'),
            print ''

        for ip in self.get_surrounding_partitions(n_partitions=n_to_print):
            if reco_info is not None and calc_missing_values == 'best' and ip == self.i_best:
                self.calculate_missing_values(reco_info, only_ip=ip)
            mark = '      '
            if ip == self.i_best:
                mark = 'best  '
            if ip == self.i_best_minus_x:
                mark = mark[:-2] + '* '
            if mark.count(' ') < len(mark):
                mark = utils.color('yellow', mark)
            right_extrastr = '' if self.n_seqs() < 400 else mark  # if line is going to be really long, put the yellow stuff also on the right side
            self.print_partition(ip, reco_info, extrastr=mark+extrastr, abbreviate=abbreviate, highlight_cluster_indices=highlight_cluster_indices, right_extrastr=right_extrastr)

    # ----------------------------------------------------------------------------------------
    def get_surrounding_partitions(self, n_partitions):
        """ return a list of partition indices centered on <self.i_best> of length <n_partitions> """
        if n_partitions is None:  # print all partitions
            ilist = range(len(self.partitions))
        else:  # print the specified number surrounding the maximum logprob
            if n_partitions < 0 or n_partitions >= len(self.partitions):
                n_partitions = len(self.partitions)
            ilist = [self.i_best, ]
            while len(ilist) < n_partitions:  # add partition numbers before and after <i_best> until we get to <n_partitions>
                if ilist[0] > 0:  # stop adding them beforehand if we've hit the first partition
                    ilist.insert(0, ilist[0] - 1)
                if len(ilist) < n_partitions and ilist[-1] < len(self.partitions) - 1:  # don't add them afterward if we already have enough, or if we're already at the end
                    ilist.append(ilist[-1] + 1)

        return ilist

    # ----------------------------------------------------------------------------------------
    def get_parent_clusters(self, ipart):
        """ Return the parent clusters that were merged to form the <ipart>th partition. """
        if ipart == 0:
            raise Exception('get_parent_clusters got ipart of zero... that don\'t make no sense yo')
        if len(self.partitions[ipart - 1]) <= len(self.partitions[ipart]):
            return None  # this step isn't a merging step -- it's a synthetic rewinding step due to multiple processes

        parents = []
        for cluster in self.partitions[ipart - 1]:  # find all clusters in the previous partition that aren't in the current one
            if cluster not in self.partitions[ipart]:
                parents.append(cluster)
        assert len(parents) == 2  # there should've been two -- those're the two that were merged to form the new cluster
        return parents

    # ----------------------------------------------------------------------------------------
    def set_synthetic_logweight_history(self, reco_info):
        # not sure if it's still relevant, but note here said: "switch clusterpath.cc back to using these"
        def potential_n_parents(partition):
            combifactor = 0
            for cluster in partition:
                n_k = len(cluster)
                combifactor += pow(2, n_k - 1) - 1
            if combifactor == 0:
                combifactor = 1
            return combifactor

        for ip in range(len(self.partitions)):
            if ip == 0:
                last_logweight = 0.
            else:
                last_logweight = self.logweights[ip-1]
            this_logweight = last_logweight + math.log(1. / potential_n_parents(self.partitions[ip]))
            self.logweights[ip] = this_logweight

    # ----------------------------------------------------------------------------------------
    def write(self, outfname, is_data, reco_info=None, true_partition=None, n_to_write=None, calc_missing_values='none', partition_lines=None):
        if utils.getsuffix(outfname) != '.csv':
            raise Exception('unhandled file extension %s' % outfname)
        if partition_lines is None:
            partition_lines = self.get_partition_lines(is_data, reco_info=reco_info, true_partition=true_partition, n_to_write=n_to_write, calc_missing_values=calc_missing_values)
        with open(outfname, 'w') as outfile:
            writer = csv.DictWriter(outfile, self.get_headers(is_data))
            writer.writeheader()
            for row in partition_lines:
                row['partition'] = ';'.join([':'.join(cluster) for cluster in row['partition']])
                if 'bad_clusters' in row:
                    row['bad_clusters'] = ';'.join(row['bad_clusters'])
                writer.writerow(row)

    # ----------------------------------------------------------------------------------------
    def get_partition_lines(self, is_data, reco_info=None, true_partition=None, n_to_write=None, calc_missing_values='none', path_index=None):  # we use this (instead of .write()) if we're writing a yaml file
        if not is_data:
            assert reco_info is not None and true_partition is not None  # we could get the true_partition if we have reco_info, but easier to just make the caller do it
        assert calc_missing_values in ['none', 'all', 'best']
        if reco_info is not None and calc_missing_values == 'all':
            self.calculate_missing_values(reco_info)

        headers = self.get_headers(is_data)
        lines = []
        for ipart in self.get_surrounding_partitions(n_partitions=n_to_write):
            part = self.partitions[ipart]

            row = {'logprob' : self.logprobs[ipart],
                   'n_clusters' : len(part),
                   'n_procs' : self.n_procs[ipart],
                   'partition' : part}
            if 'ccf_under' in headers:
                if reco_info is not None and calc_missing_values == 'best' and ipart == self.i_best:
                    self.calculate_missing_values(reco_info, only_ip=ipart)
                if self.ccfs[ipart][0] is not None and self.ccfs[ipart][1] is not None:
                    row['ccf_under'], row['ccf_over'] = self.ccfs[ipart]  # for now assume we calculated the ccfs if we did adj mi
            if 'n_true_clusters' in headers:
                row['n_true_clusters'] = len(true_partition)
            if 'bad_clusters' in headers:
                row['bad_clusters'] = self.get_bad_clusters(part, reco_info, true_partition)
            if 'path_index' in headers:
                row['path_index'] = path_index
                row['logweight'] = self.logweights[ipart]
            if 'seed_unique_id' in headers:
                row['seed_unique_id'] = self.seed_unique_id

            lines.append(row)

        return lines

    # ----------------------------------------------------------------------------------------
    def get_bad_clusters(self, partition, reco_info, true_partition):
        bad_clusters = []  # inferred clusters that aren't really all from the same event
        for ic in range(len(partition)):
            same_event = utils.from_same_event(reco_info, partition[ic])  # are all the sequences from the same event?
            entire_cluster = True  # ... and if so, are they the entire true cluster?
            if same_event:
                reco_id = reco_info[partition[ic][0]]['reco_id']  # they've all got the same reco_id then, so pick an aribtrary one
                true_clusters = [cluster for cluster in true_partition if reco_info[cluster[0]]['reco_id'] == reco_id]  # NOTE I think this doesn't work right with shm indels in the cdr3
                assert len(true_clusters) == 1
                true_cluster = true_clusters[0]
                for uid in true_cluster:
                    if uid not in partition[ic]:
                        entire_cluster = False
                        break
            else:
                entire_cluster = False
            if not same_event or not entire_cluster:
                bad_clusters.append(':'.join(partition[ic]))

        if len(bad_clusters) > 25:
            bad_clusters = ['too', 'long']

        return bad_clusters

    # ----------------------------------------------------------------------------------------
    def write_presto_partitions(self, outfname, input_info):
        with open(outfname, 'w') as outfile:
            iclust = 0
            for cluster in self.partitions[self.i_best]:
                for uid in cluster:
                    assert len(input_info[uid]['seqs']) == 1
                    outfile.write('>%s|CLONE=%d\n%s\n' % (uid, iclust, input_info[uid]['seqs'][0]))
                iclust += 1

    # ----------------------------------------------------------------------------------------
    # make tree for the single cluster in the last partition of <partitions> (which is in general *not* the last partition in self.partitions, since the calling function stops at self.i_best)
    def make_single_tree(self, partitions, annotations, uid_set, get_fasttrees=False, n_max_cons_seqs=10, debug=False):
        # NOTE don't call this externally -- if you want a single tree, call make_trees() with <i_only_cluster> set
        def getline(uidstr, uid_set=None):
            if uidstr in annotations:  # if we have this exact annotation
                return annotations[uidstr]
            else:
                if uid_set is None:
                    uid_set = set(uidstr.split(':'))  # should only get called if it's a singleton
                # note that for internal nodes in a fasttree-derived subtree, the uids will be out of order compared the the annotation keys
                for line in annotations.values():  # we may actually have the annotation for every subcluster (e.g. if --calculate-alternative-naive-seqs was set), but in case we don't, this is fine
                    if len(uid_set & set(line['unique_ids'])) > 0:  # just take the first one with any overlap. Yeah, it's not necessarily the best, but its naive sequence probably isn't that different, and for just getting the fasttree it reeeeeeaaaallly doesn't matter
                        return line
            raise Exception('couldn\'t find uid %s in annotations' % uid)
        def getseq(uid):
            line = getline(uid)
            return line['seqs'][line['unique_ids'].index(uid)]
        def lget(uid_list):
            return ':'.join(uid_list)

        default_edge_length = 999999  # it's nice to have the edges all set to something that's numeric (so the trees print), but also obvious wrong, if we forget to set somebody
        assert len(partitions[-1]) == 1
        root_label = lget(partitions[-1][0])  # we want the order of the uids in the label to correspond to the order in self.partitions
        tns = dendropy.TaxonNamespace([root_label])
        root_node = dendropy.Node(taxon=tns.get_taxon(root_label))
        root_node.uids = uid_set  # each node keeps track of the uids of its children
        dtree = dendropy.Tree(taxon_namespace=tns, seed_node=root_node)
        if debug:
            print '    starting tree with %d leaves' % len(uid_set)
        for ipart in reversed(range(len(partitions) - 1)):  # dendropy seems to only have fcns to build a tree from the root downward, so we loop starting with the last partition (- 1 is because the last partition is guaranteed to be just one cluster)
            for lnode in dtree.leaf_node_iter():  # look for leaf nodes that contain uids from two clusters in this partition, and add those as children
                tclusts = [c for c in partitions[ipart] if len(set(c) & lnode.uids) > 0]
                if len(tclusts) < 2:
                    continue
                for tclust in tclusts:
                    ttaxon = dendropy.Taxon(lget(tclust))
                    tns.add_taxon(ttaxon)
                    child = lnode.new_child(taxon=ttaxon, edge_length=default_edge_length)
                    child.uids = set(tclust)
                if debug:
                    print '      ipart %d' % ipart
                    print '        split node: %d --> %s      %s --> %s' % (len(lnode.uids), ' '.join([str(len(tc)) for tc in tclusts]), lnode.taxon.label, ' '.join([c.taxon.label for c in lnode.child_node_iter()]))

        # split existing leaves, which are probably not singletons (they're probably from the initial naive sequence collapse step) into subtrees such that each leaf is a singleton
        for lnode in dtree.leaf_node_iter():
            if len(lnode.uids) == 1:
                continue
            if get_fasttrees and len(lnode.uids) > 2:
                seqfos = [{'name' : uid, 'seq' : getseq(uid)} for uid in lnode.taxon.label.split(':')]  # may as well add them in the right order, although I don't think it matters
                subtree = treeutils.get_fasttree_tree(seqfos, getline(lnode.taxon.label, uid_set=lnode.uids)['naive_seq'], suppress_internal_node_taxa=True)  # note that the fasttree distances get ignored below (no idea if they'd be better than what we set down there, but they probably wouldn't be consistent, so I'd rather ignore them)
                for tmpnode in subtree.postorder_node_iter():
                    if tmpnode.is_leaf():
                        tmpnode.uids = set([tmpnode.taxon.label])
                    else:
                        tmpnode.uids = set([uid for c in tmpnode.child_node_iter() for uid in c.uids])
                        ttaxon = dendropy.Taxon(lget(tmpnode.uids))
                        subtree.taxon_namespace.add_taxon(ttaxon)
                        tmpnode.taxon = ttaxon  # ...and use the string of leaf nodes, even though they'll be in the wrong order (I think these get ignored when I call label_nodes() below, but it's still tidier to have them right in the meantime, and anyway since I'm suppressing internal taxa I think I need to set them to something)

                if debug:
                    print '   adding subtree with %d leaves from fastree at leaf node %s' % (len(seqfos), lnode.taxon.label)
                    print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=subtree))
                dtree.taxon_namespace.add_taxa(subtree.taxon_namespace)
                lnode.add_child(subtree.seed_node)
                assert len(lnode.child_edges()) == 1  # we're iterating over leaves, so this should always be true
                lnode.child_edges()[0].collapse()
            else:  # just add a star subtree
                for uid in lnode.taxon.label.split(':'):  # may as well add them in the right order, although I don't think it matters
                    ttaxon = dendropy.Taxon(uid)
                    tns.add_taxon(ttaxon)
                    child = lnode.new_child(taxon=ttaxon, edge_length=default_edge_length)
                    child.uids = set([uid])
                if debug:
                    print '      added %d singleton children for %s' % (len(lnode.uids), lnode.taxon.label)

        # in order to set edge lengths, we need node sequences, so first set leaf node seqs
        for lnode in dtree.leaf_node_iter():
            assert len(lnode.uids) == 1
            lnode.seq = getseq(lnode.taxon.label)
            lnode.n_descendent_leaves = 1  # keep track of how many leaf nodes contributed to each node's consensus sequence (these are leaves, so it's trivally 1). This is less accurate than keeping track of all the sequences, but also faster

        # then set internal node seqs as the consensus of their children, and set the distance as hamming distance to child seqs
        if debug:
            print '    adding edge lengths either from fasttree %s or cons seq %s' % (utils.color('blue', 'x'), utils.color('red', 'x'))
        min_edge_length = None  # setting this is nice for better debug viewing
        for node in dtree.postorder_internal_node_iter():  # includes root node
            child_cons_seq_counts = [c.n_descendent_leaves for c in node.child_node_iter()]
            total_descendent_leaves = sum(child_cons_seq_counts)
            if total_descendent_leaves > n_max_cons_seqs:  # if there's tons of descendent leaves, we don't want to pass them all to the consensus fcn since it's slow, so we choose them in proportion to their actual proportions, but scaled down to <n_max_cons_seqs>
                child_cons_seq_counts = [int(n_max_cons_seqs * csc / float(total_descendent_leaves)) for csc in child_cons_seq_counts]
                child_cons_seq_counts = [max(1, csc) for csc in child_cons_seq_counts]  # don't eliminate any sequences entirely (this makes the proportions less accurate (in some cases), but is the easy way to handle the case where there's a ton of singleton children
            if debug:
                print '  %s' % utils.color('green', node.taxon.label)
                csc_str = '  (reduced: %s)' % ' '.join([str(csc) for csc in child_cons_seq_counts]) if total_descendent_leaves > n_max_cons_seqs else ''
                print '      desc leaves per child: %s%s' % (' '.join(str(c.n_descendent_leaves) for c in node.child_node_iter()), csc_str)
            child_seqfos = [{'name' : cn.taxon.label + '-leaf-' + str(il), 'seq' : cn.seq} for cn, count in zip(node.child_node_iter(), child_cons_seq_counts) for il in range(count)]
            node.seq = utils.cons_seq(0.01, aligned_seqfos=child_seqfos, tie_resolver_seq=getline(root_label)['naive_seq'])  #, debug=debug)  # the consensus has an N at every position where the constituent sequences gave a tie. But Ns screw up the distances (especially because once we *get* an N, we can't get rid of it and it's propagated all the way up the tree), and in almost all cases the correct choice should be the naive base, so we use that
            node.n_descendent_leaves = total_descendent_leaves
            for edge in node.child_edge_iter():
                from_fasttree = False
                if edge.length == default_edge_length:  # otherwise it was set by fasttree, and it's probably better than what we'd get from this (it'd be nice to skip the cons seq stuff for the whole fasttree subtree, but then we don't have the cons seqs we need for later)
                    edge.length = utils.hamming_distance(edge.head_node.seq, node.seq) / float(len(node.seq))
                else:
                    from_fasttree = True
                if min_edge_length is not None:
                    edge.length = max(min_edge_length, edge.length)
                if debug:
                    print '       %6.3f   %s  %s' % (edge.length, utils.color('blue' if from_fasttree else 'red', 'x'), edge.head_node.taxon.label)

        if debug:
            print '        naive seq %s' % getline(root_label)['naive_seq'] # NOTE might be worthwhile to add an edge connecting seed node and the actual naive sequence (i.e. for cases where our approximate naive is off)
            print '    root cons seq %s' % utils.color_mutants(getline(root_label)['naive_seq'], dtree.seed_node.seq)

        for node in dtree.preorder_node_iter():
            del node.uids
            del node.seq
            del node.n_descendent_leaves

        treeutils.label_nodes(dtree, ignore_existing_internal_node_labels=True, ignore_existing_internal_taxon_labels=True, debug=debug)
        dtree.update_bipartitions()  # probably don't really need this
        if debug:
            print treeutils.utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree, width=250))

        return dtree

    # ----------------------------------------------------------------------------------------
    def make_trees(self, annotations, i_only_cluster=None, get_fasttrees=False, debug=False):  # makes a tree for each cluster in the most likely (not final) partition
        # i_only_cluster: only make the tree corresponding to the <i_only_cluster>th cluster in the best partition
        if self.i_best is None:
            return

        if debug:
            print '  making %d tree%s over %d partitions' % (len(self.partitions[self.i_best]) if i_only_cluster is None else 1, 's' if i_only_cluster is None else '', self.i_best + 1)
        if self.trees is None:
            self.trees = [None for _ in self.partitions[self.i_best]]
        else:
            assert len(self.trees) == len(self.partitions[self.i_best])  # presumably because we were already called with <i_only_cluster> set for a different cluster
        for i_cluster in range(len(self.partitions[self.i_best])):
            if i_only_cluster is not None and i_cluster != i_only_cluster:
                continue
            final_cluster = self.partitions[self.i_best][i_cluster]
            uid_set = set(final_cluster)  # usually the set() isn't doing anything, but sometimes I think we have uids duplicated between clusters, e.g. I think when seed partitioning (or even within a cluster, because order matters within a cluster because of bcrham caching)
            sub_partitions = [[] for _ in range(self.i_best + 1)]  # new list of partitions, but only including clusters that overlap with <final_cluster>
            for ipart in range(self.i_best + 1):
                for tmpclust in self.partitions[ipart]:
                    if len(set(tmpclust) & uid_set) == 0:
                        continue
                    sub_partitions[ipart].append(tmpclust)  # note that many of these adjacent sub-partitions can be identical, if the merges happened between clusters that correspond to a different final cluster
            self.trees[i_cluster] = self.make_single_tree(sub_partitions, annotations, uid_set, get_fasttrees=get_fasttrees, debug=debug)
