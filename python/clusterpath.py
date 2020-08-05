import dendropy
import os
import sys
import math
import csv
import copy
import collections
import itertools
import numpy

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
    def find_iparts_for_cluster(self, cluster):  # get index of partitions in which a list of uids (i.e. a cluster) appears
        return [ip for ip in range(len(self.partitions)) if cluster in self.partitions[ip]]  # NOTE just returns zero-length list if it isn't there

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
                    raise Exception('\'partition\' not among headers in %s, maybe this isn\'t a partition file? (if you\'re running \'view-output\' on a deprecated csv output file, you may need to run \'view-annotations\' instead, to tell it that this is an annotation file rather than a partition file)' % fname)
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

            # NOTE this isn't really the true partition
            true_partition = utils.get_partition_from_reco_info(reco_info, ids=[uid for cluster in self.partitions[ip] for uid in cluster])  # can't use the true partition we might already have in the calling function, since we need this to only have uids from this partition (i.e. this isn't really the true partition)
            new_vals = utils.new_ccfs_that_need_better_names(self.partitions[ip], true_partition, reco_info=reco_info, seed_unique_id=self.seed_unique_id)
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
    def print_partition(self, ip, reco_info=None, extrastr='', abbreviate=True, highlight_cluster_indices=None, print_partition_indices=False, right_extrastr='', sort_by_size=True):  # NOTE <highlight_cluster_indices> and <print_partition_indices> are quite different despite sounding similar, but I can't think of something else to call the latter that makes more sense
        if ip > 0:  # delta between this logprob and the previous one
            delta_str = '%.1f' % (self.logprobs[ip] - self.logprobs[ip-1])
        else:
            delta_str = ''
        print '      %s  %-12.2f%-7s   %s%-5d  %4d' % (extrastr, self.logprobs[ip], delta_str, ('%-5d  ' % ip) if print_partition_indices else '', len(self.partitions[ip]), self.n_procs[ip]),

        print '    ' + self.get_ccf_str(ip),

        sorted_clusters = self.partitions[ip]
        if sort_by_size:  # it's often nicer to *not* sort by cluster size here, since preserving the order frequently makes it obvious which clusters are merging as your eye scans downward through the output
            sorted_clusters = sorted(sorted_clusters, key=lambda c: len(c), reverse=True)
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
    def print_partitions(self, reco_info=None, extrastr='', abbreviate=True, print_header=True, n_to_print=None, calc_missing_values='none', highlight_cluster_indices=None, print_partition_indices=False, ipart_center=None, sort_by_size=True):
        assert calc_missing_values in ['none', 'all', 'best']
        if reco_info is not None and calc_missing_values == 'all':
            self.calculate_missing_values(reco_info)

        if print_header:
            print '    %s%7s %10s   %-7s %s%5s  %4s' % (' '*utils.len_excluding_colors(extrastr), '', 'logprob', 'delta', 'index  ' if print_partition_indices else '', 'clusters', 'n_procs'),
            if reco_info is not None or self.we_have_a_ccf:
                print ' %5s %5s' % ('purity', 'completeness'),
            print ''

        for ip in self.get_surrounding_partitions(n_to_print, i_center=ipart_center):
            if reco_info is not None and calc_missing_values == 'best' and ip == self.i_best:
                self.calculate_missing_values(reco_info, only_ip=ip)
            mark = '      '
            if ip == self.i_best:
                mark = 'best  '
            if ip == self.i_best_minus_x:
                mark = mark[:-2] + '* '
            if mark.count(' ') < len(mark):
                mark = utils.color('yellow', mark)
            right_extrastr = '' if self.n_seqs() < 200 else mark  # if line is going to be really long, put the yellow stuff also on the right side
            self.print_partition(ip, reco_info, extrastr=extrastr+mark, abbreviate=abbreviate, highlight_cluster_indices=highlight_cluster_indices, print_partition_indices=print_partition_indices, right_extrastr=right_extrastr, sort_by_size=sort_by_size)

    # ----------------------------------------------------------------------------------------
    def get_surrounding_partitions(self, n_partitions, i_center=None):
        """ return a list of partition indices centered on <self.i_best> of length <n_partitions> """
        if i_center is None:
            i_center = self.i_best
        if n_partitions is None:  # print all partitions
            ilist = range(len(self.partitions))
        else:  # print the specified number surrounding (by default) the maximum logprob
            if n_partitions < 0 or n_partitions >= len(self.partitions):
                n_partitions = len(self.partitions)
            ilist = [i_center, ]
            while len(ilist) < n_partitions:  # add partition numbers before and after <i_center> until we get to <n_partitions>
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
        for ipart in self.get_surrounding_partitions(n_to_write):
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
    def seed_cluster(self, ipart=None):
        if self.seed_unique_id is None:
            print '  %s self.seed_unique_id isn\'t set' % utils.color('red', 'error')
            return  # I don't like returning different types of things, but this lets the calling code handle things
        seed_clusters = [c for c in self.partitions[ipart if ipart is not None else self.i_best] if self.seed_unique_id in c]
        if len(seed_clusters) == 0:
            print '  %s couldn\'t find cluster containing seed %s' % (utils.color('red', 'error'), self.seed_unique_id)
        elif len(seed_clusters) > 1:
            print '  %s more than one cluster containing seed %s (this can be expected if this is the non-best partition)' % (utils.color('red', 'error'), self.seed_unique_id)
        else:
            return seed_clusters[0]

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
        print '   writing presto partition %s' % outfname
        assert utils.getsuffix(outfname) in ['.fa', '.fasta']  # already checked in processargs.py
        with open(outfname, 'w') as outfile:
            iclust = 0
            for cluster in self.partitions[self.i_best]:
                for uid in cluster:
                    assert len(input_info[uid]['seqs']) == 1
                    outfile.write('>%s|CLONE=%d\n%s\n' % (uid, iclust, input_info[uid]['seqs'][0]))
                iclust += 1

    # ----------------------------------------------------------------------------------------
    # make tree for the single cluster in the last partition of <partitions> (which is in general *not* the last partition in self.partitions, since the calling function stops at self.i_best)
    def make_single_tree(self, partitions, annotations, uid_set, naive_seq_name, get_fasttrees=False, n_max_cons_seqs=10, debug=False):
        # NOTE don't call this externally -- if you want a single tree, call make_trees() with <i_only_cluster> set
        def getline(uidstr, uid_set=None):
            if uidstr == naive_seq_name:
                assert False  # shouldn't actually happen (I think)
            elif uidstr in annotations:  # if we have this exact annotation
                return annotations[uidstr]
            else:
                if uid_set is None:
                    uid_set = set(uidstr.split(':'))  # should only get called if it's a singleton
                # note that for internal nodes in a fasttree-derived subtree, the uids will be out of order compared the the annotation keys
                for line in annotations.values():  # we may actually have the annotation for every subcluster (e.g. if --calculate-alternative-annotations was set), but in case we don't, this is fine
                    # print uid_set, len(uid_set & set(line['unique_ids']))
                    if len(uid_set & set(line['unique_ids'])) > 0:  # just take the first one with any overlap. Yeah, it's not necessarily the best, but its naive sequence probably isn't that different, and for just getting the fasttree it reeeeeeaaaallly doesn't matter
                        return line
            raise Exception('couldn\'t find uid %s in annotations' % uidstr)
        def getseq(uid):
            if uid == naive_seq_name:
                sorted_lines = sorted([l for l in annotations.values()], key=lambda l: len(l['unique_ids']), reverse=True)  # since we're making a tree, all the annotations are by definition clonal, so it doesn't really matter, but may as well get the naive sequence from the largest one
                return sorted_lines[0]['naive_seq']
            else:
                return utils.per_seq_val(getline(uid), 'seqs', uid)
        def lget(uid_list):
            return ':'.join(uid_list)

        # check for repeated uids (was only from seed uid, which shouldn't happen any more, but the code below throws an infinite loop if we do, so may as well be careful)
        for partition in partitions:
            if sum(len(c) for c in partition) > len(set(u for c in partition for u in c)):
                repeated_uids = [u for u, count in collections.Counter([u for c in partition for u in c]).items() if count > 1]
                raise Exception('found %d uid%s in more than one cluster (%s)' % (len(repeated_uids), utils.plural(len(repeated_uids)), ', '.join(repeated_uids)))

        default_edge_length = 999999  # it's nice to have the edges all set to something that's numeric (so the trees print), but also obvious wrong, if we forget to set somebody
        assert len(partitions[-1]) == 1
        root_label = lget(partitions[-1][0])  # we want the order of the uids in the label to correspond to the order in self.partitions
        tns = dendropy.TaxonNamespace([root_label])
        root_node = dendropy.Node(taxon=tns.get_taxon(root_label))
        root_node.uids = uid_set  # each node keeps track of the uids of its children
        dtree = dendropy.Tree(taxon_namespace=tns, seed_node=root_node, is_rooted=True)
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
                subtree = treeutils.get_fasttree_tree(seqfos, suppress_internal_node_taxa=True)  # note that the fasttree distances get ignored below (no idea if they'd be better than what we set down there, but they probably wouldn't be consistent, so I'd rather ignore them)
                for tmpnode in subtree.postorder_node_iter():
                    if tmpnode.is_leaf():
                        tmpnode.uids = set([tmpnode.taxon.label])
                    else:
                        tmpnode.uids = set([uid for c in tmpnode.child_node_iter() for uid in c.uids])
                        assert tmpnode.taxon is None  # make sure that we didn't do any modification to the fasttree tree to label the internal nodes
                        ttaxon = dendropy.Taxon(lget(tmpnode.uids))
                        subtree.taxon_namespace.add_taxon(ttaxon)
                        tmpnode.taxon = ttaxon  # ...and use the string of leaf nodes, even though they'll be in the wrong order (I think these get ignored when I call label_nodes() below, but it's still tidier to have them right in the meantime, and anyway since I'm suppressing internal taxa I think I need to set them to something)
                if debug:
                    print '   adding subtree with %d leaves from fastree at leaf node %s' % (len(seqfos), lnode.taxon.label)
                    print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=subtree, width=350, label_fcn=lambda x: '<combo>' if ':' in x else x))
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
        treeutils.collapse_zero_length_leaves(dtree, uid_set, debug=debug)
        # dtree.update_bipartitions(suppress_unifurcations=False)  # probably don't really need this UPDATE now i'm commenting it since it gets run by the zero length leaf collapse fcn
        if debug:
            print treeutils.utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree, width=250))

        # make sure we didn't lose any uids in the process of making the tree (this shouldn't happen any more, but it kept happening when i thought i'd fixed it before, so...)
        node_set = set(n.taxon.label for n in dtree.preorder_node_iter())
        missing_uids = uid_set - node_set
        if len(missing_uids) > 0:
            raise Exception('final clusterpath tree is missing %d input uids: %s' % (len(missing_uids), ' '.join(missing_uids)))

        return dtree

    # ----------------------------------------------------------------------------------------
    def deduplicate_seed_uid(self, debug=False):
        n_removed_list = [0 for _ in self.partitions]
        new_partitions = [[] for _ in self.partitions]
        for ipart in range(len(self.partitions)):
            found = False  # leave the seed uid in the first cluster you find it in, remove it from any others
            for cluster in self.partitions[ipart]:
                new_cluster = copy.deepcopy(cluster)
                if self.seed_unique_id in new_cluster:
                    if found:
                        while self.seed_unique_id in new_cluster:  # i think it can actually be in one cluster more than once
                            n_removed_list[ipart] += 1
                            new_cluster.remove(self.seed_unique_id)
                    else:
                        found = True
                if len(new_cluster) > 0:
                    new_partitions[ipart].append(new_cluster)
        if debug:
            print '  removed seed uid %s from %s clusters across %d partitions' % (self.seed_unique_id, ' '.join([str(n) for n in n_removed_list]), len(n_removed_list))

        return new_partitions

    # ----------------------------------------------------------------------------------------
    def get_sub_path(self, uid_set, partitions=None, annotations=None):  # get list of partitions (from 0 to self.i_best+1) restricted to clusters that overlap with uid_set
        if partitions is None:
            partitions = self.partitions
        sub_partitions = [[] for _ in range(self.i_best + 1)]  # new list of partitions, but only including clusters that overlap with uid_set
        sub_annotations = {}
        for ipart in range(self.i_best + 1):
            for tmpclust in partitions[ipart]:
                if len(set(tmpclust) & uid_set) == 0:
                    continue
                sub_partitions[ipart].append(tmpclust)  # note that many of these adjacent sub-partitions can be identical, if the merges happened between clusters that correspond to a different final cluster
                if annotations is not None and ':'.join(tmpclust) in annotations:
                    sub_annotations[':'.join(tmpclust)] = annotations[':'.join(tmpclust)]
        return sub_partitions, sub_annotations

    # ----------------------------------------------------------------------------------------
    def make_trees(self, annotations, i_only_cluster=None, get_fasttrees=False, naive_seq_name='XnaiveX', debug=False):  # makes a tree for each cluster in the most likely (not final) partition
        # i_only_cluster: only make the tree corresponding to the <i_only_cluster>th cluster in the best partition
        if self.i_best is None:
            return

        partitions = self.partitions
        if self.seed_unique_id is not None:
            partitions = self.deduplicate_seed_uid(debug=debug)

        if debug:
            print '  making %d tree%s over %d partitions' % (len(partitions[self.i_best]) if i_only_cluster is None else 1, 's' if i_only_cluster is None else '', self.i_best + 1)
        if self.trees is None:
            self.trees = [None for _ in partitions[self.i_best]]
        else:
            assert len(self.trees) == len(partitions[self.i_best])  # presumably because we were already called with <i_only_cluster> set for a different cluster
        for i_cluster in range(len(partitions[self.i_best])):
            if i_only_cluster is not None and i_cluster != i_only_cluster:
                continue
            uid_set = set(partitions[self.i_best][i_cluster])  # usually the set() isn't doing anything, but sometimes I think we have uids duplicated between clusters, e.g. I think when seed partitioning (or even within a cluster, because order matters within a cluster because of bcrham caching)
            sub_partitions, sub_annotations = self.get_sub_path(uid_set, partitions=partitions, annotations=annotations)
            self.trees[i_cluster] = self.make_single_tree(sub_partitions, sub_annotations, uid_set, naive_seq_name, get_fasttrees=get_fasttrees, debug=debug)

    # ----------------------------------------------------------------------------------------
    # cartoon explaining algorithm here https://github.com/psathyrella/partis/commit/ede140d76ff47383e0478c25fae8a9a9fa129afa#commitcomment-40981229
    def merge_light_chain(self, lcpath, heavy_annotations, light_annotations, h_ipart=None, l_ipart=None, check_partitions=False, debug=False):  # this assumes <self> is the heavy chain path, but it doesn't really need to be NOTE the clusters in the resulting partition generally have the uids in a totally different order to in either of the original partitions
        # ----------------------------------------------------------------------------------------
        def akey(klist):
            return ':'.join(klist)
        # ----------------------------------------------------------------------------------------
        def any_in_common(l1, l2):  # true if any uids in any cluster in l1 are found in any clusters in l2
            for tclust in l1:
                tset = set(tclust)
                if any(len(tset & set(tc)) > 0 for tc in l2):
                    return True
            return False
        # ----------------------------------------------------------------------------------------
        def common_clusters(tclust, tlist, return_indices=False):  # return all clusters in tlist that have uids in common with tclust
            tset = set(tclust)
            return [(i if return_indices else c) for i, c in enumerate(tlist) if len(set(c) & tset) > 0]
        # ----------------------------------------------------------------------------------------
        def is_clean_partition(putative_partition):  # make sure the list of clusters is actually disjoint
            return not any(len(set(c1) & set(c2)) > 0 for c1, c2 in itertools.combinations(putative_partition, 2))
        # ----------------------------------------------------------------------------------------
        def resolve_discordant_clusters(single_cluster, single_annotation, cluster_list, annotation_list, tdbg=False):  # reapportions all uids from single_cluster and cluster_list into return_clusts (splitting first by cdr3 and then by naive hamming distance)
            # NOTE single_cluster and cluster_list in general have quite different sets of uids, and that's fine. All that matters here is we're trying to find all the clusters that should be split from one another (without doing some all against all horror)
            if len(cluster_list) == 1:  # nothing to do
                return [single_cluster]
            adict = utils.get_annotation_dict(annotation_list)
            cdr3_groups = utils.group_seqs_by_value(cluster_list, lambda c: adict[akey(c)]['cdr3_length'])  # there's already utils.split_clusters_by_cdr3(), but it uses different inputs (e.g. sw_info) so i think it makes sense to not use it here
            if tdbg:
                print '   %s one cluster vs %d clusters' % (utils.color('blue', 'syncing'), len(cluster_list))
                print '     split into %d cdr3 groups' % len(cdr3_groups)
            lo_hbound, hi_hbound = utils.get_naive_hamming_bounds('likelihood', overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))  # these are the wider bounds, so < lo is almost certainly clonal, > hi is almost certainly not
            return_clusts = []
            for icdr, cdrgroup in enumerate(cdr3_groups):  # within each cdr3 group, split (i.e. use the cluster boundaries from cluster_list rather than single_cluster) if naive hfrac is > hi_hbound (but then there's shenanigans to adjudicate between different possibilities)
                if tdbg: print utils.color('purple', '      icdr %d' % icdr)

                # first figure out who needs to be split from whom
                clusters_to_split = {akey(c) : [] for c in cdrgroup}  # map from each cluster ('s key) to a list of clusters from which it should be split
                for c1, c2 in itertools.combinations(cdrgroup, 2):  # we could take account of the hfrac of both chains at this point, but looking at only the "split" one rather than the "merged" one, as we do here, is i think equivalent to assuming the merged one has zero hfrac, which is probably fine, since we only split if the split chain is very strongly suggesting we split
                    hfrac = utils.hamming_fraction(adict[akey(c1)]['naive_seq'], adict[akey(c2)]['naive_seq'])  # all clusters with the same cdr3 len have been padded in waterer so their naive seqs are the same length
                    if hfrac > hi_hbound:
                        clusters_to_split[akey(c1)].append(c2)
                        clusters_to_split[akey(c2)].append(c1)

                # then do the splitting
                if tdbg:
                    print '                  N to     new'
                    print '          size    split   cluster?'
                tmpclusts_for_return = []  # final (return) clusters for this cdr3 class
                for cclust in cdrgroup:
                    split_clusts = clusters_to_split[akey(cclust)]
                    if tdbg: print '         %4d    %3d' % (len(cclust), len(split_clusts)),
                    found_one = False
                    for rclust in tmpclusts_for_return:  # look for an existing return cluster to which we can merge cclust, i.e. that doesn't have any uids from which we want to split
                        if any_in_common([rclust], split_clusts):  # if any uid in rclust is in a cluster from which we want to be split, skip it, i.e. don't merge with that cluster (note that we have to do it by uid because the rclusts are already merged so don't necessarily correspond to any existing cluster)
                            continue
                        if found_one: print 'it happened!'  # TODO remove this. i'm just kind of curious if it ever happens in practice
                        if tdbg: print '     merging with size %d' % len(rclust)
                        rclust += cclust
                        found_one = True
                        break  # i.e. we just merge with the first one we find and stop looking; if there's more than one, it means we could merge all three together if we wanted (see diagram on some paper somewhere [triangle inequality]), but i doubt it'll matter either way, and this is easier
                    if not found_one:
                        if tdbg: print '      y'
                        tmpclusts_for_return.append(cclust)  # if we didn't find an existing cluster that we can add it to, add it as a new cluster

                return_clusts += tmpclusts_for_return

            if debug:
                print '      returning: %s' % ' '.join([str(len(c)) for c in return_clusts])
                # ClusterPath(partition=return_clusts).print_partitions(abbreviate=True)
            return return_clusts

        # ----------------------------------------------------------------------------------------
        if h_ipart is None:
            h_ipart = self.i_best
        if l_ipart is None:
            l_ipart = lcpath.i_best
        if debug:
            self.print_partitions(extrastr=utils.color('blue', 'heavy  '), print_partition_indices=True, ipart_center=h_ipart, n_to_print=1, sort_by_size=False)
            lcpath.print_partitions(extrastr=utils.color('blue', 'light  '), print_partition_indices=True, ipart_center=l_ipart, n_to_print=1, sort_by_size=False, print_header=False)
        if h_ipart != self.i_best or l_ipart != lcpath.i_best:
            print '  %s using non-best partition index for%s%s' % (utils.color('red', 'note'), (' heavy: %d'%h_ipart) if h_ipart != self.i_best else '', (' light: %d'%l_ipart) if l_ipart != lcpath.i_best else '')

        init_partitions = {'h' : self.partitions[h_ipart], 'l' : lcpath.partitions[l_ipart]}
        common_uids, _, _ = utils.check_intersection_and_complement(init_partitions['h'], init_partitions['l'], only_warn=True, a_label='heavy', b_label='light')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
        if len(common_uids) == 0:
            raise Exception('no uids in common between heavy and light')
        annotation_dict = {'h' : utils.get_annotation_dict(heavy_annotations), 'l' : utils.get_annotation_dict(light_annotations)}

        final_partition = []
        if debug:
            print '    N        N       hclusts     lclusts       h/l'
            print '  hclusts  lclusts    sizes       sizes      overlaps'
        for h_initclust, l_initclust in [(c, None) for c in init_partitions['h']] + [(None, c) for c in init_partitions['l']]:  # just loops over each single cluster in h and l partitions, but in a way that we know whether the single cluster is from h or l
            single_chain, list_chain = 'h' if l_initclust is None else 'l', 'l' if l_initclust is None else 'h'
            single_cluster = h_initclust if single_chain == 'h' else l_initclust
            cluster_list = common_clusters(single_cluster, init_partitions[list_chain])
            single_annotation = annotation_dict[single_chain][akey(single_cluster)]
            annotation_list = [annotation_dict[list_chain][akey(c)] for c in cluster_list]

            if debug:
                hclusts, lclusts = ([single_cluster], cluster_list) if single_chain == 'h' else (cluster_list, [single_cluster])
                overlaps = [[len(set(hc) & set(lc)) for lc in lclusts] for hc in hclusts]
                overlapstr = '   '.join([' '.join(str(ov) for ov in ovlist) for ovlist in overlaps])
                def getcstr(clist): return ' '.join(str(len(c)) for c in clist)
                hcstr, lcstr = getcstr(hclusts), getcstr(lclusts)
                cw = 10
                if len(hcstr) < cw and len(lcstr) < cw:  # fits on a single line
                    print ('    %2d      %2d         %-'+str(cw)+'s  %-'+str(cw)+'s  %s') % (len(hclusts), len(lclusts), hcstr, lcstr, overlapstr)
                else:  # split the last few columns over multiple lines
                    print ('    %2d      %2d         %-s') % (len(hclusts), len(lclusts), hcstr)
                    print ('    %2s      %2s         %-'+str(cw)+'s%-s') % ('', '', '', lcstr)
                    print ('    %2s      %2s         %-'+str(cw)+'s%-'+str(cw)+'s   %s') % ('', '', '', '', overlapstr)

            resolved_clusters = resolve_discordant_clusters(copy.deepcopy(single_cluster), single_annotation, copy.deepcopy(cluster_list), annotation_list)
            if check_partitions:
                assert is_clean_partition(resolved_clusters)
            if debug:
                print '    adding %d resolved cluster%s to %d clusters in final partition' % (len(resolved_clusters), utils.plural(len(resolved_clusters)), len(final_partition))
                print '      ifclust N rclusts'
            n_clean = 0
            for ifclust in range(len(final_partition)):  # iteration won't get as far as any clusters that we're just adding, which is what we want
                fclust = final_partition[ifclust]
                if not any_in_common([fclust], resolved_clusters):  # this is probably faster than combining it with getting the common cluster indices below, but maybe not
                    n_clean += 1
                    continue
                irclusts = common_clusters(fclust, resolved_clusters, return_indices=True)  # indices of any resolved_clusters that overlap with this fclust
                if debug: dbgstr = []
                new_fset = set(fclust)  # we'll remove uids from this, and then replace fclust with its remains
                for irclust in irclusts:  # resolve any discrepancies between these newly-resolved clusters and fclust
                    rset = set(resolved_clusters[irclust])
                    common_uids = new_fset & rset
                    if len(new_fset) > len(rset):  # remove the common ids from the larger one (effectively splitting according to the splittier one)
                        new_fset -= common_uids
                        if debug: dbgstr.append('  fclust %d --> %d' % (len(new_fset) + len(common_uids), len(new_fset)))
                    else:
                        rset -= common_uids
                        if debug: dbgstr.append('  rclust %d --> %d' % (len(rset) + len(common_uids), len(rset)))
                    resolved_clusters[irclust] = list(rset)
                if debug:
                    print '       %4d  %4d  %s' % (ifclust, len(irclusts), ''.join(dbgstr))
                final_partition[ifclust] = list(new_fset)
            if debug:
                print '       %d fclusts clean' % n_clean
            assert is_clean_partition(resolved_clusters)
            final_partition += resolved_clusters

        if debug:
            print '    removing %d/%d empty clusters' % (final_partition.count([]), len(final_partition))
        final_partition = [c for c in final_partition if len(c) > 0]
        # if debug:
        #     print '    final: %s' % ' '.join([str(len(c)) for c in final_partition])
        def chstr(n_before, n_after):
            if n_before == n_after: return ''
            else: return ' ' + utils.color('red', '%+d' % (n_after - n_before))
        print '   N clusters:\n        h %4d --> %-4d%s\n        l %4d --> %-4d%s'  % (len(init_partitions['h']), len(final_partition), chstr(len(init_partitions['h']), len(final_partition)),
                                                                                       len(init_partitions['l']), len(final_partition), chstr(len(init_partitions['l']), len(final_partition)))

        if check_partitions:
            assert is_clean_partition(final_partition)
            for initpart in init_partitions.values():
                assert len(set([u for c in initpart for u in c]) - set([u for c in final_partition for u in c])) == 0  # everybody from both initial partitions is in final_partition
            assert len(set([u for c in final_partition for u in c]) - set([u for c in init_partitions['h'] for u in c]) - set([u for c in init_partitions['l'] for u in c])) == 0  # nobody extra got added (i don't see how this could happen, but maybe it's just checking that I didnt' modify the initial partitions)

        return final_partition
