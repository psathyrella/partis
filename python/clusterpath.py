import os
import sys
import math
import csv

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
class ClusterPath(object):
    def __init__(self, initial_path_index):
        self.initial_path_index = initial_path_index

        # NOTE make *damn* sure if you add another list here that you also take care of it in remove_first_partition()
        self.n_lists = 5  # just to make sure you don't forget
        self.partitions = []
        self.logprobs = []
        self.logweights = []
        self.adj_mis = []
        self.n_procs = []

        self.best_minus = 30.  # rewind by this many units of log likelihood when merging separate processes
        self.i_best, self.i_best_minus_x = None, None
        # self.conservative_best_minus_ten_partition, self.conservative_max_minus_ten_logprob = None, None  # have to keep track of these *separately* since they don't necessarily occur in <self.partitions>, if this path is the result of merging a number of others

    def update_best_minus_x_partition(self):
        if math.isinf(self.logprobs[self.i_best]):  # if logprob is infinite, set best and best minus x to the latest one
            self.i_best_minus_x = self.i_best
            return
        for ip in range(len(self.partitions)):  # they should be in order of increasing logprob
            if self.logprobs[ip] > self.logprobs[self.i_best] - self.best_minus and self.n_procs[ip] == self.n_procs[self.i_best]:  # TODO is this exactly what I want to do with n_procs?
                self.i_best_minus_x = ip
                # self.conservative_best_minus_ten_partition = self.partitions[self.i_best_minus_ten]  # override later if necessary
                # self.conservative_max_minus_ten_logprob = self.logprobs[self.i_best_minus_ten]
                break

    def add_partition(self, partition, logprob, n_procs, logweight=None, adj_mi=None):
        # # don't add it if it's the same as the last partition
        # UPDATE we can get in trouble if we don't add duplicates, because they can have different numbers of procs... the duplicates don't really matter, anyway, as for most purposes I ignore partitions with greater than 1 proc
        # if len(self.partitions) > 0 and len(partition) == len(self.partitions[-1]) and logprob == self.logprobs[-1]:
        #     print 'NOT ADDING with n_procs %d --> %d' % (self.n_procs[-1], n_procs)
        #     return
        self.partitions.append(partition)
        self.logprobs.append(logprob)
        self.logweights.append(logweight)
        self.adj_mis.append(adj_mi)
        self.n_procs.append(n_procs)
        if math.isinf(logprob) or self.i_best is None or logprob > self.logprobs[self.i_best] or n_procs < self.n_procs[self.i_best]:  # if we haven't set i_best yet, or if this partition is more likely than i_best, or if i_best is set for a larger number of procs
            self.i_best = len(self.partitions) - 1
        self.update_best_minus_x_partition()

    def remove_first_partition(self):
        # NOTE after you do this, none of the 'best' shit is any good any more
        assert self.n_lists == 5  # make sure we didn't add another list and forget to put it in here
        self.partitions.pop(0)
        self.logprobs.pop(0)
        self.logweights.pop(0)
        self.adj_mis.pop(0)
        self.n_procs.pop(0)

    # ----------------------------------------------------------------------------------------
    def readfile(self, fname):
        with opener('r')(fname) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                partition = [cl.split(':') for cl in line['clusters'].split(';')]
                logweight = float(line['logweight']) if 'logweight' in line else None
                adj_mi = float(line['adj_mi']) if 'adj_mi' in line else None
                self.add_partition(partition, float(line['logprob']), int(line['n_procs']), logweight=logweight, adj_mi=adj_mi)

    # ----------------------------------------------------------------------------------------
    def print_partition(self, ip, reco_info=None, extrastr='', one_line=True, abbreviate=True):
        if one_line:
            if ip > 0:  # delta between this logprob and the previous one
                delta_str = '%.1f' % (self.logprobs[ip] - self.logprobs[ip-1])
            else:
                delta_str = ''
            print '      %5s  %-12.2f%-7s   %-5d  %5d' % (extrastr, self.logprobs[ip], delta_str, len(self.partitions[ip]), self.n_procs[ip]),

            # logweight (and inverse of number of potential parents)
            if self.logweights[ip] is not None:
                way_str, logweight_str = '', ''
                expon = math.exp(self.logweights[ip])
                n_ways = 0 if expon == 0. else 1. / expon
                way_str = ('%.1f' % n_ways) if n_ways < 1e7 else ('%8.1e' % n_ways)
                logweight_str = '%8.3f' % self.logweights[ip]

            # adj mi
            if reco_info is not None:
                adj_mi_str = ''
                if self.adj_mis[ip] is None:
                    adj_mi_str = 'na'
                else:
                    if self.adj_mis[ip] > 1e-3:
                        adj_mi_str = '%-8.3f' % self.adj_mis[ip]
                    else:
                        adj_mi_str = '%-8.0e' % self.adj_mis[ip]
                print '      %8s   ' % (adj_mi_str),
            if self.logweights[ip] is not None:
                print '   %10s    %8s   ' % (way_str, logweight_str),
        else:
            print '  %5s partition   %-15.2f' % (extrastr, self.logprobs[ip]),
            if reco_info is not None:
                print '    %-8.2f' % (self.adj_mis[ip]),
            print ''
            print '   clonal?   ids'

        # clusters
        for cluster in self.partitions[ip]:
            same_event = utils.from_same_event(reco_info is None, reco_info, cluster)
            if same_event is None:
                same_event = -1

            if abbreviate:
                cluster_str = ':'.join(['o' for uid in cluster])
            else:
                cluster_str = ':'.join([str(uid) for uid in cluster])
            if not same_event:
                cluster_str = utils.color('red', cluster_str)
            
            if one_line:
                if abbreviate:
                    print ' %s' % cluster_str,
                else:
                    print '   %s' % cluster_str,
            else:
                print '     %d    %s' % (int(same_event), cluster_str)
        if one_line:
            print ''

    # ----------------------------------------------------------------------------------------
    def get_partition_subset(self, n_partitions):
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
    def print_partitions(self, reco_info=None, extrastr='', one_line=True, abbreviate=True, print_header=True, n_to_print=None):
        if print_header:
            print '    %7s %10s   %-7s %5s  %5s' % ('', 'logprob', 'delta', 'clusters', 'n_procs'),
            if reco_info is not None:
                print ' %8s' % ('adj mi'),
            if self.logweights[0] is not None:
                print '  %10s  %7s' % ('pot.parents', 'logweight'),
            print ''

        for ip in self.get_partition_subset(n_partitions=n_to_print):
            mark = ''
            if ip == self.i_best:
                mark += '*'
            if ip == self.i_best_minus_x:
                mark += '*'
            self.print_partition(ip, reco_info, extrastr=mark+extrastr, one_line=one_line, abbreviate=abbreviate)

    # ----------------------------------------------------------------------------------------
    def set_synthetic_logweight_history(self, reco_info):
        # TODO switch clusterpath.cc back to using these
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
    def write_partitions(self, writer, is_data, reco_info, true_partition, smc_particles, path_index, n_to_write=None, calc_adj_mi=None):
        for ipart in self.get_partition_subset(n_partitions=n_to_write):
            part = self.partitions[ipart]
            cluster_str = ''
            bad_clusters = []  # inferred clusters that aren't really all from the same event
            for ic in range(len(part)):
                if ic > 0:
                    cluster_str += ';'
                cluster_str += ':'.join(part[ic])
                if not is_data:
                    same_event = utils.from_same_event(is_data, reco_info, part[ic])  # are all the sequences from the same event?
                    entire_cluster = True  # ... and if so, are they the entire true cluster?
                    if same_event:
                        reco_id = reco_info[part[ic][0]]['reco_id']  # they've all got the same reco_id then, so pick an aribtrary one
                        true_cluster = true_partition[reco_id]
                        for uid in true_cluster:
                            if uid not in part[ic]:
                                entire_cluster = False
                                break
                    else:
                        entire_cluster = False
                    if not same_event or not entire_cluster:
                        bad_clusters.append(':'.join(part[ic]))

            if len(bad_clusters) > 25:
                bad_clusters = ['too', 'long']
            row = {'logprob' : self.logprobs[ipart],
                   'n_clusters' : len(part),
                   'n_procs' : self.n_procs[ipart],
                   'clusters' : cluster_str}
            if smc_particles > 1:
                row['path_index'] = path_index
                row['logweight'] = self.logweights[ipart]
            if not is_data:
                if calc_adj_mi is None or self.adj_mis[ipart] != -1:  # if we don't want to write any adj mis, or if we already calculated it
                    row['adj_mi'] = self.adj_mis[ipart]
                else:
                    if calc_adj_mi == 'best' and ipart == self.i_best:  # only calculate adj_mi for the best partition
                        row['adj_mi'] = utils.mutual_information_to_true(part, reco_info)
                    else:
                        row['adj_mi'] = self.adj_mis[ipart]
                row['n_true_clusters'] = len(true_partition)
                row['bad_clusters'] = ';'.join(bad_clusters)
            writer.writerow(row)
