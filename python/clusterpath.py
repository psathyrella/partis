import os
import sys
import math

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
        for ip in range(len(self.partitions)):  # they should be in order of increasing logprob
            if self.logprobs[ip] > self.logprobs[self.i_best] - self.best_minus and self.n_procs[ip] == self.n_procs[self.i_best]:  # TODO is this exactly what I want to do with n_procs?
                self.i_best_minus_x = ip
                # self.conservative_best_minus_ten_partition = self.partitions[self.i_best_minus_ten]  # override later if necessary
                # self.conservative_max_minus_ten_logprob = self.logprobs[self.i_best_minus_ten]
                break

    def add_partition(self, partition, logprob, logweight, adj_mi, n_procs=-1):
        # don't add it if it's the same as the last partition
        if len(self.partitions) > 0 and len(partition) == len(self.partitions[-1]) and logprob == self.logprobs[-1]:
            return
        self.partitions.append(partition)
        self.logprobs.append(logprob)
        self.logweights.append(logweight)
        self.adj_mis.append(adj_mi)
        self.n_procs.append(n_procs)
        if self.i_best is None or n_procs < self.n_procs[self.i_best] or logprob > self.logprobs[self.i_best]:  # TODO is this exactly what I want to do with n_procs?
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
    def print_partition(self, ip, reco_info, extrastr='', one_line=False, abbreviate=True):
        if one_line:
            if ip > 0:
                delta_str = '%.1f' % (self.logprobs[ip] - self.logprobs[ip-1])
            else:
                delta_str = ''
            expon = math.exp(self.logweights[ip])
            n_ways = 0 if expon == 0. else 1. / expon
            way_str = ('%.1f' % n_ways) if n_ways < 1e7 else ('%8.1e' % n_ways)
            if self.adj_mis[ip] > 1e-3:
                adj_mi_str = '%-8.3f' % self.adj_mis[ip]
            else:
                adj_mi_str = '%-8.0e' % self.adj_mis[ip]
            print '      %5s  %-10.2f%-7s   %8s   %5d   %10s    %8.3f   ' % (extrastr, self.logprobs[ip], delta_str, adj_mi_str, len(self.partitions[ip]), way_str, self.logweights[ip]),
        else:
            print '  %5s partition   %-15.2f    %-8.2f' % (extrastr, self.logprobs[ip], self.adj_mis[ip])
            print '   clonal?   ids'
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
    def print_partitions(self, reco_info, extrastr='', one_line=False, abbreviate=True, header=True):
        if header:
            print '    %7s  %7s   %-7s %8s     %5s   %10s  %7s' % ('', 'logprob', 'delta', 'adj mi', 'clusters', 'pot.parents', 'logweight')
        for ip in range(len(self.partitions)):
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
