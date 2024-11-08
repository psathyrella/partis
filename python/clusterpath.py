from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import dendropy
import os
import sys
import math
import csv
import copy
import collections
import itertools
import numpy
from io import open
csv.field_size_limit(sys.maxsize)

from . import utils
from . import treeutils

# ----------------------------------------------------------------------------------------
# print a single partition without having to make a cluster path
def ptnprint(partition, **kwargs):
    ClusterPath(partition=partition).print_partitions(**kwargs)

# ----------------------------------------------------------------------------------------
class ClusterPath(object):
    def __init__(self, initial_path_index=0, seed_unique_id=None, partition=None, fname=None, partition_lines=None, add_pairwise_metrics=False):  # <partition> is a fully-formed partition, while <partition_lines> is straight from reading a file (perhaps could combine them, but I don't want to think through it now)
        # could probably remove path index since there's very little chance of doing smc in the future, but the path-merging code in glomerator was _very_ difficult to write, so I'm reluctant to nuke it
        self.initial_path_index = initial_path_index  # NOTE this is set to None if it's nonsensical, e.g. if we're merging several paths with different indices
        self.add_pairwise_metrics = add_pairwise_metrics

        # NOTE see remove_partition()
        self.list_names = ['partitions', 'logprobs', 'n_procs', 'ccfs', 'perf_metrics', 'logweights']
        # ccfs: pair of floats (not just a float) for each partition
        # perf_metrics: new, more general way to store performance metrics (could eventually remove old ccfs, but they're used in a ton of places so not doing it now)
        for lname in self.list_names:
            setattr(self, lname, [])

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
    def best(self):  # return best partition NOTE adding this very late, so there's a ton of places where i could go back and use it
        return self.partitions[self.i_best]

    # ----------------------------------------------------------------------------------------
    def last(self):
        return self.partitions[len(self.partitions) - 1]

    # ----------------------------------------------------------------------------------------
    def bmx(self):
        return self.partitions[self.i_best_minus_x]

    # ----------------------------------------------------------------------------------------
    def n_seqs(self, ip=0):  # number of sequences in each partition (shouldn't depend on which partition, but I'm not sure that I absolutely forbid the number to change)
        if len(self.partitions) == 0:
            return 0
        return len([u for c in self.partitions[ip] for u in c])

    # ----------------------------------------------------------------------------------------
    def find_iparts_for_cluster(self, cluster):  # get index of partitions in which a list of uids (i.e. a cluster) appears
        return [ip for ip in range(len(self.partitions)) if cluster in self.partitions[ip]]  # NOTE just returns zero-length list if it isn't there

    # ----------------------------------------------------------------------------------------
    def get_headers(self, is_simu):
        headers = ['logprob', 'n_clusters', 'n_procs', 'partition']
        if is_simu:
            headers += ['n_true_clusters', 'ccf_under', 'ccf_over', 'perf_metrics']
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
    # NOTE this will presumably screw up self.logprobs, self.n_procs, self.ccfs, self.perf_metrics, and self.logweights
    def add_cluster_to_all_partitions(self, cluster, allow_duplicates=False, skip_duplicates=False, debug=True):
        # ----------------------------------------------------------------------------------------
        def check(iptn, cluster, duplicate_ids):
            new_duplicates = set(cluster) & set(u for c in self.partitions[iptn] for u in c)
            if not allow_duplicates:
                raise Exception('adding cluster with uids %s already in partition' % ' '.join(new_duplicates))
            if debug:
                duplicate_ids |= set(new_duplicates)
        # ----------------------------------------------------------------------------------------
        assert len(cluster) > 0
        duplicate_ids = set()
        for iptn in range(len(self.partitions)):
            if any(len(set(cluster) & set(c)) > 0 for c in self.partitions[iptn]):
                if skip_duplicates:
                    continue
                check(iptn, cluster, duplicate_ids)
            self.partitions[iptn].append(cluster)
        if debug and len(duplicate_ids) > 0:
            print('  %s %d uid%s appeared multiple times when adding cluster to all partitions%s' % (utils.color('yellow', 'warning'), len(duplicate_ids), utils.plural(len(duplicate_ids)), ': '+' '.join(duplicate_ids) if len(duplicate_ids) < 10 else ''))

    # ----------------------------------------------------------------------------------------
    def add_partition(self, partition, logprob, n_procs, logweight=None, ccfs=None, perf_metrics=None):
        if ccfs is None:
            ccfs = [None, None]
        if perf_metrics is None:
            perf_metrics = {}
        # NOTE you typically want to allow duplicate (in terms of log prob) partitions, since they can have different n procs
        self.partitions.append(partition)  # NOTE not deep copied
        self.logprobs.append(logprob)
        self.n_procs.append(n_procs)
        self.logweights.append(logweight)
        if len(ccfs) != 2:
            raise Exception('tried to add partition with ccfs of length %d (%s)' % (len(ccfs), ccfs))
        self.ccfs.append(ccfs)
        self.perf_metrics.append(perf_metrics)
        if ccfs.count(None) != len(ccfs):
            self.we_have_a_ccf = True  # also indicates that we potentially have perf_metrics
        # set this as the best partition if 1) we haven't set i_best yet, 2) this partition is more likely than i_best, or 3) i_best is set for a partition with a larger number of procs or 4) logprob is infinite (i.e. we didn't calculate the full partitions logprob))
        if self.i_best is None or logprob > self.logprobs[self.i_best] or n_procs < self.n_procs[self.i_best] or math.isinf(logprob):
            self.i_best = len(self.partitions) - 1
        self.update_best_minus_x_partition()
        self.trees = None  # they'll be out of date if the best partition changed (so I guess in principle we could only do this if self.i_best changes, but this seems tidier)

    # ----------------------------------------------------------------------------------------
    def remove_partition(self, ip_to_remove):  # NOTE doesn't update self.we_have_a_ccf, but it probably won't change, right?
        for lname in self.list_names:
            getattr(self, lname).pop(ip_to_remove)

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
                    print('%s seed uids for each line not all the same %s %s' % (utils.color('yellow', 'warning'), line['seed_unique_id'], self.seed_unique_id))

            if process_csv:
                line['partition'] = [cluster_str.split(':') for cluster_str in line['partition'].split(';')]

            ccfs = [None, None]
            if 'ccf_under' in line and 'ccf_over' in line:  # I don't know what I want to do if there's one but not the other, but it shouldn't be possible
                if all(v not in ['', -1] for v in [line['ccf_under'], line['ccf_over']]):  # we shouldn't be writing any files any more with -1 set as a default value, but paired clustering used to do it in bin/partis.py, so we gotta handle it for old files
                    ccfs = [float(line['ccf_under']), float(line['ccf_over'])]
                self.we_have_a_ccf = True
            perf_metrics = {}
            if 'perf_metrics' in line:
                perf_metrics = line['perf_metrics']

            self.add_partition(line['partition'], float(line['logprob']), int(line.get('n_procs', 1)), logweight=float(line.get('logweight', 0)), ccfs=ccfs, perf_metrics=perf_metrics)

    # ----------------------------------------------------------------------------------------
    def calculate_missing_values(self, reco_info=None, true_partition=None, only_ip=None, fail_frac=None):  # NOTE adding true_partition argument very late, so there's probably lots of places that call the fcns that call this (printer + line getter) that would rather pass in true_partition and reco_info
        for ip in range(len(self.partitions)):
            if only_ip is not None and ip != only_ip:
                continue

            if self.ccfs[ip][0] is not None and self.ccfs[ip][1] is not None:  # already have them
                continue

            if true_partition is None and reco_info is not None:
                true_partition = utils.get_partition_from_reco_info(reco_info)  # full true partition, which we may have to modify below
            cfpart, cftrue = self.partitions[ip], true_partition
            if self.seed_unique_id is not None:  # for seed clustering, we need to add to the inferred partition any ids from the true seed cluster that we missed (as well as make sure to keep [in true partition] any non-seed-family ids that're in the inferred partition)
                true_seed_cluster = utils.get_single_entry([c for c in true_partition if self.seed_unique_id in c])
                missing_ids = list(set(true_seed_cluster) - utils.ptn_ids(cfpart))
                cfpart = copy.deepcopy(cfpart) + [missing_ids]  # add any missing true seed cluster ids to the inferred partition (as a single cluster)
                cftrue = [list(set(c) & utils.ptn_ids(cfpart)) for c in cftrue]  # then remove from true partition ids not in the inferred partition (i.e. non-seed ids), except for any non-seed ones that are in the inferred partition
                cftrue = [c for c in cftrue if len(c) > 0]  # remove empties
            else:  # whereas for non-seed, add any missing (failed) ids to the inferred partition as singletons (fail if there's too many)
                cfpart = utils.add_missing_uids_to_partition(cfpart, true_partition, warn=True, fail_frac=fail_frac, ref_label='true', miss_label='inferred')
            new_vals = utils.per_seq_correct_cluster_fractions(cfpart, cftrue, reco_info=reco_info, seed_unique_id=self.seed_unique_id)  # it doesn't need reco_info, but having it does save a step
            if None not in new_vals:  # if the function finds messed up partitions, it returns None, None (at this point, this seems to just happens when a uid was found in multiple clusters, which happens for earlier partitions (n_procs > 1) when seed_unique_id is set, since we pass seed_unique_id to all the subprocs)
                self.ccfs[ip] = new_vals
                self.we_have_a_ccf = True
            self.perf_metrics[ip] = {}
            if self.add_pairwise_metrics:
                self.perf_metrics[ip]['pairwise'] = utils.pairwise_cluster_metrics('pairwise', cfpart, cftrue, debug=True)  # eventually might be nice to also include ccf numbers in here (which shouldn't really even be called ccfs)

    # ----------------------------------------------------------------------------------------
    def get_ccf_str(self, ip):
        ccf_str_list = [('%5s' % '-') if ccf is None else ('%5.2f' % ccf) for ccf in self.ccfs[ip]]
        return ' %s ' % ' '.join(ccf_str_list)

    # ----------------------------------------------------------------------------------------
    def print_partition(self, ip, reco_info=None, extrastr='', abbreviate=True, dont_print_clusters=False, highlight_cluster_indices=None, print_partition_indices=False, right_extrastr='', sort_by_size=True, max_sizestr_len=0):  # NOTE <highlight_cluster_indices> and <print_partition_indices> are quite different despite sounding similar, but I can't think of something else to call the latter that makes more sense
        # ----------------------------------------------------------------------------------------
        def ccol(tclust, cstr, iclust):
            if reco_info is not None and not utils.from_same_event(reco_info, tclust):
                cstr = utils.color('red', cstr)
            if self.seed_unique_id is not None and self.seed_unique_id in tclust:
                cstr = utils.color('reverse_video', cstr)
            if highlight_cluster_indices is not None and iclust in highlight_cluster_indices:
                cstr = utils.color('red', cstr)
            return cstr
        # ----------------------------------------------------------------------------------------
        if ip > 0:  # delta between this logprob and the previous one
            delta_str = '%.1f' % (self.logprobs[ip] - self.logprobs[ip-1])
        else:
            delta_str = ''
        print('      %s  %-12.2f%-7s   %s%-5d  %4d' % (extrastr, self.logprobs[ip], delta_str, ('%-5d  ' % ip) if print_partition_indices else '', len(self.partitions[ip]), self.n_procs[ip]), end=' ')

        if self.we_have_a_ccf:
            print('    ' + self.get_ccf_str(ip), end=' ')

        if not dont_print_clusters:
            sorted_clusters = self.partitions[ip]
            if sort_by_size:  # it's often nicer to *not* sort by cluster size here, since preserving the order frequently makes it obvious which clusters are merging as your eye scans downward through the output
                sorted_clusters = sorted(sorted_clusters, key=lambda c: len(c), reverse=True)
            n_total_seqs = sum(len(c) for c in self.partitions[ip])
            csstr = ' '.join(ccol(c, str(len(c)), i) for i, c in enumerate(sorted_clusters) if len(c) > 1)  # would be nice to use utils.cluster_size_str(), but i need more color configurability atm
            csstr += ' (+%d)' % len([c for c in sorted_clusters if len(c)==1])
            print('   %s%s' % (csstr, ' '*(max_sizestr_len - utils.len_excluding_colors(csstr))), end=' ')
            for iclust in range(len(sorted_clusters)):
                tclust = sorted_clusters[iclust]
                cstr = ''
                if n_total_seqs < 100:
                    def ustr(u): return 'o' if len(u) > 3 and abbreviate else str(u)
                    cstr = ':'.join(ustr(u) for u in tclust)
                cstr = ccol(tclust, cstr, iclust)
                print(' %s%s' % ('' if abbreviate else '  ', cstr), end=' ')
        print('%s' % right_extrastr, end=' ')
        print('')

    # ----------------------------------------------------------------------------------------
    def print_partitions(self, reco_info=None, true_partition=None, extrastr='', abbreviate=True, dont_print_clusters=False, print_header=True, n_to_print=None, calc_missing_values='none',
                         highlight_cluster_indices=None, print_partition_indices=False, ipart_center=None, sort_by_size=True):
        assert calc_missing_values in ['none', 'all', 'best']
        if (reco_info is not None or true_partition is not None) and calc_missing_values == 'all':
            self.calculate_missing_values(reco_info=reco_info, true_partition=true_partition)

        if print_header:
            print('    %s%7s %10s   %-7s %s%5s  %4s' % (' '*utils.len_excluding_colors(extrastr), '', 'logprob', 'delta', 'index  ' if print_partition_indices else '', 'clusters', 'n_procs'), end=' ')
            if reco_info is not None or self.we_have_a_ccf:
                print(' %5s %5s' % ('purity', 'completeness'), end=' ')
            print('    sizes')
            print('')

        ptns_to_print = self.get_surrounding_partitions(n_to_print, i_center=ipart_center)
        max_sizestr_len = max(len(' '.join(str(len(c)) for c in self.partitions[i])) for i in ptns_to_print)
        for ip in ptns_to_print:
            if (reco_info is not None or true_partition is not None) and calc_missing_values == 'best' and ip == self.i_best:
                self.calculate_missing_values(reco_info=reco_info, true_partition=true_partition, only_ip=ip)
            mark = '      '
            if ip == self.i_best:
                mark = 'best  '
            if ip == self.i_best_minus_x:
                mark = mark[:-2] + '* '
            if mark.count(' ') < len(mark):
                mark = utils.color('yellow', mark)
            right_extrastr = '' if self.n_seqs() < 200 else mark  # if line is going to be really long, put the yellow stuff also on the right side
            self.print_partition(ip, reco_info, extrastr=extrastr+mark, abbreviate=abbreviate, dont_print_clusters=dont_print_clusters, highlight_cluster_indices=highlight_cluster_indices,
                                 print_partition_indices=print_partition_indices, right_extrastr=right_extrastr, sort_by_size=sort_by_size, max_sizestr_len=max_sizestr_len)

    # ----------------------------------------------------------------------------------------
    def get_surrounding_partitions(self, n_partitions, i_center=None):
        """ return a list of partition indices centered on <self.i_best> of length <n_partitions> """
        if i_center is None:
            i_center = self.i_best
        if n_partitions is None:  # print all partitions
            ilist = list(range(len(self.partitions)))
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
    def write(self, outfname, is_data, reco_info=None, true_partition=None, n_to_write=None, calc_missing_values='none', partition_lines=None):  # NOTE could now i think also remove is_data from here (like i just did for get_partition_lines())
        if utils.getsuffix(outfname) != '.csv':
            raise Exception('unhandled file extension %s' % outfname)
        if partition_lines is None:
            partition_lines = self.get_partition_lines(reco_info=reco_info, true_partition=true_partition, n_to_write=n_to_write, calc_missing_values=calc_missing_values)
        with open(outfname, utils.csv_wmode()) as outfile:
            writer = csv.DictWriter(outfile, self.get_headers(not is_data))
            writer.writeheader()
            for row in partition_lines:
                row['partition'] = ';'.join([':'.join(cluster) for cluster in row['partition']])
                if 'bad_clusters' in row:
                    row['bad_clusters'] = ';'.join(row['bad_clusters'])
                writer.writerow(row)

    # ----------------------------------------------------------------------------------------
    def get_partition_lines(self, reco_info=None, true_partition=None, n_to_write=None, calc_missing_values='none', path_index=None, fail_frac=0.05, add_pairwise_metrics=False):  # we use this (instead of .write()) if we're writing a yaml file
        if add_pairwise_metrics:
            self.add_pairwise_metrics = True  # kind of weird to reset it here (maybe overriding value from constructor) but whatevs
        assert calc_missing_values in ['none', 'all', 'best']
        is_simu = (reco_info is not None or true_partition is not None) or any(cfs!=[None, None] for cfs in self.ccfs)  # second clause is for cases where we've read a cpath with simulation info (i.e. ccfs), but don't have reco_info or true_partition
        if is_simu and calc_missing_values == 'all':
            self.calculate_missing_values(reco_info=reco_info, true_partition=true_partition, fail_frac=fail_frac)

        headers = self.get_headers(is_simu)  # NOTE it would be better (at least now, i think) to look at what is set e.g. in ccfs to determine headers, rather than use is_simu
        lines = []
        for ipart in self.get_surrounding_partitions(n_to_write):
            part = self.partitions[ipart]

            row = {'logprob' : self.logprobs[ipart],
                   'n_clusters' : len(part),
                   'n_procs' : self.n_procs[ipart],
                   'partition' : part}
            if 'ccf_under' in headers:
                if is_simu and calc_missing_values == 'best' and ipart == self.i_best:
                    self.calculate_missing_values(reco_info=reco_info, true_partition=true_partition, only_ip=ipart, fail_frac=fail_frac)
                if self.ccfs[ipart][0] is not None and self.ccfs[ipart][1] is not None:
                    row['ccf_under'], row['ccf_over'] = self.ccfs[ipart]  # for now assume we calculated the ccfs if we did adj mi
                    row['perf_metrics'] = self.perf_metrics[ipart]
            if 'n_true_clusters' in headers:
                row['n_true_clusters'] = '' if true_partition is None else len(true_partition)
            if 'bad_clusters' in headers:
                row['bad_clusters'] = '' if (reco_info is None or true_partition is None) else self.get_bad_clusters(part, reco_info, true_partition)
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
            print('  %s self.seed_unique_id isn\'t set' % utils.color('red', 'error'))
            return  # I don't like returning different types of things, but this lets the calling code handle things
        seed_clusters = [c for c in self.partitions[ipart if ipart is not None else self.i_best] if self.seed_unique_id in c]
        if len(seed_clusters) == 0:
            print('  %s couldn\'t find cluster containing seed %s' % (utils.color('red', 'error'), self.seed_unique_id))
        elif len(seed_clusters) > 1:
            print('  %s more than one cluster containing seed %s (this can be expected if this is the non-best partition)' % (utils.color('red', 'error'), self.seed_unique_id))
        else:
            return seed_clusters[0]

    # ----------------------------------------------------------------------------------------
    def print_seed_cluster_size(self, ipart=None, queries_to_include=None, lstr=''):
        seed_cluster = self.seed_cluster(ipart=ipart)
        qtistr = 'excluding seed seq%s:' % ('' if (queries_to_include is None or queries_to_include==[self.seed_unique_id]) else ' and --queries-to-include')
        if queries_to_include is None:
            queries_to_include = []
        non_qti_size = len(set(seed_cluster) - set(queries_to_include + [self.seed_unique_id]))  # NOTE don't modify queries_to_include
        qtistr = '%s %d' % (qtistr, non_qti_size)
        pstr = '%spartition' % lstr
        print('  seed cluster size in %s: %d, %s' % (('best %s'%pstr) if ipart is None else '%s with index %d (best at index %d)'%(pstr, ipart, self.i_best), len(seed_cluster), qtistr))

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
        print('   writing presto partition %s' % outfname)
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
        # ----------------------------------------------------------------------------------------
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
            raise Exception('couldn\'t find uid %s in annotations (this is likely because your partition history/cluster path doesn\'t maintain the same uids at each step, which for instance happens if you ran seed partitioning, in which case don\'t try to infer clusterpath trees)' % uidstr)
        # ----------------------------------------------------------------------------------------
        def getseq(uid):
            if uid == naive_seq_name:
                sorted_lines = sorted([l for l in annotations.values()], key=lambda l: len(l['unique_ids']), reverse=True)  # since we're making a tree, all the annotations are by definition clonal, so it doesn't really matter, but may as well get the naive sequence from the largest one
                return sorted_lines[0]['naive_seq']
            else:
                return utils.per_seq_val(getline(uid), 'seqs', uid)
        # ----------------------------------------------------------------------------------------
        def get_naive_seq(uidstr):
            return getline(uidstr)['naive_seq']
        # ----------------------------------------------------------------------------------------
        def lget(uid_list):
            return ':'.join(uid_list)
        # ----------------------------------------------------------------------------------------
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
            print('    starting tree with %d leaves and %d partition%s:' % (len(uid_set), len(partitions), utils.plural(len(partitions))))
            for ip, ptn in enumerate(partitions):
                ptnprint(ptn, extrastr='        ', print_header=ip==0)
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
                    print('      ipart %d' % ipart)
                    print('        split node: %d --> %s      %s --> %s' % (len(lnode.uids), ' '.join([str(len(tc)) for tc in tclusts]), lnode.taxon.label, ' '.join([c.taxon.label for c in lnode.child_node_iter()])))

        # split existing leaves, which are probably not singletons (they're probably from the initial naive sequence collapse step) into subtrees such that each leaf is a singleton
        for lnode in dtree.leaf_node_iter():
            if len(lnode.uids) == 1:
                continue
            if get_fasttrees and len(lnode.uids) > 2:
                seqfos = [{'name' : uid, 'seq' : getseq(uid)} for uid in lnode.taxon.label.split(':')]  # may as well add them in the right order, although I don't think it matters
                subtree, _, _ = treeutils.run_tree_inference('fasttree', input_seqfos=seqfos, suppress_internal_node_taxa=True, no_naive=True)  # note that the fasttree distances get ignored below (no idea if they'd be better than what we set down there, but they probably wouldn't be consistent, so I'd rather ignore them)
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
                    print('   adding subtree with %d leaves from fastree at leaf node %s' % (len(seqfos), lnode.taxon.label))
                    print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=subtree, width=350, label_fcn=lambda x: '<combo>' if ':' in x else x)))
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
                    print('      added %d singleton children for %s' % (len(lnode.uids), lnode.taxon.label))

        # in order to set edge lengths, we need node sequences, so first set leaf node seqs
        for lnode in dtree.leaf_node_iter():
            assert len(lnode.uids) == 1
            lnode.seq = getseq(lnode.taxon.label)
            lnode.n_descendent_leaves = 1  # keep track of how many leaf nodes contributed to each node's consensus sequence (these are leaves, so it's trivally 1). This is less accurate than keeping track of all the sequences, but also faster

        # then set internal node seqs as the consensus of their children, and set the distance as hamming distance to child seqs
        if debug:
            print('    adding edge lengths either from fasttree %s or cons seq %s' % (utils.color('blue', 'x'), utils.color('red', 'x')))
        min_edge_length = None  # setting this is nice for better debug viewing
        for node in dtree.postorder_internal_node_iter():  # includes root node
            child_cons_seq_counts = [c.n_descendent_leaves for c in node.child_node_iter()]
            total_descendent_leaves = sum(child_cons_seq_counts)
            if total_descendent_leaves > n_max_cons_seqs:  # if there's tons of descendent leaves, we don't want to pass them all to the consensus fcn since it's slow, so we choose them in proportion to their actual proportions, but scaled down to <n_max_cons_seqs>
                child_cons_seq_counts = [int(n_max_cons_seqs * csc / float(total_descendent_leaves)) for csc in child_cons_seq_counts]
                child_cons_seq_counts = [max(1, csc) for csc in child_cons_seq_counts]  # don't eliminate any sequences entirely (this makes the proportions less accurate (in some cases), but is the easy way to handle the case where there's a ton of singleton children
            if debug:
                print('  %s' % utils.color('green', node.taxon.label))
                csc_str = '  (reduced: %s)' % ' '.join([str(csc) for csc in child_cons_seq_counts]) if total_descendent_leaves > n_max_cons_seqs else ''
                print('      desc leaves per child: %s%s' % (' '.join(str(c.n_descendent_leaves) for c in node.child_node_iter()), csc_str))
            child_seqfos = [{'name' : cn.taxon.label + '-leaf-' + str(il), 'seq' : cn.seq} for cn, count in zip(node.child_node_iter(), child_cons_seq_counts) for il in range(count)]
            # node.seq = utils.old_bio_cons_seq(0.01, aligned_seqfos=child_seqfos, tie_resolver_seq=getline(root_label)['naive_seq'])  #, debug=debug)  # the consensus has an N at every position where the constituent sequences gave a tie. But Ns screw up the distances (especially because once we *get* an N, we can't get rid of it and it's propagated all the way up the tree), and in almost all cases the correct choice should be the naive base, so we use that
            node.seq = utils.cons_seq(aligned_seqfos=child_seqfos)  # UPDATE using new fcn that just picks one if there's a tie (and is ten times faster)
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
                    print('       %6.3f   %s  %s' % (edge.length, utils.color('blue' if from_fasttree else 'red', 'x'), edge.head_node.taxon.label))

        hfrac, hdist = utils.hamming_fraction(get_naive_seq(root_label), dtree.seed_node.seq, also_return_distance=True)
        if debug:
            if hdist > 0:
                print('  note: naive and root cons seq differ at %d positions (hfrac %.2f)' % (hdist, hfrac))
            print('        naive seq %s' % get_naive_seq(root_label)) # NOTE might be worthwhile to add an edge connecting seed node and the actual naive sequence (i.e. for cases where our approximate naive is off)
            print('    root cons seq %s' % utils.color_mutants(get_naive_seq(root_label), dtree.seed_node.seq))

        for node in dtree.preorder_node_iter():
            del node.uids
            del node.seq
            del node.n_descendent_leaves

        old_dtree = dtree  # this is copied from treeutils.get_tree_with_dummy_branches()
        naive_taxon = dendropy.Taxon(naive_seq_name)
        old_dtree.taxon_namespace.add_taxon(naive_taxon)
        naive_node = dendropy.Node(taxon=naive_taxon)
        dtree = dendropy.Tree(seed_node=naive_node, taxon_namespace=old_dtree.taxon_namespace, is_rooted=True)
        naive_node.add_child(old_dtree.seed_node)
        for edge in naive_node.child_edge_iter():
            edge.length = hfrac
        if debug:
            print('  added new naive root node at distance %.3f above previous root' % hfrac)

        treeutils.label_nodes(dtree, ignore_existing_internal_node_labels=True, ignore_existing_internal_taxon_labels=True, root_is_external=True, debug=debug)
        treeutils.collapse_zero_length_leaves(dtree, uid_set, debug=debug)
        # dtree.update_bipartitions(suppress_unifurcations=False)  # probably don't really need this UPDATE now i'm commenting it since it gets run by the zero length leaf collapse fcn
        if debug:
            print(treeutils.utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree, width=250)))

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
            print('  removed seed uid %s from %s clusters across %d partitions' % (self.seed_unique_id, ' '.join([str(n) for n in n_removed_list]), len(n_removed_list)))

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
    def make_trees(self, annotations, i_only_cluster=None, get_fasttrees=False, naive_seq_name='naive', debug=False):  # makes a tree for each cluster in the most likely (not final) partition
        # i_only_cluster: only make the tree corresponding to the <i_only_cluster>th cluster in the best partition
        if self.i_best is None:
            return

        partitions = self.partitions
        if self.seed_unique_id is not None:
            partitions = self.deduplicate_seed_uid(debug=debug)

        if debug:
            print('  making %d tree%s over %d partitions' % (len(partitions[self.i_best]) if i_only_cluster is None else 1, 's' if i_only_cluster is None else '', self.i_best + 1))
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
    def get_single_tree(self, line, get_fasttrees=False, debug=False):  # make tree for <line>
        i_only_cluster = self.best().index(line['unique_ids'])
        self.make_trees(annotations=utils.get_annotation_dict([line]), i_only_cluster=i_only_cluster, get_fasttrees=get_fasttrees, debug=debug)
        return self.trees[i_only_cluster]
