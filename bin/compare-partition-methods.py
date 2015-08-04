#!/usr/bin/env python
import os
import argparse
import glob
import sys
import numpy
import random
import re
from collections import OrderedDict
import time
import csv
from subprocess import check_call, Popen, check_output
import itertools
sys.path.insert(1, './python')

from humans import humans
from plotting import legends, colors, linewidths
import plotting
from hist import Hist
import seqfileopener
import utils
from clusterplot import ClusterPlot
from glomerator import Glomerator

fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output'
parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
parser.add_argument('--mutation-multipliers', default='1')
parser.add_argument('--data', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--indels', action='store_true')
parser.add_argument('--lonely-leaves', action='store_true')
parser.add_argument('--bak', action='store_true')
parser.add_argument('--count-distances', action='store_true')
parser.add_argument('--n-leaf-list', default='10')
parser.add_argument('--subset', type=int)
parser.add_argument('--n-to-partition', type=int, default=5000)
parser.add_argument('--n-data-to-cache', type=int, default=50000)
parser.add_argument('--n-sim-seqs', type=int, default=10000)
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
parser.add_argument('--startstoplist')  # list of istartstops for comparisons
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'write-plots', 'compare-sample-sizes', 'compare-subsets']
parser.add_argument('--actions', required=True)  #default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)
args.mutation_multipliers = utils.get_arg_list(args.mutation_multipliers, intify=True)
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, intify=True)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)
args.startstoplist = utils.get_arg_list(args.startstoplist)

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

if args.subset is not None:
    args.n_to_partition = 1300
if args.istartstop is not None:
    args.n_to_partition = args.istartstop[1] - args.istartstop[0]

# ----------------------------------------------------------------------------------------
procs = []
changeorandomcrapstr = '_db-pass_parse-select_clone-pass.tab'

# # ----------------------------------------------------------------------------------------
# legends['v-true'] = 'true (V indels)'
# legends['cdr3-true'] = 'true (CDR3 indels)'
# legends['v-indels'] = 'full partis (V indels)'
# legends['cdr3-indels'] = 'full partis (CDR3 indels)'
# args.is_data = False
# args.debug = False
# args.use_all_steps = False
# args.normalize_axes = []
# args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None

# input_info, reco_info = seqfileopener.get_seqfile_info('v-indels.csv', is_data=False)
# v_truehist = plotting.get_cluster_size_hist(utils.get_true_partition(reco_info).values())
# args.infnames = ['v-indels-partitions-new.csv', ]
# v_cplot = ClusterPlot(args)

# input_info, reco_info = seqfileopener.get_seqfile_info('cdr3-indels.csv', is_data=False)
# cdr3_truehist = plotting.get_cluster_size_hist(utils.get_true_partition(reco_info).values())
# args.infnames = ['cdr3-indels-partitions-new.csv', ]
# cdr3_cplot = ClusterPlot(args)
# print 'v-indels: %f' % v_cplot.adj_mi_at_max_logprob
# print 'cdr3-indels: %f' % cdr3_cplot.adj_mi_at_max_logprob
# plotting.plot_cluster_size_hists(os.getenv('www') + '/partis/tmp/foo.svg',
#                                  OrderedDict([
#                                      ['v-true', v_truehist],
#                                      ['v-indels', v_cplot.tmp_cluster_size_hist],
#                                      ['cdr3-true', cdr3_truehist],
#                                      ['cdr3-indels', cdr3_cplot.tmp_cluster_size_hist]
#                                  ]),
#                                  title='%d leaves, %dx mutation, indels' % (10, 1), xmax=10*3.01)
# sys.exit()
# # ----------------------------------------------------------------------------------------

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [ datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    # files = []
    # datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    # for datadir in datadirs:
    #     files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces
    files = [ fsdir.replace('_output', 'data') + '/' + args.dataset + '/' + human + '/shuffled.csv' for human in humans[args.dataset]]

# ----------------------------------------------------------------------------------------
def get_title(label, n_leaves, mut_mult):
    if args.data:
        title = 'data (%s %s)' % (args.dataset, label)
    else:
        title = '%d leaves, %dx mutation' % (n_leaves, mut_mult)
        if args.indels:
            title += ', indels'
    return title

# ----------------------------------------------------------------------------------------
def leafmutstr(n_leaves, mut_mult):
    return_str = 'simu-' + str(n_leaves) + '-leaves-' + str(mut_mult) + '-mutate'
    if args.indels:
        return_str += '-indels'
    if args.lonely_leaves:
        return_str += '-lonely-leaves'
    return return_str

# ----------------------------------------------------------------------------------------
def get_outdirname(label, no_subset=False):
    outdirname = fsdir + '/' + label
    if not no_subset:
        if args.subset is not None:
            outdirname += '/subset-' + str(args.subset)
        if args.istartstop is not None:
            outdirname += '/istartstop-' + '-'.join([str(i) for i in args.istartstop])
    return outdirname

# ----------------------------------------------------------------------------------------
def get_simfname(label, n_leaves, mut_mult, no_subset=False):
    return get_outdirname(label, no_subset=no_subset) + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'

# ----------------------------------------------------------------------------------------
def get_mixcr_outdir(label, n_leaves, mut_mult):
    basedir = '/fh/fast/matsen_e/dralph/work/mixcr/' + label
    if args.data:
        outdir = basedir + '/data'
    else:
        outdir = basedir + '/' + leafmutstr(n_leaves, mut_mult)
    if args.subset is not None:
        outdir += '/subset-' + str(args.subset)
    if args.istartstop is not None:
        outdir += '/istartstop-' + '-'.join([str(i) for i in args.istartstop])
    return outdir

# ----------------------------------------------------------------------------------------
def get_changeo_outdir(label, n_leaves, mut_mult):
    if args.bak:
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo.bak/' + label
    else:
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo/' + label
    if args.data:
        imgtdir = changeo_fsdir + '/data'
    else:    
        imgtdir = changeo_fsdir + '/' + leafmutstr(n_leaves, mut_mult).replace('-', '_')
    return imgtdir

# ----------------------------------------------------------------------------------------
def generate_incorrect_partition(true_partition, n_misassigned, error_type, debug=False):
    new_partition = list(true_partition)
    if debug:
        print '  before', new_partition
    for _ in range(n_misassigned):
        iclust = random.randint(0, len(new_partition) - 1)  # choose a cluster from which to remove a sequence
        iseq = random.randint(0, len(new_partition[iclust]) - 1)  # and choose the sequence to remove from this cluster
        uid = new_partition[iclust][iseq]
        new_partition[iclust].remove(uid)  # remove it
        if [] in new_partition:
            new_partition.remove([])
        if error_type == 'singletons':  # put the sequence in a cluster by itself
            new_partition.append([uid, ])
            if debug:
                print '    %s: %d --> singleton' % (uid, iclust)
        elif error_type == 'reassign':  # choose a different cluster to add it to
            inewclust = iclust
            while inewclust == iclust:
                inewclust = random.randint(0, len(new_partition) - 1)
            new_partition[inewclust].append(uid)
            if debug:
                print '    %s: %d --> %d' % (uid, iclust, inewclust)
        else:
            raise Exception('%s not among %s' % (error_type, 'singletons, reassign'))
    if debug:
        print '  after', new_partition
    return new_partition

# ----------------------------------------------------------------------------------------
def write_latex_table(adj_mis):
    for mut_mult in args.mutation_multipliers:
        for n_leaves in args.n_leaf_list:
            if n_leaves == args.n_leaf_list[0]:
                print '\\textbf{multiplier} & \\textbf{program}  ',
            print ' & %d  ' % n_leaves,
        print '\\\\'
        print '\\hline'
        iname = 0
        for name in adj_mis[args.n_leaf_list[0]][args.mutation_multipliers[0]]:
            if name == 'vollmers-0.5':
                continue
            if iname == 0:
                print '\\multirow{3}{*}{$\\times %d$} & ' % mut_mult,
            else:
                print '& ',
            print '%25s' % legends.get(name, name),
            iname += 1
            for n_leaves in args.n_leaf_list:
                val, err = adj_mis[n_leaves][mut_mult][name]
                print '  &    %5.2f $\\pm$ %.2f' % (val, err),
            print '\\\\'

# ----------------------------------------------------------------------------------------
def parse_vollmers(these_hists, these_adj_mis, seqfname, outdir, reco_info):
    vollmers_fname = seqfname.replace('.csv', '-run-viterbi.csv')
    with open(vollmers_fname) as vfile:
        vreader = csv.DictReader(vfile)
        for line in vreader:
            vhist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['clusters']))
            histfname = outdir + '/hists/vollmers-'  + line['threshold'] + '.csv'
            vhist.write(histfname)
            these_hists['vollmers-' + line['threshold']] = vhist
            if not args.data:
                these_adj_mis['vollmers-' + line['threshold']] = float(line['adj_mi']), -1.
                print '    %s: %f' % ('vollmers-' + line['threshold'], float(line['adj_mi']))
                write_adj_mi(float(line['adj_mi']), outdir + '/adj_mis/' + os.path.basename(histfname))

                # truehist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['true_clusters']))  # true partition is also written here, for lack of a better place (note that it's of course the same for all thresholds)
                truehist = plotting.get_cluster_size_hist(utils.get_true_partition(reco_info).values())
                truehist.write(outdir + '/hists/true.csv')  # will overwite itself a few times
                these_hists['true'] = truehist

# ----------------------------------------------------------------------------------------
def parse_changeo(label, n_leaves, mut_mult, these_hists, these_adj_mis, simfname, simfbase, outdir, reco_info):
    indir = get_changeo_outdir(label, n_leaves, mut_mult)  #fsdir.replace('/partis-dev/_output', '/changeo')
    if args.data:
        fbase = 'data'
    else:
        fbase = simfbase.replace('-', '_')
    if args.bak:
        infname = indir + '/' + fbase + changeorandomcrapstr
    else:
        # infname = indir + '/' + fbase + '/' + changeorandomcrapstr
        infname = indir + '/' + changeorandomcrapstr
    if args.subset is not None:
        infname = infname.replace(changeorandomcrapstr, 'subset-' + str(args.subset) + changeorandomcrapstr)
    if args.istartstop is not None:  # TODO not yet functional
        infname = infname.replace(changeorandomcrapstr, 'istartstop_' + '_'.join([str(i) for i in args.istartstop]) + changeorandomcrapstr)
    
    id_clusters = {}  # map from cluster id to list of seq ids
    with open(infname) as chfile:
        reader = csv.DictReader(chfile, delimiter='\t')
        for line in reader:
            clid = line['CLONE']
            uid = line['SEQUENCE_ID']
            if clid not in id_clusters:
                id_clusters[clid] = []
            id_clusters[clid].append(uid)
    
    partition = [ids for ids in id_clusters.values()]
    these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
    these_hists['changeo'].write(outdir + '/hists/changeo.csv')
    if not args.data:
        adj_mi_fname = infname.replace(changeorandomcrapstr, '-adj_mi.csv')
        check_call(['cp', adj_mi_fname, outdir + '/adj_mis/changeo.csv'])
        these_adj_mis['changeo'] = (read_adj_mi(adj_mi_fname), -1.)
        print '    %s: %f' % ('changeo', these_adj_mis['changeo'][0])

# ----------------------------------------------------------------------------------------
def parse_mixcr(these_hists, these_adj_mis, seqfname, outdir, reco_info):
    mixfname = seqfname.replace('.csv', '-mixcr.tsv')
    cluster_size_list = []  # put 'em in a list first so we know where to put the hist limits
    max_cluster_size = 1
    with open(mixfname) as mixfile:
        reader = csv.DictReader(mixfile, delimiter='\t')
        for line in reader:
            csize = int(line['Clone count'])
            cluster_size_list.append(csize)
            if csize > max_cluster_size:
                max_cluster_size = csize
    mixhist = Hist(max_cluster_size, 0.5, max_cluster_size + 0.5)
    for csize in cluster_size_list:
        mixhist.fill(csize)
    these_hists['mixcr'] = mixhist
    mixhist.write(outdir + '/hists/mixcr.csv')
    if not args.data:
        these_adj_mis['mixcr'] = -1., -1.

# ----------------------------------------------------------------------------------------
def parse_partis(action, these_hists, these_adj_mis, seqfname, outdir):
    args.infnames = [seqfname.replace('.csv', '-' + action + '.csv'), ]  # NOTE make sure not to add any args here that conflict with the real command line args
    args.is_data = args.data
    args.use_all_steps = False
    args.normalize_axes = []
    args.debug = False
    args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
    cplot = ClusterPlot(args)
    cplot.tmp_cluster_size_hist.write(outdir + '/hists/' + action + '.csv')
    these_hists[action + ' partis'] = cplot.tmp_cluster_size_hist
    if not args.data:
        these_adj_mis[action + ' partis'] = cplot.adj_mi_at_max_logprob, -1.
        print '    %s: %f' % (action + ' partis', cplot.adj_mi_at_max_logprob)
        write_adj_mi(cplot.adj_mi_at_max_logprob, outdir + '/adj_mis/' + action + '.csv')

# ----------------------------------------------------------------------------------------
def write_adj_mi(adj_mi, fname):
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, 'w') as outfile:
        writer = csv.DictWriter(outfile, ['adj_mi'])
        writer.writerow({'adj_mi' : adj_mi})

# ----------------------------------------------------------------------------------------
def read_adj_mi(fname):
    """ read file with a single number """
    with open(fname) as infile:
        reader = csv.DictReader(infile, fieldnames=['adj_mi'])
        for line in reader:
            return float(line['adj_mi'])
    raise Exception('something wrong with %s file' % fname)

# ----------------------------------------------------------------------------------------
def make_distance_plots(cachefname, reco_info):
    cachevals = {}
    singletons = []
    with open(cachefname) as cachefile:
        reader = csv.DictReader(cachefile)
        for line in reader:
            unique_ids = line['unique_ids'].split(':')
            cachevals[line['unique_ids']] = float(line['logprob'])
            if len(unique_ids) == 1:
                singletons.append(unique_ids[0])

    def get_joint_key(u1, u2):
        jk1, jk2 = u1 + ':' + u2, u2 + ':' + u1
        if jk1 in cachevals:
            return jk1
        elif jk2 in cachevals:
            return jk2
        else:
            return None

    def find_nearest_clonemate(uid):
        """ nearest is largest logprob, dinglebumpkus """
        maxval, maxmate = None, None
        for uidstr, val in cachevals.items():
            uids = uidstr.split(':')
            if len(uids) != 2:
                continue
            if uid not in uids:
                continue
            if maxval is None or val > maxval:
                maxval = val
                if uids[0] == uid:
                    maxmate = uids[1]
                else:
                    maxmate = uids[0]
        return maxmate

    nbins, xmin, xmax = 30, -20, 50
    nearest_clones = Hist(nbins, xmin, xmax)
    for uid in singletons:
        mate = find_nearest_clonemate(uid)
        if mate is None:
            continue
        jk = get_joint_key(uid, mate)
        if jk is None:
            continue
        lratio = cachevals[jk] - cachevals[uid] - cachevals[mate]
        nearest_clones.fill(lratio)

    hclones = Hist(nbins, xmin, xmax)
    hnot = Hist(nbins, xmin, xmax)
    for st_a, st_b in itertools.combinations(singletons, 2):
        jk1, jk2 = st_a + ':' + st_b, st_b + ':' + st_a
        if jk1 in cachevals:
            jk = jk1
        elif jk2 in cachevals:
            jk = jk2
        else:
            continue
        lratio = cachevals[jk] - cachevals[st_a] - cachevals[st_b]
        # print '%f - %f - %f = %f' % (cachevals[jk], cachevals[st_a], cachevals[st_b], lratio)
        if utils.from_same_event(args.data, reco_info, [st_a, st_b]):
            hclones.fill(lratio)
        else:
            hnot.fill(lratio)

    fig, ax = plotting.mpl_init()
    plots = {}
    nearest_clones.normalize()
    hclones.normalize()
    hnot.normalize()
    plots['clonal'] = ax.plot(hclones.get_bin_centers()[1:-1], hclones.bin_contents[1:-1], label='clonal')
    plots['not'] = ax.plot(hnot.get_bin_centers()[1:-1], hnot.bin_contents[1:-1], label='not')
    plots['nearest_clones'] = ax.plot(nearest_clones.get_bin_centers()[1:-1], nearest_clones.bin_contents[1:-1], label='nearest clones')
    plotting.mpl_finish(ax, os.getenv('www') + '/partis/tmp', 'foop')

# ----------------------------------------------------------------------------------------
def write_all_plot_csvs(label):
    hists, adj_mis = {}, {}
    for n_leaves in args.n_leaf_list:
        for mut_mult in args.mutation_multipliers:
            write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis)

    if not args.data:
        write_latex_table(adj_mis)

# ----------------------------------------------------------------------------------------
def write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = OrderedDict()
        adj_mis[n_leaves][mut_mult] = OrderedDict()
    these_hists = hists[n_leaves][mut_mult]
    these_adj_mis = adj_mis[n_leaves][mut_mult]

    plotdir = os.getenv('www') + '/partis/clustering/' + label
    if args.data:
        seqfname = get_simfname(label, n_leaves, mut_mult).replace(leafmutstr(n_leaves, mut_mult), 'data')  # hackey hackey hackey
        simfbase = None
        csvdir = os.path.dirname(seqfname) + '/data'
        plotfname = plotdir + '/plots/data.svg'
        title = get_title(label, n_leaves, mut_mult)
    else:
        seqfname = get_simfname(label, n_leaves, mut_mult)
        simfbase = leafmutstr(n_leaves, mut_mult)
        csvdir = os.path.dirname(seqfname) + '/' + simfbase
        plotfname = plotdir + '/plots/' + simfbase + '.svg'
        title = get_title(label, n_leaves, mut_mult)

    input_info, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=args.data)
    if args.count_distances:
        make_distance_plots(seqfname.replace('.csv', '-partition-cache.csv'), reco_info)
        return

    # then vollmers annotation (and true hists)
    parse_vollmers(these_hists, these_adj_mis, seqfname, csvdir, reco_info)

    # mixcr
    # parse_mixcr(these_hists, these_adj_mis, seqfname, csvdir, reco_info)

    # # then changeo
    # parse_changeo(label, n_leaves, mut_mult, these_hists, these_adj_mis, seqfname, simfbase, csvdir, reco_info)

    # partis stuff
    # for ptype in ['vsearch-', 'naive-hamming-', '']:
    for ptype in ['', ]:
        parse_partis(ptype + 'partition', these_hists, these_adj_mis, seqfname, csvdir)

    plotting.plot_cluster_size_hists(plotfname, these_hists, title=title, xmax=n_leaves*3.01)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

    # for name, adj_mi in these_adj_mis.items():
    #     with open(os.path.dirname(seqfname) + '/adj_mis/' + name.replace(' ', '_') + '.csv', 'w') as adjmifile:
    #         adjmifile.write(adj_mi + '\n')

# ----------------------------------------------------------------------------------------
def compare_all_subsets(label):
    print '\n\nfigure out why all the adj mis are screwed up\n\n'
    hists, adj_mis = {}, {}
    for n_leaves in args.n_leaf_list:
        for mut_mult in args.mutation_multipliers:
            compare_each_subsets(label, n_leaves, mut_mult, hists, adj_mis)

    if not args.data:
        write_latex_table(adj_mis)

# ----------------------------------------------------------------------------------------
def compare_each_subsets(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = OrderedDict()
        adj_mis[n_leaves][mut_mult] = OrderedDict()
    these_hists = hists[n_leaves][mut_mult]
    these_adj_mis = adj_mis[n_leaves][mut_mult]

    basedir = fsdir + '/' + label
    # expected_methods = ['vollmers-0.9', 'mixcr', 'changeo', 'vsearch-partition', 'naive-hamming-partition', 'partition']  # mostly so we can specify the order
    expected_methods = ['vollmers-0.9', 'mixcr', 'vsearch-partition', 'naive-hamming-partition', 'partition']  # mostly so we can specify the order
    if not args.data:
        expected_methods.insert(0, 'true')
    # expected_methods = ['vollmers-0.9', 'vsearch-partition', 'naive-hamming-partition']
    tmp_adj_mis = OrderedDict()
    for method in expected_methods:
        method_hists = []
        if args.n_subsets is not None:
            subdirs = [basedir + '/subset-' + str(isub) for isub in range(args.n_subsets)]
        elif args.startstoplist is not None:
            subdirs = [basedir + '/istartstop-' + istartstop.replace(',', '-') for istartstop in args.startstoplist]
        else:
            assert False
        for subdir in subdirs:
            # hists
            if args.data:
                histfname = subdir + '/data/hists/' + method + '.csv'
            else:
                histfname = subdir + '/' + leafmutstr(n_leaves, mut_mult) + '/hists/' + method + '.csv'
            method_hists.append(Hist(fname=histfname))
            if method == 'true':
                continue

            if not args.data:
                # adj mis
                fname = subdir + '/' + leafmutstr(n_leaves, mut_mult) + '/adj_mis/' + method + '.csv'
                if method not in tmp_adj_mis:
                    tmp_adj_mis[method] = []
                if method == 'mixcr':
                    tmp_adj_mis[method].append(-1.)
                else:
                    tmp_adj_mis[method].append(read_adj_mi(fname))

        these_hists[method] = plotting.make_mean_hist(method_hists)

    plotdir = os.getenv('www') + '/partis/clustering/subsets/' + label
    if args.data:
        title = get_title(label, n_leaves, mut_mult)
        plotfname = plotdir + '/plots/data.svg'
        xmax = 10
    else:
        title = get_title(label, n_leaves, mut_mult)
        plotfname = plotdir + '/plots/' + leafmutstr(n_leaves, mut_mult) + '.svg'
        xmax = n_leaves*3.01
    plotting.plot_cluster_size_hists(plotfname, these_hists, title=title, xmax=xmax)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])
    if not args.data:
        for meth, vals in tmp_adj_mis.items():
            mean = numpy.mean(vals)
            std = numpy.std(vals)
            these_adj_mis[meth] = (mean, std)
            print '  %30s %.3f +/- %.3f' % (meth, mean, std)

# ----------------------------------------------------------------------------------------
def get_misassigned_adj_mis(simfname, misassign_fraction, nseq_list, error_type):
    input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
    n_reps = 1
    uid_list = input_info.keys()
    new_partitions = {}
    for nseqs in nseq_list:
        for irep in range(n_reps):  # repeat <nreps> times
            istart = irep * nseqs
            istop = istart + nseqs
            uids = uid_list[istart : istop]
            true_partition = utils.get_true_clusters(uids, reco_info).values()
            n_misassigned = int(misassign_fraction * nseqs)
            new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type=error_type)
            # new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type='singletons')
            new_partitions[nseqs] = new_partition
    return {nseqs : utils.mutual_information(new_partitions[nseqs], reco_info) for nseqs in nseq_list}

# ----------------------------------------------------------------------------------------
def make_adj_mi_vs_sample_size_plot(label, n_leaves, mut_mult, nseq_list, adj_mis):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    fsize = 20
    mpl.rcParams.update({
        # 'font.size': fsize,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })
    fig, ax = plt.subplots()
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

    plots = {}

    colors = {'vsearch-partition' : '#008b8b',
              'naive-hamming-partition' : '#ff8c00',
              'partition' : '#cc0000',
              '0.1-true-singletons' : '#006600',
              '0.1-true-reassign' : '#32cd32'}
    adj_mis['0.1-true-singletons'] = get_misassigned_adj_mis(get_simfname(label, n_leaves, mut_mult), 0.1, nseq_list, 'singletons')
    adj_mis['0.1-true-reassign'] = get_misassigned_adj_mis(get_simfname(label, n_leaves, mut_mult), 0.1, nseq_list, 'reassign')
    for meth in adj_mis:
        linewidth = 2
        linestyle = '-'
        if 'true' in meth:
            linewidth = 4
            linestyle = '--'
        plots[meth] = ax.plot(nseq_list, [adj_mis[meth][ns] for ns in nseq_list], linewidth=linewidth, label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=1)
        if 'true' not in meth:
            plt.scatter(nseq_list, [adj_mis[meth][ns] for ns in nseq_list], color=colors.get(meth, 'grey'), linestyle='-', alpha=1, s=[30 for _ in range(len(nseq_list))])

    legend = ax.legend(loc=(0.55, 0.1))
    sns.despine(trim=True, bottom=True)
    sns.set_style("ticks")
    plt.title('%d leaves, %dx mutation' % (n_leaves, mut_mult))
    plt.xlabel('sample size')
    plt.ylabel('adjusted MI')
    plt.subplots_adjust(bottom=0.14, left=0.14)
    ax.set_xscale('log')
    xmin, xmax = nseq_list[0], nseq_list[-1] + 10
    plt.xlim(xmin, xmax)
    plt.ylim(0, 1.08)
    potential_xticks = [5, 10, 25, 50, 100, 300, 1000]
    xticks = [xt for xt in potential_xticks if xt < xmax]
    plt.xticks(xticks, [str(xt) for xt in xticks])
    plotdir = os.getenv('www') + '/partis/clustering/sample-sizes/' + label
    outfname = plotdir + '/plots/sample-sizes.svg'
    if not os.path.exists(plotdir + '/plots'):
        os.makedirs(plotdir + '/plots')
    plt.savefig(outfname)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

    return

# ----------------------------------------------------------------------------------------
def compare_sample_sizes(label, n_leaves, mut_mult):
    """ see sample-size-partition.py """
    raise Exception('It is possible that I may have broken this')
    expected_methods = ['vsearch-partition', 'naive-hamming-partition', 'partition']  # mostly so we can specify the order
    # expected_methods = ['vollmers-0.9', 'partition']  # mostly so we can specify the order
    basedir = fsdir + '/' + label
    nseq_list = []
    adj_mis, hists = OrderedDict(), OrderedDict()
    for meth in expected_methods:
        adj_mis[meth], hists = {}, {}
    # for dirname in glob.glob(basedir + '/istartstop-*'):
    #     dumstr, istart, istop = dirname.split('/')[-1].split('-')
    #     istart, istop = int(istart), int(istop)
    for istartstop in args.startstoplist:
        istart, istop = [int(i) for i in istartstop.split(',')]
        nseqs = istop - istart
        dirname = basedir + '/istartstop-' + str(istart) + '-' + str(istop)
        nseq_list.append(nseqs)
        for meth in expected_methods:
            action = meth
            if 'vollmers' in meth:
                action = 'run-viterbi'
            if args.data:
                csvfname = dirname + '/data' + '-' + action + '.csv'
            else:
                csvfname = dirname + '/' + leafmutstr(n_leaves, mut_mult) + '-' + action + '.csv'
            args.infnames = [csvfname, ]  # NOTE make sure not to add any args here that conflict with the real command line args
            args.is_data = args.data
            args.use_all_steps = False
            args.normalize_axes = []
            args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
            args.debug = False
            cplot = ClusterPlot(args)
            # cplot.tmp_cluster_size_hist.write(outdir + '/hists/' + action + '.csv')
            if args.data:
                adj_mis[meth][nseqs] = cplot.adj_mi_at_max_logprob
            hists[meth][nseqs] = cplot.tmp_cluster_size_hist

    nseq_list.sort()

    if not args.data:
        make_adj_mi_vs_sample_size_plot(label, n_leaves, mut_mult, nseq_list, adj_mis)

    mean_hists = OrderedDict()
    for meth in expected_methods:
        mean_hists[meth] = plotting.make_mean_hist(hists[meth])
    plotdir = os.getenv('www') + '/partis/clustering/' + label
    if args.data:
        plotfname = plotdir + '/plots/data.svg'
        title = get_title(label, n_leaves, mut_mult)
        xmax = None
    else:
        plotfname = plotdir + '/plots/' + simfbase + '.svg'
        title = get_title(label, n_leaves, mut_mult)
        xmax = n_leaves*3.01
    plotting.plot_cluster_size_hists(plotfname, mean_hists, title=title, xmax=xmax)

# ----------------------------------------------------------------------------------------
def get_seqfile(action, label, datafname, n_leaves=None, mut_mult=None):

    def slice_csv(csv_infname, csv_outfname):
        if os.path.exists(csv_outfname):
            print '      slicefile exists %s' % csv_outfname
            return
        print '      subsetting %d seqs with indices %d --> %d' % (args.istartstop[1] - args.istartstop[0], args.istartstop[0], args.istartstop[1])
        if not os.path.exists(os.path.dirname(csv_outfname)):
            os.makedirs(os.path.dirname(csv_outfname))
        remove_csv_infname = False
        # if '.bz2' in csv_infname:  # put ~10 more lines than we need into a tmp file that isn't bzipped
        #     tmpfname = '.tmp'.join(os.path.splitext(csv_outfname))
        #     assert '.tsv' in csv_infname  # aw, screw it, just replace tabs with commas
        #     assert '.csv' in csv_outfname
        #     check_call('bzgrep -m' + str(args.istartstop[1] + 10) + ' . ' + csv_infname + ' | sed \'s/\t/,/g\' >' + tmpfname, shell=True)
        #     csv_infname = tmpfname
        #     remove_csv_infname = True
        # check_call('head -n1 ' + csv_infname + ' | sed -e s/name/unique_id/ -e s/nucleotide/seq/ >' + csv_outfname, shell=True)
        check_call('head -n1 ' + csv_infname + ' >' + csv_outfname, shell=True)
        check_call('sed -n \'' + str(args.istartstop[0] + 2) + ',' + str(args.istartstop[1] + 1) + ' p\' ' + csv_infname + '>>' + csv_outfname, shell=True)  # NOTE conversion from standard zero indexing to sed inclusive one-indexing (and +1 for header line)
        if remove_csv_infname:
            assert '/dralph/' in csv_infname
            os.remove(csv_infname)

    if args.data:
        if args.istartstop is not None:
            subfname = fsdir + '/' + label + '/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/data.csv'
            slice_csv(datafname, subfname)
            datafname = subfname

        seqfname = datafname
    else:
        assert n_leaves is not None and mut_mult is not None
        simfname = get_simfname(label, n_leaves, mut_mult, no_subset=True)

        if args.subset is not None:
            ntot = int(check_output(['wc', '-l', simfname]).split()[0]) - 1
            subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
            if os.path.exists(subsimfname):
                print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d / %d' % (args.subset, args.n_subsets)
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                n_per_subset = int(float(ntot) / args.n_subsets)  # NOTE ignores remainders, i.e. last few sequences
                istart = args.subset * n_per_subset + 2  # NOTE sed indexing (one-indexed with inclusive bounds). Also note extra +1 to avoid header
                istop = istart + n_per_subset - 1
                check_call('sed -n \'' + str(istart) + ',' + str(istop) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)
            simfname = subsimfname

        if args.istartstop is not None:
            subsimfname = simfname.replace(label + '/', label + '/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
            slice_csv(simfname, subsimfname)
            simfname = subsimfname

        seqfname = simfname

    return seqfname

# ----------------------------------------------------------------------------------------
def execute(action, label, datafname, n_leaves=None, mut_mult=None):
    real_action = action
    if 'partition' in action:
        real_action = 'partition'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + real_action

    extras = []
    seqfname = get_seqfile(action, label, datafname, n_leaves, mut_mult)
    if args.data:
        cmd += ' --datafname ' + seqfname + ' --is-data'
        if args.dataset == 'adaptive':
            extras += ['--skip-unproductive', ]
    else:
        cmd += ' --simfname ' + seqfname

    def output_exists(outfname):
        if os.path.exists(outfname):
            if args.overwrite:
                print '                      overwriting %s' % outfname
                return False
            else:
                print '                      output exists, skipping (%s)' % outfname
                return True
        else:
            return False

    def get_outputname():
        if args.data:
            return get_outdirname(label) + '/data-' + action + '.csv'
        else:
            return ('-' + action).join(os.path.splitext(seqfname))

    if action == 'cache-data-parameters':
        if output_exists(fsdir + '/' + label + '/data'):
            return
        extras += ['--n-max-queries', + args.n_data_to_cache]
        n_procs = args.n_data_to_cache / 500
    elif action == 'simulate':
        if output_exists(seqfname):
            return
        extras += ['--n-sim-events', int(float(args.n_sim_seqs) / n_leaves)]
        extras += ['--n-leaves', n_leaves, '--mutation-multiplier', mut_mult]
        if args.indels:
            extras += ['--indel-frequency', 0.5]
        if args.lonely_leaves:
            extras += ['--constant-number-of-leaves', ]
        n_procs = 10
    elif action == 'cache-simu-parameters':
        if output_exists(seqfname.replace('.csv', '')):
            return
        n_procs = 20
    elif action == 'partition':
        outfname = get_outputname()
        if output_exists(outfname):
            return
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition]
        if args.count_distances:
            extras += ['--persistent-cachefname', ('-cache').join(os.path.splitext(outfname))]  # '--n-partition-steps', 1, 
        n_procs = max(1, args.n_to_partition / 15)  # something like 15 seqs/process to start with
    elif action == 'naive-hamming-partition':
        outfname = get_outputname()
        if output_exists(outfname):
            return
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition, '--auto-hamming-fraction-bounds']
        n_procs = max(1, args.n_to_partition / 30)
    elif action == 'vsearch-partition':
        outfname = get_outputname()
        # outfname = '-vsearch-partition'.join(os.path.splitext(seqfname))
        if output_exists(outfname):
            return
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition, '--naive-vsearch']
        n_procs = max(1, args.n_to_partition / 100)  # only used for ighutil step
    elif action == 'run-viterbi':
        outfname = get_outputname()
        # outfname = '-run-viterbi'.join(os.path.splitext(seqfname))
        if output_exists(outfname):
            return
        cmd += ' --outfname ' + outfname
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        extras += ['--n-max-queries', args.n_to_partition]
        n_procs = max(1, args.n_to_partition / 50)
    elif action == 'run-changeo':

        def untar_imgt(imgtdir):
            tar_cmd = 'mkdir ' + imgtdir + ';'
            tar_cmd += ' tar Jxvf ' + imgtdir + '.txz --exclude=\'IMGT_HighV-QUEST_individual_files_folder/*\' -C ' + imgtdir
            check_call(tar_cmd, shell=True)

        imgtdir = get_changeo_outdir(label, n_leaves, mut_mult)
        if os.path.isdir(imgtdir):
            print '                      already untar\'d into %s' % imgtdir
        else:
            if os.path.exists(imgtdir + '.txz'):
                untar_imgt(imgtdir)
            else:
                print '   hmm... imgtdir not there... maybe we only have the subsets'

        if args.subset is not None:
            subset_dir = imgtdir + '/subset-' + str(args.subset)
            if not os.path.exists(subset_dir):
                os.makedirs(subset_dir)
                tsvfnames = glob.glob(imgtdir + '/*.txt')
                check_call(['cp', '-v', imgtdir + '/11_Parameters.txt', subset_dir + '/'])
                tsvfnames.remove(imgtdir + '/11_Parameters.txt')
                tsvfnames.remove(imgtdir + '/README.txt')
                input_info, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=False)
                subset_ids = input_info.keys()
                utils.subset_files(subset_ids, tsvfnames, subset_dir)
            imgtdir = subset_dir
        if args.istartstop is not None:
            subset_dir = imgtdir + '/istartstop_' + '_'.join([str(i) for i in args.istartstop])
            if os.path.exists(subset_dir + '.txz'):
                untar_imgt(subset_dir)
            elif not os.path.exists(subset_dir):
                os.makedirs(subset_dir)
                tsvfnames = glob.glob(imgtdir + '/*.txt')
                check_call(['cp', '-v', imgtdir + '/11_Parameters.txt', subset_dir + '/'])
                tsvfnames.remove(imgtdir + '/11_Parameters.txt')
                tsvfnames.remove(imgtdir + '/README.txt')
                input_info, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=args.data)
                subset_ids = input_info.keys()
                utils.subset_files(subset_ids, tsvfnames, subset_dir)
            imgtdir = subset_dir

        def run(cmdstr):
            print 'RUN %s' % cmdstr
            check_call(cmdstr.split())
        
        resultfname = imgtdir + changeorandomcrapstr
        if os.path.exists(resultfname):
            print '                         changeo already finished (%s)' % resultfname
            return

        # check_call(['./bin/csv2fasta', seqfname])
        fastafname = os.path.splitext(seqfname)[0] + '.fasta'
        utils.csv_to_fasta(seqfname, outfname=fastafname)  #, name_column='name' if args.data else 'unique_id', seq_column='nucleotide' if args.data else 'seq')
        bindir = '/home/dralph/work/changeo/changeo'
        start = time.time()
        # check_call(['mv', seqfname.replace('.csv', '.fa'), seqfname.replace('.csv', '.fasta')])
        cmd = bindir + '/MakeDb.py imgt -i ' + imgtdir + ' -s ' + fastafname
        run(cmd)
        cmd = bindir + '/ParseDb.py select -d ' + imgtdir + '_db-pass.tab -f FUNCTIONAL -u T'
        run(cmd)
        cmd = bindir + '/DefineClones.py bygroup -d ' + imgtdir + '_db-pass_parse-select.tab --act first --model m1n --dist 7'
        run(cmd)
        print '        changeo time: %.3f' % (time.time()-start)

        # read changeo's output and toss it into a csv
        input_info, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=args.data)
        id_clusters = {}  # map from cluster id to list of seq ids
        with open(resultfname) as chfile:
            reader = csv.DictReader(chfile, delimiter='\t')
            for line in reader:
                clid = line['CLONE']
                uid = line['SEQUENCE_ID']
                if clid not in id_clusters:
                    id_clusters[clid] = []
                id_clusters[clid].append(uid)
        
        partition = [ids for ids in id_clusters.values()]
        # these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
        if not args.data:
            write_adj_mi(utils.mutual_information(partition, reco_info, debug=True), imgtdir + '-adj_mi.csv')
        return
    elif action == 'run-mixcr':
        binary = '/home/dralph/work/mixcr/mixcr-1.2/mixcr'
        mixcr_workdir = get_mixcr_outdir(label, n_leaves, mut_mult)
        if not os.path.exists(mixcr_workdir):
            os.makedirs(mixcr_workdir)

        # fastafname = os.path.splitext(seqfname)[0] + '.fasta'
        infname = mixcr_workdir + '/' + os.path.basename(os.path.splitext(seqfname)[0] + '.fasta')
        outfname = os.path.splitext(seqfname)[0] + '-mixcr.tsv'
        if os.path.exists(outfname):
            print '                      mixcr output exists, skipping (%s)' % outfname
            return

        # check_call(['./bin/csv2fasta', seqfname])
        utils.csv_to_fasta(seqfname, outfname=infname, n_max_lines=args.n_to_partition)  #, name_column='name' if args.data else 'unique_id', seq_column='nucleotide' if args.data else 'seq'
        # check_call('head -n' + str(2*args.n_to_partition) + ' ' + fastafname + ' >' + infname, shell=True)
        # os.remove(seqfname.replace('.csv', '.fa'))

        def run(cmdstr):
            print 'RUN %s' % cmdstr
            check_call(cmdstr.split())

        start = time.time()
        cmd = binary + ' align -f --loci IGH ' + infname + ' ' + infname.replace('.fasta', '.vdjca')
        run(cmd)
        cmd = binary + ' assemble -f ' + infname.replace('.fasta', '.vdjca') + ' ' + infname.replace('.fasta', '.clns')
        run(cmd)
        cmd = binary + ' exportClones ' + infname.replace('.fasta', '.clns') + ' ' + infname.replace('.fasta', '.txt')
        run(cmd)
        print '        mixcr time: %.3f' % (time.time()-start)
        check_call(['cp', '-v', infname.replace('.fasta', '.txt'), outfname])

        return
    elif action == 'compare-sample-sizes':
        compare_sample_sizes(label, n_leaves, mut_mult)
        return
    elif action == 'compare-subsets':
        assert False
    else:
        raise Exception('bad action %s' % action)

    cmd +=  ' --plotdir ' + os.getenv('www') + '/partis'
    n_proc_str = str(n_procs)
    if n_procs > 10:
        n_proc_str += ':10'
        extras += ['--slurm', '--workdir', fsdir + '/_tmp/' + str(random.randint(0,99999))]
        
    extras += ['--n-procs', n_proc_str]

    cmd += utils.get_extra_str(extras)
    print '   ' + cmd
    # return
    # check_call(cmd.split())
    # return
    if args.data:
        logbase = fsdir + '/' + label + '/_logs/data-' + action
    else:
        logbase = fsdir + '/' + label + '/_logs/' + leafmutstr(n_leaves, mut_mult) + '-' + action
    if args.subset is not None:
        logbase = logbase.replace('_logs/', '_logs/subset-' + str(args.subset) + '/')
    if args.istartstop is not None:
        logbase = logbase.replace('_logs/', '_logs/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
    if not os.path.exists(os.path.dirname(logbase)):
        os.makedirs(os.path.dirname(logbase))
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    time.sleep(5)

# ----------------------------------------------------------------------------------------
for datafname in files:
    # if '/B/' not in datafname:
    #     continue
    if args.dataset == 'stanford':
	    human = os.path.basename(datafname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
	    human = re.findall('[ABC]', datafname)[0]
    print 'run', human
    label = human
    if args.bak:
        label += '.bak'

    if args.only_run is not None and human not in args.only_run:
        continue

    performance_fname = fsdir + '/' + label + '/partition-performance.csv'
    performance_file_header = ['method', 'mut_mult', 'n_leaves', 'adj_mi', 'n_clusters', 'n_true_clusters']

    for action in args.actions:
        if action == 'write-plots' or action == 'compare-subsets':
            continue
        print '      ----> ', action
        if action == 'cache-data-parameters':
            execute(action, label, datafname)
            continue
        for n_leaves in args.n_leaf_list:
            print '  ----> ', n_leaves, ' leaves'
            for mut_mult in args.mutation_multipliers:
                print '         ----> mutate', mut_mult
                execute(action, label, datafname, n_leaves, mut_mult)
                # sys.exit()

    if 'write-plots' in args.actions:
        write_all_plot_csvs(label)
    if 'compare-subsets' in args.actions:
        compare_all_subsets(label)

    break
