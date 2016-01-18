import os
import glob
import sys
import numpy
import random
from collections import OrderedDict
import time
import csv
from subprocess import check_call, Popen, check_output, PIPE
import itertools
sys.path.insert(1, './python')
csv.field_size_limit(sys.maxsize)
from hist import Hist
import seqfileopener
import utils
import baseutils
from humans import humans
# from clusterplot import ClusterPlot
from clusterpath import ClusterPath
import plotting

changeorandomcrapstr = '_db-pass_parse-select_clone-pass.tab'
metrics = ['adj_mi', 'ccf_under', 'ccf_over']

# ----------------------------------------------------------------------------------------
def FOOP():
    glfo = utils.read_germline_set('data/imgt')
    _, reco_info = seqfileopener.get_seqfile_info('/fh/fast/matsen_e/dralph/work/partis-dev/_output/A/simu-3-leaves-1.0-mutate.csv', is_data=False, n_max_queries=5)
    true_partition = utils.get_true_partition(reco_info)
    utils.generate_distance_based_incorrect_partition(glfo, reco_info, true_partition, hfrac_error=0.03, merge_fraction=0.1, debug=True)
    sys.exit()
    tholds = [0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.3, 0.9]
    # tholds = [0.005, 0.04, 0.9]
    adj_mis, ccf_unders, ccf_overs, nccf_a, nccf_b = [], [], [], [], []
    nccf_product = []
    for thold in tholds:
        cpath = ClusterPath()
        cpath.readfile(str(thold) + '.csv')
        partition = cpath.partitions[cpath.i_best]
        # cpath.print_partition(cpath.i_best)
        ccf = utils.correct_cluster_fractions(partition, true_partition)
        nccf = utils.new_ccfs_that_need_better_names(partition, true_partition, reco_info)
        print '%5.3f   %5.4f    %5.4f  %5.4f    %5.4f  %5.4f' % (thold, cpath.adj_mis[cpath.i_best], ccf[0], ccf[1], nccf[1], nccf[0])
        adj_mis.append(cpath.adj_mis[cpath.i_best])
        ccf_unders.append(ccf[0])
        ccf_overs.append(ccf[1])
        nccf_a.append(nccf[1])
        nccf_b.append(nccf[0])
        nccf_product.append(nccf[0] * nccf[1])

    fig, ax = plotting.mpl_init()
    ax.plot(tholds, adj_mis, label='adj mi')
    ax.plot(tholds, ccf_unders, label='ccf under', color='#cc0000')
    ax.plot(tholds, ccf_overs, label='ccf over', color='#cc0000', linestyle='--')
    ax.plot(tholds, nccf_a, label='new ccf a', color='#4e8975')
    ax.plot(tholds, nccf_b, label='new ccf b', color='#4e8975', linestyle='--')
    ax.plot(tholds, nccf_product, label='new ccf product', color='#006600')
    plotting.mpl_finish(ax, os.getenv('www') + '/partis/tmp', 'foop', log='x', xticks=tholds, xticklabels=tholds)  #, leg_loc=())
    # hist = plotting.get_cluster_size_hist(cpath.partitions[cpath.i_best])

# ----------------------------------------------------------------------------------------
def mut_str(args, label, n_leaves, mut_mult):
    # if float(mut_mult).is_integer():  # TODO fix this
    if mut_mult - int(mut_mult) == 0.:  # TODO fix this
        return '%dx mutation' % mut_mult
    else:
        return '%.1fx mutation' % mut_mult

# ----------------------------------------------------------------------------------------
def get_title(args, label, n_leaves, mut_mult, hfrac_bounds=None):
    if args.data:
        title = 'data (%s %s)' % ('fixme', label)
    else:
        title = '%d leaves, %s' % (n_leaves, mut_str(args, label, n_leaves, mut_mult))
        if hfrac_bounds is not None:
            title += ', %.2f-%.2f hfrac' % tuple(hfrac_bounds)
        if args.istartstop is not None:
            title += ', %d seqs' % (args.istartstop[1] - args.istartstop[0])
        if args.indels:
            title += ', indels'
    return title

# ----------------------------------------------------------------------------------------
def get_str(a_list, delimiter='-'):
    return delimiter.join([str(item) for item in a_list])

# ----------------------------------------------------------------------------------------
def leafmutstr(args, n_leaves, mut_mult, hfrac_bounds=None):
    return_str = 'simu-' + str(n_leaves) + '-leaves-' + str(mut_mult) + '-mutate'
    if hfrac_bounds is not None:
        return_str += '-hfrac-bounds-' + get_str(hfrac_bounds)
    if args.indels:
        return_str += '-indels'
    if args.lonely_leaves:
        return_str += '-lonely-leaves'
    if args.mimic:
        return_str += '-mimic'
    return return_str

# ----------------------------------------------------------------------------------------
def get_outdirname(args, label, no_subset=False):
    outdirname = args.fsdir + '/' + label
    if not no_subset:
        if args.subset is not None:
            outdirname += '/subset-' + str(args.subset)
        if args.istartstop is not None:
            outdirname += '/istartstop-' + get_str(args.istartstop)
    return outdirname

# ----------------------------------------------------------------------------------------
def get_simfname(args, label, n_leaves, mut_mult, no_subset=False):
    return get_outdirname(args, label, no_subset=no_subset) + '/' + leafmutstr(args, n_leaves, mut_mult) + '.csv'

# ----------------------------------------------------------------------------------------
def get_program_workdir(args, program_name, label, n_leaves, mut_mult):
    basedir = '/fh/fast/matsen_e/dralph/work/' + program_name + '/' + label
    if args.data:
        outdir = basedir + '/data'
    else:
        outdir = basedir + '/' + leafmutstr(args, n_leaves, mut_mult)
    if args.subset is not None:
        outdir += '/subset-' + str(args.subset)
    if args.istartstop is not None:
        outdir += '/istartstop-' + get_str(args.istartstop)
    return outdir

# ----------------------------------------------------------------------------------------
def get_changeo_outdir(args, label, n_leaves, mut_mult):
    if args.bak:
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo.bak/' + label
    else:
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo/' + label
    if args.data:
        imgtdir = changeo_fsdir + '/data'
    else:
        imgtdir = changeo_fsdir + '/' + leafmutstr(args, n_leaves, mut_mult).replace('-', '_')
    return imgtdir

# ----------------------------------------------------------------------------------------
def deal_with_parse_results(info, outdir, vname, partition, hist, metrics=None, info_vname=None):
    if info_vname is None:
        info_vname = vname
    if partition is not None:
        info['partitions'][info_vname] = partition
    info['hists'][info_vname] = hist
    hist.write(outdir + '/hists/' + vname + '.csv')
    if metrics is not None and vname != 'true':
        for mname, val in metrics.items():
            info[mname][info_vname] = val
            write_float_val(outdir + '/' + mname + '/' + vname + '.csv', val, mname)

# ----------------------------------------------------------------------------------------
def parse_true(args, info, outdir, true_partition):
    # well, not really parse per se
    truehist = plotting.get_cluster_size_hist(true_partition)
    deal_with_parse_results(info, outdir, 'true', true_partition, truehist, metrics=None)

# ----------------------------------------------------------------------------------------
def parse_vollmers(args, info, vollmers_fname, outdir, reco_info, true_partition):
    n_lines = 0
    with open(vollmers_fname) as vfile:
        vreader = csv.DictReader(vfile)
        for line in vreader:
            n_lines += 1
            partitionstr = line['partition'] if 'partition' in line else line['clusters']  # backwards compatibility -- used to be 'clusters' and there's still a few old files floating around
            partition = utils.get_partition_from_str(partitionstr)
            metrics = None
            if not args.data:
                utils.check_intersection_and_complement(partition, true_partition)
                ccfs = utils.new_ccfs_that_need_better_names(partition, true_partition, reco_info)
                metrics = {'adj_mi' : float(line['adj_mi']), 'ccf_under' : ccfs[0], 'ccf_over' : ccfs[1]}

            deal_with_parse_results(info, outdir, 'vollmers-' + line['threshold'], partition, plotting.get_cluster_size_hist(partition), metrics)

    if n_lines < 1:
        raise Exception('zero partition lines read from %s' % vollmers_fname)

# ----------------------------------------------------------------------------------------
def parse_changeo(args, label, n_leaves, mut_mult, info, simfbase, outdir):
    print 'TODO update for all sorts of shit (including hfrac_bounds)'
    raise Exception('rerun all the changeo stuff accounting for missing uids')
    indir = get_changeo_outdir(args, label, n_leaves, mut_mult)  #fsdir.replace('/partis-dev/_output', '/changeo')
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
    if args.istartstop is not None:
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
    metrics = None
    if not args.data:
        adj_mi_fname = infname.replace(changeorandomcrapstr, '-adj_mi.csv')
        metrics = {'adj_mi' : read_float_val(adj_mi_fname, 'adj_mi')}
        for etype in ['under', 'over']:
            ccf_fname = infname.replace(changeorandomcrapstr, '-ccf_' + etype+ '.csv')
            metrics['ccf_' + etype] = read_float_val(ccf_fname, 'ccf_' + etype)
    deal_with_parse_results(info, outdir, 'changeo', partition, plotting.get_cluster_size_hist(partition), metrics)

# ----------------------------------------------------------------------------------------
def parse_mixcr(args, info, seqfname, outdir):
    print 'TODO update for all sorts of shit (including hfrac_bounds)'
    raise Exception('rerun all the mixcr stuff accounting for missing uids')
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
    deal_with_parse_results(info, outdir, 'mixcr', None, mixhist, None)

# ----------------------------------------------------------------------------------------
def parse_partis(args, action, info, outfname, outdir, reco_info, true_partition):
    cpath = ClusterPath()
    cpath.readfile(outfname)
    hist = plotting.get_cluster_size_hist(cpath.partitions[cpath.i_best])
    partition = cpath.partitions[cpath.i_best]
    vname = action
    info_vname = vname + ' partis'  # arg, shouldn't have done it that way
    metrics = None
    if not args.data:
        # ccfs = utils.new_ccfs_that_need_better_names(cpath.partitions[cpath.i_best], true_partition, reco_info)
        ccfs = cpath.ccfs[cpath.i_best]
        metrics = {'adj_mi' : cpath.adj_mis[cpath.i_best], 'ccf_under' : ccfs[0], 'ccf_over' : ccfs[1]}
    deal_with_parse_results(info, outdir, action, partition, hist, metrics, info_vname)

# ----------------------------------------------------------------------------------------
def get_synthetic_partition_type(stype):
    misfrac, mistype, threshold = None, None, None
    if 'distance' in stype:
        mistype, threshold = stype.split('-')
        threshold = float(threshold)
    else:
        misfrac, mistype = stype.split('-')
        misfrac = float(misfrac)
    return misfrac, mistype, threshold

# ----------------------------------------------------------------------------------------
def generate_synthetic_partitions(args, label, n_leaves, mut_mult, seqfname, base_outfname):
    _, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=False)
    glfo = utils.read_germline_set(args.datadir)
    true_partition = utils.get_true_partition(reco_info)
    for stype in args.synthetic_partitions:
        misfrac, mistype, threshold = get_synthetic_partition_type(stype)
        vname = 'misassign-' + stype
        outfname = base_outfname.replace('.csv', '-' + vname + '.csv')
        if output_exists(args, outfname):
            continue
        new_partition = utils.generate_incorrect_partition(true_partition, misassign_fraction=misfrac, error_type=mistype, threshold=threshold, reco_info=reco_info, glfo=glfo)
        cpath = ClusterPath()
        adj_mi = utils.adjusted_mutual_information(true_partition, new_partition)
        ccfs = utils.new_ccfs_that_need_better_names(new_partition, true_partition, reco_info)
        cpath.add_partition(new_partition, logprob=float('-inf'), n_procs=1, adj_mi=adj_mi, ccfs=ccfs)
        cpath.write(outfname, is_data=False, reco_info=reco_info, true_partition=true_partition)

# ----------------------------------------------------------------------------------------
def parse_synthetic(args, info, outdir, true_partition, base_outfname):
    for stype in args.synthetic_partitions:
        misfrac, mistype, threshold = get_synthetic_partition_type(stype)
        assert False  # update for threshold
        vname = 'misassign-' + stype
        cpath = ClusterPath()
        outfname = base_outfname.replace('.csv', '-' + vname + '.csv')
        cpath.readfile(outfname)
        partition = cpath.partitions[cpath.i_best]
        hist = plotting.get_cluster_size_hist(partition)
        ccfs = cpath.ccfs[cpath.i_best]
        metrics = {'adj_mi' : cpath.adj_mis[cpath.i_best], 'ccf_under' : ccfs[0], 'ccf_over' : ccfs[1]}
        deal_with_parse_results(info, outdir, vname, partition, hist, metrics)

# ----------------------------------------------------------------------------------------
def write_float_val(fname, val, valname):
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, 'w') as outfile:
        writer = csv.DictWriter(outfile, [valname])
        writer.writerow({valname : val})

# ----------------------------------------------------------------------------------------
def read_float_val(fname, valname):
    """ read file with a single number """
    with open(fname) as infile:
        reader = csv.DictReader(infile, fieldnames=[valname])
        for line in reader:
            return float(line[valname])
    raise Exception('something wrong with %s file' % fname)

# ----------------------------------------------------------------------------------------
def make_a_distance_plot(args, metric, combinations, reco_info, cachevals, plotdir, plotname, plottitle):
    def get_joint_key(k1, k2):
        """ figure out which order we have <k1>, <k2> in the cache (if neither, return None) """
        jk1 = k1 + ':' + k2
        jk2 = k2 + ':' + k1
        if jk1 in cachevals:
            return jk1
        elif jk2 in cachevals:
            return jk2
        else:
            return None

    hstyles = ['plain', 'zoom-logy', 'zoom']
    if metric == 'logprob':
        nbins, xmin, xmax = 40, -55, 70
    elif metric == 'naive_hfrac':
        nbins, xmin, xmax = 70, 0., 0.65
    else:
        assert False
    hists = OrderedDict()
    htypes = ['nearest-clones', 'farthest-clones', 'all-clones', 'not']
    for ht in htypes:
        hists[ht] = Hist(nbins, xmin, xmax)
    bigvals, smallvals = {}, {}
    for key_a, key_b in combinations:  # <key_[ab]> is colon-separated string (not a list of keys)
        a_ids, b_ids = key_a.split(':'), key_b.split(':')
        # if not utils.from_same_event(reco_info, a_ids) or not utils.from_same_event(reco_info, b_ids):  # skip clusters that were erroneously merged -- i.e., in effect, assume the previous step didn't screw up at all
        #     raise Exception('woop')
        jk = get_joint_key(key_a, key_b)
        if jk is None:  # if we don't have the joint logprob cached
            continue

        if metric == 'logprob':
            # jk = get_joint_key(key_a, key_b)
            # if jk is None:  # if we don't have the joint logprob cached
            #     continue
            lratio = cachevals[jk] - cachevals[key_a] - cachevals[key_b]
            # print '%f - %f - %f = %f' % (cachevals[jk], cachevals[key_a], cachevals[key_b], lratio),
            mval = lratio
        elif metric == 'naive_hfrac':
            # mval = utils.hamming_fraction(cachevals[key_a], cachevals[key_b])
            mval = cachevals[jk]
        else:
            assert False

        if utils.from_same_event(reco_info, a_ids + b_ids):  # NOTE if a_ids and b_ids are already all from same event, you can optimize this some more
            hists['all-clones'].fill(mval)
            for key in (key_a, key_b):
                if key not in bigvals:
                    bigvals[key] = mval
                if mval > bigvals[key]:
                    bigvals[key] = mval
                if key not in smallvals:
                    smallvals[key] = mval
                if mval < smallvals[key]:
                    smallvals[key] = mval
        else:
            hists['not'].fill(mval)

    if metric == 'logprob':
        bigkey = 'nearest'  # i.e. for logprob ratio, big values mean sequences are nearby
        smallkey = 'farthest'
    elif metric == 'naive_hfrac':
        bigkey = 'farthest'
        smallkey = 'nearest'
    for val in bigvals.values():
        hists[bigkey + '-clones'].fill(val)
    for val in smallvals.values():
        hists[smallkey + '-clones'].fill(val)

    ignore = False
    for hist in hists.values():
        hist.normalize(include_overflows=not ignore, expect_empty=True)

    print ' ', metric, '----------------'
    for hs in hstyles:
        print '   ', hs
        fig, ax = plotting.mpl_init()

        plots = {}
        plots['clonal'] = hists['all-clones'].mpl_plot(ax, ignore_overflows=ignore, label='clonal', alpha=0.7, linewidth=4, color='#6495ed')
        plots['not'] = hists['not'].mpl_plot(ax, ignore_overflows=ignore, label='non-clonal', linewidth=3, color='#2e8b57')  #linewidth=7, alpha=0.5)
        plots['nearest'] = hists['nearest-clones'].mpl_plot(ax, ignore_overflows=ignore, label='nearest clones', linewidth=3)
        plots['farthest'] = hists['farthest-clones'].mpl_plot(ax, ignore_overflows=ignore, label='farthest clones', linewidth=3, linestyle='--')
        if 'log' in hs:
            ax.set_yscale('log')

        leg_loc = (0.1, 0.6)
        xmin = hists['not'].xmin
        xmax = hists['not'].xmax
        ybounds = None
        if 'zoom' in hs:
            if metric == 'logprob':
                nbins, xmin, xmax = 40, 0, 30
                ybounds = (0., 0.1)
            elif metric == 'naive_hfrac':
                nbins, xmin, xmax = 30, 0.02, 0.2
                ybounds = (0., 0.1)
            else:
                assert False
        else:
            if metric == 'logprob':
                pass  #ybounds = (0., 0.)
            elif metric == 'naive_hfrac':
                ybounds = (0., 0.8)

        delta = xmax - xmin
        if metric == 'logprob':
            xlabel = 'log likelihood ratio'
        elif metric == 'naive_hfrac':
            leg_loc = (0.5, 0.6)
            xlabel = 'naive hamming fraction'
        plotting.mpl_finish(ax, plotdir + '/' + hs, plotname, title=plottitle, xlabel=xlabel, ylabel='counts' if 'un-normed' in hs else 'frequency', xbounds=[xmin - 0.03*delta, xmax + 0.03*delta], ybounds=ybounds, leg_loc=leg_loc)
        plotting.make_html(plotdir + '/' + hs)  # this'll overwrite itself a few times

# ----------------------------------------------------------------------------------------
def make_distance_plots(args, baseplotdir, label, n_leaves, mut_mult, cachefname, reco_info, metric):
    cachevals = {}
    everybody, singletons, pairs, triplets, quads = [], [], [], [], []
    with open(cachefname) as cachefile:
        reader = csv.DictReader(cachefile)
        # iline = 0
        for line in reader:
            if metric == 'logprob':
                if line[metric] == '':
                    continue
                cachevals[line['unique_ids']] = float(line['logprob'])
            elif metric == 'naive_hfrac':
                cachevals[line['unique_ids']] = -1. if line['naive_hfrac'] == '' else float(line['naive_hfrac'])  # we need the singletons, even if they don't have hfracs
            else:
                assert False

            unique_ids = line['unique_ids'].split(':')

            # TODO hmm, do I really want to just skip non-clonal ones?
            if not utils.from_same_event(reco_info, unique_ids):
                continue

            everybody.append(line['unique_ids'])
            
            if len(unique_ids) == 1:
                singletons.append(line['unique_ids'])
            elif len(unique_ids) == 2:
                pairs.append(line['unique_ids'])  # use the string so it's obvious which order to use when looking in the cache
            elif len(unique_ids) == 3:
                triplets.append(line['unique_ids'])  # use the string so it's obvious which order to use when looking in the cache
            elif len(unique_ids) == 4:
                quads.append(line['unique_ids'])  # use the string so it's obvious which order to use when looking in the cache
            # iline += 1
            # if iline > 10:
            #     break

    plotdir = baseplotdir + '/distances'
    plotname = leafmutstr(args, n_leaves, mut_mult)

    print 'everybody'
    make_a_distance_plot(args, metric, itertools.combinations(everybody, 2), reco_info, cachevals, plotdir=plotdir + '/' + metric + '/everybody', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (everybody)')

    print 'singletons'
    make_a_distance_plot(args, metric, itertools.combinations(singletons, 2), reco_info, cachevals, plotdir=plotdir + '/' + metric + '/singletons', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (singletons)')

    # print 'one pair one singleton'
    # one_pair_one_singleton = []
    # for ipair in range(len(pairs)):
    #     for ising in range(len(singletons)):
    #         one_pair_one_singleton.append((pairs[ipair], singletons[ising]))
    # make_a_distance_plot(args, metric, one_pair_one_singleton, reco_info, cachevals, plotdir=plotdir + '/' + metric + '/one-pair-one-singleton', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (pair + single)')

    # print 'one triplet one singleton'
    # one_triplet_one_singleton = []
    # for itriplet in range(len(triplets)):
    #     for ising in range(len(singletons)):
    #         one_triplet_one_singleton.append((triplets[itriplet], singletons[ising]))
    # make_a_distance_plot(args, metric, one_triplet_one_singleton, reco_info, cachevals, plotdir=plotdir + '/' + metric + '/one-triplet-one-singleton', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (triple + single)')

    # print 'one quad one singleton'
    # one_quad_one_singleton = []
    # for iquad in range(len(quads)):
    #     for ising in range(len(singletons)):
    #         one_quad_one_singleton.append((quads[iquad], singletons[ising]))
    # make_a_distance_plot(args, metric, one_quad_one_singleton, reco_info, cachevals, plotdir=plotdir + '/' + metric + '/one-quad-one-singleton', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (quad + single)')



    # print 'two pairs'
    # two_pairs = []
    # for ipair in range(len(pairs)):
    #     for jpair in range(ipair + 1, len(pairs)):
    #         two_pairs.append((pairs[ipair], pairs[jpair]))
    # make_a_distance_plot(args, metric, two_pairs, reco_info, cachevals, plotdir=plotdir + '/' + metric + '/two-pairs', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (pair + pair)')

    # print 'two triplets'
    # two_triplets = []
    # for itriplet in range(len(triplets)):
    #     for jtriplet in range(itriplet + 1, len(triplets)):
    #         two_triplets.append((triplets[itriplet], triplets[jtriplet]))
    # make_a_distance_plot(args, metric, two_triplets, reco_info, cachevals, plotdir=plotdir + '/' + metric + '/two-triplets', plotname=plotname, plottitle=get_title(args, label, n_leaves, mut_mult) + ' (triplet + triplet)')

# ----------------------------------------------------------------------------------------
def write_all_plot_csvs(args, label, parameterlist, datafname):
    baseplotdir = os.getenv('www') + '/partis/clustering/' + label
    info = {k : {} for k in metrics + ['hists', 'partitions']}  # NOTE I think I'm not using <info> for anything. But it shouldn't hurt to keep it around
    for params in parameterlist:
        write_each_plot_csvs(args, baseplotdir, label, params['n_leaves'], params['mut_mult'], info, params['hfrac_bounds'], datafname)

    check_call(['./bin/permissify-www', baseplotdir])
    print 'finished!'

# ----------------------------------------------------------------------------------------
def write_each_plot_csvs(args, baseplotdir, label, n_leaves, mut_mult, all_info, hfrac_bounds, datafname):
    for k in all_info:
        if n_leaves not in all_info[k]:
            all_info[k][n_leaves] = {}
        if mut_mult not in all_info[k][n_leaves]:
            all_info[k][n_leaves][mut_mult] = OrderedDict()
        if hfrac_bounds is not None:
            if get_str(hfrac_bounds) not in all_info[k][n_leaves][mut_mult]:
                all_info[k][n_leaves][mut_mult][get_str(hfrac_bounds)] = OrderedDict()

    if hfrac_bounds is None:
        this_info = {k : all_info[k][n_leaves][mut_mult] for k in all_info}
    else:
        this_info = {k : all_info[k][n_leaves][mut_mult][get_str(hfrac_bounds)] for k in all_info}

    plotdir = baseplotdir + '/subsets'
    if args.subset is not None:
        plotdir += '/subset-' + str(args.subset)
    if args.istartstop is not None:
        plotdir += '/istartstop-' + get_str(args.istartstop)

    seqfname = get_seqfile(args, datafname, label, n_leaves, mut_mult)
    if args.data:
        simfbase = None
        plotname = 'data'
    else:
        simfbase = leafmutstr(args, n_leaves, mut_mult, hfrac_bounds)
        plotname = simfbase

    _, reco_info = seqfileopener.get_seqfile_info(seqfname, is_data=args.data)
    if args.count_distances:
        for metric in ['logprob', 'naive_hfrac']:
            make_distance_plots(args, plotdir, label, n_leaves, mut_mult, seqfname.replace('.csv', '-partition-cache.csv'), reco_info, metric)
        print '\nTODO why am I returning here?'
        return

    csvdir = seqfname.replace('.csv', '')
    if hfrac_bounds is not None:
        csvdir += '-hfrac-bounds-' + get_str(hfrac_bounds)

    true_partition = utils.get_true_partition(reco_info)
    if not args.data:
        parse_true(args, this_info, csvdir, true_partition)
        parse_synthetic(args, this_info, csvdir, true_partition, get_outputname(args, label, 'synthetic', seqfname, hfrac_bounds))

    if 'run-viterbi' in args.expected_methods:
        parse_vollmers(args, this_info, get_outputname(args, label, 'run-viterbi', seqfname, hfrac_bounds), csvdir, reco_info, true_partition)
    if 'changeo' in args.expected_methods:
        parse_changeo(args, label, n_leaves, mut_mult, this_info, simfbase, csvdir)
    if 'mixcr' in args.expected_methods:
        parse_mixcr(args, this_info, seqfname, csvdir)

    # partis stuff
    for action in [a for a in args.expected_methods if 'partition' in a]:
        parse_partis(args, action, this_info, get_outputname(args, label, action, seqfname, hfrac_bounds), csvdir, reco_info, true_partition)

    log = 'xy'
    if not args.data and n_leaves <= 10:
        log = 'x'
    title = get_title(args, label, n_leaves, mut_mult, hfrac_bounds)
    plotting.plot_cluster_size_hists(plotdir + '/cluster-size-distributions/' + plotname + '.svg', this_info['hists'], title=title, log=log)  #, xmax=n_leaves*6.01
    plotting.make_html(plotdir + '/cluster-size-distributions')  # this runs a bunch more times than it should
    if not args.no_similarity_matrices:  # they're kinda slow is all
        for meth1, meth2 in itertools.combinations(this_info['partitions'].keys(), 2):
            if '0.5' in meth1 or '0.5' in meth2:  # skip vollmers 0.5
                continue
            n_biggest_clusters = 40  # if args.data else 30)
            plotting.plot_cluster_similarity_matrix(plotdir + '/similarity-matrices/' + (meth1 + '-' + meth2).replace('partition ', ''), plotname, meth1, this_info['partitions'][meth1], meth2, this_info['partitions'][meth2], n_biggest_clusters=n_biggest_clusters, title=get_title(args, label, n_leaves, mut_mult))

# ----------------------------------------------------------------------------------------
def rearrange_metrics_vs_n_leaves(args, valdict, mut_mult_to_use):
    plotvals = OrderedDict()
    for n_leaves in args.n_leaf_list:
        for meth, (val, err) in valdict[n_leaves][mut_mult_to_use].items():
            if meth not in plotvals:
                plotvals[meth] = OrderedDict()
            plotvals[meth][n_leaves] = val, err
    return plotvals

# ----------------------------------------------------------------------------------------
def compare_subsets(args, label):
    baseplotdir = os.getenv('www') + '/partis/clustering/' + label
    if args.hfrac_bound_list is not None and args.istartstop is not None:
        baseplotdir += '/subsets/istartstop-' + get_str(args.istartstop)
    info = {k : {} for k in metrics + ['hists', ]}
    print 'TODO rationalize all these different plotdirs'
    for n_leaves in args.n_leaf_list:
        print '%d leaves' % n_leaves
        for mut_mult in args.mutation_multipliers:
            print '  %.1f mutation' % mut_mult
            compare_subsets_for_each_leafmut(args, baseplotdir, label, n_leaves, mut_mult, info)

    if not args.data:
        if args.plot_mean_of_subsets:  # plot adj mi and stuff for means over subsets, as a function of n_leaves
            for mut_mult in args.mutation_multipliers:
                for metric in metrics:
                    plotvals = rearrange_metrics_vs_n_leaves(args, info[metric], mut_mult)
                    plotting.plot_adj_mi_and_co(plotvals, mut_mult, baseplotdir + '/means-over-subsets/metrics', metric, xvar='n_leaves', title='%dx mutation' % mut_mult)
            plotting.make_html(baseplotdir + '/means-over-subsets/metrics')
        elif args.hfrac_bound_list is not None:
            if len(args.expected_methods) != 1:
                print args.expected_methods
                raise Exception('needs updating if not')
            if args.expected_methods == ['partition']:
                plotting.make_html(baseplotdir + '/plots-vs-thresholds/metrics' + '/logprobs', n_columns=4)
            elif args.expected_methods == ['naive-hamming-partition']:
                plotting.make_html(baseplotdir + '/plots-vs-thresholds/metrics' + '/naive-hfracs', n_columns=4)
            elif args.expected_methods == ['vsearch-partition']:
                plotting.make_html(baseplotdir + '/plots-vs-thresholds/metrics' + '/vsearch-naive-hfracs', n_columns=2)
            else:
                assert False

    check_call(['./bin/permissify-www', baseplotdir])

# ----------------------------------------------------------------------------------------
def compare_subsets_for_each_leafmut(args, baseplotdir, label, n_leaves, mut_mult, all_info):
    """ Where "subsets" means any of a number of different things, or more accurately, this function name means almost nothing. """

    per_subset_info = read_histfiles_and_co(args, label, n_leaves, mut_mult)
    this_info = get_this_info(all_info, n_leaves, mut_mult)

    if args.hfrac_bound_list is not None:  # plot things as a function of hfrac bounds
        plotdir = baseplotdir + '/plots-vs-thresholds/metrics'
        if len(args.expected_methods) != 1:
            print args.expected_methods
            raise Exception('needs updating if not')
        if args.expected_methods == ['partition']:
            plotdir += '/logprobs'
        elif args.expected_methods == ['naive-hamming-partition']:
            plotdir += '/naive-hfracs'
        elif args.expected_methods == ['vsearch-partition']:
            plotdir += '/vsearch-naive-hfracs'
        else:
            assert False
        print 'TODO add mutation level to plot, not just multiplier'
        plotting.plot_metrics_vs_thresholds(args.expected_methods[0], [b[0] for b in args.hfrac_bound_list], per_subset_info, plotdir, plotfname=leafmutstr(args, n_leaves, mut_mult), title=get_title(args, label, n_leaves, mut_mult))
        # for metric in metrics:
        #     for method in get_expected_methods_to_plot(args, metric):
        #         this_info[metric][method] = per_subset_info[metric][method]
    elif args.plot_mean_of_subsets:  # fill this_info with hists of mean over subsets, and plot them
        plot_means_over_subsets(args, label, n_leaves, mut_mult, this_info, per_subset_info, baseplotdir)
    else:  # if we're not averaging over the subsets for each leafmut (and we're not plotting vs thresholds), then we want to plot adj_mi (and whatnot) as a function of subset (presumably each subset is a different size)
        if not args.data:
            for metric in metrics:
                plotvals = OrderedDict()
                for method, values in per_subset_info[metric].items():
                    plotvals[method] = OrderedDict([(nseqs , (val, 0.)) for nseqs, val in zip(nseq_list, values)])
                metric_plotdir = os.getenv('www') + '/partis/clustering/' + label + '/plots-vs-subsets/metrics'
                plotting.plot_adj_mi_and_co(plotvals, mut_mult, metric_plotdir, metric, xvar='nseqs', title=get_title(args, label, n_leaves, mut_mult))
                plotting.make_html(metric_plotdir)

# ----------------------------------------------------------------------------------------
def get_expected_methods_to_plot(args, metric=None):
    # expected_methods = ['vollmers-0.9', 'mixcr', 'changeo', 'vsearch-partition', 'naive-hamming-partition', 'partition']
    expected_methods = list(args.expected_methods)
    print 'TODO make sure expected_methods is ok'
    if len(args.synthetic_partitions) > 0:
        expected_methods += list(args.synthetic_partitions)
    if not args.data and metric == 'hists':
        expected_methods.insert(0, 'true')

    return expected_methods

# ----------------------------------------------------------------------------------------
def read_histfiles_and_co(args, label, n_leaves, mut_mult):
    """ where "subsets" means any of a number of different things """

    basedir = args.fsdir + '/' + label
    if args.n_subsets is not None:
        subdirs = [basedir + '/subset-' + str(isub) for isub in range(args.n_subsets)]
    elif args.istartstoplist is not None:
        subdirs = [basedir + '/istartstop-' + str(istartstop[0]) + '-' + str(istartstop[1]) for istartstop in args.istartstoplist]
        nseq_list = [istartstop[1] - istartstop[0] for istartstop in args.istartstoplist]
    else:
        if args.istartstop is not None:
            basedir += '/istartstop-' + get_str(args.istartstop)
        subdirs = [basedir, ]

    # ----------------------------------------------------------------------------------------
    def read_values(subdir, method, per_subset_info, hfrac_bounds=None):
        if metric == 'hists':
            if args.data:
                raise Exception('I think this needs updating for data')
                # return = subdir + '/data/hists/' + method + '.csv'
            else:
                histfname = subdir + '/' + leafmutstr(args, n_leaves, mut_mult, hfrac_bounds) + '/hists/' + method + '.csv'
            hist = Hist(fname=histfname)
            per_subset_info[metric][method].append(hist)
        else:
            fname = subdir + '/' + leafmutstr(args, n_leaves, mut_mult, hfrac_bounds) + '/' + metric+ '/' + method + '.csv'
            value = read_float_val(fname, metric)
            per_subset_info[metric][method].append(value)

    per_subset_info = {k : OrderedDict() for k in metrics + ['hists', ]}
    for metric in per_subset_info:
        if n_leaves == 1 and metric == 'adj_mi':
            continue
        for method in get_expected_methods_to_plot(args, metric):
            if metric != 'hists' and (args.data or method == 'true'):
                continue
            if method not in per_subset_info[metric]:
                per_subset_info[metric][method] = []

            if len(subdirs) > 1:
                for subdir in subdirs:
                    read_values(subdir, method, per_subset_info)
            else:
                for hfrac_bounds in args.hfrac_bound_list:
                    read_values(subdirs[0], method, per_subset_info, hfrac_bounds)

    return per_subset_info

# ----------------------------------------------------------------------------------------
def get_this_info(all_info, n_leaves, mut_mult):
    for k in all_info:
        if n_leaves not in all_info[k]:
            all_info[k][n_leaves] = {}
        if mut_mult not in all_info[k][n_leaves]:
            all_info[k][n_leaves][mut_mult] = OrderedDict()
    return {k : all_info[k][n_leaves][mut_mult] for k in all_info}

# ----------------------------------------------------------------------------------------
def plot_means_over_subsets(args, label, n_leaves, mut_mult, this_info, per_subset_info, baseplotdir):
    print 'TODO needs testing'
    for method in get_expected_methods_to_plot(args):
        this_info['hists'][method] = plotting.make_mean_hist(per_subset_info['hists'][method])
    cluster_size_plotdir = baseplotdir + '/means-over-subsets/cluster-size-distributions'
    log = 'xy'
    if args.data:
        title = get_title(args, label, n_leaves, mut_mult)
        plotfname = cluster_size_plotdir + '/data.svg'
        xmax = 10
    else:
        title = get_title(args, label, n_leaves, mut_mult)
        plotfname = cluster_size_plotdir + '/' + leafmutstr(args, n_leaves, mut_mult) + '.svg'
        xmax = n_leaves*6.01
        if n_leaves <= 10:
            log = 'x'
    plotting.plot_cluster_size_hists(plotfname, this_info['hists'], title=title, xmax=xmax, log=log)
    plotting.make_html(cluster_size_plotdir)

    if not args.data:
        for metric in metrics:
            print '   ', metric
            for meth, vals in per_subset_info[metric].items():  # add the mean over subsets to <this_info>
                mean = numpy.mean(vals)
                if mean == -1.:  # e.g. mixcr
                    continue
                std = numpy.std(vals)
                this_info[metric][meth] = (mean, std)
                print '        %30s %.3f +/- %.3f' % (meth, mean, std)

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
            true_partition = utils.get_true_partition(reco_info, ids=uids)
            new_partition = utils.generate_incorrect_partition(true_partition, misassign_fraction, error_type=error_type)
            # new_partition = utils.generate_incorrect_partition(true_partition, misassign_fraction, error_type='singletons')
            new_partitions[nseqs] = new_partition
    return {nseqs : utils.adjusted_mutual_information(new_partitions[nseqs], utils.get_true_partition(reco_info, ids=new_partitions[nseqs].keys())) for nseqs in nseq_list}

# ----------------------------------------------------------------------------------------
def output_exists(args, outfname):
    if os.path.exists(outfname):
        if os.stat(outfname).st_size == 0:
            print '                      deleting zero length %s' % outfname
            os.remove(outfname)
            return False
        elif args.overwrite:
            print '                      overwriting %s' % outfname
            if os.path.isdir(outfname):
                raise Exception('output %s is a directory, rm it by hand' % outfname)
            else:
                os.remove(outfname)
            return False
        else:
            print '                      output exists, skipping (%s)' % outfname
            return True
    else:
        return False

# ----------------------------------------------------------------------------------------
def run_changeo(args, label, n_leaves, mut_mult, seqfname):
    def untar_imgt(imgtdir):
        tar_cmd = 'mkdir ' + imgtdir + ';'
        tar_cmd += ' tar Jxvf ' + imgtdir + '.txz --exclude=\'IMGT_HighV-QUEST_individual_files_folder/*\' -C ' + imgtdir
        check_call(tar_cmd, shell=True)

    imgtdir = get_changeo_outdir(args, label, n_leaves, mut_mult)
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
        check_call(cmdstr.split(), env=os.environ)

    resultfname = imgtdir + changeorandomcrapstr
    if output_exists(args, resultfname):
        return

    fastafname = os.path.splitext(seqfname)[0] + '.fasta'
    utils.csv_to_fasta(seqfname, outfname=fastafname)  #, name_column='name' if args.data else 'unique_id', seq_column='nucleotide' if args.data else 'seq')
    bindir = '/home/dralph/work/changeo/changeo/bin'
    os.environ['PYTHONPATH'] = bindir.replace('/bin', '')
    start = time.time()
    cmd = bindir + '/MakeDb.py imgt -i ' + imgtdir + ' -s ' + fastafname + ' --failed'
    # cmd = bindir + '/MakeDb.py imgt -h'
    run(cmd)
    cmd = bindir + '/ParseDb.py select -d ' + imgtdir + '_db-pass.tab'
    if args.data:
        cmd += ' -f FUNCTIONAL -u T'
    else:  # on simulation we don't want to skip any (I'm not forbidding stop codons in simulation)
        cmd += ' -f FUNCTIONAL -u T F'
    run(cmd)
    # cmd = bindir + '/DefineClones.py bygroup -d ' + imgtdir + '_db-pass_parse-select.tab --act first --model m1n --dist 7'
    cmd = bindir + '/DefineClones.py bygroup -d ' + imgtdir + '_db-pass_parse-select.tab --model hs1f --norm len --act set --dist 0.2'
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
        true_partition = utils.get_true_partition(reco_info)
        subset_of_true_partition = utils.remove_missing_uids_from_true_partition(true_partition, partition)
        print 'removed from true: %.3f' % utils.adjusted_mutual_information(subset_of_true_partition, partition)

        partition_with_uids_added = utils.add_missing_uids_as_singletons_to_inferred_partition(utils.get_true_partition(reco_info), partition)
        print 'added to inferred: %.3f' % utils.adjusted_mutual_information(true_partition, partition_with_uids_added)
        assert False  # need to work out why changeo is filtering out so many seqs, and decide how to treat them if I can't fix it

        write_float_val(imgtdir + '-adj_mi.csv', adj_mi, 'adj_mi')
        raise Exception('dont do this here, do it somewhere else')
        ccfs = utils.new_ccfs_that_need_better_names(partition, true_partition, reco_info)
        write_float_val(imgtdir + '-ccf_under.csv', ccfs[0], 'ccf_under')
        write_float_val(imgtdir + '-ccf_over.csv', ccfs[1], 'ccf_over')

# ----------------------------------------------------------------------------------------
def run_mixcr(args, label, n_leaves, mut_mult, seqfname):
    binary = '/home/dralph/work/mixcr/mixcr-1.2/mixcr'
    mixcr_workdir = get_program_workdir(args, 'mixcr', label, n_leaves, mut_mult)
    if not os.path.exists(mixcr_workdir):
        os.makedirs(mixcr_workdir)

    # fastafname = os.path.splitext(seqfname)[0] + '.fasta'
    infname = mixcr_workdir + '/' + os.path.basename(os.path.splitext(seqfname)[0] + '.fasta')
    outfname = os.path.splitext(seqfname)[0] + '-mixcr.tsv'
    if output_exists(args, outfname):
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


# ----------------------------------------------------------------------------------------
def run_igscueal(args, label, n_leaves, mut_mult, seqfname):
    igscueal_dir = '/home/dralph/work/IgSCUEAL'
    # outfname = os.path.splitext(seqfname)[0] + '-igscueal.tsv'
    # if output_exists(args, outfname):
    #     return
    workdir = get_program_workdir(args, 'igscueal', label, n_leaves, mut_mult)

    infname = workdir + '/' + os.path.basename(os.path.splitext(seqfname)[0] + '.fasta')

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    utils.csv_to_fasta(seqfname, outfname=infname, n_max_lines=args.n_to_partition)  #, name_column='name' if args.data else 'unique_id', seq_column='nucleotide' if args.data else 'seq'
    # write cfg file (.bf)
    sed_cmd = 'sed'
    replacements = {'igscueal_dir' : igscueal_dir,
                    'input_fname' : infname,
                    'results_fname' : workdir + '/results.tsv',
                    'rearrangement_fname' : workdir + '/rearrangement.tsv',
                    'tree_assignment_fname' : workdir + '/tree_assignment.tsv'}
    for pattern, replacement in replacements.items():
        sed_cmd += ' -e \'s@xxx-' + pattern + '-xxx@' + replacement + '@\''
    template_cfgfname = igscueal_dir + '/TopLevel/MPIScreenFASTA.bf'
    cfgfname = workdir + '/cfg.bf'
    sed_cmd += ' ' + template_cfgfname + ' >' + cfgfname
    check_call(sed_cmd, shell=True)

    # cmd = 'salloc -N 3 mpirun -np 3 /home/dralph/work/hyphy/hyphy-master/HYPHYMPI ' + cfgfname
    # srun --mpi=openmpi
    cmd = 'srun --exclude=data/gizmod.txt mpirun -np 2 /home/dralph/work/hyphy/hyphy-master/HYPHYMPI ' + cfgfname

    ntot = int(check_output(['wc', '-l', infname]).split()[0]) / 2
    n_procs = max(1, int(float(ntot) / 10))
    n_per_proc = int(float(ntot) / n_procs)  # NOTE ignores remainders, i.e. last few sequences
    workdirs = []
    start = time.time()
    procs = []
    for iproc in range(n_procs):
        workdirs.append(workdir + '/igs-' + str(iproc))
        if not os.path.exists(workdirs[-1]):
            os.makedirs(workdirs[-1])

        if len(procs) - procs.count(None) > 500:  # can't have more open files than something like this
            print '        too many procs (len %d    none %d)' % (len(procs), procs.count(None))
            procs.append(None)
            continue

        suboutfname = replacements['results_fname'].replace(workdir, workdirs[-1])
        if os.path.exists(suboutfname) and os.stat(suboutfname).st_size != 0:
            print '    %d already there (%s)' % (iproc, suboutfname)
            procs.append(None)
            continue

        check_call(['cp', cfgfname, workdirs[-1] + '/'])
        check_call(['sed', '-i', 's@' + workdir + '@' + workdirs[-1] + '@', workdirs[-1] + '/' + os.path.basename(cfgfname)])

        subinfname = workdirs[-1] + '/' + os.path.basename(infname)
        istart = 2 * iproc * n_per_proc + 1  # NOTE sed indexing (one-indexed with inclusive bounds), and factor of two for fasta file
        istop = istart + 2 * n_per_proc - 1
        check_call('sed -n \'' + str(istart) + ',' + str(istop) + ' p\' ' + infname + '>' + subinfname, shell=True)

        print '     starting %d' % iproc
        procs.append(Popen(cmd.replace(workdir, workdirs[-1]).split(), stdout=PIPE, stderr=PIPE))
        # procs.append(Popen(['sleep', '10']))

    while procs.count(None) < len(procs):
        for iproc in range(n_procs):
            if procs[iproc] is not None and procs[iproc].poll() is not None:  # it's finished
                stdout, stderr = procs[iproc].communicate()
                print '\nproc %d' % iproc
                print 'out----\n', stdout, '\n-----'
                print 'err----\n', stderr, '\n-----'
                procs[iproc] = None
            time.sleep(0.1)
    print '      igscueal time: %.3f' % (time.time()-start)

# ----------------------------------------------------------------------------------------
def slice_file(args, csv_infname, csv_outfname):  # not necessarily csv
    if os.path.exists(csv_outfname):
        if not os.path.exists(csv_outfname.replace('.csv', '.fa')):
            utils.csv_to_fasta(csv_outfname)
        return
    print '      subsetting %d seqs with indices %d --> %d' % (args.istartstop[1] - args.istartstop[0], args.istartstop[0], args.istartstop[1])
    if not os.path.exists(os.path.dirname(csv_outfname)):
        os.makedirs(os.path.dirname(csv_outfname))
    if '.csv' in csv_infname:  # if it's actually a csv
        remove_csv_infname = False
        check_call('head -n1 ' + csv_infname + ' >' + csv_outfname, shell=True)
        check_call('sed -n \'' + str(args.istartstop[0] + 2) + ',' + str(args.istartstop[1] + 1) + ' p\' ' + csv_infname + '>>' + csv_outfname, shell=True)  # NOTE conversion from standard zero indexing to sed inclusive one-indexing (and +1 for header line)
        if remove_csv_infname:
            assert '/dralph/' in csv_infname
            os.remove(csv_infname)
    elif '.fa' in csv_infname:
        input_info, _ = seqfileopener.get_seqfile_info(csv_infname, is_data=True)
        with open(csv_outfname, 'w') as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'seq'))
            writer.writeheader()
            iseq = -1
            for line in input_info.values():  # hackey, but it's an ordered dict so it should be ok
                iseq += 1
                if iseq < args.istartstop[0]:
                    continue
                if iseq >= args.istartstop[1]:
                    break
                writer.writerow({'unique_id' : line['unique_id'], 'seq' : line['seq']})
        # print 'sed -n \'' + str(2*args.istartstop[0] + 1) + ',' + str(2*args.istartstop[1] + 1) + ' p\' ' + csv_infname + '>>' + csv_outfname
        # check_call('sed -n \'' + str(2*args.istartstop[0] + 1) + ',' + str(2*args.istartstop[1] + 1) + ' p\' ' + csv_infname + '>>' + csv_outfname, shell=True)  # NOTE conversion from standard zero indexing to sed inclusive one-indexing (and multiply by two for fasta file)

# ----------------------------------------------------------------------------------------
def get_seqfile(args, datafname, label, n_leaves, mut_mult):

    if args.data:
        if args.istartstop is None:
            seqfname = datafname
        else:
            subfname = args.fsdir + '/' + label + '/istartstop-' + get_str(args.istartstop) + '/data.csv'
            slice_file(args, datafname, subfname)
            seqfname = subfname
    else:
        if not args.data:
            assert n_leaves is not None and mut_mult is not None
        simfname = get_simfname(args, label, n_leaves, mut_mult, no_subset=True)

        if args.subset is not None:
            ntot = int(check_output(['wc', '-l', simfname]).split()[0]) - 1
            n_per_subset = int(float(ntot) / args.n_subsets)  # NOTE ignores remainders, i.e. last few sequences
            args.n_to_partition = n_per_subset  # we use this to set the number of procs (and maybe for other things as well)
            subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
            if os.path.exists(subsimfname):
                pass
                # print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d / %d' % (args.subset, args.n_subsets)
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                istart = args.subset * n_per_subset + 2  # NOTE sed indexing (one-indexed with inclusive bounds). Also note extra +1 to avoid header
                istop = istart + n_per_subset - 1
                check_call('sed -n \'' + str(istart) + ',' + str(istop) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)
            simfname = subsimfname

        if args.istartstop is not None:
            subsimfname = simfname.replace(label + '/', label + '/istartstop-' + get_str(args.istartstop) + '/')
            slice_file(args, simfname, subsimfname)
            simfname = subsimfname

        seqfname = simfname

    return seqfname

# ----------------------------------------------------------------------------------------
def get_outputname(args, label, action, seqfname, hfrac_bounds):
    if args.data:
        outputname = get_outdirname(args, label) + '/data-' + action + '.csv'
    else:
        outputname = ('-' + action).join(os.path.splitext(seqfname))
    if hfrac_bounds is not None:
        outputname = outputname.replace('.csv', '-hfrac-bounds-' + get_str(hfrac_bounds) + '.csv')
    return outputname

# ----------------------------------------------------------------------------------------
def execute(args, action, datafname, label, n_leaves, mut_mult, procs, hfrac_bounds=None):
    cmd = './bin/run-driver.py --label ' + label + ' --action '
    if 'partition' in action:
        cmd += ' partition'
    else:
        cmd += ' ' + action
    cmd += ' --stashdir ' + args.fsdir + ' --old-style-dir-structure'

    extras = []
    seqfname = get_seqfile(args, datafname, label, n_leaves, mut_mult)
    if args.data:
        cmd += ' --datafname ' + seqfname + ' --is-data'
        # if args.dataset == 'adaptive':
        extras += ['--skip-unproductive', ]
    else:
        cmd += ' --simfname ' + seqfname

    n_procs, n_fewer_procs = 1, 1
    if action == 'cache-data-parameters':
        outfname = args.fsdir + '/' + label + '/data'
        extras += ['--n-max-queries', + args.n_data_to_cache]
        n_procs = max(1, args.n_data_to_cache / 500)
        n_fewer_procs = min(500, args.n_data_to_cache / 2000)
    elif action == 'simulate':
        outfname = seqfname
        extras += ['--n-sim-events', int(float(args.n_sim_seqs) / n_leaves)]
        extras += ['--n-leaves', n_leaves, '--mutation-multiplier', mut_mult]
        if args.indels:
            extras += ['--indel-frequency', 0.5]
        if args.lonely_leaves:
            extras += ['--constant-number-of-leaves', ]
        if args.mimic:
            extras += ['--mimic-data-read-length', ]
        n_procs = 10
    elif action == 'cache-simu-parameters':
        outfname = seqfname.replace('.csv', '')
        n_procs = 20
        n_fewer_procs = min(500, args.n_sim_seqs / 2000)
    elif action == 'partition':
        outfname = get_outputname(args, label, action, seqfname, hfrac_bounds)
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition]
        if args.count_distances:
            extras += ['--cache-naive-hfracs', '--persistent-cachefname', ('-cache').join(os.path.splitext(outfname))]  # '--n-partition-steps', 1,

        if hfrac_bounds is not None:
            print 'TODO change name from hfrac_bounds'
            assert hfrac_bounds[0] == hfrac_bounds[1]
            extras += ['--logprob-ratio-threshold', hfrac_bounds[0]]

        n_procs = max(1, args.n_to_partition / 100)
        n_fewer_procs = min(500, args.n_to_partition / 2000)
    elif action == 'naive-hamming-partition':
        outfname = get_outputname(args, label, action, seqfname, hfrac_bounds)
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition, '--naive-hamming']

        if hfrac_bounds is not None:
            extras += ['--naive-hamming-bounds', get_str(hfrac_bounds, delimiter=':')]

        n_procs = max(1, args.n_to_partition / 200)
    elif action == 'vsearch-partition':
        outfname = get_outputname(args, label, action, seqfname, hfrac_bounds)
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', args.n_to_partition, '--naive-vsearch']
        # # ----------------------------------------------------------------------------------------
        # print 'TODO remove this'
        # extras += ['--persistent-cachefname', seqfname.replace('.csv', '-naive-seq-cache.csv')]
        # # ----------------------------------------------------------------------------------------
        if hfrac_bounds is not None:
            extras += ['--naive-hamming-bounds', get_str(hfrac_bounds, delimiter=':')]
        n_procs = max(1, args.n_to_partition / 100)  # only used for ighutil step
        # n_procs = 10
        # n_fewer_procs = 5
    elif action == 'run-viterbi':
        outfname = get_outputname(args, label, action, seqfname, hfrac_bounds)
        cmd += ' --outfname ' + outfname
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        extras += ['--n-max-queries', args.n_to_partition]
        n_procs = max(1, args.n_to_partition / 50)
    elif action == 'run-changeo':
        run_changeo(args, label, n_leaves, mut_mult, seqfname)
        return
    elif action == 'run-mixcr':
        run_mixcr(args, label, n_leaves, mut_mult, seqfname)
        return
    elif action == 'run-igscueal':
        run_igscueal(args, label, n_leaves, mut_mult, seqfname)
        return
    elif action == 'synthetic':
        outfname = get_outputname(args, label, action, seqfname, hfrac_bounds)
        generate_synthetic_partitions(args, label, n_leaves, mut_mult, seqfname, outfname)
        return
    else:
        raise Exception('bad action %s' % action)

    if output_exists(args, outfname):
        return

    # cmd += ' --plotdir ' + os.getenv('www') + '/partis'
    if n_procs > 500:
        print 'reducing n_procs %d --> %d' % (n_procs, 500)
        n_procs = 500
    n_proc_str = str(n_procs)
    extras += ['--workdir', args.fsdir.replace('_output', '_tmp') + '/' + str(random.randint(0, 99999))]
    if action != 'simulate':
        extras += ['--slurm', ]
    if n_procs > 10:
        n_fewer_procs = max(1, n_fewer_procs)
        n_proc_str += ':' + str(n_fewer_procs)

    extras += ['--n-procs', n_proc_str]

    cmd += baseutils.get_extra_str(extras)
    print '   ' + cmd
    # return

    logbase = os.path.dirname(outfname) + '/_logs/' + os.path.basename(outfname).replace('.csv', '')
    if action not in logbase:
        logbase += '-' + action

    if not os.path.exists(os.path.dirname(logbase)):
        os.makedirs(os.path.dirname(logbase))
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    time.sleep(30)  # 300sec = 5min
