import os
import sys
import utils
import plotting
import lbplotting
import treeutils
import collections
import json
import numpy
import math
from hist import Hist
import copy

# ----------------------------------------------------------------------------------------
vlabels = {
    'obs_frac' : 'fraction sampled',
    'n-sim-seqs-per-gen' : 'N/gen',
    'obs-times' : 't obs',
    'carry-cap' : 'carry cap',
}
linestyles = {'lbi' : '-', 'lbr' : '-', 'dtr' : '--'}
linestyles.update(plotting.linestyles)
linewidths = {'lbi' : 2.5, 'lbr' : 2.5, 'dtr' : 3}
hard_colors = {'dtr' : '#626262',
               'aa-lbi' : '#e043b9',
               'aa-lbr' : '#e043b9'}  # don't like the cycle colors these end up with
hard_colors.update(plotting.colors)
def metric_color(metric):  # as a fcn to avoid import if we're not plotting
    if metric in hard_colors:
        return hard_colors[metric]
    if 'synth-' in metric:
        return 'grey'
    mstrlist = ['shm:lbi:cons-dist-aa:cons-dist-nuc:dtr:aa-lbi', 'delta-lbi:lbr:dtr:aa-lbr']
    metric_colors = {m : plotting.frozen_pltcolors[i % len(plotting.frozen_pltcolors)] for mstrs in mstrlist for i, m in enumerate(mstrs.split(':'))}
    return metric_colors.get(metric, 'red')
ltexts = {
    'single' : 'single chain',
    'joint' : 'joint',
}

# ----------------------------------------------------------------------------------------
# NOTE it's kind of arbitrary what gets a subdir (e.g. per_x and perf_metric) vs what doesn't (e.g. locus, choice_grouping)
def get_comparison_plotdir(args, mtmp, per_x=None, extra_str='', perf_metric=None): # NOTE that in make_plots() per_x is None indicates we're doing paired clustering, but here it does not (can also mean we want the parent dir)
    plotdir = '%s/%s%s/plots' % (args.base_outdir, args.label, '/'+args.version if hasattr(args, 'version') else '')
    if mtmp is not None:  # for tree metrics, this is the metric, whereas for paired this is the method and perf_metric is the metric
        plotdir += '/' + mtmp
        if mtmp == 'combined' and args.combo_extra_str is not None:
            plotdir += '-' + args.combo_extra_str
    if extra_str != '':
        assert mtmp is not None
        plotdir += '_' + extra_str
    if per_x is not None:
        plotdir += '/' + per_x
    if perf_metric is not None:
        plotdir += '/' + perf_metric
    return plotdir

# ----------------------------------------------------------------------------------------
def cp_val(cpath, ptilestr, yfname):
    if ptilestr == 'precision':
        rval = cpath.ccfs[cpath.i_best][0]
    elif ptilestr == 'sensitivity':
        rval = cpath.ccfs[cpath.i_best][1]
    elif ptilestr == 'f1':
        rval = sys.modules['scipy.stats'].hmean(cpath.ccfs[cpath.i_best])
    elif ptilestr == 'cln-frac':
        rval = utils.collision_fraction(cpath.best())
    else:
        assert False
    if rval is None:
        print '  %s read none type val from %s' % (utils.color('yellow', 'warning'), yfname)
    return rval

# ----------------------------------------------------------------------------------------
def read_hist_csv(fname, ptilestr):
    hist = Hist(fname=fname)
    if 'pcfrac-' in ptilestr:
        blabels = {'pcfrac-corr' : 'correct', 'pcfrac-mis' : 'mispaired', 'pcfrac-un' : 'unpaired'}
        pval = hist.bin_contents[hist.find_bin(None, label=blabels[ptilestr])]
    else:
        pval = hist.get_mean()
    return {ptilestr : pval}

# ----------------------------------------------------------------------------------------
def readlog(fname, metric, locus, ptntype):
    # ----------------------------------------------------------------------------------------
    def timestr():
        if ptntype == 'joint':
            return 'merge time'
        elif metric == 'partition':
            return 'loop time:'
        elif metric == 'vsearch-partition':
            return 'vsearch time:'
        else:
            return '%s time:' % metric
    # ----------------------------------------------------------------------------------------
    with open(fname) as lfile:
        flines = lfile.readlines()
    if 'partition' in metric:  # partis methods
        lstrs = ['./bin/partis', timestr()] if ptntype=='single' else [timestr()]
        tlines = []
        for tln in flines:
            for lstr in lstrs:
                if lstr in tln:
                    tlines.append((lstr, tln))
        if ptntype == 'single':
            assert len(tlines) == 5
            tlines = tlines[1:]  # remove the paired cmd
            tvals = {}
            for ltmp, itmp in zip(['igh', 'igk'], [0, 2]):
                lstr, cmdline = tlines[itmp]
                assert ltmp == utils.get_val_from_arglist(cmdline.split(), '--locus')
                _, timeline = tlines[itmp + 1]
                _, _, timestr = timeline.split()
                tvals[ltmp] = float(timestr)
            return {'time-reqd' : tvals[locus]}
        else:
            _, timeline = utils.get_single_entry(tlines)
            _, _, timestr = timeline.split()
            return {'time-reqd' : float(timestr)}
    else:  # other methods
        assert ptntype == 'single'  # at least for now
        tlines = [l for l in flines if timestr() in l]
        tloci = ['igh', 'igk'] if metric=='scoper' else [locus]
        if len(tlines) != len(tloci):
            raise Exception('couldn\'t find exactly %d time lines for loci %s (got %d) in %s' % (len(tloci), tloci, len(tlines), fname))
        tvals = {}
        for ltmp, tline in zip(tloci, tlines):
            if metric == 'scoper':
                tltmp, _, _, timestr = tline.split()
                assert tltmp == ltmp
            else:
                _, _, timestr = tline.split()
            tvals[ltmp] = float(timestr)
        return {'time-reqd' : tvals[locus]}

# ----------------------------------------------------------------------------------------
# <metric>: for tree metrics this is the metric (e.g. lbi), for paired clustering this is the method (e.g. partis) and <ptilestr> is the metric
# <ptilestr>: x var in ptile plots (y var in final plots), i.e. whatever's analagous to [var derived from] 'affinity' or 'n-ancestor' (<ptilelabel> is its label), i.e. for paired this is f1, precision, sensitivity
# <xvar>: x var in *final* plot
def make_plots(args, svars, action, metric, ptilestr, ptilelabel, xvar, fnfcn=None, per_x=None, choice_grouping=None, min_ptile_to_plot=75., use_relative_affy=False, metric_extra_str='',
               locus=None, ptntype=None, xdelim='_XTRA_', pdirfcn=None, fnames=None, make_legend=False, debug=False):  # NOTE I started trying to split fcns out of here, but you have to pass around too many variables it's just not worth it
    # ----------------------------------------------------------------------------------------
    def legstr(label, title=False):
        if label is None: return None
        jstr = '\n' if title else '; '
        tmplist = [legdict.get(l, l.replace('-', ' ')) for l in label.split('; ')]
        if title and args.pvks_to_plot is not None:  # if we're only plotting specific values, put them in the legend str (typically we're just plotting one value)
            assert isinstance(args.pvks_to_plot, list)  # don't really need this
            for il in range(len(tmplist)):
                subpvks = [pvk.split('; ')[il] for pvk in args.pvks_to_plot]
                tmplist[il] += ': %s' % ' '.join(legdict.get(spvk, spvk) for spvk in subpvks)
        lstr = jstr.join(tmplist)
        # if per_x is None and ptntype is not None and label in args.plot_metrics:  # need to add single/joint to the method
        #     lstr += ' %s' % ltexts[ptntype]
        if title:
            lstr = lstr.replace(' fraction (nuc)', '')
        return lstr
    # ----------------------------------------------------------------------------------------
    def nsimevts():
        if per_x is None:  # for paired, we should always be accessing n sim events as a scan var
            assert False
        else:  # see hasattr assertions below
            return args.n_sim_events_per_proc
    # ----------------------------------------------------------------------------------------
    def get_obs_frac(vlists, varnames, return_dbg=False):
        if per_x is None:  # really just so when you're looking at this fcn you don't need to worry what happens for paired clustering
            assert False
        obs_times = utils.vlval(args, vlists, varnames, 'obs-times')
        n_per_gen_vals = utils.vlval(args, vlists, varnames, 'n-sim-seqs-per-gen')
        if len(obs_times) == len(n_per_gen_vals):  # note that this duplicates logic in bcr-phylo simulator.py
            n_sampled = sum(n_per_gen_vals)
        elif len(n_per_gen_vals) == 1:
            n_sampled = len(obs_times) * n_per_gen_vals[0]
        else:
            assert False
        n_total = utils.vlval(args, vlists, varnames, 'carry-cap')  # note that this is of course the number alive at a given time, and very different from the total number that ever lived
        obs_frac = n_sampled / float(n_total)
        if return_dbg:
            dbgstr = '    %-12s %-12s   %-5d     %8s / %-4d = %.3f' % (' '.join(str(o) for o in obs_times), ' '.join(str(n) for n in n_per_gen_vals), n_total,
                                                                       ('(%s)' % ' + '.join(str(n) for n in n_per_gen_vals)) if len(obs_times) == len(n_per_gen_vals) else ('%d * %d' % (len(obs_times), n_per_gen_vals[0])),
                                                                       n_total, obs_frac)
            return dbgstr
        else:
            return obs_frac
    # ----------------------------------------------------------------------------------------
    def get_n_seqs(vlists, varnames):
        n_leaves, n_events = [utils.vlval(args, vlists, varnames, vstr) for vstr in ['n-leaves', 'n-sim-events']]
        return int(n_leaves) * int(n_events)
    # ----------------------------------------------------------------------------------------
    def pvkeystr(vlists, varnames):
        # ----------------------------------------------------------------------------------------
        def valstr(vname):
            if vname == 'obs_frac':  # obs_frac has to be treated diffferently since it's a derived value from several parameters, whereas others are just one parameter (i.e. it only matters if we're going to set obs_frac as final_plot_xvar)
                return '%.4f' % get_obs_frac(vlists, varnames)  # shouldn't happen for paired
            vval = utils.vlval(args, vlists, varnames, vname)
            def strfcn(x):
                return str(x)  # TODO
            if isinstance(vval, list):
                return ', '.join(strfcn(v) for v in vval)
            else:
                return strfcn(vval)
        # ----------------------------------------------------------------------------------------
        pvnames = sorted(set(varnames) - set(['seed', xvar]))
        if args.legend_var is not None:  # pvnames == ['n-sim-seqs-per-gen']:  # if this is the only thing that's different between different runs (except for the x variable and seed/replicate) then we want to use obs_frac
            pvnames = [args.legend_var]  # ['obs_frac']
        pvkey = '; '.join(valstr(vn) for vn in pvnames)  # key identifying each line of a different color
        pvlabel[0] = '; '.join(vlabels.get(vn, vn) for vn in pvnames)
        return pvkey
    # ----------------------------------------------------------------------------------------
    def get_ytmpfo(yamlfo, iclust=None):
        if 'percentiles' in yamlfo:  # new-style files
            ytmpfo = yamlfo['distr-hists' if distr_hists else 'percentiles']
            if per_x == 'per-seq':
                ytmpfo = ytmpfo.get('per-seq', ytmpfo)['all-clusters' if iclust is None else 'iclust-%d' % iclust]  # distr-hists don't have 'per-seq'/'per-cluster' level
            else:
                ytmpfo = ytmpfo.get('per-cluster', ytmpfo)[choice_grouping]
        else:  # old-style files
            ytmpfo = yamlfo
            if iclust is not None:
                if 'iclust-%d' % iclust not in ytmpfo:
                    print '    %s requested per-cluster ptile vals, but they\'re not in the yaml file (probably just an old file)' % utils.color('yellow', 'warning')  # I think it's just going to crash on the next line anyway
                ytmpfo = ytmpfo['iclust-%d' % iclust]
        return ytmpfo
    # ----------------------------------------------------------------------------------------
    def yval_key(ytmpfo):
        if ptilestr == 'affinity' and 'mean_affy_ptiles' in ytmpfo:  # old-style files used shortened version
            return 'mean_affy_ptiles'
        else:
            return 'mean_%s_ptiles' % ptilestr
    # ----------------------------------------------------------------------------------------
    def get_ptile_diff_vals(ytmpfo, iclust=None):  # the perfect line is higher for lbi, but lower for lbr, hence the abs(). Occasional values can go past/better than perfect, so maybe it would make sense to reverse sign for lbi/lbr rather than taking abs(), but I think this is better
        ytmpfo = get_ytmpfo(ytmpfo, iclust=iclust)
        return [abs(pafp - afp) for lbp, afp, pafp in zip(ytmpfo['lb_ptiles'], ytmpfo[yval_key(ytmpfo)], ytmpfo['perfect_vals']) if lbp > min_ptile_to_plot]
    # ----------------------------------------------------------------------------------------
    def get_dhist_vals(ytmpfo, iclust=None):  # specify args.distr_hist_limit in one of three forms: (N, 3) (take 3 from each family), (frac, 0.01) (take 1% of family/total), or (val, 7) (take those with lb value greater than 7)
        # ----------------------------------------------------------------------------------------
        def gethist(tkey):
            yfo = ytmpfo[k]
            htmp = Hist(n_bins=yfo['n_bins'], xmin=yfo['xmin'], xmax=yfo['xmax'])
            for ibin in htmp.ibiniter(True):
                htmp.bin_contents[ibin] = yfo['bin_contents'][ibin]
            return htmp
        # ----------------------------------------------------------------------------------------
        ytmpfo = get_ytmpfo(ytmpfo, iclust=iclust)
        hzero, hother = [gethist(k) for k in ['zero', 'not 0']]
        for tattr in ['n_bins', 'xmin', 'xmax']:
            assert getattr(hzero, tattr) == getattr(hother, tattr)
        if args.distr_hist_limit[0] == 'val':
            ibin = hzero.find_bin(args.distr_hist_limit[1])
        elif args.distr_hist_limit[0] in ['frac', 'N']:
            hsum = copy.deepcopy(hzero)
            hsum.add(hother)
            tsum, total = 0., hsum.integral(True)
            tfrac = args.distr_hist_limit[1]
            if args.distr_hist_limit[0] == 'N':
                tfrac = args.distr_hist_limit[1] / total
                if iclust is None:
                    tfrac *= nsimevts()
            for ibin in hsum.ibiniter(True, reverse=True):
                tsum += hsum.bin_contents[ibin]
                if tsum >= tfrac * total:
                    break
        else:
            assert False
        n_zero, n_other = [h.integral(False, ibounds=(ibin, h.n_bins+2)) for h in [hzero, hother]]
        return n_zero, n_other, hzero.low_edges[ibin]
    # ----------------------------------------------------------------------------------------
    def get_varname_str():
        return ''.join(utils.wfmt(vlabels.get(v, v), 25 if per_x is None else 10) for v in varnames)
    def get_varval_str(vstrs):
        return ''.join(utils.wfmt(v, 24 if per_x is None else 9) for v in vstrs)
    # ----------------------------------------------------------------------------------------
    def read_plot_info():
        # ----------------------------------------------------------------------------------------
        def add_plot_vals(ytmpfo, vlists, varnames, iclust=None):
            # ----------------------------------------------------------------------------------------
            def getikey():
                if args.n_replicates == 1:
                    raise Exception('this should be checked, especially as regards nsimevts()')
                if args.n_replicates == 1 and treat_clusters_together:  # could add 'or per_x is None' next to treat_clusters_together in all these places, but whatever
                    ikey = None
                    def initfcn(): return []  # i swear it initially made more sense for this to be such a special case
                elif args.n_replicates == 1:  # but more than one event per proc
                    ikey = iclust
                    def initfcn(): return {i : [] for i in range(nsimevts())}
                elif treat_clusters_together:  # but more than one replicate/seed (paired should always use this one)
                    ikey = vlists[varnames.index('seed')]
                    def initfcn(): return {i : [] for i in utils.sargval(args, 'seed')}
                else:  # both of 'em non-trivial
                    ikey = '%d-%d' % (vlists[varnames.index('seed')], iclust)
                    def initfcn(): return {('%d-%d' % (i, j)) : [] for i in utils.sargval(args, 'seed') for j in range(nsimevts())}
                return ikey, initfcn
            # ----------------------------------------------------------------------------------------
            if per_x is None:
                fval = ytmpfo[ptilestr]
                if debug:
                    print ' %.2f' % fval,
            elif distr_hists:
                n_zero, n_other, lo_edge = get_dhist_vals(ytmpfo, iclust=iclust)
                if n_zero + n_other == 0:
                    # print '  %s zero total in distr hists' % utils.color('yellow', 'warning')
                    return
                fval =  n_zero / float(n_zero + n_other)
                if debug:
                    nd = str(0 if lo_edge >= 10. else 1)
                    print '%s%d/%d=%.2f (%s)' % ('    ' if iclust in [0, None] else ' ', n_zero, n_zero + n_other, fval, ('%3.'+nd+'f')%lo_edge),
            else:
                diff_vals = get_ptile_diff_vals(ytmpfo, iclust=iclust)
                if len(diff_vals) == 0:
                    missing_vstrs['empty'].append((iclust, vstrs))  # empty may be from empty list in yaml file, or may be from none of them being above <min_ptile_to_plot>
                    return
                fval = numpy.mean(diff_vals)  # diff_to_perfect
                if debug:
                    print ' %.2f' % fval,
            tau = utils.vlval(args, vlists, varnames, xvar)  # not necessarily tau anymore (i think it's just final_plot_xvar?)
            ikey, initfcn = getikey()
            pvkey = pvkeystr(vlists, varnames)  # key identifying each line in the plot, each with a different color, (it's kind of ugly to get the label here but not use it til we plot, but oh well)
            if pvkey not in plotvals:
                plotvals[pvkey] = initfcn()
            plotlist = plotvals[pvkey][ikey] if ikey is not None else plotvals[pvkey]  # it would be nice if the no-replicate-families-together case wasn't treated so differently
            plotlist.append((tau, fval))  # NOTE this appends to plotvals, the previous line is just to make sure we append to the right place
            if args.x_legend_var is not None:
                if 'mfreq' in args.x_legend_var:
                    mfreq = utils.get_mean_mfreq(pdirfcn(varnames, vstrs) + '/hmm')
                    mstr = ('%.0f' % (100*mfreq)) if '-pct' in args.x_legend_var else '%.3f' % mfreq
                elif args.x_legend_var == 'n-seqs':
                    assert ptilestr == 'time-reqd'  # might require changing something otherwise (?)
                    mstr = '%d' % get_n_seqs(vlists, varnames)
                else:
                    assert False
                if mstr not in xleg_vals:
                    xleg_vals[tau] = mstr
                else:
                    if mstr != xleg_vals[tau]:
                        print '  %s different values for derived var %s: %s vs %s' % (utils.color('yellow', 'warning'), args.x_legend_var, mstr, xleg_vals[tau])
        # ----------------------------------------------------------------------------------------
        def get_iclusts(yamlfo, yfname):
            assert per_x is not None  # just to make clear we don't get here for paired
            iclusts_in_file = []
            if 'percentiles' in yamlfo:  # if there's info for each cluster, it's in sub-dicts under 'iclust-N' (older files won't have it)
                iclusts_in_file = sorted([int(k.split('-')[1]) for k in yamlfo['percentiles']['per-seq'] if 'iclust-' in k])
            else:
                iclusts_in_file = sorted([int(k.split('-')[1]) for k in yamlfo if 'iclust-' in k])
            missing_iclusts = [i for i in range(nsimevts()) if i not in iclusts_in_file]
            if len(missing_iclusts) > 0:
                print '  %s missing %d/%d iclusts (i = %s) from file %s' % (utils.color('red', 'error'), len(missing_iclusts), nsimevts(), ' '.join(str(i) for i in missing_iclusts), yfname)
            # assert iclusts_in_file == list(range(nsimevts()))  # I'm not sure why I added this (presumably because I thought I might not see missing ones any more), but I'm seeing missing ones now (because clusters were smaller than min_selection_metric_cluster_size)
            return iclusts_in_file
        # ----------------------------------------------------------------------------------------
        def read_smetric_file(vlists, vstrs):
            if debug:
                dbgstr = get_obs_frac(vlists, varnames, return_dbg=True)
                print '%s   | %s' % (get_varval_str(vstrs), dbgstr),
            yfname = fnfcn(varnames, vstrs, metric, ptilestr, cg=choice_grouping, tv=lbplotting.ungetptvar(ptilestr), use_relative_affy=use_relative_affy, extra_str=metric_extra_str)
            try:
                with open(yfname) as yfile:
                    yamlfo = json.load(yfile)  # too slow with yaml
            except IOError:  # os.path.exists() is too slow with this many files
                missing_vstrs['missing'].append((None, vstrs))
                return
            if treat_clusters_together:
                add_plot_vals(yamlfo, vlists, varnames)
            else:
                for iclust in get_iclusts(yamlfo, yfname):
                    add_plot_vals(yamlfo, vlists, varnames, iclust=iclust)
        # ----------------------------------------------------------------------------------------
        def read_pairclust_file(vlists, vstrs):
            # ----------------------------------------------------------------------------------------
            if debug:
                print '%s   | %s' % (get_varval_str(vstrs), ''),
            yfname = fnfcn(varnames, vstrs)
            try:
                if ptilestr == 'time-reqd':
                    ytmpfo = readlog(yfname, metric, locus, ptntype)
                elif 'pcfrac-' in ptilestr or ptilestr == 'naive-hdist':
                    ytmpfo = read_hist_csv(yfname, ptilestr)
                else:
                    _, _, cpath = utils.read_output(yfname, skip_annotations=True)
                    ytmpfo = {ptilestr : cp_val(cpath, ptilestr, yfname)}
            except IOError:  # os.path.exists() is too slow with this many files
                missing_vstrs['missing'].append((None, vstrs))
                return
            add_plot_vals(ytmpfo, vlists, varnames)
        # ----------------------------------------------------------------------------------------
        if debug:
            if per_x is None:
                print '%s              values' % get_varname_str()
            else:
                print '%s   | obs times    N/gen        carry cap       fraction sampled           values' % get_varname_str()
        missing_vstrs = {'missing' : [], 'empty' : []}
        for vlists, vstrs in zip(val_lists, valstrs):  # why is this called vstrs rather than vstr?
            if per_x is None:
                read_pairclust_file(vlists, vstrs)
            else:
                read_smetric_file(vlists, vstrs)
            if debug:
                print ''

        # print info about missing and empty results
        n_printed, n_max_print = 0, 5
        for mkey, vstrs_list in missing_vstrs.items():  # ok now it's iclust and vstrs list, but what tf am I going to name that
            if len(vstrs_list) == 0:
                continue
            print '        %s: %d families' % (mkey, len(vstrs_list))
            print '     %s   iclust' % get_varname_str()
            for iclust, vstrs in vstrs_list:
                tfn = fnfcn(varnames, vstrs) if per_x is None else fnfcn(varnames, vstrs, metric, ptilestr, cg=choice_grouping, tv=lbplotting.ungetptvar(ptilestr), use_relative_affy=use_relative_affy, extra_str=metric_extra_str)  # arg this is ugly
                print '      %s    %4s    %s' % (get_varval_str(vstrs), iclust, tfn)
                n_printed += 1
                if n_printed >= n_max_print:
                    print '             [...]'
                    print '      skipping %d more lines' % (len(vstrs_list) - n_max_print)
                    break

        # average over the replicates/clusters
        if (args.n_replicates > 1 or not treat_clusters_together) and len(plotvals) > 0:
            if debug:
                print '  averaging over %d replicates' % args.n_replicates,
                if per_x is not None and nsimevts() is not None:
                    if treat_clusters_together:
                        print '(treating %d clusters per proc together)' % nsimevts(),
                    else:
                        print 'times %d clusters per proc:' % nsimevts(),
                print ''
                tmplen = str(max(len(pvkey) for pvkey in plotvals))
                print ('    %'+tmplen+'s   N used  N expected') % 'pvkey'
                dbgvals = []
            for pvkey, ofvals in plotvals.items():
                mean_vals, err_vals = [], []
                ofvals = {i : vals for i, vals in ofvals.items() if len(vals) > 0}  # remove zero-length ones (which should [edit: maybe?] correspond to 'missing'). Note that this only removes one where *all* the vals are missing, whereas if they're partially missing they values they do have will get added as usual below
                n_used = []  # just for dbg
                tmpvaldict = collections.OrderedDict()  # rearrange them into a dict keyed by the appropriate tau/xval
                for ikey in ofvals:  # <ikey> is an amalgamation of iseeds and icluster, e.g. '20-0'
                    for pairvals in ofvals[ikey]:
                        tau, tval = pairvals  # reminder: tau is not in general (any more) tau, but is the variable values fulfilling the original purpose of tau (i think x values?) in the plot
                        tkey = tuple(tau) if isinstance(tau, list) else tau  # if it's actually tau, it will be a single value, but if xvar is set to, say, n-sim-seqs-per-gen then it will be a list
                        if tkey not in tmpvaldict:  # these will usually get added in order, except when there's missing ones in some ikeys
                            tmpvaldict[tkey] = []
                        tmpvaldict[tkey].append(tval)
                use_sort = False # UPDATE i think we can just not sort? it means you have to have the order right on the command line # xvar != 'parameter-variances' and 'None' not in tmpvaldict.keys()  # we want None to be first
                tvd_keys = sorted(tmpvaldict) if use_sort else tmpvaldict.keys()  # for parameter-variances we want to to keep the original ordering from the command line
                for tau in tvd_keys:  # note that the <ltmp> for each <tau> are in general different if some replicates/clusters are missing or empty
                    ltmp = tmpvaldict[tau]  # len of <ltmp> is N seeds (i.e. procs) times N clusters per seed
                    mean_vals.append((tau, numpy.mean(ltmp)))
                    err_vals.append((tau, numpy.std(ltmp, ddof=1) / math.sqrt(len(ltmp))))  # standard error on mean (for standard deviation, comment out denominator)
                    n_used.append(len(ltmp))
                    if debug:
                        dbgvals.append((tau, mean_vals[-1][1], err_vals[-1][1]))
                plotvals[pvkey] = mean_vals
                errvals[pvkey] = err_vals
                if debug:
                    n_expected = args.n_replicates
                    if not treat_clusters_together:
                        n_expected *= nsimevts()
                    print ('     %'+tmplen+'s     %s     %4d%s') % (pvkey if pvkey!='' else '-', ('%4d' % n_used[0]) if len(set(n_used)) == 1 else utils.color('red', ' '.join(str(n) for n in set(n_used))), n_expected, '' if n_used[0] == n_expected else utils.color('red', ' <--'))
            if debug:
                print '    final values:'
                for tau, meanval, errval in dbgvals:
                    print '     %6s  %5.2f +/- %.1f' % (tau, meanval, errval)
    # ----------------------------------------------------------------------------------------
    def plotcall(pvkey, xticks, diffs_to_perfect, yerrs, mtmp, ipv=None, imtmp=None, label=None, add_to_leg=False, alpha=0.5, estr=''):
        markersize = 15  # 1 if len(xticks) > 1 else 15
        linestyle = linestyles.get(mtmp, '-')
        if args.plot_metrics.count(mtmp) > 1 and estr != '':
            linestyle = 'dotted'
        linewidth = linewidths.get(mtmp, 3)
        color = None
        if ipv is not None:
            color = plotting.frozen_pltcolors[ipv % len(plotting.frozen_pltcolors)]
        elif imtmp is not None:  # used to us <imtmp> to directly get color, but now we want to get the same colors no matter the matplotlib version and order on the command line, so now it just indicates that we should add the metric str
            # color = plotting.frozen_pltcolors[imtmp % len(plotting.pltcolors)]
            color = metric_color(mtmp)
        if yerrs is not None:
            ax.errorbar(xticks, diffs_to_perfect, yerr=yerrs, color=color, alpha=alpha, linewidth=linewidth, markersize=markersize, marker='.', linestyle=linestyle)  #, title='position ' + str(position)) #, label=legstr(label)
        else:
            ax.plot(xticks, diffs_to_perfect, color=color, alpha=alpha, linewidth=linewidth)  # , label=legstr(label)
        if add_to_leg:
            dlabel = mtmp if len(pvk_list) == 1 else label
            if not args.dont_plot_extra_strs and estr != '':
                dlabel += ' %s' % estr
            # ax.plot([], [], label=legstr(dlabel), alpha=alpha, linewidth=linewidth, linestyle=linestyle, color='grey' if ipv is not None else color, marker='.', markersize=0)
            leg_entries[legstr(dlabel)] = {'alpha' : alpha, 'linewidth' : linewidth, 'linestyle' : linestyle, 'color' : 'grey' if (ipv is not None and len(pvk_list)==1) else color} #, marker='.', markersize=0)
        for dtp in diffs_to_perfect:  # NOTE can't add all at once bc of namespace issues
            all_yvals.add(dtp)
    # ----------------------------------------------------------------------------------------
    def getplotname(mtmp):
        if per_x is None:
            return '%s-%s-vs-%s-%s-%s' % (ptilestr, mtmp if mtmp is not None else 'combined', xvar, ptntype, locus)
        elif per_x == 'per-seq':
            return '%s%s-%s-ptiles-vs-%s-%s' % (affy_key_str, ptilestr, mtmp if mtmp is not None else 'combined', xvar, choice_grouping)
        else:
            return '%s-ptiles-vs-%s' % (choice_grouping.replace('-vs', ''), xvar)
    # ----------------------------------------------------------------------------------------
    def getkwargs(estr):
        if per_x is None:
            return {'perf_metric' : ptilestr}
        else:
            return {'per_x' : per_x, 'extra_str' : estr}
    # ----------------------------------------------------------------------------------------
    def get_outfname(mtmp, estr):
        return '%s/%s.yaml' % (get_comparison_plotdir(args, mtmp, **getkwargs(estr)), getplotname(mtmp))
    # ----------------------------------------------------------------------------------------
    def titlestr(metric):
        if per_x is None:
            return '%s' % plotting.legends.get(metric, metric)
        else:
            return '%s:' % lbplotting.mtitlestr(per_x, metric, short=True, max_len=7)
    # ----------------------------------------------------------------------------------------
    def getxticks(xvals):
        l_xvar = xvar
        if args.x_legend_var is not None:
            if action == 'plot':  # if combining plots, we're reading files in which the conversion already happened before writing
                xvals = [xleg_vals[x] for x in xvals]
            l_xvar = args.x_legend_var
        xlabel = legdict.get(l_xvar, l_xvar.replace('-', ' '))
        if l_xvar == 'parameter-variances':  # special case cause we don't parse this into lists and whatnot here
            xticks, xticklabels = [], []
            global_pv_vars = None
            for ipv, pv_cft_str in enumerate(xvals):  # <pv_cft_str> corresponds to one bcr-phylo run, but can contain more than one parameter variance specification
                xticks.append(ipv)
                pv_vars, xtl_strs = [], []
                for pvar_str in pv_cft_str.split('_c_'):
                    assert '..' in pvar_str  # don't handle the uniform-distribution-with-variance method a.t.m.
                    pvar, pvals = pvar_str.split(',')
                    def fmt(v, single=False):
                        if pvar == 'selection-strength':
                            if v == 1.: fstr = '%.0f'
                            else: fstr = '%.2f' if single else '%.1f'
                        else:
                            fstr = '%d'
                        return fstr % v
                    pv_vars.append(pvar)
                    pvlist = [float(pv) for pv in pvals.split('..')]
                    xtlstr = '%s-%s'%(fmt(min(pvlist)), fmt(max(pvlist))) if min(pvlist) != max(pvlist) else fmt(pvlist[0], single=True)
                    xtl_strs.append(xtlstr)
                xticklabels.append('\n'.join(xtl_strs))
                if global_pv_vars is None:
                    global_pv_vars = pv_vars
                if pv_vars != global_pv_vars:
                    raise Exception('each bcr-phylo run has to have the same variables with parameter variances, but got %s and %s' % (global_pv_vars, pv_vars))
            xlabel = ', '.join(legdict.get(p, p.replace('-', ' ')) for p in global_pv_vars)
        elif isinstance(xvals[0], tuple) or isinstance(xvals[0], list):  # if it's a tuple/list (not sure why it's sometimes one vs other times the other), use (more or less arbitrary) integer x axis values
            def tickstr(t):
                if len(t) < 4:
                    return ', '.join(str(v) for v in t)
                else:
                    return '%s -\n %s\n(%d)' % (t[0], t[-1], len(t)) #, t[1] - t[0])
            xticks = list(range(len(xvals)))
            xticklabels = [tickstr(t) for t in xvals]
        elif None in xvals or 'None' in xvals:  # e.g. mean-cells-per-droplet has to handle mixed none/str + float values
            xticks = list(range(len(xvals)))  # arbitrarily put the ticks at integer vals starting with 0
            xticklabels = ['n/a' if t=='None' else str(t) for t in xvals]  # labels are the actual xvals tho
        else:
            xticks = [float(x) for x in xvals]  # i guess we can just float() all of them (and ignore svartypes)?
            xticklabels = [str(t) for t in xvals]
        if l_xvar in ['n-seqs', 'n-sim-events']:
            # ax.ticklabel_format(style='sci')
            for itick, (xt, xtl) in enumerate(zip(xticks, xticklabels)):
                # xticklabels[itick] = '10^{%d}' % math.log(xt)
                if xtl[-3:] == '000':
                    xticklabels[itick] = xticklabels[itick][:-3] + 'k'
        return xticks, xticklabels, xlabel

    # ----------------------------------------------------------------------------------------
    if per_x is None:  # paired clustering
        assert hasattr(args, 'n_sim_events_list') and not hasattr(args, 'n_sim_events_per_proc')
        distr_hists, affy_key_str, treat_clusters_together = False, '', True  # it's important that treat_clusters_together is True, but the others i think aren't used for paired clustering
        legdict = plotting.legends
        if args.x_legend_var is not None:
            xleg_vals = {}
    else:  # tree metrics
        assert not hasattr(args, 'n_sim_events_list') and hasattr(args, 'n_sim_events_per_proc')
        assert choice_grouping is not None
        distr_hists = ptilestr == 'n-ancestor'  # hackey, but ok for now
        if metric == 'lbr' and args.dont_observe_common_ancestors:
            print '    skipping lbr when only observing leaves'
            return
        affy_key_str = 'relative-' if (use_relative_affy and ptilestr=='affinity') else ''  # NOTE somewhat duplicates lbplotting.rel_affy_str()
        treat_clusters_together = nsimevts() is None or (per_x == 'per-seq' and choice_grouping == 'among-families')  # if either there's only one family per proc, or we're choosing cells among all the clusters in a proc together, then things here generally work as if there were only one family per proc (note that I think I don't need the 'per-seq' since it shouldn't be relevant for 'per-cluster', but it makes it clearer what's going on)
        legdict = treeutils.legtexts
        legdict.update(lbplotting.metric_for_target_distance_labels)
    pvlabel = ['?']  # arg, this is ugly (but it does work...)
    _, varnames, val_lists, valstrs = utils.get_var_info(args, svars)
    plotvals, errvals = collections.OrderedDict(), collections.OrderedDict()
    fig, ax = plotting.mpl_init()
    all_xtks, all_xtls, xlabel, all_yvals = [], [], None, set()
    leg_entries = collections.OrderedDict()
    if action == 'plot':
        read_plot_info()
        outfo = []
        if len(plotvals) == 0:
            print '  %s no plotvals for %s %s %s' % (utils.color('yellow', 'warning'), metric, per_x, choice_grouping)
            return
        for ipv, pvkey in enumerate(plotvals):
            xvals, diffs_to_perfect = zip(*plotvals[pvkey])
            xticks, xticklabels, xlabel = getxticks(xvals)
            all_xtks += [x for x in xticks if x not in all_xtks]
            all_xtls += [l for l in xticklabels if l not in all_xtls]
            # assert xvals == tuple(sorted(xvals))  # this definitely can happen, but maybe not atm? and maybe not a big deal if it does. So maybe should remove this
            yerrs = zip(*errvals[pvkey])[1] if pvkey in errvals else None  # each is pairs tau, err
            if args.x_legend_var is not None:
                xvals = [xleg_vals[x] for x in xvals]
            plotcall(pvkey, xticks, diffs_to_perfect, yerrs, metric, ipv=ipv, label=pvkey, estr=metric_extra_str)
            outfo.append((pvkey, {'xvals' : xvals, 'yvals' : diffs_to_perfect, 'yerrs' : yerrs}))
        with open(get_outfname(metric, metric_extra_str), 'w') as yfile:  # write json file to be read by 'combine-plots'
            json.dump(outfo, yfile)
        title = titlestr(metric)
        plotdir = get_comparison_plotdir(args, metric, **getkwargs(metric_extra_str))  # per_x=per_x, extra_str=metric_extra_str
    elif action == 'combine-plots':
        pvks_from_args = set([pvkeystr(vlists, varnames) for vlists in val_lists])  # have to call this fcn at least once just to set pvlabel (see above) [but now we're also using the results below UPDATE nvmd didn't end up doing it that way, but I'm leaving the return value there in case I want it later]
        plotfos = collections.OrderedDict()
        for mtmp, estr in zip(args.plot_metrics, args.plot_metric_extra_strs):  # <mtmp>: metric for trees, method for paired
            if per_x is not None and ptilestr not in [v for v, l in lbplotting.single_lbma_cfg_vars(mtmp, final_plots=True)]:  # i.e. if the <ptilestr> (variable name) isn't in any of the (variable name, label) pairs (e.g. n-ancestor for lbi; we need this here because of the set() in the calling block)
                continue
            ofn = get_outfname(mtmp, estr)
            if not os.path.exists(ofn):
                print '    %s missing %s' % (utils.color('yellow', 'warning'), ofn)
                continue
            with open(ofn) as yfile:
                mkey = mtmp
                if estr != '':
                    mkey = '%s%s%s' % (mtmp, xdelim, estr)  # this is ugly, but we need to be able to split it apart in the loop just below here
                plotfos[mkey] = collections.OrderedDict(json.load(yfile))
                if len(plotfos[mkey]) == 0:
                    raise Exception('read zero length info from %s' % ofn)  # if this happens when we're writing the file (above), we can skip it, but  I think we have to crash here (just rerun without this metric/extra_str). It probably means you were missing the dtr files for this per_x/cgroup
        if len(plotfos) == 0:
            print '  nothing to plot'
            return
        pvks_from_file = set([tuple(pfo.keys()) for pfo in plotfos.values()])  # list of lists of pv keys (to make sure the ones from each metric's file are the same)
        if len(pvks_from_file) > 1:  # eh, they can be different now if I ran different metrics with different argument lists
            print '  %s different lists of pv keys for different metrics: %s' % (utils.color('yellow', 'warning'), pvks_from_file)
            pvk_list = sorted(pvks_from_file, key=len)[0]  # use the shortest one
        else:
            pvk_list = list(pvks_from_file)[0]
        if args.pvks_to_plot is not None:
            # pvk_list = [p for p in list(pvks_from_file)[0] if p in pvks_from_args]  # don't do it this way since if you only ask it to plot one value it'll get the wrong file path (since it'll no longer make a subdir level for that variable)
            ptmp = [p for p in pvk_list if p in args.pvks_to_plot]
            if len(ptmp) == 0:
                raise Exception('requirement in --pvks-to-plot \'%s\' didn\'t give us any from the list %s' % (args.pvks_to_plot, pvk_list))
            pvk_list = ptmp
        for ipv, pvkey in enumerate(pvk_list):
            for imtmp, (mkey, pfo) in enumerate(plotfos.items()):
                mtmp, estr = (mkey, '') if xdelim not in mkey else mkey.split(xdelim)
                xticks, xticklabels, xlabel = getxticks(pfo[pvkey]['xvals'])
                plotcall(pvkey, xticks, pfo[pvkey]['yvals'], pfo[pvkey]['yerrs'], mtmp, label=pvkey if (imtmp == 0 and len(pvk_list) > 1) else None, ipv=ipv if len(pvk_list) > 1 else None, imtmp=imtmp, add_to_leg=ipv==0 or len(pvk_list)>1, estr=estr)
                all_xtks += [x for x in xticks if x not in all_xtks]
                all_xtls += [l for l in xticklabels if l not in all_xtls]
        # if ''.join(args.plot_metric_extra_strs) == '':  # no extra strs
        #     title = '+'.join(plotfos) + ': '
        # else:
        #     title = '+'.join(set(args.plot_metrics)) + ': '
        title = ''
        plotdir = get_comparison_plotdir(args, 'combined', per_x=per_x)
    else:
        assert False

    xmin, xmax = 0.95 * min(all_xtks), 1.05 * max(all_xtks)
    ymin, ymax = ax.get_ylim()
    # if ptilestr == 'affinity':
    #     ymin, ymax = 0, max(ymax, 25)
    # elif ptilestr == 'n-ancestors':
    #     ymin, ymax = 0, max(ymax, 1.5)

    log, adjust = '', {}
    if xvar == 'lb-tau' and len(all_xtks) > 1:
        ax.plot([1./args.seq_len, 1./args.seq_len], (ymin, ymax), linewidth=3, alpha=0.7, color='darkred', linestyle='--') #, label='1/seq len')
    if xvar in ['carry-cap', 'n-sim-events', 'n-leaves']:
        log = 'x'

    if ax.get_ylim()[1] < 1:
        adjust['left'] = 0.21
    if ax.get_ylim()[1] < 0.01:
        adjust['left'] = 0.26
    adjust['bottom'] = 0.25
    adjust['top'] = 0.9
    if all_xtls is not None and '\n' in all_xtls[0]:
        adjust['bottom'] = 0.3
        import matplotlib.pyplot as plt
        plt.xlabel('xlabel', fontsize=14)

    n_ticks = 4
    ylabel, leg_loc, yticks, yticklabels = '', None, None, None
    if per_x is None:
        ylabel = ptilelabel
        if ptilestr == 'time-reqd':
            if ltexts[ptntype] == 'joint':
                title += 'joint merge time'
            else:
                title += ' %s: single chain partition time' % locus
            ymin, ymax = [mfcn(sorted(all_yvals)) for mfcn in [min, max]]
            if 'y' not in log: log += 'y'
            yticks, yticklabels = plotting.timeticks, plotting.timeticklabels
            in_ticks = [y for y in yticks if y > ymin and y < ymax]  # ticks that will actually show up
            if ymin < in_ticks[0]:  # if the lowest tick isn't right at the min, add another tick below
                ymin = yticks[yticks.index(in_ticks[0]) - 1]  # don't think there's any reason to add it to in_ticks
            if ymax > in_ticks[-1]:  # same for max
                ymax = yticks[yticks.index(in_ticks[-1]) + 1]
        else:
            if 'pcfrac-' in ptilestr:
                title = ''
            else:
                title += ' %s: %s %s' % (locus, ltexts[ptntype], legdict.get(ptilestr, ptilestr))
            ymin = 0
            if ptilestr != 'naive-hdist':
                ymax = 1.05
        leg_loc = [0.7, 0.15]
        if args.final_plot_xvar == 'n-leaves' and '--constant-number-of-leaves' in args.simu_extra_args:
            xlabel += ' (constant)'
    else:
        # dy = (ymax - ymin) / float(n_ticks - 1)
        # if ptilestr != 'affinity':
        #     yticks = [int(y) if ptilestr == 'affinity' else utils.round_to_n_digits(y, 3) for y in numpy.arange(ymin, ymax + 0.5*dy, dy)]
        #     yticklabels = ['%s'%y for y in yticks]
        if per_x == 'per-seq':
            title += 'choosing %s' % (choice_grouping.replace('within-families', 'within each family').replace('among-', 'among all '))
        if use_relative_affy:
            fig.text(0.5, 0.87, 'relative %s' % ptilestr, fontsize=15, color='red', fontweight='bold')
        leg_loc = [0.04, 0.6]
        # if metric != 'lbi' and len(title) < 17:
        #     leg_loc[0] = 0.7
        if distr_hists:
            ylabel = 'frac. correct'
            if args.distr_hist_limit[0] == 'val':
                ylabel += ' (%s > %.1f)' % (titlestr(metric), args.distr_hist_limit[1])
            elif args.distr_hist_limit[0] == 'frac':
                ylabel += ' (top %.0f%%)' % (100*args.distr_hist_limit[1])
            elif args.distr_hist_limit[0] == 'N':
                ylabel += ' (top %d per family)' % args.distr_hist_limit[1]
            else:
                assert False
        else:
            ylabel = '%s to perfect' % ('percentile' if ptilelabel == 'affinity' else ptilelabel)

    ffn = plotting.mpl_finish(ax, plotdir, getplotname(metric),
                              xlabel=xlabel,
                              # ylabel='%s to perfect\nfor %s ptiles in [%.0f, 100]' % ('percentile' if ptilelabel == 'affinity' else ptilelabel, ylabelstr, min_ptile_to_plot),
                              ylabel=ylabel,
                              title=title,  # leg_title=legstr(pvlabel[0], title=True), leg_prop={'size' : 12}, leg_loc=leg_loc,
                              xticks=all_xtks, xticklabels=all_xtls, xticklabelsize=12 if all_xtls is not None and '\n' in all_xtls[0] else 16,
                              yticks=yticks, yticklabels=yticklabels,
                              xbounds=(xmin, xmax), ybounds=(ymin, ymax), log=log, adjust=adjust,
    )
    if make_legend:
        leg_entries = collections.OrderedDict(reversed(list(leg_entries.items())))  # the last things to be *plotted* cover earlier things, but we also want them to be at the top of the legend
        plotting.plot_legend_only(leg_entries, plotdir, 'legend', title=legstr(pvlabel[0], title=True))  #[(l, leg_entries['color']) for l, lfo in leg_entries], )
    if fnames is not None:
        fnames.append(ffn)
