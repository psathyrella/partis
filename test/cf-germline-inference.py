#!/usr/bin/env python
import yaml
import glob
import csv
import math
import copy
import numpy
from collections import OrderedDict
import os
import random
import argparse
import sys
import subprocess
import colored_traceback.always
sys.path.insert(1, './python')

import utils
import glutils
from hist import Hist
sys.path.insert(1, './datascripts')
import heads

sim_locus = 'igh'
region = 'v'

study_translations = {  # for zenodo
    'jason-mg' : 'myasthenia-gravis',
    'jason-influenza' : 'influenza',
}
varval_titles = {
    'mfreq' : 'SHM rate',
    'nsnp' : 'SNPs to new allele',
    'multi-nsnp' : 'SNPs to multiple new alleles',
    'prevalence' : 'new-allele prevalence',
    'n-leaves' : 'mean leaves per clonal family',
    'weibull' : 'tree balance',
}

varval_legend_titles = {
    'mfreq' : 'SHM',
    'nsnp' : 'N SNPs',
    'multi-nsnp' : 'N SNPs',
    'prevalence' : 'prevalence',
    'n-leaves' : 'N leaves',
    'weibull' : 'Weibull parameter',
}

all_methods = ['tigger-default', 'igdiscover', 'partis']
characters = ['subject', 'isotype', 'timepoint']
def methstr(meth):
    if meth == 'tigger-default':
        return 'TIgGER'
    if meth == 'igdiscover':
        return 'IgDiscover'
    elif meth == 'full':
        return 'full IMGT'
    else:
        return meth
def diffstr(difficulty):
    if difficulty == 'easy':
        return 'low SHM'
    elif difficulty == 'hard':
        return 'high SHM'
    else:
        assert False
def gls_sim_str(diff, iproc):
    return '%s%s' % (diffstr(diff), (' (%s)' % str(iproc)) if iproc != '' else '')  # cast to string 'cause sometimes it isn't a string (e.g. 'total')

# ----------------------------------------------------------------------------------------
def varvalstr(name, val):
    if name == 'multi-nsnp':
        valstr = ':'.join([str(v) for v in val])
    elif name == 'alcluster':
        valstr = '-'.join(['n%s-%s' % (k, v) for k, v in val.items()])
    elif name == 'gls-gen':
        valstr = 'simu'
    else:
        valstr = str(val)
    return valstr

# ----------------------------------------------------------------------------------------
def legend_str(args, val):
    if args.action == 'mfreq':
        # lstr = '%.1fx' % val
        if val < 0.15:
            return 'low'
        elif val > 0.5 and val < 1.5:
            return 'typical'
        elif val > 1.5:
            return 'high'
        else:
            assert False
    elif args.action == 'nsnp':
        lstr = '%d' % val
    elif args.action == 'multi-nsnp':
        lstr = '%s' % '+'.join([str(v) for v in val])
    elif args.action == 'prevalence':
        lstr = '%d%%' % (100*val)
    elif args.action == 'n-leaves':
        lstr = '%.1f' % val
    elif args.action == 'weibull':
        lstr = '%.1f' % val
    elif args.action == 'alcluster':
        lstr = 'er...'
    else:
        assert False
    return lstr

# ----------------------------------------------------------------------------------------
def get_outdir(args, baseoutdir, varname, varval, n_events=None):
    outdir = baseoutdir
    if args.action == 'gls-gen':
        outdir += '/' + varvalstr(varname, varval) + '/' + args.gls_gen_difficulty
    if args.action == 'data':
        outdir += '/' + varvalstr(varname, varval)
    else:
        outdir += '/' + varname + '-' + varvalstr(varname, varval) + '/n-events-' + str(n_events)
    return outdir

# ----------------------------------------------------------------------------------------
def get_single_performance(outdir, method, debug=False):
    sglfo = glutils.read_glfo(outdir + '/germlines/simulation', locus=sim_locus)
    iglfo = glutils.read_glfo(outdir + '/' + method + '/sw/germline-sets', locus=sim_locus)
    glutils.synchronize_glfos(ref_glfo=sglfo, new_glfo=iglfo, region=region)
    missing_alleles = set(sglfo['seqs'][region]) - set(iglfo['seqs'][region])
    spurious_alleles = set(iglfo['seqs'][region]) - set(sglfo['seqs'][region])
    if debug:
        if len(missing_alleles) > 0:
            print '    %2d  missing %s' % (len(missing_alleles), ' '.join([utils.color_gene(g) for g in missing_alleles]))
        if len(spurious_alleles) > 0:
            print '    %2d spurious %s' % (len(spurious_alleles), ' '.join([utils.color_gene(g) for g in spurious_alleles]))
        if len(missing_alleles) == 0 and len(spurious_alleles) == 0:
            print '    none missing'
    return {
        'missing' : len(missing_alleles),
        'spurious' : len(spurious_alleles),
        'total' : len([g for g in sglfo['seqs'][region] if '+' in g]),  # anybody with a '+' should be a new allele
    }

# ----------------------------------------------------------------------------------------
def get_gls_fname(outdir, method, locus, sim_truth=False, data=False, annotation_performance_plots=False):  # NOTE duplicates/depends on code in test-germline-inference.py
    if annotation_performance_plots:
        return outdir + '/' + method + '/annotation-performance-plots/sw/mutation'
    gls_dir = get_gls_dir(outdir, method, sim_truth=sim_truth, data=data, annotation_performance_plots=annotation_performance_plots)
    return glutils.get_fname(gls_dir, locus, region)

# ----------------------------------------------------------------------------------------
def get_gls_dir(outdir, method, sim_truth=False, data=False, annotation_performance_plots=False):  # NOTE duplicates/depends on code in test-germline-inference.py
    if data:
        if method == 'partis' or method == 'full':
            outdir += '/hmm/germline-sets'  # NOTE this is inside the datascripts output dir, also NOTE doesn't use <method> (since we only have partis for a method a.t.m., although could use --label or --extra-str to differentiate)
        else:
            outdir += '/' + method
    elif sim_truth:
        outdir += '/germlines/simulation'
    elif method == 'partis' or method == 'full':
        outdir += '/' + method + '/sw/germline-sets'
    elif 'tigger' in method or method == 'igdiscover':
        outdir += '/' + method
    else:
        assert False
    return outdir

# ----------------------------------------------------------------------------------------
def make_gls_tree_plot(args, plotdir, plotname, glsfnames, glslabels, locus, ref_label=None, title=None, title_color=None, legends=None, legend_title=None, pie_chart_faces=False):
    # ete3 requires its own python version, so we run as a subprocess
    cmdstr = 'export PATH=%s:$PATH && xvfb-run -a ./bin/plot-gl-set-trees.py' % args.ete_path
    cmdstr += ' --plotdir ' + plotdir
    cmdstr += ' --plotname ' + plotname
    cmdstr += ' --glsfnames ' + ':'.join(glsfnames)
    cmdstr += ' --glslabels ' + ':'.join(glslabels)
    if ref_label is not None:
        cmdstr += ' --ref-label ' + ref_label
    if title is not None:
        cmdstr += ' --title="%s"' % title
    if title_color is not None:
        cmdstr += ' --title-color %s' % title_color
    if legends is not None:
        cmdstr += ' --legends=' + ':'.join('"%s"' % l for l in legends)
    if legend_title is not None:
        cmdstr += ' --legend-title="%s"' % legend_title
    if pie_chart_faces:
        cmdstr += ' --pie-chart-faces'
    cmdstr += ' --locus ' + locus
    if args.plotcache:
        cmdstr += ' --use-cache'
    if args.only_print:
        cmdstr += ' --only-print'
    utils.simplerun(cmdstr, shell=True, debug=args.dry_run, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def print_gls_gen_summary_table(args, baseoutdir):
    latex = True
    varname = args.action
    varval = 'simu'
    statuses = ['missing', 'spurious', 'ok']
    meanvals = {m : {st : [] for st in statuses} for m in args.methods}
    for method in args.methods:
        for iproc in range(args.iteststart, args.n_tests):
            resultfile = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/' + str(iproc) + '/' + method + '/gls-gen-plots/results.yaml'
            with open(resultfile) as yamlfile:
                yamlfo = yaml.safe_load(yamlfile)
            for status in statuses:
                meanvals[method][status].append(len(yamlfo[status]))

    print '%20s     %s' % ('', ' '.join([('%-16s' % st) for st in statuses]))
    for meth in args.methods:
        print '%20s' % methstr(meth),
        for status in statuses:
            mean = float(sum(meanvals[meth][status])) / len(meanvals[meth][status])
            err = numpy.std(meanvals[meth][status], ddof=1) / math.sqrt(len(meanvals[meth][status]))
            print '  %s   %s%4.1f %s %-4.1f%s' % ('&' if latex else '', '$' if latex else '', mean, '\\pm' if latex else '+/-', err, '$' if latex else ''),
        print '%s' % ('\\\\' if latex else '')

# ----------------------------------------------------------------------------------------
def get_gls_gen_annotation_performance_plots(args, baseoutdir):
    import plotting
    import plotconfig
    methcolors = {  # NOTE started from scolors in bin/plot-gl-set-trees.py htmlcolorcods.com, and slide each one a little rightward
        'tigger-default' : '#dd4d39',
        'igdiscover' : '#55ab7a', #60ac84',
        'partis' : '#6b83ca', #758bcd',
        'full' : '#858585',
    }
    lstyledict = {}  # 'tigger-default' : '--'}
    lwdict = {'full' : 9, 'igdiscover' : 8, 'partis' : 5, 'tigger-default' : 2}  # methods are sorted below, so it's always [full, igdiscover, partis, tigger]
    linewidths = [lwdict[m] for m in args.methods]
    colors = [methcolors[meth] for meth in args.methods]
    linestyles = [lstyledict.get(m, '-') for m in args.methods]
    alphas = [0.8 if m in ['full', 'igdiscover'] else 1 for m in args.methods]

    varname = args.action
    varval = 'simu'
    plotnames = ['v_hamming_to_true_naive', 'v_muted_bases']
    xtitles = ['V distance to true naive', 'inferred - true']
    meanvals = {pn : {m : [] for m in args.methods} for pn in plotnames}
    print '  annotations: %s' % get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events)
    all_hists = {pn : [] for pn in plotnames}
    for iproc in range(args.iteststart, args.n_tests):
        outdir = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/' + str(iproc)  # duplicates code in bin/test-germline-inference.py
        plotdir = outdir + '/annotation-performance-plots'
        print '    %s' % plotdir
        if not args.only_print:
            utils.prep_dir(plotdir, wildlings=['*.png', '*.svg', '*.csv'])

        # shenanigans for the six (three easy and thre hard) of 'em that go in the paper pdf
        make_legend = (iproc > 2) or (iproc == 0) # and args.gls_gen_difficulty == 'easy')
        make_xtitle = (iproc > 2) or (iproc == 2)
        make_ytitle = (iproc > 2) or (args.gls_gen_difficulty == 'easy')

        for plotname in plotnames:
            hists = {meth : Hist(fname=get_gls_fname(outdir, meth, sim_locus, annotation_performance_plots=True) + '/' + plotname + '.csv', title=methstr(meth) if make_legend else None) for meth in args.methods}
            for meth in args.methods:
                if hists[meth].overflow_contents() != 0.0:
                    print '  %s %s non-zero under/overflow %f' % (utils.color('red', 'error'), methstr(meth), hists[meth].overflow_contents())
                meanvals[plotname][meth].append(hists[meth].get_mean())
            if args.only_print:
                continue
            plotting.draw_no_root(hists[args.methods[0]], log='y', plotdir=plotdir, plotname=plotname, more_hists=[hists[m] for m in args.methods[1:]], colors=colors, ytitle='sequences' if make_ytitle else None,
                                  xtitle=xtitles[plotnames.index(plotname)] if make_xtitle else '',
                                  plottitle=gls_sim_str(args.gls_gen_difficulty, iproc),
                                  linewidths=linewidths, linestyles=linestyles, alphas=alphas, remove_empty_bins=True, square_bins=True)
            all_hists[plotname].append(hists)

    print '  total plots'
    plotdir = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/annotation-performance-plots'
    print '    %s' % plotdir
    if not args.only_print:
        utils.prep_dir(plotdir, wildlings=['*.png', '*.svg', '*.csv'])
    for plotname in plotnames:
        total_hists = {}
        for meth in args.methods:
            xmin = min([hdict[meth].xmin for hdict in all_hists[plotname]])
            xmax = max([hdict[meth].xmax for hdict in all_hists[plotname]])
            total_hists[meth] = Hist(xmax - xmin, xmin, xmax, title=all_hists[plotname][0][meth].title)
            for hdict in all_hists[plotname]:
                assert hdict[meth].integral(include_overflows=True) > 100  # make sure it isn't normalized (this is a shitty way to do this)
                bin_centers = hdict[meth].get_bin_centers()
                for ibin in range(len(hdict[meth].low_edges)):
                    xval = bin_centers[ibin]
                    for _ in range(int(hdict[meth].bin_contents[ibin])):
                        total_hists[meth].fill(xval)
        plotting.draw_no_root(total_hists[args.methods[0]], log='y', plotdir=plotdir, plotname='total-' + plotname, more_hists=[total_hists[m] for m in args.methods[1:]], colors=colors, ytitle='sequences' if make_ytitle else None,
                              xtitle=xtitles[plotnames.index(plotname)],
                              plottitle=gls_sim_str(args.gls_gen_difficulty, iproc=''),
                              linewidths=linewidths, linestyles=linestyles, alphas=alphas, remove_empty_bins=True, square_bins=True)

    for plotname in plotnames:
        if 'muted_bases' in plotname:  #  mean value isn't meaningful
            continue
        print plotname
        for meth in args.methods:
            mean = float(sum(meanvals[plotname][meth])) / len(meanvals[plotname][meth])
            err = numpy.std(meanvals[plotname][meth], ddof=1) / math.sqrt(len(meanvals[plotname][meth]))
            print '   %15s  %6.3f / %d = %6.2f +/- %6.2f' % (methstr(meth), sum(meanvals[plotname][meth]), len(meanvals[plotname][meth]), mean, err)

# ----------------------------------------------------------------------------------------
def get_gls_gen_tree_plots(args, baseoutdir, method):
    varname = args.action
    varval = 'simu'
    for iproc in range(args.iteststart, args.n_tests):
        outdir = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/' + str(iproc)
        print '%-2d                            %s' % (iproc, outdir)
        simfname = get_gls_fname(outdir, method=None, locus=sim_locus, sim_truth=True)
        inffname = get_gls_fname(outdir, method, sim_locus)
        make_gls_tree_plot(args, outdir + '/' + method + '/gls-gen-plots', varvalstr(varname, varval),
                           glsfnames=[simfname, inffname],
                           glslabels=['sim', 'inf'],
                           locus=sim_locus, ref_label='sim', legend_title=methstr(method), title=gls_sim_str(args.gls_gen_difficulty, iproc))  #, title_color=method)

# ----------------------------------------------------------------------------------------
def get_character_str(character, charval):
    if character == 'subject':
        return charval
    elif character == 'timepoint':
        if charval == 'merged':
            return ''
        elif charval[0] == 'W':
            return 'week %d' % int(charval[1:])
        elif charval[0] == 'M':
            return 'month %d' % int(charval[1:])
        elif 'dpi' in charval:
            return charval.replace('dpi', ' dpi')
        elif charval[0] in ['p', 'm'] and charval[-1] in ['h', 'd']:
            plusminusstr = charval[0].replace('p', '+').replace('m', '-')
            number = int(charval[1:-1])
            unitstr = charval[-1].replace('h', 'hour').replace('d', 'day')
            return '%s%d %s%s' % (plusminusstr, number, unitstr, utils.plural(number))
        else:
            raise Exception('not sure what to do with %s' % charval)
    elif character == 'isotype':
        return 'Ig%s' % charval.upper()
    else:
        assert False

# ----------------------------------------------------------------------------------------
def all_the_same(chstr, mfolist):
    if chstr not in mfolist[0]:
        raise Exception('column \'%s\' not available for %s (choices: %s)' % (chstr, nstr, ' '.join(mfolist[0].keys())))
    return len(set([mfo[chstr] for mfo in mfolist])) == 1

# ----------------------------------------------------------------------------------------
def get_dset_title(mfolist):
    return_strs = []
    for character in characters:
        if mfolist[0][character] == '':
            continue
        if len(mfolist) == 1 or all_the_same(character, mfolist):
            return_strs.append('%s' % get_character_str(character, mfolist[0][character]))
    return ' '.join(return_strs)

# ----------------------------------------------------------------------------------------
def get_dset_legends(mfolist):
    assert len(mfolist) > 1
    legends = []
    for mfo in mfolist:
        return_strs = []
        for character in characters:
            if mfo[character] == '':
                continue
            if not all_the_same(character, mfolist):
                return_strs.append('%s' % get_character_str(character, mfo[character]))
        legends.append(' '.join(return_strs) + ' only')
    return legends

# ----------------------------------------------------------------------------------------
def get_data_plots(args, baseoutdir, methods, study, dsets):
    metafos = heads.read_metadata(study)
    assert len(set([metafos[ds]['locus'] for ds in dsets]))  # make sure everybody has the same locus
    mfo = metafos[dsets[0]]
    data_outdirs = [heads.get_datadir(study, 'processed', extra_str='gls-gen-paper-' + args.label) + '/' + ds for ds in dsets]
    outdir = get_outdir(args, baseoutdir, varname='data', varval=study + '/' + '-vs-'.join(dsets))  # for data, only the plots go here, since datascripts puts its output somewhere else
    if len(dsets) > 1 and len(methods) == 1:  # sample vs sample
        glslabels = dsets
        title = get_dset_title([metafos[ds] for ds in dsets])
        if study != 'kate-qrs':
            title += '  %s' % methstr(methods[0])
        title_color = methods[0]
        legends = get_dset_legends([metafos[ds] for ds in dsets])
        legend_title = methstr(methods[0]) if study == 'kate-qrs' else None  # for kate-qrs we need to put the subject _and_ the isotype in the title, so there's no room for the method
        pie_chart_faces = False
        print '%s:' % utils.color('green', methods[0]),
    elif len(methods) > 1 and len(dsets) == 1:  # method vs method
        glslabels = methods
        title = get_dset_title([mfo])
        title_color = None
        legends = [methstr(m) + ' only' for m in methods]
        legend_title = None
        pie_chart_faces = len(methods) > 2  # True
        print '%s:' % utils.color('green', dsets[0]),
    else:
        raise Exception('one of \'em has to be length 1: %d %d' % (len(methods), len(dsets)))
    print '%s' % (' %s ' % utils.color('light_blue', 'vs')).join(glslabels)
    make_gls_tree_plot(args, outdir + '/' + '-vs-'.join(methods) + '/gls-gen-plots', study + '-' + '-vs-'.join(dsets),
                       glsfnames=[get_gls_fname(ddir, meth, locus=mfo['locus'], data=True) for ddir in data_outdirs for meth in methods],
                       glslabels=glslabels,
                       locus=mfo['locus'],
                       title=title,
                       title_color=title_color,
                       legends=legends,
                       legend_title=legend_title,
                       pie_chart_faces=pie_chart_faces)

# ----------------------------------------------------------------------------------------
def plot_single_test(args, baseoutdir, method):
    import plotting
    plot_types = ['missing', 'spurious']

    def get_performance(varname, varval):
        perf_vals = {pt : [] for pt in plot_types + ['total']}
        for iproc in range(args.iteststart, args.n_tests):
            single_vals = get_single_performance(get_outdir(args, baseoutdir, varname, varval, n_events=n_events) + '/' + str(iproc), method=method)
            for ptype in plot_types + ['total']:
                perf_vals[ptype].append(single_vals[ptype])
        return perf_vals

    plotvals = []
    for varval in args.varvals:
        print '%s %s' % (args.action, varvalstr(args.action, varval))
        plotvals.append({pt : {k : [] for k in ['xvals', 'ycounts', 'ytotals']} for pt in plot_types})
        for n_events in args.n_event_list:
            perf_vals = get_performance(varname=args.action, varval=varval)
            print '  %d' % n_events
            print '    iproc    %s' % ' '.join([str(i) for i in range(args.iteststart, args.n_tests)])
            print '    missing  %s' % ' '.join([str(v) for v in perf_vals['missing']]).replace('0', ' ')
            print '    spurious %s' % ' '.join([str(v) for v in perf_vals['spurious']]).replace('0', ' ')
            for ptype in plot_types:
                count = sum(perf_vals[ptype])
                plotvals[-1][ptype]['xvals'].append(n_events)
                plotvals[-1][ptype]['ycounts'].append(count)
                plotvals[-1][ptype]['ytotals'].append(sum(perf_vals['total']))
    for ptype in plot_types:
        print '    %s' % baseoutdir + '/' + ptype + '.svg'
        plotting.plot_gl_inference_fractions(baseoutdir, ptype,
                                             [pv[ptype] for pv in plotvals],
                                             labels=[legend_str(args, v) for v in args.varvals],
                                             xlabel='sample size',
                                             ylabel='%s alleles / total' % ptype,
                                             leg_title=varval_legend_titles.get(args.action, None),
                                             title=varval_titles.get(args.action, None)
        )

# ----------------------------------------------------------------------------------------
def write_single_zenodo_subdir(zenodo_dir, args, study, dset, method, mfo):
    method_outdir = heads.get_datadir(study, 'processed', extra_str='gls-gen-paper-' + args.label) + '/' + dset
    gls_dir = get_gls_dir(method_outdir, method, data=True)
    print '            %s --> %s' % (gls_dir, zenodo_dir)
    glfo = glutils.read_glfo(gls_dir, mfo['locus'], remove_orfs='partis' in method)
    glutils.write_glfo(zenodo_dir, glfo)
    if method == 'partis':
        # allele finding plots
        plotdir = gls_dir.replace('hmm/germline-sets', 'plots/sw/allele-finding')
        if not os.path.exists(zenodo_dir + '/fits'):
            os.makedirs(zenodo_dir + '/fits')
        for genedir in glob.glob(plotdir + '/try-0/*'):  # would be nice to copy html, but links will be wrong
            subprocess.check_call(['cp', '-r', genedir, zenodo_dir + '/fits/'])

        # csv prevalence files
        for region in utils.regions:
            with open(gls_dir.replace('/germline-sets', '/%s_gene-probs.csv' % region)) as infile:
                reader = csv.DictReader(infile)
                countfo = {line['%s_gene' % region] : int(line['count']) for line in reader}
                old_total = sum(countfo.values())
                orf_genes = [g for g in countfo if g not in glfo['seqs'][region]]  # this is kind of dangerous... but the genes are read from the same parameter dir that we're reading this prevalence file, so the only way it's gonna be missing is if we just removed it with the read_glfo() line above
                for ogene in orf_genes:
                    # if region == 'v':
                    #     _, nearest_gene, _ = glutils.find_nearest_gene_with_same_cpos(glfo, glfo['seqs'][region][ogene])  # oops, that's dumb... of course it isn't there
                    # else:
                    nearest_gene = glutils.find_nearest_gene_using_names(glfo, ogene)
                    # print '  adding %d to %s from %s' % (countfo[ogene], utils.color_gene(nearest_gene), utils.color_gene(ogene))
                    countfo[nearest_gene] += countfo[ogene]
                for ogene in orf_genes:
                    del countfo[ogene]
                assert old_total == sum(countfo.values())
                with open('%s/%s_gene-probs.csv' % (zenodo_dir, region), 'w') as outfile:
                    writer = csv.DictWriter(outfile, ('%s_gene' % region, 'count'))
                    writer.writeheader()
                    for gene in countfo:
                        writer.writerow({'%s_gene' % region : gene, 'count' : countfo[gene]})
    elif method == 'tigger-default':
        # doesn't seem to have written anything
        pass
    elif method == 'igdiscover':
        # for fname in ['errorhistograms.pdf', 'V_usage.pdf', 'V_usage.tab']:
        #     subprocess.check_call(['cp', '%s/work/final/%s' % (gls_dir, fname), zenodo_dir + '/'])
        subprocess.check_call(['cp', '-r', '%s/work/final' % gls_dir, zenodo_dir + '/'])  # aw, screw it, just write everything. The simulation stuff is already huge, anyway
    else:
        assert False

# ----------------------------------------------------------------------------------------
def write_zenodo_files(args, baseoutdir):
    for study, dset in [v.split('/') for v in args.varvals]:
        print '%-10s  %15s   %s' % (study, dset, baseoutdir)
        metafos = heads.read_metadata(study)
        for method in args.methods:
            outdir = get_outdir(args, baseoutdir, varname='data', varval='zenodo/%s/%s/%s' % (study_translations.get(study, study), dset, method.replace('-default', '')))
            print '  %-15s' % method
            write_single_zenodo_subdir(outdir, args, study, dset, method, metafos[dset])

# ----------------------------------------------------------------------------------------
def plot_tests(args, baseoutdir, method, method_vs_method=False, annotation_performance_plots=False, print_summary_table=False):
    if args.action == 'gls-gen':
        if annotation_performance_plots:
            assert method is None
            get_gls_gen_annotation_performance_plots(args, baseoutdir)
        elif print_summary_table:
            assert method is None
            print_gls_gen_summary_table(args, baseoutdir)
        else:
            get_gls_gen_tree_plots(args, baseoutdir, method)
    elif args.action == 'data':
        dsetfos = [v.split('/') for v in args.varvals]  # (study, dset)
        if method_vs_method:
            assert method is None
            for study, dset in dsetfos:
                get_data_plots(args, baseoutdir, args.methods, study, [dset])
        else:
            sample_groups = []
            for study in all_data_groups:
                for samples in all_data_groups[study]:
                    if len([s for s in samples if [study, s] in dsetfos]) == len(samples):  # if they're all in <dsetfos>
                        sample_groups.append((study, samples))
            for study, samples in sample_groups:
                get_data_plots(args, baseoutdir, [method], study, samples)
                for sample in samples:
                    dsetfos.remove([study, sample])
            for study, dset in dsetfos:
                print 'hmmmm %s (probably need to set --method-vs-method' % dset  # crashes below since both method and dset lists are of length one
                # get_data_plots(args, baseoutdir, [method], study, [dset])
    else:
        plot_single_test(args, baseoutdir, method)

# ----------------------------------------------------------------------------------------
def get_base_cmd(args, n_events, method):
    cmd = './bin/test-germline-inference.py'
    cmd += ' --n-procs ' + str(args.n_procs_per_test) + ' --n-tests ' + str(args.n_tests)
    if args.iteststart != 0:
        cmd += ' --iteststart ' + str(args.iteststart)
    cmd += ' --methods ' + method
    cmd += ' --n-sim-events ' + str(n_events)
    if not args.no_slurm:
        cmd += ' --slurm'
    if args.action != 'gls-gen':
        cmd += ' --inf-v-genes ' + args.v_genes[0]
    return cmd

# ----------------------------------------------------------------------------------------
def run_single_test(args, baseoutdir, val, n_events, method):
    cmd = get_base_cmd(args, n_events, method)
    outdir = get_outdir(args, baseoutdir, args.action, val, n_events=n_events)
    sim_v_genes = [args.v_genes[0]]
    mut_mult = None
    nsnpstr, nindelstr = '1', ''
    if args.action == 'mfreq':
        mut_mult = val
    elif args.action == 'nsnp':
        nsnpstr = str(val)
    elif args.action == 'multi-nsnp':
        nsnpstr = ':'.join([str(n) for n in val])
        sim_v_genes *= len(val)
    elif args.action == 'prevalence':
        cmd += ' --allele-prevalence-freqs ' + str(1. - val) + ':' + str(val)  # i.e. previously-known allele has 1 - p, and new allele has p
    elif args.action == 'n-leaves':
        cmd += ' --n-leaves ' + str(val)  # NOTE default of 1 (for other tests) is set in test-germline-inference.py
        cmd += ' --n-leaf-distribution geometric'
        cmd += ' --n-max-queries ' + str(n_events)  # i.e. we simulate <n_events> rearrangement events, but then only use <n_events> sequences for inference
    elif args.action == 'weibull':
        cmd += ' --n-leaves 5'  # NOTE default of 1 (for other tests) is set in test-germline-inference.py
        cmd += ' --n-leaf-distribution geometric'
        cmd += ' --root-mrca-weibull-parameter ' + str(val)
        cmd += ' --n-max-queries ' + str(n_events)  # i.e. we simulate <n_events> rearrangement events, but then only use <n_events> sequences for inference
    elif args.action == 'alcluster':
        nsnpstr = val['snp']
        nindelstr = val['indel']
        sim_v_genes *= len(val['snp'].split(':'))
    elif args.action == 'gls-gen':
        nsnpstr = '1:1:2:3:5'
        nindelstr = '' # '0:0:0:0:0:0:0:0'
        if args.gls_gen_difficulty == 'easy':
            # genes_per_region_str = '20:5:3'
            # n_sim_alleles_per_gene_str = '1,2:1,2:1,2'
            min_allele_prevalence_freq = 0.15
            mut_mult = 0.3
        elif args.gls_gen_difficulty == 'hard':
            # genes_per_region_str = '25:5:3'
            # n_sim_alleles_per_gene_str = '1,2,3:1,2:1,2'
            min_allele_prevalence_freq = 0.05
            mut_mult = 1.
        else:
            assert False
        # cmd += ' --n-genes-per-region ' + genes_per_region_str
        # cmd += ' --n-sim-alleles-per-gene ' + n_sim_alleles_per_gene_str
        cmd += ' --min-allele-prevalence-freq ' + str(min_allele_prevalence_freq)
        cmd += ' --gls-gen'
    else:
        assert False

    if mut_mult is not None:
        cmd += ' --mut-mult ' + str(mut_mult)

    if args.action != 'gls-gen':
        cmd += ' --sim-v-genes ' + ':'.join(sim_v_genes)
    if '--nosim' not in cmd:
        if nsnpstr != '':
            cmd += ' --nsnp-list ' + nsnpstr
        if nindelstr != '':
            cmd += ' --nindel-list ' + nindelstr
    if args.plot_annotation_performance:
        cmd += ' --plot-annotation-performance'
    cmd += ' --outdir ' + outdir
    utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_data(args, baseoutdir, study, dset, method):
    cmd = './datascripts/run.py cache-parameters'
    cmd += ' --study ' + study
    cmd += ' --samples ' + dset
    assert args.label is not None  # it's got a default now, so it shouldn't anymore be None
    cmd += ' --extra-str gls-gen-paper-' + args.label
    if args.no_slurm:
        cmd += ' --no-slurm'
    cmd += ' --n-procs ' + str(args.n_procs_per_test)
    if args.n_random_queries is not None:
        assert method == 'partis' or method == 'tigger-default' or method == 'igdiscover'  # I don't think it works for any others a.t.m.
        cmd += ' --n-random-queries ' + str(args.n_random_queries)
    if args.check:
        cmd += ' --check'
    cmd += ' --force-default-initial-germline-dir'
    if method != 'partis':
        cmd += ' --other-method ' + method
    if args.n_max_jobs is not None:
        cmd += ' --n-max-jobs ' + str(args.n_max_jobs)

    utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_tests(args, baseoutdir, method):
    if args.action == 'gls-gen':
        n_events = args.gls_gen_events
        val = 'simu'
        run_single_test(args, baseoutdir, val, n_events, method)
    elif args.action == 'data':
        for var in args.varvals:
            study, dset = var.split('/')
            run_data(args, baseoutdir, study, dset, method)
    else:
        for val in args.varvals:
            for n_events in args.n_event_list:
                run_single_test(args, baseoutdir, val, n_events, method)

# ----------------------------------------------------------------------------------------
default_varvals = {
    'mfreq' : '0.1:1.0:2.0',
    'nsnp' : '1:2:3:4',
    'multi-nsnp' : '1,1:1,3:2,3',
    'prevalence' : '0.1:0.2:0.3',
    'n-leaves' : '1.5:3:10',
    'weibull' : '0.3:0.5:1.3',
    'alcluster' : [
        # {'snp' : 25, 'indel' : 3},
        {'snp' : '75:100:125', 'indel' : '3:4:5'},
    ],
    'gls-gen' : None,
    'data' : {
        # 'jason-mg' : [
        #     'AR02-igh', 'AR03-igh', 'AR04-igh', 'AR05-igh', 'HD07-igh', 'HD09-igh', 'HD10-igh', 'HD13-igh', 'MK02-igh', 'MK03-igh', 'MK04-igh', 'MK05-igh', 'MK08-igh',
        # ],
        # 'sheng-gssp' : [
        #     'lp23810-m-pool',  'lp23810-g-pool', 'lp08248-m-pool', 'lp08248-g-pool',
        #     # 'lp23810-m-pool', 'lp08248-m-pool',  # igm only, for method-vs-method including igdiscover
        # ],
        # 'three-finger' : ['3ftx-1-igh'], #, 'pla2-1-igh'],
        # 'kate-qrs' : ['1g', '4g', '1k', '1l', '4k', '4l'],
        # 'laura-mb-2' : ['BF520-m-W1', 'BF520-m-M9', 'BF520-g-W1', 'BF520-g-M9'], #, 'BF520-k-W1', 'BF520-l-W1', 'BF520-k-M9', 'BF520-l-M9']
        # 'jason-influenza' : ['FV-igh-m2d', 'FV-igh-p3d', 'FV-igh-p7d'],
        # 'jason-influenza' : [
        #     # 'FV-igh-m8d', 'FV-igh-m2d', 'FV-igh-m1h', 'FV-igh-p1h', 'FV-igh-p1d', 'FV-igh-p3d', 'FV-igh-p7d', 'FV-igh-p14d', 'FV-igh-p21d', 'FV-igh-p28d',
        #     # 'GMC-igh-m8d', 'GMC-igh-m2d', 'GMC-igh-m1h', 'GMC-igh-p1h', 'GMC-igh-p1d', 'GMC-igh-p3d', 'GMC-igh-p7d', 'GMC-igh-p14d', 'GMC-igh-p21d', 'GMC-igh-p28d',
        #     # 'IB-igh-m8d', 'IB-igh-m2d', 'IB-igh-m1h', 'IB-igh-p1h', 'IB-igh-p1d', 'IB-igh-p3d', 'IB-igh-p7d', 'IB-igh-p14d', 'IB-igh-p21d', 'IB-igh-p28d',
        #     'FV-igh', 'GMC-igh', 'IB-igh',  # merged
        # ],
        # 'davide-gl-valid' : ['B10', 'B11', 'B12', 'B13', 'B14', 'B16', 'B17', 'B18', 'B19', 'B20', 'B21'],
        'crotty-fna' : ['A01',],
    }
}
all_data_groups = {
    'kate-qrs' : [
        ['1g', '4g'],
        ['1k', '4k'],
        ['1l', '4l'],
    ],
    'laura-mb-2' : [
        ['BF520-m-W1', 'BF520-m-M9'],
        ['BF520-g-W1', 'BF520-g-M9'],
        ['BF520-k-W1', 'BF520-k-M9'],
        ['BF520-l-W1', 'BF520-l-M9'],
    ],
    'sheng-gssp' : [
        ['lp23810-m-pool', 'lp23810-g-pool'],
        ['lp08248-m-pool', 'lp08248-g-pool'],
    ],
    'jason-influenza' : [
        ['FV-igh-m2d', 'FV-igh-p3d', 'FV-igh-p7d'],
        # ['FV-igh-m8d', 'FV-igh-m1h', 'FV-igh-p28d'],
        # ['FV-igh-p21d', 'FV-igh-p1h', 'FV-igh-p1d'],
        # # 'FV-igh-p14d'  # odd one out (also the one where tigger crashes)
        ['GMC-igh-m2d', 'GMC-igh-p3d', 'GMC-igh-p7d'],
        # ['GMC-igh-m8d', 'GMC-igh-m1h', 'GMC-igh-p28d'],
        # ['GMC-igh-p21d', 'GMC-igh-p1h', 'GMC-igh-p1d'],
        # # 'GMC-igh-p14d'  # odd one out (also the one where tigger crashes)
        ['IB-igh-m2d', 'IB-igh-p3d', 'IB-igh-p7d'],
        # ['IB-igh-m8d', 'IB-igh-m1h', 'IB-igh-p28d'],
        # ['IB-igh-p21d', 'IB-igh-p1h', 'IB-igh-p1d'],
        # # 'IB-igh-p14d'  # odd one out (also the one where tigger crashes)
    ],
    # 'davide-gl-valid' : [
    #     ['B12', 'B16'],
    # ],
    # for comparing different individuals:
    # 'jason-influenza' : [
    #     ['FV-igh', 'GMC-igh'],
    # ],
    # 'jason-mg' : [
    #     ['HD07-igh', 'HD09-igh'],
    #     ['HD10-igh', 'HD13-igh'],
    #     ['AR02-igh', 'AR03-igh'],
    #     ['AR04-igh', 'AR05-igh'],
    #     ['MK02-igh', 'MK03-igh'],
    #     ['MK04-igh', 'MK05-igh'],
    #     # , 'MK08-igh'],
    # ],

}
default_varvals['data'] = ':'.join([study + '/' + heads.full_dataset(heads.read_metadata(study), dset) for study in default_varvals['data'] for dset in default_varvals['data'][study]])
for study in all_data_groups:
    for idp in range(len(all_data_groups[study])):
        all_data_groups[study][idp] = [heads.full_dataset(heads.read_metadata(study), ds) for ds in all_data_groups[study][idp]]
# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['mfreq', 'nsnp', 'multi-nsnp', 'prevalence', 'n-leaves', 'weibull', 'alcluster', 'gls-gen', 'data'])
parser.add_argument('--methods', default='partis') # not using <choices> 'cause it's harder since it's a list
parser.add_argument('--method-vs-method', action='store_true')
parser.add_argument('--v-genes', default='IGHV4-39*01')
parser.add_argument('--varvals')
parser.add_argument('--n-event-list', default='1000:2000:4000:8000')  # NOTE modified later for multi-nsnp also NOTE not used for gen-gset
parser.add_argument('--gls-gen-events', type=int, default=50000)
parser.add_argument('--gls-gen-difficulty', default='easy', choices=['easy', 'hard'])
parser.add_argument('--n-random-queries', type=int)
parser.add_argument('--n-max-jobs', type=int)
parser.add_argument('--n-tests', type=int, default=3)
parser.add_argument('--iteststart', type=int, default=0)
parser.add_argument('--n-procs-per-test', type=int, default=5)
parser.add_argument('--plot', action='store_true')
parser.add_argument('--write-zenodo-files', action='store_true')
parser.add_argument('--plot-annotation-performance', action='store_true')
parser.add_argument('--print-table', action='store_true')
parser.add_argument('--no-slurm', action='store_true')
parser.add_argument('--plotcache', action='store_true')
parser.add_argument('--only-print', action='store_true')
parser.add_argument('--check', action='store_true')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--label', default='xxx')
parser.add_argument('--ete-path', default='/home/' + os.getenv('USER') + '/anaconda_ete/bin')
args = parser.parse_args()

args.methods = sorted(utils.get_arg_list(args.methods))
args.v_genes = utils.get_arg_list(args.v_genes)
args.n_event_list = utils.get_arg_list(args.n_event_list, intify=True)

# ----------------------------------------------------------------------------------------
alfdir = utils.fsdir() + '/partis/allele-finder'
baseoutdir = alfdir
if args.label is not None:
    baseoutdir += '/' + args.label
baseoutdir += '/' + args.action

if args.varvals is None:
    args.varvals = default_varvals[args.action]
kwargs = {}
if args.action == 'mfreq' or args.action == 'prevalence' or args.action == 'n-leaves' or args.action == 'weibull':
    kwargs['floatify'] = True
if args.action == 'nsnp':
    kwargs['intify'] = True
if args.action != 'alcluster':  # could also do this for data i think, if i remove that line up there ^
    args.varvals = utils.get_arg_list(args.varvals, **kwargs)
if args.action == 'multi-nsnp':
    args.varvals = [[int(n) for n in gstr.split(',')] for gstr in args.varvals]  # list of nsnps for each test, e.g. '1,1:2,2' runs two tests: 1) two new alleles, each with one snp and 2) two new alleles each with 2 snps
    factor = numpy.median([(len(nl) + 1) / 2. for nl in args.varvals])  # i.e. the ratio of (how many alleles we'll be dividing the events among), to (how many we'd be dividing them among for the other [single-nsnp] tests)
    args.n_event_list = [int(factor * n) for n in args.n_event_list]

if args.write_zenodo_files:
    assert args.action == 'data'  # would need to implement it
    write_zenodo_files(args, baseoutdir)
elif args.plot:
    if args.method_vs_method:
        plot_tests(args, baseoutdir, method=None, method_vs_method=True)
    elif args.plot_annotation_performance:
        plot_tests(args, baseoutdir, method=None, annotation_performance_plots=True)
    elif args.print_table:
        plot_tests(args, baseoutdir, method=None, print_summary_table=True)
    else:
        for method in [m for m in args.methods if m != 'simu']:
            plot_tests(args, baseoutdir, method)
else:
    for method in args.methods:
        run_tests(args, baseoutdir, method)
