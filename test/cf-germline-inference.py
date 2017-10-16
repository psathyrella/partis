#!/usr/bin/env python
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

legend_titles = {
    'mfreq' : 'mutation',
    'nsnp' : 'N SNPs',
    'multi-nsnp' : 'N SNPs',
    'prevalence' : 'prevalence',
    'n-leaves' : 'mean N leaves',
}

all_methods = ['tigger-default', 'igdiscover', 'partis']
characters = ['subject', 'isotype', 'timepoint']
def methstr(meth):
    if meth == 'tigger-default':
        return 'tigger'
    elif meth == 'full':
        return 'full IMGT'
    else:
        return meth
def diffstr(difficulty):
    if difficulty == 'easy':
        return 'simple'
    elif difficulty == 'hard':
        return 'complex'
    else:
        assert False

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
    return glutils.get_fname(outdir, locus, region)

# ----------------------------------------------------------------------------------------
def make_gls_tree_plot(args, plotdir, plotname, glsfnames, glslabels, locus, ref_label=None, leaf_names=False, title=None, title_color=None, legends=None, legend_title=None, pie_chart_faces=False):
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
    if leaf_names:
        cmdstr += ' --leaf-names'
    cmdstr += ' --locus ' + locus
    if args.plotcache:
        cmdstr += ' --use-cache'
    if args.only_print:
        cmdstr += ' --only-print'
    utils.simplerun(cmdstr, shell=True, debug=args.dry_run, dryrun=args.dry_run)

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
    varname = args.action
    varval = 'simu'
    plotnames = ['v_hamming_to_true_naive', 'v_muted_bases']
    meanvals = {pn : {methstr(m) : [] for m in args.methods} for pn in plotnames}
    for iproc in range(args.iteststart, args.n_tests):
        outdir = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/' + str(iproc)  # duplicates code in bin/test-germline-inference.py
        plotdir = outdir + '/annotation-performance-plots'
        print '    %s' % plotdir
        utils.prep_dir(plotdir, wildlings=['*.png', '*.svg', '*.csv'])
        for plotname in plotnames:
            print '      %s/%s.svg' % (plotdir, plotname)
            hists = [Hist(fname=get_gls_fname(outdir, meth, sim_locus, annotation_performance_plots=True) + '/' + plotname + '.csv', title=methstr(meth)) for meth in args.methods]
            for hist in hists:
                if hist.overflow_contents() != 0.0:
                    print '  %s %s non-zero under/overflow %f' % (utils.color('red', 'error'), hist.title, hist.overflow_contents())
                meanvals[plotname][hist.title].append(hist.get_mean())
                # print '%10s  %6.3f' % (hist.title, hist.get_mean())
            colors = [methcolors[meth] for meth in args.methods]
            plotting.draw_no_root(hists[0], log='y', plotdir=plotdir, plotname=plotname, more_hists=hists[1:], colors=colors, ytitle='sequences')

    for plotname in plotnames:
        if 'muted_bases' in plotname:  #  mean value isn't meaningful
            continue
        print plotname
        for method in args.methods:
            methtitle = methstr(method)
            mean = float(sum(meanvals[plotname][methtitle])) / len(meanvals[plotname][methtitle])
            print '   %15s  %6.3f / %d = %6.3f' % (methstr(methtitle), sum(meanvals[plotname][methtitle]), len(meanvals[plotname][methtitle]), mean)

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
                           locus=sim_locus, ref_label='sim', legend_title=methstr(method), title=diffstr(args.gls_gen_difficulty))  #, title_color=method)

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
        pie_chart_faces = True
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
        plotting.plot_gl_inference_fractions(baseoutdir, ptype, [pv[ptype] for pv in plotvals], labels=[legend_str(args, v) for v in args.varvals], xlabel='sample size', ylabel='fraction %s' % ptype, leg_title=legend_titles.get(args.action, None), title=ptype + ' alleles')

# ----------------------------------------------------------------------------------------
def plot_tests(args, baseoutdir, method, method_vs_method=False, annotation_performance_plots=False):
    if args.action == 'gls-gen':
        if annotation_performance_plots:
            assert method is None
            get_gls_gen_annotation_performance_plots(args, baseoutdir)
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
    if args.plot_performance:
        cmd += ' --plot-performance'
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
        # #     'HD07-igk', 'HD07-igl', 'AR03-igk', 'AR03-igl',
        # ],
        # 'sheng-gssp' : [
        #     'lp23810-m-pool',  'lp23810-g-pool', 'lp08248-m-pool', 'lp08248-g-pool',
        #     # 'lp23810-m-pool', 'lp08248-m-pool',
        # ],
        # 'three-finger' : ['3ftx-1-igh'], #, 'pla2-1-igh'],
        # 'kate-qrs' : ['1g', '4g', '1k', '1l', '4k', '4l'],
        # 'laura-mb-2' : ['BF520-m-W1', 'BF520-m-M9', 'BF520-g-W1', 'BF520-g-M9'], #, 'BF520-k-W1', 'BF520-l-W1', 'BF520-k-M9', 'BF520-l-M9']
        # 'jason-influenza' : ['GMC-igh-m8d', 'GMC-igh-m1h', 'GMC-igh-p28d'],
        # 'jason-influenza' : [
        #     # 'FV-igh-m8d', 'FV-igh-m2d', 'FV-igh-m1h', 'FV-igh-p1h', 'FV-igh-p1d', 'FV-igh-p3d', 'FV-igh-p7d', 'FV-igh-p14d', 'FV-igh-p21d', 'FV-igh-p28d',
        #     # 'GMC-igh-m8d', 'GMC-igh-m2d', 'GMC-igh-m1h', 'GMC-igh-p1h', 'GMC-igh-p1d', 'GMC-igh-p3d', 'GMC-igh-p7d', 'GMC-igh-p14d', 'GMC-igh-p21d', 'GMC-igh-p28d',
        #     # 'IB-igh-m8d', 'IB-igh-m2d', 'IB-igh-m1h', 'IB-igh-p1h', 'IB-igh-p1d', 'IB-igh-p3d', 'IB-igh-p7d', 'IB-igh-p14d', 'IB-igh-p21d', 'IB-igh-p28d',
        #     'FV-igh', 'GMC-igh', 'IB-igh',  # merged
        # ],
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
parser.add_argument('--plot-performance', action='store_true')
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

if args.plot:
    if args.method_vs_method:
        plot_tests(args, baseoutdir, method=None, method_vs_method=True)
    elif args.plot_performance:
        plot_tests(args, baseoutdir, method=None, annotation_performance_plots=True)
    else:
        for method in [m for m in args.methods if m != 'simu']:
            plot_tests(args, baseoutdir, method)
else:
    for method in args.methods:
        run_tests(args, baseoutdir, method)
