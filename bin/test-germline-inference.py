#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy
import copy
import random
import argparse
import time
import sys
import os
import glob
import colored_traceback.always

from subprocess import check_call

import partis.utils as utils
import partis.glutils as glutils
import partis.processargs as processargs

# ----------------------------------------------------------------------------------------
def cov_cmd():
    return 'coverage3 run --append'

# ----------------------------------------------------------------------------------------
def get_outfname(args, method, annotation_performance_plots=False, return_parent_gl_dir=False):
    outdir = args.outdir + '/' + method
    if not annotation_performance_plots:  # default: output is igh/ighv.fasta
        if method == 'partis' or method == 'full':  # parameter directory, not regular file (although, could change it to the gls .fa in sw/)
            outdir += '/sw/germline-sets'
        if not return_parent_gl_dir:
            return glutils.get_fname(outdir, args.locus, 'v')
        else:
            return outdir
    else:  # product of running partis annotation with --plot-annotation-performance
        return outdir + '/annotation-performance-plots'

# ----------------------------------------------------------------------------------------
def simulate(args):
    if utils.output_exists(args, args.simfname):
        return
    cmd_str = args.partis_path + ' simulate --n-sim-events ' + str(args.n_sim_events) + ' --outfname ' + args.simfname + ' --n-leaves ' + str(args.n_leaves) + ' --rearrange-from-scratch --force-dont-generate-germline-set --shm-parameter-dir ' + utils.get_partis_dir() + '/data/recombinator/scratch-parameters'
    cmd_str += ' --allow-nonfunctional-scratch-seqs'
    if args.n_leaf_distribution is None:
        cmd_str += ' --constant-number-of-leaves'
    else:
        cmd_str += ' --n-leaf-distribution ' + args.n_leaf_distribution
    if args.mut_mult is not None:
        cmd_str += ' --mutation-multiplier ' + str(args.mut_mult)
    if args.root_mrca_weibull_parameter is not None:
        cmd_str += ' --root-mrca-weibull-parameter ' + str(args.root_mrca_weibull_parameter)

    if args.n_procs is not None:
        cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd_str += ' --batch-system slurm'

    args.allele_prevalence_fname = args.workdir + '/allele-prevalence-freqs.csv'

    # figure what genes we're using
    if args.gls_gen:
        sglfo = glutils.read_glfo(args.default_germline_dir, locus=args.locus)
        glutils.remove_v_genes_with_bad_cysteines(sglfo)
        glutils.generate_germline_set(sglfo, args, new_allele_info=args.new_allele_info, dont_remove_template_genes=args.dont_remove_template_genes, debug=True)
        cmd_str += ' --allele-prevalence-fname ' + args.allele_prevalence_fname
    else:
        sglfo = glutils.read_glfo(args.default_germline_dir, locus=args.locus, only_genes=(args.sim_v_genes + args.dj_genes))
        added_snp_names = glutils.generate_new_alleles(sglfo, args.new_allele_info, debug=True, remove_template_genes=(not args.dont_remove_template_genes))  # NOTE template gene removal is the default for glutils.generate_germline_set

        if args.allele_prevalence_freqs is not None:
            if not utils.is_normed(args.allele_prevalence_freqs):
                raise Exception('--allele-prevalence-freqs %s not normalized' % args.allele_prevalence_freqs)
            if len(args.allele_prevalence_freqs) != len(sglfo['seqs']['v']):  # already checked when parsing args, but, you know...
                raise Exception('--allele-prevalence-freqs %d not the same length as sglfo %d' % (len(args.allele_prevalence_freqs), len(sglfo['seqs']['v'])))
            gene_list = sorted(sglfo['seqs']['v']) if len(added_snp_names) == 0 else list(set(args.sim_v_genes)) + added_snp_names
            prevalence_freqs = {'v' : {g : f for g, f in zip(gene_list, args.allele_prevalence_freqs)}, 'd' : {}, 'j' : {}}
            glutils.write_allele_prevalence_freqs(prevalence_freqs, args.allele_prevalence_fname)
            cmd_str += ' --allele-prevalence-fname ' + args.allele_prevalence_fname

    glutils.write_glfo(args.outdir + '/germlines/simulation', sglfo)
    cmd_str += ' --initial-germline-dir ' + args.outdir + '/germlines/simulation'
    # glutils.print_glfo(sglfo)

    # run simulation
    if args.seed is not None:
        cmd_str += ' --random-seed ' + str(args.seed)
    utils.simplerun(cmd_str, dryrun=args.dryrun)

# ----------------------------------------------------------------------------------------
def run_other_method(args, method):
    if method not in ['tigger-default', 'tigger-tuned', 'igdiscover']:  # really just to make it easier to search for this fcn
        assert False
    assert args.n_max_queries is None
    if utils.output_exists(args, get_outfname(args, method)):
        return
    simfasta = utils.getprefix(args.simfname) + '.fa'
    utils.csv_to_fasta(args.simfname, outfname=simfasta, overwrite=False, remove_duplicates=True)
    cmd = './test/%s-run.py' % method.split('-')[0]
    if method == 'tigger-tuned':
        cmd += ' --tuned-tigger-params'
    cmd += ' --infname ' + simfasta
    cmd += ' --outfname ' + get_outfname(args, method)
    if args.species != 'human':
        cmd += ' --species %s' % args.species
    if args.overwrite:
        cmd += ' --overwrite'
    if args.gls_gen:
        cmd += ' --gls-gen'
        cmd += ' --glfo-dir ' + utils.get_partis_dir() + '/' + args.default_germline_dir  # the partis mehods have this as the default internally, but we want/have to set it explicitly here
    else:
        cmd += ' --glfo-dir ' + args.inf_glfo_dir
    cmd += ' --simulation-germline-dir ' + args.outdir + '/germlines/simulation'  # alleleclusterer is the only one that really uses this, but for now I want its dbg output to have the sim info
    if method != 'igdiscover':  # for now we're saving all the igdiscover output/intermediate files, so we write them to an output dir
        cmd += ' --workdir ' + args.workdir + '/' + method
    if args.n_procs is not None:
        cmd += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd += ' --slurm'

    utils.simplerun(cmd, dryrun=args.dryrun)

# ----------------------------------------------------------------------------------------
def run_performance_plot(args, method):
    perf_outdir = get_outfname(args, method, annotation_performance_plots=True)
    if utils.output_exists(args, perf_outdir):
        return

    cmd_str = args.partis_path + ' cache-parameters --infname ' + args.simfname + ' --plot-annotation-performance'
    cmd_str += ' --is-simu --simulation-germline-dir ' + args.outdir + '/germlines/simulation'
    cmd_str += ' --initial-germline-dir ' + get_outfname(args, method, return_parent_gl_dir=True)  # i.e. use the inferred glfo from <method>
    cmd_str += ' --parameter-dir ' + perf_outdir + '/dummy-parameter-dir'
    cmd_str += ' --plotdir ' + perf_outdir
    cmd_str += ' --only-smith-waterman --leave-default-germline --dont-write-parameters'  # i.e. we really want to annotate, not cache parameters, but then it'd look for a parameter dir
    if args.n_procs is not None:
        cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.n_max_queries is not None:
        cmd_str += ' --n-max-queries ' + str(args.n_max_queries)  # NOTE do *not* use --n-random-queries, since it'll change the cluster size distribution
    if args.slurm:
        cmd_str += ' --batch-system slurm'
    if args.seed is not None:
        cmd_str += ' --random-seed ' + str(args.seed)
    utils.simplerun(cmd_str, dryrun=args.dryrun)

# ----------------------------------------------------------------------------------------
def run_partis_parameter_cache(args, method):
    if utils.output_exists(args, get_outfname(args, method)):
        return

    paramdir = args.outdir + '/' + method
    plotdir = args.outdir + '/' + method + '/plots'

    # remove any old sw cache files
    sw_cachefiles = glob.glob(paramdir + '/sw-cache-*.csv')
    if len(sw_cachefiles) > 0:
        for cachefname in sw_cachefiles:
            check_call(['rm', '-v', cachefname])
            sw_cache_gldir = cachefname.replace('.csv', '-glfo')
            if os.path.exists(sw_cache_gldir):  # if stuff fails halfway through, you can get one but not the other
                glutils.remove_glfo_files(sw_cache_gldir, args.locus)
                # os.rmdir(sw_cache_gldir)

    # generate germline set and cache parameters
    cmd_str = args.partis_path + ' cache-parameters --infname ' + args.simfname + ' --only-smith-waterman'
    cmd_str += ' --initial-germline-dir %s' % (args.default_germline_dir if args.gls_gen else args.inf_glfo_dir)
    if method == 'partis':
        cmd_str += ' --debug-allele-finding' # --always-find-new-alleles'
        cmd_str += ' --is-simu --simulation-germline-dir ' + args.outdir + '/germlines/simulation'  # alleleclusterer is the only one that really uses this, but for now I want its dbg output to have the sim info
        if args.allele_cluster:
            cmd_str += ' --allele-cluster'
            if args.kmeans_allele_cluster:
                cmd_str += ' --kmeans-allele-cluster'
    elif method == 'full':
        cmd_str += ' --leave-default-germline'
    else:
        assert False

    if args.species != 'human':
        cmd_str += ' --species %s' % args.species

    if args.n_procs is not None:
        cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.n_max_queries is not None:
        cmd_str += ' --n-max-queries ' + str(args.n_max_queries)  # NOTE do *not* use --n-random-queries, since it'll change the cluster size distribution
    if args.slurm:
        cmd_str += ' --batch-system slurm'

    cmd_str += ' --parameter-dir ' + paramdir
    cmd_str += ' --plotdir ' + plotdir
    if args.seed is not None:
        cmd_str += ' --random-seed ' + str(args.seed)
    if args.plot_and_fit_absolutely_everything is not None:
        cmd_str += ' --plot-and-fit-absolutely-everything ' + str(args.plot_and_fit_absolutely_everything)
    utils.simplerun(cmd_str, dryrun=args.dryrun)

# ----------------------------------------------------------------------------------------
def write_inf_glfo(args):  # read default glfo, restrict it to the specified alleles, and write to somewhere where all the methods can read it
    # NOTE this dir should *not* be modified by any of the methods
    inf_glfo = glutils.read_glfo(args.default_germline_dir, locus=args.locus, only_genes=args.inf_v_genes + args.dj_genes)
    print('  writing initial inference glfo with %d v: %s' % (len(inf_glfo['seqs']['v']), ' '.join([utils.color_gene(g) for g in inf_glfo['seqs']['v']])))
    glutils.write_glfo(args.inf_glfo_dir, inf_glfo)

# ----------------------------------------------------------------------------------------
def run_tests(args):
    print('seed %d' % args.seed)
    # all fcns return immediately if output already exists

    if 'simu' in args.methods:
        simulate(args)
        args.methods.remove('simu')

    if not args.gls_gen:
        write_inf_glfo(args)
    for method in args.methods:
        if args.plot_annotation_performance:
            run_performance_plot(args, method)
        elif method == 'partis' or method == 'full':
            run_partis_parameter_cache(args, method)
        else:
            run_other_method(args, method)

# ----------------------------------------------------------------------------------------
def multiple_tests(args):
    def getlogdir(iproc):
        logdir = args.outdir + '/' + str(iproc) + '/logs'
        if args.plot_annotation_performance:
            logdir += '/annotation-performance-plots'
        return logdir + '/' + '-'.join(args.methods)
    def cmd_str(iproc):
        clist = copy.deepcopy(sys.argv)
        utils.remove_from_arglist(clist, '--n-tests', has_arg=True)
        utils.remove_from_arglist(clist, '--iteststart', has_arg=True)
        utils.replace_in_arglist(clist, '--outdir', args.outdir + '/' + str(iproc))
        utils.replace_in_arglist(clist, '--seed', str(args.seed + iproc))
        # clist.append('--slurm')
        return ' '.join(clist)

    for iproc in range(args.iteststart, args.n_tests):  # don't overwrite old log files... need to eventually fix this so it isn't necessary
        def lfn(iproc, ilog):
            logfname =  args.outdir + '/' + str(iproc) + '/log'
            if ilog > 0:
                logfname += '.' + str(ilog)
            return logfname

    cmdfos = [{'cmd_str' : cmd_str(iproc),
               'workdir' : args.workdir + '/' + str(iproc),
               'logdir' : getlogdir(iproc),
               'outfname' : args.outdir + '/' + str(iproc)}
              for iproc in range(args.iteststart, args.n_tests)]
    if args.dryrun:
        for iproc in range(args.iteststart, args.n_tests):
            utils.simplerun(cmdfos[iproc - args.iteststart]['cmd_str'], dryrun=True)
        return
    for iproc in range(args.iteststart, args.n_tests):
        logd = getlogdir(iproc)
        if os.path.exists(logd + '/log'):
            ilog = 0
            while os.path.exists(logd + '/log.' + str(ilog)):
                ilog += 1
            check_call(['mv', '-v', logd + '/log', logd + '/log.' + str(ilog)])
    print('  look for logs in %s' % args.outdir)
    utils.run_cmds(cmdfos, debug='write')

# ----------------------------------------------------------------------------------------

# # ----------------------------------------------------------------------------------------
# from hist import Hist
# import plotting
# fig, ax = plotting.mpl_init()

# ntrees = 1000
# distrs = [
#     # (1.5, 'geo'),
#     # (3, 'geo'),
#     (10, 'geo'),
#     # (25, 'geo'),
#     # (2.3, 'zipf'),
#     # (1.8, 'zipf'),
#     # (1.3, 'zipf'),
# ]

# # ----------------------------------------------------------------------------------------
# def getsubsample(vals):
#     print vals
#     iclust = 0
#     seqs = []
#     for v in vals:
#         seqs += [iclust for _ in range(v)]
#         iclust += 1
#     print seqs
#     subseqs = numpy.random.choice(seqs, size=ntrees, replace=False)
#     # subseqs = seqs[:ntrees]
#     print subseqs
#     import itertools
#     subvals = []
#     for _, group in itertools.groupby(sorted(subseqs)):
#         for what in set(group):
#             subg = [s for s in subseqs if s == what]
#             print what, len(subg)
#             subvals.append(len(subg))
#     print subvals
#     return subvals

# ih = 0
# for n_leaves, fcn in distrs:
#     if fcn == 'zipf':
#         vals = numpy.random.zipf(n_leaves, size=ntrees)  # NOTE <n_leaves> is not the mean here
#     elif fcn == 'geo':
#         vals = numpy.random.geometric(1. / n_leaves, size=ntrees)
#     else:
#         assert False
#     nbins = 100
#     htmp = Hist(nbins, -0.5, nbins - 0.5)
#     for v in vals:
#         htmp.fill(v)
#     htmp.mpl_plot(ax, color=plotting.default_colors[ih], errors=False, label='%s %.1f' % (fcn, numpy.mean(vals)))
# # ----------------------------------------------------------------------------------------
#     hsub = Hist(nbins, -0.5, nbins - 0.5)
#     subvals = getsubsample(vals)
#     for v in subvals:
#         hsub.fill(v)
#     hsub.mpl_plot(ax, color=plotting.default_colors[ih], errors=False, label='%s %.1f' % (fcn, numpy.mean(subvals)), linestyle='--')
# # ----------------------------------------------------------------------------------------
#     ih += 1

# plotting.mpl_finish(ax, utils.fsdir() + '/partis/tmp/tmp', 'baz', xbounds=(0.9, nbins), log='y')
# sys.exit()
# # ----------------------------------------------------------------------------------------

example_str = '\n    '.join(['example usage:',
                             'one new allele separated by 3 snps from existing allele:',
                             '    ./bin/test-germline-inference.py --n-sim-events 2000 --n-procs 10 --sim-v-genes=IGHV1-18*01 --inf-v-genes=IGHV1-18*01 --snp-positions 27,55,88',
                             'one new allele [i.e. that the inference doesn\'t know about, but that in this case is in IMGT] separated by 1 snp from existing allele:',
                             '    ./bin/test-germline-inference.py --n-sim-events 2000 --n-procs 10 --sim-v-genes=IGHV4-39*01:IGHV4-39*02 --inf-v-genes=IGHV4-39*01',
                             'generate a full germline set for simulation, and then try to infer it:',
                             '    ./bin/test-germline-inference.py --n-sim-events 2000 --n-procs 10 --gls-gen'])
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, epilog=example_str)
parser.add_argument('--n-sim-events', type=int, default=20, help='number of simulated rearrangement events')
parser.add_argument('--n-max-queries', type=int, help='number of queries to use for inference from the simulation sample')
parser.add_argument('--n-leaves', type=float, default=1., help='see bin/partis --help')
parser.add_argument('--n-leaf-distribution', help='see bin/partis --help')
parser.add_argument('--root-mrca-weibull-parameter', type=float, help='see bin/partis --help')
parser.add_argument('--n-procs', type=int)
parser.add_argument('--seed', type=int, default=int(time.time()), help='random seed')
parser.add_argument('--gls-gen', action='store_true', help='generate a random germline set from scratch (parameters specified above), and infer a germline set from scratch, instead of using --sim-v-genes, --dj-genes, --inf-v-genes.')
parser.add_argument('--sim-v-genes', default='IGHV4-39*01:IGHV4-39*08', help='V genes to use for simulation')
parser.add_argument('--inf-v-genes', default='IGHV4-39*01', help='V genes to use for inference')
parser.add_argument('--dj-genes', default='IGHD6-19*01:IGHJ4*02', help='D and J genes to use for both simulation and inference')
parser.add_argument('--snp-positions', help='colon-separated list (length must equal length of <--sim-v-genes>) of comma-separated snp positions for each gene, e.g. for two genes you might have \'3,71:45\'')
parser.add_argument('--nsnp-list', help='colon-separated list (length must equal length of <--sim-v-genes> unless --gls-gen) of the number of snps to generate for each gene (each snp at a random position). If --gls-gen, then this still gives the number of snpd genes, but it isn\'t assumed to be the same length as anything [i.e. we don\'t yet know how many v genes there\'ll be]')
parser.add_argument('--indel-positions', help='see --snp-positions (a.t.m. the indel length distributions are hardcoded)')
parser.add_argument('--nindel-list', help='see --nsnp-list')
parser.add_argument('--n-genes-per-region', default='::', help='see bin/partis --help')
parser.add_argument('--n-sim-alleles-per-gene', default='::', help='see bin/partis --help')
parser.add_argument('--min-sim-allele-prevalence-freq', default=glutils.default_min_allele_prevalence_freq, type=float, help='see bin/partis --help')
parser.add_argument('--allele-prevalence-freqs', help='colon-separated list of allele prevalence frequencies, including newly-generated snpd genes (ordered alphabetically)')
parser.add_argument('--dont-remove-template-genes', action='store_true', help='when generating snps, *don\'t* remove the original gene before simulation')  # NOTE template gene removal is the default for glutils.generate_germline_set
parser.add_argument('--mut-mult', type=float, help='DO NOT USE use --mutation-multiplier (see below)')
parser.add_argument('--mutation-multiplier', type=float, help='see bin/partis --help')  # see note below
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dryrun', action='store_true')
parser.add_argument('--allele-cluster', action='store_true', help='see bin/partis --help')
parser.add_argument('--kmeans-allele-cluster', action='store_true', help='see bin/partis --help')
parser.add_argument('--plot-annotation-performance', action='store_true', help='see bin/partis --help')
parser.add_argument('--methods', default='simu:partis', help='colon-separated list of methods to run. By default runs simulation, and then partis inference (igdiscover and tigger, if installed, are the other options)')
parser.add_argument('--outdir', default=utils.fsdir() + '/partis/allele-finder')
parser.add_argument('--inf-glfo-dir', help='default set below')
parser.add_argument('--simfname', help='default set below')
parser.add_argument('--workdir', default=utils.fsdir() + '/_tmp/hmms/' + str(random.randint(0, 999999)))
parser.add_argument('--n-tests', type=int, help='instead of just running once, run <N> independent tests simultaneously')
parser.add_argument('--iteststart', type=int, default=0, help='for use with --n-tests, if you want to add more tests on')
parser.add_argument('--plot-and-fit-absolutely-everything', type=int, help='fit every single position for this <istart> and write every single corresponding plot (slow as hell, and only for debugging/making plots for paper)')
parser.add_argument('--partis-path', default='./bin/partis')
parser.add_argument('--prepend-coverage-command', action='store_true', help='see bin/partis --help')
parser.add_argument('--species', default='human', choices=('human', 'macaque'))
parser.add_argument('--locus', default='igh')
parser.add_argument('--allele-prevalence-fname', help='for internal use only (set above)')

args = parser.parse_args()
assert args.locus == 'igh'  # would just need to update some things, e.g. propagate through to the various methods
args.methods = utils.get_arg_list(args.methods)
available_methods = set(['simu', 'partis', 'full', 'tigger-default', 'tigger-tuned', 'igdiscover'])
if len(set(args.methods) - available_methods) > 0:
    raise Exception('unexpected --methods: %s' % ' '.join(set(args.methods) - available_methods))
# args.default_germline_dir = 'old-glfo/%s' % args.species  # 'data/germlines/%s' % args.species  # NOTE gad damnit, I just deleted old-glfo, had no idea what it was for
print('  %s hopefully old-glfo/ isn\'t needed to recreate old results (see comment)' % utils.color('yellow', 'note:'))
args.default_germline_dir = 'data/germlines/%s' % args.species  # 'data/germlines/%s' % args.species

args.generate_germline_set = args.gls_gen  # for compatibility with bin/partis (i.e. so they can both use the fcn in processargs, but I don't have to rewrite either)
args.mut_mult = args.mutation_multiplier  # for compatibility with bin/partis (i.e. so they can both use the fcn in processargs, but I don't have to rewrite either)
if args.generate_germline_set:  # if we're generating/inferring a whole germline set these are either set automatically or not used
    delattr(args, 'sim_v_genes')
    delattr(args, 'inf_v_genes')
    delattr(args, 'dj_genes')
    args.allele_prevalence_freqs = None
    args.inf_glfo_dir = None
else:
    args.dj_genes = utils.get_arg_list(args.dj_genes)
    args.sim_v_genes = utils.get_arg_list(args.sim_v_genes)
    args.inf_v_genes = utils.get_arg_list(args.inf_v_genes)
    args.allele_prevalence_freqs = utils.get_arg_list(args.allele_prevalence_freqs, floatify=True)

processargs.process_gls_gen_args(args)  # well, also does stuff with non-gls-gen new allele args

if args.inf_glfo_dir is None:
    args.inf_glfo_dir = args.outdir + '/germlines/inference'
if args.simfname is None:
    args.simfname = args.outdir + '/simu.yaml'

if args.prepend_coverage_command:
    args.partis_path = '%s %s' % (cov_cmd(), args.partis_path)

if args.seed is not None:
    random.seed(args.seed)
    numpy.random.seed(args.seed)

if args.n_tests is not None:
    multiple_tests(args)
else:
    run_tests(args)
