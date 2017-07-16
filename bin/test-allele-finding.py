#!/usr/bin/env python
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
sys.path.insert(1, './python')
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')

import utils
import glutils

# ----------------------------------------------------------------------------------------
def get_outfname(args, method):
    outdir = args.outdir + '/' + method
    if method == 'partis' or method == 'full':  # parameter directory, not regular file (although, could change it to the gls .fa in sw/)
        outdir += '/sw/germline-sets'
    return glutils.get_fname(outdir, args.locus, 'v')

# ----------------------------------------------------------------------------------------
def simulate(args):
    if utils.output_exists(args, args.simfname):
        return
    cmd_str = args.partis_path + ' simulate --n-sim-events ' + str(args.n_sim_events) + ' --n-leaves ' + str(args.n_leaves) + ' --rearrange-from-scratch --outfname ' + args.simfname
    if args.n_leaf_distribution is None:
        cmd_str += ' --constant-number-of-leaves'
    else:
        cmd_str += ' --n-leaf-distribution ' + args.n_leaf_distribution
    if args.mut_mult is not None:
        cmd_str += ' --mutation-multiplier ' + str(args.mut_mult)

    cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd_str += ' --batch-system slurm --subsimproc'

    allele_prevalence_fname = args.workdir + '/allele-prevalence-freqs.csv'

    # figure what genes we're using
    if args.gls_gen:
        assert args.sim_v_genes is None and args.allele_prevalence_freqs is None

        remove_template_genes = False

        sglfo = glutils.read_glfo('data/germlines/human', locus=args.locus)
        glutils.remove_v_genes_with_bad_cysteines(sglfo)
        glutils.generate_germline_set(sglfo, args.n_genes_per_region, args.n_sim_alleles_per_gene, args.min_allele_prevalence_freq, allele_prevalence_fname, new_allele_info=args.new_allele_info, remove_template_genes=remove_template_genes)
        cmd_str += ' --allele-prevalence-fname ' + allele_prevalence_fname
    else:
        sglfo = glutils.read_glfo('data/germlines/human', locus=args.locus, only_genes=(args.sim_v_genes + args.dj_genes))
        added_snp_names = glutils.generate_new_alleles(sglfo, args.new_allele_info, debug=True, remove_template_genes=args.remove_template_genes)

        if args.allele_prevalence_freqs is not None:
            if not utils.is_normed(args.allele_prevalence_freqs):
                raise Exception('--allele-prevalence-freqs %s not normalized' % args.allele_prevalence_freqs)
            if len(args.allele_prevalence_freqs) != len(sglfo['seqs']['v']):  # already checked when parsing args, but, you know...
                raise Exception('--allele-prevalence-freqs %d not the same length as sglfo %d' % (len(args.allele_prevalence_freqs), len(sglfo['seqs']['v'])))
            gene_list = sorted(sglfo['seqs']['v']) if len(added_snp_names) == 0 else list(set(args.sim_v_genes)) + added_snp_names
            prevalence_freqs = {'v' : {g : f for g, f in zip(gene_list, args.allele_prevalence_freqs)}, 'd' : {}, 'j' : {}}
            glutils.write_allele_prevalence_freqs(prevalence_freqs, allele_prevalence_fname)
            cmd_str += ' --allele-prevalence-fname ' + allele_prevalence_fname

    print '  simulating with %d v: %s' % (len(sglfo['seqs']['v']), ' '.join([utils.color_gene(g) for g in sglfo['seqs']['v']]))
    glutils.write_glfo(args.outdir + '/germlines/simulation', sglfo)
    cmd_str += ' --initial-germline-dir ' + args.outdir + '/germlines/simulation'

    # run simulation
    if args.seed is not None:
        cmd_str += ' --seed ' + str(args.seed)
    utils.simplerun(cmd_str, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_other_method(args, method):
    if method not in ['tigger-default', 'tigger-tuned', 'igdiscover']:  # really just to make it easier to search for this fcn
        assert False
    if utils.output_exists(args, get_outfname(args, method)):
        return
    simfasta = utils.getprefix(args.simfname) + '.fa'
    utils.csv_to_fasta(args.simfname, outfname=simfasta, overwrite=False, remove_duplicates=True)
    cmd = './test/%s-run.py' % method.split('-')[0]
    if method == 'tigger-tuned':
        cmd += ' --tuned-tigger-params'
    cmd += ' --infname ' + simfasta
    cmd += ' --outfname ' + get_outfname(args, method)
    if args.overwrite:
        cmd += ' --overwrite'
    if args.gls_gen:
        cmd += ' --gls-gen'
        cmd += ' --glfo-dir ' + partis_dir + '/data/germlines/human'  # the partis mehods have this as the default internally, but we want/have to set it explicitly here
    else:
        cmd += ' --glfo-dir ' + args.inf_glfo_dir
    cmd += ' --simulation-germline-dir ' + args.outdir + '/germlines/simulation'  # alleleclusterer is the only one that really uses this, but for now I want its dbg output to have the sim info
    if method != 'igdiscover':  # for now we're saving all the igdiscover output/intermediate files, so we write them to an output dir
        cmd += ' --workdir ' + args.workdir + '/' + method
    cmd += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd += ' --slurm'

    utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_partis(args, method):
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
    if method == 'partis':
        cmd_str += ' --debug-allele-finding' # --always-find-new-alleles'
        cmd_str += ' --is-simu --simulation-germline-dir ' + args.outdir + '/germlines/simulation'  # alleleclusterer is the only one that really uses this, but for now I want its dbg output to have the sim info
    elif method == 'full':
        cmd_str += ' --leave-default-germline'
    else:
        assert False

    cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.n_max_queries is not None:
        cmd_str += ' --n-max-queries ' + str(args.n_max_queries)  # NOTE do *not* use --n-random-queries, since it'll change the cluster size distribution
    if args.slurm:
        cmd_str += ' --batch-system slurm'

    if not args.gls_gen:  # otherwise it uses the default (full) germline dir
        cmd_str += ' --initial-germline-dir ' + args.inf_glfo_dir  # --dont-remove-unlikely-alleles

    cmd_str += ' --parameter-dir ' + paramdir
    cmd_str += ' --only-overall-plots --plotdir ' + plotdir
    if args.seed is not None:
        cmd_str += ' --seed ' + str(args.seed)
    if args.plot_and_fit_absolutely_everything is not None:
        cmd_str += ' --plot-and-fit-absolutely-everything ' + str(args.plot_and_fit_absolutely_everything)
    utils.simplerun(cmd_str, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def write_inf_glfo(args):  # read default glfo, restrict it to the specified alleles, and write to somewhere where all the methods can read it
    # NOTE this dir should *not* be modified by any of the methods
    inf_glfo = glutils.read_glfo('data/germlines/human', locus=args.locus, only_genes=args.inf_v_genes + args.dj_genes)
    print '  writing initial inference glfo with %d v: %s' % (len(inf_glfo['seqs']['v']), ' '.join([utils.color_gene(g) for g in inf_glfo['seqs']['v']]))
    glutils.write_glfo(args.inf_glfo_dir, inf_glfo)

# ----------------------------------------------------------------------------------------
def run_tests(args):
    print 'seed %d' % args.seed
    # all fcns return immediately if output already exists

    if 'simu' in args.methods:
        simulate(args)
        args.methods.remove('simu')

    if not args.gls_gen:
        write_inf_glfo(args)
    for method in args.methods:
        if method == 'partis' or method == 'full':
            run_partis(args, method)
        else:
            run_other_method(args, method)

# ----------------------------------------------------------------------------------------
def multiple_tests(args):
    def getlogdir(iproc):
        return args.outdir + '/' + str(iproc) + '/logs/' + '-'.join(args.methods)
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
    if args.dry_run:
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
    print '  look for logs in %s' % args.outdir
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
#     subseqs = numpy.random.choice(seqs, size=ntrees)
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
                             './bin/test-allele-finding.py --n-sim-events 2000 --n-procs 10 --sim-v-genes=IGHV1-18*01 --inf-v-genes=IGHV1-18*01 --snp-positions 27,55,88',
                             './bin/test-allele-finding.py --n-sim-events 2000 --n-procs 10 --sim-v-genes=IGHV4-39*01:IGHV4-39*02 --inf-v-genes=IGHV4-39*01'])
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=example_str)
parser.add_argument('--n-sim-events', type=int, default=20, help='number of simulated rearrangement events')
parser.add_argument('--n-max-queries', type=int, help='number of queries to use for inference from the simulation sample')
parser.add_argument('--n-leaves', type=float, default=1.)
parser.add_argument('--n-leaf-distribution')
parser.add_argument('--n-procs', type=int, default=2)
parser.add_argument('--seed', type=int, default=int(time.time()))
parser.add_argument('--gls-gen', action='store_true', help='generate a random germline set from scratch (parameters specified above), and infer a germline set from scratch, instead of using --sim-v-genes, --dj-genes, --inf-v-genes, and --snp-positions.')
parser.add_argument('--sim-v-genes', default='IGHV4-39*01:IGHV4-39*05', help='.')
parser.add_argument('--inf-v-genes', default='IGHV4-39*01', help='.')
parser.add_argument('--dj-genes', default='IGHD6-19*01:IGHJ4*02', help='.')
parser.add_argument('--snp-positions', help='colon-separated list (length must equal length of <--sim-v-genes>) of comma-separated snp positions for each gene, e.g. for two genes you might have \'3,71:45\'')
parser.add_argument('--nsnp-list', help='colon-separated list (length must equal length of <--sim-v-genes> unless --gls-gen) of the number of snps to generate for each gene (each snp at a random position). If --gls-gen, then this still gives the number of snpd genes, but it isn\'t assumed to be the same length as anything [i.e. we don\'t yet know how many v genes there\'ll be]')
parser.add_argument('--indel-positions', help='see --snp-positions (a.t.m. the indel length distributions are hardcoded)')
parser.add_argument('--nindel-list', help='see --nsnp-list')
parser.add_argument('--n-genes-per-region', default='20:5:3')
parser.add_argument('--n-sim-alleles-per-gene', default='1,2:1,2:1,2')
parser.add_argument('--min-allele-prevalence-freq', type=float, default=0.1)
parser.add_argument('--allele-prevalence-freqs', help='colon-separated list of allele prevalence frequencies, including newly-generated snpd genes (ordered alphabetically)')
parser.add_argument('--remove-template-genes', action='store_true', help='when generating snps, remove the original gene before simulation')
parser.add_argument('--mut-mult', type=float)
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--methods', default='simu:partis')
parser.add_argument('--outdir', default=utils.fsdir() + '/partis/allele-finder')
parser.add_argument('--inf-glfo-dir', help='default set below')
parser.add_argument('--simfname', help='default set below')
parser.add_argument('--workdir', default=utils.fsdir() + '/_tmp/hmms/' + str(random.randint(0, 999999)))
parser.add_argument('--n-tests', type=int)
parser.add_argument('--iteststart', type=int, default=0)
parser.add_argument('--plot-and-fit-absolutely-everything', type=int, help='fit every single position for this <istart> and write every single corresponding plot (slow as hell, and only for debugging/making plots for paper)')
parser.add_argument('--partis-path', default='./bin/partis')
parser.add_argument('--locus', default='igh')

args = parser.parse_args()
args.dj_genes = utils.get_arg_list(args.dj_genes)
args.sim_v_genes = utils.get_arg_list(args.sim_v_genes)
args.inf_v_genes = utils.get_arg_list(args.inf_v_genes)
args.allele_prevalence_freqs = utils.get_arg_list(args.allele_prevalence_freqs, floatify=True)
args.methods = utils.get_arg_list(args.methods)
available_methods = set(['simu', 'partis', 'full', 'tigger-default', 'tigger-tuned', 'igdiscover'])
if len(set(args.methods) - available_methods) > 0:
    raise Exception('unexpected --methods: %s' % ' '.join(set(args.methods) - available_methods))

positions = {
    'snp' : utils.get_arg_list(args.snp_positions),
    'indel' : utils.get_arg_list(args.indel_positions),
}
numbers = {
    'snp' : utils.get_arg_list(args.nsnp_list, intify=True),
    'indel' : utils.get_arg_list(args.nindel_list, intify=True),
}
args.snp_positions = None  # just to make sure you don't accidentally use them
args.indel_positions = None
args.nsnp_list = None
args.nindel_list = None
# if we're generating/inferring a whole germline set these are either set automatically or not used
if args.gls_gen:
    args.sim_v_genes = None
    args.inf_v_genes = None
    args.dj_genes = None
    args.allele_prevalence_freqs = None
    args.inf_glfo_dir = None

n_new_alleles = None
mtypes = ['snp', 'indel']
for mtype in mtypes:
    if positions[mtype] is not None:  # if specific positions were specified on the command line
        positions[mtype] = [[int(p) for p in pos_str.split(',')] for pos_str in positions[mtype]]
        if len(positions[mtype]) != len(args.sim_v_genes):
            raise Exception('--%s-positions %s and --sim-v-genes %s not the same length (%d vs %d)' % (mtype, positions[mtype], args.sim_v_genes, len(positions[mtype]), len(args.sim_v_genes)))
    if numbers[mtype] is not None:
        if not args.gls_gen and len(numbers[mtype]) != len(args.sim_v_genes):
            raise Exception('--n%s-list %s and --sim-v-genes %s not the same length (%d vs %d)' % (mtype, numbers[mtype], args.sim_v_genes, len(numbers[mtype]), len(args.sim_v_genes)))
        if positions[mtype] is not None:
            raise Exception('can\'t specify both --n%s-list and --%s-positions' % (mtype, mtype))
        positions[mtype] = [[None for _ in range(number)] for number in numbers[mtype]]  # the <None> tells glutils to choose a position at random
    if positions[mtype] is not None:
        if n_new_alleles is None:
            n_new_alleles = len(positions[mtype])
        if len(positions[mtype]) != n_new_alleles:
            raise Exception('mismatched number of new alleles for %s' % ' vs '.join(mtypes))
if n_new_alleles is None:
    n_new_alleles = 0
for mtype in mtypes:
    if positions[mtype] is None:  # if it wasn't specified at all, i.e. we don't want to generate any new alleles
        positions[mtype] = [[] for _ in range(n_new_alleles)]
args.new_allele_info = [{'gene' : args.sim_v_genes[igene] if not args.gls_gen else None,
                         'snp-positions' : positions['snp'][igene],
                         'indel-positions' : positions['indel'][igene]}
                        for igene in range(n_new_alleles)]

if args.allele_prevalence_freqs is not None:
    # easier to check the length after we've generated snpd genes (above)
    if not utils.is_normed(args.allele_prevalence_freqs):
        raise Exception('--allele-prevalence-freqs %s not normalized' % args.allele_prevalence_freqs)
if args.inf_glfo_dir is None:
    args.inf_glfo_dir = args.outdir + '/germlines/inference'
if args.simfname is None:
    args.simfname = args.outdir + '/simu.csv'

if args.seed is not None:
    random.seed(args.seed)
    numpy.random.seed(args.seed)

if args.n_tests is not None:
    multiple_tests(args)
else:
    run_tests(args)
