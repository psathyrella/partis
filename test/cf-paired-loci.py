#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math
import collections
import multiprocessing
import yaml
import operator
import random
from io import open

import partis.utils as utils
import partis.paircluster as paircluster
import partis.scanplot as scanplot
import partis.plotting as plotting
import partis.clusterpath as clusterpath
import partis.glutils as glutils

# ----------------------------------------------------------------------------------------
script_base = os.path.basename(__file__).replace('cf-', '').replace('.py', '')
partition_types = ['single', 'joint']
all_perf_metrics = ['precision', 'sensitivity', 'f1', 'time-reqd', 'naive-hdist', 'cln-frac', 'n-clusters', 'pairwise-prec', 'pairwise-sens', 'pairwise-f1', 'coar', 'rf', 'mrca', 'n-pars-trees']  # pcfrac-*: pair info cleaning correct fraction, cln-frac: collision fraction
pcfrac_metrics = ['pcfrac-%s%s'%(t, s) for s in ['', '-ns'] for t in ['correct', 'mispaired', 'unpaired', 'correct-family', 'near-family']]  # '-ns': non-singleton
all_perf_metrics += pcfrac_metrics
synth_actions = ['synth-%s'%a for a in ['distance-0.00', 'distance-0.005', 'distance-0.01', 'distance-0.02', 'distance-0.03', 'reassign-0.10', 'singletons-0.40', 'singletons-0.20']]
ptn_actions = ['partition', 'partition-lthresh', 'star-partition', 'vsearch-partition', 'subset-partition', 'vsearch-subset-partition', 'annotate', 'vjcdr3-0.9', 'vjcdr3-0.8', 'scoper', 'mobille', 'igblast', 'linearham', 'enclone'] + synth_actions  # using the likelihood (rather than hamming-fraction) threshold makes basically zero difference
print('%s default iqtree in treeutils.py is now 3 (rather than 1.6.12), which may require changing things here to get identical results as in the past)' % utils.wrnstr())
phylo_actions = ['iqtree', 'raxml', 'gctree', 'gctree-mut-mult', 'gctree-no-dag', 'igphyml']  # , 'iqtree-1.6.beta3', 'iqtree-2.3.1'  # , 'gctree-base'
tree_perf_actions = ['%s-tree-perf'%a for a in phylo_actions]  # it would be really nice to run tree perf during the phylo action, but i can't figure out a good way to do that (main problem is getting access to both true and inferred annotations in a sensible way)
phylo_naive_actions = ['phylo-naive-iqtree', 'phylo-naive-iqtree-fuzz-0', 'phylo-naive-iqtree-fuzz-1', 'phylo-naive-iqtree-fuzz-2']
after_actions = ['replay-plot', 'cache-parameters', 'merge-paired-partitions', 'get-selection-metrics', 'parse-linearham-trees', 'write-fake-paired-annotations', 'tree-perf']  + ptn_actions + phylo_actions + tree_perf_actions + phylo_naive_actions  # actions that come after simulation (e.g. partition)
plot_actions = ['single-chain-partis', 'single-chain-scoper']
def is_single_chain(action):
    return 'synth-' in action or 'vjcdr3-' in action or 'single-chain-' in action or action in ['mobille', 'igblast', 'linearham']
def is_joint_method(action):
    return action in ['enclone', 'single-chain-partis']  # yes, 'single-chain-partis' is a "joint method": this is just a hackey way to get the single chain line on the joint plot, we don't actually want it [twice] on the single chain plot

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simu:cache-parameters:partition:plot')  # can also be merge-paired-partitions and get-selection-metrics
parser.add_argument('--merge-paired-partitions', action='store_true', help='for partis partition actions, don\'t re-partition, just merge paired partitions')
parser.add_argument('--n-sim-events-list', default='10', help='N sim events in each repertoire/"proc"/partis simulate run')
parser.add_argument('--n-leaves-list', help='NOTE can use either this or \'n-leaf-distribution\' for \'hist\' n leaf distr (and depending on zip vars you need to use one or the other)') #'2:3:4:10') #1 5; do10)
parser.add_argument('--n-sim-seqs-per-generation-list')  # only for bcr-phylo
parser.add_argument('--carry-cap-list')  # only for bcr-phylo
parser.add_argument('--constant-number-of-leaves-list')
parser.add_argument('--n-leaf-distribution-list', help='NOTE can use either this or \'n-leaves\' for \'hist\' n leaf distr (and depending on zip vars you need to use one or the other)')
parser.add_argument('--mean-cells-per-droplet-list') #, default='None')
parser.add_argument('--fraction-of-reads-to-remove-list')
parser.add_argument('--bulk-data-fraction-list')
parser.add_argument('--allowed-cdr3-lengths-list') # NOTE need to replace , with M and : with L in each element e.g. ighM15-34:ighM36-43:ighM45-52:ighM54-61:ighM63-73
parser.add_argument('--n-genes-per-region-list')
parser.add_argument('--n-sim-alleles-per-gene-list')
parser.add_argument('--scratch-mute-freq-list') #, type=float, default=1)
parser.add_argument('--mutation-multiplier-list') #, type=float, default=1)
parser.add_argument('--obs-times-list')  # only for bcr-phylo
parser.add_argument('--tree-imbalance-list')
parser.add_argument('--context-depend-list')
parser.add_argument('--target-distance-list')
parser.add_argument('--aa-paratope-positions-list')
parser.add_argument('--aa-struct-positions-list')
parser.add_argument('--leaf-sampling-scheme-list')
parser.add_argument('--n-naive-seq-copies-list')
parser.add_argument('--dont-observe-common-ancestors-list')
parser.add_argument('--biggest-naive-seq-cluster-to-calculate-list')
parser.add_argument('--biggest-logprob-cluster-to-calculate-list')
parser.add_argument('--prep', action='store_true', help='only for mobille run script atm')
parser.add_argument('--antn-perf', action='store_true', help='calculate annotation performance values')
parser.add_argument('--calc-antns', action='store_true', help='calculate annotations')
parser.add_argument('--simu-type', default='partis', choices=['partis', 'bcr-phylo', 'gcdyn'], help='which simulation package to use: partis (i.e. TreeSim/bpp), bcr-phylo, or gcdyn')
parser.add_argument('--dont-observe-common-ancestors', action='store_true')
parser.add_argument('--data-cluster-size-hist-fname', default='/fh/fast/matsen_e/processed-data/partis/goo-dengue-10x/count-params-v0/d-14/parameters/igh+igk/igh/hmm/cluster_size.csv') #/fh/fast/matsen_e/processed-data/partis/10x-examples/v1/hs-1-postvax/parameters/igh+igk/igh/hmm/cluster_size.csv')  # ick ick ick
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--bcrham-time', action='store_true')
utils.add_scanvar_args(parser, script_base, all_perf_metrics)
args = parser.parse_args()
args.scan_vars = {
    'simu' : ['seed', 'n-leaves', 'n-sim-seqs-per-generation', 'carry-cap', 'constant-number-of-leaves', 'n-leaf-distribution', 'scratch-mute-freq', 'mutation-multiplier', 'obs-times', 'tree-imbalance', 'context-depend', 'target-distance', 'aa-paratope-positions', 'aa-struct-positions', 'leaf-sampling-scheme', 'n-naive-seq-copies', 'dont-observe-common-ancestors', 'mean-cells-per-droplet', 'fraction-of-reads-to-remove', 'bulk-data-fraction', 'allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene', 'n-sim-events', 'dataset-in'],
    'cache-parameters' : ['biggest-naive-seq-cluster-to-calculate', 'biggest-logprob-cluster-to-calculate'],  # only really want these in 'partition', but this makes it easier to point at the right parameter dir
    'partition' : ['biggest-naive-seq-cluster-to-calculate', 'biggest-logprob-cluster-to-calculate'],
}
args.str_list_vars = ['allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene', 'n-sim-seqs-per-generation', 'obs-times', 'aa-paratope-positions', 'aa-struct-positions']
args.recurse_replace_vars = ['allowed-cdr3-lengths']  # ick ick ick
args.bool_args = ['constant-number-of-leaves', 'dont-observe-common-ancestors']  # NOTE different purpose to svartypes below (this isn't to convert all the values to the proper type, it's just to handle flag-type args
if 'all-pcfrac' in args.perf_metrics:
    args.perf_metrics = args.perf_metrics.replace('all-pcfrac', ':'.join(pcfrac_metrics))
args.paired_loci = True
for targ in ['actions', 'plot_metrics']:
    tlist = getattr(args, targ).split(':')
    if 'tree-perf' in tlist:
        tpi = tlist.index('tree-perf')
        setattr(args, targ, ':'.join(tlist[: tpi] + tree_perf_actions + tlist[tpi + 1 :]))
utils.process_scanvar_args(args, after_actions, plot_actions, all_perf_metrics)
if args.antn_perf:
    if any(a in phylo_actions for a in args.actions):
        raise Exception('can\'t set --antn-perf for phylo actions, since --antn-perf requires --is-simu and thus reads true annotations, but phylo actions infer trees on inferred annotations')
    print('  --antn-perf: turning on --make-plots')
    args.make_plots = True

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    actstr = action
    if action in ['cache-parameters', 'partition', 'single-chain-partis'] + phylo_actions:
        actstr = 'partis'
    elif action == 'single-chain-scoper':
        actstr = 'scoper'
    return '%s/%s' % (utils.svoutdir(args, varnames, vstrs, action), actstr)

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action, locus=None, single_chain=False, single_file=False, logfile=False, pmetr=None):
    outdir = odir(args, varnames, vstrs, action)
    if action == 'cache-parameters' and not logfile:
        outdir += '/parameters'
        if locus is None and not single_file:
            return outdir
    if single_file:
        assert locus is None
        locus = 'igk'
    assert locus is not None or action in phylo_actions + tree_perf_actions + ['write-fake-paired-annotations', 'replay-plot']
    if logfile:
        ofn = '%s/%s%s.log' % (outdir, 'work/%s/'%locus if action=='mobille' else '', action)
    elif pmetr is not None and 'pcfrac-' in pmetr:
        ofn = '%s/true-pair-clean-performance.csv' % outdir #, pmetr.replace('pcfrac', '').replace('-ns', '') if 'pcfrac-' in pmetr else '')
    elif pmetr is not None and pmetr == 'naive-hdist':
        if 'phylo-naive' in action:
            ofn = '%s/%s/mutation/hamming_to_true_naive.csv' % (outdir, locus)
        else:
            ofn = '%s/single-chain/plots/%s/hmm/mutation/hamming_to_true_naive.csv' % (outdir, locus)
    elif action == 'cache-parameters':
        ofn = '%s/%s' % (outdir, locus)
        if single_file:
            # ofn += '/hmm/germline-sets/%s/%sv.fasta' % (locus, locus)
            ofn += '/hmm/all-mean-mute-freqs.csv'
    elif action in phylo_actions:
        ofn = '%s/%s-annotations.yaml' % (outdir, action)  # there's also the selection-metrics.yaml file in the same dir, but i think atm I'm only using this one so don't need an option here to choose
    elif action in tree_perf_actions:
        ofn = '%s/tree-perf-vals.yaml' % outdir
    elif action == 'write-fake-paired-annotations':
        ofn = '%s/fake-paired-annotations.yaml' % outdir
    elif action == 'replay-plot':
        ofn = '%s/diff-vals.yaml' % outdir
    elif 'phylo-naive' in action:  # duplicates naive-hdist bit above
        ofn = '%s/%s/mutation/hamming_to_true_naive.csv' % (outdir, locus)
    else:
        ofn = paircluster.paired_fn(outdir, locus, suffix='.yaml', actstr=None if action=='simu' else 'partition', single_chain=single_chain or is_single_chain(action))
    return ofn

# ----------------------------------------------------------------------------------------
# distance-based partitions get made by running partis, but then we make the other types here
def make_synthetic_partition(action, varnames, vstrs):
    for ltmp in plot_loci():
        _, _, true_cpath = utils.read_output(paircluster.paired_fn(odir(args, varnames, vstrs, 'simu'), ltmp, suffix='.yaml'))
        _, mistype, misfrac = action.split('-')
        new_partition = utils.generate_incorrect_partition(true_cpath.best(), float(misfrac), mistype)
        new_cpath = clusterpath.ClusterPath(partition=new_partition)
        new_cpath.calculate_missing_values(true_partition=true_cpath.best())
        ofn = ofname(args, varnames, vstrs, action, locus=ltmp, single_chain=True)
        utils.write_annotations(ofn, None, [], utils.annotation_headers, partition_lines=new_cpath.get_partition_lines())  # could now use utils.write_only_partition()
        print('    %s: wrote synthetic partition to %s' % (ltmp, ofn))

# ----------------------------------------------------------------------------------------
def get_replacefo():  # ick
    if args.tree_imbalance_list is None:
        return None
    rfo = {}
    for imb in args.tree_imbalance_list:
        if imb in ['None', None]:  # use None to leave default/not set
            rfo[imb] = None
        else:
            rfo[imb] = ['input-simulation-treefname', imbalfname(imb)]  # NOTE atm this has to be length 2
    return {'tree-imbalance' : rfo}

# ----------------------------------------------------------------------------------------
def get_cmd(action, base_args, varnames, vlists, vstrs, synth_frac=None):
    if action == 'scoper':
        cmd = './test/scoper-run.py --indir %s --outdir %s --simdir %s' % (ofname(args, varnames, vstrs, 'cache-parameters'), odir(args, varnames, vstrs, action), odir(args, varnames, vstrs, 'simu'))
        return cmd
    if action in ['mobille', 'igblast', 'linearham', 'enclone'] + tree_perf_actions + phylo_naive_actions:
        binstr = 'mobille-igblast' if action in ['mobille', 'igblast'] else ('tree-perf' if action in tree_perf_actions else action)
        if action in phylo_naive_actions:
            binstr = 'phylo-naive'
        if action in tree_perf_actions:
            cmd = './bin/%s-run.py --true-tree-file %s --inferred-tree-file %s' % (binstr, ofname(args, varnames, vstrs, 'write-fake-paired-annotations', single_file=True), ofname(args, varnames, vstrs, action.replace('-tree-perf', ''), single_file=True))
        else:
            cmd = './test/%s-run.py --simdir %s' % (binstr, odir(args, varnames, vstrs, 'simu'))
        cmd += ' --outdir %s' % odir(args, varnames, vstrs, action)
        if action in ['mobille', 'igblast']:  # i don't think both of them need all these
            cmd += ' --id-str %s --base-imgt-outdir %s' % ('_'.join('%s-%s'%(n, s) for n, s in zip(varnames, vstrs)), '%s/%s/%s/imgt-output' % (args.base_outdir, args.label, args.version))
        if action == 'linearham':
            if not args.antn_perf:
                raise Exception('running linearham action without --antn-perf set, which means you likely also didn\'t set it for the partition action (although it has on direct effect on the linearham action)')
            cmd += ' --partis-outdir %s --n-sim-events %d' % (odir(args, varnames, vstrs, 'partition'), int(utils.vlval(args, vlists, varnames, 'n-sim-events')))
            cmd += ' --docker --local-docker-image'
        if 'phylo-naive' in action:
            tree_inf_method = action.replace('phylo-naive-', '')
            if 'fuzz-' in action:
                tree_inf_method = '-'.join(tree_inf_method.split('-')[:-2])
            import partis.treeutils as treeutils
            assert tree_inf_method in treeutils.all_phylo_methods
            cmd += ' --input-partition-dir %s --tree-inference-method %s' % (odir(args, varnames, vstrs, 'partition'), tree_inf_method)
            if 'fuzz' in action:
                cmd += ' --phylo-naive-inference-fuzz %s' % action.split('-')[-1]
        if args.n_sub_procs > 1:
            cmd += ' --n-procs %d' % args.n_sub_procs
        if args.prep:
            cmd += ' --prep'
            # then do this by hand, and submit to imgt/high vquest by hand, then download results and put them in the proper dir (run mobille run script to get dir)
            # tar czf /path/somewhere/to/rsync/imgt-input.tgz /fh/local/dralph/partis/paired-loci/vs-shm/v2/seed-*/scratch-mute-freq-*/mobille/work/*/imgt-input/*.fa
        if args.overwrite:
            cmd += ' --overwrite'
        return cmd
    actstr = action
    if 'synth-distance-' in action or action in ['vsearch-partition', 'partition-lthresh', 'star-partition']:
        actstr = 'partition'
    if action == 'vsearch-subset-partition':
        actstr = 'subset-partition'
    if 'vjcdr3-' in action:
        actstr = 'annotate'
    if args.merge_paired_partitions:
        actstr = 'merge-paired-partitions'
    if action in phylo_actions:
        actstr = 'infer-trees'
    if action == 'simu' and args.simu_type == 'gcdyn':
        cmd = './bin/gcdyn-simu-run.py --rm-last-ighj-base --outdir %s --n-sub-procs %d' % (odir(args, varnames, vstrs, action), args.n_sub_procs)
    else:
        binstr, actstr, odstr = ('bcr-phylo-run.py', '--actions %s'%actstr, 'base') if args.simu_type=='bcr-phylo' and action=='simu' else ('partis', actstr.replace('simu', 'simulate'), 'paired')
        cmd = './bin/%s %s --paired-loci --%s-outdir %s --n-procs %d' % (binstr, actstr, odstr, odir(args, varnames, vstrs, action), args.n_sub_procs)
    if action == 'replay-plot':
        cmd = './projects/replay-plot.py --min-seqs-per-gc 70 --normalize --n-max-simu-trees 61 --simu-like-dir %s --outdir %s --bcr-phylo' % (os.path.dirname(ofname(args, varnames, vstrs, 'write-fake-paired-annotations')), os.path.dirname(ofname(args, varnames, vstrs, action)))
        if args.simu_type == 'gcdyn':
            cmd += ' --default-naive-affinity 0'
    elif action == 'simu':
        if args.simu_type == 'partis':
            cmd += ' --simulate-from-scratch --no-per-base-mutation'
        cmd += ' %s' % ' '.join(base_args)
        if args.simu_type in ['partis', 'bcr-phylo'] and args.single_light_locus is not None:
            cmd += ' --single-light-locus %s' % args.single_light_locus
        if args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
        for vname, vstr in zip(varnames, vstrs):
            cmd = utils.add_to_scan_cmd(args, vname, vstr, cmd, replacefo=get_replacefo())
        tmp_astrs = [a for a in ['--n-leaves', '--n-leaf-distribution'] if a in cmd.split() and utils.get_val_from_arglist(cmd.split(), a) == 'hist' ]
        if len(tmp_astrs) > 0:
            astr = utils.get_single_entry(tmp_astrs)  # i think it'll break/doesn't make sense if there's more than one
            if astr == '--n-leaves':
                cmd = ' '.join(utils.remove_from_arglist(cmd.split(), '--n-leaves', has_arg=True))
                cmd += ' --n-leaf-distribution hist'
            cmd += ' --n-leaf-hist-fname %s' %  args.data_cluster_size_hist_fname
        if args.simu_type == 'bcr-phylo':
            # raise Exception('need to fix duplicate uids coming from bcr-phylo (they get modified in seqfileopener, which is ok, but then the uids in the final partition don\'t match the uids in the true partition')
            cmd += ' --dont-get-tree-metrics --only-csv-plots --mutated-outpath --min-ustr-len 20'  # NOTE don't increase the mutation rate it makes everything terminate early  --base-mutation-rate 1'  # it's nice to jack up the mutation rate so we get more mutations in less time (higher than this kills off all leaves, not sure why, altho i'm sure it's obvious if i thought about it)
            if args.dont_observe_common_ancestors:
                cmd += ' --dont-observe-common-ancestors'
            if args.overwrite:
                cmd += ' --overwrite'
            # cmd += ' --debug 1'
    else:
        cmd += ' --paired-indir %s' % odir(args, varnames, vstrs, 'simu')
        if action in ['vsearch-partition', 'vsearch-subset-partition']:
            cmd += ' --naive-vsearch'  # same as --fast
        if 'synth-distance-' in action:
            synth_hfrac = float(action.replace('synth-distance-', ''))
            cmd += ' --synthetic-distance-based-partition --naive-hamming-bounds %.2f:%.2f' % (synth_hfrac, synth_hfrac)
        if 'vjcdr3-' in action:
            cmd += ' --annotation-clustering --annotation-clustering-threshold %.2f' % float(action.split('-')[1])
        if action == 'partition-lthresh':
            cmd += ' --paired-naive-hfrac-threshold-type likelihood'
        if action == 'star-partition':
            cmd += ' --subcluster-annotation-size None'
        if action not in phylo_actions:  # we want to infer the trees, so need to ignore true annotations
            cmd += ' --is-simu'
        if action not in ['cache-parameters', 'subset-partition', 'vsearch-subset-partition']:
            cmd += ' --refuse-to-cache-parameters'
        if 'synth-distance-' in action or 'vjcdr3-' in action or action in ['vsearch-partition', 'partition-lthresh', 'star-partition', 'annotate']:
            cmd += ' --parameter-dir %s' % ofname(args, varnames, vstrs, 'cache-parameters')
        if action in ptn_actions and 'vjcdr3-' not in action and action not in ['subset-partition', 'vsearch-subset-partition'] and not args.make_plots and not args.antn_perf and not args.calc_antns:
            cmd += ' --dont-calculate-annotations'
        if action == 'write-fake-paired-annotations':
            cmd += ' --extra-annotation-columns tree'
        if args.simu_type == 'gcdyn':
            cmd += ' --species mouse --initial-germline-dir datascripts/meta/taraki-gctree-2021-10/germlines --no-insertions-or-deletions'
        for vname, vstr in zip(varnames, vstrs):
            if vname in args.scan_vars['simu']:
                continue
            if action == 'cache-parameters' and vname in ['biggest-naive-seq-cluster-to-calculate', 'biggest-logprob-cluster-to-calculate']:
                continue  # ick ick ick
            cmd = utils.add_to_scan_cmd(args, vname, vstr, cmd, replacefo=get_replacefo())
        if args.make_plots and action != 'cache-parameters':
            cmd += ' --plotdir paired-outdir'
            if action in ptn_actions:
                cmd += ' --no-partition-plots'
            if args.antn_perf:
                cmd += ' --plot-annotation-performance --only-csv-plots'
        if args.inference_extra_args is not None:
            cmd += ' %s' % args.inference_extra_args
        if action in phylo_actions:  # should maybe remove plotdir and annotation performance?
            cmd = ' '.join(utils.remove_from_arglist(cmd.split(), '--simultaneous-true-clonal-seqs'))  # --is-simu can't be set for this, but if it isn't set then this causes it to crash
            cmd += ' --tree-inference-method %s --infer-trees-with-only-leaves' % action  # some simulation writes seqs for internal nodes (which we need), but we don't want to pass them to the tree inference methods

    return cmd

# ----------------------------------------------------------------------------------------
# TODO combine this also with fcns in cf-tree-metrics.py (and put it in scanplot)
def run_scan(action):
    base_args, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars[action])
    cmdfos = []
    print('  %s: running %d combinations of: %s' % (utils.color('blue_bkg', action), len(valstrs), ' '.join(varnames)))
    if args.debug:
        print('   %s' % ' '.join(varnames))
    n_already_there, existing_ofn = 0, None
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print('   %s' % ' '.join(vstrs))

        single_file, locus = (True, None) if args.data_in_cfg is None else (False, args.data_in_cfg[vstrs[0]]['locus'])
        ofn = ofname(args, varnames, vstrs, action, single_file=single_file, locus=locus)
        if args.merge_paired_partitions:
            assert action in ['partition', 'vsearch-partition']
        else:
            if utils.output_exists(args, ofn, debug=False):
                n_already_there += 1
                existing_ofn = ofn
                continue

        if 'synth-reassign-' in action or 'synth-singletons-' in action:
            make_synthetic_partition(action, varnames, vstrs)
            continue

        if action == 'simu' and args.data_in_cfg is not None:
            utils.make_scanvar_data_links(args, varnames, vstrs, action, ofn, ofname, odir)

        cmd = get_cmd(action, base_args, varnames, val_lists, vstrs)
        # utils.simplerun(cmd, logfname='%s-%s.log'%(odir(args, varnames, vstrs, action), action), dryrun=args.dry)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, varnames, vstrs, action),
            'workdir' : '%s/partis-work/%d' % (args.workdir, icombo),
        }]

    utils.run_scan_cmds(args, cmdfos, '%s.log'%(action if not args.merge_paired_partitions else 'merge-paired-partitions'), len(valstrs), n_already_there, existing_ofn)

# ----------------------------------------------------------------------------------------
def imbalfname(ibval):
    return '%s/linearham-simulation/processed-trees/imbal-%s.trees' % (args.base_outdir, ibval)

# ----------------------------------------------------------------------------------------
def parse_linearham_trees():
    # downloaded from here https://zenodo.org/record/3746832, then cat'd together all the trees
    # here we sort by the imbalance, then put in separated files rounded to second decimal place
    import partis.treeutils as treeutils
    ibtrees, ibvals = {}, []
    with open('/fh/local/dralph/partis/paired-loci/linearham-simulation/simulation_trees/all.trees') as tfile:
        for line in tfile:
            dtr = treeutils.get_dendro_tree(treestr=line.strip(), no_warn=True) #debug=True)
            imbal = treeutils.get_imbalance(dtr)
            rval = '%.2f' % imbal #utils.round_to_n_digits(imbal, 1)
            if rval not in ibtrees:
                ibtrees[rval] = []
            ibtrees[rval].append(dtr)
            ibvals.append(imbal)
    print(' '.join(['%.3f'%v for v in sorted(ibvals)]))
    utils.mkdir(imbalfname('xxx'), isfile=True)
    print('    writing trees to %s' % os.path.dirname(imbalfname('xxx')))
    for ibval, treelist in sorted(list(ibtrees.items()), key=operator.itemgetter(0)):
        print('      %3d trees to %s' % (len(treelist), os.path.basename(imbalfname(ibval))))
        with open(imbalfname(ibval), 'w') as ofile:
            for dtr in treelist:
                ofile.write(dtr.as_string(schema='newick'))

# ----------------------------------------------------------------------------------------
def plot_loci(methods=None):
    if methods is not None and any(m in tree_perf_actions for m in methods):  # maybe should be all(), i dunno
        return [None]
    if args.single_light_locus is None:
        return utils.sub_loci('ig')
    else:
        return [utils.heavy_locus('ig'), args.single_light_locus]

# ----------------------------------------------------------------------------------------
def get_fnfcn(method, locus, ptntype, pmetr):
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, method, locus=locus, single_chain=ptntype=='single', logfile=pmetr=='time-reqd', pmetr=pmetr)
    return tmpfcn

# ----------------------------------------------------------------------------------------
def get_pdirfcn(locus):
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, 'cache-parameters', locus=locus)
    return tmpfcn

# ----------------------------------------------------------------------------------------
def get_ptntypes(method):
    ptypes = partition_types
    # if is_single_chain(method):
    #     ptypes = [t for t in ptypes if t != 'joint']  # this is probably needlessly general
    return ptypes

# ----------------------------------------------------------------------------------------
def skip_this(pmetr, ptntype, method, ltmp):
    if 'pcfrac-' in pmetr and (ptntype != 'joint' or ltmp != 'igh'):  # only plot pair info cleaning fractions for joint ptntype
        return True
    if pmetr == 'naive-hdist' and ptntype != 'single':  # only do annotation performance for single chain (at least for now)
        return True
    if pmetr == 'time-reqd' and ptntype == 'joint' and is_single_chain(method): # or method=='vsearch-partition')):
        return True
    if ptntype != 'joint' and is_joint_method(method):
        return True
    if args.bcrham_time and ptntype == 'joint':
        return True
    if method in tree_perf_actions and ptntype == 'single':
        return True
    return False

# ----------------------------------------------------------------------------------------
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

for action in args.actions:
    if action == 'parse-linearham-trees':
        parse_linearham_trees()
    elif action in ['simu'] + after_actions + plot_actions:
        run_scan(action)
    elif action in ['plot', 'combine-plots']:
        if args.dry:
            print('  --dry: not plotting')
            continue
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars['partition'])
        if action == 'plot':
            print('plotting %d combinations of %d variable%s (%s) to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), scanplot.get_comparison_plotdir(args, None)))
            fnames = {meth : {pmetr : [[] for _ in partition_types] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            procs = []
            for method in args.plot_metrics:  # NOTE in cf-tree-metrics.py these are [selection] metrics, but here they're [clustering] methods
                for pmetr in args.perf_metrics:
                    utils.prep_dir(scanplot.get_comparison_plotdir(args, method) + '/' + pmetr, wildlings=['*.html', '*.svg', '*.yaml'])  # , subdirs=args.perf_metrics
                for ipt, ptntype in enumerate(get_ptntypes(method)):
                    for ltmp in plot_loci(methods=[method]):
                        for pmetr in args.perf_metrics:
                            if skip_this(pmetr, ptntype, method, ltmp):
                                continue
                            print('  %12s  %6s partition: %3s %s' % (method, ptntype.replace('single', 'single chain'), ltmp, pmetr))
                            arglist, kwargs = (args, args.scan_vars['partition'], action, method, pmetr, args.final_plot_xvar), {'fnfcn' : get_fnfcn(method, ltmp, ptntype, pmetr), 'locus' : ltmp, 'ptntype' : ptntype, 'fnames' : fnames[method][pmetr][ipt], 'pdirfcn' : get_pdirfcn(ltmp), 'script_base' : script_base, 'debug' : args.debug}
                            if args.test:
                                scanplot.make_plots(*arglist, **kwargs)
                            else:
                                procs.append(multiprocessing.Process(target=scanplot.make_plots, args=arglist, kwargs=kwargs))
            if not args.test:
                utils.run_proc_functions(procs)
            # NOTE fnames doesn't get filled by the multiprocessing stuff, so if you want the html you have to set --test
            for method in args.plot_metrics:
                for pmetr in args.perf_metrics:
                    pmcdir = scanplot.get_comparison_plotdir(args, method) + '/' + pmetr
                    fnames[method][pmetr] = [[f.replace(pmcdir+'/', '') for f in flist] for flist in fnames[method][pmetr]]
                    plotting.make_html(pmcdir, n_columns=3, fnames=fnames[method][pmetr])  # this doesn't work unless --test is set since multiprocessing uses copies of <fnames>, but whatever, just run combine-plots
        elif action == 'combine-plots':
            cfpdir = scanplot.get_comparison_plotdir(args, 'combined')
            utils.prep_dir(cfpdir, wildlings=['*.html', '*.svg'])
            fnames, iplot = [[] for _ in args.perf_metrics], 0
            for ipm, pmetr in enumerate([m for m in args.perf_metrics if 'pcfrac-correct-family' not in m]):  # see note in read_hist_csv()
                print('    ', pmetr)
                for ptntype in partition_types:
                    for ltmp in plot_loci(methods=args.plot_metrics):
                        if skip_this(pmetr, ptntype, args.plot_metrics[0], ltmp):  # kind of sucks to use [0], but if i'm skipping one i'm probably skipping all of them
                            continue
                        scanplot.make_plots(args, args.scan_vars['partition'], action, None, pmetr, args.final_plot_xvar, locus=ltmp, ptntype=ptntype, fnames=fnames[int(ipm//3) if 'pcfrac-' in pmetr else ipm], make_legend=ltmp in ['igh', None], leg_label='-'+ptntype, script_base=script_base, debug=args.debug)
                        # iplot += 1
            fnames += [[os.path.dirname(fnames[0][0]) + '/legend-%s.svg'%ptntype] for ptntype in partition_types]
            plotting.make_html(cfpdir, n_columns=3 if len(plot_loci())==3 else 4, fnames=fnames)  # NOTE the pcfrac ones have to be first in the list for the ipm/3 thing to work
        else:
            raise Exception('unsupported action %s' % action)
    else:
        raise Exception('unsupported action %s' % action)

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
