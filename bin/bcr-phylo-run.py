#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import csv
import colored_traceback.always
import collections
import copy
import os
import sys
import numpy
import math
import time
import traceback
import itertools
import yaml
from io import open

import partis.utils as utils
import partis.indelutils as indelutils
import partis.treeutils as treeutils
from partis.event import RecombinationEvent
import partis.paircluster as paircluster

bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'
ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def simdir():
    return '%s/selection/simu' % args.base_outdir
def infdir():
    return '%s/selection/partis' % args.base_outdir
def evtdir(i, igcr=None):
    return '%s/event-%d%s' % (simdir(), i, '' if igcr is None else '/round-%d'%igcr)
def spath(tstr, unsampled=False):  # use spath() for building command line args, whereas use get_simfn() to get [an] inidivudal file e.g. for utils.output_exists(), as well as for passing to fcns in paircluster.py
    if args.mutated_outpath and tstr == 'mutated':
        opth = args.base_outdir
    else:
        opth = '%s/%s-simu' % (simdir(), tstr)
    return '%s%s%s' % (opth, '-unsampled' if (args.tpsample and unsampled) else '', '' if args.paired_loci else '.yaml')
def sfname(tstr, ltmp, lpair=None, unsampled=False):
    if ltmp is None: assert not args.paired_loci
    return paircluster.paired_fn(spath(tstr, unsampled=unsampled), ltmp, lpair=lpair, suffix='.yaml') if args.paired_loci else spath(tstr, unsampled=unsampled)
def naive_fname(ltmp, lpair=None):  # args are only used for paired loci (but we pass the whole fcn to another fcn, so we need the signature like this)
    return sfname('naive', ltmp, lpair=lpair) #paircluster.paired_fn(spath('naive'), ltmp, lpair=lpair, suffix='.yaml') if args.paired_loci else spath('naive')
def bcr_phylo_fasta_fname(outdir):
    return '%s/%s.fasta' % (outdir, args.extrastr)
def get_simfn(ltmp, lpair=None, joint=False):  # NOTE joint has no effect, but is needed for passing to paircluster.write_concatd_output_files()
    return sfname('mutated', ltmp, lpair=lpair)
def get_unsampled_simfn(ltmp, lpair=None, joint=False):  # ugly to have this, but signature has to be this so we can pass it to fcns in paircluster
    return sfname('mutated', ltmp, lpair=lpair, unsampled=True)
# ----------------------------------------------------------------------------------------
def ipath(stype):  # path/file name for building command line args
    rpath = infdir()
    if args.paired_loci:
        return rpath
    assert stype in ['params', 'partition', 'plots']
    rpath = '%s/%s' % (rpath, stype)
    if stype == 'partition':
        rpath += '.yaml'
    return rpath
# ----------------------------------------------------------------------------------------
def ifname(stype, ltmp='igh'):  # path/files for utils.output_exists()
    rpath = ipath(stype)
    if args.paired_loci:
        if stype == 'partition':
            rpath = paircluster.paired_fn(rpath, ltmp, suffix='.yaml', actstr=stype)
        else:
            rpath += '/parameters/%s' % ltmp
    if stype == 'params':
        rpath += '/hmm/hmms'
    return rpath
# ----------------------------------------------------------------------------------------
def lpairs():
    return utils.locus_pairs[ig_or_tr]
# ----------------------------------------------------------------------------------------
def rearrange():
    if utils.output_exists(args, naive_fname('igh'), outlabel='naive simu', offset=4):  # just look for the merged igh file, since it's about the last to be written (and both paired subdirs may not be there)
        return
    cmd = './bin/partis simulate --simulate-from-scratch --mutation-multiplier 0.0001 --n-leaves 1 --constant-number-of-leaves'  # tends to get in infinite loop if you actually pass 0. (yes, I should fix this)
    cmd += ' --debug %d --random-seed %d --n-sim-events %d' % (int(args.debug), args.seed, args.n_sim_events if not args.restrict_to_single_naive_seq else 1)
    if args.paired_loci:
        cmd += ' --paired-loci --paired-outdir %s' % spath('naive')
    else:
        cmd += ' --outfname %s' % spath('naive')
    if args.restrict_available_genes:
        assert not args.paired_loci
        cmd += ' --only-genes IGHV1-18*01:IGHJ1*01'
    if args.rearr_extra_args is not None:
        cmd += ' %s' % args.rearr_extra_args
    if args.single_light_locus is not None:
        cmd += ' --single-light-locus %s' % args.single_light_locus
    if args.n_procs > 1:
        cmd += ' --n-procs %d' % args.n_procs
    if args.slurm:
        cmd += ' --batch-system slurm'
    utils.simplerun(cmd, dryrun=args.dry_run, debug=True)

# ----------------------------------------------------------------------------------------
def get_vpar_val(parg, pval, debug=False):  # get value of parameter/command line arg that is allowed to (but may not at the moment) be drawn from a variable distribution (note we have to pass in <pval> for args that are lists)
    if args.parameter_variances is None or parg not in args.parameter_variances:  # default: just use the single, fixed value from the command line
        return pval
    if args.n_gc_rounds is not None and parg in ['obs-times', 'n-sim-seqs-per-generation']:
        raise Exception('shouldn\'t get here (see exception elsewhere)')
    def sfcn(x):  # just for dbg/exceptions
        return str(int(x)) if parg != 'selection-strength' else ('%.2f' % x)
    pvar = args.parameter_variances[parg]
    if '..' in pvar:  # list of allowed values NOTE pval is *not* used if we're choosing from several choices (ick, but not sure what else to do)
        dbgstr = '[%s]' % pvar.replace('..', ', ')
        return_val = numpy.random.choice([float(pv) for pv in pvar.split('..')])
    else:  # actual parameter variance (i know, this is ugly)
        parg_bounds = {'min' : {'n-sim-seqs-per-generation' : 1}, 'max' : {}}
        pmean = pval
        pvar = float(pvar)
        pmin, pmax = pmean - 0.5 * pvar, pmean + 0.5 * pvar
        if pmin < 0:
            raise Exception('min parameter value for %s less than 0 (from mean %s and half width %s)' % (parg, sfcn(pmean), sfcn(pvar)))
        if parg == 'selection-strength' and pmax > 1:
            raise Exception('max parameter value for %s greater than 1 (from mean %s and half width %s)' % (parg, sfcn(pmean), sfcn(pvar)))
        if parg in parg_bounds['min'] and pmin < parg_bounds['min'][parg]:
            raise Exception('min value too small for %s: %f < %f' % (parg, pmin, parg_bounds['min'][parg]))
        if parg in parg_bounds['max'] and pmax > parg_bounds['max'][parg]:
            raise Exception('max value too large for %s: %f > %f' % (parg, pmax, parg_bounds['max'][parg]))
        dbgstr = '[%6s, %6s]' % (sfcn(pmin), sfcn(pmax))
        return_val = numpy.random.uniform(pmin, pmax)
    if parg != 'selection-strength':
        return_val = int(return_val)
    if debug:
        print('  %30s --> %-7s  %s' % (dbgstr, sfcn(return_val), parg))
    return return_val

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(naive_seq, outdir, ievent, uid_str_len=None, igcr=None):
    if utils.output_exists(args, bcr_phylo_fasta_fname(outdir), outlabel='bcr-phylo', offset=4):
        return None

    cmd = '%s/bin/simulator.py' % bcr_phylo_path
    if args.run_help:
        cmd += ' --help'
    else:
        cmd += ' --lambda0 %f' % args.base_mutation_rate
        if args.no_selection:
            cmd += ' --no_selection'
        else:
            cmd += ' --selection_strength %f' % get_vpar_val('selection-strength', args.selection_strength)
        for astr in ['obs-times', 'n-sim-seqs-per-generation']:  # for search: obs_times n_sim_seqs_per_generation
            aval = getattr(args, astr.replace('-', '_'))
            tstr = ' '.join('%d' % get_vpar_val(astr, t) for t in (aval if args.n_gc_rounds is None else aval[igcr]))
            cmd += ' --%s %s' % (astr.replace('n-sim-seqs-per-generation', 'n-to-sample').replace('-', '_'), tstr)  # ick
        cmd += ' --metric_for_target_dist %s' % args.metric_for_target_distance
        if args.multifurcating_tree:
            cmd += ' --multifurcating_tree'
        if args.aa_paratope_positions is not None:
            cmd += ' --aa_paratope_positions %s' % args.aa_paratope_positions
        if args.aa_struct_positions is not None:
            cmd += ' --aa_struct_positions %s' % args.aa_struct_positions
        if args.dont_mutate_struct_positions:
            cmd += ' --dont_mutate_struct_positions'
        if args.skip_stops:
            cmd += ' --skip_stops_when_mutating'
        if args.allow_stops:
            cmd += ' --allow_stops_in_functional_seqs'
        cmd += ' --target_dist %d' % args.target_distance
        cmd += ' --target_count %d' % args.target_count
        cmd += ' --carry_cap %d' % get_vpar_val('carry-cap', args.carry_cap)
        if not args.dont_observe_common_ancestors:
            cmd += ' --observe_common_ancestors'
        if args.leaf_sampling_scheme is not None:
            cmd += ' --leaf_sampling_scheme %s' % args.leaf_sampling_scheme
        if args.n_target_clusters is not None:
            cmd += ' --n_target_clusters %d' % args.n_target_clusters
        # cmd += ' --target_cluster_distance 1'
        if args.tdist_scale is not None:
            cmd += ' --tdist_scale %d' % args.tdist_scale
        if args.tdist_weights is not None:
            cmd += ' --tdist_weights %s' % args.tdist_weights
        if args.min_target_distance is not None:
            cmd += ' --min_target_distance %d' % args.min_target_distance
        if args.min_effective_kd is not None:
            cmd += ' --min_effective_kd %d' % args.min_effective_kd
        if args.n_naive_seq_copies is not None:
            cmd += ' --n_naive_seq_copies %d' % args.n_naive_seq_copies
        if args.n_gc_rounds is not None and igcr > 0:
            init_fn = '%s/init-seqs.fa' % outdir
            if not args.dry_run:
                isfos = utils.read_fastx(bcr_phylo_fasta_fname(evtdir(ievent, igcr=igcr - 1)))
                if args.n_reentry_seqs is not None:
                    if args.n_reentry_seqs > len(isfos):
                        print('  %s --n-reentry-seqs %d greater than number of observed seqs %d in %s' % (utils.wrnstr(), args.n_reentry_seqs, len(isfos), bcr_phylo_fasta_fname(evtdir(ievent, igcr=igcr - 1))))
                    isfos = numpy.random.choice(isfos, size=args.n_reentry_seqs, replace=False)
                utils.write_fasta(init_fn, isfos)
            cmd += ' --initial_seq_file %s' % init_fn

    cmd += ' --debug %d' % args.debug
    cmd += ' --n_tries 1000'
    if args.context_depend == 0:
        cmd += ' --no_context'
    cmd += ' --no_plot'
    if args.only_csv_plots:
        cmd += ' --dont_write_hists'
    cmd += ' --outbase %s/%s' % (outdir, args.extrastr)
    cmd += ' --random_seed %d' % (args.seed + ievent)  # NOTE if args.n_gc_rounds is set, it's *really* important that this is the same for each round since it ensures we have the same target sequence
    if uid_str_len is not None:
        cmd += ' --uid_str_len %d' % uid_str_len
    cmd += ' --naive_seq %s' % naive_seq

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cfo = None
    if args.n_procs == 1:
        utils.run_ete_script(cmd, dryrun=args.dry_run)
    else:
        cmd = utils.run_ete_script(cmd, return_for_cmdfos=True, dryrun=args.dry_run)
        cfo = {'cmd_str' : cmd, 'workdir' : outdir, 'outfname' : bcr_phylo_fasta_fname(outdir)}
    sys.stdout.flush()
    return cfo

# ----------------------------------------------------------------------------------------
def parse_bcr_phylo_output(glfos, naive_events, outdir, ievent, uid_info):
    # ----------------------------------------------------------------------------------------
    def split_seqfos(seqfos):
        hline, lline = naive_events[ievent]
        hseqfos, lseqfos = [], []
        for sfo in seqfos:
            padseq = utils.pad_nuc_seq(hline['naive_seq'])
            assert len(sfo['seq']) == len(padseq) + len(lline['naive_seq'])
            hseqfos.append({'name' : sfo['name'], 'seq' : sfo['seq'][ : len(hline['naive_seq'])]})
            lseqfos.append({'name' : sfo['name'], 'seq' : sfo['seq'][len(padseq) : ]})
        return hseqfos, lseqfos
    # ----------------------------------------------------------------------------------------
    def read_kdvals(kdfname):
        nodefo = {}
        with open(kdfname) as kdfile:
            reader = csv.DictReader(kdfile)
            for line in reader:
                nodefo[line['uid']] = {
                    'kd' : float(line['kd']),
                    'relative_kd' : float(line['relative_kd']),
                    'lambda' : float(line['lambda']) if line['lambda'] != '' else None,  # bcr-phylo used to not run the lambda update fcn after last iteratio, which resulted in empty lambda values, but it shouldn't happen any more (but leaving here for backwards compatibility)
                    'target_index' : int(line['target_index']),
                    'target_distance' : float(line['target_distance']),
                    'time' : int(line['time']),
                }
        return nodefo
    # ----------------------------------------------------------------------------------------
    def get_mature_line(sfos, naive_line, glfo, nodefo, dtree, target_sfos, dup_translations=None, locus=None):
        assert len(naive_line['unique_ids']) == 1  # enforces that we ran naive-only, 1-leaf partis simulation above
        assert not indelutils.has_indels(naive_line['indelfos'][0])  # would have to handle this below
        if args.debug:
            utils.print_reco_event(naive_line)
        reco_info = collections.OrderedDict()
        for sfo in sfos:
            mline = utils.get_non_implicit_copy(naive_line)
            del mline['tree']
            mline['unique_ids'] = [sfo['name']]
            mline['seqs'] = [sfo['seq']]
            mline['input_seqs'] = [sfo['seq']]  # it's really important to set both the seqs (since they're both already in there from the naive line)
            mline['duplicates'] = [[]]
            reco_info[sfo['name']] = mline
            try:
                utils.add_implicit_info(glfo, mline)
            except:  # TODO not sure if I really want to leave this in long term, but it shouldn't hurt anything (it's crashing on unequal naive/mature sequence lengths, and I need this to track down which event it is) UPDATE: yeah it was just because something crashed in the middle of writing a .fa file
                print('implicit info adding failed for ievent %d in %s' % (ievent, outdir))
                lines = traceback.format_exception(*sys.exc_info())
                print(utils.pad_lines(''.join(lines)))  # NOTE this will still crash on the next line if implicit info adding failed
        final_line = utils.synthesize_multi_seq_line_from_reco_info([sfo['name'] for sfo in sfos], reco_info)

        ftree = copy.deepcopy(dtree)
        if locus is not None:
            def ltr(u): return u + '-' + locus
            new_nodefo = {}
            for u_old in nodefo:
                new_nodefo[ltr(u_old)] = nodefo[u_old]
            nodefo = new_nodefo
            treeutils.translate_labels(ftree, [(u, ltr(u)) for u in final_line['unique_ids']])
            final_line['unique_ids'] = [ltr(u) for u in final_line['unique_ids']]
            assert len(sfos) == len(final_line['unique_ids'])
            for iseq, sfo in enumerate(sfos):
                naive_id = naive_line['unique_ids'][0]
                assert naive_id.count('-') == 1
                bstr = naive_id.replace('-'+locus, '')
                pids = final_line['paired-uids'][iseq]
                assert len(pids) == 1 and pids[0].find(bstr) == 0 and pids[0].count('-') == 1 and pids[0].split('-')[1] in utils.loci  # if uid is xxx-igh, paired id shoud be e.g. xxx-igk
                final_line['paired-uids'][iseq] = [p.replace(bstr, sfo['name']) for p in pids]

        tmp_trns = {}
        for iu, old_id in enumerate(final_line['unique_ids']):  # NOTE this only translates the uids, for paired h/l we still need to go back through and translate paired uids
            if old_id in uid_info['all_uids']:
                new_id, uid_info['n_duplicate_uids'] = utils.choose_non_dup_id(old_id, uid_info['n_duplicate_uids'], uid_info['all_uids'])
                tmp_trns[old_id] = new_id
                final_line['unique_ids'][iu] = new_id
            uid_info['all_uids'].add(final_line['unique_ids'][iu])
        if len(tmp_trns) > 0:
            for old_id, new_id in tmp_trns.items():
                nodefo[new_id] = nodefo[old_id]
                del nodefo[old_id]
            treeutils.translate_labels(ftree, [(o, n) for o, n in tmp_trns.items()], expect_missing=True)
        if dup_translations is not None:
            dup_translations.update(tmp_trns)

        if len(set(nodefo) - set(final_line['unique_ids'])) > 0:  # uids in the kd file but not the <line> (i.e. not in the newick/fasta files) are probably just bcr-phylo discarding internal nodes
            print('        in kd file, but missing from final_line (probably just internal nodes that bcr-phylo wrote to the tree without names): %s' % (set(nodefo) - set(final_line['unique_ids'])))
        if len(set(final_line['unique_ids']) - set(nodefo)) > 0:
            print('        in final_line, but missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(nodefo)))
        final_line['affinities'] = [1. / nodefo[u]['kd'] for u in final_line['unique_ids']]
        if args.affinity_measurement_error is not None:
            # before = final_line['affinities']
            final_line['affinities'] = [numpy.random.normal(a, args.affinity_measurement_error * a) for a in final_line['affinities']]
            # print '  '.join('%.4f'%v for v in before)
            # print '  '.join('%.4f'%v for v in final_line['affinities'])
        final_line['relative_affinities'] = [1. / nodefo[u]['relative_kd'] for u in final_line['unique_ids']]
        final_line['lambdas'] = [nodefo[u]['lambda'] for u in final_line['unique_ids']]
        final_line['nearest_target_indices'] = [nodefo[u]['target_index'] for u in final_line['unique_ids']]
        final_line['min_target_distances'] = [nodefo[u]['target_distance'] for u in final_line['unique_ids']]
        final_line['generation-times'] = [nodefo[u]['time'] for u in final_line['unique_ids']]
        ftree.scale_edges(1. / numpy.mean([len(s) for s in final_line['seqs']]))  # note that if --paired-loci is set then most edges will still be the wrong length (compared to the mutations in the single-locus sequences), i.e. best not to use this much until treeutils.combine_selection_metrics(), where we rescale to the full h+l length
        # treeutils.compare_tree_distance_to_shm(ftree, final_line, debug=True)
        if args.debug:
            print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ftree), padwidth=12))
        final_line['tree'] = ftree.as_string(schema='newick')
        if args.debug:
            utils.print_reco_event(final_line) #, extra_print_keys=['lambdas'])

        tmp_event = RecombinationEvent(glfo)  # I don't want to move the function out of event.py right now
        tmp_event.set_reco_id(final_line, irandom=ievent)  # not sure that setting <irandom> here actually does anything
        final_line['target_seqs'] = [[tfo['seq'] for tfo in target_sfos] for _ in final_line['unique_ids']]  # NOTE it would be nice to put this in a per-family key, but it ends up we want it to behave like an input meta info key, and input meta keys need to be per-seq since they're inherently a property of each sequence. So instead we duplicate this info across all seqs to which it applies
        return final_line
    # ----------------------------------------------------------------------------------------
    def translate_duplicate_pids(mpair, dup_translations):
        if len(dup_translations) == 0:
            return
        assert len(set(len(l['unique_ids']) for l in mpair)) == 1  # make sure h and l annotations have the same length
        for atn1, atn2 in itertools.permutations(mpair, 2):
            # print ':'.join([utils.color('red' if atn1['paired-uids'][i]!=[u] else None, u) for i, u in enumerate(atn2['unique_ids'])])
            for pids, uid in zip(atn1['paired-uids'], atn2['unique_ids']):  # this is just to double check things, so could be removed
                assert len(pids) == 1
                if pids[0] != uid:
                    assert pids[0] in dup_translations and dup_translations[pids[0]] == uid
                    del dup_translations[pids[0]]
            atn1['paired-uids'] = [[u] for u in atn2['unique_ids']]  # seems a bit hackey to reset all of them, not just the translated one, but whatever

    # ----------------------------------------------------------------------------------------
    # extract kd values from pickle file (used to need a separate script since ete3 needed python 3, but now probably doesn't)
    kdfname, nwkfname = '%s/kd-vals.csv' % outdir, '%s/simu.nwk' % outdir
    if not utils.output_exists(args, kdfname, outlabel='kd/nwk conversion', offset=4):  # eh, don't really need to check for both kd and nwk file, chances of only one being missing are really small, and it'll just crash when it looks for it a couple lines later
        cmd = './bin/read-bcr-phylo-trees.py --pickle-tree-file %s/%s_lineage_tree.p --kdfile %s --newick-tree-file %s' % (outdir, args.extrastr, kdfname, nwkfname)
        utils.run_ete_script(cmd, debug=args.n_procs==1)
    nodefo = read_kdvals(kdfname)
    dtree = treeutils.get_dendro_tree(treefname=nwkfname)
    seqfos = utils.read_fastx(bcr_phylo_fasta_fname(outdir))  # output mutated sequences from bcr-phylo
    target_seqfos = utils.read_fastx('%s/%s_targets.fa' % (outdir, args.extrastr))
    if args.paired_loci:
        mpair = []
        dup_translations = {}
        for tline, sfos, tsfos in zip(naive_events[ievent], split_seqfos(seqfos), split_seqfos(target_seqfos)):
            mpair.append(get_mature_line(sfos, tline, glfos[tline['loci'][0]], nodefo, dtree, target_seqfos, dup_translations=dup_translations, locus=tline['loci'][0]))
        translate_duplicate_pids(mpair, dup_translations)
        return mpair
    else:
        return get_mature_line(seqfos, naive_events[ievent], glfos[0], nodefo, dtree, target_seqfos)

# ----------------------------------------------------------------------------------------
def read_rearrangements():
    if args.paired_loci:
        lp_infos = paircluster.read_lpair_output_files(lpairs(), naive_fname, dbgstr='naive simulation')
        naive_events = paircluster.get_all_antn_pairs(lp_infos)
        glfos, _, _ = paircluster.concat_heavy_chain(lpairs(), lp_infos)  # per-locus glfos with concat'd heavy chain
    else:
        glfo, naive_events, _ = utils.read_output(naive_fname(None))
        glfos = [glfo]
    return glfos, naive_events

# ----------------------------------------------------------------------------------------
def write_simulation(glfos, mutated_events, unsampled=False):
    opath = spath('mutated', unsampled=unsampled)
    print('  writing%s annotations to %s' % ('' if not args.tpsample else (' unsampled' if unsampled else ' timepoint sampled'), opath))
    mheads = []
    if args.n_gc_rounds is not None:
        mheads += ['gc-rounds', 'generation-times']
    if args.tpsample:
        mheads += ['timepoints']
    headers = utils.simulation_headers + mheads
    if args.paired_loci:
        lp_infos = {}
        for lpair in lpairs():
            lpfos = {k : {} for k in ['glfos', 'antn_lists', 'cpaths']}  # cpaths i think don't get used
            mevents = [(hl, ll) for hl, ll in mutated_events if [hl['loci'][0], ll['loci'][0]] == lpair]  # grab the events for this h/l pair
            for ltmp, levents in zip(lpair, zip(*mevents)):
                lpfos['antn_lists'][ltmp] = levents
                lpfos['glfos'][ltmp] = glfos[ltmp]
            lp_infos[tuple(lpair)] = lpfos
        sfcn = get_unsampled_simfn if (args.tpsample and unsampled) else get_simfn
        paircluster.write_lpair_output_files(lpairs(), lp_infos, sfcn, headers)
        glfos, antn_lists, _ = paircluster.concat_heavy_chain(lpairs(), lp_infos)  # per-locus glfos with concat'd heavy chain
        paircluster.write_concatd_output_files(glfos, antn_lists, sfcn, headers)
        outfos, metafos = paircluster.get_combined_outmetafos(antn_lists, extra_meta_headers=mheads)
        paircluster.write_combined_fasta_and_meta(opath+'/all-seqs.fa', opath+'/meta.yaml', outfos, metafos)
    else:
        utils.write_annotations(opath, glfos[0], mutated_events, headers)

# ----------------------------------------------------------------------------------------
def combine_gc_rounds(glfos, mevt_lists):
    # ----------------------------------------------------------------------------------------
    def fix_evt(igcr, sum_time, evt, all_gtimes, ltmp=None):
        evt['generation-times'] = [t + sum_time for t in evt['generation-times']]
        evt['gc-rounds'] = [igcr for _ in evt['unique_ids']]
        if ltmp is None:
            def tfcn(u, t): return '%s-%s'%(u, t)
        else:
            def tfcn(u, t): lstr = u.split('-')[-1]; assert lstr in utils.loci; return u.replace('-'+lstr, '-%s-%s'%(t, lstr))
        trns = {u : tfcn(u, t) for u, t in zip(evt['unique_ids'], evt['generation-times'])}
        if ltmp is not None:
            trns.update({p : tfcn(p, t) for pids, t in zip(evt['paired-uids'], evt['generation-times']) for p in pids})
        utils.translate_uids([evt], trns=trns, translate_pids=args.paired_loci)  # kind of annoying to add the timepoint to the uid, but otherwise we get duplicate uids in different rounds (and we can't change the random seed atm, or else the target seqs will be out of whack)
        all_gtimes |= set(evt['generation-times'])
    # ----------------------------------------------------------------------------------------
    if utils.output_exists(args, get_simfn('igh'), outlabel='mutated simu', offset=4):
        return None
    assert len(mevt_lists) == args.n_gc_rounds
    assert len(set(len(l) for l in mevt_lists)) == 1  # all rounds should have the same number of events
    sum_time = 0
    assert len(args.obs_times) == args.n_gc_rounds  # also checked elsewhere
    all_gtimes, stlist = set(), []
    for igcr in range(args.n_gc_rounds):
        if args.paired_loci:
            for epair in mevt_lists[igcr]:
                for evt in epair:
                    fix_evt(igcr, sum_time, evt, all_gtimes, ltmp=evt['loci'][0])
        else:
            for evt in mevt_lists[igcr]:
                fix_evt(igcr, sum_time, evt, all_gtimes)
        sum_time += args.obs_times[igcr][-1]
        stlist.append(sum_time)
    print('    merging %d events over %d gc rounds with final generation times%s: %s' % (len(mevt_lists[0]), args.n_gc_rounds, '' if args.dont_observe_common_ancestors else ' (including observed common ancestors)', ' '.join(utils.color('blue' if t in stlist else None, str(t)) for t in sorted(all_gtimes))))
    merged_events = []
    for ievt in range(len(mevt_lists[0])):
        if args.paired_loci:
            mpair = []
            lpair = [l['loci'][0] for l in mevt_lists[0][0]]
            for ilocus, ltmp in enumerate(lpair):
                mgevt = utils.combine_events(glfos[ltmp], [evts[ievt][ilocus] for evts in mevt_lists], meta_keys=['gc-rounds', 'generation-times'])
                mpair.append(mgevt)
            merged_events.append(mpair)
        else:
            mgevt = utils.combine_events(glfos[0], [evts[ievt] for evts in mevt_lists], meta_keys=['gc-rounds', 'generation-times'])
            merged_events.append(mgevt)

    write_simulation(glfos, merged_events, unsampled=args.tpsample)

    return merged_events

# ----------------------------------------------------------------------------------------
def sample_tp_seqs(glfos, evt_list, l_evts=None, ltmp=None):
    id_list = [u for l in evt_list for u in l['unique_ids']]
    if len(set(id_list)) != len(id_list):
        raise Exception('duplicate ids in final events')  # shouldn't be able to get to here, but if it does it'll break stuff below
    for fevt in evt_list:
        if 'timepoints' in fevt:  # shouldn't be in there
            print('  %s \'timepoints\' already in event (overwriting)' % utils.wrnstr())
        fevt['timepoints'] = [None for _ in fevt['unique_ids']]
    all_gtimes = set(t for l in evt_list for t in l['generation-times'])
    gt_ids = {t : [] for t in all_gtimes}  # map from each generation time to list of all remaining uids with that time
    for tline in evt_list:
        for tid, gtime in zip(tline['unique_ids'], tline['generation-times']):
            gt_ids[gtime].append(tid)
    print('                       N      generation     N          N')
    print('         timepoint   total      time       chosen   remaining')
    for tpfo in args.sequence_sample_times:
        if any(t not in all_gtimes for t in tpfo['times']):
            raise Exception('generation time %s not among actual final times: %s' % (' '.join(str(t) for t in tpfo['times'] if t not in all_gtimes), ' '.join(str(t) for t in sorted(all_gtimes))))
        allowed_uids = [u for gt in tpfo['times'] for u in gt_ids[gt]]
        if tpfo['n'] > len(allowed_uids):
            print('  %s not enough allowed seqs remain (%d > %d, probably didn\'t sample enough sequences at allowed times %s)' % (utils.wrnstr(), tpfo['n'], len(allowed_uids), ' '.join(str(t) for t in tpfo['times'])))
        chosen_ids = numpy.random.choice(allowed_uids, size=tpfo['n'], replace=False)
        if len(chosen_ids) != tpfo['n']:
            print('  %s couldn\'t choose enough seqs (only got %d)' % (utils.wrnstr(), len(chosen_ids)))
        n_chosen = {}
        for gtime in tpfo['times']:
            n_before = len(gt_ids[gtime])
            gt_ids[gtime] = [u for u in gt_ids[gtime] if u not in chosen_ids]
            n_chosen[gtime] = n_before - len(gt_ids[gtime])
        for igt, gtime in enumerate(tpfo['times']):
            print('    %12s      %3s      %3d         %4d       %4d' % (tpfo['name'] if igt==0 else '', '%d'%tpfo['n'] if igt==0 else '', gtime, n_chosen[gtime], len(gt_ids[gtime])))
        for fevt in evt_list:
            fevt['timepoints'] = [tpfo['name'] if u in chosen_ids else t for u, t in zip(fevt['unique_ids'], fevt['timepoints'])]
    for ievt, fevt in enumerate(evt_list):
        iseqs_to_keep = [i for i, t in enumerate(fevt['timepoints']) if t is not None]
        if l_evts is not None:
            hevt, levt = fevt, l_evts[ievt]
            hloc, lloc = [e['loci'][0] for e in [hevt, levt]]
            def htrans(u, hloc, lloc): return u.replace('-'+hloc, '-'+lloc)
            assert [htrans(u, hloc, lloc) for u in hevt['unique_ids']] == levt['unique_ids']
            levt['timepoints'] = [t for t in hevt['timepoints']]  # NOTE *has* to happen before restrict_to_iseqs() (duh)
            utils.restrict_to_iseqs(levt, iseqs_to_keep, glfos[levt['loci'][0]], remove_tree=True)
        utils.restrict_to_iseqs(fevt, iseqs_to_keep, glfos[0 if ltmp is None else ltmp], remove_tree=True)
    # utils.print_reco_event(levt, extra_print_keys=['timepoints', 'gc-rounds', 'generation-times'])

# ----------------------------------------------------------------------------------------
def write_timepoint_sampled_sequences(glfos, final_events):
    if utils.output_exists(args, get_simfn('igh'), outlabel='mutated simu', offset=4):
        return None
    if args.paired_loci:
        h_evts, l_evts = list(zip(*final_events))
        sample_tp_seqs(glfos, h_evts, l_evts=l_evts, ltmp=h_evts[0]['loci'][0])
    else:
        sample_tp_seqs(glfos, final_events)

    write_simulation(glfos, final_events)

# ----------------------------------------------------------------------------------------
def simulate(igcr=None):

    if igcr in [None, 0]:
        rearrange()

    glfos, naive_events = read_rearrangements()
    if args.dry_run:
        for ievent in range(args.n_sim_events):
            _ = run_bcr_phylo('<NAIVE_SEQ>', evtdir(ievent, igcr=igcr), ievent, igcr=igcr)
        return None, None
    if args.restrict_to_single_naive_seq:
        print('  --restrict-to-single-naive-seq: using same naive event for all mutation simulations')
        assert len(naive_events) == 1
        naive_events = [naive_events[0] for _ in range(args.n_sim_events)]
    else:
        assert len(naive_events) == args.n_sim_events

    outdirs = [evtdir(i, igcr=igcr) for i in range(len(naive_events))]

    start = time.time()
    cmdfos = []
    if args.n_procs > 1:
        print('    starting %d events%s' % (len(naive_events), '' if args.n_procs==1 else ' with N max simultaneous procs %d'%args.n_procs))
    uid_str_len = args.min_ustr_len  # UPDATE don't need to increase this any more since I'm handling duplicates when above + int(math.log(len(naive_events), 7))  # if the final sample's going to contain many trees, it's worth making the uids longer so there's fewer collisions/duplicates (note that this starts getting pretty slow if it's bigger than 7 or so)
    for ievent, outdir in enumerate(outdirs):
        if args.n_sim_events > 1 and args.n_procs == 1:
            print('  %s %d' % (utils.color('blue', 'ievent'), ievent))
        if args.paired_loci:
            hnseq, lnseq = [l['naive_seq'] for l in naive_events[ievent]]
            naive_seq = utils.pad_nuc_seq(hnseq) + lnseq
        else:
            naive_seq = naive_events[ievent]['naive_seq']
        cfo = run_bcr_phylo(naive_seq, outdir, ievent, uid_str_len=uid_str_len, igcr=igcr)  # if n_procs > 1, doesn't run, just returns cfo
        if cfo is not None:
            print('      %s %s' % (utils.color('red', 'run'), cfo['cmd_str']))
            cmdfos.append(cfo)
    if args.n_procs > 1 and len(cmdfos) > 0:
        utils.run_cmds(cmdfos, shell=True, n_max_procs=args.n_procs, batch_system='slurm' if args.slurm else None, allow_failure=True, debug='print')
    print('  bcr-phylo run time: %.1fs' % (time.time() - start))

    if utils.output_exists(args, get_simfn('igh'), outlabel='mutated simu', offset=4):  # i guess if it crashes during the plotting just below, this'll get confused
        return None, None

    start = time.time()
    uid_info = {'all_uids' : set(), 'n_duplicate_uids' : 0}  # stuff just for dealing with duplicate uids
    mutated_events = []
    for ievent, outdir in enumerate(outdirs):
        mutated_events.append(parse_bcr_phylo_output(glfos, naive_events, outdir, ievent, uid_info))
    if uid_info['n_duplicate_uids'] > 0:
        print('  %s renamed %d duplicate uids from %d bcr-phylo events' % (utils.color('yellow', 'warning'), uid_info['n_duplicate_uids'], len(mutated_events)))
    print('  parsing time: %.1fs' % (time.time() - start))

    if igcr is None:
        write_simulation(glfos, mutated_events, unsampled=args.tpsample)

    if not args.only_csv_plots:
        import partis.lbplotting as lbplotting
        for ievent, outdir in enumerate(outdirs):
            if args.paired_loci:
                lpair = [l['loci'][0] for l in mutated_events[ievent]]
                evtlist = mutated_events[ievent]
            else:
                lpair = None
                evtlist = [mutated_events[ievent]]
            lbplotting.plot_bcr_phylo_simulation(outdir + '/plots', outdir, evtlist, args.extrastr, lbplotting.metric_for_target_distance_labels[args.metric_for_target_distance], lpair=lpair)
        # utils.simplerun('cp -v %s/simu_collapsed_runstat_color_tree.svg %s/plots/' % (outdir, outdir))

    return glfos, mutated_events

# ----------------------------------------------------------------------------------------
def simulseq_args():
    cstr = ''
    if args.restrict_to_single_naive_seq:
        print('  note: using --all-seqs-simultaneous because --restrict-to-single-naive-seq was set')
        cstr += ' --all-seqs-simultaneous'
    if args.n_gc_rounds is None and not args.tpsample:
        cstr += ' --is-simu'
        if '--all-seqs-simultaneous' not in cstr:
            cstr += ' --simultaneous-true-clonal-seqs'
    elif args.n_sim_events == 1:
        print('  %s not using --is-simu since --n-gc-rounds or --sequence-sample-time-fname are set, so e.g. plots won\'t use true info, and true tree won\'t be set' % utils.wrnstr())
        if '--all-seqs-simultaneous' not in cstr:
            cstr += ' --all-seqs-simultaneous'
    else:
        print('  %s not using any of --is-simu or --simultaneous-true-clonal-seqs since either --n-gc-rounds or --sequence-sample-time-fname are set with more than one event, so e.g. plots won\'t use true info, and true tree won\'t be set' % utils.wrnstr())
    return cstr

# ----------------------------------------------------------------------------------------
def cache_parameters():
    if utils.output_exists(args, ifname('params'), outlabel='parameters', offset=4):
        return
    cmd = './bin/partis cache-parameters --random-seed %d --no-indels' % args.seed  # forbid indels because in the very rare cases when we call them, they're always wrong, and then they screw up the simultaneous true clonal seqs option
    cmd += simulseq_args()
    fstr = ' --paired-loci --paired-indir %s --paired-outdir %s' if args.paired_loci else ' --infname %s --parameter-dir %s'
    cmd += fstr % (spath('mutated'), ipath('params'))
    if args.all_inference_plots:
        cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else ipath('plots'))
    if args.meta_info_key_to_color is not None:
        cmd += ' --meta-info-key-to-color %s' % args.meta_info_key_to_color
    if args.inf_extra_args is not None:
        cmd += ' %s' % args.inf_extra_args
    if args.n_procs > 1:
        cmd += ' --n-procs %d' % args.n_procs
    if args.slurm:
        cmd += ' --batch-system slurm'
    if args.n_max_queries is not None:
        cmd += ' --n-max-queries %d' % args.n_max_queries
    utils.simplerun(cmd, debug=True, dryrun=args.dry_run)
    sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def partition():
    if utils.output_exists(args, ifname('partition'), outlabel='partition', offset=4):
        return
    cmd = './bin/partis partition --random-seed %d' % args.seed
    cmd += simulseq_args()
    fstr = ' --paired-loci --paired-indir %s --paired-outdir %s' if args.paired_loci else (' --infname %%s --parameter-dir %s --outfname %%s' % ipath('params'))
    cmd += fstr % (spath('mutated'), ipath('partition'))
    #  --write-additional-cluster-annotations 0:5  # I don't think there was really a good reason for having this
    if not args.dont_get_tree_metrics:
        cmd += ' --get-selection-metrics'
    if args.tree_inference_method is not None:
        cmd += ' --tree-inference-method %s' % args.tree_inference_method
    if not args.dont_get_tree_metrics or args.all_inference_plots:
        cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else ipath('plots'))
    if not args.all_inference_plots:
        cmd += ' --no-partition-plots'
    if args.meta_info_key_to_color is not None:
        cmd += ' --meta-info-key-to-color %s' % args.meta_info_key_to_color
    if args.inf_extra_args is not None:
        cmd += ' %s' % args.inf_extra_args
    if args.lb_tau is not None:
        cmd += ' --lb-tau %f' % args.lb_tau
    if args.n_procs > 1:
        cmd += ' --n-procs %d' % args.n_procs
    if args.slurm:
        cmd += ' --batch-system slurm'
    if args.n_max_queries is not None:
        cmd += ' --n-max-queries %d' % args.n_max_queries
    if args.extra_smetric_plots is not None:
        cmd += ' --selection-metric-plot-cfg %s' % ':'.join(treeutils.default_plot_cfg + args.extra_smetric_plots + ['distr'])
    utils.simplerun(cmd, debug=True, dryrun=args.dry_run)
    # cmd = './bin/partis get-selection-metrics --outfname %s/partition.yaml' % infdir()
    # utils.simplerun(cmd, debug=True) #, dryrun=True)
    sys.stdout.flush()

# ----------------------------------------------------------------------------------------
all_actions = ('simu', 'cache-parameters', 'partition')
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter)
parser.add_argument('--actions', default=':'.join(all_actions), help='which actions to run')
parser.add_argument('--base-outdir', default='%s/partis/bcr-phylo/test' % os.getenv('fs', default=os.getenv('HOME')), help='base output dir')
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--run-help', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--only-csv-plots', action='store_true')
parser.add_argument('--dont-get-tree-metrics', action='store_true', help='Partition without getting tree metrics, presumably because you want to run them yourself later')
parser.add_argument('--tree-inference-method')
parser.add_argument('--seed', type=int, default=1, help='random seed (note that bcr-phylo doesn\'t seem to support setting its random seed)')
parser.add_argument('--n-procs', type=int, default=1)
parser.add_argument('--extrastr', default='simu', help='doesn\'t really do anything, but it\'s required by bcr-phylo EDIT ok it actually doesn\'t need it, it\'s just that the output files look weird without it because they start with \'_\' if it\'s empty')
parser.add_argument('--n-sim-seqs-per-generation', default='100', help='Number of sequences to sample at each time in --obs-times.')
parser.add_argument('--n-sim-events', type=int, default=1, help='number of simulated rearrangement events')
parser.add_argument('--n-max-queries', type=int, help='during parameter caching and partitioning, stop after reading this many queries from simulation file (useful for dtr training samples where we need massive samples to actually train the dtr, but for testing various metrics need far smaller samples).')
parser.add_argument('--obs-times', default='100:120', help='Times (generations) at which to select sequences for observation. Note that this is the time that the sequences existed in/exited the gc, not necessarily the time at which we "sequenced" them. If --n-gc-rounds is set, this must be a colon-separated list of comma-separated lists (see help under that arg).')
parser.add_argument('--sequence-sample-time-fname', help='Times at which we "sample" for sequencing, i.e. at which we draw blood (as opposed to --obs-times, which is the generation time at which a cell leaves the gc). Specified in a yaml file as a list of key : val pairs, with key the timepoint label and values a dict with keys \'n\' (total number of sequences) and \'times\' (dict keyed by gc round time (if --n-gc-rounds is set, otherwise just a list) of generation times from which to sample those \'n\' sequences uniformly at random). If set, a new output file/dir is created by inserting \'timepoint-sampled\' at end of regular output file. See example in data/sample-seqs.yaml.')
parser.add_argument('--n-naive-seq-copies', type=int, help='see bcr-phylo docs')
parser.add_argument('--n-gc-rounds', type=int, help='number of rounds of gc entry, i.e. if set, upon gc completion we choose --n-reentry-seqs sampled seqs with which to seed a new (otherwise identical) gc reaction. Results for each gc round N are written to subdirs round-N/ for each event, then all sampled sequences from all reactions are collected into the normal output file locations, with input meta info key \'gc-rounds\' specifying their gc round. If this arg is set, then --obs-times must be a colon-separated list (of length --n-gc-rounds) of comma-separated lists, where each sample time is relative to the *start* of that round.')
parser.add_argument('--n-reentry-seqs', type=int, help='number of sampled seqs from previous round (chosen randomly) to inject into the next gc round (if not set, we take all of them).')
parser.add_argument('--carry-cap', type=int, default=1000, help='carrying capacity of germinal center')
parser.add_argument('--target-distance', type=int, default=15, help='Desired distance (number of non-synonymous mutations) between the naive sequence and the target sequences.')
parser.add_argument('--tdist-scale', type=int, help='see bcr-phylo docs')
parser.add_argument('--tdist-weights', help='see bcr-phylo docs')
parser.add_argument('--metric-for-target-distance', default='aa', choices=['aa', 'nuc', 'aa-sim-ascii', 'aa-sim-blosum'], help='see bcr-phylo docs')
parser.add_argument('--target-count', type=int, default=1, help='Number of target sequences to generate.')
parser.add_argument('--n-target-clusters', type=int, help='number of cluster into which to divide the --target-count target seqs (see bcr-phylo docs)')
parser.add_argument('--min-target-distance', type=int, help='see bcr-phylo docs')
parser.add_argument('--min-effective-kd', type=float, help='see bcr-phylo docs')
parser.add_argument('--affinity-measurement-error', type=float, help='Fractional measurement error on affinity: if set, replace \'affinities\' in final partis line with new values smeared with a normal distribution with this fractional width, i.e. <a> is replaced with a value drawn from a normal distribution with mean <a> and width <f>*<a> for this fraction <f>.')
parser.add_argument('--base-mutation-rate', type=float, default=0.365, help='see bcr-phylo docs')
parser.add_argument('--selection-strength', type=float, default=1., help='see bcr-phylo docs')
parser.add_argument('--context-depend', type=int, default=0, choices=[0, 1])  # i wish this could be a boolean, but having it int makes it much much easier to interface with the scan infrastructure in cf-tree-metrics.py
parser.add_argument('--aa-paratope-positions', help='see bcr-phylo docs')
parser.add_argument('--aa-struct-positions', help='see bcr-phylo docs')
parser.add_argument('--dont-mutate-struct-positions', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--skip-stops', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--allow-stops', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--no-selection', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--multifurcating-tree', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--restrict-available-genes', action='store_true', help='restrict v and j gene choice to one each (so context dependence is easier to plot)')
parser.add_argument('--restrict-to-single-naive-seq', action='store_true', help='restrict all events to use the same naive sequence')
parser.add_argument('--lb-tau', type=float, help='')
parser.add_argument('--dont-observe-common-ancestors', action='store_true')
parser.add_argument('--leaf-sampling-scheme', help='see bcr-phylo help')
parser.add_argument('--parameter-variances', help='if set, parameters vary from family to family in one of two ways 1) the specified parameters are drawn from a uniform distribution of the specified width (with mean from the regular argument) for each family. Format example: n-sim-seqs-per-generation,10:carry-cap,150 would give --n-sim-seqs-per-generation +/-5 and --carry-cap +/-75, or 2) parameters for each family are chosen from a \'..\'-separated list, e.g. obs-times,75..100..150')
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--paired-loci', action='store_true')
parser.add_argument('--parameter-plots', action='store_true', help='DEPRECATED')
parser.add_argument('--all-inference-plots', action='store_true')
parser.add_argument('--meta-info-key-to-color')
parser.add_argument('--single-light-locus', help='set to igk or igl if you want only that one; otherwise each event is chosen at random (see partis help)')
parser.add_argument('--rearr-extra-args', help='')
parser.add_argument('--inf-extra-args', help='')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--mutated-outpath', action='store_true', help='write final (mutated) output file[s] to --base-outdir, rather than the default of burying them in subdirs with intermediate files')
parser.add_argument('--extra-smetric-plots', default=':'.join(treeutils.default_plot_cfg))
parser.add_argument('--min-ustr-len', type=int, default=5, help='min length of hashed uid strs (longer makes collisions less likely, but it\'s hard to avoid them entirely since independent bcr-phylo procs choose the uids for each family)')

args = parser.parse_args()

if args.parameter_plots:
    print('  %s transferring deprecated arg --parameter-plots to --all-inference-plots' % utils.wrnstr())
    args.all_inference_plots = True
    delattr(args, 'parameter_plots')
if args.seed is not None:
    numpy.random.seed(args.seed)
args.obs_times = utils.get_arg_list(args.obs_times, intify=True, list_of_lists=args.n_gc_rounds is not None)
args.n_sim_seqs_per_generation = utils.get_arg_list(args.n_sim_seqs_per_generation, intify=True, list_of_lists=args.n_gc_rounds is not None)
args.actions = utils.get_arg_list(args.actions, choices=all_actions)
args.parameter_variances = utils.get_arg_list(args.parameter_variances, key_val_pairs=True, choices=['selection-strength', 'obs-times', 'n-sim-seqs-per-generation', 'carry-cap', 'metric-for-target-distance'])  # if you add more, make sure the bounds enforcement and conversion stuff in get_vpar_val() are still ok
args.extra_smetric_plots = utils.get_arg_list(args.extra_smetric_plots, choices=treeutils.all_plot_cfg)
if args.rearr_extra_args is not None:
    args.rearr_extra_args = args.rearr_extra_args.replace('@', ' ')  # ick this sucks
if args.inf_extra_args is not None:
    args.inf_extra_args = args.inf_extra_args.replace('@', ' ')  # ick this sucks
if args.affinity_measurement_error is not None:
    assert args.affinity_measurement_error >= 0
    if args.affinity_measurement_error > 1:
        print('  note: --affinity-measurement-error %.2f is greater than 1 -- this is fine as long as it\'s on purpose, but will result in smearing by a normal with width larger than each affinity value (and probably result in some negative values).' % args.affinity_measurement_error)
if args.n_gc_rounds is not None:
    assert len(args.obs_times) == args.n_gc_rounds
    for otlist in args.obs_times:
        if otlist != sorted(otlist):  # various things assume it's sorted
            raise Exception('obs times within each gc round must be sorted')
    otstrs = ['%s' % ' '.join(str(t) for t in otlist) for otlist in args.obs_times]
    def fgt(i, t): return t + sum(args.obs_times[j][-1] for j in range(i))
    fgstrs = ['%s' % ' '.join(str(fgt(i, t)) for t in otlist) for i, otlist in enumerate(args.obs_times)]
    print('    --obs-times at each of %d gc rounds: %s  (final generation times: %s)' % (args.n_gc_rounds, ', '.join(otstrs), ', '.join(fgstrs)))
    if len(args.n_sim_seqs_per_generation) != args.n_gc_rounds and len(args.n_sim_seqs_per_generation) == 1:
        args.n_sim_seqs_per_generation = [args.n_sim_seqs_per_generation[0] for _ in range(args.n_gc_rounds)]
    assert len(args.n_sim_seqs_per_generation) == args.n_gc_rounds
    if args.parameter_variances is not None:  # don't feel like implementing this atm
        if any(a in args.parameter_variances for a in ['obs-times', 'n-sim-seqs-per-generation']):
            raise Exception('haven\'t implemented parameter variances for --obs-times/--n-sim-seqs-per-generation with multiple gc rounds')
setattr(args, 'sequence_sample_times', None)
setattr(args, 'tpsample', False)
if args.sequence_sample_time_fname is not None:
    print('  reading timepoint sample times from %s' % args.sequence_sample_time_fname)
    with open(args.sequence_sample_time_fname) as sfile:
        yamlfo = yaml.load(sfile, Loader=yaml.CLoader)
    if args.n_gc_rounds is not None:  # have to translate the "local" gc round times to final times, then also concatenate them into ones list for all gc rounds
        sum_time = 0
        for igcr in range(args.n_gc_rounds):
            for tpfo in yamlfo:
                if igcr not in tpfo['times']:
                    continue
                tpfo['times'][igcr] = [t + sum_time for t in tpfo['times'][igcr]]
            sum_time += args.obs_times[igcr][-1]
        for tpfo in yamlfo:
            tpfo['times'] = sorted(t for tlist in tpfo['times'].values() for t in tlist)  # should already be sorted, but whatever
    args.sequence_sample_times = yamlfo
    args.tpsample = True  # just a shorthand
    delattr(args, 'sequence_sample_time_fname')

assert args.extrastr == 'simu'  # I think at this point this actually can't be changed without changing some other things

# ----------------------------------------------------------------------------------------
if 'simu' in args.actions:
    if args.n_gc_rounds is None:
        glfos, final_events = simulate()
    else:
        mevt_lists = []  # list (for each gc round) of [sub]lists, where each sublist is the sampled seqs from that round for each event
        for igcr in range(args.n_gc_rounds):
            glfos, mevts = simulate(igcr=igcr)
            mevt_lists.append(mevts)
        if not args.dry_run:
            final_events = combine_gc_rounds(glfos, mevt_lists)
    if not args.dry_run and args.tpsample:
        write_timepoint_sampled_sequences(glfos, final_events)
if 'cache-parameters' in args.actions:
    cache_parameters()
if 'partition' in args.actions:
    partition()
