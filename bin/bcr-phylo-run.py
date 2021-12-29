#!/usr/bin/env python
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

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
import utils
import indelutils
import treeutils
from event import RecombinationEvent
import paircluster

ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'
ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def simdir():
    return '%s/selection/simu' % args.base_outdir
def infdir():
    return '%s/selection/partis' % args.base_outdir
def evtdir(i):
    return '%s/event-%d' % (simdir(), i)
def spath(tstr):  # use spath() for building command line args, whereas use simfname() to get [an] inidivudal file e.g. for utils.output_exists(), as well as for passing to fcns in paircluster.py
    if args.mutated_outpath and tstr == 'mutated':
        opth = args.base_outdir
    else:
        opth = '%s/%s-simu' % (simdir(), tstr)
    return '%s%s' % (opth, '' if args.paired_loci else '.yaml')
def sfname(tstr, ltmp, lpair=None):
    if ltmp is None: assert not args.paired_loci
    return paircluster.paired_fn(spath(tstr), ltmp, lpair=lpair, suffix='.yaml') if args.paired_loci else spath(tstr)
def naive_fname(ltmp, lpair=None):  # args are only used for paired loci (but we pass the whole fcn to another fcn, so we need the signature like this)
    return sfname('naive', ltmp, lpair=lpair) #paircluster.paired_fn(spath('naive'), ltmp, lpair=lpair, suffix='.yaml') if args.paired_loci else spath('naive')
def bcr_phylo_fasta_fname(outdir):
    return '%s/%s.fasta' % (outdir, args.extrastr)
def simfname(ltmp, lpair=None, joint=False):  # NOTE joint has no effect, but is needed for passing to paircluster.write_concatd_output_files()
    return sfname('mutated', ltmp, lpair=lpair)
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
    cmd += ' --debug %d --random-seed %d --n-sim-events %d' % (int(args.debug), args.seed, args.n_sim_events)
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
        print '  %30s --> %-7s  %s' % (dbgstr, sfcn(return_val), parg)
    return return_val

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(naive_seq, outdir, ievent, uid_str_len=None):
    if utils.output_exists(args, bcr_phylo_fasta_fname(outdir), outlabel='bcr-phylo', offset=4):
        return None

    cmd = '%s/bin/simulator.py' % bcr_phylo_path
    if args.run_help:
        cmd += ' --help'
    else:
        cmd += ' --lambda0 %f' % args.base_mutation_rate
        cmd += ' --selection_strength %f' % get_vpar_val('selection-strength', args.selection_strength)
        cmd += ' --obs_times %s' % ' '.join(['%d' % get_vpar_val('obs-times', t) for t in args.obs_times])
        cmd += ' --n_to_sample %s' % ' '.join('%d' % get_vpar_val('n-sim-seqs-per-generation', n) for n in args.n_sim_seqs_per_generation)
        cmd += ' --metric_for_target_dist %s' % args.metric_for_target_distance
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
        if args.no_selection:
            cmd += ' --no_selection'
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
        if args.n_initial_seqs is not None:
            cmd += ' --n_initial_seqs %d' % args.n_initial_seqs

    cmd += ' --debug %d' % args.debug
    cmd += ' --n_tries 1000'
    if args.context_depend == 0:
        cmd += ' --no_context'
    cmd += ' --no_plot'
    if args.only_csv_plots:
        cmd += ' --dont_write_hists'
    cmd += ' --outbase %s/%s' % (outdir, args.extrastr)
    cmd += ' --random_seed %d' % (args.seed + ievent)
    if uid_str_len is not None:
        cmd += ' --uid_str_len %d' % uid_str_len
    cmd += ' --naive_seq %s' % naive_seq

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cfo = None
    if args.n_procs == 1:
        utils.run_ete_script(cmd, ete_path, dryrun=args.dry_run)
    else:
        cmd, _ = utils.run_ete_script(cmd, ete_path, return_for_cmdfos=True, tmpdir=outdir, dryrun=args.dry_run)
        cfo = {'cmd_str' : cmd, 'workdir' : outdir, 'outfname' : bcr_phylo_fasta_fname(outdir)}
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
                }
        return nodefo
    # ----------------------------------------------------------------------------------------
    def get_mature_line(sfos, naive_line, glfo, nodefo, dtree, target_sfos, locus=None):
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
                print 'implicit info adding failed for ievent %d in %s' % (ievent, outdir)
                lines = traceback.format_exception(*sys.exc_info())
                print utils.pad_lines(''.join(lines))  # NOTE this will still crash on the next line if implicit info adding failed
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

        dup_translations = {}
        for iu, old_id in enumerate(final_line['unique_ids']):
            if old_id in uid_info['all_uids']:
                new_id, uid_info['n_duplicate_uids'] = utils.choose_non_dup_id(old_id, uid_info['n_duplicate_uids'], uid_info['all_uids'])
                dup_translations[old_id] = new_id
                final_line['unique_ids'][iu] = new_id
            uid_info['all_uids'].add(final_line['unique_ids'][iu])
        if len(dup_translations) > 0:
            for old_id, new_id in dup_translations.items():
                nodefo[new_id] = nodefo[old_id]
                del nodefo[old_id]
            treeutils.translate_labels(ftree, [(o, n) for o, n in dup_translations.items()])

        # extract kd values from pickle file (use a separate script since it requires ete/anaconda to read)
        if len(set(nodefo) - set(final_line['unique_ids'])) > 0:  # uids in the kd file but not the <line> (i.e. not in the newick/fasta files) are probably just bcr-phylo discarding internal nodes
            print '        in kd file, but missing from final_line (probably just internal nodes that bcr-phylo wrote to the tree without names): %s' % (set(nodefo) - set(final_line['unique_ids']))
        if len(set(final_line['unique_ids']) - set(nodefo)) > 0:
            print '        in final_line, but missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(nodefo))
        final_line['affinities'] = [1. / nodefo[u]['kd'] for u in final_line['unique_ids']]
        final_line['relative_affinities'] = [1. / nodefo[u]['relative_kd'] for u in final_line['unique_ids']]
        final_line['lambdas'] = [nodefo[u]['lambda'] for u in final_line['unique_ids']]
        final_line['nearest_target_indices'] = [nodefo[u]['target_index'] for u in final_line['unique_ids']]
        final_line['min_target_distances'] = [nodefo[u]['target_distance'] for u in final_line['unique_ids']]
        ftree.scale_edges(1. / numpy.mean([len(s) for s in final_line['seqs']]))
        if args.debug:
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ftree), padwidth=12)
        final_line['tree'] = ftree.as_string(schema='newick')
        if args.debug:
            utils.print_reco_event(final_line) #, extra_print_keys=['lambdas'])

        tmp_event = RecombinationEvent(glfo)  # I don't want to move the function out of event.py right now
        tmp_event.set_reco_id(final_line, irandom=ievent)  # not sure that setting <irandom> here actually does anything
        final_line['target_seqs'] = [tfo['seq'] for tfo in target_sfos]
        return final_line

    # ----------------------------------------------------------------------------------------
    kdfname, nwkfname = '%s/kd-vals.csv' % outdir, '%s/simu.nwk' % outdir
    if not utils.output_exists(args, kdfname, outlabel='kd/nwk conversion', offset=4):  # eh, don't really need to check for both kd and nwk file, chances of only one being missing are really small, and it'll just crash when it looks for it a couple lines later
        cmd = './bin/read-bcr-phylo-trees.py --pickle-tree-file %s/%s_lineage_tree.p --kdfile %s --newick-tree-file %s' % (outdir, args.extrastr, kdfname, nwkfname)
        utils.run_ete_script(cmd, ete_path, debug=args.n_procs==1)
    nodefo = read_kdvals(kdfname)
    dtree = treeutils.get_dendro_tree(treefname=nwkfname)
    seqfos = utils.read_fastx(bcr_phylo_fasta_fname(outdir))  # output mutated sequences from bcr-phylo
    target_seqfos = utils.read_fastx('%s/%s_targets.fa' % (outdir, args.extrastr))
    if args.paired_loci:
        mevents = []
        for tline, sfos, tsfos in zip(naive_events[ievent], split_seqfos(seqfos), split_seqfos(target_seqfos)):
            mevents.append(get_mature_line(sfos, tline, glfos[tline['loci'][0]], nodefo, dtree, target_seqfos, locus=tline['loci'][0]))
        return mevents
    else:
        return get_mature_line(seqfos, naive_events[ievent], glfos[0], nodefo, dtree, target_seqfos)

# ----------------------------------------------------------------------------------------
def read_rearrangements():
    if args.paired_loci:
        lp_infos = paircluster.read_lpair_output_files(lpairs(), naive_fname, dbgstr='naive simulation')
        naive_events = paircluster.get_both_lpair_antn_pairs(lpairs(), lp_infos)
        glfos, _, _ = paircluster.concat_heavy_chain(lpairs(), lp_infos)  # per-locus glfos with concat'd heavy chain
    else:
        glfo, naive_events, _ = utils.read_output(naive_fname(None))
        glfos = [glfo]
    return glfos, naive_events

# ----------------------------------------------------------------------------------------
def write_simulation(glfos, mutated_events):
    if args.paired_loci:
        lp_infos = {}
        for lpair in lpairs():
            lpfos = {k : {} for k in ['glfos', 'antn_lists', 'cpaths']}  # cpaths i think don't get used
            mevents = [(hl, ll) for hl, ll in mutated_events if [hl['loci'][0], ll['loci'][0]] == lpair]  # grab the events for this h/l pair
            for ltmp, levents in zip(lpair, zip(*mevents)):
                lpfos['antn_lists'][ltmp] = levents
                lpfos['glfos'][ltmp] = glfos[ltmp]
            lp_infos[tuple(lpair)] = lpfos
        paircluster.write_lpair_output_files(lpairs(), lp_infos, simfname, utils.simulation_headers)
        glfos, antn_lists, _ = paircluster.concat_heavy_chain(lpairs(), lp_infos)  # per-locus glfos with concat'd heavy chain
        paircluster.write_concatd_output_files(glfos, antn_lists, simfname, utils.simulation_headers)
        paircluster.write_merged_simu(antn_lists, spath('mutated')+'/all-seqs.fa', spath('mutated')+'/meta.yaml')
    else:
        utils.write_annotations(spath('mutated'), glfos[0], mutated_events, utils.simulation_headers)

# ----------------------------------------------------------------------------------------
def simulate():

    rearrange()

    glfos, naive_events = read_rearrangements()
    if args.dry_run:
        for ievent in range(args.n_sim_events):
            _ = run_bcr_phylo('<NAIVE_SEQ>', evtdir(ievent), ievent)
        return
    assert len(naive_events) == args.n_sim_events

    outdirs = [evtdir(i) for i in range(len(naive_events))]

    start = time.time()
    cmdfos = []
    if args.n_procs > 1:
        print '    starting %d events' % len(naive_events)
    uid_str_len = args.min_ustr_len  # UPDATE don't need to increase this any more since I'm handling duplicates when above + int(math.log(len(naive_events), 7))  # if the final sample's going to contain many trees, it's worth making the uids longer so there's fewer collisions/duplicates (note that this starts getting pretty slow if it's bigger than 7 or so)
    for ievent, outdir in enumerate(outdirs):
        if args.n_sim_events > 1 and args.n_procs == 1:
            print '  %s %d' % (utils.color('blue', 'ievent'), ievent)
        if args.paired_loci:
            hnseq, lnseq = [l['naive_seq'] for l in naive_events[ievent]]
            naive_seq = utils.pad_nuc_seq(hnseq) + lnseq
        else:
            naive_seq = naive_events[ievent]['naive_seq']
        cfo = run_bcr_phylo(naive_seq, outdir, ievent, uid_str_len=uid_str_len)  # if n_procs > 1, doesn't run, just returns cfo
        if cfo is not None:
            print '      %s %s' % (utils.color('red', 'run'), cfo['cmd_str'])
            cmdfos.append(cfo)
    if args.n_procs > 1 and len(cmdfos) > 0:
        utils.run_cmds(cmdfos, shell=True, n_max_procs=args.n_procs, batch_system='slurm' if args.slurm else None, allow_failure=True, debug='print')
    print '  bcr-phylo run time: %.1fs' % (time.time() - start)

    if utils.output_exists(args, simfname('igh'), outlabel='mutated simu', offset=4):  # i guess if it crashes during the plotting just below, this'll get confused
        return

    start = time.time()
    uid_info = {'all_uids' : set(), 'n_duplicate_uids' : 0}  # stuff just for dealing with duplicate uids
    mutated_events = []
    for ievent, outdir in enumerate(outdirs):
        mutated_events.append(parse_bcr_phylo_output(glfos, naive_events, outdir, ievent, uid_info))
    if uid_info['n_duplicate_uids'] > 0:
        print '  %s renamed %d duplicate uids from %d bcr-phylo events' % (utils.color('yellow', 'warning'), uid_info['n_duplicate_uids'], len(mutated_events))
    print '  parsing time: %.1fs' % (time.time() - start)

    print '  writing annotations to %s' % spath('mutated')
    write_simulation(glfos, mutated_events)

    if not args.only_csv_plots:
        import lbplotting
        for ievent, outdir in enumerate(outdirs):
            if args.paired_loci:
                lpair = [l['loci'][0] for l in mutated_events[ievent]]
                evtlist = mutated_events[ievent]
            else:
                lpair = None
                evtlist = [mutated_events[ievent]]
            lbplotting.plot_bcr_phylo_simulation(outdir + '/plots', outdir, evtlist, args.extrastr, lbplotting.metric_for_target_distance_labels[args.metric_for_target_distance], lpair=lpair)
        # utils.simplerun('cp -v %s/simu_collapsed_runstat_color_tree.svg %s/plots/' % (outdir, outdir))

# ----------------------------------------------------------------------------------------
def cache_parameters():
    if utils.output_exists(args, ifname('params'), outlabel='parameters', offset=4):
        return
    cmd = './bin/partis cache-parameters --random-seed %d --no-indels' % args.seed  # forbid indels because in the very rare cases when we call them, they're always wrong, and then they screw up the simultaneous true clonal seqs option
    fstr = ' --paired-loci --paired-indir %s --paired-outdir %s' if args.paired_loci else ' --infname %s --parameter-dir %s'
    cmd += fstr % (spath('mutated'), ipath('params'))
    if args.parameter_plots:
        cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else ipath('plots'))
    if args.inf_extra_args is not None:
        cmd += ' %s' % args.inf_extra_args
    if args.n_procs > 1:
        cmd += ' --n-procs %d' % args.n_procs
    if args.slurm:
        cmd += ' --batch-system slurm'
    if args.n_max_queries is not None:
        cmd += ' --n-max-queries %d' % args.n_max_queries
    utils.simplerun(cmd, debug=True, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def partition():
    if utils.output_exists(args, ifname('partition'), outlabel='partition', offset=4):
        return
    cmd = './bin/partis partition --simultaneous-true-clonal-seqs --is-simu --random-seed %d' % args.seed
    fstr = ' --paired-loci --paired-indir %s --paired-outdir %s' if args.paired_loci else (' --infname %%s --parameter-dir %s --outfname %%s' % ipath('params'))
    cmd += fstr % (spath('mutated'), ipath('partition'))
    #  --write-additional-cluster-annotations 0:5  # I don't think there was really a good reason for having this
    if not args.dont_get_tree_metrics:
        cmd += ' --get-selection-metrics --plotdir %s --no-partition-plots' % ('paired-outdir' if args.paired_loci else ipath('plots'))
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
    if args.distr_hists:
        cmd += ' --selection-metric-plot-cfg %s' % ':'.join(treeutils.default_plot_cfg + ['distr'])
    utils.simplerun(cmd, debug=True, dryrun=args.dry_run)
    # cmd = './bin/partis get-selection-metrics --outfname %s/partition.yaml' % infdir()
    # utils.simplerun(cmd, debug=True) #, dryrun=True)

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
parser.add_argument('--seed', type=int, default=1, help='random seed (note that bcr-phylo doesn\'t seem to support setting its random seed)')
parser.add_argument('--n-procs', type=int, default=1)
parser.add_argument('--extrastr', default='simu', help='doesn\'t really do anything, but it\'s required by bcr-phylo EDIT ok it actually doesn\'t need it, it\'s just that the output files look weird without it because they start with \'_\' if it\'s empty')
parser.add_argument('--n-sim-seqs-per-generation', default='100', help='Number of sequences to sample at each time in --obs-times.')
parser.add_argument('--n-sim-events', type=int, default=1, help='number of simulated rearrangement events')
parser.add_argument('--n-max-queries', type=int, help='during parameter caching and partitioning, stop after reading this many queries from simulation file (useful for dtr training samples where we need massive samples to actually train the dtr, but for testing various metrics need far smaller samples).')
parser.add_argument('--obs-times', default='100:120', help='Times (reproductive rounds) at which to selection sequences for observation.')
parser.add_argument('--carry-cap', type=int, default=1000, help='carrying capacity of germinal center')
parser.add_argument('--target-distance', type=int, default=15, help='Desired distance (number of non-synonymous mutations) between the naive sequence and the target sequences.')
parser.add_argument('--tdist-scale', type=int, help='see bcr-phylo docs')
parser.add_argument('--tdist-weights', help='see bcr-phylo docs')
parser.add_argument('--metric-for-target-distance', default='aa', choices=['aa', 'nuc', 'aa-sim-ascii', 'aa-sim-blosum'], help='see bcr-phylo docs')
parser.add_argument('--target-count', type=int, default=1, help='Number of target sequences to generate.')
parser.add_argument('--n-target-clusters', type=int, help='number of cluster into which to divide the --target-count target seqs (see bcr-phylo docs)')
parser.add_argument('--min-target-distance', type=int, help='see bcr-phylo docs')
parser.add_argument('--min-effective-kd', type=float, help='see bcr-phylo docs')
parser.add_argument('--n-initial-seqs', type=int, help='see bcr-phylo docs')
parser.add_argument('--base-mutation-rate', type=float, default=0.365, help='see bcr-phylo docs')
parser.add_argument('--selection-strength', type=float, default=1., help='see bcr-phylo docs')
parser.add_argument('--context-depend', type=int, default=0, choices=[0, 1])  # i wish this could be a boolean, but having it int makes it much much easier to interface with the scan infrastructure in cf-tree-metrics.py
parser.add_argument('--aa-paratope-positions', help='see bcr-phylo docs')
parser.add_argument('--aa-struct-positions', help='see bcr-phylo docs')
parser.add_argument('--dont-mutate-struct-positions', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--skip-stops', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--allow-stops', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--no-selection', action='store_true', help='see bcr-phylo docs')
parser.add_argument('--restrict-available-genes', action='store_true', help='restrict v and j gene choice to one each (so context dependence is easier to plot)')
parser.add_argument('--lb-tau', type=float, help='')
parser.add_argument('--dont-observe-common-ancestors', action='store_true')
parser.add_argument('--leaf-sampling-scheme', help='see bcr-phylo help')
parser.add_argument('--parameter-variances', help='if set, parameters vary from family to family in one of two ways 1) the specified parameters are drawn from a uniform distribution of the specified width (with mean from the regular argument) for each family. Format example: n-sim-seqs-per-generation,10:carry-cap,150 would give --n-sim-seqs-per-generation +/-5 and --carry-cap +/-75, or 2) parameters for each family are chosen from a \'..\'-separated list, e.g. obs-times,75..100..150')
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--paired-loci', action='store_true')
parser.add_argument('--parameter-plots', action='store_true')
parser.add_argument('--single-light-locus', help='set to igk or igl if you want only that one; otherwise each event is chosen at random (see partis help)')
parser.add_argument('--rearr-extra-args', help='')
parser.add_argument('--inf-extra-args', help='')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--mutated-outpath', action='store_true', help='write final (mutated) output file[s] to --base-outdir, rather than the default of burying them in subdirs with intermediate files')
parser.add_argument('--distr-hists', action='store_true', help='include lb distribution hists in plotting')
parser.add_argument('--min-ustr-len', type=int, default=10, help='min length of hashed uid strs (longer makes collisions less likely, but it\'s hard to avoid them entirely since independent bcr-phylo procs choose the uids for each family)')

args = parser.parse_args()

if args.seed is not None:
    numpy.random.seed(args.seed)
args.obs_times = utils.get_arg_list(args.obs_times, intify=True)
args.n_sim_seqs_per_generation = utils.get_arg_list(args.n_sim_seqs_per_generation, intify=True)
args.actions = utils.get_arg_list(args.actions, choices=all_actions)
args.parameter_variances = utils.get_arg_list(args.parameter_variances, key_val_pairs=True, choices=['selection-strength', 'obs-times', 'n-sim-seqs-per-generation', 'carry-cap', 'metric-for-target-distance'])  # if you add more, make sure the bounds enforcement and conversion stuff in get_vpar_val() are still ok

assert args.extrastr == 'simu'  # I think at this point this actually can't be changed without changing some other things

# ----------------------------------------------------------------------------------------
if 'simu' in args.actions:
    simulate()
if 'cache-parameters' in args.actions:
    cache_parameters()
if 'partition' in args.actions:
    partition()
