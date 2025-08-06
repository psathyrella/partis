#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import copy
import time
import random
import sys
import subprocess
# import multiprocessing
import numpy
import scipy
import math
import itertools
import collections
import traceback
from io import open
import colored_traceback.always
import os
import json
import csv
import operator
import glob
from collections import defaultdict
from pathlib import Path
import shutil
partis_dir = str(Path(__file__).parent.parent)
if not os.path.exists(partis_dir):
    print('%s partis dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir)
sys.path.insert(1, partis_dir) # + '/python')

import python.utils as utils
import python.glutils as glutils
import python.treeutils as treeutils
import python.processargs as processargs
import python.seqfileopener as seqfileopener
from python.partitiondriver import PartitionDriver
from python.clusterpath import ClusterPath
import python.paircluster as paircluster
from python.parametercounter import ParameterCounter
from python.corrcounter import CorrCounter
from python.waterer import Waterer

# ----------------------------------------------------------------------------------------
def run_simulation(args):
    # ----------------------------------------------------------------------------------------
    def run_sub_cmds(working_gldir):
        # ----------------------------------------------------------------------------------------
        def get_sub_events(iproc):
            n_sub_events = args.n_sim_events // args.n_procs
            if iproc == args.n_procs - 1 and args.n_sim_events % args.n_procs > 0:  # do any extra ones in the last proc (has to match sub file writer below)
                n_sub_events += args.n_sim_events % args.n_procs
            n_per_proc_list.append(n_sub_events)
            return n_sub_events
        # ----------------------------------------------------------------------------------------
        def get_workdir(iproc):
            return args.workdir + '/sub-' + str(iproc)
        # ----------------------------------------------------------------------------------------
        def get_outfname(iproc):
            return get_workdir(iproc) + '/' + os.path.basename(args.outfname)
        # ----------------------------------------------------------------------------------------
        def write_sub_file(subtype, iproc, n_sub_events, objlist):  # split a list of either trees or heavy chain events into sub lists to send to each sub process
            assert subtype in ['trees', 'heavy-chain-events']
            subfname = '%s/%s-sub-%d.%s'%(args.workdir, subtype, iproc, 'nwk' if subtype=='trees' else 'yaml')
            utils.replace_in_arglist(clist, '--input-simulation-treefname' if subtype=='trees' else '--heavy-chain-event-fname', subfname)
            n_sub_objs = len(objlist) // args.n_procs
            sub_obj_list = objlist[iproc * n_sub_objs : (iproc + 1) * n_sub_objs]
            if iproc == args.n_procs - 1 and len(objlist) % args.n_procs > 0:  # add any extra ones to the last proc (has to match n_sub_events above) NOTE it's important that we add them to the *last* proc (it used to be the first one), since otherwise the h/l paired correlations are wrong, since things end up in a different order after you read files
                sub_obj_list += objlist[args.n_procs * n_sub_objs : ]
            if len(sub_obj_list) == 0:
                raise Exception('couldn\'t split up %d %s among %d procs' % (len(objlist), subtype, args.n_procs))
            if len(sub_obj_list) < n_sub_events:
                print('  note: number of %s %d smaller than number of events %d for sub proc %d (of %d)' % (subtype, len(sub_obj_list), n_sub_events, iproc, args.n_procs))
            work_fnames.append(subfname)
            if subtype == 'trees':
                with open(subfname, 'w') as sfile:
                    sfile.writelines(sub_obj_list)
            else:
                utils.write_annotations(subfname, heavy_chain_glfo, sub_obj_list, utils.simulation_headers)  # there shouldn't be any reason to add the extra headers here, since it's just for input to light chain correlations
        # ----------------------------------------------------------------------------------------
        if args.input_simulation_treefname is not None:  # split up the input trees among the sub procs
            with open(args.input_simulation_treefname) as tfile:
                treelines = tfile.readlines()
            if len(treelines) < args.n_sim_events:
                print('  note: total number of trees %d less than --n-sim-events %d' % (len(treelines), args.n_sim_events))

        cmdfos = []
        for iproc in range(args.n_procs):
            n_sub_events = get_sub_events(iproc)
            clist = copy.deepcopy(sys.argv)
            if shutil.which('partis'):  # use version in PATH if it's there (pipx seems to leave two incompatible versions lying around)
                clist[0] = 'partis'
            utils.replace_in_arglist(clist, '--n-procs', '1')
            clist.append('--im-a-subproc')
            utils.replace_in_arglist(clist, '--random-seed', str(args.random_seed + iproc))
            utils.replace_in_arglist(clist, '--workdir', get_workdir(iproc))
            utils.replace_in_arglist(clist, '--outfname', get_outfname(iproc))
            utils.replace_in_arglist(clist, '--n-sim-events', str(n_sub_events))
            if working_gldir is not None:
                utils.replace_in_arglist(clist, '--initial-germline-dir', working_gldir)
            if args.allele_prevalence_fname is not None:
                utils.replace_in_arglist(clist, '--allele-prevalence-fname', args.allele_prevalence_fname)
            if args.input_simulation_treefname is not None:
                write_sub_file('trees', iproc, n_sub_events, treelines)
            if heavy_chain_events is not None:
                write_sub_file('heavy-chain-events', iproc, n_sub_events, heavy_chain_events)
            cmdstr = ' '.join(clist)
            if args.debug:
                print('  %s %s' % (utils.color('red', 'run'), cmdstr))
            cmdfos.append({'cmd_str' : cmdstr, 'workdir' : get_workdir(iproc), 'logdir' : args.workdir+'/log-'+str(iproc), 'outfname' : get_outfname(iproc)})  # logdirs have to be different than <workdirs> since ./bin/partis (rightfully) barfs if its workdir already exists

        utils.run_cmds(cmdfos, batch_system=args.batch_system, batch_options=args.batch_options, batch_config_fname=args.batch_config_fname, debug='print')
        file_list = [cmdfos[i]['outfname'] for i in range(args.n_procs)]
        utils.merge_simulation_files(args.outfname, file_list, utils.add_lists(list(utils.simulation_headers), args.extra_annotation_columns),  # list() cast is terrible, but somehow I've ended up with some of the headers as lists and some as tuples, and I can't track down all the stuff necessary to synchronize them a.t.m.
                                      n_total_expected=args.n_sim_events, n_per_proc_expected=n_per_proc_list, use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)

        return cmdfos
    # ----------------------------------------------------------------------------------------
    if args.paired_loci:
        run_all_loci(args)
        return

    from python.recombinator import Recombinator

    utils.prep_dir(args.workdir)

    if args.generate_trees:  # note that if we're not simulating from scratch we kind of don't need to generate trees in a separate step (we can just set --choose-trees-in-order and heavy and light will use the same trees in the same order), but for scratch simulation they end up different (I think just from generating the germline set below). It anyway seems safer to really specify the same trees.
        import python.treegenerator as treegenerator
        tmp_shm_parameter_dir = utils.parameter_type_subdir(args, args.shm_parameter_dir) if args.shm_parameter_dir is not None else None
        treegen = treegenerator.TreeGenerator(args, tmp_shm_parameter_dir)
        treegen.generate_trees(args.random_seed, args.outfname, args.workdir)
        return

    default_prevalence_fname = args.workdir + '/allele-prevalence.csv'
    if args.generate_germline_set and args.allele_prevalence_fname is None:
        args.allele_prevalence_fname = default_prevalence_fname

    if args.initial_germline_dir is not None:  # the command line explicitly told us where to get the glfo
        if args.reco_parameter_dir is not None:
            print('  note: getting germline sets from --initial-germline-dir, even though --reco-parameter-dir was also set')
        input_gldir = args.initial_germline_dir
    else:  # if it wasn't explicitly set, we have to decide what to use
        if args.rearrange_from_scratch:  # just use the default
            input_gldir = args.default_initial_germline_dir
        else:  # otherwise assume we're supposed to use the glfo in --reco-parameter-dir
            tpdir = utils.parameter_type_subdir(args, args.reco_parameter_dir)
            if not os.path.exists(tpdir):
                raise Exception('--parameter-dir (with --parameter-type \'%s\') doesn\'t exist: %s' % (utils.get_parameter_type(args, args.reco_parameter_dir), tpdir))
            input_gldir = tpdir + '/' + glutils.glfo_dir
    glfo = glutils.read_glfo(input_gldir, args.locus, only_genes=args.only_genes)

    working_gldir = None  # if we generate a germline set and we're also going to run subprocesses, we need to write that glfo to disk so the subprocs can see it
    if not args.im_a_subproc and args.generate_germline_set:
        glutils.generate_germline_set(glfo, args, new_allele_info=args.new_allele_info, debug=args.debug>1)  # NOTE removes unwanted genes from <glfo>
        if utils.getsuffix(args.outfname) == '.csv':
            print('  writing generated germline set to %s/' % args.outfname.replace('.csv', '-glfo'))
            glutils.write_glfo(args.outfname.replace('.csv', '-glfo'), glfo)
        if args.n_procs > 1:
            working_gldir = args.workdir + '/' + glutils.glfo_dir
            glutils.write_glfo(working_gldir, glfo)

    heavy_chain_events, heavy_chain_glfo = None, None
    if args.paired_correlation_values is not None and not utils.has_d_gene(args.locus):  # if we're doing paired correlations and this is light chain, we need to read the heavy chain annotations so we can pass the proper igh event into each light chain rearrangement
        heavy_chain_glfo, heavy_chain_events, _ = utils.read_output(args.heavy_chain_event_fname, dont_add_implicit_info=True)
        print('    read %d heavy chain events from %s for paired correlations: %s' % (len(heavy_chain_events), args.heavy_chain_event_fname, ' '.join(':'.join(l['unique_ids']) for l in heavy_chain_events)))

    # ----------------------------------------------------------------------------------------
    def make_events(n_events, workdir, outfname, random_ints, heavy_chain_events=None):
        reco = Recombinator(args, glfo, seed=args.random_seed, workdir=workdir, heavy_chain_events=heavy_chain_events)
        start = time.time()
        events = []
        for ievt in range(n_events):
            event = reco.combine(random_ints[ievt], i_choose_tree=ievt if args.choose_trees_in_order else None, i_heavy_event=ievt if heavy_chain_events is not None else None)
            events.append(event)
        if args.check_tree_depths or args.debug:
            reco.print_validation_values()
        utils.write_annotations(outfname, glfo, events, utils.add_lists(list(utils.simulation_headers), args.extra_annotation_columns), synth_single_seqs=utils.getsuffix(outfname) == '.csv', use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)  # keep writing the csv as single-sequence lines, just for backwards compatibility (now trying to switch to synthesizing single seq lines when reading) NOTE list() cast is terrible, but somehow I've ended up with some of the headers as lists and some as tuples, and I can't track down all the stuff necessary to synchronize them a.t.m.
        print('    made %d event%s with %d seqs in %.1fs (%.1fs of which was running bppseqgen)' % (len(events), utils.plural(len(events)), sum(len(l['unique_ids']) for l in events), time.time()-start, sum(reco.validation_values['bpp-times'])))

    if not args.im_a_subproc:
        print('simulating')

    n_per_proc_list, work_fnames = [], []
    if args.n_procs == 1:
        make_events(args.n_sim_events, args.workdir, args.outfname, [random.randint(0, numpy.iinfo(numpy.int32).max) for _ in range(args.n_sim_events)], heavy_chain_events=heavy_chain_events)
    else:
        cmdfos = run_sub_cmds(working_gldir)

    if not args.im_a_subproc and args.allele_prevalence_fname is not None:  # check final prevalence freqs
        glutils.check_allele_prevalence_freqs(args.outfname, glfo, args.allele_prevalence_fname, only_region='v', debug=args.debug>1)
        if args.allele_prevalence_fname == default_prevalence_fname:
            os.remove(default_prevalence_fname)
    if not args.im_a_subproc and args.allowed_cdr3_lengths is not None:  # check final cdr3 lengths
        _, alist, _ = utils.read_output(args.outfname)
        clens = [l['cdr3_length'] for l in alist]
        clcounts = {l : clens.count(l) for l in set(clens)}
        print('  final cdr3 lengths: %s (with counts: %s)' % (' '.join(str(l) for l in sorted(set(clens))), ',  '.join('%d %d'%(l, c) for l, c in sorted(list(clcounts.items()), key=operator.itemgetter(0)))))

    if working_gldir is not None:
        glutils.remove_glfo_files(working_gldir, args.locus)
    if args.n_procs > 1:
        for iproc in range(args.n_procs):
            os.rmdir(cmdfos[iproc]['logdir'])
    if not args.im_a_subproc:  # remove the dummy output file if necessary (if --outfname isn't set on the command line, we still write to a temporary output file so we can do some checks)
        if args.outfname == processargs.get_dummy_outfname(args.workdir):  # note that we *don't* want to use the light_chain arg here, that's only used in the parent process that calls the heavy and light subprocesses. In the parent process we don't want anything removed here since we still need to read both files, and in the subprocesses we want light_chain to be False
            work_fnames.append(args.outfname)
        utils.rmdir(args.workdir, fnames=work_fnames)

# ----------------------------------------------------------------------------------------
def read_qti_file(args):
    if args.queries_to_include_fname is None:
        return None
    uid_only_sfos = []
    if utils.getsuffix(args.queries_to_include_fname) == '.yaml':  # this is pretty ugly, but i'm not sure there's a better way given current constraints
        seqfos = utils.read_seqfos(args.queries_to_include_fname)
        if args.locus is not None and any('locus' in s for s in seqfos):  # remove any that have locus info that doesn't match args.locus (if args.locus is None, we're combining loci)
            seqfos = [s for s in seqfos if 'locus' not in s or s['locus'] == args.locus]
        full_sfos = [s for s in seqfos if s['seq'] is not None]
        uid_only_sfos = [s for s in seqfos if s['seq'] is None]
        more_input_info = collections.OrderedDict([(s['name'], {'unique_ids' : [s['name']], 'seqs' : [s['seq']]}) for s in full_sfos])
        utils.add_input_meta_keys(set(k for sfo in seqfos for k in sfo))
        for sfo in full_sfos:
            for mkey in set(utils.input_metafile_keys) & set(sfo):
                lkey = utils.input_metafile_keys[mkey]
                more_input_info[sfo['name']][lkey] = [sfo[mkey]]
    else:
        more_input_info, _, _ = seqfileopener.read_sequence_file(args.queries_to_include_fname, args.is_data, args=args)
    print('  --queries-to-include-fname: adding %d extra sequence%s%s from %s' % (len(more_input_info), utils.plural(len(more_input_info)), '' if len(uid_only_sfos)==0 else ' (plus %d uid-only queries)'%len(uid_only_sfos), args.queries_to_include_fname))
    queries_to_add = list(more_input_info.keys()) + [s['name'] for s in uid_only_sfos]
    if args.queries_to_include is None:
        args.queries_to_include = []
    if len(set(queries_to_add) & set(args.queries_to_include)) > 0:  # could be actual overlap, or could be we've called this fcn twice
        print('  %s overlap between --queries-to-include and extra input meta info in %s (maybe accidentally called this fcn more than once?): %s' % (utils.wrnstr(), args.queries_to_include_fname, ' '.join(set(queries_to_add) & set(args.queries_to_include))))
    args.queries_to_include += [k for k in queries_to_add if k not in args.queries_to_include]  # whatever, just remove any duplicates
    return more_input_info

# ----------------------------------------------------------------------------------------
def read_inputs(args, actions):
    if actions[0] == 'cache-parameters':  # for parameter caching, use the default in data/germlines/human unless something else was set on the command line
        gldir = args.default_initial_germline_dir if args.initial_germline_dir is None else args.initial_germline_dir
    elif runs_on_existing_output(actions[0]):
        if utils.getsuffix(args.outfname) == '.csv':  # old way
            if args.initial_germline_dir is not None:
                gldir = args.initial_germline_dir
            elif args.parameter_dir is not None and os.path.exists(utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir):
                gldir = utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir
            else:
                raise Exception('couldn\'t guess germline info location with deprecated .csv output file: either set it with --intitial-germline-dir or --parameter-dir, or use .yaml output files so germline info is written to the same file as the rest of the output')
        elif utils.getsuffix(args.outfname) == '.yaml':  # new way
            gldir = None  # gets set when we read the glfo from the yaml in partitiondriver
        else:
            raise Exception('unhandled annotation file suffix %s' % args.outfname)
    else:
        if args.initial_germline_dir is not None:
            gldir = args.initial_germline_dir
        elif args.parameter_dir is not None and os.path.exists(utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir):
            gldir = utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir
        else:
            tstr = ''
            if args.parameter_dir is not None:
                if not os.path.exists(utils.parameter_type_subdir(args, args.parameter_dir)):
                    raise Exception('--parameter-dir (with --parameter-type \'%s\') doesn\'t exist: %s' % (utils.get_parameter_type(args, args.parameter_dir), utils.parameter_type_subdir(args, args.parameter_dir)))
                if not os.path.exists(utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir):
                    tstr = ' (--parameter-dir was set to %s, but the corresponding glfo dir %s doesn\'t exist)' % (args.parameter_dir, utils.parameter_type_subdir(args, args.parameter_dir) + '/' + glutils.glfo_dir)
            raise Exception('couldn\'t guess germline info location: set it with either --intitial-germline-dir or --parameter-dir%s' % tstr)

    glfo = None
    if gldir is not None:
        template_glfo = None
        if args.sanitize_input_germlines:
            print('   using default germline dir %s for template glfo (e.g. for codon positions)' % args.default_initial_germline_dir)
            template_glfo = glutils.read_glfo(args.default_initial_germline_dir, args.locus)
        glfo = glutils.read_glfo(gldir, args.locus, only_genes=args.only_genes, template_glfo=template_glfo, add_dummy_name_components=args.sanitize_input_germlines, debug=2 if args.sanitize_input_germlines else False)

    simglfo = None
    if not args.is_data:
        if args.infname is not None and utils.getsuffix(args.infname) == '.yaml':
            simglfo = None  # cause we'll read it from the input file
        elif args.simulation_germline_dir is not None:  # if an explicit dir was set on the command line
            simglfo = glutils.read_glfo(args.simulation_germline_dir, locus=args.locus)  # probably don't apply <args.only_genes> (?)
        elif args.outfname is not None and args.infname is None:
            print('  note: setting --infname to value from --outfname since --is-simu was set and we need input simulation info')
            args.infname = args.outfname
        else:
            raise Exception('couldn\'t find simulation germline info: either set --infname to a .yaml simulation file, or set --simulation-germline-dir')

    more_input_info = read_qti_file(args)  # also set args.queries_to_include

    if args.infname is None:  # put this *after* setting queries_to_include from a file, since we need that set
        return None, None, glfo, simglfo

    input_info, reco_info, yaml_glfo = seqfileopener.read_sequence_file(args.infname, args.is_data, n_max_queries=args.n_max_queries, args=args, simglfo=simglfo, more_input_info=more_input_info, dont_add_implicit_info=args.dont_add_simu_implicit_info)
    if len(input_info) == 0:
        sys.exit(0)
    if not args.is_data and yaml_glfo is not None:  # NOTE is is extremely important that <glfo> doesn't get set to the true info in a simulation yaml
        simglfo = yaml_glfo

    if len(input_info) > 1000 and args.n_procs == 1:
        print('  note: running on %d sequences with only %d process%s. This may be kinda slow, so it might be a good idea to set --n-procs N to the number of processors on your local machine, or look into non-local parallelization with --batch-system.' % (len(input_info), args.n_procs, utils.plural(args.n_procs, prefix='e')))
    if len(input_info) > 1000 and args.outfname is None and actions != ['cache-parameters']:
        print('  note: running on a lot of sequences (%d) without setting --outfname. Which is ok, but there will be no persistent record of the results (except the parameter directory).' % len(input_info))
    return input_info, reco_info, glfo, simglfo

# ----------------------------------------------------------------------------------------
def default_parameter_dir(args):
    if args.paired_loci and args.paired_indir is not None:
        instr = args.paired_indir
    elif args.infname is not None:
        instr = args.infname
    else:
        assert args.action in processargs.actions_not_requiring_input
        return utils.dummy_str  # this is a shitty convention, but code further on crashes if I let the parameter dir be None
    return '_output/%s' % utils.getprefix(instr).replace('/', '_')

# ----------------------------------------------------------------------------------------
def is_subset_action(action):
    return action in ['subset-partition', 'subset-annotate']

# ----------------------------------------------------------------------------------------
def subset_partition(args):
    # ----------------------------------------------------------------------------------------
    def subdir(igroup):
        return '%s/isub-%d' % (args.paired_outdir, igroup)
    # ----------------------------------------------------------------------------------------
    def input_fn(indir):
        return '%s/input-seqs.fa' % indir
    # ----------------------------------------------------------------------------------------
    def ptnfn(bdir, ltmp, lpair=None):
        return paircluster.paired_fn(bdir, ltmp, lpair=lpair, actstr='partition', suffix='.yaml')
    # ----------------------------------------------------------------------------------------
    def get_cmd(infname, outdir, merged_odir=None):
        clist = copy.deepcopy(sys.argv)
        if shutil.which('partis'):  # use version in PATH if it's there (pipx seems to leave two incompatible versions lying around)
            clist[0] = 'partis'
        assert is_subset_action(clist[1])  # ugh, but utils.replace_in_arglist() only handles -- style arg strs
        clist[1] = args.action.split('-')[1]
        utils.remove_from_arglist(clist, '--n-subsets', has_arg=True)
        utils.replace_in_arglist(clist, '--paired-outdir', outdir)
        # utils.replace_in_arglist(clist, '--persistent-cachefname', 'paired-outdir')  # probably not worth copying hmm cache files since the padding won't be right (also i'm not sure they're being read properly? would need to check more carefully) UDPATE now i have a padding fcn that might work for this (utils.re_pad_hmm_seqs())
        if args.input_metafnames is not None:
            utils.replace_in_arglist(clist, '--input-metafnames', '%s/meta.yaml'%(os.path.dirname(infname) if merged_odir is None else merged_odir))
        if merged_odir is None:  # subset runs
            utils.remove_from_arglist(clist, '--paired-indir', has_arg=True)
            utils.replace_in_arglist(clist, '--infname', infname)
            utils.remove_from_arglist(clist, '--is-simu')  # don't want to deal with subsetting input simulation file (and dual-use of --paired-indir), so just ignore --is-simu until merged run
        else:  # merged run
            utils.remove_from_arglist(clist, '--infname', has_arg=True)
            if args.paired_indir is None:  # if --paired-indir wasn't set (i.e. --infname *was* set), we want to use the extract-pairing-info/split-loci output in --paired-outdir
                clist += ['--paired-indir', merged_odir]
            clist += ['--input-partition-fname', merged_odir, '--continue-from-input-partition', '--parameter-dir', '%s/parameters'%merged_odir,
                      '--ignore-sw-pair-info', '--refuse-to-cache-parameters', '--ignore-default-input-metafile']  # (don't need to ignore default input metafile for subset runs since they use --infname)
        return ' '.join(clist)
    # ----------------------------------------------------------------------------------------
    def run_subsets(infname):
        tmp_seed_id = args.seed_unique_id  # UGH (have to unset seed id before reading input file since a) we don't need it here and b) it'll still be a list of length two here, which'll break when it tries to treat it as a single sequence id)
        args.seed_unique_id = None
        input_info, _, _ = seqfileopener.read_sequence_file(infname, True, n_max_queries=args.n_max_queries, args=args, dont_add_implicit_info=True)  # i think setting is_data=True makes sense (we don't want reco_info)
        args.seed_unique_id = tmp_seed_id
        print('%s %d seqs with %d subsets from %s' % (utils.color('blue_bkg', 'running %s'%args.action), len(input_info), args.n_subsets, infname))
        input_seqfos = [s for l in input_info.values() for s in utils.seqfos_from_line(l)]
        sfo_lists = utils.subset_paired_queries(input_seqfos, args.droplet_id_separators, args.droplet_id_indices, n_subsets=args.n_subsets, input_info=None if args.input_metafnames is None else input_info,
                                                seed_unique_ids=args.seed_unique_id, queries_to_include=args.queries_to_include, debug=True)  # don't want to re-split droplet ids etc if pair info is already there (i.e. use input_info if it has pair info, which it probably does if input meta file is set)
        for isub, sfos in enumerate(sfo_lists):
            ofn = ptnfn(subdir(isub), 'igh')
            if os.path.exists(ofn):
                print('    isub-%s output exists, skipping: %s' % (isub, ofn))
                continue
            utils.write_fasta(input_fn(subdir(isub)), sfos)
            if args.input_metafnames is not None:
                metafos = {}
                for sfo in sfos:
                    assert len(input_info[sfo['name']]['unique_ids']) == 1
                    metafos[sfo['name']] = {mk : input_info[sfo['name']][lk][0] for mk, lk in utils.input_metafile_keys.items() if lk in input_info[sfo['name']]}
                utils.jsdump('%s/meta.yaml'%subdir(isub), metafos)
            # NOTE to get args.seed_unique_id working, I'd need to run parameter caching as a separate step here (probably calling remove_action_specific_args()) (plus other stuff below)
            utils.simplerun(get_cmd(input_fn(subdir(isub)), subdir(isub)), logfname='%s/partition.log'%subdir(isub), dryrun=args.dry_run, extra_str='  %s '%utils.color('blue', 'isub %d:'%isub))
    # ----------------------------------------------------------------------------------------
    def merge_subset_files():
        # NOTE it would be nice to merge duplicates in here (from --collapse-duplicate-sequences, but jeez it seems like it'll be a bit fiddly to implement it
        # ----------------------------------------------------------------------------------------
        def print_dup_warning(ptn, wstr):
            n_tot, n_uniq = sum(len(c) for c in ptn), len(set(u for c in ptn for u in c))
            if n_tot != n_uniq:
                print('    %s %s has duplicates in partition: sum of cluster sizes %d vs %d unique queries' % (utils.wrnstr(), wstr, n_tot, n_uniq))
        # ----------------------------------------------------------------------------------------
        def handle_overlap_merge(lpfos, tmp):
            lpfos['cpaths'][ltmp].partitions, lpfos['antn_lists'][ltmp] = utils.merge_overlapping_clusters(lpfos['cpaths'][ltmp].partitions, glfo=lpfos['glfos'][ltmp], antn_list=lpfos['antn_lists'][ltmp], dbgstr='%s: '%utils.locstr(ltmp), debug=args.debug_paired_clustering)
            print_dup_warning(lpfos['cpaths'][ltmp].best(), 'after')
        # ----------------------------------------------------------------------------------------
        lpair_lpfos, locus_lpfos, i_empty = [], [], []
        for isub in range(args.n_subsets):
            tmp_lpfos = paircluster.read_paired_dir(subdir(isub))
            if all(tmp_lpfos[tuple(lpair)]['glfos'] is None for lpair in utils.locus_pairs[args.ig_or_tr]):  # this should mean that we're keeping all unpaired-seqs, and there were no paired annotations for this subset proc, so it copied/linked the single chain files to the joint location
                i_empty.append(isub)
                tmp_lpfos = paircluster.read_paired_dir(subdir(isub), joint=True)
                locus_lpfos.append(tmp_lpfos)
            else:
                lpair_lpfos.append(tmp_lpfos)
        if len(i_empty) > 0:
            print('    --keep-all-unpaired-seqs: using joint partition files for %d/%d empty paired subdir annotations (indices %s)' % (len(i_empty), args.n_subsets, ' '.join(str(i) for i in i_empty)))

        # <subset_merged_lp_infos>: just has subsets merged (i.e. still split into lpairs)
        # <locus_merged_lp_infos>: also has merged heavy chain
        subset_merged_lp_infos = paircluster.concat_lpair_chains(utils.locus_pairs[args.ig_or_tr], lpair_lpfos, dont_deep_copy=True)  # i think it's ok to not deep copy here, since this fcn doesn't deduplicate anything, and once we merge these subset lpfos, we don't use the originals for anything
        if args.seed_unique_id is not None or args.queries_to_include is not None:  # these will have been passed to every subproc, so need to be deduplicated here so they don't show up multiple times in the merged proc
            print('  merging overlaps in subset merged lp_infos')
            for lpair in subset_merged_lp_infos:
                for ltmp in lpair:
                    handle_overlap_merge(subset_merged_lp_infos[lpair], ltmp)
        def seedid(x): assert args.seed_unique_id is None
        # if seed id or queries to include is set, we *need* to keep the duplicates so we can merge overlapping clusters below
        locus_merged_lp_infos = paircluster.handle_concatd_heavy_chain(utils.locus_pairs[args.ig_or_tr], subset_merged_lp_infos, dont_calculate_annotations=args.dont_calculate_annotations,
                                                                       seed_unique_id=seedid('s'), dont_deduplicate=args.seed_unique_id is not None or args.queries_to_include is not None, debug=True) #args.debug_paired_clustering)  # note that we need to deep copy here, since e.g. if we deduplicate the merged igh antn list, we don't want that to change the annotations in the original lists

        if len(locus_lpfos) > 0:
            locus_merged_lp_infos = paircluster.concat_locus_chains(utils.sub_loci(args.ig_or_tr), [locus_merged_lp_infos] + locus_lpfos, dont_deep_copy=True, dbgstr=' (from single-chain results due to --keep-all-unpaired-seqs)', debug=args.debug_paired_clustering)

        if args.seed_unique_id is not None or args.queries_to_include is not None:  # these will have been passed to every subproc, so need to be deduplicated here so they don't show up multiple times in the merged proc
            print('  merging overlaps in locus merged lp_infos')
            for ltmp in utils.sub_loci(args.ig_or_tr):
                handle_overlap_merge(locus_merged_lp_infos, ltmp)

        print('   writing merged outputs from %d subsets to %s' % (args.n_subsets, merged_odir))
        def ofn_fcn(ltmp, lpair=None, joint=None): return ptnfn(utils.fpath(merged_odir), ltmp, lpair=lpair)  # I *think* i don't need joint
        headers = utils.add_lists(list(utils.annotation_headers), args.extra_annotation_columns)
        paircluster.write_lpair_output_files(utils.locus_pairs[args.ig_or_tr], subset_merged_lp_infos, ofn_fcn, headers, use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)
        paircluster.write_concatd_output_files(locus_merged_lp_infos['glfos'], locus_merged_lp_infos['antn_lists'], ofn_fcn, headers, use_pyyaml=args.write_full_yaml_output, cpaths=locus_merged_lp_infos['cpaths'], dont_write_git_info=args.dont_write_git_info, write_light_chain_files=len(locus_lpfos) > 0)
        outfos, metafos = paircluster.get_combined_outmetafos(locus_merged_lp_infos['antn_lists']) #, extra_meta_headers=[h for h in headers if h in utils.reversed_input_metafile_keys]) hmm, i can't just do this, since this probably has a bunch of keys that aren't in the annotation, but otoh i don't want to completely forget about maybe adding the option
        paircluster.write_combined_fasta_and_meta('%s/all-seqs.fa'%merged_odir, '%s/meta.yaml'%merged_odir, outfos, metafos, write_locus_files=True)  # need this meta file (not the original input one) to pick up pair cleaning (maybe for other reasons as well)
        utils.merge_parameter_dirs(merged_odir, subdir, args.n_subsets, ig_or_tr=args.ig_or_tr)  # , include_hmm_cache_files=True (see note about persistent cache files above)
        # maybe don't need this now I've added the locus_lpfos stuff above?
        # if args.keep_all_unpaired_seqs and len(locus_merged_lp_infos['antn_lists']) == 0:
        #     print('    --keep-all-unpaired-seqs: zero length paired results, so merging single chain outputs')
        #     for ltmp in utils.sub_loci(args.ig_or_tr):
        #         utils.merge_yamls(ptnfn(merged_odir, ltmp), [ptnfn(subdir(i), ltmp) for i in range(args.n_subsets)], headers, use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)
    # ----------------------------------------------------------------------------------------
    if args.input_metafnames is None:
        if args.paired_indir is None:
            args.input_metafnames = ['%s/input/meta.yaml' % args.paired_outdir]  # this is where prep_inputs() will put them
        else:
            args.input_metafnames = ['%s/meta.yaml' % args.paired_indir]
        print('  note: --input-metafnames not specified, so setting to default location %s' % args.input_metafnames[0])
    infname = args.infname
    if infname is None and os.path.exists('%s/all-seqs.fa'%args.paired_indir):  # simulation has a single file
        infname = '%s/all-seqs.fa'%args.paired_indir
    if infname is None and any(os.path.exists('%s/%s.fa'%(args.paired_indir, l)) for l in utils.sub_loci(args.ig_or_tr)):  # UGH (data paired indir, need to cat together the per-locus .fa files into a file in the paired outdir)
        infname = '%s/input-seqs.fa' % args.paired_outdir
        if not os.path.exists(infname):
            utils.prep_dir(args.paired_outdir, allow_other_files=True)
            utils.simplerun('cat %s >%s/input-seqs.fa' % (' '.join('%s/%s.fa'%(args.paired_indir, l) for l in utils.sub_loci(args.ig_or_tr)), args.paired_outdir), shell=True, extra_str='    ')
    if infname is None:
        raise Exception('couldn\'t figure out an input file')
    if all(any(os.path.exists(ptnfn(subdir(i), 'ig'+l)) for l in 'hkl') for i in range(args.n_subsets)):  # if e.g. keeping all unpaired seqs, we may only have an output file for one locus, hence the any() (yes this means it won't rerun if it crashes partway through writing the three files)
        print('  all subset partitions exist (e.g. %s)' % ptnfn(subdir(0), 'igh'))
    else:
        run_subsets(infname)

    merged_odir = '%s/merged-subsets' % args.paired_outdir
    if os.path.exists(ptnfn(merged_odir, 'igh')):
        print('  merged subset files exist (e.g. %s)' % ptnfn(merged_odir, 'igh'))
    else:
        merge_subset_files()

    final_odir = '%s' % args.paired_outdir
    if args.action == 'subset-partition':
        utils.simplerun(get_cmd(infname, final_odir, merged_odir=merged_odir), logfname='%s/merged-partition.log'%final_odir, dryrun=args.dry_run)
    else:  # for 'subset-annotate' just link the merged dir
        print('  linking merged outputs to final dir')
        for ltmp in utils.sub_loci(args.ig_or_tr):
            utils.makelink(final_odir, ptnfn(merged_odir, ltmp), ptnfn(final_odir, ltmp), dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_partitiondriver(args):
    if args.paired_loci:
        run_all_loci(args)
        return

    if args.parameter_dir is None:
        args.parameter_dir = default_parameter_dir(args)
        print('  note: --parameter-dir not set, so using default: %s' % args.parameter_dir)

    actions = [args.action]  # do *not* use <args.action> after this (well, for anything other than checking what was actually set on the command line)
    if args.action in ['annotate', 'partition'] and not os.path.exists(args.parameter_dir):
        if args.refuse_to_cache_parameters:
            raise Exception('--parameter-dir %s doesn\'t exist, and --refuse-to-cache-parameters was set' % args.parameter_dir)
        actions = ['cache-parameters'] + actions
        print('  parameter dir does not exist, so caching a new set of parameters before running action \'%s\': %s' % (actions[1], args.parameter_dir))
        if args.seed_unique_id is not None:  # if we're auto parameter caching for/before seed partitioning, we *don't* (yet) want to remove non-clonal sequences, since we need all the non-clonal sequences to get better parameters (maybe at some point we want to be able to count parameters just on this lineage, but for now let's keep it simple)
            raise Exception('if setting --seed-unique-id for \'partition\', you must first explicitly run \'cache-parameters\' in order to ensure that parameters are cached on all sequences, not just clonally related sequences.')

    input_info, reco_info, glfo, simglfo = read_inputs(args, actions)
    if input_info is not None and len(input_info) < args.n_procs:
        print('  note: reducing N procs to the number of seqs %d --> %d' % (args.n_procs, len(input_info)))
        args.n_procs = len(input_info)
    if args.max_n_seqs_to_likelihood_cluster is not None and not args.naive_vsearch and not args.no_naive_vsearch and args.seed_unique_id is None and len(input_info) > args.max_n_seqs_to_likelihood_cluster:
        print('  note: found %d (>%d) seqs in input, so turning on --naive-vsearch (i.e. --fast) to speed things up. You can turn this behavior off by either setting --no-naive-vsearch or increasing --max-n-seqs-to-likelihood-cluster' % (len(input_info), args.max_n_seqs_to_likelihood_cluster))
        args.naive_vsearch = True
    parter = PartitionDriver(args, glfo, input_info, simglfo, reco_info)
    parter.run(actions)
    if not runs_on_existing_output(args.action):  # mostly wanted to avoid rewriting the persistent hmm cache file
        parter.clean()

# ----------------------------------------------------------------------------------------
# it would be really nice if this didn't need to exist, it's complicated (but it's necessitated by needing to auto-cache parameters on the full un-split igh seqs when splitting loci)
# note that i could also just stop having action-specific arguments, but oh well this seems better for now at least
def remove_action_specific_args(clist, argsaction, tmpaction, debug=True):  # remove args specific to command line action <argsaction> that are not also specific to current action <tmpaction>
    if argsaction not in subargs:  # nothing to do
        return
    def anames(argdicts): return set(afo['name'] for afo in argdicts)
    pargs = anames(parent_args)
    act_args = anames(subargs[argsaction]) - anames(subargs[tmpaction])
    iarg = 0
    removed_strs = []
    while iarg < len(clist):
        if clist[iarg][:2] != '--':  # skip arg values
            iarg += 1
            continue
        if clist[iarg] not in pargs and any(astr.find(clist[iarg]) == 0 for astr in act_args):  # have to allow for incompletely-written args, since argparser handles those
            removed_strs.append(clist[iarg])
            clist.pop(iarg)
            if iarg < len(clist) and clist[iarg][:2] != '--':  # if it has an arg value, also remove that
                removed_strs.append(clist[iarg])
                clist.pop(iarg)
        else:
            iarg += 1
    if debug and len(removed_strs) > 0:
        print('    removed %d arg strs that were specific to action \'%s\': %s' % (len(removed_strs), argsaction, ' '.join(removed_strs)))

# ----------------------------------------------------------------------------------------
def run_all_loci(args, ig_or_tr='ig'):
    # ----------------------------------------------------------------------------------------
    def getpdir(ltmp, bpdir=None, lpair=None, write=False):  # pass in <ltmp> of '' to get the base dir (note: i'm adding <lpair> late, and i'm not really sure when we should use it -- i think not always?)
        if bpdir is None:
            if write and args.parameter_out_dir is not None:
                base_dir = args.parameter_out_dir
            else:
                base_dir = utils.non_none([args.parameter_dir, args.paired_outdir, default_parameter_dir(args)])
            bpdir = '%s%s' % (base_dir, '/parameters' if args.parameter_dir is None else '')
        if lpair is not None and (write or os.path.exists('%s/%s'%(bpdir, '+'.join(lpair)))):
            bpdir = '%s/%s'%(bpdir, '+'.join(lpair))  # NOTE this is pretty similar to utils.parameter_type_subdir(), but i think it's better to have it separate?
        return '%s/%s' % (bpdir, ltmp)
    # ----------------------------------------------------------------------------------------
    def getifn(ltmp, lpair=None):
        if args.paired_indir is not None:
            yfn, ffn = [paircluster.paired_fn(args.paired_indir, ltmp, lpair=lpair, suffix=sx) for sx in ['.yaml', '.fa']]
            if [os.path.exists(f) for f in [yfn, ffn]].count(True) == 2:
                raise Exception('both %s and %s exist, not sure which to use' % (yfn, ffn))
            return yfn if os.path.exists(yfn) else ffn
        elif args.infname is not None:  # if we ran split-input.py on --infname
            return '%s/%s.fa' % (utils.non_none([args.paired_outdir, args.workdir]), ltmp)
        else:
            return 'XXX-DUMMY-XXX'  # ok this kind of sucks, but i think it's better than rewriting all the stuff for merge-paired-partitions (which is the one that doesn't need input stuff specified)
    # ----------------------------------------------------------------------------------------
    def getodir(lpair=None, single_chain=False, is_plotting=False):
        odir = utils.non_none([args.paired_outdir, args.workdir])
        if is_plotting and args.plotdir is not None:
            odir = args.plotdir
        return odir + paircluster.subd(seed_unique_id=args.seed_unique_id, lpair=lpair, single_chain=single_chain)
    # ----------------------------------------------------------------------------------------
    def getplotdir(ltmp, lpair=None, single_chain=False, tmpaction=None):
        return '%s/plots%s%s' % (getodir(lpair=lpair, single_chain=single_chain, is_plotting=True), ('/'+ltmp) if ltmp is not None else '', '/parameters' if tmpaction=='cache-parameters' else '')
    # ----------------------------------------------------------------------------------------
    def getofn(ltmp, joint=False, partition_only=False, trees=False, lpair=None, input_meta=False, fasta=False, persistent_cache=False, recursed=False, gctrees=False):  # NOTE <joint> and <lpair> do *not* always have the same effect (e.g. differs with single-chain/)
        if trees:
            bname = '%strees.nwk' % (ltmp+'-' if gctrees else '')
        elif args.action == 'simulate':
            return paircluster.paired_fn(getodir(), ltmp, lpair=lpair, suffix='.yaml')  # kind of weird that this getodir() call doesn't need any args
        else:
            sfx = 'fa' if fasta else 'yaml'
            if args.action == 'simulate':
                bstr = 'simu'
            elif input_meta:
                bstr = 'meta'
            elif fasta:
                bstr = 'input-seqs'
            elif persistent_cache:
                bstr = 'persistent-cache'
                sfx = 'csv'
            elif args.action in ['view-output', 'get-selection-metrics'] and not recursed and not os.path.exists(getofn(ltmp, joint=joint, partition_only=partition_only, trees=trees, lpair=lpair, recursed=True)):  # ok this is ugly, but it's nice to have 'view-output' automatically figure out if it's simulation or inference output (and yes it's also nice to have different file names for them (partition-igh.yaml vs igh.yaml))
                bstr = ''
            else:
                bstr ='partition'
            bname = '%s%s%s.%s' % (bstr, '' if ltmp is None or bstr=='' else '-', '' if ltmp is None else ltmp, sfx)
        ofn = '%s/%s' % (getodir(lpair=lpair, single_chain=not joint and lpair is None), bname)
        if partition_only:
            ofn = utils.insert_before_suffix('-only-partition', ofn)
        return ofn
    # ----------------------------------------------------------------------------------------
    def get_meta_fns(ltmp, joint=False, lpair=None, tdbg=False):
        mfnames = []
        if args.input_metafnames is not None and ltmp is None:  # if ltmp is set, we're rewriting the meta files for each locus
            mfnames += args.input_metafnames
        elif os.path.exists(getofn(ltmp, joint=joint, lpair=lpair, input_meta=True)):
            mfnames.append(getofn(ltmp, joint=joint, lpair=lpair, input_meta=True))
        if not args.ignore_default_input_metafile and getifn('meta') not in mfnames and os.path.exists(getifn('meta')):  # look in a default-ish location to which e.g. paircluster.write_combined_fasta_and_meta() writes it
            if tdbg:
                print('  note: adding input metafname from default location (%s)' % getifn('meta'))
            mfnames.append(getifn('meta'))
        return mfnames
    # ----------------------------------------------------------------------------------------
    def get_merged_fn(ftype):  # ok this name sucks, but this fcn is for in/out files with all loci together
        return '%s/%s.%s' % (getodir(), ftype, 'fa' if 'seqs' in ftype else 'yaml')
    # ----------------------------------------------------------------------------------------
    def seedid(ltmp):
        if args.seed_unique_id is None or ltmp not in args.seed_loci:
            return None
        assert len(args.seed_unique_id) == 2 and len(args.seed_loci) == 2  # yes, this is also tested elsewhere
        return utils.get_single_entry([u for u, l in zip(args.seed_unique_id, args.seed_loci) if l == ltmp])
    # ----------------------------------------------------------------------------------------
    def sloci(tdict=None):  # if tdict is set, only return those that're in tdict
        rloci = utils.sub_loci(ig_or_tr) # if args.seed_unique_id is None else args.seed_loci
        if tdict is not None:
            rloci = [l for l in rloci if l in tdict]
        return rloci
    # ----------------------------------------------------------------------------------------
    def spairs():
        return utils.locus_pairs[ig_or_tr]  # if args.seed_unique_id is None else [args.seed_loci]
    # ----------------------------------------------------------------------------------------
    def cov_cmd():
        return 'coverage3 run --append'
    # ----------------------------------------------------------------------------------------
    def prep_inputs():
        has_input_meta_pair_info = False
        if args.input_metafnames is not None:  # need to see if we were passed any pair info
            has_input_meta_pair_info = any('paired-uids' in mfo for mfo in utils.read_json_yamls(args.input_metafnames).values())
            if has_input_meta_pair_info:
                print('  note: found \'paired-uids\' in --input-metafnames, so not extracting pairing info')
        mfname = getofn(None, joint=True, input_meta=True)
        if is_subset_action(args.action):  # UGH
            mfname = '%s/input/%s' % (os.path.dirname(mfname), os.path.basename(mfname))
        if not args.guess_pairing_info and not args.no_pairing_info and not has_input_meta_pair_info:
            if utils.getsuffix(args.infname) == '.yaml':
                new_fn = '%s/input-seqs.fa' % utils.non_none([args.paired_outdir, args.workdir])
                print('  note: converting input .yaml file to .fa so we can extract pairing info + split loci (new --infname: %s)' % new_fn)
                utils.simplerun('%s/bin/parse-output.py %s %s' % (utils.get_partis_dir(), args.infname, new_fn), extra_str='    ')
                args.infname = new_fn
            cmd = '%s%s/bin/extract-pairing-info.py %s %s' % (cov_cmd()+' ' if args.prepend_coverage_command else '', args.partis_dir, args.infname, mfname)
            if args.n_max_queries > 0:
                cmd += ' --n-max-queries %d' % args.n_max_queries
            if args.n_random_queries is not None:
                cmd += ' --n-random-queries %d --random-seed %d' % (args.n_random_queries, args.random_seed)
            if args.droplet_id_separators is not None:
                cmd += ' --droplet-id-separators %s' % args.droplet_id_separators
            if args.droplet_id_indices is not None:
                cmd += ' --droplet-id-indices %s' % ':'.join(str(i) for i in args.droplet_id_indices)
            try:
                utils.simplerun(cmd, dryrun=args.dry_run)
            except:
                elines = traceback.format_exception(*sys.exc_info())
                print(utils.pad_lines(''.join(elines)))
                print('  note: couldn\'t extract 10x-style droplet ids from sequence ids (see above), continuing without pairing info')
        if is_subset_action(args.action):  # don't want loci to be split, since we pass all loci to each subset process in one file
            return
        cmd = '%s%s/bin/split-loci.py %s --outdir %s' % (cov_cmd()+' ' if args.prepend_coverage_command else '', args.partis_dir, args.infname, getodir())
        if args.reverse_negative_strands:
            cmd += ' --reverse-negative-strands'
        if os.path.exists(mfname):  # NOTE this won't show up up in dry_run since it won't exist yet
            cmd += ' --input-metafname %s' % mfname
        if args.guess_pairing_info:
            cmd += ' --guess-pairing-info'
        if args.droplet_id_separators is not None:
            cmd += ' --droplet-id-separators %s' % args.droplet_id_separators
        if args.droplet_id_indices is not None:
            cmd += ' --droplet-id-indices %s' % ':'.join(str(i) for i in args.droplet_id_indices)
        if args.n_max_queries > 0:
            cmd += ' --n-max-queries %d' % args.n_max_queries
        if args.n_random_queries is not None:
            cmd += ' --n-random-queries %d --random-seed %d' % (args.n_random_queries, args.random_seed)
        utils.simplerun(cmd, dryrun=args.dry_run)
        if args.paired_outdir is None:
            no_pairing_info = not os.path.exists(mfname)
            work_fnames.extend(paircluster.paired_dir_fnames(getodir(), no_pairing_info=no_pairing_info and args.is_data, include_failed=True, include_meta=True))  # simulation input files should have pairing info (which is then overriden if an input meta file is specified)
    # ----------------------------------------------------------------------------------------
    def auto_cache_params():  # note that this somewhat duplicates code in run_partitiondriver() above (but they need to be separate because of the all-igh-seqs vs only-paired-with-igk/l thing)
        if runs_on_existing_output(args.action):  # don't need a parameter dir to merge existing paired partitions
            return False
        existing_pdir_loci = [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(getpdir(l))]
        if len(existing_pdir_loci) == 0:
            if args.refuse_to_cache_parameters:
                raise Exception('parameter dirs don\'t exist, but --refuse-to-cache-parameters was set: %s' % ' '.join(getpdir(l) for l in utils.sub_loci(ig_or_tr)))
            return True  # all missing, so auto cache parameters
        elif len(existing_pdir_loci) == len(utils.sub_loci(ig_or_tr)):
            return False  # all there, don't auto-cache
        else:
            print('  %s parameters exist for some (%s) but not all (%s) loci, so not auto-caching parameters into %s' % (utils.color('yellow', 'warning'), ' '.join(existing_pdir_loci), ' '.join(set(utils.sub_loci(ig_or_tr)) - set(existing_pdir_loci)), getpdir('')))
            return False
    # ----------------------------------------------------------------------------------------
    def get_n_events(l_locus, arg_events, tdbg=False):
        lfrac = args.light_chain_fractions[l_locus]
        def round_up(): return l_locus=='igk' if lfrac==0.5 else lfrac>0.5  # have to pick one of 'em if they're both exactly 0.5
        n_events = (math.ceil if round_up() else math.floor)(arg_events * lfrac)  # round one of them up, the other one down
        if tdbg:
            print('  %s+%s: %d * %.2f = %d events' % (h_locus, l_locus, arg_events, lfrac, n_events))
        return int(n_events)
    # ----------------------------------------------------------------------------------------
    def translate_locus_names(clist, lpair, ltmp):
        if args.treefname is not None:  # have to translate uids (i.e. add locus)
            tfn = getofn(ltmp, trees=True, lpair=lpair, gctrees=True)
            treeutils.write_translated_trees(tfn, translation_fcn=lambda x: '%s-%s'%(x, ltmp), infname=args.treefname)
            utils.replace_in_arglist(clist, '--treefname', tfn)
            work_fnames.append(tfn)
        if args.queries_to_include_fname is not None: # NOTE this isn't really right, since the file we read here has uids repeated, once for h and once for l, so when the file we write here gets read, the duplicate uid translator adds "-2", and the seqs are messed up then, but it's actually ok for now since we only need the uids in there
            seqfos = utils.read_fastx(args.queries_to_include_fname)
            new_sfos = [{'name' : '%s-%s'%(s['name'], ltmp), 'seq' : s['seq']} for s in seqfos]
            qfn = utils.insert_before_suffix('-qti', getofn(ltmp, lpair=lpair, fasta=True))
            utils.write_fasta(qfn, new_sfos)
            utils.replace_in_arglist(clist, '--queries-to-include-fname', qfn)
            work_fnames.append(qfn)
        if args.input_partition_fname is not None:
            _, _, cpath = utils.read_output(args.input_partition_fname)
            optn = [['%s-%s'%(u, ltmp) for u in c] for c in cpath.best()]
            pfn = '%s/input-partitions-%s.yaml' % (getodir(lpair=lpair), ltmp)
            utils.write_only_partition(pfn, optn)
            utils.replace_in_arglist(clist, '--input-partition-fname', pfn)
            work_fnames.append(pfn)
    # ----------------------------------------------------------------------------------------
    def run_step(tmpaction, ltmp, auto_cache=False, skip_missing_input=False, skip_missing_output=False, lpair=None, joint=False):
        # ----------------------------------------------------------------------------------------
        def prep_args(ltmp):
            clist = copy.deepcopy(sys.argv)
            if shutil.which('partis'):  # use version in PATH if it's there (pipx seems to leave two incompatible versions lying around)
                clist[0] = 'partis'
            if any('=' in a for a in clist):  # arg (have to handle --arg=val style syntax)
                new_clist = []
                for cstr in clist:
                    if '=' in cstr and hasattr(args, cstr.split('=')[0].lstrip('-').replace('-', '_')):
                        print('    note: replacing \'%s\' --> \'%s\' in clist' % (cstr, cstr.split('=')))
                        new_clist += cstr.split('=')
                    else:
                        new_clist.append(cstr)
                clist = new_clist
            utils.remove_from_arglist(clist, '--paired-loci')
            utils.remove_from_arglist(clist, '--single-light-locus', has_arg=True)
            utils.remove_from_arglist(clist, '--dry-run')
            utils.remove_from_arglist(clist, '--reverse-negative-strands')
            utils.remove_from_arglist(clist, '--paired-indir', has_arg=True)
            utils.remove_from_arglist(clist, '--paired-outdir', has_arg=True)
            utils.remove_from_arglist(clist, '--light-chain-fractions', has_arg=True)
            utils.remove_from_arglist(clist, '--seed-loci', has_arg=True)
            utils.remove_from_arglist(clist, '--random-seed-seq')
            utils.remove_from_arglist(clist, '--dont-calculate-annotations')  # we need the single-chain partition steps to write annotations, although maybe it's just for pairing info so maybe we could write that separately, but whatever
            utils.remove_from_arglist(clist, '--guess-pairing-info')
            utils.remove_from_arglist(clist, '--no-pairing-info')
            utils.remove_from_arglist(clist, '--n-random-queries', has_arg=True)
            utils.remove_from_arglist(clist, '--n-max-queries', has_arg=True)
            utils.remove_from_arglist(clist, '--existing-output-run-cfg', has_arg=True)
            utils.remove_from_arglist(clist, '--remove-nonfunctional-seqs')  # have to do this after we have both chains, so we can remove both seqs if either one is nonfunctional
            if tmpaction == 'view-output':
                utils.remove_from_arglist(clist, '--input-partition-fname', has_arg=True)
            if args.prepend_coverage_command:
                clist = cov_cmd().split() + clist
            utils.remove_from_arglist(clist, '--prepend-coverage-command')
            if 'single' not in args.existing_output_run_cfg:
                utils.remove_from_arglist(clist, '--get-selection-metrics')
                utils.remove_from_arglist(clist, '--plot-partitions')
            if tmpaction in ['generate-trees', 'simulate']:
                rseedstr = utils.get_val_from_arglist(clist, '--random-seed') if '--random-seed' in clist else str(args.random_seed)  # could probably just use args.random_seed for both cases?
                # need to send a different seed both to h+k and h+l, and to each of h/l within each locus pair (even if --random-seed wasn't set on the command line, since otherwise the subprocs will use the time, which may be the same). Note that sending the same seed to different locus simulations leaves them correlated, which screws up locus correlations [omfg UGH that bug was hard to find]  NOTE can't just increment it (e.g. have 1 for h+k and 2 for h+l), since then if there's multiple subprocs they'll potentially mostly have the same seeds [like this: utils.locus_pairs[ig_or_tr].index(lpair)]
                utils.replace_in_arglist(clist, '--random-seed', str(utils.hashint(ltmp + '+'.join(lpair) + rseedstr)))  # NOTE sends same seed to tree generation and igh, which i think is fine?
                for tstr in ['', 'reco-', 'shm-']:  # deal with simulation {,reco,shm}-parameter dir stuff
                    argstr = '--%sparameter-dir' % tstr
                    if utils.is_in_arglist(clist, argstr) > 0:
                        utils.replace_in_arglist(clist, argstr, getpdir(ltmp, bpdir=utils.get_val_from_arglist(clist, argstr), lpair=lpair))  # the trees aren't really associated with the locus if we're generating trees, but we just use it to get the mutation rate (and maybe that doesn't really get used?)  # NOTE this will use the heavy chain parameter dir and ignore the light chain one, which I think is ok (it's just used for tree scaling, and it would be better to eventually have different scaling for light chain, but it should be fine)
            if tmpaction == 'generate-trees':
                utils.remove_from_arglist(clist, '--n-procs', has_arg=True)
                utils.remove_from_arglist(clist, '--n-sim-events', has_arg=True)
                utils.remove_from_arglist(clist, '--allowed-cdr3-lengths', has_arg=True)  # format will confuse it
                clist += ['--generate-trees']
                utils.replace_in_arglist(clist, '--outfname', getofn(ltmp, trees=True, lpair=lpair))
                work_fnames.append(getofn(ltmp, trees=True, lpair=lpair))
                utils.replace_in_arglist(clist, '--n-trees', str(get_n_events(lpair[1], utils.non_none([args.n_trees, args.n_sim_events]))), insert_after='--generate-trees')
            else:
                utils.insert_in_arglist(clist, ['--locus', ltmp], args.action)
                if args.action == 'simulate':
                    clist += ['--choose-trees-in-order', '--outfname', getofn(ltmp, lpair=lpair)]
                    if args.input_simulation_treefname is None:
                        clist += ['--input-simulation-treefname', getofn(ltmp, trees=True, lpair=lpair)]
                    if args.paired_correlation_values is not None and not utils.has_d_gene(ltmp):
                        clist += ['--heavy-chain-event-fname', getofn('igh', lpair=lpair)]
                    if args.only_genes is not None:
                        utils.replace_in_arglist(clist, '--only-genes', ':'.join([g for g in args.only_genes if utils.get_locus(g)==ltmp]))
                    utils.replace_in_arglist(clist, '--n-sim-events', str(get_n_events(lpair[1], args.n_sim_events)), insert_after='--choose-trees-in-order')
                    if args.allowed_cdr3_lengths is not None:
                        if ltmp in args.allowed_cdr3_lengths:
                            utils.replace_in_arglist(clist, '--allowed-cdr3-lengths', ':'.join(args.allowed_cdr3_lengths[ltmp].split('c')))  # ick, i hate adding some new split rule with 'c' here (but not sure what else to do)
                        else:
                            utils.remove_from_arglist(clist, '--allowed-cdr3-lengths', has_arg=True)
                else:
                    if auto_cache:
                        clist[clist.index(args.action)] = 'cache-parameters'
                    elif args.action != 'cache-parameters':
                        if '--refuse-to-cache-parameters' not in clist:
                            clist += ['--refuse-to-cache-parameters']
                    if auto_cache or tmpaction in ['annotate', 'plot-partitions', 'view-output']:
                        remove_action_specific_args(clist, args.action, tmpaction)
                    if tmpaction != 'view-output' or not args.is_data:
                        lpt = lpair if tmpaction != 'update-meta-info' else None  # ok this is weird, but for updating meta info we need to read also unpaired seqs (so meta info is set for them). Maybe other actions also? but not sure
                        utils.replace_in_arglist(clist, '--infname', getifn(ltmp, lpair=lpt))  # NOTE this is in some cases replaced below
                    utils.replace_in_arglist(clist, '--parameter-dir', getpdir(ltmp))
                    if args.seed_unique_id is not None and tmpaction in ['plot-partitions', 'view-output', 'get-selection-metrics']:  # 'partition',   # for 'partition' we no longer want to run with seed unique id set, since we're now dropping unseeded clusters beforehand, and don't want the single-chain procs to drop them (since that leaves us with seemingly-unpaired seqs later on)
                        utils.replace_in_arglist(clist, '--seed-unique-id', seedid(ltmp))
                    if args.sw_cachefname is None:
                        clist += ['--sw-cachefname', '%s/sw-cache.yaml'%getpdir(ltmp)]
                    if tmpaction != 'cache-parameters':
                        utils.insert_in_arglist(clist, ['--outfname', getofn(ltmp, joint=joint, lpair=lpair)], '--infname' if '--infname' in clist else None, has_arg=True)
                        if args.paired_outdir is None:
                            work_fnames.append(getofn(ltmp, joint=joint, lpair=lpair))
                            if not joint and lpair is None and getodir(single_chain=True) not in work_fnames:
                                work_fnames.append(getodir(single_chain=True))

            mfnames = get_meta_fns(None, joint=True, tdbg=True)  # , lpair=lpair
            if args.infname is not None and os.path.exists(getofn(None, joint=True, input_meta=True)):  # if we automatically generated an input meta file with extract-pairing-info.py, we need to add it (but this should be the one in the *parent* dir with original pairing info, whereas below in 'annotate' we want the sub-loci ones (after cleaning)
                mfnames.insert(0, getofn(None, joint=True, lpair=lpair, input_meta=True))  # NOTE this is *not* the same file as in the if statement (also note that the insert(0, ) makes sure that any values from the explicitly-passed faile will overwrite those from the file generated by prep_inputs()
            if len(mfnames) > 0:
                utils.replace_in_arglist(clist, '--input-metafnames', ':'.join(sorted(set(mfnames), key=mfnames.index)))  # not sure how i get duplicates, but oh well this gets rid of em

            if args.persistent_cachefname is not None:
                assert args.persistent_cachefname == 'paired-outdir'
                if tmpaction == 'cache-parameters':
                    utils.remove_from_arglist(clist, '--persistent-cachefname', has_arg=True)
                else:
                    utils.replace_in_arglist(clist, '--persistent-cachefname', getofn(ltmp, lpair=lpair, persistent_cache=True))

            if joint and tmpaction in ['annotate', 'view-output', 'get-selection-metrics', 'plot-partitions']:  # used to have: 'args.plot_partitions or' which I'm pretty sure was wrong but too chicken to remove completely
                clist[clist.index(args.action)] = tmpaction
            if args.seed_unique_id is not None and tmpaction == 'partition':
                utils.replace_in_arglist(clist, '--infname', getofn(ltmp, lpair=lpair, fasta=True))  # this does almost the same thing as the line under 'annotate', but note that the tmp .fa files end up in different dirs, since lpair is unset here, whereas it's set for 'annotate'
                utils.remove_from_arglist(clist, '--is-simu')  # at least for now, doing this since we're using fasta input file (which maybe is ok permanently?)
                utils.remove_from_arglist(clist, '--seed-unique-id', has_arg=True)
                utils.remove_from_arglist(clist, '--plotdir', has_arg=True)
            if joint and tmpaction == 'annotate':  # <joint> here (and just above) is just a proxy for if this is running annotations on paired clustering output (as opposed to if we're just running the plain annotate action)
                clist += ['--ignore-sw-pair-info']  # maybe other actions besides 'annotate' should have this, but not sure
                utils.remove_from_arglist(clist, '--queries-to-include-fname', has_arg=True)  # we end up with inconsistent sequences (well, it's just that one is N-padded and the other isn't)
                utils.remove_from_arglist(clist, '--plotdir', has_arg=True)  # i think there isn't anything we want to plot when we're getting the joint partition annotations?
                utils.remove_from_arglist(clist, '--plot-annotation-performance')  # i think there isn't anything we want to plot when we're getting the joint partition annotations?
                utils.remove_from_arglist(clist, '--persistent-cachefname', has_arg=True)  # i think there isn't anything we want to plot when we're getting the joint partition annotations?
                utils.remove_from_arglist(clist, '--annotation-clustering')
                utils.remove_from_arglist(clist, '--annotation-clustering-threshold', has_arg=True)
                utils.replace_in_arglist(clist, '--infname', getofn(ltmp, lpair=lpair, fasta=True))
                utils.replace_in_arglist(clist, '--input-partition-fname', getofn(ltmp, joint=True, partition_only=True, lpair=lpair))
                utils.remove_from_arglist(clist, '--debug', has_arg=True)  # it crashes in utils.print_true_events() if debug is turned on, and I don't care about this debug output anyway
                utils.remove_from_arglist(clist, '--is-simu')  # now that i'm making a tmp .fa file for this step's input, if i want to set --is-simu i'd need to pass in a simulation germline dir, but i think i don't have any reason to set --is-simu anyway?
                utils.remove_from_arglist(clist, '--simultaneous-true-clonal-seqs')  # if this is set, it'll already be taken care of in --input-partition-fname
                if os.path.exists(getofn(ltmp, joint=True, lpair=lpair, input_meta=True)):
                    utils.replace_in_arglist(clist, '--input-metafnames', getofn(ltmp, joint=True, lpair=lpair, input_meta=True))  # have to use the new, rewritten input meta file (with corrected pairing info)
            if not joint and tmpaction in ['annotate', 'partition']:  # <joint> here (and just above) is just a proxy for if this is running annotations on paired clustering output (as opposed to if we're just running the plain action)
                if args.input_partition_fname is not None:
                    utils.replace_in_arglist(clist, '--input-partition-fname', paircluster.paired_fn(utils.get_val_from_arglist(clist, '--input-partition-fname'), ltmp, suffix='.yaml', actstr='partition'))
            if args.dont_calculate_annotations and tmpaction == 'partition':
                clist += ['--use-sw-annotations']  # we need the annotations just for propagating pair info (yes it'd be better to arrange it differently but whatever)
            if tmpaction in ['annotate', 'partition']:  # eh, it's better to do this stuff after reading all the annotations
                utils.remove_from_arglist(clist, '--count-parameters')  # we want to only write these for the annotation step, with the merged partition
                utils.remove_from_arglist(clist, '--count-correlations')

            if args.plotdir is not None and not (joint and tmpaction=='annotate' or args.seed_unique_id is not None and tmpaction=='partition'):  # this adds it if it isn't there, and the latter two cases we just removed it (see above) and don't want it there
                utils.replace_in_arglist(clist, '--plotdir', getplotdir(ltmp, lpair=lpair, single_chain=not joint, tmpaction=tmpaction))
            if args.plot_partitions or tmpaction == 'plot-partitions':  # partitiondriver turns on partition plots if plotdir and input partition fname are both specified, so setting --plot-partitions isn't necessary... but they result in potentially different plotdirs, which is potentially important if you want plots to go somewhere specific
                if tmpaction != 'annotate':  # don't need it
                    utils.remove_from_arglist(clist, '--plot-partitions')
                if tmpaction == 'plot-partitions' and args.is_data:  # for simulation we still need it
                    utils.remove_from_arglist(clist, '--infname', has_arg=True)
                if '--plotdir' not in clist:
                    clist += ['--plotdir', getplotdir(ltmp, lpair=lpair, single_chain=not joint, tmpaction=tmpaction)]
            if tmpaction == 'get-selection-metrics':
                if args.tree_inference_method != 'linearham':
                    utils.remove_from_arglist(clist, '--parameter-dir', has_arg=True)  # remove these just to make the cmd line simpler/easier to read
                utils.remove_from_arglist(clist, '--sw-cachefname', has_arg=True)
                utils.remove_from_arglist(clist, '--input-metafnames', has_arg=True)
                if args.is_data:  # for simulation we still need it
                    utils.remove_from_arglist(clist, '--infname', has_arg=True)
                utils.remove_from_arglist(clist, '--get-selection-metrics')
                utils.remove_from_arglist(clist, '--debug', has_arg=True)
            if args.guess_pairing_info and (tmpaction in ['get-selection-metrics', 'view-output'] or tmpaction in ['cache-parameters', 'annotate', 'partition'] and not joint):  # maybe should do it for other actions? UPDATE adding 'partition' here, but not totally certain i should, and don't have a quick way to test it
                translate_locus_names(clist, lpair, ltmp)
            if tmpaction == 'view-output' and (args.extra_print_keys is None or 'paired-uids' not in args.extra_print_keys):
                utils.replace_in_arglist(clist, '--extra-print-keys', ':'.join(['paired-uids'] + utils.non_none([args.extra_print_keys, []])))

            return clist
        # ----------------------------------------------------------------------------------------
        lpstr = utils.color('blue', ltmp if lpair is None else '+'.join(lpair))  # just for dbg
        ltstr = '' if lpair is None else utils.color('blue', ' '+ltmp)  # just for dbg
        def sdbgstr(dstr, fn): return '%s:%s %s, skipping %s (%s)' % (lpstr, ltstr, dstr, tmpaction, fn)
        if skip_missing_input:
            ifn = getifn(ltmp)
            if args.seed_unique_id is not None and tmpaction == 'partition':  # maybe should also do for 'annotate'?
                ifn = getofn(ltmp, lpair=lpair, fasta=True)
            if not os.path.exists(ifn):
                print(sdbgstr('input file missing', ifn))
                if not args.is_data and '-DUMMY-' in ifn:
                    print('      maybe need to seed --paired-indir (in addition to --paired-outdir)?')
                return
            if os.stat(ifn).st_size == 0:
                print(sdbgstr('zero length input file', ifn))
                return
            if tmpaction in ['annotate', 'partition'] and not os.path.exists(getpdir(ltmp)):  # this probably means that no sequences had annotations (i.e. all failed) for this locus when parameter caching
                print(sdbgstr('parameter dir missing', getpdir(ltmp)))
                return
        if skip_missing_output:
            if not os.path.exists(getofn(ltmp, joint=joint, lpair=lpair)):
                print(sdbgstr('output file missing', getofn(ltmp, joint=joint, lpair=lpair)))
                return
        if args.seed_unique_id is not None and tmpaction=='plot-partitions' and not joint:
            print('    %s: not plotting single chain partitions with seed id set' % ltmp)
            return
        print('%s %s:%s%s' % (utils.color('blue_bkg', tmpaction), lpstr, ltstr, utils.color('blue', ' merged') if joint and lpair is None else ''))
        sys.stdout.flush()
        utils.simplerun(' '.join(prep_args(ltmp)), dryrun=args.dry_run)
    # ----------------------------------------------------------------------------------------
    def rewrite_input_metafo(ltmp, lpair, joint_partition, antn_dict, unpaired_seqs, single_antn_list):  # replace old paired uids with new, fixed ones (also writes tmp input meta file, even if there wasn't an original input meta file)
        old_metafos = {}
        if any(os.path.exists(f) for f in get_meta_fns(None, joint=True)):
            old_metafos = utils.read_json_yamls(get_meta_fns(None, joint=True))
        if args.queries_to_include_fname is not None and utils.getsuffix(args.queries_to_include_fname) == '.yaml':  # this kind of sucks, but if you want to add weird extra edge case functionality, it's gonna make your code messier
            for sfo in utils.read_seqfos(args.queries_to_include_fname):
                if sfo['name'] not in old_metafos:
                    old_metafos[sfo['name']] = {}
                for mkey in set(utils.input_metafile_keys) & set(sfo):
                    old_metafos[sfo['name']][mkey] = sfo[mkey]
        new_metafos = {}
        single_pids = {u : pids for l in antn_dict.values() for u, pids in zip(l['unique_ids'], l['paired-uids'])}
        if args.pair_unpaired_seqs_with_paired_family:
            single_antn_dict = {u : utils.get_single_entry([l for l in single_antn_list if u in l['unique_ids']]) for u in unpaired_seqs[ltmp]}  # NOTE keys are single seqs, i.e. it's *not* equivalent to getting the dict corresponding to single_antn_list
        for jclust in joint_partition:
            for uid in jclust:
                new_metafos[uid] = old_metafos.get(uid, {})  # not sure it's possible for it to be missing, but maybe
                if uid in unpaired_seqs[ltmp]:  # if we determined in paircluster.remove_badly_paired_seqs() that we weren't sure who <uid> was paired with, we need to set <pids> to []
                    if args.pair_unpaired_seqs_with_paired_family and uid in single_antn_dict:
                        pids = utils.per_seq_val(single_antn_dict[uid], 'paired-uids', uid)
                    else:
                        pids = utils.input_metafile_defaults('paired-uids')
                else:
                    pids = single_pids.get(uid, utils.input_metafile_defaults('paired-uids'))
                new_metafos[uid]['paired-uids'] = pids
        if not os.path.exists(getodir(lpair=lpair)):
            os.makedirs(getodir(lpair=lpair))
        utils.jsdump(getofn(ltmp, joint=True, lpair=lpair, input_meta=True), new_metafos)
        work_fnames.append(getofn(ltmp, joint=True, lpair=lpair, input_meta=True))
    # ----------------------------------------------------------------------------------------
    def write_joint_locus(ltmp, lpair, joint_partition, antn_dict, glfo, single_antn_list, unpaired_seqs, ccfs=None):  # have to use <single_antn_list> since <antn_dict> is missing input seqs, i think mostly/entirely from removing badly paired seqs
        # ----------------------------------------------------------------------------------------
        def write_fasta_input():
            single_seqs = {u : s for l in single_antn_list for u, s in zip(l['unique_ids'], l['input_seqs'])}
            missing_seqs = []
            with open(getofn(ltmp, lpair=lpair, fasta=True), 'w') as sfile:
                for uid in [u for c in joint_partition for u in c]:
                    if uid in single_seqs:
                        sfile.write('>%s\n%s\n' % (uid, single_seqs[uid]))
                    else:
                        missing_seqs.append(uid)
            if len(missing_seqs) > 0:
                print('    %s missing %d/%d seqs when writing annotation input: %s' % (utils.color('yellow', 'warning'), len(missing_seqs), sum(len(c) for c in joint_partition), ' '.join(missing_seqs)))  # not sure if this will actually happen or not
            work_fnames.append(getofn(ltmp, lpair=lpair, fasta=True))
            if args.paired_outdir is None:
                work_fnames.append(os.path.dirname(getofn(ltmp, lpair=lpair, fasta=True)))  # there's probably a better place to do this, but oh well
        # ----------------------------------------------------------------------------------------
        if len(joint_partition) == 0:
            return
        rewrite_input_metafo(ltmp, lpair, joint_partition, antn_dict, unpaired_seqs, single_antn_list)
        missing_clusters = []
        for jclust in joint_partition:
            if ':'.join(jclust) not in antn_dict:
                missing_clusters.append(jclust)
                # overlap_clusters = [l['unique_ids'] for l in annotation_list if len(set(l['unique_ids']) & set(jclust)) > 0]
                # print '      %3d: %s' % (len(jclust), ' '.join('%3d/%3d'%(len(set(c) & set(jclust)), len(c)) for c in overlap_clusters))
        if len(missing_clusters) > 0:  # we also get the annotations for the non-missing ones, which should end up the same so it's a bit of a waste, but it's nice to not have to re-read the output file and rewrite it with the non-missing ones
            print('  need to get annotations for %d/%d joint %s clusters' % (len(missing_clusters), len(joint_partition), ltmp))
        jcp = ClusterPath(partition=joint_partition, seed_unique_id=seedid(ltmp))
        if args.seed_unique_id is not None:
            jcp.print_seed_cluster_size(queries_to_include=args.queries_to_include, lstr='joint ')
        partition_lines = jcp.get_partition_lines()
        for pln in partition_lines:
            if ccfs is not None:
                pln['ccf_under'], pln['ccf_over'] = ccfs[ltmp]['joint']
        headers = utils.add_lists(list(utils.annotation_headers), args.extra_annotation_columns)
        if not os.path.exists(getpdir(ltmp)) or args.dont_calculate_annotations:  # no parameter dir, so can't get new annotations (or we were told not to get annotations) NOTE in the former case, could maybe just subset the existing annotations according to the joint partition?
            if not os.path.exists(getpdir(ltmp)):
                print('  %s no %s --parameter-dir was set (or it doesn\'t exist), so we can\'t get annotations for the new joint partition (the joint partition, and the non-joint/non-merged annotations, are still in the output dir)' % (utils.color('yellow', 'warning'), ltmp))
            if args.dont_calculate_annotations:
                print('    not calculating annotations')
            utils.write_annotations(getofn(ltmp, joint=True, lpair=lpair), glfo, [], headers, partition_lines=partition_lines)  # only actually write the partition
            if args.airr_output:
                utils.write_airr_output(utils.replace_suffix(getofn(ltmp, joint=True, lpair=lpair), '.tsv'), [], cpath=jcp, args=args)
        else:
            utils.write_annotations(getofn(ltmp, joint=True, partition_only=True, lpair=lpair), glfo, [], headers, partition_lines=partition_lines)  # write just the partition (no annotations) to a tmp file for use by the annotation step below
            work_fnames.append(getofn(ltmp, joint=True, partition_only=True, lpair=lpair))
            write_fasta_input()
            run_step('annotate', ltmp, lpair=lpair, joint=True)  # this runs partis 'annotate' (it would be nice if i could somehow use partitiondriver.get_cluster_annotations() here, but i don't think it really makes sense/is workable)
    # ----------------------------------------------------------------------------------------
    def combine_simu_chains():
        # ----------------------------------------------------------------------------------------
        def rm_nonfunc_seqs(antns_to_remove, hline, lline):  # remove both seqs (h+l) if either is nonfunctional
            assert len(hline['unique_ids']) == len(lline['unique_ids'])
            functional_iseqs = [iseq for iseq in range(len(hline['unique_ids'])) if utils.is_functional(hline, iseq) and utils.is_functional(lline, iseq)]
            if len(functional_iseqs) == 0:  # none are functional
                print('      %s no functional seqs in family with %d seqs' % (utils.wrnstr(), len(hline['unique_ids'])))
                antns_to_remove.append([hline, lline])
            elif len(functional_iseqs) < len(hline['unique_ids']):  # it's generally very rare for them to all be functional
                print('      removing %d nonfunctional seqs (of %d)' % (len(hline['unique_ids']) - len(functional_iseqs), len(hline['unique_ids'])))
                for ltmp, tln in zip(lpair, [hline, lline]):
                    utils.restrict_to_iseqs(tln, functional_iseqs, lpfos['glfos'][ltmp])
            else:
                print('    all seqs functional, so not removing any: %d / %d' % (len(functional_iseqs), len(hline['unique_ids'])))
        # ----------------------------------------------------------------------------------------
        def rm_empty_annotations(antns_to_remove):
            print('    %s removing %d families with zero functional seqs' % (utils.wrnstr(), len(antns_to_remove)))
            for atn_pair in antns_to_remove:
                for ltmp, tln in zip(lpair, atn_pair):
                    alist = lpfos['antn_lists'][ltmp]
                    i_to_remove = utils.get_single_entry([i for i, l in enumerate(alist) if l['unique_ids']==tln['unique_ids']])  # do it in two steps to make sure we remove exactly one of them
                    lpfos['antn_lists'][ltmp] = [l for i, l in enumerate(alist) if i!=i_to_remove]  # NOTE don't need to update the cpaths since they already get replaced below
                    if len(lpfos['antn_lists'][ltmp]) == 0:
                        print('  %s removed all %s families' % (utils.wrnstr(), ltmp))
        # ----------------------------------------------------------------------------------------
        lp_infos = {}
        for lpair in spairs():
            lpfos = paircluster.read_locus_output_files(lpair, getofn, lpair=lpair, dont_add_implicit_info=not args.debug_paired_clustering, dbgstr='simulation', debug=args.debug_paired_clustering)
            lp_infos[tuple(lpair)] = lpfos
            if None in list(lpfos.values()):
                continue
            print('%s: synchronizing heavy and light chain simulation trees and rewriting output files in %s/' % (utils.color('blue', '+'.join(lpair)), getodir(lpair=lpair)))
            antns_to_remove = []
            for hline, lline in paircluster.get_antn_pairs(lpair, lpfos):
                treeutils.merge_heavy_light_trees(hline, lline, use_identical_uids=False)
                if args.remove_nonfunctional_seqs:
                    rm_nonfunc_seqs(antns_to_remove, hline, lline)
            if len(antns_to_remove) > 0:
                rm_empty_annotations(antns_to_remove)
            for ltmp in lpair:  # replace the old cpaths with new ones using the new uids (this should be faster than translating individual uids as we go)
                lpfos['cpaths'][ltmp] = ClusterPath(partition=utils.get_partition_from_annotation_list(lpfos['antn_lists'][ltmp]), seed_unique_id=seedid(ltmp))
        concat_lpfos = paircluster.handle_concatd_heavy_chain(spairs(), lp_infos, dont_calculate_annotations=args.dont_calculate_annotations, seed_unique_id=seedid(utils.heavy_locus(ig_or_tr)), debug=args.debug_paired_clustering)
        headers = utils.add_lists(list(utils.simulation_headers), args.extra_annotation_columns)
        outfos, metafos = paircluster.get_combined_outmetafos(concat_lpfos['antn_lists']) #, extra_meta_headers=[h for h in headers if h in utils.reversed_input_metafile_keys]) hmm, i can't just do this, since this probably has a bunch of keys that aren't in the annotation, but otoh i don't want to completely forget about maybe adding the option
        outfos = paircluster.modify_simu_pair_info(args, outfos, metafos, lp_infos, concat_lpfos)
        paircluster.write_lpair_output_files(spairs(), lp_infos, getofn, headers, use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)  # NOTE this rewrites the files we just read, which isn't really great in terms of resiliency to crashing partway through things
        paircluster.write_concatd_output_files(concat_lpfos['glfos'], concat_lpfos['antn_lists'], getofn, headers, use_pyyaml=args.write_full_yaml_output, work_fnames=work_fnames if args.paired_outdir is None else None, cpaths=concat_lpfos['cpaths'], dont_write_git_info=args.dont_write_git_info)
        paircluster.write_combined_fasta_and_meta(get_merged_fn('all-seqs'), get_merged_fn('meta'), outfos, metafos)
    # ----------------------------------------------------------------------------------------
    def read_true_files(add_selection_metrics=None, no_plot=False, merged=False, debug=False):
        assert args.paired_indir is not None
        def ofn_fcn(ltmp, lpair=None): return paircluster.paired_fn(args.paired_indir, ltmp, lpair=lpair, suffix='.yaml')
        kwargs = {'add_selection_metrics' : add_selection_metrics, 'plotdir_fcn' : None if (args.plotdir is None or no_plot) else getplotdir, 'dont_add_implicit_info' : not debug, 'lb_tau' : args.lb_tau, 'debug' : debug}
        if merged:  # three of 'em, one for each locus (i.e. with merged igh)
            outfos = paircluster.read_locus_output_files(sloci(), ofn_fcn, **kwargs)
            return outfos
        else:  # separate igh+igk and igh+igl
            lp_infos = paircluster.read_lpair_output_files(spairs(), ofn_fcn, **kwargs)
            return lp_infos
    # ----------------------------------------------------------------------------------------
    def print_joint_vs_input(input_cpaths, input_antn_lists, concat_cpaths, debug=False):
        # ----------------------------------------------------------------------------------------
        def getjcstr(u, jfamilies):  # str for length of <u>'s cluster in <jfamilies>
            jfs = [f for f in jfamilies if u in f]
            if len(jfs) == 0:
                return utils.color('blue', '?')
            else:
                jf = utils.get_single_entry(jfs)
                return ('%d: %s'%(len(jf), ' '.join(jf)))
        # ----------------------------------------------------------------------------------------
        tmpstrs = ['   N clusters with all seqs:'] \
                  + ['%s %4d --> %-4d%s'  % (utils.locstr(l), len(input_cpaths[l].best()), len(concat_cpaths[l].best()), paircluster.chstr(len(input_cpaths[l].best()), len(concat_cpaths[l].best()))) for l in sloci(tdict=concat_cpaths)]
        print('\n        '.join(tmpstrs))
        if debug:
            print('annotations for initial partition clusters, showing (at right) the joint partition clusters into which they were split') # (%s: probably paired with other chain):' % utils.color('blue', '?')
            for ltmp in sloci(tdict=concat_cpaths):
                print('%s' % utils.color('green', ltmp))
                print('  initial:')
                input_cpaths[ltmp].print_partitions(n_to_print=1, extrastr='      ')
                print('  joint:')
                concat_cpaths[ltmp].print_partitions(n_to_print=1, extrastr='      ')
                print('  annotations:')
                input_antn_dict = utils.get_annotation_dict(input_antn_lists[ltmp])
                for tclust in input_cpaths[ltmp].best():  # loop over clusters in the initial partition
                    jfamilies = [c for c in concat_cpaths[ltmp].best() if len(set(tclust) & set(c)) > 0]  # clusters in the joint partition that overlap with this cluster
                    uid_extra_strs = [getjcstr(u, jfamilies) for u in tclust]
                    utils.print_reco_event(input_antn_dict[':'.join(tclust)], uid_extra_strs=uid_extra_strs, extra_str='      ', queries_to_emphasize=args.seed_unique_id)
    # ----------------------------------------------------------------------------------------
    def count_parameter_things(lp_infos, lpair):  # NOTE this duplicates a ton of code in partitiondriver.process_annotation_output()
        # TODO should add perf plotter here as well perfplotter = PerformancePlotter('hmm') if self.args.plot_annotation_performance else None
        lpk = tuple(lpair)
        if lp_infos[lpk]['glfos'] is None:
            return
        if args.count_parameters:
            for ltmp in lpair:
                pcounter = ParameterCounter(lp_infos[lpk]['glfos'][ltmp], args, count_correlations=args.count_correlations)  # count within-locus correlations
                for tline in lp_infos[lpk]['antn_lists'][ltmp]:
                    pcounter.increment(tline)
                if args.plotdir is not None:
                    pcounter.plot(getplotdir(ltmp, lpair=lpair), only_csv=args.only_csv_plots, only_overall=not args.make_per_gene_plots, make_per_base_plots=args.make_per_gene_per_base_plots)
                if not args.dont_write_parameters:
                    pcounter.write('%s/hmm'%getpdir(ltmp, lpair=lpair, write=True))  # note: shouldn't need true_pcounter here, we just need it written once somewhere
                    # self.write_hmms(XXX)  # TODO would require copying a bit from partitiondriver.write_hmms()
        if args.count_correlations:
            ccounter = CorrCounter(paired_loci=lpair)  # count between-locus correlations
            ccounter.incr_cluster_pairs(lp_infos, lpair)
            ccounter.plot('%s/correlations'%getplotdir(None, lpair=lpair), only_mi=True, only_csv=args.only_csv_plots)
    # ----------------------------------------------------------------------------------------
    def remove_unseeded_seqs(debug=True):
        # ----------------------------------------------------------------------------------------
        def read_sw(ltmp):
            if not os.path.exists('%s/sw-cache.yaml'%getpdir(ltmp)):
                return []
            wargs = copy.deepcopy(args)
            wargs.seed_unique_id = seedid(ltmp) if ltmp in args.seed_loci else None
            wargs.locus = ltmp
            wargs.is_data = True  # don't need the simu info, and it would complicate things to make it work
            waterer = Waterer(wargs, None, None, None, None, locus=ltmp)
            waterer.read_cachefile('%s/sw-cache.yaml'%getpdir(ltmp), ignore_seed_unique_id=True, quiet=True)
            return [waterer.info[q] for q in waterer.info['queries']]
        # ----------------------------------------------------------------------------------------
        print('   removing seqs very different from seed seq')
        sw_antns = {l : read_sw(l) for l in utils.sub_loci(ig_or_tr)}
        seed_lines = [utils.get_single_entry([l for l in sw_antns[ltmp] if l['unique_ids'][0]==seedid(ltmp)]) for ltmp in args.seed_loci]
        def hbfcn(ltmp): return utils.get_naive_hamming_bounds('likelihood', overall_mute_freq=numpy.mean([f for l in sw_antns[ltmp] for f in l['mut_freqs']]))  # with 'likelihood', these are the wider bounds, so < lo is almost certainly clonal, > hi is almost certainly not
        _, hi_hbounds = list(zip(*[hbfcn(l) for l in args.seed_loci]))
        def hfrc(tl, sl): return utils.hamming_fraction(tl['naive_seq'], sl['naive_seq'], align_if_necessary=True)
        def hfracs(tlines): return [hfrc(tl, sl) for tl, sl in zip(tlines, seed_lines)]  # shouldn't need to align, but if it does it'll print a warning, which is probably better than crashing? but if it's aligning i think something is probably pretty wrong

        # first add to <final_antns> all the correctly/uniquely-paired seqs with the same cdr3/similar naive to seed seq
        final_antns = {l : [] for l in utils.sub_loci(ig_or_tr)}
        antn_pairs, antn_pair_ids = [], set()  # NOTE this is quite similar to paircluster.find_cluster_pairs(), although we don't have lp_infos here, so had to pass in just the antn lists (but even then it was too slow)
        antn_dicts = {l : utils.get_annotation_dict(sw_antns[l]) for l in args.seed_loci}
        hloc, lloc = args.seed_loci
        for hline in sw_antns[hloc]:
            pids = hline['paired-uids']
            if len(pids) != 1 or ':'.join(pids[0]) not in antn_dicts[lloc]:
                continue
            lline = antn_dicts[lloc][':'.join(pids[0])]
            antn_pairs.append([hline, lline])
            antn_pair_ids |= set([l['unique_ids'][0] for l in [hline, lline]])
        n_removed = {'cdr3' : 0, 'hfrac' : 0}
        for tlines in antn_pairs:
            if any(tl['cdr3_length'] != sl['cdr3_length'] for tl, sl in zip(tlines, seed_lines)):
                n_removed['cdr3'] += 1
                continue
            if args.n_final_clusters is None and args.min_largest_cluster_size is None and any(hf > hb for hf, hb in zip(hfracs(tlines), hi_hbounds)):
                n_removed['hfrac'] += 1
                continue
            for ltmp, tline in zip(args.seed_loci, tlines):
                final_antns[ltmp].append(tline)
        # then add seqs with bad (un/multi) pairing, since we may want the former, and the latter may get fixed by the pairing info cleaner
        n_single_added = {l : 0 for l in args.seed_loci}
        for il, ltmp in enumerate(args.seed_loci):
            other_antn_dict = antn_dicts[args.seed_loci[(il+1)%2]]
            for tline in sw_antns[ltmp]:
                if tline['unique_ids'][0] in antn_pair_ids:  # already decided whether or not to keep it
                    continue
                pids = utils.get_single_entry(tline['paired-uids'])
                if len(pids) == 1:  # if it's correctly/uniquely paired, we should have already evaluated it
                    continue
                if all(u not in other_antn_dict for u in pids):  # if none of the pids are in the <other_antn_dict>, there's no hope that it's correclty paired
                    continue
                if tline['cdr3_length'] != seed_lines[il]['cdr3_length']:
                    continue
                if args.n_final_clusters is None and args.min_largest_cluster_size is None and hfrc(tline, seed_lines[il]) > hi_hbounds[il]:
                    continue
                final_antns[ltmp].append(tline)
                n_single_added[ltmp] += 1
        if debug:
            mtruestr = ''
            if not args.is_data:
                for ltmp in args.seed_loci:
                    true_partition = utils.read_cpath(getifn(ltmp), n_max_queries=args.n_max_queries).best()
                    true_seed_cluster = utils.get_single_entry([c for c in true_partition if seedid(ltmp) in c])
                    final_ids = [u for l in final_antns[ltmp] for u in l['unique_ids']]
                    missing_true_ids = [u for u in true_seed_cluster if u not in final_ids]
                    if len(missing_true_ids) == 0:
                        mtruestr = '    ' + utils.color('green', '-')
                    else:
                        mtruestr = '    ' + utils.color('red', str(len(missing_true_ids)))
            print('                           wrong/no     removed      singly         kept%s' % ('' if args.is_data else '   missing'))
            print('           locus  before   pairing    cdr3  hfrac    added   after  frac%s' % ('' if args.is_data else '   true ids'))
            for ltmp in args.seed_loci:
                n_init, n_final = len(sw_antns[ltmp]), len(final_antns[ltmp])
                print('             %s    %5d    %5d     %5d %5d    %5d   %5d   %.3f%s' % (utils.locstr(ltmp), n_init, n_init - len(antn_pairs), n_removed['cdr3'], n_removed['hfrac'], n_single_added[ltmp], n_final, n_final / float(n_init), mtruestr))

        taken_uids = [u for alist in final_antns.values() for l in alist for u in l['unique_ids']]  # uids that we actually want, that might be related to the seed ab
        all_paired_uids = set([u for alist in final_antns.values() for l in alist for pids in l['paired-uids'] for u in pids])  # all uids that're paired with a uid that we took (includes ones that we took)
        missing_paired_uids = set(all_paired_uids) - set(taken_uids)  # uids that we need to get in the next step
        if len(missing_paired_uids) > 0:  # shouldn't happen any more unless pairing info is messed up, i.e. we added single seqs (we used to remove by cdr3/hfrac on single chains), but it's really bad if we miss one, since then it looks like its partner is unpaired, which screws up paired clustering
            if debug:
                print('     adding seqs paired with uids we kept')
                print('           locus  before  added  after  kept frac')
            for ltmp in utils.sub_loci(ig_or_tr):
                n_before = len(final_antns[ltmp])
                final_antns[ltmp] += [l for l in sw_antns[ltmp] if utils.get_single_entry(l['unique_ids']) in missing_paired_uids]
                if debug:
                    print('             %s   %5d   %5d  %5d    %.3f' % (utils.locstr(ltmp), n_before, len(final_antns[ltmp]) - n_before, len(final_antns[ltmp]), 0 if len(sw_antns[ltmp])==0 else len(final_antns[ltmp]) / float(len(sw_antns[ltmp]))))
            final_uids = [u for alist in final_antns.values() for l in alist for u in l['unique_ids']]
            missing_uids = missing_paired_uids - set(final_uids)
            if debug:
                print('      taken %d, paired with taken %d (total %d)' % (len(taken_uids), len(missing_paired_uids), len(taken_uids) + len(missing_paired_uids)))
                # print '             final %d%s' % (len(final_uids), ' (missing %d)'%len(missing_uids) if len(missing_uids)>0 else '')
            if len(missing_uids) > 0:
                print('  %s missing %d expected uids paired with seqs maybe clonal to seed ab' % (utils.color('red', 'error'), len(missing_uids)))

        qti_queries = []
        if args.queries_to_include_fname is not None:  # need to exclude these, so the single-chain partition procs don't then read them from both infname and queries to include file
            rfcn = utils.read_seqfos if utils.getsuffix(args.queries_to_include_fname) in ['.yaml', '.json'] else utils.read_fastx
            qti_queries = [s['name'] for s in rfcn(args.queries_to_include_fname)]
        for ltmp in utils.sub_loci(ig_or_tr):  # NOTE that e.g. --n-random-queries and --n-max-queries get applied *after* this, in the single chain partition steps (which may or may not be what you want, since that means they're being applied after you've removed sequences in the previous lines)
            if len(final_antns[ltmp]) == 0:
                continue
            utils.write_fasta(getofn(ltmp, fasta=True), [utils.get_single_entry(utils.seqfos_from_line(l)) for l in final_antns[ltmp] if l['unique_ids'][0] not in qti_queries])
            work_fnames.append(getofn(ltmp, fasta=True))
#     # ----------------------------------------------------------------------------------------
# # NOTE god damnit, this doesn't work -- you can't just add them as singletons, since then paired clustering thinks they're supposed to be split, which fucks up the paired partitions
#     # when seed clustering, each single-chain process starts by removing all seqs with cdr3 len different to the seed seqs, but here we have to put some of them back in.
#     # Otherwise the seqs they're paired with would show up as unpaired, so we wouldn't be able to properly split during paired clustering, and they'd end up overmerged
#     def reinstate_missing_uids(single_outfos):
#         print '    reinstating uids from sw info that were missing from single chain partition output (i.e. that were removed with cdr3 length different to seed seq)'
#         print '      reading cached sw results '
#         input_antns = {}
#         for ltmp in utils.sub_loci(ig_or_tr):
#             waterer = Waterer(args, None, None, None, None, locus=ltmp)
#             waterer.read_cachefile('%s/sw-cache.yaml'%getpdir(ltmp), ignore_seed_unique_id=True)
#             input_antns[ltmp] = [waterer.info[q] for q in waterer.info['queries']]
#         print '           locus  before     after  added'
#         uids_to_add = set([u for alist in single_outfos['antn_lists'].values() for l in alist for pids in l['paired-uids'] for u in pids])  # only add ids that are paired with an id that's in the existing output
#         for ltmp in utils.sub_loci(ig_or_tr):  # sloci(single_outfos['glfos']):
#             existing_uids = []  # uids that *weren't* removed during seed clustering
#             if ltmp in single_outfos['glfos']:
#                 existing_uids += [u for c in single_outfos['cpaths'][ltmp].best() for u in c]  # NOTE there are in general lots of duplicate uids in antn_lists, i.e. *don't* use antn_lists here, i.e. *don't* ever use the antn_lists as a partition (also note that different partitions have different uids, especially for seed partitioning)
#             else:
#                 single_outfos['cpaths'][ltmp] = ClusterPath(partition=[])  # this should be the light chain locus that *wasn't* the seed seq's
#                 single_outfos['antn_lists'][ltmp] = []
#             n_added = 0
#             for line in input_antns[ltmp]:
#                 uid = utils.get_single_entry(line['unique_ids'])
#                 if uid in existing_uids or uid not in uids_to_add:
#                     continue
#                 # NOTE i'm not going to bother merging the glfos, which maybe is ok?
#                 single_outfos['cpaths'][ltmp].add_cluster_to_all_partitions([uid], skip_duplicates=args.seed_unique_id is not None)
#                 single_outfos['antn_lists'][ltmp].append(line)
#                 n_added += 1
#             print '             %s   %5d     %5d   %5d' % (utils.locstr(ltmp), len(existing_uids), sum(len(c) for c in single_outfos['cpaths'][ltmp].best()), n_added)
    # ----------------------------------------------------------------------------------------
    def restrict_to_single_partition(single_outfos):  # replace the cpaths with new ones that have only a single partition (the paired clustering atm only handles one partition, and the simplest way to control/enforce which is to restrict to a single one to start with)
        pfcn = 'best'
        if args.n_final_clusters is not None or args.min_largest_cluster_size is not None:
            pfcn = 'last'
            print('    %s%s: using last partition (rather than best) from single-locus output cpaths' % ('--n-final-clusters' if args.n_final_clusters else '', '--min-largest-cluster-size' if args.min_largest_cluster_size else ''))
            print('                     N clusters     N uids')
            print('                   best   last   best   last')
            for ltmp in sloci():
                bptn, lptn = [getattr(single_outfos['cpaths'][ltmp], pf)() for pf in ('best', 'last')]
                print('            %s   %5d %6d   %5d %6d' % (utils.locstr(ltmp), len(bptn), len(lptn), sum(len(c) for c in bptn), sum(len(c) for c in lptn)))
        for ltmp in sloci():
            single_outfos['cpaths'][ltmp] = ClusterPath(partition=getattr(single_outfos['cpaths'][ltmp], pfcn)(), seed_unique_id=seedid(ltmp))
    # ----------------------------------------------------------------------------------------
    def combine_inf_chains():
        print('%s chains' % utils.color('blue_bkg', 'combining'))
        sys.stdout.flush()
        single_outfos = paircluster.read_locus_output_files(sloci(), getofn, dont_add_implicit_info=not args.debug_paired_clustering, dbgstr='partition', debug=args.debug_paired_clustering)  # single chain cpaths (and other info)
        restrict_to_single_partition(single_outfos)
        # doesn't work: (see note above)
        # if args.seed_unique_id is not None:
        #     reinstate_missing_uids(single_outfos)
        if None in list(single_outfos.values()):
            return None
        if not any('paired-uids' in l for alist in single_outfos['antn_lists'].values() for l in alist):
            print('  no paired uids in any annotations, not merging single chain partitions')
            return None
        true_outfos = None if args.is_data else {k : read_true_files(merged=mgd, no_plot=True) for k, mgd in [('merged', True), ('lpairs', False)]}
        paircluster.clean_pair_info(args, single_outfos['cpaths'], single_outfos['antn_lists'], plotdir=getplotdir(None) if args.plotdir is not None else None, performance_outdir=None if args.is_data else getodir(), true_outfos=true_outfos, debug=args.debug_paired_clustering)
        lp_infos = {}
        for lpair in spairs():
            print('%s: synchronizing heavy and light chain cluster paths' % utils.color('blue', '+'.join(lpair)))
            sys.stdout.flush()
            ploci = {ch : l for ch, l in zip(utils.chains, lpair)}
            if any(l not in single_outfos['antn_lists'] or len(single_outfos['antn_lists'][l])==0 for l in lpair):
                print('    no single chain output')
                continue
            lp_cpaths, lp_antn_lists, unpaired_seqs = paircluster.remove_badly_paired_seqs(ploci, single_outfos, keep_all_unpaired=args.keep_all_unpaired_seqs, debug=args.debug_paired_clustering)  # NOTE that if we re-annotated after doing this, we'd in general get different annotations (i.e. this is really the first step in paired clustering, while the rest of the steps happen in the next fcn call)
            if all(len(alist)==0 for alist in lp_antn_lists.values()):
                print('    %s zero sequences remain after removing un/badly paired seqs. Maybe you\'re missing pairing info?' % utils.color('yellow', 'warning'))
                continue
            seedids = None if args.seed_unique_id is None else {l : u for l, u in zip(args.seed_loci, args.seed_unique_id)}
            if args.n_final_clusters is not None:
                print('    %s: --n-final-clusters is set, but there isn\'t yet a way to specify a particular number of clusters during paired clustering. So currently this means we over-merge as far as possible (which is correct only if you set --n-final-clusters to 1)' % utils.color('yellow', 'warning'))
            if args.min_largest_cluster_size is not None:
                print('    %s: --min-largest-cluster-size is set, but there isn\'t yet a way to specify a particular number of clusters during paired clustering. So currently we don\'t do anything specific for --min-largest-cluster-size during paired clustering, so you\'ll probably get smaller clusters than you want. If this is a problem, you can try using --n-final-clusters 1 instead.' % utils.color('yellow', 'warning'))
            mstart = time.time()
            joint_partitions, cluster_pairs, ccfs = paircluster.merge_chains(ploci, lp_cpaths, lp_antn_lists, unpaired_seqs=unpaired_seqs, true_outfos=true_outfos, input_cpaths=single_outfos['cpaths'], input_antn_lists=single_outfos['antn_lists'],
                                                                             seed_unique_ids=seedids, overmerge=args.n_final_clusters is not None, naive_hamming_bound_type=args.paired_naive_hfrac_threshold_type, debug=args.debug_paired_clustering, fail_frac=args.max_ccf_fail_frac if args.seed_unique_id is None else None)  # , iparts={'igl' : 7}
            if args.pair_unpaired_seqs_with_paired_family:
                paircluster.pair_unpaired_seqs_with_paired_family(ploci, unpaired_seqs, cluster_pairs, single_outfos['antn_lists'], debug=args.debug_paired_clustering)
            print('  merge time %.1f' % (time.time() - mstart))
            for ltmp in lpair:  # 'writes' them by running annotate with the joint partition as input partition
                if args.seed_unique_id is not None:  # NOTE it's kind of dumb/shitty to keep the non-seeded h/l pair all the way through this (i.e. to do paired clustering on igl when the seed seq is igh+igk), but we really need to not lose the seqs in there, since then we lose track of their chain, which means it looks like the seqs they're paired with that we *do* care about are unpaired, which screws up their joint partition
                    n_before = len(joint_partitions[ltmp])
                    joint_partitions[ltmp], _ = utils.split_partition_with_criterion(joint_partitions[ltmp], lambda cluster: seedid(ltmp) in cluster)
                    print('     removed %d unseeded clusters (of %d total, leaving %d) from %s joint partition' % (n_before - len(joint_partitions[ltmp]), n_before, len(joint_partitions[ltmp]), ltmp))
                write_joint_locus(ltmp, lpair, joint_partitions[ltmp], utils.get_annotation_dict(lp_antn_lists[ltmp]), single_outfos['glfos'][ltmp], single_outfos['antn_lists'][ltmp], unpaired_seqs, ccfs=ccfs)
            lp_infos[tuple(lpair)] = paircluster.read_locus_output_files(lpair, getofn, lpair=lpair, dont_add_implicit_info=not args.debug_paired_clustering and not args.count_parameters, dbgstr='partition')
            count_parameter_things(lp_infos, lpair)
        if args.pair_unpaired_seqs_with_paired_family and not args.is_data:  # note that this will overwrite the plots/csvs that were written during normal pair info cleaning, but i think that's more or less ok
            paircluster.plot_fraction_correctly_paired(single_outfos['cpaths'], None, true_outfos=true_outfos, antn_lists=single_outfos['antn_lists'], calc_near_fams=True, performance_outdir=None if args.is_data else getodir(), plotdir=getplotdir(None) if args.plotdir is not None else None)
        if args.input_metafnames is not None:  # need to add any new input meta keys to headers
            seqfileopener.read_input_metafo(args.input_metafnames, [])
        concat_lpfos = paircluster.handle_concatd_heavy_chain(spairs(), lp_infos, dont_calculate_annotations=args.dont_calculate_annotations, seed_unique_id=seedid(utils.heavy_locus(ig_or_tr)), debug=args.debug_paired_clustering)
        paircluster.write_concatd_output_files(concat_lpfos['glfos'], concat_lpfos['antn_lists'], getofn, utils.add_lists(list(utils.annotation_headers), args.extra_annotation_columns), use_pyyaml=args.write_full_yaml_output,
                                               work_fnames=work_fnames if args.paired_outdir is None else None, cpaths=concat_lpfos['cpaths'], true_outfos=true_outfos, dont_write_git_info=args.dont_write_git_info, airr_output=args.airr_output, args=args, fail_frac=args.max_ccf_fail_frac if args.seed_unique_id is None else None)
        print_joint_vs_input(single_outfos['cpaths'], single_outfos['antn_lists'], concat_lpfos['cpaths'], debug=args.debug_paired_clustering)
        if args.keep_all_unpaired_seqs and len(lp_infos) == 0:  # ugh
            print('    --keep-all-unpaired-seqs: zero length paired results, so linking single chain output to final/merged results dir, since otherwise there\'d be no final output files')
            for ltmp in utils.sub_loci(ig_or_tr):
                # maybe there's some reason to rewrite this? I don't think so, though (maybe restrict_to_single_partition() does something, or there's some other modification)
                # utils.write_annotations(getofn(ltmp, joint=True), single_outfos['glfos'][ltmp], single_outfos['antn_lists'][ltmp], utils.add_lists(list(utils.annotation_headers), args.extra_annotation_columns), use_pyyaml=args.write_full_yaml_output, dont_write_git_info=args.dont_write_git_info)
                if os.path.exists(getofn(ltmp)):
                    utils.makelink(getodir(), getofn(ltmp), getofn(ltmp, joint=True), extra_str='        ', debug=True)
        # ccfs = None
        # if len(single_outfos['cpaths']['igh'].best()) > 5000:
        #     print '    not comparing single vs joint partitions to save time'
        # else:
        #     ccfs = {}
        #     for ltmp in sloci(tdict=concat_lpfos['cpaths']):
        #         dbg_str = '%s single vs joint concat\'d '%utils.locstr(ltmp)
        #         ccfs[ltmp] = paircluster.compare_partition_pair(single_outfos['cpaths'][ltmp].best(), concat_lpfos['cpaths'][ltmp].best(), dbg_str=dbg_str, cf_label='single', ref_label='joint', add_to_ref=True, debug=True) # antn_list=XXX
        # return ccfs
    # ----------------------------------------------------------------------------------------
    def run_cf_plotdirs(cfpdir, ltmp):
        if not os.path.exists(getplotdir(ltmp, single_chain=True)) or not os.path.exists(getplotdir(ltmp)) or 'sizes' not in utils.non_none([args.partition_plot_cfg, utils.default_ptn_plot_cfg]):
            return
        subd = 'partitions/sizes'
        titlestr = ltmp
        # if ccfs is not None and ltmp in ccfs:
        #     titlestr += ':@specif.%.2f'%ccfs[ltmp][0]
        cmdstr = '%s/bin/compare-plotdirs.py --outdir %s/%s --names single:joint --plotdirs %s/%s:%s/%s' % (args.partis_dir, cfpdir, ltmp, getplotdir(ltmp, single_chain=True), subd, getplotdir(ltmp), subd)
        cmdstr += ' --make-parent-html --add-to-title %s --log xy --translegend=0.07:0.07' % titlestr
        utils.simplerun(cmdstr, dryrun=args.dry_run)

        # NOTE this has to happen *after* running compare-plotdirs.py, since it cleans the dirs
        if not args.dry_run:
            import python.plotting as plotting
            plotting.plot_csim_matrix_from_files(cfpdir + '/' + ltmp, 'csim-matrix', 'single partition', getofn(ltmp, joint=False), 'joint partition', getofn(ltmp, joint=True), 35)
    # ----------------------------------------------------------------------------------------
    def process_all_outputs(tmpaction):
        assert tmpaction in ['update-meta-info'] + ex_out_cfg_actions
        if tmpaction=='update-meta-info' or 'paired' in args.existing_output_run_cfg:
            for lpair in spairs():  # run plots/meta update for each locus in each paired (e.g. igh+igk/) subdir (single-chain plots for each locus are run above)
                for ltmp in lpair:
                    run_step(tmpaction, ltmp, lpair=lpair, joint=True, skip_missing_output=True)
        if tmpaction=='update-meta-info' or 'merged' in args.existing_output_run_cfg:
            for ltmp in sloci():  # heavy chain plots with the joint partition combining igh paired with both light loci (make light chain links just so people know where to find them/consistency)
                if utils.has_d_gene(ltmp):
                    run_step(tmpaction, ltmp, joint=True)
                else:
                    # adding 'plot-partitions' late because something's wrong with the link call (it's crashing), and maybe i just shouldn't make the link for 'plot-partitions'?
                    if tmpaction not in ['update-meta-info', 'plot-partitions'] and os.path.exists(getplotdir(ltmp, lpair=utils.getlpair(ltmp))):  # light chain in the paired subdirs is the same, so for plotting just link the results
                        utils.makelink(getodir(), getplotdir(ltmp, lpair=utils.getlpair(ltmp)), getplotdir(ltmp), dryrun=args.dry_run)
        if tmpaction=='plot-partitions' and args.plotdir is not None and args.sub_plotdir is None:  # if --sub-plotdir is set, we presumably already plotted these, and atm they don't get put in the sub plotdir anyway so we'd just be overwriting them with the same thing
            cfpdir = '%s/cf-plots' % getodir(is_plotting=True)
            if not args.no_partition_plots:
                for ltmp in sloci():  # heavy chain plots with the joint partition combining igh paired with both light loci (make light chain links just so people know where to find them/consistency)
                    run_cf_plotdirs(cfpdir, ltmp)
            if os.path.exists(cfpdir):
                fnoutstr, _ = utils.simplerun('find %s -type f -name *.svg' % cfpdir, return_out_err=True, dryrun=args.dry_run)
                import python.plotting as plotting
                fnlist = sorted(fnoutstr.strip().split('\n'))
                plotting.make_html(cfpdir, fnames=[[f for f in fnlist if 'cluster-sizes' in f], [f for f in fnlist if 'cluster-sizes' not in f]])
    # ----------------------------------------------------------------------------------------
    def get_pntns(smetrics=False):  # name kind of sucks
        if not args.is_data:
            if args.paired_indir is None:
                raise Exception('can\'t read simulation files without setting --paired-indir')
            lp_infos = read_true_files(add_selection_metrics=args.selection_metrics_to_calculate if smetrics and len(args.existing_output_run_cfg) > 0 else None, no_plot=True, debug=args.debug or args.debug_paired_clustering)
            antn_pairs = paircluster.get_all_antn_pairs(lp_infos)  # NOTE this will cause a bunch of warnings (mostly later) if this simulation has multiple cells per droplet/missing pair info
        else:
            lp_infos = paircluster.read_lpair_output_files(spairs(), getofn, dont_add_implicit_info=not args.debug_paired_clustering, read_selection_metrics=smetrics and len(args.existing_output_run_cfg) > 0, dbgstr='partition', debug=args.debug or args.debug_paired_clustering)
            _ = read_qti_file(args)  # sets <args.queries_to_include>
            antn_pairs = paircluster.find_all_cluster_pairs(lp_infos, min_cluster_size=args.min_paired_cluster_size_to_read if smetrics else None, min_cluster_arg_str=' (controlled by --min-paired-cluster-size-to-read)' if smetrics else '')  # , required_keys=['tree-info']
        antn_pairs = sorted(antn_pairs, key=lambda x: sum(len(l['unique_ids']) for l in x), reverse=True)  # sort by the sum of h+l ids (if i could start over i might sort by the number of common ids)
        fake_pntns, mpfo_lists, mtpys = paircluster.make_fake_hl_pair_antns(args, antn_pairs)
        return antn_pairs, fake_pntns, mpfo_lists, mtpys
    # ----------------------------------------------------------------------------------------
    paired_loci_actions = ['cache-parameters', 'annotate', 'partition', 'subset-partition', 'subset-annotate', 'merge-paired-partitions', 'get-selection-metrics', 'plot-partitions', 'view-output', 'simulate', 'update-meta-info', 'write-fake-paired-annotations']
    if args.action not in paired_loci_actions:
        raise Exception('action %s not supported for paired loci (choose from %s)' % (args.action, paired_loci_actions))
    if args.paired_outdir is None:
        print('  note: --paired-outdir is not set, so there will be no persistent record of the results (except the parameter directory).')
    utils.prep_dir(args.workdir)  # NOTE it's kind of weird to have workdir and paired_outdir, but workdir is for stuff we throw out by default, whereas paired_outdir is for stuff we save
    work_fnames = []
    if args.infname is not None:  # if --infname is set (i.e. if there's no --paired-indir) we need to extract pair info and split loci
        prep_inputs()
    if is_subset_action(args.action):
        subset_partition(args)
        return
    args.tree_inference_outdir = getodir()
    if args.tree_inference_subdir is not None:
        args.tree_inference_outdir = '%s/%s' % (args.tree_inference_outdir, args.tree_inference_subdir)

    if args.action == 'simulate':
        if not args.dry_run:
            paircluster.prep_paired_dir(getodir(), clean=True, suffix='.yaml', extra_files=[getofn(None, trees=True, lpair=lp) for lp in utils.locus_pairs[ig_or_tr]])
        for lpair in utils.locus_pairs[ig_or_tr]:
            h_locus, l_locus = lpair
            print('%s: simulating %d event%s' % (utils.color('blue', '+'.join(lpair)), get_n_events(l_locus, args.n_sim_events), utils.plural(get_n_events(l_locus, args.n_sim_events))))
            if get_n_events(l_locus, args.n_sim_events) == 0:
                continue
            if args.input_simulation_treefname is None:
                run_step('generate-trees', h_locus, lpair=lpair)
            else:
                print('    %s drawing input simulation trees from --input-simulation-treefname *in order* for both heavy and light (rather than randomly, as for single chain)' % utils.color('yellow', 'warning'))
            for ltmp in lpair:
                run_step(args.action, ltmp, lpair=lpair)
        combine_simu_chains()
        if args.paired_outdir is None:
            work_fnames.extend(paircluster.paired_dir_fnames(getodir(), suffix='.yaml') + [get_merged_fn('all-seqs'), get_merged_fn('meta')])  # using extend() rather than += in a few places because of local variable/scoping behavior
    else:
        if args.action != 'cache-parameters' and auto_cache_params():  # ick. we can't just let the normal auto-parameter caching happen within each locus pair below, since it'll infer separately on the igh seqs paired with igk vs igl
            if args.seed_unique_id is not None:  # if we're auto parameter caching for/before seed partitioning, we *don't* (yet) want to remove non-clonal sequences, since we need all the non-clonal sequences to get better parameters (maybe at some point we want to be able to count parameters just on this lineage, but for now let's keep it simple)
                raise Exception('if setting --seed-unique-id for \'partition\', you must first explicitly run \'cache-parameters\' in order to ensure that parameters are cached on all sequences, not just clonally related sequences.')
            missing_pdir_loci = [l for l in sloci() if not os.path.exists(getpdir(l))]
            print('  missing %d locus parameter dirs (%s), so caching a new set of parameters before running action \'%s\':  %s' % (len(missing_pdir_loci), ' '.join(missing_pdir_loci), args.action, ' '.join(getpdir(l) for l in missing_pdir_loci)))
            for ltmp in sloci():
                run_step('cache-parameters', ltmp, auto_cache=True, skip_missing_input=True)
        if args.random_seed_seq:
            (args.seed_unique_id, args.seed_loci), _ = utils.choose_seed_unique_id(os.path.dirname(getifn(utils.heavy_locus)), None, None, paired=True, choose_random=True)
        ex_out_cfg_actions = ['plot-partitions', 'get-selection-metrics', 'write-fake-paired-annotations']  # actions whose behavior is determined by args.existing_output_run_cfg (I'm not sure why I didn't use runs_on_existing_output(), but maybe because these are a subset? would be nice to clean this up eventually)
        if args.action not in ['merge-paired-partitions', 'view-output']+ex_out_cfg_actions or (args.action in ex_out_cfg_actions and 'single' in args.existing_output_run_cfg):
            if args.seed_unique_id is not None and args.action == 'partition':
                remove_unseeded_seqs()
            for ltmp in sloci():
                run_step(args.action, ltmp, skip_missing_input=args.action != 'plot-partitions', skip_missing_output=runs_on_existing_output(args.action))
        if args.action in ['partition', 'merge-paired-partitions', 'annotate'] and not args.dry_run:  # ok it's weird to combine inf chains for 'annotate', but really we just want the merged/joint heavy chain output file in the normal location
            combine_inf_chains()
        ext_actions = []
        if args.action=='update-meta-info':
            ext_actions.append(args.action)
        ptnplot = args.action=='plot-partitions' or args.plot_partitions or (args.action in ['partition', 'merge-paired-partitions'] and args.plotdir is not None and not args.no_partition_plots)  # ick
        if ptnplot:
            ext_actions.append('plot-partitions')
        if args.action=='get-selection-metrics' or args.get_selection_metrics:
            ext_actions.append('get-selection-metrics')
        for atmp in ext_actions:
            process_all_outputs(atmp)  # includes several run_step() calls
        if args.action == 'view-output' or args.debug_paired_clustering or args.debug:
            for ltmp in sloci():
                run_step('view-output', ltmp, joint=True, skip_missing_input=not args.is_data, skip_missing_output=True)
        antn_pairs, fake_pntns, mpfo_lists, mtpys = None, None, None, None
        if args.action=='get-selection-metrics' or args.get_selection_metrics:
            print('%s' % utils.color('blue_bkg', 'getting combined h/l selection metrics'))
            antn_pairs, fake_pntns, mpfo_lists, mtpys = get_pntns(smetrics=True)
            treeutils.combine_selection_metrics(antn_pairs, fake_pntns, mpfo_lists, mtpys, plotdir=None if args.plotdir is None else getplotdir(None), args=args, tree_inference_outdir=getodir())
        if args.action == 'write-fake-paired-annotations':
            antn_pairs, fake_pntns, mpfo_lists, mtpys = get_pntns()
            utils.write_annotations('%s/fake-paired-annotations.yaml'%args.paired_outdir, None, fake_pntns, utils.add_lists(list(utils.annotation_headers), args.extra_annotation_columns) + utils.fake_paired_columns)
        if ptnplot and 'no-fake-paired' not in args.existing_output_run_cfg:
            from python.partitionplotter import PartitionPlotter
            print('%s' % utils.color('blue_bkg', 'plotting combined h/l partitions'))
            partplotter = PartitionPlotter(args)
            if fake_pntns is None:
                antn_pairs, fake_pntns, mpfo_lists, mtpys = get_pntns()  # smetrics=True
            inf_lines, true_lines = (None, fake_pntns) if not args.is_data else (utils.get_annotation_dict(fake_pntns), None)
            partplotter.plot('%s/partitions/%s'%(getplotdir(None), 'true' if not args.is_data else 'inferred'),
                             [l['unique_ids'] for l in (true_lines if not args.is_data else list(inf_lines.values())) if len(l['unique_ids']) > 0],
                             utils.get_annotation_dict(true_lines) if not args.is_data else inf_lines,
                             args=args)  # you get len 0 uids e.g. if reading simulation with multiple cells per droplet/missing reads

    if not args.dry_run:  # this will crash if you set plotdir, but not paired outdir, but who cares
        if args.paired_outdir is None and args.seed_unique_id is not None:
            work_fnames.extend([getodir(), os.path.dirname(getodir())])  # ick
        utils.rmdir(args.workdir, fnames=work_fnames)

# ----------------------------------------------------------------------------------------
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
class SaneArgumentParser(argparse.ArgumentParser):
  """Disables prefix matching in ArgumentParser."""  # yes, it's convenient to only need to type part of the arg, but since we have to parse through the command line strings e.g. for --paired-loci, it's extremely dangerous, and has led to bad bugs
  def _get_option_tuples(self, option_string):
    """Prevent argument parsing from looking for prefix matches."""
    return []

usage_msg = """partis <action> [options]
list available actions: partis --help
details on the options for each individual action: partis <action> --help
"""
formatter_class = MultiplyInheritedFormatter
parser = SaneArgumentParser(formatter_class=MultiplyInheritedFormatter, usage=usage_msg)
subparsers = parser.add_subparsers(dest='action')
parent_parser = SaneArgumentParser(add_help=False)

parent_args = []
parent_args.append({'name' : '--locus', 'kwargs' : {'default' : 'igh', 'choices' : utils.loci, 'help' : 'which immunoglobulin or t-cell receptor locus? This is automatically unset if --paired-loci is set.'}})
parent_args.append({'name' : '--ig-or-tr', 'kwargs' : {'default' : 'ig', 'choices' : ['ig', 'tr'], 'help' : 'if --locus is not set (i.e. if --paired-loci is set), this specifies whether we\'re running on ig or tr data (if --locus *is* set, then --ig-or-tr is set automatically). TODO probably needs a bit of testing to work for tcrs'}})  # TODO also there's probably lots of places that should use this new ig_or_tr that's in <args>
parent_args.append({'name' : '--paired-loci', 'kwargs' : {'action' : 'store_true', 'help' : 'Set this if input contains sequences from more than one locus (igh+igk+igl all together). Input can be specified either with --infname (in which case it will be automatically split apart by loci), or with --paired-indir (whose files must conform to the same conventions). It will then run the specified action on each of the single locus input files, and (if specified) merge the resulting partitions.'}})
parent_args.append({'name' : '--reverse-negative-strands', 'kwargs' : {'action' : 'store_true', 'help' : 'If --paired-loci is set, align every sequence both forwards and revcomp\'d, then for each sequence keep the sense with better alignment. If *not* running with --paired-loci, then first run bin/split-loci.py separately with --reverse-negative-strands.'}})
parent_args.append({'name' : '--dry-run', 'kwargs' : {'action' : 'store_true', 'help' : 'Just print subprocess commands that would be run without actually running them (only implemented for --paired-loci).'}})
parent_args.append({'name' : '--species', 'kwargs' : {'default' : 'human', 'choices' : ('human', 'macaque', 'mouse', 'c57bl', 'balbc'), 'help' : 'Which species?'}})
parent_args.append({'name' : '--queries', 'kwargs' : {'help' : 'Colon-separated list of query names to which to restrict the analysis'}})
parent_args.append({'name' : '--queries-to-include', 'kwargs' : {'help' : 'When reading the input file, look for and include these additional uids when --n-random-queries is set (i.e. when they otherwise might be skipped). *Not* compatible with --n-max-queries. Contrast with --queries, which includes *only* the indicated uids. Additionally, these queries are treated differently e.g. when partition plotting or removing duplicate sequences.'}})
parent_args.append({'name' : '--queries-to-include-fname', 'kwargs' : {'help' : 'In cases where you want certain sequences to be included in --queries-to-include or --seed-unique-id, but these sequences are not in --infname, you can put them in this file. Can either be fasta, or if you need to also include input meta info (e.g. timepoint, subject), use yaml (list of dicts, e.g. [{"name" : "x", "seq" : "ACGT", "affinity" : 1.3}]).'}})
parent_args.append({'name' : '--reco-ids', 'kwargs' : {'help' : 'Colon-separated list of rearrangement-event IDs to which we restrict ourselves'}})  # or recombination events
parent_args.append({'name' : '--n-max-queries', 'kwargs' : {'type' : int, 'default' : -1, 'help' : 'Maximum number of query sequences to read from input file, starting from beginning of file'}})
parent_args.append({'name' : '--n-random-queries', 'kwargs' : {'type' : int, 'help' : 'choose this many queries at random from entire input file'}})
parent_args.append({'name' : '--istartstop', 'kwargs' : {'help' : 'colon-separated start:stop line indices for input sequence file (with python slice conventions, e.g. if set to \'2:4\' will skip the zeroth and first sequences, and then take the following two sequences, and then skip all subsequence sequences). Applied before any other input filters, e.g. --n-max-queries, --queries, --reco-ids, etc.'}})

parent_args.append({'name' : '--debug', 'kwargs' : {'type' : int, 'default' : 0, 'choices' : [0, 1, 2], 'help' : 'Overall debug verbosity level. See also --sw-debug (which defaults to this), and --debug-allele-finding and --debug-paired-clustering (which are independent).'}})
parent_args.append({'name' : '--sw-debug', 'kwargs' : {'type' : int, 'choices' : [0, 1, 2], 'help' : 'Debug level for Smith-Waterman. If not set, defaults to level given by --debug.'}})
parent_args.append({'name' : '--debug-paired-clustering', 'kwargs' : {'action' : 'store_true', 'help' : 'print lots of ascii debug info while merging paired partitions'}})
parent_args.append({'name' : '--print-chosen-abs', 'kwargs' : {'action' : 'store_true', 'help' : 'print detailed ascii art information on chosen abs (see --ab-choice-cfg etc)'}})
parent_args.append({'name' : '--abbreviate', 'kwargs' : {'action' : 'store_true', 'help' : 'Abbreviate/translate sequence ids to improve readability of partition debug output. Uses a, b, c, ..., aa, ab, ...'}})
parent_args.append({'name' : '--print-git-commit', 'kwargs' : {'action' : 'store_true', 'help' : 'print sys.argv, git commit hash, and tag info'}})
parent_args.append({'name' : '--dont-write-git-info', 'kwargs' : {'action' : 'store_true', 'help' : 'Don\'t write git tag/commit info to yaml output files (used to enable diffing of test results)'}})
parent_args.append({'name' : '--n-final-clusters', 'kwargs' : {'type' : int, 'help' : 'If you reach the maximum likelihood partition and there are still more than this many clusters, attempt to keep merging until there aren\'t.  If --min-largest-cluster-size is also set, we stop if either of their criteria are satisfied. Set this also when reading existing output (e.g. plot-partitions) in order to select the proper partition.'}})
parent_args.append({'name' : '--min-largest-cluster-size', 'kwargs' : {'type' : int, 'help' : 'If you reach the maximum likelihood partition and the largest cluster isn\'t this big, attempt to keep merging until it is. If --n-final-clusters is also set, we stop if either of their criteria are satisfied. Set this also when reading existing output (e.g. plot-partitions) in order to select the proper partition.'}})
parent_args.append({'name' : '--only-print-best-partition', 'kwargs' : {'action' : 'store_true', 'help' : 'When printing annotations or reading existing output (e.g. view-output, plot-partitions), instead of the default of printing/reading the annotation for every cluster for which we calculated one (see e.g. --calculate-alternative-annotations and --n-partitions-to-write), only print annotations for clusters in the best partition.'}})
parent_args.append({'name' : '--only-print-seed-clusters', 'kwargs' : {'action' : 'store_true', 'help' : 'same as --only-print-best-partition, but in addition, only print/read the seed cluster(s). Note that if --only-print-best-partition is *not* set, then there will be more than one seed cluster.'}})
parent_args.append({'name' : '--only-print-queries-to-include-clusters', 'kwargs' : {'action' : 'store_true', 'help' : 'same as --only-print-best-partition, but in addition, only print/read the cluster(s) corresponding to --queries-to-include/--queries-to-include-fname. Note that if --only-print-best-partition is *not* set, then there will be more than one cluster for each such query.'}})  # NOTE we could just assume this if queries to include are specified, but that's not the behavior when we're plotting (in which case we include everybody but highlight the queries to include}})
parent_args.append({'name' : '--print-trees', 'kwargs' : {'action' : 'store_true', 'help' : 'When printing ascii annotations (e.g. with \'view-output\'), print any tree that\'s in each annotation.'}})
parent_args.append({'name' : '--cluster-indices', 'kwargs' : {'help' : 'indices of clusters (when sorted largest to smallest) to process for actions that read existing output, e.g. \'view-output\' and \'get-selection-metrics\'. Specified as a colon-separated list, where each item can be either a single integer or a python slice-style range of integers, e.g. 0:3-6:50 --> 0:3:4:5:50'}})
parent_args.append({'name' : '--partition-index-to-print', 'kwargs' : {'type' : int, 'help' : 'like --cluster-indices, but restricts to the partition with this index (rather than default of printing best partition).'}})
parent_args.append({'name' : '--extra-print-keys', 'kwargs' : {'help' : 'add the value of these keys (colon-separated list) as a new column to the ascii art annotation'}})

parent_args.append({'name' : '--n-procs', 'kwargs' : {'type' : int, 'default' : utils.auto_n_procs(), 'help' : 'Number of processes over which to parallelize (defaults to the number of cpus on the machine). This is usually the maximum that will be initialized at any given time, but for internal reasons, certain steps (e.g. smith waterman) sometimes use slightly more.'}})
parent_args.append({'name' : '--n-max-to-calc-per-process', 'kwargs' : {'default' : 250, 'help' : 'if a bcrham process calc\'d more than this many fwd + vtb values (and this is the first time with this number of procs), don\'t decrease the number of processes in the next step (default %(default)d)'}})
parent_args.append({'name' : '--min-hmm-step-time', 'kwargs' : {'default' : 2., 'help' : 'if a clustering step takes fewer than this many seconds, always reduce n_procs'}})
parent_args.append({'name' : '--batch-system', 'kwargs' : {'choices' : ['slurm', 'sge'], 'help' : 'batch system with which to attempt paralellization'}})
parent_args.append({'name' : '--batch-options', 'kwargs' : {'help' : 'additional options to apply to --batch-system (e.g. --batch-options : "--foo bar")'}})
parent_args.append({'name' : '--batch-config-fname', 'kwargs' : {'default' : '/etc/slurm-llnl/slurm.conf', 'help' : 'system-wide batch system configuration file name'}})  # for when you're running the whole thing within one slurm allocation, i.e. with  % salloc --nodes N ./bin/partis [...]

parent_args.append({'name' : '--only-smith-waterman', 'kwargs' : {'action' : 'store_true', 'help' : 'Exit after finishing smith-waterman.'}})
parent_args.append({'name' : '--count-parameters', 'kwargs' : {'action' : 'store_true', 'help' : 'force parameter counting when action is not cache-parameters (presumably so that you can plot them, or to get multi-hmm parameters from partitioning)'}})
parent_args.append({'name' : '--count-correlations', 'kwargs' : {'action' : 'store_true', 'help' : 'during parameter counting, also count correlations between parameters (and plot them, if --plotdir is set). NOTE that this should only be run on clustered/partitioned sequences, since otherwise it will count families more than once and won\'t mean much.'}})
parent_args.append({'name' : '--dont-write-parameters', 'kwargs' : {'action' : 'store_true', 'help' : 'don\'t write parameters to disk even if you\'ve counted them (mostly for use in conjunction with --only-smith-waterman, so you can avoid cluttering up your file system)'}})
parent_args.append({'name' : '--partis-dir', 'kwargs' : {'default' : partis_dir, 'help' : 'for internal use only'}})
parent_args.append({'name' : '--ig-sw-binary', 'kwargs' : {'default' : partis_dir + '/packages/ig-sw/src/ig_align/ig-sw', 'help' : 'Path to ig-sw executable.'}})
parent_args.append({'name' : '--vsearch-binary', 'kwargs' : { 'help' : 'Path to vsearch binary (vsearch binaries for linux and darwin are pre-installed in bin/, but for other systems you need to get your own)'}})
parent_args.append({'name' : '--is-simu', 'kwargs' : {'action' : 'store_true', 'help' : 'Set if running on simulated sequences'}})
parent_args.append({'name' : '--skip-unproductive', 'kwargs' : {'action' : 'store_true', 'help' : 'Skip sequences which Smith-Waterman determines to be unproductive (i.e. if they have stop codons, out of frame cdr3, or mutated cyst/tryp/phen)'}})
parent_args.append({'name' : '--skip-in-frame-rearrangements', 'kwargs' : {'action' : 'store_true', 'help' : 'keep only sequences that had out-of-frame rearrangements, i.e. with cdr3 length not a multiple of 3 (note that this is *different* to the in_frames key, which tells if the conserved codons are in frame with respect to the start of V, and thus which depends also on shm indels)'}})
parent_args.append({'name' : '--collapse-duplicate-sequences', 'kwargs' : {'action' : 'store_true', 'help' : 'Collapse all sequences that are identical from 5\' end of v to 3\' end of j under a single uid. The duplicate uids then appear in the output key \'duplicates\'. This collapse happens after any framework insertion/constant region trimming (see --dont-remove-framework-insertions).'}})
parent_args.append({'name' : '--also-remove-duplicate-sequences-with-different-lengths', 'kwargs' : {'action' : 'store_true', 'help' : 'If --collapse-duplicate-sequences is set, by default we collapse together only queries that have exactly the same (coding/vdj) sequence. If this argument is also set, we also collapse sequences that are sub/super strings of each other (we keep/index by the longest one). So, e.g., this would also collapse sequences that code for the same bcr but have different read lengths/coverage. NOTE this is not the same as collapsing all sequences that are identical in the parts that they share, which would be much more computationally expensive (e.g. if one is missing start of V, the other end of J).'}})
parent_args.append({'name' : '--dont-remove-framework-insertions', 'kwargs' : {'action' : 'store_true', 'help' : 'DEPRECATED see leader_seqs and c_gene_seqs columns, which are always available. Old help msg: By default we trim anything to the 5\' of V (which partis calls fv_insertion) and 3\' of J (jf_insertion), since partis ignores these regions (i.e. doesn\'t at the moment do any alignment to the leader or to constant genes), and it simplifies the output. This argument turns that off, so coordinates start from the very start of the read (even if 5\' of V), and fv and jf insertions will be set (rather than empty)'}})
parent_args.append({'name' : '--align-constant-regions', 'kwargs' : {'action' : 'store_true', 'help' : 'if your reads extend 5\' of V and/or 3\' of J, you can set this arg to align these against reference sequences for leader a C genes, which are stored in the \'leaders\' and \'c_genes\' keys.'}})
parent_args.append({'name' : '--dont-rescale-emissions', 'kwargs' : {'action' : 'store_true', 'help' : 'Don\'t scale each hmm\'s emission probabilities to account for the branch length of each individual sequence.'}})
parent_args.append({'name' : '--no-indels', 'kwargs' : {'action' : 'store_true', 'help' : 'Tell smith-waterman not to look for indels, by drastically increasing the gap-open penalty (you can also set the penalty directly).'}})
parent_args.append({'name' : '--no-indel-gap-open-penalty', 'kwargs' : {'type' : int, 'default' : 1000, 'help' : 'set --gap-open-penalty to this when --no-indels is set (also used in python/waterer.py'}})
parent_args.append({'name' : '--no-insertions-or-deletions', 'kwargs' : {'action' : 'store_true', 'help' : 'Forbid (real, not effective) insertions and deletions in simulation and in hmm (but not sw). Removes all (real) insertions and deletions when reading hmm output (note that they\'re still allowed in smith waterman annotations). Useful as an ad hoc way to enforce a particular naive sequence: split your naive sequence into v/d/j parts, put them into a germline set directory, and set --initial-germline-dir and --no-insertions-or-deletions). To forbid indels in smith-waterman, see --no-indels'}})
parent_args.append({'name' : '--random-seed', 'kwargs' : {'type' : int, 'default' : int(time.time()), 'help' : 'Random seed used by many different things, but especially when reshuffling sequences between partitioning steps, and during simulation. Set this if you want to get exactly the same result when rerunning.'}})
parent_args.append({'name' : '--seed', 'kwargs' : {'type' : int, 'help' : 'DEPRECATED use --random-seed'}})
parent_args.append({'name' : '--min-observations-per-gene', 'kwargs' : {'type' : int, 'default' : 20, 'help' : 'If a gene is observed fewer times than this, we average over other alleles/primary versions/etc. (e.g. in recombinator and hmmwriter). Also used as a more general "this isn\'t very many observations" value.'}})
parent_args.append({'name' : '--no-per-base-mfreqs', 'kwargs' : {'action' : 'store_true', 'help' : 'When making the HMM model files, instead of the default of accounting for different rates to different bases (i.e. A->T vs A->G), do *not* account for the different rates to different bases. This is only really useful for testing the new simulation option --per-base-mutation.'}})
parent_args.append({'name' : '--region-end-exclusion-length', 'kwargs' : {'type' : int, 'default' : 0, 'help' : 'when counting/writing parameters, ignore this many bases abutting non-templated insertions for calculating mutation frequencies (note: doesn\'t make a difference (hence set to 0 by default) probably because we\'re setting a strongish prior on these bases when writing hmms anyway'}})
parent_args.append({'name' : '--allow-conserved-codon-deletion', 'kwargs' : {'action' : 'store_true', 'help' : 'When building hmm yaml model files during parameter caching, allow deletions that extend through the conserved codons (cyst in V and tryp/phen in J) (by default such deletions are forbidden; see https://github.com/psathyrella/partis/issues/308). NOTE that this has *no* effect if you\'ve already cached parameters.'}})
parent_args.append({'name' : '--subcluster-annotation-size', 'kwargs' : {'default' : 3, 'help' : 'when running the bcrham viterbi algorithm, instead of running the multi-hmm on the entire family, split the family into subclusters of (about) this size, then replace each subcluster with its naive sequence and iterate until running on one cluster of (about) this size consisting entirely of inferred naive/ancestor sequences (see https://github.com/psathyrella/partis/issues/308). Set to the string \'None\' to turn off.'}})

parent_args.append({'name' : '--only-genes', 'kwargs' : {'help' : 'Colon-separated list of genes to which to restrict the analysis. If any regions (V/D/J) are not represented among these genes, these regions are left unrestricted. If running \'simulate\', you probably also want to set --force-dont-generate-germline-set.'}})
parent_args.append({'name' : '--n-max-per-region', 'kwargs' : {'default' : '3:5:2', 'help' : 'Number of best smith-waterman matches (per region, in the format v:d:j) to pass on to the hmm. Note that regardless of these values, s-w only writes hmm model files corresponding to the best match for each sequence.'}})
parent_args.append({'name' : '--gap-open-penalty', 'kwargs' : {'type' : int, 'default' : 30, 'help' : 'Penalty for indel creation in Smith-Waterman step.'}})
parent_args.append({'name' : '--max-vj-mut-freq', 'kwargs' : {'type' : float, 'default' : 0.4, 'help' : 'skip sequences whose mutation rates in V and J are greater than this (it\'s really not possible to get meaningful smith-waterman matches above this)'}})
parent_args.append({'name' : '--max-logprob-drop', 'kwargs' : {'type' : float, 'default' : 5., 'help' : 'stop glomerating when the total logprob has dropped by this much'}})
parent_args.append({'name' : '--n-simultaneous-seqs', 'kwargs' : {'type' : int, 'help' : 'Number of simultaneous sequences on which to run the multi-HMM (e.g. 2 for a pair hmm)'}})
parent_args.append({'name' : '--all-seqs-simultaneous', 'kwargs' : {'action' : 'store_true', 'help' : 'Run all input sequences simultaneously, i.e. equivalent to setting --n-simultaneous-seqs to the number of input sequences.'}})
parent_args.append({'name' : '--simultaneous-true-clonal-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'If action is annotate/cache-parameters, run true clonal sequences together simultaneously with the multi-HMM. If actions is partition, don\'t actually run clustering and instead use the true partition (useful for e.g. validating selection metrics, where you don\'t want to be conflating partition performance with selection metric performance). For search: use-true-partition'}})
parent_args.append({'name' : '--mimic-data-read-length', 'kwargs' : {'action' : 'store_true', 'help' : 'In simulation, trim V 5\' and D 3\' to mimic read lengths seen in data (must also be set when caching parameters)'}})

parent_args.append({'name' : '--infname', 'kwargs' : {'help' : 'input sequence file in .fa, .fq, .csv, or partis output .yaml (if .csv, specify id string and sequence headers with --name-column and --seq-column)'}})
parent_args.append({'name' : '--paired-indir', 'kwargs' : {'help' : 'Directory with input files for use with --paired-loci. Must conform to file naming conventions from bin/split-loci.py (really paircluster.paired_dir_fnames()), i.e. the files generated when --infname and --paired-loci are set.'}})
parent_args.append({'name' : '--guess-pairing-info', 'kwargs' : {'action' : 'store_true', 'help' : utils.did_help['guess']}})
parent_args.append({'name' : '--no-pairing-info', 'kwargs' : {'action' : 'store_true', 'help' : 'don\'t try to extract pairing info even though --paired-loci is set (useful if the sequence ids will confuse the pair info extraction code)'}})
parent_args.append({'name' : '--droplet-id-separators', 'kwargs' : {'help' : utils.did_help['seps']}})
parent_args.append({'name' : '--droplet-id-indices', 'kwargs' : {'help' : utils.did_help['indices']}})
parent_args.append({'name' : '--name-column', 'kwargs' : {'help' : 'column/key name for sequence ids in input csv/yaml file (default: \'unique_ids\'). If set to \'fasta-info-index-N\' for an integer N, it will take the Nth (zero-indexed) value from a fasta files uid line.'}})  # for search: includes what used to be --fasta-info-index
parent_args.append({'name' : '--seq-column', 'kwargs' : {'help' : 'column/key name for nucleotide sequences in input csv/yaml file (default: \'input_seqs\')'}})
parent_args.append({'name' : '--sanitize-input-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'If set, replaces unexpected characters in input sequences with supported alternatives, e.g. R --> N (only from fasta though).'}})
parent_args.append({'name' : '--input-metafnames', 'kwargs' : {'help' : 'colon-separated list of yaml/json files, each of which has meta information for the sequences in --infname (and --queries-to-include-fname, although that file can also include its own input meta info), keyed by sequence id. If running multiple steps (e.g. cache-parameters and partition), this must be set for all steps. See https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#input-meta-info for an example.'}})
parent_args.append({'name' : '--input-metafname', 'kwargs' : {'help' : 'DEPRECATED use --input-metafnames'}})
parent_args.append({'name' : '--ignore-default-input-metafile', 'kwargs' : {'action' : 'store_true', 'help' : 'If set, do not add input meta file from default location (e.g. in --paired-indir) to --input-metafnames. Necessary e.g. for subset-partition to avoid re-reading uncleaned pair info during final merged partition step.'}})
parent_args.append({'name' : '--input-partition-fname', 'kwargs' : {'help' : 'partis-style json/yaml file from which to read a partition for input to specified action. For \'annotate\' action, annotates sequences in --infname using the partition in this file (rather than the default of annotating each sequence individually). For \'partition\' action, by default uses the input partition instead of doing any actual partitioning (to continue partitioning starting from the input partition, set --continue-from-input-partition). Note that if you\'re using this with --paired-loci, you should specify the entire directory, rather than a single file.'}})
parent_args.append({'name' : '--input-partition-index', 'kwargs' : {'type' : int, 'help' : 'Index of the partition to be read from --input-partition-fname (if unset, defaults to the best partition). To figure out which index you want, you probably want to run the view-output action on the file.'}})
parent_args.append({'name' : '--outfname', 'kwargs' : {'help' : 'output file name'}})
parent_args.append({'name' : '--paired-outdir', 'kwargs' : {'help' : 'Directory for all output files when --paired-loci is set, i.e. involving multiple loci in input and/or paired heavy/light information.'}})
parent_args.append({'name' : '--dont-calculate-annotations', 'kwargs' : {'action' : 'store_true', 'help' : 'Don\'t calculate annotations for the final partition (so just the partition is written to output).'}})
parent_args.append({'name' : '--use-sw-annotations', 'kwargs' : {'action' : 'store_true', 'help' : 'Instead of running hmm annotation on each cluster in the final partition, use the smith-waterman annotations to synthesize a multi-sequence annotation for each final cluster.'}})
parent_args.append({'name' : '--write-full-yaml-output', 'kwargs' : {'action' : 'store_true', 'help' : 'By default, we write yaml output files using the json subset of yaml, since it\'s much faster. If this is set, we instead write full yaml, which is more human-readable (but also much slower).'}})
parent_args.append({'name' : '--presto-output', 'kwargs' : {'action' : 'store_true', 'help' : 'Write output file(s) in presto/changeo format. Since this format depends on a particular IMGT alignment, this depends on a fasta file with imgt-gapped alignments for all the V, D, and J germline genes. The default in data/germlines/<species>/imgt-alignments/, is probably fine for most cases. For the \'annotate\' action, a single .tsv file is written with annotations (so --outfname suffix must be .tsv). For the \'partition\' action, a fasta file is written with cluster information (so --outfname suffix must be .fa or .fasta), as well as a .tsv in the same directory with the corresponding annotations.'}})
parent_args.append({'name' : '--airr-output', 'kwargs' : {'action' : 'store_true', 'help' : 'Write output file(s) in AIRR-C format (if --outfname has suffix .tsv, only the airr .tsv is written; however if --outfname has suffix .yaml, both the standard partis .yaml file and an airr .tsv are written). A description of the airr columns can be found here https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields.'}})
parent_args.append({'name' : '--airr-input', 'kwargs' : {'action' : 'store_true', 'help' : 'Read --infname in airr format. Equivalent to setting \'--seq-column sequence --name-column sequence_id\'.'}})
parent_args.append({'name' : '--extra-annotation-columns', 'kwargs' : {'help' : 'Extra columns to add to the (fairly minimal) set of information written by default to annotation output files (choose from: %s)' % ' '.join(utils.extra_annotation_headers)}})  # NOTE '-columns' in command line arg, but '-headers' in utils (it's more consistent that way, I swear}})
parent_args.append({'name' : '--cluster-annotation-fname', 'kwargs' : {'help' : 'output file for cluster annotations (default is <--outfname>-cluster-annotations.<suffix>)'}})
parent_args.append({'name' : '--parameter-dir', 'kwargs' : {'help' : 'Directory to/from which to write/read sample-specific parameters. If not specified, a default location is used (and printed to std out). If it does not exist, we infer parameters before proceeding to the desired action.'}})
parent_args.append({'name' : '--parameter-type', 'kwargs' : {'default' : None, 'choices' : utils.parameter_type_choices, 'help' : 'Use parameters from Smith-Waterman (sw), single-sequence HMM (hmm), or multi-sequence HMM (multi-hmm) subdirectories for inference and simulation? \'multi-hmm\' is best, but requires running \'partition\' with --count-parameters set; \'hmm\' is much better than \'sw\', but \'sw\' is occasionally useful for debugging. If not set, looks first for \'multi-hmm\', then \'hmm\', then \'sw\', if none are found defaults to %s. You shouldn\'t generally need to set this unless you want to do something somewhat weird.'%utils.default_parameter_type}})
parent_args.append({'name' : '--parameter-out-dir', 'kwargs' : {'help' : 'Special parameter dir for writing multi-hmm parameters, i.e. when running annotate or partition with --count-parameters set (if not set, defaults to <--parameter-dir>/multi-hmm).'}})
parent_args.append({'name' : '--refuse-to-cache-parameters', 'kwargs' : {'action' : 'store_true', 'help' : 'Disables auto parameter caching, i.e. if --parameter-dir doesn\'t exist, instead of inferring parameters, raise an exception. Useful for batch/production use where you want to make sure you\'re caching parameters in a separate step.'}})
parent_args.append({'name' : '--persistent-cachefname', 'kwargs' : {'help' : 'Name of file which will be used as an initial cache file (if it exists), and to which all cached info will be written out before exiting. Must be set to \'paired-outdir\' if --paired-loci is set.'}})
parent_args.append({'name' : '--sw-cachefname', 'kwargs' : {'help' : 'Smith-Waterman cache file name. Default is set using a hash of all the input sequence ids (in partitiondriver, since we have to read the input file first).'}})
parent_args.append({'name' : '--write-sw-cachefile', 'kwargs' : {'action' : 'store_true', 'help' : 'Write sw results to the sw cache file during actions for which we\'d normally only look for an existing one (i.e annotate and partition).'}})
parent_args.append({'name' : '--workdir', 'kwargs' : {'help' : 'Temporary working directory (default is set below)'}})

parent_args.append({'name' : '--plotdir', 'kwargs' : {'help' : 'Base directory to which to write plots (by default this is not set, and consequently no plots are written. If --paired-loci is set, you can set --plotdir to \'paired-outdir\' to copy over the value from that argument.'}})
parent_args.append({'name' : '--plot-annotation-performance', 'kwargs' : {'action' : 'store_true', 'help' : 'Compare inferred and true annotation values (deletion lengths, gene calls, etc.). If --plotdir is set, write corresponding plots to disk, and see also --print-n-worst-annotations.'}})
parent_args.append({'name' : '--plot-partitions', 'kwargs' : {'action' : 'store_true', 'help' : 'after running annotation, plot partitions on resulting cluster annotations. This only makes sense if partitions were specified with --input-partition-fname, and only with \'annotate\' and \'merge-paired-partitions\' actions. Can also use \'plot-partitions\' action'}})
parent_args.append({'name' : '--partition-plot-cfg', 'kwargs' : {'default' : ':'.join(utils.default_ptn_plot_cfg), 'help' : 'colon-separated list of plot types to make for partition plots. Choices: %s.' % ' '.join(utils.all_ptn_plot_cfg)}})
parent_args.append({'name' : '--print-n-worst-annotations', 'kwargs' : {'type' : int, 'help' : 'For use with --plot-annotation-performance: print ascii annotations for the N least accurate annotations according to several annotation values (e.g. VD insertion length, distance to true naive sequence). One of two ways to visualize annotation performance -- the other is by setting --plotdir to look at summary plots over all families. View with "less -RS".'}})
parent_args.append({'name' : '--no-partition-plots', 'kwargs' : {'action' : 'store_true', 'help' : 'don\'t make paritition plots, even if --plotdir is set (presumably because you want other plots, e.g. for selection metrics -- the partition plots have a lot more depenencies)'}})
parent_args.append({'name' : '--only-csv-plots', 'kwargs' : {'action' : 'store_true', 'help' : 'skip writing actual image files, which can quite be slow, and only write the csv/yaml summaries (where implemented)'}})
parent_args.append({'name' : '--make-per-gene-plots', 'kwargs' : {'action' : 'store_true', 'help' : 'in addition to plots aggregating over genes, write plots displaying info for each gene of, e.g., per position shm rate, deletion frequencies'}})
parent_args.append({'name' : '--make-per-gene-per-base-plots', 'kwargs' : {'action' : 'store_true', 'help' : 'in addition to the plots made by --make-per-gene-plots, also make the per-gene, per-base plots (i.e. showing A->T vs A->G (this is quite slow, like a few seconds per gene plot).'}})
parent_args.append({'name' : '--linearham-dir', 'kwargs' : {'default' : ('%s/work/linearham' % os.getenv('HOME')) if os.getenv('HOME') is not None else None, 'help' : 'path to linearham main dir (necessary if you want to use linearham without docker)'}})
parent_args.append({'name' : '--meta-info-to-emphasize', 'kwargs' : {'help' : 'Input meta info (or regular annotation) key to emphasize (highlight in red) in various plots, similar to --queries-to-include. Specify as comma-separated key-value pair, for instance \'timepoints,+8d\' would highlight all sequences with timepoint \'+8d\'. Can be any annotation key or input meta key. For now only supports one key-val pair, but in future should support colon-separated list.'}})
parent_args.append({'name' : '--meta-info-key-to-color', 'kwargs' : {'help' : 'Like --meta-info-to-emphasize, except for this key we choose a different color for each value.'}})
parent_args.append({'name' : '--meta-info-bkg-vals', 'kwargs' : {'help' : 'If set, these values for --meta-info-key-to-color will be deemphasized (e.g. smaller points are used). Use e.g. if the vast majority of your points have these values, so you want to highlight other values.'}})
parent_args.append({'name' : '--meta-emph-formats', 'kwargs' : {'help' : 'Formatting instructions for keys used in --meta-info-to-emphasize and --meta-info-key-to-color, specified as a colon-separated list of comma-separated key/val pairs. Currently you can do three things: specify a different (presumably prettier) name for a key, specify a different name for a value, and specify to emphasize/color based on length rather than value. For instance \'--meta-info-key-to-color multiplicities --meta-info-to-emphasize paired-uids,1 --meta-emph-formats multiplicities,multi:paired-uids,len\' will color by multiplicities, but label this in the plots with \'multi\', and emphasize seqs with exactly one paired uid (i.e. their \'paired-uids\' list is length 1). To specify a better name for a value, use =, e.g. \'specificitys,yes=HSV+\' will use \'HSV+\' in place of \'yes\' for the key \'specificitys\''}})
parent_args.append({'name' : '--existing-output-run-cfg', 'kwargs' : {'help' : 'colon-separated list specifying the paired clustering results on which to run "existing output" actions, e.g. get-selection-metrics and plot-partitions. \'single\': single chain results, \'paired\': paired results in e.g. h+k/ subdirs, \'merged\': merged igh results (does *not* run on link to e.g. h+k/ results). If not set, only runs on the \'fake paired\' annotations resulting from smashing together each pair of h+l seqs (which are also run, unless \'no-fake-paired\' is in this arg).'}})
parent_args.append({'name' : '--cluster-size-bins', 'kwargs' : {'help' : 'colon-separated list of bin low edges for cluster size plots when partition plotting'}})
parent_args.append({'name' : '--dont-label-queries-to-include', 'kwargs' : {'action' : 'store_true', 'help' : 'when plotting --queries-to-include queries (and --seed-unique-id), neither highlight nor label them. I.e. these queries are handled as normal during partitioning, but treated like non-emphasized queries when plotting (the idea is that if you have a very large number of queries that you want to include, you don\'t necessarily want them emphasized in plots. Not yet implemented for all plot types.'}})

parent_args.append({'name' : '--default-initial-germline-dir', 'kwargs' : {'default' : partis_dir + '/data/germlines', 'help' : 'For internal use only (used to pass around the default location). To specify your own germline directory from which to start, use --initial-germline-dir instead.'}})
parent_args.append({'name' : '--initial-germline-dir', 'kwargs' : {'help' : 'Directory from which to read initial germline info (either from which to begin inference, or with which to simulate). For inference, it is only actually read when caching parameters, while subsequent steps read the inferred germline set from --parameter-dir. NOTE default is set below, because control flow is too complicated for argparse. See also --simulation-germline-dir.'}})
parent_args.append({'name' : '--sanitize-input-germlines', 'kwargs' : {'action' : 'store_true', 'help' : 'By default we require all gene names to look like IGHVx*y. If you set this option when you\'re also setting --initial-germline-dir, then we instead allow arbitrary strings as gene names in this input germline directory, and we sanitize them by adding the correct locus and region info at the start (so make sure that you have --locus set correctly), and adding an arbitrary allele string at the end (e.g. <stuff> --> IGHV<stuff>*x).'}})
parent_args.append({'name' : '--simulation-germline-dir', 'kwargs' : {'help' : 'Location of germline info that was used for simulation. Not generally necessary, since the germline set that was actually used for simulation is written to the output yaml file, but if using deprecated csv output files, or sometimes when evaluating germline inference, it\'s necessary to have the original, initial germline info that was passed to simulation.'}})
parent_args.append({'name' : '--dont-add-simu-implicit-info', 'kwargs' : {'action' : 'store_true', 'help' : 'when running on simulation without complete annotation information, this lets you turn off implicit info adding, which means things may still work using only the annotation info that\'s there'}})
parent_args.append({'name' : '--aligned-germline-fname', 'kwargs' : {'help' : 'fasta file with imgt-gapped alignments for each V, D, and J gene (used by --presto-output). The defaults are in data/germlines/<species>/imgt-alignments/ (set in processargs.py). The existing alignments are automatically extended with any missing genes, but no guarantees are made as to the similarity of this to what imgt would do.'}})
parent_args.append({'name' : '--dtr-path', 'kwargs' : {'help' : 'path to directory with decision tree regression model files'}})  # , 'default' : partis_dir + '/data/selection-metrics/dtr-models'

parent_args.append({'name' : '--n-max-alleles-per-gene', 'kwargs' : {'type' : int, 'default' : None, 'help' : 'if set, after allele inference the germline set is reduced such that each imgt-name-defined gene has at most this many alleles, with the alleles assigned to the fewest sequences removed first.'}})
parent_args.append({'name' : '--typical-genes-per-region-per-subject', 'kwargs' : {'default' : '55:25:6', 'help' : 'typical number of alleles per subject for the v, d, and j regions (used by --min-allele-prevalence-fraction)'}})
parent_args.append({'name' : '--min-allele-prevalence-fraction', 'kwargs' : {'type' : float, 'default' : 0.0005, 'help' : 'Remove any V alleles that represent less than this fraction of the repertoire (rescaled using --typical-genes-per-region-per-subject for D and J). Converted to --min-allele-prevalence-fractions (note plural!) internally.'}})

parent_args.append({'name' : '--leave-default-germline', 'kwargs' : {'action' : 'store_true', 'help' : 'Turns off all germline inference, i.e. just assigns each sequence to its closest match in the initial germline set. Generally this will result in many spurious/misassigned genes, and usually makes smith-waterman much slower (since it has to align against way more genes), but it does reduce the number of smith-waterman steps that have to be run. Equivalent to setting --dont-remove-unlikely-alleles and --dont-find-new-alleles, as well as ensuring that --allele-cluster is off.'}})
parent_args.append({'name' : '--dont-remove-unlikely-alleles', 'kwargs' : {'action' : 'store_true', 'help' : 'Turn off the allele-removal step of germline inference (see --leave-default-germline).'}})
parent_args.append({'name' : '--allele-cluster', 'kwargs' : {'action' : 'store_true', 'help' : 'Turn on clustering-based germline inference. This is automatically turned on for non-human species. To turn this off, you either have to set --leave-default-germline, or edit python/processargs.py.'}})
parent_args.append({'name' : '--kmeans-allele-cluster', 'kwargs' : {'action' : 'store_true', 'help' : 'it\'s possible (even likely!) that a kmeans-style clustering approach would work better for clustering-based germline inference, but it isn\'t fully implemented at the moment. Nonetheless, this will turn it on.'}})
parent_args.append({'name' : '--dont-find-new-alleles', 'kwargs' : {'action' : 'store_true', 'help' : 'Turn off fit-based germline inference.'}})
# parent_args.append({'name' : '--always-find-new-alleles', 'kwargs' : {'action' : 'store_true', 'help' : 'By default we only look for new alleles if a repertoire\'s mutation rate is amenable to reasonable new-allele sensitivity (i.e. if it\'s not crazy high). This overrides that.'}})
parent_args.append({'name' : '--debug-allele-finding', 'kwargs' : {'action' : 'store_true', 'help' : 'print lots of debug info on new-allele fits'}})
parent_args.append({'name' : '--new-allele-fname', 'kwargs' : {'help' : 'fasta fname to which to write any new alleles (they are also written, together with previously-known alleles that are also in the sample, to --parameter-dir)'}})
parent_args.append({'name' : '--n-max-snps', 'kwargs' : {'type' : int, 'default' : 8, 'help' : 'when new-allele finding, look for new V alleles separated from existing V alleles by up to this many SNPs. Also used for allele removal (corresponding numbers for D and J are set automatically)'}})
parent_args.append({'name' : '--n-max-mutations-per-segment', 'kwargs' : {'type' : int, 'default' : 23, 'help' : 'when new-allele finding, exclude sequences which have more than this many mutations in the V segment'}})
parent_args.append({'name' : '--min-allele-finding-gene-length', 'kwargs' : {'type' : int, 'default' : 150, 'help' : 'if (after excluding particularly short reads) the reads for a gene are shorter than this, then don\'t look for new alleles with/on this gene'}})
parent_args.append({'name' : '--plot-and-fit-absolutely-everything', 'kwargs' : {'type' : int, 'help' : 'fit every single position for this <istart> and write every single corresponding plot (slow as hell, and only for debugging/making plots for paper)'}})
parent_args.append({'name' : '--get-selection-metrics', 'kwargs' : {'action' : 'store_true', 'help' : 'calculate selection metrics for each cluster (can also be run on existing output with \'get-selection-metrics\' action). By default it calculates fast, but very approximate, trees (see manual for details).'}})
parent_args.append({'name' : '--min-selection-metric-cluster-size', 'kwargs' : {'type' : int, 'default' : treeutils.default_min_selection_metric_cluster_size, 'help' : 'When calculating selection metrics, ignore clusters smaller than this. See also --min-paired-cluster-size-to-read, which is similar but applies earlier, when reading clusters from files.'}})
parent_args.append({'name' : '--min-paired-cluster-size-to-read', 'kwargs' : {'type' : int, 'default' : treeutils.default_min_selection_metric_cluster_size, 'help' : 'When reading paired annotations when getting selection metrics or plotting partitions, ignore clusters with either N h or l ids smaller than this. See also --min-selection-metric-cluster-size, which is similar but only skips selection metric calculation.'}})
parent_args.append({'name' : '--treefname', 'kwargs' : {'help' : 'newick-formatted file with a tree corresponding to the sequences either in --infname (if making new output, i.e. action is annotate or partition) or --outfname (if reading existing output, i.e. action is get-selection-metrics) (unrelated to --input-simulation-treefname).'}})
parent_args.append({'name' : '--tree-inference-method', 'kwargs' : {'choices' : ['fasttree', 'iqtree', 'iqtree-1.6.beta3', 'iqtree-2.3.1', 'raxml', 'gctree', 'gctree-base', 'gctree-mut-mult', 'gctree-no-dag', 'linearham', 'igphyml', 'cpath'], 'help' : 'Method to use when inferring trees (default: fasttree)'}})
parent_args.append({'name' : '--tree-inference-subdir', 'kwargs' : {'help' : 'Subdirectory of the (automatically-set) tree inference workdir to which to write tree inference output files. By default, these are written to a subdir of --outfname/--paired-outdir with the name of the inference method, but this argument allows multiple versions, as with --sub-plotdir.'}})
parent_args.append({'name' : '--infer-trees-with-only-leaves', 'kwargs' : {'action' : 'store_true', 'help' : 'Discard internal nodes in true trees when inferring phylogenetic trees.'}})
parent_args.append({'name' : '--infer-trees-with-collapsed-duplicate-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'collapse identical sequences when inferring trees (and then again when plotting them, since unfortunately it\'s really hard to pass the information between those two scripts). Yes, it sucks having both this and --collapse-duplicate-sequences (which acts during smith-waterman), but sometimes we can\'t use that, e.g. with queries to include or subset partitioning (at least without more work).'}})
parent_args.append({'name' : '--phylo-naive-inference-fuzz', 'kwargs' : {'type' : int, 'help' : 'when inferring trees, replace the original inferred naive seq with the inferred root from the phylo program. When passing the naive seq/root to the phylo program, we replace uncertain positions (e.g. non-templated insertions) with Ns in the naive sequence before passing to tree inference. The phylo-inferred naive sequence is added to the annotation as an observed sequence named \'phylo-naive\'.'}})
parent_args.append({'name' : '--selection-metrics-to-calculate', 'kwargs' : {'default' : 'lbi:lbr:aa-lbi:aa-lbr:cons-dist-aa', 'help' : 'colon-separated list of selection metrics to calculate.'}})
parent_args.append({'name' : '--selection-metric-plot-cfg', 'kwargs' : {'default' : ':'.join(treeutils.default_plot_cfg), 'help' : 'colon-separated list of plot types to make for selection metrics (see treeutils.plot_tree_metrics()). Choices: %s.' % ' '.join(treeutils.all_plot_cfg)}})
parent_args.append({'name' : '--extra-daffy-metrics', 'kwargs' : {'help' : 'colon-separated list of extra metrics (beyond defaults of %s) for which to make delta-affinity performance plots'%':'.join(treeutils.daffy_metrics)}})
parent_args.append({'name' : '--sub-plotdir', 'kwargs' : {'help' : 'A subdir to which to write partition and selection metric plots (e.g. if plotting several different tree inference methods or --affinity-key[s]). Useful for instance with paired loci where you can\'t always control plot dirs for individual commands.'}})
parent_args.append({'name' : '--slice-bin-fname', 'kwargs' : {'default' : '%s/data/selection-metrics/slice-bins.yaml'%partis_dir, 'help' : 'yaml file with N bins for each slice variable in the selection metric sliced specificity correlation plots (see default for example)'}})
parent_args.append({'name' : '--selection-metric-fname', 'kwargs' : {'help' : 'yaml file to which to write selection metrics. If not set, and if --add-metrics-to-outfname is not set, this defaults to <--outfname>.replace(<suffix>, \'-selection-metrics\' + <suffix>)'}})
parent_args.append({'name' : '--add-selection-metrics-to-outfname', 'kwargs' : {'action' : 'store_true', 'help' : 'If set, instead of writing a separate file with selection metrics, we include them in <--outfname> under the key \'tree-info\'.'}})
parent_args.append({'name' : '--lb-tau', 'kwargs' : {'type' : float, 'help' : 'exponential decay length for local branching index (lbi).'}})
parent_args.append({'name' : '--dont-normalize-lbi', 'kwargs' : {'action' : 'store_true', 'help' : 'By default, we normalize tau so that 0. is a rough minimum, and 1. is a fairly large value. But since you can\'t normalize for tau larger than 1/seq_len, it\'s useful to turn off normalization.'}})
parent_args.append({'name' : '--include-relative-affy-plots', 'kwargs' : {'action' : 'store_true', 'help' : 'in addition to validation plots using actual affinity from simulation, also make plots using \'relative\' affinity (see bcr-phylo), i.e. a cell\'s affinity relative only to other cells alive at the same point in time. Useful because the selection metrics attempt to predict the effects of fitness on branchiness only at a given time (or: if you sample a lot of intermediate ancestors, your selection metric performance will look artificially bad because the ancestors, while having higher affinity than other cells alive at the same time, have affinity lower than the leaves).'}})
parent_args.append({'name' : '--affinity-key', 'kwargs' : {'help' : 'Name of output key/column to use instead of \'affinities\' in selection metric plotting.'}})
parent_args.append({'name' : '--invert-affinity', 'kwargs' : {'action' : 'store_true', 'help' : 'take the multiplicative inverse of --affinity-key (which must be set if using this arg) before making selection metric plots (e.g. for ic50)'}})
parent_args.append({'name' : '--use-droplet-id-for-combo-id', 'kwargs' : {'action' : 'store_true', 'help' : 'when combining h/l selection metrics, this forces to only use the droplet id (rather than also adding contig ids).'}})
parent_args.append({'name' : '--ab-choice-cfg', 'kwargs' : {'default' : partis_dir + '/data/selection-metrics/ab-choice.yaml', 'help' : 'yaml file with config info for how to choose abs.'}})
parent_args.append({'name' : '--choose-all-abs', 'kwargs' : {'action' : 'store_true', 'help' : 'instead of reading --ab-choice-cfg to decide which abs to choose, take them all (presumably because you just want info for all seqs dumped to --chosen-ab-fname).'}})
parent_args.append({'name' : '--label-tree-nodes', 'kwargs' : {'action' : 'store_true', 'help' : 'label all nodes in the tree plots (NOTE that you may need to also set a tree plot color args, e.g. --meta-info-key-to-color, in order for the node labels to display)'}})
parent_args.append({'name' : '--label-leaf-nodes', 'kwargs' : {'action' : 'store_true', 'help' : 'label only leaf nodes in the tree plots (NOTE that you may need to also set a tree plot color args, e.g. --meta-info-key-to-color, in order for the node labels to display)'}})
parent_args.append({'name' : '--label-root-node', 'kwargs' : {'action' : 'store_true', 'help' : 'label root node in the tree plots (NOTE that you may need to also set a tree plot color args, e.g. --meta-info-key-to-color, in order for the node labels to display)'}})
parent_args.append({'name' : '--node-size-key', 'kwargs' : {'help' : 'similar to --meta-info-key-to-color, but this controls the size of the node in tree plots'}})
parent_args.append({'name' : '--mutation-label-cfg', 'kwargs' : {'help' : 'colon-separated list of strings specifying how to label nodes in tree plots. \'all\': label all branches with mutations along that branch (like \'3 nuc, 1 aa\'), \'leaf\': only label terminal branches (with total mutations from root), \'mut-strs\': include str for each individual AA mutation (e.g. N54M), \'short\': add only the number of nuc and aa mutations (like \'3 (1)\').'}})
parent_args.append({'name' : '--branch-color-key', 'kwargs' : {'help' : 'see bin/plot-lb-trees.py help'}})
parent_args.append({'name' : '--node-label-regex', 'kwargs' : {'help' : 'see bin/plot-lb-trees.py help'}})
parent_args.append({'name' : '--min-face-size', 'kwargs' : {'type' : int, 'help' : 'see bin/plot-lb-trees.py help'}})
parent_args.append({'name' : '--max-face-size', 'kwargs' : {'type' : int, 'help' : 'see bin/plot-lb-trees.py help'}})

parent_args.append({'name' : '--no-sw-vsearch', 'kwargs' : {'action' : 'store_true', 'help' : 'By default, we get preliminary V annotations from vsearch that we pass to sw. This improves accuracy because the vsearch shm rate estimates let us separate sequences into groups by optimal sw match/mismatch scores, and improves speed by letting us reverse shm V indels before running sw (so fewer sw iterations; also vsearch is faster than sw). Setting this option skips the preliminary vsearch step (although if caching parameters, we still use vsearch to remove less-likely genes; but no vsearch info is passed to sw)'}})

parent_args.append({'name' : '--collapse-similar-paired-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'if set, when cleaning heavy/light pair info we collapse highly similar sequences within droplets (i.e. extra seqs that likely result from sequencing error).'}})
parent_args.append({'name' : '--pair-unpaired-seqs-with-paired-family', 'kwargs' : {'action' : 'store_true', 'help' : 'After finishing normal paired clustering, continue by pairing any unpaired sequences in each heavy/light family pair with the partner of their nearest paired family member. This is designed for use when you have both single cell and bulk samples, and you want to use the single seq pairing info to arrive at some approximate pair info for the bulk sample.'}})
parent_args.append({'name' : '--add-unpaired-seqs-to-fake-paired-annotations', 'kwargs' : {'action' : 'store_true', 'help' : 'when making the fake h+l "paired" annotation (i.e. smashing h+l seqs together) for use in selection metrics and/or partition plotting, by default we ignore unpaired seqs. This option includes them, with Ns for the missing opposite-chain sequence.'}})
parent_args.append({'name' : '--keep-all-unpaired-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'By default when paired clustering, seqs with no pair info are kept only if they\'re in a (single-chain) family with at least one paired seq (i.e. families consisting entirely of unpaired seqs are discarded). If this is set, instead we keep all unpaired seqs.'}})
parent_args.append({'name' : '--ignore-sw-pair-info', 'kwargs' : {'action' : 'store_true', 'help' : 'If we already have paired clustering results and we\'re reading them for a subsequent step (for instance to get annotations) that is also reading sw cache info, we usually want to ignore the pair info in the sw cache file, since it is pre-cleaning/uncorrected. This argument instructs it to do that. Note that even if this argument isn\'t set, it should correctly figure out to overwrite the sw pair info with the correct pair info, but it\'s better/safer to also not read the sw pair info to begin with.'}})

parent_args.append({'name' : '--max-ccf-fail-frac', 'kwargs' : {'type' : float, 'default' : 0.05, 'help' : 'when calculating clustering performance metrics (correct cluster fractions, purity/completeness), crash if more than this fraction of sequences are missing from the inferred partition'}})

parent_args.append({'name' : '--prepend-coverage-command', 'kwargs' : {'action' : 'store_true', 'help' : 'when running subprocesses, prepend the python coverage append command to their command lines'}})

# ----------------------------------------------------------------------------------------
subconfig = collections.OrderedDict((
    # NOTE when adding a new action, also update utils.existing_output_action below
    ('version'          , {'func' : int, 'help' : 'print version information and exit'}),  # int doesn't do anything, it's just because I have to put something here
    ('cache-parameters' , {'func' : run_partitiondriver, 'help' : 'Cache parameter values and write hmm model files.'}),
    ('annotate'         , {'func' : run_partitiondriver, 'help' : 'Annotate sequences in input file, i.e. run the viterbi algorithm, using pre-existing parameter directory.'}),
    ('partition'        , {'func' : run_partitiondriver, 'help' : 'Partition sequences in input file into clonally-related families using pre-existing parameter directory.'}),
    ('subset-partition' , {'func' : run_partitiondriver, 'help' : 'Split input into --n-subsets subsets, cache parameters and partition each subset individually, then partition the merged results. When combined e.g. with --small-clusters-to-ignore, this will speed up (and reduce memory usage of) partitioning on large samples.'}),
    ('subset-annotate'  , {'func' : run_partitiondriver, 'help' : 'Same as \'subset-partition\', but annotate instead of partition.'}),
    ('simulate'         , {'func' : run_simulation,      'help' : 'Generate simulated sequences based on information in pre-existing parameter directory.'}),
    ('view-output'      , {'func' : run_partitiondriver, 'help' : 'Print partitions and/or annotations from an existing output yaml file.'}),
    ('merge-paired-partitions', {'func' : run_partitiondriver, 'help' : 'Merge the heavy and light chain partition paths from two existing output files.'}),
    ('view-alternative-annotations' , {'func' : run_partitiondriver, 'help' : 'Print (to std out) a comparison of the naive sequences and v/d/j gene calls corresponding to sub- and super-clusters of the cluster specified with --queries. You must have specified --calculate-alternative-annotations in a previous partition step so that this information was saved.'}),
    ('plot-partitions'  , {'func' : run_partitiondriver, 'help' : 'Plot existing partitions and cluster annotations.'}),
    ('get-selection-metrics' , {'func' : run_partitiondriver, 'help' : 'Calculate selection metrics using existing output.'}),
    ('infer-trees' , {'func' : run_partitiondriver, 'help' : 'Infer phylogenetic trees on annotations in existing output file.'}),
    ('get-linearham-info', {'func' : run_partitiondriver, 'help' : 'Write input file for linearham (to --linearham-info-fname), using a previous partis output (--outfname) file as input.'}),
    ('update-meta-info',   {'func' : run_partitiondriver, 'help' : 'Read existing output files and update their info from --input-metafnames, i.e. if your input meta info changed but you don\'t want to rerun everything. Also updates --queries-to-include info.'}),
    ('write-fake-paired-annotations',   {'func' : run_partitiondriver, 'help' : 'Read paired output from --paired-outdir, and write \'fake\' annotations by concatenating h and l seqs to file in the same dir (useful e.g. for phylo method input).'}),
    # deprecated actions:
    ('view-annotations' , {'func' : run_partitiondriver, 'help' : 'Mostly deprecated: Print annotations from an existing old-style annotation output csv (for current yaml output files, use \'view-output\').'}),
    ('view-partitions'  , {'func' : run_partitiondriver, 'help' : 'Mostly deprecated: Print partitions from an existing old-style partition output csv (for current yaml output files, use \'view-output\').'}),
    ('view-alternative-naive-seqs'  , {'func' : run_partitiondriver, 'help' : 'DEPRECATED use \'view-alternative-annotations\''}),
    ('run-viterbi'      , {'func' : run_partitiondriver, 'help' : 'DEPRECATED use \'annotate\''}),
))

utils.existing_output_actions = [
    'view-output',
    'merge-paired-partitions',
    'view-alternative-annotations',
    'plot-partitions',
    'get-selection-metrics',
    'infer-trees',
    'get-linearham-info',
    'update-meta-info',
    'write-fake-paired-annotations',
    'view-annotations',
    'view-partitions',
    'view-alternative-naive-seqs',
]

# ----------------------------------------------------------------------------------------
def runs_on_existing_output(actionstr):
    return actionstr in utils.existing_output_actions

# ----------------------------------------------------------------------------------------
subargs = {subname : [] for subname in subconfig}

# ----------------------------------------------------------------------------------------
subargs['annotate'].append({'name' : '--annotation-clustering', 'kwargs' : {'action' : 'store_true', 'help' : 'Perform annotation-based clustering: group together sequences with the same V and J, same CDR3 length, and 90%% cdr identity. Very, very inaccurate.'}})
subargs['annotate'].append({'name' : '--annotation-clustering-threshold', 'kwargs' : {'default' : 0.9, 'type' : float, 'help' : 'threshold used by --annotation-clustering (i.e. a cluster is defined as the same v and j gene plus this threshold for single-linkage clustering on cdr3 nucleotide hamming distance)'}})
subargs['annotate'].append({'name' : '--paired-naive-hfrac-threshold-type', 'kwargs' : {'default' : 'naive-hamming', 'choices' : ['naive-hamming', 'likelihood'], 'help' : 'see help for this option under \'partition\' action'}})

# ----------------------------------------------------------------------------------------
subargs['partition'].append({'name' : '--naive-vsearch', 'kwargs' : {'action' : 'store_true', 'help' : 'Very fast clustering: infer naive (unmutated ancestor) for each input sequence, then cluster using vsearch (in --cluster_fast mode, see https://doi.org/10.7717/peerj.2584). This is less accurate than default clustering, since it knows nothing about vdj rearrangement, but since it clusters on naive sequences it\'s still much more accurate than any clustering method on mutated sequences.'}})
subargs['partition'].append({'name' : '--fast', 'kwargs' : {'action' : 'store_true', 'help' : 'Synonym for --naive-vsearch (see that help msg)'}})
subargs['partition'].append({'name' : '--no-naive-vsearch', 'kwargs' : {'action' : 'store_true', 'help' : 'Use this to turn off --naive-vsearch (--fast) in cases where it\'s been turned on automatically.'}})
subargs['partition'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'Throw out all sequences that are not clonally related to this sequence id. Much much much faster than partitioning the entire sample (well, unless your whole sample is one family). If --paired-loci is set, this must be a colon-separated list of length two with the heavy:light sequence ids, and you must also set --seed-loci.'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['partition'].append({'name' : '--seed-seq', 'kwargs' : {'help' : 'DEPRECATED use --seed-unique-id and --queries-to-include-fname'}})
subargs['partition'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'Only used when both --seed-unique-id and --paired-loci are set: colon-separated list of length two with the heavy and light chain loci, e.g. \'igh:igk\''}})
subargs['partition'].append({'name' : '--random-seed-seq', 'kwargs' : {'action' : 'store_true', 'help' : 'choose a sequence at random from the input file, and use it as the seed for seed partitioning (as if it had been set as the --seed-unique-id)'}})
subargs['partition'].append({'name' : '--naive-hamming-bounds', 'kwargs' : {'help' : 'Clustering bounds (lo:hi colon-separated pair) on naive sequence hamming distance. If not specified, the bounds are set based on the per-dataset mutation levels. For most purposes should be left at the defaults.'}})
subargs['partition'].append({'name' : '--logprob-ratio-threshold', 'kwargs' : {'type' : float, 'default' : 18., 'help' : 'Reaches a min value of <this> minus five for large clusters. Do *not* change this.'}})
subargs['partition'].append({'name' : '--paired-naive-hfrac-threshold-type', 'kwargs' : {'default' : 'naive-hamming', 'choices' : ['naive-hamming', 'likelihood'], 'help' : 'type of naive hamming fraction threshold to use for paired clustering splitting. Roughly: \'naive-hamming\' is the threshold at which seqs have equal probability of being clonal vs non-clonal, whereas, \'likelihood\' is a higher threshold above which almost all sequences are non-clonal (both are calculated per-dataset based on shm frequency).'}})
subargs['partition'].append({'name' : '--synthetic-distance-based-partition', 'kwargs' : {'action' : 'store_true', 'help' : 'Use simulation truth info to create a synthetic distance-based partition (for validation). Clusters true naive sequences with vsearch, i.e. effectively removes the component of uncertainty that comes from naive sequence inference, leaving only the component from clustering itself.'}})
subargs['partition'].append({'name' : '--cache-naive-hfracs', 'kwargs' : {'action' : 'store_true', 'help' : 'In addition to naive sequences and log probabilities, also cache naive hamming fractions between cluster pairs. Only really useful for plotting or testing.'}})
subargs['partition'].append({'name' : '--biggest-naive-seq-cluster-to-calculate', 'kwargs' : {'type' : int, 'default' : 15, 'help' : 'start thinking about subsampling before you calculate anything if cluster is bigger than this'}})
subargs['partition'].append({'name' : '--biggest-logprob-cluster-to-calculate', 'kwargs' : {'type' : int, 'default' : 15, 'help' : 'start thinking about subsampling before you calculate anything if cluster is bigger than this'}})
subargs['partition'].append({'name' : '--n-partitions-to-write', 'kwargs' : {'type' : int, 'default' : 10, 'help' : 'Number of partitions (surrounding the best partition) to write to output file.'}})
# subargs['partition'].append({'name' : '--naive-swarm', 'kwargs' : {'action' : 'store_true', 'help' : 'Use swarm instead of vsearch, which the developer recommends. Didn\'t seem to help much, and needs more work to optimize threshold, so DO NOT USE.'}})
subargs['partition'].append({'name' : '--small-clusters-to-ignore', 'kwargs' : {'help' : 'colon-separated list (or dash-separated inclusive range) of cluster sizes to throw out after several partition steps. E.g. \'1:2\' will, after <--n-steps-at-which-to-ignore-small-clusters> partition steps, throw out all singletons and pairs. Alternatively, \'1-10\' will ignore all clusters with size less than 11.'}})
subargs['partition'].append({'name' : '--n-steps-after-which-to-ignore-small-clusters', 'kwargs' : {'type' : int, 'default' : 3, 'help' : 'number of partition steps after which to throw out small clusters (where "small" is controlled by <--small-clusters-to-ignore>). (They\'re thrown out before this if we get to n_procs one before this).'}})
subargs['partition'].append({'name' : '--calculate-alternative-annotations', 'kwargs' : {'action' : 'store_true', 'help' : 'write to disk all the information necessary to, in a later step (\'view-alternative-annotations\'), print alternative inferred naive sequences (i.e. visualize uncertainty in the inferred naive sequence). This is largely equivalent to setting --write-additional-cluster-annotations to \'sys.maxint:sys.maxint\'.'}})
subargs['partition'].append({'name' : '--max-cluster-size', 'kwargs' : {'type' : int, 'help' : 'stop clustering immediately if any cluster grows larger than this (useful for limiting memory usage, which can become a problem when the final partition contains very large clusters)'}})
subargs['partition'].append({'name' : '--write-additional-cluster-annotations', 'kwargs' : {'help' : 'in addition to writing annotations for each cluster in the best partition, also write annotations for all the clusters in several partitions on either side of the best partition. Specified as a pair of numbers \'m:n\' for m partitions before, and n partitions after, the best partition.'}})
subargs['partition'].append({'name' : '--add-pairwise-clustering-metrics', 'kwargs' : {'action' : 'store_true', 'help' : 'In addition to the default purity/completeness clustering performance metrics, also write pairwise ones (as defined e.g. in mobille paper).'}})
subargs['partition'].append({'name' : '--naive-hamming-cluster', 'kwargs' : {'action' : 'store_true', 'help' : 'agglomerate purely with naive hamming distance, i.e. set the low and high preclustering bounds to the same value. Not recommended at this point: if you want fast use --naive-vsearch. (Default clustering uses naive hamming distance enough that it\'s not much slower than this, anyway).'}})
subargs['partition'].append({'name' : '--max-n-seqs-to-likelihood-cluster', 'kwargs' : {'type' : int, 'default' : 50000, 'help' : 'If the repertoire has more than this many seqs, --naive-vsearch (--fast) is turned on, i.e. naive vsearch is used for clustering without any likelihoods. Turn off with --no-naive-vsearch.'}})
subargs['partition'].append({'name' : '--continue-from-input-partition', 'kwargs' : {'action' : 'store_true', 'help' : 'When --input-partition-fname is set, continue partitoning starting from the input partition (instead of the default of keeping the input partition without any additional clustering). Arguably this should be the default, but I can\'t change it without breaking backward compatibility.'}})

subargs['subset-partition'].append({'name' : '--n-subsets', 'kwargs' : {'type' : int, 'default' : 5, 'help' : 'Number of subsets into which to divide the input for \'subset-partition\''}})
subargs['subset-annotate'].append({'name' : '--n-subsets', 'kwargs' : {'type' : int, 'default' : 5, 'help' : 'Number of subsets into which to divide the input for \'subset-annotate\''}})

# ----------------------------------------------------------------------------------------
# basic simulation:
subargs['simulate'].append({'name' : '--mutation-multiplier', 'kwargs' : {'type' : float, 'help' : 'Multiply observed branch lengths by some factor when simulating, e.g. if in data it was 0.05, but you want closer to ten percent in your simulation, set this to 2'}})
subargs['simulate'].append({'name' : '--n-sim-events', 'kwargs' : {'type' : int, 'default' : 1, 'help' : 'Number of rearrangement events to simulate'}})
subargs['simulate'].append({'name' : '--n-trees', 'kwargs' : {'type' : int, 'help' : 'Number of phylogenetic trees from which to choose during simulation (we pre-generate this many trees before starting a simulation run, then for each rearrangement event choose one at random -- so this should be at least of order the number of simulated events, so your clonal families don\'t all have the same tree).'}})
subargs['simulate'].append({'name' : '--n-leaf-distribution', 'kwargs' : {'default' : None, 'choices' : ['geometric', 'box', 'zipf', 'hist'], 'help' : 'Distribution from which to draw the number of leaves for each tree. Except for \'hist\' (which is a histogram inferred from data) the distribution\'s parameter comes from --n-leaves. If not set, defaults to --default-scratch-n-leaf-distribution (if rearranging from scratch), or \'hist\' (if rearranging with inferred parameters). Note that in the latter case if the hist file doesn\'t exist, it\'ll fall back to one of the others (with warning). If planning to use \'hist\', you\'ll need to run \'partition\' with --count-parameters set; otherwise the only parameters available will have been inferred from the singleton partition (before partitioning), so the cluster size hist will only have entries in the 1-bin (it will also warn about this).'}})
subargs['simulate'].append({'name' : '--n-leaf-hist-fname', 'kwargs' : {'help' : 'file with histogram of cluster sizes (i.e. from a parameter dir) from which to sample simulation family sizes'}})
subargs['simulate'].append({'name' : '--default-scratch-n-leaf-distribution', 'kwargs' : {'default' : 'geometric', 'choices' : ['geometric', 'box', 'zipf'], 'help' : 'Default for --n-leaf-distribution. Do not set, this is for internal use only.'}})
subargs['simulate'].append({'name' : '--n-leaves', 'kwargs' : {'type' : float, 'default' : 5., 'help' : 'Parameter determining the shape of the distribution from which the number of leaves per tree is drawn (--n-leaf-distribution). For \'geometric\' and \'box\', this is the mean; while for \'zipf\' it is the exponent, i.e. *not* the mean -- see docs for numpy.random.zipf for details. Has *no* effect if --n-leaf-distribution is \'hist\'.'}})
subargs['simulate'].append({'name' : '--constant-number-of-leaves', 'kwargs' : {'action' : 'store_true', 'help' : 'Give all trees the same number of leaves (i.e. override --n-leaf-distribution).'}})
subargs['simulate'].append({'name' : '--allowed-cdr3-lengths', 'kwargs' : {'help' : 'For single chain simulation, a colon-separated list of nucleotide cdr3 lengths to which to restrict the simulation (ranges also allowed, e.g. 36-40 [python slice conventions]). NOTE that our cdr3 definition includes both conserved codons, which differs from the imgt definition (sorry). If --paired-loci is set, specify different lengths for each loci as key/val pairs e.g. \'igh,39c42-46:igk,30c33c36\' (each \'c\' will be replaced by a colon when passing to each locus\'s simulation process). If a locus is not present in the dict, no restrictions will be placed on its lengths. If a locus is present more than once, the lengths for each additional entry are appended. Typical human lengths: igh (36-72), igk (30-39), igl (30-42).'}})
subargs['simulate'].append({'name' : '--remove-nonfunctional-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'Remove non-functional sequences from simulated rearrangement events. Note that because this happens after generating SHM (since we have no way to tell bppseqgen to only generate functional sequences), you will in general need to specify a (potentialy much) larger value for --n-leaves if you set --remove-nonfunctional-seqs. Typically, the vast majority of nonfunctional simulated sequences are due to stop codons generated by SHM. See also --mutate-stop-codons.'}})
subargs['simulate'].append({'name' : '--gtrfname', 'kwargs' : {'default' : partis_dir + '/data/recombinator/gtr.txt', 'help' : 'File with list of GTR parameters. Fed into bppseqgen along with the chosen tree. Corresponds to an arbitrary dataset at the moment, but eventually will be inferred per-dataset.'}})  # NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt)
subargs['simulate'].append({'name' : '--root-mrca-weibull-parameter', 'kwargs' : {'type' : float, 'help' : 'if set, uses TreeSimGM (instead of TreeSim), and passes this value as the weibull parameter (e.g. 0.1: long root-mrca distance, lots of shared mutation; 5: short, little) NOTE requires installation of TreeSimGM'}})
subargs['simulate'].append({'name' : '--input-simulation-treefname', 'kwargs' : {'help' : 'file with newick-formatted lines corresponding to trees to use for simulation. Note that a) the tree depths are rescaled according to the shm rates requested by other command line arguments, i.e. the depths in the tree file are ignored, and b) the resulting sequences do not use the leaf names from the trees (unrelated to --treefname).'}})
subargs['simulate'].append({'name' : '--generate-trees', 'kwargs' : {'action' : 'store_true', 'help' : 'Run the initial tree-generation step of simulation (writing to --outfname), without then proceeding to actually generate the sequences. Used for paired heavy/light simulation so we can pass the same list of trees to both.'}})
subargs['simulate'].append({'name' : '--choose-trees-in-order', 'kwargs' : {'action' : 'store_true', 'help' : 'Instead of the default of choosing a tree at random from the list of trees (which is either generated at the start of the simulation run, or passed in with --input-simulation-treefname), instead choose trees sequentially. If there\'s more events than trees, it cycles through the list again. Used for paired heavy/light simulation.'}})
subargs['simulate'].append({'name' : '--check-tree-depths', 'kwargs' : {'action' : 'store_true', 'help' : 'check how close the fraction of mutated bases that we get from bppseqgen is to the depth of the trees that we passed in.'}})
subargs['simulate'].append({'name' : '--no-per-base-mutation', 'kwargs' : {'action' : 'store_true', 'help' : 'By default the simulation uses both a different SHM rate for each position in each gene, and also a different rate to each base for each of these positions (e.g. position 37 in IGHV1-2*02 might mutate 3x as much as neighboring bases, and it might mutate to T twice as often as G). This option turns off the per-base part of that, so while each position in each gene still has a different overall rate, the rate to each base is the same (e.g. A->T is the same as A->G). Run time with this option set is about 10 times faster than the default. NOTE also --no-per-base-mfreqs'}})
subargs['simulate'].append({'name' : '--mutate-conserved-codons', 'kwargs' : {'action' : 'store_true', 'help' : 'By default we don\'t let mutations occur in conserved codons (cyst in V and tryp/phen in J), but if this option is set we let them mutate.'}})
subargs['simulate'].append({'name' : '--mutate-stop-codons', 'kwargs' : {'action' : 'store_true', 'help' : 'bpp-seqgen doesn\'t have a way to forbid stop codon mutations, and to keep the trees more correct, by default we leave these in the final sequences (and with high shm most sequences will have stop codons). If this arg is set, we instead add additional mutations to each stop codon until it is no longer a stop. This generally means that almost all output sequences will be productive.'}})
subargs['simulate'].append({'name' : '--light-chain-fractions', 'kwargs' : {'default' : 'igk,0.67:igl,0.33', 'help' : 'fraction of events with igk vs igl in paired simulation, in form \'igk,f1:igl,f2\', where f1+f2 must equal to 1. See also --single-light-locus.'}})
subargs['simulate'].append({'name' : '--single-light-locus', 'kwargs' : {'help' : 'Instead of choosing the light chain locus for each rearrangement at random, use this locus for all simulated rearrangements. See also --light-chain-fractions (which this option overrides).'}})
# shm indels:
subargs['simulate'].append({'name' : '--indel-frequency', 'kwargs' : {'default' : 0., 'type' : float, 'help' : 'fraction of simulated sequences with SHM indels'}})
subargs['simulate'].append({'name' : '--n-indels-per-indeld-seq', 'kwargs' : {'default' : '1:2', 'help' : 'list of integers from which to choose the number of SHM indels in each sequence that we\'ve already decided has indels'}})
subargs['simulate'].append({'name' : '--mean-indel-length', 'kwargs' : {'type' : float, 'default' : 5, 'help' : 'mean length of each SHM indel (geometric distribution)'}})
subargs['simulate'].append({'name' : '--indel-location', 'kwargs' : {'help' : 'Where to put simulated SHM indels. Set either to a single integer position at which all indels will occur, or choose from one of three regions (each excluding the first and last five positions in the sequence): None (anywhere in sequence), \'v\' (before cysteine), \'cdr3\' (within cdr3). If set to a single position, any entries in --n-indels-per-indeld-seq greater than 1 are removed, since we don\'t want a bunch of them at the same point.'}})
# from-scratch (rather than mimicking a particular data sample):
subargs['simulate'].append({'name' : '--rearrange-from-scratch', 'kwargs' : {'action' : 'store_true', 'help' : 'Don\'t use an existing parameter directory for rearrangement-level parameters, and instead make up some plausible stuff from scratch. Have to also set --shm-parameter-dir.'}})
subargs['simulate'].append({'name' : '--mutate-from-scratch', 'kwargs' : {'action' : 'store_true', 'help' : 'Don\'t use an existing parameter directory for shm-level (mutation) parameters, and instead make up stuff from scratch (by default this means shm rate varies over positions and sequences, but is the same for all regions). Have to also set --reco-parameter-dir.'}})
subargs['simulate'].append({'name' : '--simulate-from-scratch', 'kwargs' : {'action' : 'store_true', 'help' : 'same as setting both --rearrange-from-scratch and --mutate-from-scratch'}})
subargs['simulate'].append({'name' : '--allow-nonfunctional-scratch-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'By default, scratch rearrangement skips potential rearrangements that are either out of frame or have stop codons. Setting this arg instead allows them, which can be useful to e.g. to speed up simulation if you don\'t care about getting only productive sequences.'}})
subargs['simulate'].append({'name' : '--shm-parameter-dir', 'kwargs' : {'help' : 'parameter directory from which to retrieve shm-level info when --rearrange-from-scratch is set (to set germline info, use --initial-germline-dir).'}})
subargs['simulate'].append({'name' : '--reco-parameter-dir', 'kwargs' : {'help' : 'parameter directory from which to retrieve rearrangement-level info when --mutate-from-scratch is set (to set germline info, use --initial-germline-dir).'}})
subargs['simulate'].append({'name' : '--scratch-mute-freq', 'kwargs' : {'type' : float, 'default' : 0.05, 'help' : 'shm fraction/rate (tree depth passed to bpp) used by --mutate-from-scratch. Note that to get very large average values you need to also set --same-mute-freq-for-all-seqs (and maybe --flat-mute-freq), since saturation means that the high end of the spectrum ends up getting pushed down while the low end doesn\'t.'}})
subargs['simulate'].append({'name' : '--flat-mute-freq', 'kwargs' : {'action' : 'store_true', 'help' : 'use the same shm fraction/rate (--scratch-mute-freq) for all positions (in practice it\'s not that much flatter than the Gamma that is used by default --mutate-from-scratch). For use with --mutate-from-scratch.'}})
subargs['simulate'].append({'name' : '--same-mute-freq-for-all-seqs', 'kwargs' : {'action' : 'store_true', 'help' : 'use the same shm rate (--scratch-mute-freq) for all families and all sequences. This only actually directly determines the mean family shm (so "-for-all-families" might be a better name), but since for the options we need TreeSim can only give ultrametric trees, this is the only lever we have. For use with --mutate-from-scratch. Note that this means we tell bppseqgen to use the same expected rate for every sequence -- there\'s still variance in the resulting number of output mutation per sequence. See note in treegenerator.choose_full_sequence_branch_length().'}})
subargs['simulate'].append({'name' : '--correlation-values', 'kwargs' : {'help' : 'By default, scratch rearrangement has no correlations between parameters, but this arg lets you turn them on. Each correlation is specified as a \'.\'-separated pair of parameter names, separated by \',\' from a float in [0, 1] specifying the magnitude of the correlation, with multiple correlations separated by \':\', e.g. \'v_gene.d_gene,0.1:d_gene.j_gene,0.7\'. Currently allowed parameter pairs are %s. Note that the correlation model is just one (fairly arbitrary) choice for how to set things up, for details consult the code. Note also that non-scratch simulation, i.e. from inferred parameters, includes correlations among all parameters.' % ', '.join(str(p) for p in utils.available_simu_correlations)}})
subargs['simulate'].append({'name' : '--paired-correlation-values', 'kwargs' : {'help' : 'like --correlation-values, except that the first parameter is in the heavy chain rearrangmeent, and the second is in the light chain rearrangement.'}})
subargs['simulate'].append({'name' : '--heavy-chain-event-fname', 'kwargs' : {'help' : 'for internal use by --paired-correlation-values: file with heavy chain events so the light chain rearrangement knows what was in the heavy chain rearrangement.'}})
# new allele/germline set generation (this is also for from-scratch)
subargs['simulate'].append({'name' : '--generate-germline-set', 'kwargs' : {'action' : 'store_true', 'help' : 'Choose a realistic germline set from the available genes. Turned on automatically when --rearrange-from-scratch or --simulate-from-scratch are set.'}})
subargs['simulate'].append({'name' : '--force-dont-generate-germline-set', 'kwargs' : {'action' : 'store_true', 'help' : 'Force --generate-germline-set to False in situations in which it would otherwise be turned on (I know, this is kind of messy, but I think it\'s the best available solution).'}})
subargs['simulate'].append({'name' : '--n-genes-per-region', 'kwargs' : {'default' : '::', 'help' : 'colon-separated list specifying the number of genes (not alleles) for each region with order v:d:j (for use with --generate-germline-set). So the *total* number of alleles is this times the number from --n-sim-alleles-per-gene). Leave a region\'s entry blank to use default: e.g. \'5::\' would choose 5 v genes, while choosing the default numbers of d and j genes. Default: %s' % glutils.default_n_genes_per_region}})
subargs['simulate'].append({'name' : '--n-sim-alleles-per-gene', 'kwargs' : {'default' : '::', 'help' : 'colon-separated list of mean alleles per gene for each region with order v:d:j (for use with --generate-germline-set). Leave a region\'s entry blank to use the default, e.g. \'1.2::\' uses default values for d and j. Default: %s' % glutils.default_n_alleles_per_gene}})
subargs['simulate'].append({'name' : '--min-sim-allele-prevalence-freq', 'kwargs' : {'default' : glutils.default_min_allele_prevalence_freq, 'type' : float, 'help' : 'minimum frequency at which alleles are allowed to occur, e.g. if it\'s 0.01 then each pair of V alleles will have a prevalence ratio between 0.01 and 1'}})
subargs['simulate'].append({'name' : '--allele-prevalence-fname', 'kwargs' : {'help' : 'mostly for internal use (specifies allele prevalence freqs, as written by glutils.generate_germline_set())'}})
subargs['simulate'].append({'name' : '--nsnp-list', 'kwargs' : {'help' : 'When --generate-germline-set is set, this is a colon-separated list whose length gives the number of novel alleles to add to the germline set. Each entry in the list is the number of SNPs to generate for the corresponding new allele. If both --nsnp-list and --nindel-list are set, they must be the same length; if one is unset, it is assumed to be all zeroes.'}})
subargs['simulate'].append({'name' : '--nindel-list', 'kwargs' : {'help' : 'When --generate-germline-set is set, this is a colon-separated list whose length gives the number of novel alleles to add to the germline set. Each entry in the list is the number of indels to generate for the corresponding new allele. If both --nsnp-list and --nindel-list are set, they must be the same length; if one is unset, it is assumed to be all zeroes.'}})
# misc
subargs['simulate'].append({'name' : '--im-a-subproc', 'kwargs' : {'action' : 'store_true', 'help' : 'Internal use only. This is set to true if this is a subprocess, i.e. set when --n-procs > 1.'}})
subargs['simulate'].append({'name' : '--mean-cells-per-droplet', 'kwargs' : {'help' : 'If set, in paired simulation add extra cells per droplet (to mimic ambiguous pairing info, as in "overloaded" 10x data). Cells are apportioned uniform-randomly among the correct number of droplets so as to arrive at the specified mean cells. Correct pairing info is left in original simulation file, while new, incorrect info is written to a new input meta file, which if passed to partis inference with --input-metafnames will override pairing info from simulation file. See also --fraction-of-reads-to-remove. If not set (or set to the string \'None\'), all droplets have exactly one cell (which is *not* the same as setting this arg to 1).'}})
subargs['simulate'].append({'name' : '--constant-cells-per-droplet', 'kwargs' : {'action' : 'store_true', 'help' : 'instead of the default of choosing the droplet for each cell uniform randomly from all droplets, apportion them as evenly as possible.'}})
subargs['simulate'].append({'name' : '--fraction-of-reads-to-remove', 'kwargs' : {'type' : float, 'help' : 'If set, in paired simulation we throw out this fraction of sequences (so that their paired seq, which remains, has no paired seq). See also --mean-cells-per-droplet.'}})
subargs['simulate'].append({'name' : '--bulk-data-fraction', 'kwargs' : {'type' : float, 'help' : 'If set, in paired simulation we remove pairing info for this fraction of sequence pairs (i.e. as if this fraction of sample is bulk data, while the rest is single cell).'}})

subargs['view-alternative-annotations'].append({'name' : '--print-all-annotations', 'kwargs' : {'action' : 'store_true', 'help' : 'In addition to printing each alternative naive seq and gene call, with information on which clusters support its correctness, also print the annotations for all of these clusters below each naive sequence and gene call.'}})
subargs['view-alternative-annotations'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['view-alternative-annotations'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['view-output'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['view-output'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['get-selection-metrics'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['get-selection-metrics'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['get-selection-metrics'].append({'name' : '--chosen-ab-fname', 'kwargs' : {'help' : 'csv file to which to write information on chosen antibody sequences'}})
subargs['infer-trees'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['infer-trees'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['merge-paired-partitions'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['merge-paired-partitions'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['merge-paired-partitions'].append({'name' : '--random-seed-seq', 'kwargs' : {'action' : 'store_true', 'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['merge-paired-partitions'].append({'name' : '--paired-naive-hfrac-threshold-type', 'kwargs' : {'default' : 'naive-hamming', 'choices' : ['naive-hamming', 'likelihood'], 'help' : 'see help for this option under \'partition\' action'}})
subargs['plot-partitions'].append({'name' : '--seed-unique-id', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)
subargs['plot-partitions'].append({'name' : '--seed-loci', 'kwargs' : {'help' : 'see help for this option under \'partition\' action'}})  # NOTE do *not* move these up above -- we forbid people to set them for auto parameter caching (see exception above)

subargs['get-linearham-info'].append({'name' : '--linearham-info-fname', 'kwargs' : {'help' : 'yaml file to which to write linearhmam input information'}})

sub_arg_groups = {'subset-partition' : ['partition'], 'subset-annotate' : ['annotate']}  # actions that use the args of other actions (i.e. their <subargs> need to be merged, e.g. 'subset-partition' needs access to all the args of 'partition')

def get_arg_names(actions):  # return set of all arg names (in the form they appear in args.__dict__) for the specified actions
    if actions == 'all':
        actions = list(subconfig.keys())
    return set([actionconf['name'][2:].replace('-', '_') for action in actions for actionconf in subargs[action]])

# above we just made a dict with lists of args, here we actually make the sub parsers
for argconf in parent_args:
    parent_parser.add_argument(argconf['name'], **argconf['kwargs'])
subparsermap = {}
for name, vals in subconfig.items():
    subparsermap[name] = subparsers.add_parser(name, parents=[parent_parser], help=vals['help'], formatter_class=MultiplyInheritedFormatter)
    subparsermap[name].set_defaults(func=vals['func'])  # set the default fcn to run for <name> (i.e. action)
    if name in sub_arg_groups:
        subargs[name] += [c for n in sub_arg_groups[name] for c in subargs[n]]
    for argconf in subargs[name]:
        subparsermap[name].add_argument(argconf['name'], **argconf['kwargs'])

# ----------------------------------------------------------------------------------------
args = parser.parse_args()

# add OR of all arguments to all subparsers to <args>, as None (to avoid having to rewrite a *##!(%ton of other code)
for missing_arg in get_arg_names('all') - set(args.__dict__):  # NOTE see also remove_action_specific_args() above, which kind of does similar things, at least if i rewrote this all from scratch only one (or neither...) would exist
    args.__dict__[missing_arg] = None

processargs.process(args)
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
start = time.time()
args.func(args)
print('      total time: %.1f' % (time.time()-start))
