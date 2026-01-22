from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import collections
import copy
import itertools
import numpy
import sys
import operator
import os
import math
import json

from . import hutils
from . import utils
from . import prutils
from .clusterpath import ptnprint, ClusterPath
from .hist import Hist
from . import glutils
from . import treeutils
from io import open

# ----------------------------------------------------------------------------------------
def chstr(n_before, n_after):  # empty if same number before and after, otherwise red +/-N
    if n_before == n_after: return ''
    else: return ' ' + utils.color('red', '%+d' % (n_after - n_before))

# ----------------------------------------------------------------------------------------
def subd(seed_unique_id=None, lpair=None, single_chain=None):
    subd = ''
    if seed_unique_id is not None:
        subd += '/seeds/%s' % '+'.join(seed_unique_id)
    if lpair is not None:
        subd += '/%s' % '+'.join(lpair)
    if single_chain:
        subd += '/single-chain'
    return subd

# ----------------------------------------------------------------------------------------
# return standardized file name (including subdirs) in directory structure that we use for paired heavy/light i/o
def paired_fn(bdir, locus, lpair=None, suffix='.fa', ig_or_tr='ig', actstr=None, seed_unique_id=None, single_chain=None):  # if set, only file(s) for this <locus>, and/or only files for this <lpair> of loci. If <lpair> is set but <locus> is None, returns subdir name
    bdir += subd(seed_unique_id=seed_unique_id, lpair=lpair, single_chain=single_chain)  # i think you only ever want to specify one of these
    if lpair is not None:
        if locus is None:
            return bdir
    if actstr == 'auto':  # return the first one that exists
        for astr in [None, 'partition']:  # i think these are the only ones i use atm
            fn = paired_fn(bdir, locus, lpair=lpair, suffix=suffix, ig_or_tr=ig_or_tr, actstr=astr, seed_unique_id=seed_unique_id, single_chain=single_chain)
            if os.path.exists(fn):
                return fn
        raise Exception('couldn\'t find paired fn with auto actstr in %s' % bdir)
    return '%s/%s%s%s' % (bdir, '' if actstr is None else actstr+'-', locus, suffix)

# ----------------------------------------------------------------------------------------
# NOTE this doesn't include [e.g.] the single-chain files, since it was originally just written for simulation + split-loci.py files
def paired_dir_fnames(bdir, no_pairing_info=False, only_paired=False, suffix='.fa', ig_or_tr='ig', include_failed=False, include_meta=False, no_dirs=False):  # return all files + dirs from previous fcn
    fnames = []
    if not only_paired:
        fnames += [paired_fn(bdir, l, suffix=suffix) for l in utils.sub_loci(ig_or_tr)]  # single-locus files
    if include_failed:
        fnames += [paired_fn(bdir, 'failed', suffix=suffix)]  # kind of hackey, but bin/split-loci.py does this, so we need to be able to clean it up
    if include_meta:
        fnames += [paired_fn(bdir, 'meta', suffix='.yaml')]  # this is also kind of hackey
    if not no_pairing_info:
        fnames += [paired_fn(bdir, l, lpair=lp, suffix=suffix) for lp in utils.locus_pairs[ig_or_tr] for l in lp]  # paired single-locus files
        if not no_dirs:
            fnames += [paired_fn(bdir, None, lpair=lp, suffix=suffix) for lp in utils.locus_pairs[ig_or_tr]]  # paired subdirs
    return fnames  # return dirs after files for easier cleaning (below)

# ----------------------------------------------------------------------------------------
def prep_paired_dir(bdir, clean=False, suffix='.fa', extra_files=None, ig_or_tr='ig'):
    if clean and os.path.exists(bdir):
        clean_paired_dir(bdir, suffix=suffix, extra_files=extra_files, expect_missing=True)
    for lpair in utils.locus_pairs[ig_or_tr]:
        utils.prep_dir(paired_fn(bdir, None, lpair=lpair))

# ----------------------------------------------------------------------------------------
def clean_paired_dir(bdir, suffix='.fa', extra_files=None, expect_missing=False, ig_or_tr='ig'):
    print('    cleaning paired files in %s/' % bdir)
    fnames = paired_dir_fnames(bdir, suffix=suffix, ig_or_tr=ig_or_tr)
    if extra_files is not None:
        fnames = extra_files + fnames  # put 'em at the start, since presumably they're actual files, not dirs
    utils.clean_files(fnames, expect_missing=expect_missing)

# ----------------------------------------------------------------------------------------
def check_for_cross_locus_duplicates(uid_lists_by_locus, extra_str=''):
    """Check that no UID appears in multiple loci. Raises exception if any duplicates found.

    Args:
        uid_lists_by_locus: dict mapping locus name to list/set of UIDs for that locus
        extra_str: optional string to prepend to error message
    """
    uid_to_loci = {}
    for locus, uids in uid_lists_by_locus.items():
        for uid in uids:
            if uid not in uid_to_loci:
                uid_to_loci[uid] = []
            uid_to_loci[uid].append(locus)

    cross_locus_dups = {u : loci for u, loci in uid_to_loci.items() if len(loci) > 1}
    if len(cross_locus_dups) > 0:
        err_str = '%sfound %d uid%s in multiple loci (which will cause hard to track crashes):\n' % (extra_str, len(cross_locus_dups), utils.plural(len(cross_locus_dups)))
        for uid, loci in sorted(cross_locus_dups.items()):
            err_str += '  %s: in %s\n' % (uid, ', '.join(loci))
        raise Exception(err_str)

# ----------------------------------------------------------------------------------------
# add_selection_metrics: list of selection metrics to add (plotdir_fcn is atm only if adding selection metrics)
def read_locus_output_files(tmploci, ofn_fcn, lpair=None, read_selection_metrics=False, add_selection_metrics=None, lb_tau=None, plotdir_fcn=None, seed_unique_id=None, dont_add_implicit_info=False, dbgstr='', debug=False):
    # ----------------------------------------------------------------------------------------
    def parse_pairing_info(ltmp, atnlist):
        if seed_unique_id is not None:  # I'm really not sure this is the best place to do this, but we have to add the seed id pairing info at some point
            assert isinstance(seed_unique_id, list) and len(seed_unique_id) == 2  # just a reminder
            sid = seedid(ltmp)  # maybe this would be too slow to run for every line?
            other_sid = utils.get_single_entry([u for u in seed_unique_id if u != sid])
            for sline in [l for l in atnlist if sid in l['unique_ids']]:
                if 'paired-uids' in sline:  # if the seed id is in a cluster with non-seed seqs, the key'll be there but it'll be None
                    pids = utils.per_seq_val(sline, 'paired-uids', sid)
                    if pids is not None and pids != [other_sid]:
                        print('  %s unexpected seq(s) paired with seed id %s: %s' % (utils.color('yellow', 'warning'), sid, pids))
                else:  # but if it's by itself the key won't be there
                    sline['paired-uids'] = [[] for _ in sline['unique_ids']]
                sline['paired-uids'][sline['unique_ids'].index(sid)] = [other_sid]
        for tline in atnlist:  # unfortunately we *also* need to replace any None values with [], since somehow when writing/reading sw cache info missing ones end up as Nones. Should really fix that rather than doing this, but oh well
            if 'paired-uids' in tline:
                inones = [i for i, pids in enumerate(tline['paired-uids']) if pids is None]
                for itmp in inones:
                    tline['paired-uids'][itmp] = utils.input_metafile_defaults('paired-uids')
            else:
                tline['paired-uids'] = [utils.input_metafile_defaults('paired-uids') for _ in tline['unique_ids']]
    # ----------------------------------------------------------------------------------------
    def read_smetrics(ofn, atnlist):
        with open(treeutils.smetric_fname(ofn)) as sfile:
            sminfos = json.load(sfile)
        sminfos = {':'.join(smfo['unique_ids']) : smfo for smfo in sminfos}  # convert from list to dict
        for atn in atnlist:
            if ':'.join(atn['unique_ids']) in sminfos:  # all the ones that were too small won't be there, for instance
                atn['tree-info'] = sminfos[':'.join(atn['unique_ids'])]
    # ----------------------------------------------------------------------------------------
    nulldict = {'glfos' : None, 'antn_lists' : [], 'cpaths' : ClusterPath(partition=[])}
    lpfos = {k : {} for k in ['glfos', 'antn_lists', 'cpaths']}
    for ltmp in tmploci:  # read single-locus output files
        ofn = ofn_fcn(ltmp, lpair=lpair)
        if not os.path.exists(ofn):
            for k in lpfos:
                lpfos[k][ltmp] = copy.deepcopy(nulldict[k])
            if debug:
                print('%s: no %s %s output file, skipping: %s' % (utils.color('blue', '+'.join(lpair) if lpair is not None else ltmp), ltmp, dbgstr, ofn))
            continue
        lpfos['glfos'][ltmp], lpfos['antn_lists'][ltmp], lpfos['cpaths'][ltmp] = utils.read_output(ofn, dont_add_implicit_info=dont_add_implicit_info, skip_failed_queries=True)
        if debug:
            print('    read %d %s annotations with %d seqs from %s' % (len(lpfos['antn_lists'][ltmp]), ltmp, sum(len(l['unique_ids']) for l in lpfos['antn_lists'][ltmp]), ofn))
        if read_selection_metrics and os.path.exists(treeutils.smetric_fname(ofn)):  # if it doesn't exist, the info should be in the regular output file
            read_smetrics(ofn, lpfos['antn_lists'][ltmp])
        if add_selection_metrics is not None:
            print('            adding %d selection metric%s to annotations read from %s%s: %s' % (len(add_selection_metrics), utils.plural(len(add_selection_metrics)), ofn_fcn(ltmp, lpair=lpair), '' if plotdir_fcn is None else ' and plotting to %s' % plotdir_fcn(ltmp, lpair=lpair), ' '.join(add_selection_metrics)))
            for smetric in add_selection_metrics:
                treeutils.calculate_individual_tree_metrics(smetric, lpfos['antn_lists'][ltmp], base_plotdir=None if plotdir_fcn is None else plotdir_fcn(ltmp, lpair=lpair), lb_tau=lb_tau)
        parse_pairing_info(ltmp, lpfos['antn_lists'][ltmp])

    # Check for UIDs that appear in multiple loci (NOTE this is also checked earlier in split-loci.py, but we keep it here as a safety check)
    uid_lists = {ltmp : [uid for atn in lpfos['antn_lists'][ltmp] for uid in atn['unique_ids']]
                 for ltmp in tmploci if lpfos['antn_lists'][ltmp] is not None}
    check_for_cross_locus_duplicates(uid_lists)

    if all(lpfos['glfos'][l] is None for l in tmploci):  # if there was no info for *any* of the loci, set Nones one level up (it's just easier to have the Nones there)
        lpfos = {k : None for k in lpfos}
    return lpfos

# ----------------------------------------------------------------------------------------
def read_lpair_output_files(lpairs, ofn_fcn, read_selection_metrics=False, add_selection_metrics=None, lb_tau=None, plotdir_fcn=None, seed_unique_id=None, dont_add_implicit_info=False, dbgstr='', debug=False):
    lp_infos = {}
    for lpair in lpairs:
        lp_infos[tuple(lpair)] = read_locus_output_files(lpair, ofn_fcn, lpair=lpair, read_selection_metrics=read_selection_metrics, add_selection_metrics=add_selection_metrics, lb_tau=lb_tau, plotdir_fcn=plotdir_fcn,
                                                         seed_unique_id=seed_unique_id, dont_add_implicit_info=dont_add_implicit_info, dbgstr=dbgstr, debug=debug)
    return lp_infos

# ----------------------------------------------------------------------------------------
# similar to read_lpair_output_files(), except you don't need to pass in the file name fcn (adding this late, so could probably use it in a lot more places)
# if <joint>, reads joint/merged h, k, l files (rather than the files in each paired (igh+ig[kl]) subdir (NOTE output/return structure is different in the two cases)
def read_paired_dir(pdir, joint=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def getofn(ltmp, lpair=None):
        ofn = paired_fn(pdir, ltmp, lpair=lpair, suffix='.yaml')
        if not os.path.exists(ofn):  # first look for simu file (e.g. igh.yaml), if it's not there look for the partition output file
            ofn = paired_fn(pdir, ltmp, lpair=lpair, actstr='partition', suffix='.yaml')
        return ofn
    # ----------------------------------------------------------------------------------------
    if joint:  # NOTE different return dict structure for joint vs not
        return read_locus_output_files(utils.sub_loci('ig'), getofn, debug=debug)
    else:
        return read_lpair_output_files(utils.locus_pairs['ig'], getofn, debug=debug)  # NOTE both these fcns ignore a bunch of kw args, which I should probably use at some point

# ----------------------------------------------------------------------------------------
def write_lpair_output_files(lpairs, lp_infos, ofn_fcn, headers, use_pyyaml=False, dont_write_git_info=False):
    def glpf(p, k, l):  # NOTE duplicates code in concat_heavy_chain()
        if lp_infos[tuple(p)][k] is None:
            return None
        return lp_infos[tuple(p)][k].get(l)
    for lpair in lpairs:
        for ltmp in lpair:
            if glpf(lpair, 'glfos', ltmp) is None:
                continue
            utils.write_annotations(ofn_fcn(ltmp, lpair=lpair), glpf(lpair, 'glfos', ltmp), glpf(lpair, 'antn_lists', ltmp), headers, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info)

# ----------------------------------------------------------------------------------------
# write concatd heavy chain <glfos> and <antn_lists>, and use <ofn_fcn> to link to existing light chain files in paired subdirs
def write_concatd_output_files(glfos, antn_lists, ofn_fcn, headers, use_pyyaml=False, work_fnames=None, cpaths=None, true_outfos=None, dont_write_git_info=False, airr_output=False, args=None, write_light_chain_files=False, fail_frac=None):
    for ltmp in sorted(glfos):  # don't really need to write igh first, but i guess it's nice to be consistent
        ofn = ofn_fcn(ltmp, joint=True)
        if utils.has_d_gene(ltmp) or write_light_chain_files:
            if not utils.has_d_gene(ltmp) and write_light_chain_files:
                print('       writing (rather than linking) light chain files: %s' % ofn)
            cp = ClusterPath(partition=utils.get_partition_from_annotation_list(antn_lists[ltmp])) if cpaths is None else cpaths[ltmp]
            partition_lines = cp.get_partition_lines(true_partition=None if true_outfos is None else true_outfos['merged']['cpaths'][ltmp].best(), calc_missing_values='best', fail_frac=fail_frac)
            utils.write_annotations(ofn, glfos[ltmp], antn_lists[ltmp], headers, partition_lines=partition_lines, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info)
            if airr_output:
                utils.write_airr_output(utils.replace_suffix(ofn, '.tsv'), antn_lists[ltmp], cpath=cp, args=args, glfo=glfos[ltmp])
        else:
            utils.makelink(os.path.dirname(ofn), ofn_fcn(ltmp, lpair=utils.getlpair(ltmp)), ofn)
            if airr_output:
                utils.makelink(os.path.dirname(ofn), utils.replace_suffix(ofn_fcn(ltmp, lpair=utils.getlpair(ltmp)), '.tsv'), utils.replace_suffix(ofn, '.tsv'))
        if work_fnames is not None:
            work_fnames.append(ofn)

# ----------------------------------------------------------------------------------------
# merges info (keys 'glfos', 'antn_lists', and 'cpaths') for ltmp (of lpair) from <lp_infos> into <glfos>, <antn_lists>, and <joint_cpaths>
# NOTE the deep copies can be really important, since we later deduplicate the concatd partitions + annotations (they didn't used to be there, which caused a bug/crash) [well even if we didn't deduplicate, it was fucking stupid to not deep copy them)
# NOTE also that you'll probably have duplicates in both the partitions and annotations after this
def merge_locus_lpfo(glfos, antn_lists, joint_cpaths, lpair, ltmp, lp_infos, dont_deep_copy=False, dbgstr='', debug=False):
    # ----------------------------------------------------------------------------------------
    def glpf(p, k, l):  # short for "get key (k) from lp_infos for lpair (p) and locus (l) NOTE duplicates code in write_lpair_output_files
        if p is None:  # single-locus (i.e. not paired)
            if lp_infos[k] is None:
                return None
            return lp_infos[k].get(l)
        if tuple(p) not in lp_infos or lp_infos[tuple(p)][k] is None:
            return None
        return lp_infos[tuple(p)][k].get(l)
    # ----------------------------------------------------------------------------------------
    def dfn(val):
        return val if dont_deep_copy else copy.deepcopy(val)
    # ----------------------------------------------------------------------------------------
    def dbgprint(ststr, add_tot=False):
        if debug:
            tlist = glpf(lpair, 'antn_lists', ltmp)
            def nseqs(tl): return sum(len(l['unique_ids']) for l in tl)
            totstr = ''
            if add_tot:
                totstr = ' for total %d (%d)' % (len(antn_lists[ltmp]), nseqs(antn_lists[ltmp]))
            print('      %s%s %s %d annotations (%d seqs)%s' % (dbgstr, utils.locstr(ltmp), ststr, len(tlist), nseqs(tlist), totstr))
    # ----------------------------------------------------------------------------------------
    if glpf(lpair, 'glfos', ltmp) is None:  # this lpair's output files were empty
        if debug:
            print('      %s%s empty files' % (dbgstr, utils.locstr(ltmp)))
        return
    if ltmp in glfos:  # originally: heavy chain second time through: merge chain glfos from those paired with igk and with igl (now more general though)
        dbgprint('adding', add_tot=True)
        glfos[ltmp] = glutils.get_merged_glfo(glfos[ltmp], glpf(lpair, 'glfos', ltmp))
        antn_lists[ltmp] += dfn(glpf(lpair, 'antn_lists', ltmp))
        if glpf(lpair, 'cpaths', ltmp) is not None:
            if len(joint_cpaths[ltmp].partitions) != 1:  # they should always be 1 anyway, and if they weren't, it'd make it more complicated to concatenate them
                print('        %s multiple partitions in cpath when merging loci' % utils.wrnstr())  # ugh, maybe it isn't a big deal? I'm getting longer ones now, probably from the keep all unpaired/subsetpartition stuff
            joint_cpaths[ltmp] = ClusterPath(partition=joint_cpaths[ltmp].best() + dfn(glpf(lpair, 'cpaths', ltmp).best()))
    else:  # light + heavy chain first time through
        dbgprint('starting with')
        glfos[ltmp] = dfn(glpf(lpair, 'glfos', ltmp))
        antn_lists[ltmp] = dfn(glpf(lpair, 'antn_lists', ltmp))
        if glpf(lpair, 'cpaths', ltmp) is not None:
            joint_cpaths[ltmp] = dfn(glpf(lpair, 'cpaths', ltmp))

# ----------------------------------------------------------------------------------------
def lp_merge_final_dbg(tloci, antn_lists):
    print('    finished with:%s' % (' zero annotations' if len(antn_lists)==0 else ''))
    for ltmp in [l for l in tloci if l in antn_lists]:
        print('     %s: %d annotations (%d seqs)' % (utils.locstr(ltmp), len(antn_lists[ltmp]), sum(len(l['unique_ids']) for l in antn_lists[ltmp])))

# ----------------------------------------------------------------------------------------
# similar to concat_heavy_chain(), except this merges loci only within each lpair (for 'subset-partition' action), whereas that merges heavy chains that're paired with either igk or igl
def concat_lpair_chains(lpairs, lpfo_list, dont_deep_copy=False, debug=False):
    # ----------------------------------------------------------------------------------------
    if debug:
        print('  concatenating locus pair loci among %d lp_infos' % len(lpfo_list))
    final_lp_infos = {}
    for lpair in lpairs:
        if debug:
            print('   %s' % utils.lpstr(lpair))
        glfos, antn_lists, joint_cpaths = {}, {}, {}
        for lp_infos in lpfo_list:
            for ltmp in lpair:
                merge_locus_lpfo(glfos, antn_lists, joint_cpaths, lpair, ltmp, lp_infos, dont_deep_copy=dont_deep_copy, debug=debug)
        final_lp_infos[tuple(lpair)] = {'glfos' : glfos, 'antn_lists' : antn_lists, 'cpaths' : joint_cpaths}
        if debug:
            lp_merge_final_dbg(lpair, antn_lists)
    return final_lp_infos

# ----------------------------------------------------------------------------------------
# like concat_lpair_chains(), except doesn't subdivide by locus pairs (but unlike concat_heavy_chain(), doesn't just merge the heavy chain)
def concat_locus_chains(tmploci, lpfo_list, dont_deep_copy=False, dbgstr='', debug=False):  # NTOE <lpfo_list> has different structure to that in concat_lpair_chains()
    if debug:
        print('  concatenating single-locus info among %d dicts%s' % (len(lpfo_list), dbgstr))
    glfos, antn_lists, cpaths = {}, {}, {}
    for ilp, lpfo in enumerate(lpfo_list):
        for iloc, ltmp in enumerate(tmploci):
            merge_locus_lpfo(glfos, antn_lists, cpaths, None, ltmp, lpfo, dont_deep_copy=dont_deep_copy, dbgstr=('%s '%ilp) if iloc==0 else '  ', debug=debug)
    if debug:
        lp_merge_final_dbg(tmploci, antn_lists)
    return {'glfos' : glfos, 'antn_lists' : antn_lists, 'cpaths' : cpaths}  # I *think* I don't need the deduplication stuff from heavy chain concat fcn below, although I think that's only because atm I'm using this fcn for totally disjoint subsets from subset-partition (so, ugh, should add it)

# ----------------------------------------------------------------------------------------
# similar to concat_lpair_chains() above (see also handle_concatd_heavy_chain() below, which runs this and then does some more processing)
def concat_heavy_chain(lpairs, lp_infos, dont_deep_copy=False, debug=False):  # yeah yeah this name sucks but i want it to be different to the one in the calling scripts
    debug = True
    if debug:
        print('  concatenating heavy chain loci:')
    glfos, antn_lists, joint_cpaths = {}, {}, {}
    for lpair in lpairs:
        for ltmp in lpair:
            merge_locus_lpfo(glfos, antn_lists, joint_cpaths, lpair, ltmp, lp_infos, dont_deep_copy=dont_deep_copy, debug=debug)
    return glfos, antn_lists, joint_cpaths

# ----------------------------------------------------------------------------------------
# concatenate partitions + annotations for the heavy chains paired with k and l
def handle_concatd_heavy_chain(lpairs, lp_infos, dont_deep_copy=False, ig_or_tr='ig', dont_calculate_annotations=False, seed_unique_id=None, dont_deduplicate=False, debug=False):  # name sucks, but can't think of anything better
    glfos, antn_lists, cpaths = concat_heavy_chain(lpairs, lp_infos, dont_deep_copy=dont_deep_copy)  # NOTE this is a pretty arbitrary way to combine the partitions for the seqs with uncertain pairing info, but whatever
    lpfos = {'glfos' : glfos, 'antn_lists' : antn_lists, 'cpaths' : cpaths}
    hloc = utils.heavy_locus(ig_or_tr)
    if len(lpfos['glfos']) == 0:
        return lpfos
    if not dont_deduplicate:
        tmp_ptn = utils.get_deduplicated_partitions([lpfos['cpaths'][hloc].best()], antn_list=lpfos['antn_lists'][hloc] if not dont_calculate_annotations else None, glfo=lpfos['glfos'][hloc], debug=debug)[0]  # have to remove duplicates from heavy partitions and annotations (since seqs that we don't have good pairing info for get put in both light chain dirs, so appear twice in concat'd heavy chain)
        lpfos['cpaths'][hloc] = ClusterPath(partition=tmp_ptn, seed_unique_id=seed_unique_id)
    lp_merge_final_dbg([hloc], antn_lists)
    return lpfos

# ----------------------------------------------------------------------------------------
# somewhat similar to get_antn_pairs() and find_cluster_pairs() below, but operates on single sequences
def find_seq_pairs(antn_lists, ig_or_tr='ig'):
    # ----------------------------------------------------------------------------------------
    def handle_atntn(ltmp, antn):
        for iseq, (tid, pids, iseq) in enumerate(zip(antn['unique_ids'], antn['paired-uids'], antn['input_seqs'])):
            tstr, ostr = ('h', 'l') if utils.has_d_gene(ltmp) else ('l', 'h')
            if len(pids) == 0:  # add all the unpaired seqs
                ofo = {'%s_id'%tstr : tid, '%s_locus'%tstr : ltmp, '%s_seq'%tstr : iseq, '%s_id'%ostr : '', '%s_locus'%ostr : '', '%s_seq'%ostr : ''}
                ofo['l_antn'] = None
            elif len(pids) == 1 and tstr == 'h':  # write h/l pairs when <ltmp> is the h locus
                if pids[0] in all_seqs:
                    lfo = all_seqs[pids[0]]
                    ofo = {'h_id' : tid, 'h_locus' : ltmp, 'h_seq' : iseq, 'l_id' : pids[0], 'l_locus' : lfo['locus'], 'l_seq' : lfo['seq']}
                    ofo['l_antn'] = all_antns[ofo['l_id']]
                else:
                    print('  missing %s' % pids[0])
                    continue
            else:
                continue
            ofo['h_antn'] = antn
            outfos.append(ofo)
    # ----------------------------------------------------------------------------------------
    all_seqs = {u : {'seq' : s, 'locus' : l} for l, alist in antn_lists.items() for a in alist for u, s in zip(a['unique_ids'], a['input_seqs'])}
    all_antns = {u : l for llist in antn_lists.values() for l in llist for u in l['unique_ids']}
    outfos = []
    for ltmp in sorted(antn_lists):
        for antn in antn_lists[ltmp]:
            handle_atntn(ltmp, antn)
    return outfos

# ----------------------------------------------------------------------------------------
def get_all_antn_pairs(lp_infos, size_sort=False, ig_or_tr='ig'):  # same as previous fcn, but lumps together both h+k and h+l, rather than just giving you one or the other
    antn_pairs = []
    for lpair in [lpk for lpk in utils.locus_pairs[ig_or_tr] if tuple(lpk) in lp_infos]:
        antn_pairs += get_antn_pairs(lpair, lp_infos[tuple(lpair)])
    if size_sort:
        antn_pairs = sorted(antn_pairs, key=lambda x: sum(len(l['unique_ids']) for l in x), reverse=True)  # sort by sum of number of heavy and light uids
    return antn_pairs

# ----------------------------------------------------------------------------------------
# similar to find_cluster_pairs(), but this just assumes the annotations are in order (i.e. probably only really safe for simulation)
def get_antn_pairs(lpair, lpfos):  # return list of (hline, lline) pairs
    if None in list(lpfos.values()):
        return []
    if len(set(len(lpfos['antn_lists'][l]) for l in lpair)) != 1:  # if the lists for both loci aren't the same length, you should probably use find_cluster_pairs() (although in reality length differences are probably from unpaired singletons getting added, so it's probably ok)
        print('  %s different length annotation lists (probably just from unpaired singletons) among %s: %s' % (utils.wrnstr(), lpair, [len(lpfos['antn_lists'][l]) for l in lpair]))
    return list(zip(*[lpfos['antn_lists'][l] for l in lpair]))

# ----------------------------------------------------------------------------------------
def find_all_cluster_pairs(lp_infos, required_keys=None, quiet=False, min_cluster_size=None, min_cluster_arg_str='', ig_or_tr='ig', debug=False):
    antn_pairs = []
    for lpair in [lpk for lpk in utils.locus_pairs[ig_or_tr] if tuple(lpk) in lp_infos]:
        antn_pairs += find_cluster_pairs(lp_infos, lpair, required_keys=required_keys, quiet=quiet, min_cluster_size=min_cluster_size, min_cluster_arg_str=min_cluster_arg_str, debug=debug)
    return antn_pairs

# ----------------------------------------------------------------------------------------
# at least most of the time i shouldn't really need this fcn since the h/l annotation lists are in the same order (and can use the previous fcns), but I dunno sometimes maybe they're not so better to be safe
# if you're sure they're ordered correctly, you can use the previuos two fcns (atm I think i'm only doing that for simulation, where i'm confident the order is correct)
def find_cluster_pairs(lp_infos, lpair, antn_lists=None, required_keys=None, quiet=False, min_cluster_size=None, min_cluster_arg_str='', debug=False):  # the annotation lists should just be in the same order, but after adding back in all the unpaired sequences to each chain they could be a bit wonky
    # ----------------------------------------------------------------------------------------
    def getpids(line):  # return uids of all seqs paired with any seq in <line>
        all_ids = []
        for ip, pids in enumerate(line['paired-uids']):
            if pids is None or len(pids) == 0:
                continue
            elif len(pids) == 1:
                # assert pids[0] not in all_ids  # this is kind of slow, and maybe it's ok to comment it?
                all_ids.append(pids[0])
            else:
                raise Exception('too many paired ids (%d) for %s: %s' % (len(pids), line['unique_ids'][ip], ' '.join(pids)))
        return set(all_ids)
    # ----------------------------------------------------------------------------------------
    if antn_lists is not None:  # can set <antn_lists> if you don't already have lp_infos, and don't need glfos (added this for use in partis.remove_unseeded_seqs(), but it was too slow)
        assert lp_infos is None
        lp_infos = {tuple(lpair) : {'antn_lists' : {l : antn_lists[l] for l in lpair}, 'cpaths' : {ltmp : ClusterPath(partition=[l['unique_ids'] for l in antn_lists[ltmp]]) for ltmp in lpair}}}
    if required_keys is None:
        required_keys = []
    if 'paired-uids' not in required_keys:
        required_keys.append('paired-uids')

    lp_antn_pairs = []
    lpk = tuple(lpair)
    if None in list(lp_infos[lpk].values()):
        return lp_antn_pairs
    h_part, l_part = [sorted(lp_infos[lpk]['cpaths'][l].best(), key=len, reverse=True) for l in lpair]
    h_atn_dict, l_atn_dict = [utils.get_annotation_dict(lp_infos[lpk]['antn_lists'][l], cpath=lp_infos[lpk]['cpaths'][l]) for l in lpair]
    if debug:
        print('  finding cluster pairs for %s partitions with cluster sizes:\n     %s: %s\n     %s: %s' % ('+'.join(lpair), lpair[0], ' '.join(str(len(c)) for c in h_part), lpair[1], ' '.join(str(len(c)) for c in l_part)))
        print('          sizes')
        print('          h   l    l index')
    n_skipped = {k : 0 for k in required_keys + ['zero-len-paired-uids', 'too-small']}
    unpaired_l_clusts = [c for c in l_part]
    for h_clust in h_part:
        if ':'.join(h_clust) not in h_atn_dict:  # note sure if this really avoids it, and i wish i knew how it could be missing, but don't have time right now to track it down
            print('    %s cluster %s missing from h_atn_dict, skipping' % (utils.wrnstr(), ':'.join(h_clust)))
            continue
        h_atn = h_atn_dict[':'.join(h_clust)]
        h_pids = getpids(h_atn)

        if any(k not in h_atn for k in required_keys):  # skip any annotations that are missing any of these keys (atm, only used to skip ones without 'tree-info', which usually means clusters that were smaller than min selection metric cluster size
            for rk in set(required_keys) - set(h_atn):
                n_skipped[rk] += 1
            continue
        if len(h_pids) == 0:
            n_skipped['zero-len-paired-uids'] += 1
            continue

        l_clusts = [c for c in l_part if len(h_pids & set(c)) > 0]
        if len(l_clusts) != 1:
            if not quiet:
                l_overlaps = [h_pids & set(c) for c in l_clusts]
                print('  %s found %d light clusters (rather than 1) for heavy cluster with size %d (%d pids) (overlaps: %s)' % (utils.color('yellow', 'warning'), len(l_clusts), len(h_clust), len(h_pids), ' '.join(str(len(o)) for o in l_overlaps)))
                print('         h clust %s' % ':'.join(utils.color('blue_bkg', u) for u in h_clust))
                for il, lct in enumerate(l_clusts):
                    print('        %s %s' % ('l clusts' if il==0 else '        ', ':'.join(utils.color('blue_bkg' if u in h_pids else None, u) for u in lct)))
            continue
        assert len(l_clusts) == 1
        if ':'.join(l_clusts[0]) not in l_atn_dict:
            print('      %s missing annotation for light chain (size %d, paired with size %d) when finding cluster pairs%s%s' % (utils.color('yellow', 'warning'), len(l_clusts[0]), len(h_clust), ' '+':'.join(l_clusts[0]) if len(l_clusts[0])<30 else '', ' '+':'.join(h_clust) if len(h_clust) < 30 else ''))
            unpaired_l_clusts.remove(l_clusts[0])  # i guess i want to remove it from here? i guess we know who it's paired with, but there's no annotation so we can't do anything with it
            continue
        l_atn = l_atn_dict[':'.join(l_clusts[0])]
        if min_cluster_size is not None and any(len(l['unique_ids']) < min_cluster_size for l in [h_atn, l_atn]):
            n_skipped['too-small'] += 1
            continue
        h_atn['loci'] = [lpair[0] for _ in h_atn['unique_ids']]  # this kind of sucks, but it seems like the best option a.t.m. (see note in event.py)
        l_atn['loci'] = [lpair[1] for _ in l_atn['unique_ids']]
        lp_antn_pairs.append((h_atn, l_atn))
        unpaired_l_clusts.remove(l_clusts[0])
        if debug:
            print('        %3d %3d   %3d' % (len(h_clust), len(l_clusts[0]), l_part.index(l_clusts[0])))
    if len(unpaired_l_clusts) > 0:
        print('    %s: %d unpaired light cluster%s after finding h/l cluster pairs' % ('+'.join(lpair), len(unpaired_l_clusts), utils.plural(len(unpaired_l_clusts))))
        # this is just too verbose atm (and hopefully not necessary?)
        # for lc in unpaired_l_clusts:
        #     if ':'.join(lc) not in l_atn_dict:
        #         print '        %s missing annotation for unpaired light chain when finding cluster pairs: %s' % (utils.color('yellow', 'warning'), ':'.join(lc))
        #         continue
        #     lpids = getpids(l_atn_dict[':'.join(lc)])
        #     hpclusts = [c for c in h_part if len(set(lpids) & set(c)) > 0]0
        #     if len(hpclusts) > 0:  # i think this would mean that the pairing info was non-reciprocal, which probably isn't really possible?
        #         print '       %s unpaired light cluster with size %d overlaps with heavy cluster(s): %s' % (utils.color('yellow', 'warning'), len(lc), ' '.join(str(len(c)) for c in hpclusts))
    keptstr = ' (keeping %d annotation pairs)' % len(lp_antn_pairs)
    if any(n > 0 for k, n in n_skipped.items() if k not in ['too-small', 'zero-len-paired-uids']):
        print('    %s: skipped %d annotations missing required keys (%s)%s' % ('+'.join(lpair), sum(n_skipped.values()), '  '.join('%s: %d'%(k, n) for k, n in sorted(n_skipped.items()) if n>0 and k!='zero-len-paired-uids'), keptstr))
    if n_skipped['zero-len-paired-uids'] > 0:
            print('    %s: skipped %d annotations with zero length paired uids%s' % ('+'.join(lpair), n_skipped['zero-len-paired-uids'], keptstr))
    if n_skipped['too-small'] > 0:
            print('    %s: skipped %d annotations with N h or l ids < %d%s%s' % ('+'.join(lpair), n_skipped['too-small'], min_cluster_size, min_cluster_arg_str, keptstr))
    if debug:
        print('  ')
    return lp_antn_pairs

# ----------------------------------------------------------------------------------------
def gsval(mfo, tch, vname, no_fail=False):
    if tch+'_iseq' not in mfo:  # ick
        return None
    assert tch in 'hlp'  # heavy, light, or fake paired
    cln, iseq = mfo[tch if tch in 'hl' else tch+'_atn'], mfo[tch+'_iseq']
    return utils.antnval(cln, vname, iseq, use_default=no_fail)
# ----------------------------------------------------------------------------------------
def get_did(args, uid, return_contigs=False):
    return utils.get_droplet_id(uid, args.droplet_id_separators, args.droplet_id_indices, return_contigs=return_contigs)

# ----------------------------------------------------------------------------------------
def both_dids(args, mfo):
    return [get_did(args, gsval(mfo, c, 'unique_ids')) for c in 'hl']

# ----------------------------------------------------------------------------------------
def sumv(mfo, kstr, imtp):
    if kstr == 'seq_mtps':  # NOTE this is the sum of utils.get_multiplicity() over identical sequences
        def vfcn(c): return imtp[c][gsval(mfo, c, 'input_seqs_aa')]
    else:
        if kstr in ['seqs', 'naive_seq']:
            def vfcn(c): return utils.pad_seq_for_translation(mfo[c], gsval(mfo, c, kstr))  # maybe don't need this, but safer to have it
        else:
            def vfcn(c): return gsval(mfo, c, kstr)
    kvals = [vfcn(c) for c in 'hl']
    return None if None in kvals else kvals[0] + kvals[1]  # needs to work for both ints and strings

# ----------------------------------------------------------------------------------------
# make fake annotations with h/l seqs smashed together (with N padding between them so l is in frame)
def make_fake_hl_pair_antns(args, antn_pairs):  # maybe better to not require <args>, but whatever, easy to fix later
    # ----------------------------------------------------------------------------------------
    def combid(mfo):  # new uid that combines h+l ids
        _, cids = zip(*[get_did(args, gsval(mfo, c, 'unique_ids'), return_contigs=True) for c in 'hl'])
        dids = both_dids(args, mfo)  # the vast majority of the time they have the same did, so this is just the did, but in simulation, if they're mispaired, they can be different
        if len(set(dids)) == 1:  # if they have the same droplet id (data or correctly paired simulation)
            if not args.is_data or args.use_droplet_id_for_combo_id:  # in simulation the droplet ids should be unique, so we can just use the droplet id as the combined id
                cmbid = dids[0]
            else:  # but in data we can get multiple cells per droplet id
                cmbid = '%s_contig_%s+%s' % (dids[0], cids[0], cids[1])
        else:  # but if they're mispaired in simulation (i.e. have different "droplet ids") then keep all the info
            assert len(set(dids)) == 2
            cmbid = '%s-%s+%s-%s' % (dids[0], mfo['h']['loci'][0], dids[1], mfo['l']['loci'][0])
        return cmbid
    # ----------------------------------------------------------------------------------------
    def update_qti(combo_ids, metric_pairs):
        for cmbid, mfo in zip(combo_ids, metric_pairs):
            hid, lid = [gsval(mfo, c, 'unique_ids') for c in 'hl']
            if any(u in args.queries_to_include for u in (hid, lid)):  # translate uids in args.queries_to_include (should maybe untranslate afterwards?)
                for u in (hid, lid):
                    if u in args.queries_to_include:
                        args.queries_to_include.remove(u)
                args.queries_to_include.append(cmbid)
    # ----------------------------------------------------------------------------------------
    def get_mtpys(metric_pairs):  # NOTE this is the sum of utils.get_multiplicity() over identical sequences
        icl_mtpys = {}
        for c in 'hl':
            seqlist = [gsval(m, c, 'input_seqs_aa') for m in metric_pairs for _ in range(gsval(m, c, 'multipy'))]
            icl_mtpys[c] = {s : seqlist.count(s) for s in set(seqlist)}
        return icl_mtpys
    # ----------------------------------------------------------------------------------------
    def get_pantn(metric_pairs, h_atn, l_atn, all_pair_ids, imtp, tdbg=False):  # return a fake annotation <p_atn> with the sum/joint metrics in it
        # ----------------------------------------------------------------------------------------
        def translate_heavy_tree(htree):
            trns = [(gsval(m, 'h', 'unique_ids'), c) for m, c in zip(metric_pairs, p_atn['unique_ids'])]  # translation from hid to the new combined h+l id we just made
            treeutils.translate_labels(htree, trns, expect_missing=True)
            htree.scale_edges(len(h_atn['seqs'][0]) / float(len(p_atn['seqs'][0])))
            return htree, htree.as_string(schema='newick')
        # ----------------------------------------------------------------------------------------
        def add_unp_seqs(cpkeys):
            # ----------------------------------------------------------------------------------------
            def ambig_seq(tch, t_atn, iseq, aa=False):
                # ----------------------------------------------------------------------------------------
                def aseq():
                    achar = utils.ambiguous_amino_acids[0] if aa else utils.ambig_base
                    return achar * len((l_atn if tch=='h' else h_atn)[tkey][0])  # NOTE doesn't account for indels
                # ----------------------------------------------------------------------------------------
                tkey = 'seqs'
                if aa:
                    tkey += '_aa'
                tseq = t_atn[tkey][iseq]
                hseq, lseq = (tseq, aseq()) if tch=='h' else (aseq(), tseq)
                if not aa:
                    hseq, lseq = [utils.pad_seq_for_translation(l, s) for l, s in zip((h_atn, l_atn), (hseq, lseq))]  # NOTE kind of duplicates code in sumv()
                return hseq + lseq
            # ----------------------------------------------------------------------------------------
            def print_amb_seqs(all_unp_ids):  # print all seqs after adding the unpaired ones (which will have ambiguous sections)
                for tstr, aa in zip(('', '_aa'), (False, True)):
                    print('       naive %s' % p_atn['naive_seq%s'%tstr])
                    for u, s in zip(p_atn['unique_ids'], p_atn['seqs%s'%tstr]):
                        utils.color_mutants(p_atn['naive_seq%s'%tstr], s, amino_acid=aa, extra_str='      %s' % (utils.color('purple', '  unp. ') if u in all_unp_ids else '       '), post_str=' '+u, print_result=True, only_print_seq=True)
            # ----------------------------------------------------------------------------------------
            mfo_ids = [gsval(m, c, 'unique_ids') for m in metric_pairs for c in 'hl']
            all_unp_ids, n_unp_added = [], {c : 0 for c in 'hl'}
            if tdbg:
                print('  iclust %d: adding unpaired seqs for paired selection metrics' % iclust)
            for tch, och, t_atn in zip('hl', 'lh', (h_atn, l_atn)):
                unp_ids = [u for u in t_atn['unique_ids'] if u not in mfo_ids]
                if tdbg:
                    print('    %s: %s' % (tch, unp_ids))
                for tid in unp_ids:
                    icseq = t_atn['unique_ids'].index(tid)  # index in h/l antn
                    ipseq = len(p_atn['seqs'])  # index in fake paired annotation (to which we're adding the unpaired seq)
                    p_atn['seqs'].append(ambig_seq(tch, t_atn, icseq))
                    p_atn['input_seqs'].append(p_atn['seqs'][ipseq])
                    p_atn['seqs_aa'].append(ambig_seq(tch, t_atn, icseq, aa=True))
                    p_atn['shm_aa'].append(utils.antnval(t_atn, 'shm_aa', icseq))
                    tmpkeys = ['unique_ids', 'n_mutations', 'mut_freqs'] + cpkeys
                    if 'multiplicities' in t_atn:
                        tmpkeys.append('multiplicities')
                    for tk in [k for k in tmpkeys if k in t_atn]:
                        assert len(p_atn[tk]) == ipseq  # if somehow duplicates sneak into <tmpkeys> again this'll trigger
                        p_atn[tk].append(t_atn[tk][icseq])
                    p_atn['has_shm_indels'].append(False)  # ick ick ick
                all_unp_ids += unp_ids
                n_unp_added[tch] += len(unp_ids)
            utils.check_per_seq_lengths(p_atn)
            if tdbg:
                print_amb_seqs(all_unp_ids)
            print('    added unpaired seqs to fake paired annotation: %s %d  %s %d' % (utils.locstr(h_atn['loci'][0]), n_unp_added['h'], utils.locstr(l_atn['loci'][0]), n_unp_added['l']))
        # ----------------------------------------------------------------------------------------
        p_atn = {'is_fake_paired' : True, 'invalid' : True}  # make a new fake annotation for the sequences that are in both h+l (they're not really 'invalid', but they *are* fake, and e.g. the indel info is wrong, so seems safer to call the 'invalid')
        p_atn['unique_ids'] = [combid(m) for m in metric_pairs]
        if args.queries_to_include is not None:
            update_qti(p_atn['unique_ids'], metric_pairs)
        if any(u in all_pair_ids for u in p_atn['unique_ids']):
            raise Exception('tried to add duplicate uid(s) %s when making paired annotation' % [u for u in p_atn['unique_ids'] if u in all_pair_ids])
        all_pair_ids |= set(p_atn['unique_ids'])
        p_atn['seqs'] = [sumv(m, 'seqs', imtp) for m in metric_pairs]
        p_atn['input_seqs'] = [s for s in p_atn['seqs']]  # NOTE do *not* let 'seqs' and 'input_seqs' point to the same list (we only need 'input_seqs' since they're what gets written to the output file)
        for tch in 'hl':
            p_atn['%s_seqs'%tch] = [gsval(m, tch, 'seqs') for m in metric_pairs]
        for mfo, pseq in zip(metric_pairs, p_atn['seqs']):
            total_hlen = 0  # total len of igh seq
            for tch in 'hl':
                _, npads = utils.pad_seq_for_translation(mfo[tch], gsval(mfo, tch, 'seqs'), return_n_padded=True)  # maybe don't need this, but safer to have it
                cbounds = [mfo[tch]['codon_positions'][r] for r in 'vj']
                offset = npads[0] + total_hlen
                p_atn['%s_offset'%tch] = offset  # index of first position in <tch>
                p_atn['%s_cdr3_bounds'%tch] = [b + offset for b in cbounds]
                p_atn['%s_frame'%tch] = utils.get_frame(mfo[tch])
                assert utils.get_cdr3_seq(mfo[tch], mfo[tch+'_iseq']) == pseq[cbounds[0] + offset : cbounds[1] + offset + 3]  # make sure the cdr3 in the padded seq (with the offset) gives the same bit as from the h/l annotation
                if tch == 'h':
                    total_hlen += npads[0] + len(gsval(mfo, tch, 'seqs')) + npads[1]
        if len(metric_pairs) == 0:
            return p_atn
        p_atn['seqs_aa'] = [sumv(m, 'seqs_aa', imtp) for m in metric_pairs]
        p_atn['naive_seq'] = sumv(metric_pairs[0], 'naive_seq', imtp)
        p_atn['naive_seq_aa'] = sumv(metric_pairs[0], 'naive_seq_aa', imtp)  # NOTE it's *really* important you don't end up translating the sum'd naive seq since i don't think they necessarily get concat'd in frame
        p_atn['n_mutations'] = [sumv(m, 'n_mutations', imtp) for m in metric_pairs]
        p_atn['shm_aa'] = [sumv(m, 'shm_aa', imtp) for m in metric_pairs]
        p_atn['mut_freqs'] = [n / float(len(s)) for n, s in zip(p_atn['n_mutations'], p_atn['seqs'])]
        p_atn['has_shm_indels'] = [False for _ in metric_pairs]  # ick ick ick
        # NOTE if you add a key here, it also has to be added below in the args.add_unpaired_seqs_to_fake_paired_annotations block
        if 'multiplicities' in h_atn:  # <h_atn> is the same as m['h'], i should really settle on one of them
            h_mults, l_mults = [[utils.get_multiplicity(m[c], uid=gsval(m, c, 'unique_ids')) for m in metric_pairs] for c in 'hl']
            if h_mults != l_mults:
                raise Exception('h and l multiplicities not the same:\n    %s\n    %s' % (h_mults, l_mults))
            p_atn['multiplicities'] = h_mults
        cpkeys = ['affinities' if args.affinity_key is None else args.affinity_key]  # per-seq keys to copy from h_atn (NOTE ignores l_atn)
        if not args.is_data:
            assert not args.add_unpaired_seqs_to_fake_paired_annotations  # not sure if it makes sense? in any case i'm pretty sure the tree wouldn't be right, and some other things would probably have to change
            _, p_atn['tree'] = translate_heavy_tree(treeutils.get_dendro_tree(treestr=h_atn['tree']))
            cpkeys.append('min_target_distances')
        if args.meta_info_key_to_color is not None and args.meta_info_key_to_color:
            cpkeys.append(args.meta_info_key_to_color)
        if args.meta_info_to_emphasize is not None:
            cpkeys += [k for k in args.meta_info_to_emphasize.keys()]
        def addme(k): return all(k in m['h'] and k in m['l'] and gsval(m, 'h', k)==gsval(m, 'l', k) for m in metric_pairs)
        cpkeys += [k for k in utils.input_metafile_keys.values() if addme(k)]  # input meta keys that are in both h and l annotations and equal in value for all mfos
        cpkeys = list(set(cpkeys))
        cpkeys = [k for k in cpkeys if k not in p_atn]
        assert len(cpkeys) == len(set(cpkeys))  # any duplicates will fuck up add_unp_seqs()
        assert all(k not in p_atn for k in cpkeys)  # copying any key that's already there (e.g. n_mutations) is also *really* bad
        for tk in [k for k in cpkeys if k in h_atn]:
            p_atn[tk] = [h_atn[tk][m['h_iseq']] for m in metric_pairs]
        if args.add_unpaired_seqs_to_fake_paired_annotations:
            add_unp_seqs(cpkeys)
        p_atn['fv_insertion'] = ''  # this stuff for left side of v is only needed so when we add aa seqs corresponding to any inferred ancestral seqs it doesn't crash trying to pad (we don't want it to pad, since it's already padded above)
        p_atn['v_5p_del'] = 0
        for iseq, mfo in enumerate(metric_pairs):
            mfo['p_atn'] = p_atn
            mfo['p_iseq'] = iseq
        return p_atn
    # ----------------------------------------------------------------------------------------
    pair_antns, mtpys, all_pair_ids = [], {}, set()
    mpfo_lists, pair_antns = [[None for _ in antn_pairs] for _ in range(2)]
    for iclust, (h_atn, l_atn) in enumerate(antn_pairs):
        for ltmp in (h_atn, l_atn):
            utils.add_seqs_aa(ltmp)
            utils.add_naive_seq_aa(ltmp)
        metric_pairs = []
        for hid, pids in zip(h_atn['unique_ids'], h_atn['paired-uids']):
            if pids is None or len(pids) == 0:  # should only have the latter now (set with .get() call in rewrite_input_metafo())
                continue
            lid = pids[0]
            if lid not in l_atn['unique_ids']:
                print('  paired light id %s missing' % lid)
                continue
            mpfo = {'iclust' : iclust, 'seqtype' : 'observed'}
            for tch, uid, ltmp in zip(('h', 'l'), (hid, lid), (h_atn, l_atn)):
                mpfo[tch] = ltmp
                mpfo[tch+'_iseq'] = ltmp['unique_ids'].index(uid)
            metric_pairs.append(mpfo)
        mpfo_lists[iclust] = metric_pairs
        mtpys[iclust] = get_mtpys(metric_pairs)
        pair_antns[iclust] = get_pantn(metric_pairs, h_atn, l_atn, all_pair_ids, mtpys[iclust])

    return pair_antns, mpfo_lists, mtpys

# ----------------------------------------------------------------------------------------
def remove_pair_info_from_bulk_data(outfos, metafos, bulk_data_fraction):
    if len(outfos) % 2 != 0:
        print('  %s outfos len not a multiple of 2' % utils.wrnstr())  # i think everybody's partner should still be in there
    n_to_remove = int(bulk_data_fraction * len(outfos) / 2.)
    n_removed, single_cell_ids = 0, list(metafos)  # <single_cell_ids> is the list of all uids that still have pair info
    while n_removed < n_to_remove:
        uid_to_remove = numpy.random.choice(single_cell_ids)
        pid_to_remove = utils.get_single_entry(metafos[uid_to_remove]['paired-uids'])
        for utmp in [uid_to_remove, pid_to_remove]:
            metafos[utmp]['paired-uids'] = []
            single_cell_ids.remove(utmp)
        n_removed += 1
    print('  removed pair info for %d / %d = %.2f sequence pairs' % (n_to_remove, len(outfos) / 2., n_to_remove / float(len(outfos)/2.)))

# ----------------------------------------------------------------------------------------
def apportion_cells_to_droplets(outfos, metafos, mean_cells_per_droplet, constant_n_cells=False, debug=False):
    n_droplets = max(1, int(0.5 * float(len(outfos)) / mean_cells_per_droplet))  # (randomly) apportion cells among this many droplets (0.5 is because <outfos> includes both heavy and light sequences)
    droplet_ids = [[] for _ in range(n_droplets)]  # list of sequence ids for each droplet
    sfo_dict = {s['name'] : s for s in outfos}  # temp, to keep track of who still needs apportioning (but we do modify its sfos, which are shared with <outfos>)
    while len(sfo_dict) > 0:
        tid = next(iter(sfo_dict))
        if constant_n_cells:
            for glen, igroup in itertools.groupby(sorted(range(len(droplet_ids)), key=lambda i: len(droplet_ids[i])), key=lambda i: len(droplet_ids[i])):  # sort droplets by number of cells, then choose any droplet that is among those with the fewest number of cells (e.g. if 7 droplets have 0 cells, choose randomly from among those 7)
                ichoices = list(igroup)
                break
            if len(ichoices) == 0:
                ichoices = list(range(len(droplet_ids)))
            idrop = numpy.random.choice(ichoices)
        else:
            idrop = numpy.random.choice(range(len(droplet_ids)))
        droplet_ids[idrop] += [tid] + metafos[tid]['paired-uids']  # add <tid> plus its paired ids to this drop (note that these are the original/correct paired ids, which is what we want)
        for uid in [tid] + metafos[tid]['paired-uids']:
            sfo_dict[uid]['droplet-ids'] = droplet_ids[idrop]
            del sfo_dict[uid]
    for sfo in outfos:
        metafos[sfo['name']]['paired-uids'] = sorted([u for u in sfo['droplet-ids'] if u != sfo['name']])
    if debug:
        print('    droplets:')
        for idrop, uids in enumerate(droplet_ids):
            print('    %5d  %s' % (idrop, ' '.join(sorted(uids))))
    print('  apportioned %d seqs among %d droplets (mean/2 %.1f): %s' % (len(outfos), n_droplets, numpy.mean([len(d) for d in droplet_ids]) / 2., ' '.join(str(len(d)) for d in droplet_ids)))

# ----------------------------------------------------------------------------------------
def remove_reads_from_droplets(outfos, metafos, fraction_of_reads_to_remove, debug=False):
    n_to_remove = int(fraction_of_reads_to_remove * len(outfos))
    ifos_to_remove = numpy.random.choice(range(len(outfos)), size=n_to_remove, replace=False)
    uids_to_remove = set(outfos[ifo]['name'] for ifo in ifos_to_remove)
    pids_of_removed_ifos = set()  # keep track of all pids of each removed id, so we can (quickly) remove the removed ids from each pid's pair info
    for rid in uids_to_remove:
        pids_of_removed_ifos |= set(metafos[rid]['paired-uids'])
        del metafos[rid]
    pids_of_removed_ifos -= uids_to_remove
    for pid in pids_of_removed_ifos:
        metafos[pid]['paired-uids'] = sorted(set(metafos[pid]['paired-uids']) - uids_to_remove)
    outfos = [outfos[ifo] for ifo in range(len(outfos)) if ifo not in ifos_to_remove]
    print('  removed %d / %d = %.2f seqs from outfos (leaving %d / %d unpaired)' % (n_to_remove, len(outfos) + n_to_remove, n_to_remove / float(len(outfos) + n_to_remove), len([ofo for ofo in outfos if len(metafos[ofo['name']]['paired-uids'])==0]), len(outfos)))
    if debug:
        print('     removed ids: %s' % ' '.join(sorted(uids_to_remove)))
        print('     left unpaired: %s' % ' '.join([ofo['name'] for ofo in outfos if len(metafos[ofo['name']]['paired-uids'])==0]))
    return outfos, sorted(uids_to_remove)  # need to sort in order to get the same order on different runs

# ----------------------------------------------------------------------------------------
# merge together info from all loci into <outfos> and <metafos>
def get_combined_outmetafos(antn_lists, extra_meta_headers=None):  # <extra_meta_headers>: pass in the annotation (usually plural) version
    outfos, metafos = [], {}
    for ltmp in antn_lists:
        for tline in antn_lists[ltmp]:
            for iseq, (uid, seq, pids) in enumerate(zip(tline['unique_ids'], tline['input_seqs'], tline['paired-uids'])):
                outfos.append({'name' : uid, 'seq' : seq})
                metafos[uid] = {'locus' : ltmp, 'paired-uids' : pids}
                if extra_meta_headers is not None:
                    if any(h not in list(utils.input_metafile_keys.values()) for h in extra_meta_headers):  # they have to be in utils.input_metafile_keys so we know what the "meta file" (i.e. usually singular) version is
                        raise Exception('extra meta headers have to be in utils.input_metafile_keys, but got unknown ones: %s' % list((h for h in extra_meta_headers if h not in list(utils.input_metafile_keys.values()))))
                    metafos[uid].update({utils.reversed_input_metafile_keys[h] : tline[h][iseq] for h in extra_meta_headers if h in tline})
    return outfos, metafos

# ----------------------------------------------------------------------------------------
# write fasta and meta file with all simulation loci together
def write_combined_fasta_and_meta(fastafname, metafname, outfos, metafos, write_locus_files=False):
    print('    writing combined h/l info to fasta and meta files in %s' % os.path.dirname(fastafname))
    utils.mkdir(fastafname, isfile=True)
    utils.write_fasta(fastafname, outfos)
    utils.jsdump(metafname, metafos)
    if write_locus_files:  # write also single-locus fasta fiels
        def kfcn(s): return metafos[s['name']]['locus']
        for ltmp, lfos in itertools.groupby(sorted(outfos, key=kfcn), key=kfcn):
            utils.write_fasta('%s/%s.fa' % (os.path.dirname(fastafname), ltmp), lfos)

# ----------------------------------------------------------------------------------------
def modify_simu_pair_info(args, outfos, metafos, lp_infos, concat_lpfos):
    # ----------------------------------------------------------------------------------------
    def update_lpf(lpfos, ltmp, uids_to_remove):
        if args.fraction_of_reads_to_remove is not None:
            partition = [[u for u in c if u not in uids_to_remove] for c in lpfos['cpaths'][ltmp].best()]
            lpfos['cpaths'][ltmp] = ClusterPath(partition=partition, seed_unique_id=lpfos['cpaths'][ltmp].seed_unique_id)
            i_atns_to_remove = []
            for iatn, atn in enumerate(lpfos['antn_lists'][ltmp]):
                atn_ids_to_remove = sorted(set(atn['unique_ids']) & set(uids_to_remove))
                if len(atn_ids_to_remove) == len(atn['unique_ids']):
                    i_atns_to_remove.append(iatn)
                elif len(atn_ids_to_remove) > 0:
                    utils.restrict_to_iseqs(atn, [i for i, u in enumerate(atn['unique_ids']) if u not in atn_ids_to_remove], lpfos['glfos'][ltmp])
            if len(i_atns_to_remove) > 0:
                lpfos['antn_lists'][ltmp] = [a for i, a in enumerate(lpfos['antn_lists'][ltmp]) if i not in i_atns_to_remove]
        for atn in lpfos['antn_lists'][ltmp]:
            atn['paired-uids'] = [metafos[u]['paired-uids'] for u in atn['unique_ids']]
    # ----------------------------------------------------------------------------------------
    if args.bulk_data_fraction is not None:
        remove_pair_info_from_bulk_data(outfos, metafos, args.bulk_data_fraction)
    if args.mean_cells_per_droplet is not None:
        apportion_cells_to_droplets(outfos, metafos, args.mean_cells_per_droplet, constant_n_cells=args.constant_cells_per_droplet)
    uids_to_remove = None
    if args.fraction_of_reads_to_remove is not None:
        outfos, uids_to_remove = remove_reads_from_droplets(outfos, metafos, args.fraction_of_reads_to_remove)
    if any(a is not None for a in [args.bulk_data_fraction, args.mean_cells_per_droplet, args.fraction_of_reads_to_remove]):
        for lpair, lpfos in lp_infos.items():
            if lpfos['glfos'] is None:  # no simulation events for this locus pair
                continue
            for ltmp in lpair:
                update_lpf(lpfos, ltmp, uids_to_remove)
        for ltmp in concat_lpfos['glfos']:
            update_lpf(concat_lpfos, ltmp, uids_to_remove)
    return outfos

# ----------------------------------------------------------------------------------------
# rename all seqs in the light chain partition to their paired heavy chain uid (also change uids in the light chain annotations). Note that pairings must, at this stage, be unique.
def translate_paired_uids(ploci, init_partitions, antn_lists):
    # first make a map from each light chain id to its paired heavy chain id
    h_paired_uids = {}
    for hline in antn_lists[ploci['h']]:  # looping over the partitions here was identical, at least when I checked it
        for h_id, pids in zip(hline['unique_ids'], hline['paired-uids']):
            if len(pids) == 0:
                raise Exception('no paired uids for %s' % h_id)  # everybody has to have exactly one paired id at this point
            elif len(pids) > 1:
                raise Exception('multiple paired uids %s for %s sequence %s' % (' '.join(pids), ploci['h'], h_id))
            h_paired_uids[pids[0]] = h_id
    cpath = ClusterPath(partition=init_partitions['l'])
    l_translations = utils.translate_uids(antn_lists[ploci['l']], trns=h_paired_uids, cpath=cpath, failstr='paired uids')
    init_partitions['l'] = cpath.best()
    # i'm too chicken to delete the old way:
    # # then go through the light chain annotations + partition swapping names
    # l_translations = {}
    # for lline in antn_lists[ploci['l']]:
    #     for iseq, l_id in enumerate(lline['unique_ids']):
    #         if l_id not in h_paired_uids:
    #             raise Exception('no paired uids for %s' % l_id)  # everybody has to have exactly one paired id at this point
    #         lline['unique_ids'][iseq] = h_paired_uids[l_id]
    #         l_translations[h_paired_uids[l_id]] = l_id  # so we can go back to <l_id> afterwards
    # if len(h_paired_uids) > 0:
    #     init_partitions['l'] = [[h_paired_uids.get(u, u) for u in c] for c in init_partitions['l']]
    return l_translations

# ----------------------------------------------------------------------------------------
# reverse action of previous fcn
# TODO should switch to utils.translate_uids() here
def untranslate_pids(ploci, init_partitions, antn_lists, l_translations, joint_partitions, antn_dict):
    for lline in antn_lists[ploci['l']]:
        lline['unique_ids'] = [l_translations.get(u, u) for u in lline['unique_ids']]
    init_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in init_partitions['l']]
    joint_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in joint_partitions['l']]

    # remove h ids from the l joint partition and vice versa
    all_loci = {u : l for ants in antn_lists.values() for antn in ants for u, l in zip(antn['unique_ids'], antn['loci'])}  # it seems like i should have this info somewhere already, but maybe not?
    for tch in joint_partitions:
        joint_partitions[tch] = [[u for u in c if all_loci[u]==ploci[tch]] for c in joint_partitions[tch]]  # i think they should all be in <all_loci>

    for ch in antn_dict:  # maybe this this still only gets used for dbg? but still should properly fix it
        antn_dict[ch] = utils.get_annotation_dict(antn_lists[ploci[ch]])

# ----------------------------------------------------------------------------------------
def pair_unpaired_seqs_with_paired_family(ploci, unpaired_seqs, cluster_pairs, antn_lists, debug=False):
    # note that at this point these are *single* chain annotations, and/since we've only just made the joint partition
    # ----------------------------------------------------------------------------------------
    def get_pids(atn, uid):
        pds = atn['paired-uids'][atn['unique_ids'].index(uid)]
        assert len(pds) in [0, 1]  # just to make sure
        return pds
    # ----------------------------------------------------------------------------------------
    def tcol(pds, oclust):
        if len(pds) == 0:
            return utils.color('blue', '-')
        elif pds[0] not in oclust:  # paired id isn't in opposite chain cluster? not sure but maybe could happen
            return utils.color('yellow', '?')
        else:
            return utils.color('green', 'x')
    # ----------------------------------------------------------------------------------------
    def cstrs(cpair, pd_dict):
        rstrs = [' '.join(tcol(pd_dict[u], oc) for u in tc) for tc, oc in zip(cpair, reversed(cpair))]
        rstrs = ['%s%s' % (tstr, ' '*(mlen - 2*len(tclust))) for mlen, tstr, tclust in zip(max_lens, rstrs, cpair)]
        return rstrs
    # ----------------------------------------------------------------------------------------
    def zstr(val):
        return '' if val==0 else str(val)
    # ----------------------------------------------------------------------------------------
    n_no_paired = 0
    if debug:
        max_lens = [2*max(len(c) for c in clist) for clist in zip(*cluster_pairs)]
        print('  pairing unpaired seqs from %d cluster pairs (N fixed \'both\': paired with unique unpaired seq from other chain; \'extra\' same, but not unique; %s paired with already-paired seq' % (len(cluster_pairs), utils.color('blue', '*')))
        print('        size     recipr    N fixed         %s   %s' % tuple(utils.wfmt(s, 2*mlen, jfmt='-') for s, mlen in zip(('heavy', 'light'), max_lens)))
        print('       h     l   paired    h     l       %s%s   %s%s' % tuple(utils.wfmt(s, mlen, jfmt='-') for mlen in max_lens for s in ['before', 'after']))
    # ----------------------------------------------------------------------------------------
    all_atns = {u : utils.get_single_entry([l for l in antn_lists[ploci[tch]] if u in l['unique_ids']]) for cp in cluster_pairs for tch, c in zip('hl', cp) for u in c}
    for cpair in cluster_pairs:
        h_atns, l_atns = [[all_atns[u] for u in tclust] for c, tclust in zip('hl', cpair)]
        all_pids = {u : get_pids(l, u) for tclust, alist in zip(cpair, (h_atns, l_atns)) for u, l in zip(tclust, alist)}
        ok_ids = [[u for u in tc if len(all_pids[u])==1 and all_pids[u][0] in oc] for tc, oc in zip(cpair, reversed(cpair))]  # list of correctly paired uids for each chain (note that this is *not* mutually exclusive to <unp_ids>) [tc: this cluster, oc: other cluster]
        n_recip_paired = utils.get_single_entry(list(set(len(ulist) for ulist in ok_ids)))
        if n_recip_paired == 0:
            n_no_paired += 1
            continue
        unp_ids = [[u for u in tc if len(all_pids[u])==0] for tc, oc in zip(cpair, reversed(cpair))]  # could i guess just use <unpaired_seqs> as well
        for tch, unlist, oklist in zip('hl', unp_ids, ok_ids):
            for unid in unlist:
                nearid = unpaired_seqs[ploci[tch]][unid]['nearest-paired']  # pair all unpaired seqs with their nearest paired seq (which i think makes sense since the paired seqs are the ones that are providing the actual information about what family the unpaired seqs are paired with)
                all_atns[unid]['paired-uids'][all_atns[unid]['unique_ids'].index(unid)] = [all_pids[nearid][0]]
        if debug:
            after_pids = {u : get_pids(l, u) for tclust, alist in zip(cpair, (h_atns, l_atns)) for u, l in zip(tclust, alist)}
            bcstrs, acstrs = [cstrs(cpair, pdlist) for pdlist in (all_pids, after_pids)]
            print('    %4d  %4d  %4s    %4s  %4s       %s %s    %s %s'  % (len(cpair[0]), len(cpair[1]), zstr(n_recip_paired), zstr(len(unp_ids[0])), zstr(len(unp_ids[1])), bcstrs[0], acstrs[0], bcstrs[1], acstrs[1]))
# ----------------------------------------------------------------------------------------
# old version that first pairs unpaired seqs in order, then apportions the leftovers (i think it's worse, but it was a fair bit of work so don't want to delete yet):
        # unp_ids = [[u for u in tc if len(all_pids[u])==0] for tc, oc in zip(cpair, reversed(cpair))]  # keep the order the same, since it seems like similar seqs end up in similar order, probably due to seq similarity from single chain clustering?
        # for uid_pair in zip(*unp_ids):  # first do the ones that're matched up
        #     for tch, tid, oid in zip('hl', uid_pair, reversed(uid_pair)):
        #         all_atns[tid]['paired-uids'][all_atns[tid]['unique_ids'].index(tid)] = [oid]
        # (ish, iln), (shorter_unids, longer_unids) = zip(*sorted(enumerate(unp_ids), key=lambda x: len(x[1])))
        # shid = shorter_unids[-1] if len(shorter_unids) > 0 else ok_ids[ish][0]  # if there are some unpaired uids from the shorter chain, use those; otherwise use the first paired one (note that this means in the latter case that we're pairing with a seq from single cell data, which is not going to be a correct pairing, but it *will* likely result in a functional ab)
        # for eid in longer_unids[len(shorter_unids) : ]:
        #     all_atns[eid]['paired-uids'][all_atns[eid]['unique_ids'].index(eid)] = [shid]
        # if debug:
        #     n_bt_fixed, n_ex_fixed = len(shorter_unids), len(longer_unids) - len(shorter_unids)
        #     after_pids = {u : get_pids(l, u) for tclust, alist in zip(cpair, (h_atns, l_atns)) for u, l in zip(tclust, alist)}
        #     bcstrs, acstrs = [cstrs(cpair, pdlist) for pdlist in (all_pids, after_pids)]
        #     print '    %4d  %4d  %4s    %4s  %4s%s       %s %s    %s %s      %s    %s'  % (len(cpair[0]), len(cpair[1]), zstr(n_recip_paired), zstr(n_bt_fixed), zstr(n_ex_fixed), ' ' if len(longer_unids)==0 or len(shorter_unids)>0 else utils.color('blue', '*'), bcstrs[0], acstrs[0], bcstrs[1], acstrs[1], ' '.join(cpair[0]), ' '.join(cpair[1]))
# ----------------------------------------------------------------------------------------

    if n_no_paired > 0:
        print('    skipped %d family pairs with no reciprocally paired seqs' % n_no_paired)

# ----------------------------------------------------------------------------------------
def remove_badly_paired_seqs(ploci, outfos, keep_all_unpaired=False, debug=False):  # remove seqs paired with the other/wrong light chain, as well as those with no pairing info (the latter we keep track of so we can insert them later into the right final cluster)
    # ----------------------------------------------------------------------------------------
    def add_unpaired(cline, iseq, uid, paired_iseqs):
        nearest_uid = None
        if len(paired_iseqs) > 0:
            nearest_uid = 'NEAR'  # now that we have process_unpaired(), we don't use the actual distance that was calculated here, we just need to know if there were any paired seqs in the family (note that this could probably be cleaned up some more, e.g. synchronize definitions of 'paired' vs 'to keep' seqs, but I don't want to change things atm) 
            # sorted_hdists = sorted([(cline['unique_ids'][i], utils.hamming_distance(cline['seqs'][i], cline['seqs'][iseq])) for i in paired_iseqs if i != iseq], key=operator.itemgetter(1))
            # nearest_uid = sorted_hdists[0][0] if len(sorted_hdists) > 0 else None
        unpaired_seqs[cline['loci'][iseq]][uid] = {'nearest' : nearest_uid, 'nearest-paired' : None, 'single-chain-family' : None}  # unless there's no other seqs in the cluster, attach it to the nearest seq by hamming distance NOTE this is 'nearest' while the below is 'nearest-paired', although that isn't quite accurate names
    # ----------------------------------------------------------------------------------------
    def process_unpaired(cline, locus, iseqs_to_keep, paired_iseqs, compare_hamming=False):
        # ----------------------------------------------------------------------------------------
        def print_cf_ham(iun, sorted_hdists):  # this is just to print some dbg to compare to the old hamming distance, so doesn't really belong here long term
            old_sorted_hdists = sorted([(cline['unique_ids'][i], utils.hamming_distance(cline['seqs'][i], cline['seqs'][iun])) for i in paired_itokeep], key=operator.itemgetter(1))
            print(old_sorted_hdists)
            print(sorted_hdists)
            if [u for u, d in old_sorted_hdists][0] != [u for u, d in sorted_hdists][0]:
                print('  %s' % utils.color('red', 'YEP'))
            utils.print_reco_event(cline, extra_print_keys=['paired-uids', 'mut_positions'])
            print('')
        # ----------------------------------------------------------------------------------------
        mut_positions = utils.get_mut_positions(cline)
        cline['mut_positions'] = mut_positions
        for iun, unid in [(i, u) for i, u in enumerate(cline['unique_ids']) if u in unpaired_seqs[locus]]:  # loop over unpaired seqs in <cline>
            if unpaired_seqs[locus][unid]['nearest'] is None:  # no paired seqs in family (maybe it's a singleton, maybe not)
                unpaired_seqs[locus][unid]['single-chain-family'] = cline['unique_ids']
                continue
            nearest_pid = None
            if len(iseqs_to_keep) > 0:
                # used to sort by hamming distance, but it was too slow, so now keep track of muted positions in each seq, and do this approximate hamming calculation (only difference i think is when the two sequences mutate at the same position but to different bases, which is rare enough idgaf)
                sorted_hdists = sorted([(cline['unique_ids'][i], len(mut_positions[i]) + len(mut_positions[iun]) - 2 * len(set(mut_positions[i]) & set(mut_positions[iun]))) for i in iseqs_to_keep], key=operator.itemgetter(1))
                nearest_pid = sorted_hdists[0][0] if len(sorted_hdists) > 0 else None
                if compare_hamming:  # compare to old way of actually calculating hamming distance
                    print_cf_ham(iun, sorted_hdists)
            unpaired_seqs[locus][unid]['nearest-paired'] = nearest_pid
            if nearest_pid is None:
                unpaired_seqs[locus][unid]['single-chain-family'] = cline['unique_ids']
    # ----------------------------------------------------------------------------------------
    cpaths, antn_lists, glfos = [outfos[k] for k in ['cpaths', 'antn_lists', 'glfos']]
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l]) for l in antn_lists}
    all_loci = {u : l for l, ants in antn_lists.items() for antn in ants for u in antn['unique_ids']}  # this includes the heavy ones, which we don't need, but oh well
    all_pids = {u : pids[0] for alist in antn_lists.values() for l in alist for u, pids in zip(l['unique_ids'], l['paired-uids']) if len(pids)==1}  # uid : pid for all uid's that have a single unique pid (which should be all of them, since we just ran pair cleaning -- otherwise we crash below) (I'm pretty sure that the partition implied by the annotations is identical to the one in <cpaths>, and it's nice to loop over annotations for this)
    unpaired_seqs = {l : {} for l in ploci.values()}  # map for each locus from the uid of each seq with no (or non-reciprocal) pairing info to the nearest sequence in its family (after merging partitions we'll insert it into the family that this nearest seq ended up in)
    lp_cpaths, lp_antn_lists = {}, {}
    total_unpaired_discarded = {c : 0 for c in ploci}  # just for debug
    print('    removing badly + unpaired seqs%s' % ('\n' if debug else ': '), end=' ')
    sys.stdout.flush()
    if debug:
        print('          N     N      N      no   other  non-  discarding    original')
        print('        before kept removed  info  light recip   unpaired      cluster')
    for tch in sorted(ploci):
        new_partition, new_antn_list = [], []
        for iclust, cluster in enumerate(cpaths[ploci[tch]].best()):
            cline = antn_dicts[ploci[tch]][':'.join(cluster)]
            paired_iseqs = [i for i, pds in enumerate(cline['paired-uids']) if len(pds) > 0]
            iseqs_to_remove = []
            n_no_info, n_other_light, n_non_reciprocal = 0, 0, 0  # just for dbg NOTE n_other_light are the only ones we *really* want to remove, since they're h seqs paired with the other light chain (e.g. we're doing h+k and they're paired with igl), whereas the other two categories we eventually want to re-add since we're not sure who they're paired with
            unpaired_to_add = []
            for iseq, uid in enumerate(cline['unique_ids']):
                pids = cline['paired-uids'][iseq]
                if len(pids) == 0:  # no pairing info
                    iseqs_to_remove.append(iseq)
                    unpaired_to_add.append((iseq, uid))
                    n_no_info += 1
                elif len(pids) > 1:  # should've all been removed by pair info cleaning
                    raise Exception('multiple paired uids for \'%s\': %s' % (uid, pids))
                else:
                    # print('        ', uid, tch, pids, all_loci[utils.get_single_entry(pids)], ploci['l'])
                    if tch == 'h' and all_loci[utils.get_single_entry(pids)] != ploci['l']:  # if it's the other light chain (e.g. igl when we're doing h+k)
                        iseqs_to_remove.append(iseq)
                        n_other_light += 1
                    else:  # also remove any non-reciprocal pairings (I think this will still miss any whose partner was removed) NOTE it would be nice to enforce reciprocal pairings in pair info cleaning, but atm i think we can't look at both chains at once in that fcn UPDATE i think we do this now
                        if all_pids[uid] not in all_pids or all_pids[all_pids[uid]] != uid:  # if uid's pid isn't in all_pids, or if it is but it's a different uid
                            iseqs_to_remove.append(iseq)
                            unpaired_to_add.append((iseq, uid))
                            n_non_reciprocal += 1
            # NOTE this means (more or less) we discard any unpaired seqs that are in families with no paired seqs, which i think is what we want to do?
            if len(unpaired_to_add) > 0:
                discarding_unpaired = True  # just for debug
                # first clause: if --keep-all-unpaired was set, and if all seqs in the cluster are actually unpaired (as opposed to a mix of unpaired and other chain, in which case keeping them here would screw them up when we later do the other chain)
                # second clause: if there's paired seqs in the cluster (i.e. if the whole cluster isn't n_no_info or n_other_light), we want to add the unpaired ones as unpaired seqs (i.e. if n_no_info + n_other_light *equals* len(cluster), we *don't* want to add them, since either they're all unpaired, or they should get added in/during the other light chain, and if they're added to both, they'll get deduplicated poorly later)
                if (keep_all_unpaired and n_no_info==len(cluster)) or n_no_info + n_other_light < len(cluster):
                    discarding_unpaired = False
                    for upi, upid in unpaired_to_add:  # note that n_other_light seqs aren't in <unpaired_to_add>
                        add_unpaired(cline, upi, upid, paired_iseqs)
                if discarding_unpaired:
                    total_unpaired_discarded[tch] += len(unpaired_to_add)
            iseqs_to_keep = [i for i in range(len(cline['unique_ids'])) if i not in iseqs_to_remove]
            process_unpaired(cline, ploci[tch], iseqs_to_keep, paired_iseqs)  # have to go back after finishing cluster since only now do we know who we ended up keeping
            if len(iseqs_to_keep) > 0:
                new_partition.append([cluster[i] for i in iseqs_to_keep])
                new_cline = utils.get_non_implicit_copy(cline)
                utils.restrict_to_iseqs(new_cline, iseqs_to_keep, glfos[ploci[tch]])  # note that this calls utils.add_implicit_info()
                new_antn_list.append(new_cline)
            if debug:
                def fstr(v): return '' if v==0 else '%d'%v
                nrstr = utils.color('red' if len(iseqs_to_keep)==0 else None, fstr(len(iseqs_to_remove)), width=3)
                print('    %s   %3d   %3s    %s     %3s   %3s    %3s      %s     %s' % (utils.locstr(ploci[tch]) if iclust==0 else ' ', len(cline['unique_ids']), fstr(len(iseqs_to_keep)), nrstr, fstr(n_no_info), fstr(n_other_light), fstr(n_non_reciprocal), utils.color('red', 'x') if len(unpaired_to_add)>0 and discarding_unpaired else ' ', ' '.join(cluster)))
        lp_cpaths[ploci[tch]] = ClusterPath(seed_unique_id=cpaths[ploci[tch]].seed_unique_id, partition=new_partition)
        lp_antn_lists[ploci[tch]] = new_antn_list

    if debug:
        print('    totals before: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in cpaths[ploci[tch]].best())) for tch in sorted(ploci)))
        print('    totals after: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in lp_cpaths[ploci[tch]].best())) for tch in sorted(ploci)))
    print('%skept %d total unpaired (i.e. setting aside in order to reintegrate after paired clustering): %s' % ('        ' if debug else '', sum(len(s) for s in unpaired_seqs.values()), '  '.join('%s %d'%(utils.locstr(l), len(unpaired_seqs[l])) for l in sorted(unpaired_seqs))))
    if keep_all_unpaired:
        print('            --keep-all-unpaired-seqs: keeping even unpaired seqs in families with no (well-)paired seqs')
    else:
        print('            note: discarding unpaired seqs in families with no (well-)paired seqs (i.e. families that are all either unpaired or paired to the other/wrong light chain), rather than keeping track of them to re-add later')
    if sum(total_unpaired_discarded.values()) > 0:
        print('        discarded %d total unpaired (in families with no (well-)paired seqs): %s' % (sum(total_unpaired_discarded.values()), '  '.join('%s %s'%(utils.locstr(ploci[c]), total_unpaired_discarded[c]) for c in sorted(ploci))))

    return lp_cpaths, lp_antn_lists, unpaired_seqs

# ----------------------------------------------------------------------------------------
# default: return True if <pids> is len 1 and is in the same family as <uid>'s correct partner
# near=True: same, but only has to be a family whose naive seq is within <max_hdist> hamming distance of the true family's naive seq (true_antn_dict must be set)
def crct_fam(true_partitions, uid, pids, near=False, max_hdist=3, true_antn_dict=None, tdbg=False):  # not mutually exclusive to the other categories (correct, mispaired, unpaired)
    # ----------------------------------------------------------------------------------------
    def getatn(utmp):
        ltmp = utmp.split('-')[-1]
        tclust = utils.get_single_entry([c for c in true_partitions[ltmp] if utmp in c])
        return true_antn_dict[ltmp][':'.join(tclust)]
    # ----------------------------------------------------------------------------------------
    if len(pids) != 1:
        if tdbg:
            print('    %s N pids %d' % (uid, len(pids)))
        return False
    if utils.is_correctly_paired(uid, pids[0]):
        if tdbg:
            print('    %s correct' % uid)
        return True
    plocus = pids[0].split('-')[-1]
    assert plocus in utils.loci
    # pntn = utils.get_single_entry([l for l in antn_dicts[plocus].values() if pids[0] in l['unique_ids']])  # paired annotation
    true_pclust = utils.get_single_entry([c for c in true_partitions[plocus] if pids[0] in c])  # cluster in true partition with paired id
    is_corr_fam = any(utils.is_correctly_paired(uid, p) for p in true_pclust)  # true if the correct paired id is in the paired cluster (i.e. we just got the wrong family member) NOTE it's weird to kind of guess if the pairing is correct, but this is the how we do it when we don't have the annotation (yes it would be better to just always have the annotation [which we only pass in if <near> is set)
    if near:  # NOTE this is pretty ugly and could really stand to be cleaned up
        if true_antn_dict is None:
            raise Exception('must set true_antn_dict if <near> is set')
        if is_corr_fam:
            return True
        inf_patn = getatn(pids[0])  # annotation for the inferred paired id
        true_uatn = getatn(uid)  # annotation for <uid>
        # true_pid = utils.get_single_entry(utils.per_seq_val(true_uatn, 'paired-uids', uid))
        # true_patn = getatn(true_pid)  # annotation for true paired id (NOTE this is quite different to inf_patn -- and that's fine, since our inference said that <uid> was associated with <pids[0]> (i.e. <uid> and the true partner of <pids[0]> are interchangeable), not that <true_pid> and <pids[0]> are interchangeable (the latter would mean that we had pair info for <true_pid>, whereas since <uid> has collided with something it is likely that its true partner was left unpaired)
        true_inf_id = utils.get_single_entry(utils.per_seq_val(inf_patn, 'paired-uids', pids[0]))  # true partner of inferred pid
        true_inf_atn = getatn(true_inf_id)
        if len(true_uatn['naive_seq']) == len(true_inf_atn['naive_seq']):
            hdist = utils.hamming_distance(true_uatn['naive_seq'], true_inf_atn['naive_seq'])  # naive hdist between the true naive seqs of <uid> and the true partner of <pids[0]>
        else:
            hdist = 9999
        return hdist <= max_hdist
    else:
        if tdbg:
            print('    %s corr fam: %6s   pid: %s   paired annotation uids: %s' % (uid, is_corr_fam, pids[0], ' '.join(true_pclust)))
        return is_corr_fam

# ----------------------------------------------------------------------------------------
def plot_fraction_correctly_paired(cpaths, antn_dicts, true_outfos, antn_lists=None, performance_outdir=None, plotdir=None, fnames=None, calc_near_fams=False):
    # ----------------------------------------------------------------------------------------
    def gpt(uid, pids):  # these are all mutually exclusive (as opposed to 'correct-family')
        if len(pids) == 0:
            return 'unpaired'
        elif len(pids) > 1:
            return 'multiple'
        elif utils.is_correctly_paired(uid, pids[0]):
            return 'correct'
        else:
            return 'mispaired'
    # ----------------------------------------------------------------------------------------
    if antn_dicts is None:
        assert antn_lists is not None
        antn_dicts = {l : utils.get_annotation_dict(antn_lists[l], cpath=cpaths[l]) for l in antn_lists}
    true_partitions = {l : true_outfos['merged']['cpaths'][l].best() for l in antn_dicts}
    true_antn_dict = {l : utils.get_annotation_dict(true_outfos['merged']['antn_lists'][l]) for l in antn_dicts} if calc_near_fams else None
    mcodes = ['correct', 'mispaired', 'unpaired', 'multiple']  # these are mutually exclusive
    kcodes = mcodes + ['correct-family', 'near-family', 'total']
    fcinfo, afo = {k : {} for k in kcodes}, {c : 0 for c in mcodes}
    for ltmp in sorted(cpaths):
        for cluster in cpaths[ltmp].best():
            atn = antn_dicts[ltmp][':'.join(cluster)]
            for uid, pids in zip(atn['unique_ids'], atn['paired-uids']):
                rcode = gpt(uid, pids)  # mutually exclusive result code
                fsize = len(utils.get_single_entry([c for c in true_partitions[ltmp] if uid in c]))
                if fsize not in fcinfo['total']:
                    for k in kcodes:  # add it to all of them so they all have the same bins
                        fcinfo[k][fsize] = 0
                if fsize not in fcinfo[rcode]:
                    fcinfo[rcode][fsize] = 0
                fcinfo['total'][fsize] += 1
                fcinfo[rcode][fsize] += 1
                afo[rcode] += 1
                if crct_fam(true_partitions, uid, pids):  # correct family
                    fcinfo['correct-family'][fsize] += 1
                if calc_near_fams and crct_fam(true_partitions, uid, pids, true_antn_dict=true_antn_dict, near=True):
                        fcinfo['near-family'][fsize] += 1
    fcinfo['all'] = afo
    for kcd in kcodes + ['all']:
        if kcd not in ['total', 'all']:
            for fsize in fcinfo[kcd]:
                fcinfo[kcd][fsize] = fcinfo[kcd][fsize] / float(fcinfo['total'][fsize])
        xbins = None if kcd=='all' else [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 30.5, 75.5, 100.5, 200.5, 500.5, 1000.5, 5000.5, 10000.5]
        fchist = hutils.make_hist_from_dict_of_counts(fcinfo[kcd], 'string' if kcd=='all' else 'int', 'pair cleaning performance', is_log_x=kcd!='all', arg_bins=xbins) #, no_sort=True)
        fcplname = 'true-pair-clean-performance%s' % ('' if kcd=='all' else '-'+kcd)
        if performance_outdir is not None:  # need to write it here so you have it even if you're not plotting
            fchist.write('%s/%s.csv'%(performance_outdir, fcplname))
        if plotdir is not None:
            fn = fchist.fullplot(plotdir, fcplname, pargs={'label' : kcd, 'ignore_overflows' : True}, fargs={'xticklabelsize' : 15, 'ylabel' : 'number of seqs'})  # , 'ybounds' : (0., 1.05)'xbounds' : (0.95, 1.05*len(fcinfo)), 
            fnames.append([fn])

# ----------------------------------------------------------------------------------------
def clean_pair_info(args, cpaths, antn_lists, plotdir=None, performance_outdir=None, max_hdist=4, true_outfos=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def check_droplet_id_groups(pid_groups, all_uids, tdbg=False):
        if not args.is_data:
            # print '  note: couldn\'t get droplet id from \'%s\', so assuming this isn\'t 10x data' % next(iter(all_uids))  # NOTE i'm not sure that this gives the same one as the previous line
            return False
        if len(set(len(g) for g in pid_groups)) == 1:
            print('    %s all pid groups are length 1 in check_droplet_id_groups(). Maybe you\'re missing pairing info?' % utils.color('yellow', 'warning'))
            return False
        # check against the droplet id method (we could just do it this way, but it would only work for 10x, and only until they change their naming convention)
        pgroup_strs = set(':'.join(sorted(pg)) for pg in pid_groups)  # <pid_groups>: list of pid groups, i.e. each element is the uids from a single droplet (for 10x), so <pgroup_strs> converts each group into an indexable key
        n_not_found, n_total = 0, 0
        if tdbg:
            print('              found?   drop id           contigs     overlaps (with any non-identical groups)')
        for dropid, drop_queries in utils.get_droplet_groups(all_uids, args.droplet_id_separators, args.droplet_id_indices):
            dqlist = list(drop_queries)
            found = ':'.join(sorted(dqlist)) in pgroup_strs  # was this exact combination of queries in pid_groups?
            n_total += 1
            if not found:  # if not, see if any pid groups have some of these queries
                n_not_found += 1
            if tdbg or not found:
                overlaps = [g for g in pgroup_strs if dropid in g]  # this should essentially always be length 1, except when we're missing pairing info in simulation (in which case i know we don't even want to be running this simulation, but whatever I'm testing edge cases)
                ostr = ' '.join(sorted(utils.get_contig_id(q, args.droplet_id_separators, args.droplet_id_indices) for q in overlaps[0].split(':'))) if len(overlaps)==1 else 'multiple'
                print('  %25s    %s   %-8s   %s' % (utils.color('green', '-') if found else utils.color('red', 'x'), dropid, ' '.join(sorted(utils.get_contig_id(q, args.droplet_id_separators, args.droplet_id_indices) for q in dqlist)),
                                                    utils.color('red', ostr if not found else '')))
        if n_not_found > 0:
            print('  %s droplet id group check failed for %d / %d groups, i.e. droplet ids parsed from uids don\'t match pair info: either pairing info is messed up or missing, or this is simulation and you didn\'t set --is-simu (if the latter, ignore this)' % (utils.color('red', 'error'), n_not_found, n_total))
        return True
    # ----------------------------------------------------------------------------------------
    def plot_uids_before(plotdir, pid_groups, all_antns):
        # ----------------------------------------------------------------------------------------
        def fnfplot(logstr, fhists, n_max_bins=15):
            from . import plotting
            fklabels = {'func' : 'all func.', 'nonfunc' : 'any non.'}
            fig, ax = plotting.mpl_init()
            for fk, fcolor in zip(fhists, plotting.default_colors):
                fhists[fk].mpl_plot(ax, label=fklabels[fk], color=fcolor, remove_empty_bins=True)
                if logstr == '':
                    fhists[fk].write('%s/%s.csv'%(plotdir, fk + '-per-drop'))
                xticks = fhists[fk].get_bin_centers()
                xticklabels = fhists[fk].bin_labels
                xbounds = None if fhists[fk].n_bins<n_max_bins else (0, n_max_bins)
            fn = plotting.mpl_finish(ax, plotdir, 'func-non-func-per-drop'+logstr, xlabel='N seqs per droplet', ylabel='counts', title='before', log='' if logstr=='' else 'y', leg_loc=(0.6, 0.7), xticks=xticks, xticklabels=xticklabels, xbounds=xbounds)
            fnames[0].append(fn)
        # ----------------------------------------------------------------------------------------
        bhist = Hist(value_list=[len(pg) for pg in pid_groups], init_int_bins=True)
        fn = bhist.fullplot(plotdir, 'seqs-per-droplet', pargs={'remove_empty_bins' : True}, fargs={'xlabel' : 'seqs per droplet', 'ylabel' : 'counts', 'title' : 'before', 'xticks' : [i for i in range(0, int(bhist.xmax+1), 2)]})
        fnames[0].append(fn)
        # fhists = {f : Hist(bhist.n_bins, bhist.xmin, bhist.xmax) for f in ['func', 'nonfunc']}
        flcounts = {f : {} for f in ['func', 'nonfunc']}
        per_seq_flcounts = {}
        for pgroup in pid_groups:
            funclist = [utils.is_functional(all_antns[p], all_antns[p]['unique_ids'].index(p)) for p in pgroup]
            # fhists['func'].fill(funclist.count(True))
            # fhists['nonfunc'].fill(funclist.count(False))
            fkey = 'nonfunc' if any(not f for f in funclist) else 'func'  # fill the nonfunc hist if *any* seq in this droplet is nonfunctional
            # lcounts = {l : lgstr(pgroup, for_plot=True).count(l) for l in utils.loci}
            lckey = lgstr(pgroup, for_plot=True)
            if lckey not in flcounts[fkey]:
                flcounts[fkey][lckey] = 0
            flcounts[fkey][lckey] += 1
            per_seq_flcounts.update({u : lckey for u in pgroup})

        # shenanigans to get bins sorted by the average fractional counts (summed over func and nonfunc individually, i.e. by the average fractional count, so bins will look out of order for func/nonfunc individually)
        ctotals = {fk : sum(flcounts[fk].values()) for fk in flcounts if len(flcounts[fk])>0}  # total counts for func and nonfunc
        flfracs = {}
        for fk in flcounts:
            for ck in flcounts[fk]:
                if ck not in flfracs:
                    flfracs[ck] = 0
                flfracs[ck] += flcounts[fk][ck] / (float(ctotals[fk]) * len(ctotals))
        assert utils.is_normed(flfracs)
        binlabels, _ = zip(*sorted(list(flfracs.items()), key=operator.itemgetter(1), reverse=True))

        fhists = {f : Hist(len(binlabels), -0.5, len(binlabels) - 0.5) for f in ['func', 'nonfunc']}
        for fstr in ['func', 'nonfunc']:
            for ibin, blabel in enumerate(binlabels):
                fhists[fstr].set_ibin(ibin + 1, flcounts[fstr].get(blabel, 0), math.sqrt(flcounts[fstr].get(blabel, 0)))
                fhists[fstr].bin_labels[ibin + 1] = blabel
        for logstr in ['', '-log']:
            fnfplot(logstr, fhists)
        return per_seq_flcounts
    # ----------------------------------------------------------------------------------------
    def plot_n_pseqs_per_seq(pstr):
        from . import plotconfig
        pidlengths = {}
        for ltmp in sorted(cpaths):
            for cluster in cpaths[ltmp].best():
                atn = antn_dicts[ltmp][':'.join(cluster)]
                pidlengths.update({u : len(set(pids) - set([u])) for u, pids in zip(atn['unique_ids'], atn['paired-uids'])})
        ahist = Hist(value_list=list(pidlengths.values()), init_int_bins=True)
        xvar = 'paired-seqs-per-seq'
        fn = ahist.fullplot(plotdir, '%s-%s'%(xvar, pstr), pargs={'remove_empty_bins' : True}, fargs={'xlabel' : plotconfig.xtitles.get(xvar), 'ylabel' : 'counts', 'title' : pstr, 'xbounds' : (-0.05, 1.05*ahist.xmax), 'xticks' : [i for i in range(0, int(ahist.xmax+1), 1 if 'after' in pstr else 2)]})
        fnames[1].append(fn)
        return pidlengths
    # ----------------------------------------------------------------------------------------
    def make_final_plots(initial_seqs_per_seq, initial_flcounts):
        final_seqs_per_seq = plot_n_pseqs_per_seq('after')
        from . import plotting
        fn = plotting.plot_smatrix(plotdir, 'pseq-matrix', xydicts=(final_seqs_per_seq, initial_seqs_per_seq), n_max_bins=12, xlabel='after', ylabel='before', lfcn=lambda x: 'miss.' if x==-1 else str(x), title='N paired seqs per seq')
        fnames[2].append(fn)
        final_flcounts = {}  # note that this has to be per seq (even though that kind of double counts) since otherwise we wouldn't have a way to determine correspondence between initial and final
        for ltmp in sorted(cpaths):
            for cluster in cpaths[ltmp].best():
                atn = antn_dicts[ltmp][':'.join(cluster)]
                final_flcounts.update({u : lgstr(set([u] + pids), for_plot=True) for u, pids in zip(atn['unique_ids'], atn['paired-uids'])})  # have to make sure <u> is included in <pids> (as well as that there's no duplicates)
        fn = plotting.plot_smatrix(plotdir, 'flcount-matrix', xydicts=(final_flcounts, initial_flcounts), kfcn=len, n_max_bins=15,
                                   lfcn=lambda x: 'miss.' if x==-1 else ('none' if x=='' else str(x)), xlabel='after', ylabel='before', title='pair combo (per seq)', tdbg=2 if debug else False)
        fnames[2].append(fn)
        plotting.make_html(plotdir, fnames=fnames)
    # ----------------------------------------------------------------------------------------
    def getloc(uid):
        return all_loci.get(uid, '?')
    # ----------------------------------------------------------------------------------------
    def gval(uid, key):  # get per-seq val for <uid>
        if uid not in all_antns:
            return None
        return utils.per_seq_val(all_antns[uid], key, uid)
    # ----------------------------------------------------------------------------------------
    def lgstr(lgroup, dont_sort=False, for_plot=False):  # return str representing locus of each uid in <lgroup>, e.g. [a, b, c] --> 'h  k  k' (or ['h', 'k', 'k'])
        sfcn = utils.pass_fcn if dont_sort else sorted
        lfcn = (lambda x: x[2]) if for_plot else utils.locstr
        lgstrs = [lfcn(l) for l in sfcn([getloc(u) for u in lgroup])]
        if for_plot and len(lgstrs) > 4:
            lcounts = [(l, lgstrs.count(l)) for l in sorted(set(lgstrs), key=lambda x: lgstrs.index(x))]  # make a set but then order it same as in the label
            lgstrs = ['%d%s'%(c, l) for l, c in lcounts]
        return ' '.join(lgstrs)
    # ----------------------------------------------------------------------------------------
    def choose_seqs_to_remove(chain_ids, remove_unproductive=False, tdbg=False):  # choose one of <chain_ids> to eliminate (based on identical/similar seq collapse and productivity)
        ids_to_remove = set(u for u in chain_ids if getloc(u)=='?')  # remove any with missing annotations
        if tdbg and len(ids_to_remove) > 0:  # i think this actually can't happen a.t.m.
            print('      removed %d with missing annotations' % len(ids_to_remove))

        # among any pairs of sequences that are [almost] identical at all non-ambiguous position, keep only the longest one (note that this is really preprocessing/error correction, so probably shouldn't really be here)
        if args.collapse_similar_paired_seqs:
            dbgstr = []
            n_equivalent = 0
            for tpair in itertools.combinations(chain_ids, 2):
                if len(set(getloc(u) for u in tpair)) > 1:  # can't be equivalent if they have different loci
                    continue
                if len(set(len(gval(u, 'seqs')) for u in tpair)) > 1:  # or if their (N-padded) sequences are different lengths
                    continue
                hdist = utils.hamming_distance(*[gval(u, 'seqs') for u in tpair])
                if tdbg:
                    dbgstr.append(utils.color('blue' if hdist==0 else 'yellow', '%d'%hdist))
                if hdist <= max_hdist:  # it would be nice to be able to combine their sequences (since they should have coverage only over different parts of vdj), but I think propagating the resulting annotation modifications would be hard
                    better_id, worse_id = sorted(tpair, key=lambda q: utils.ambig_frac(gval(q, 'seqs')))  # if we're tossing one with hdist > 0, it might make more sense to keep the lower-shm one if they're the same length, but I don't think it really matters (in the end we're still just guessing which is the right pairing)
                    ids_to_remove.add(worse_id)
                    n_equivalent += 1
            if tdbg and len(dbgstr) > 0:
                print('        %d pair%s equivalent with hdists %s' % (n_equivalent, utils.plural(n_equivalent), ' '.join(dbgstr)))

        # if specified, remove unproductive (only on real data, since simulation usually has lots of stop codons)
        if args.is_data and remove_unproductive:
            dbgstr = []
            unproductive_ids = []
            for uid in chain_ids:
                if not utils.is_functional(all_antns[uid], all_antns[uid]['unique_ids'].index(uid)):
                    unproductive_ids.append(uid)
                    if tdbg:
                        dbgstr.append(utils.is_functional_dbg_str(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), sep='+'))
            ids_to_remove |= set(unproductive_ids)
            if tdbg and len(unproductive_ids) > 0:
                print('        %d unproductive  %s' % (len(unproductive_ids), ',  '.join(dbgstr)))

        return ids_to_remove

    # ----------------------------------------------------------------------------------------
    def update_all_pid_info(finished_ids, tdbg=False):  # update pids in all annotations for all loci after finishing cleaning a cluster
        finished_ipgs = [pid_ids[u] for u in finished_ids]
        finished_pgroups = [pid_groups[i] for i in finished_ipgs]
        for pgroup in finished_pgroups:
            for uid in [u for u in pgroup if u not in finished_ids]:
                antn = all_antns[uid]
                iseq = antn['unique_ids'].index(uid)
                antn['paired-uids'][iseq] = [p for p in antn['paired-uids'][iseq] if p not in finished_ids]
    # ----------------------------------------------------------------------------------------
    def update_pid_info(cpids, tdbg=False):  # remove the two ids <cpids> from everybody's paired ids, since we've just decided they're properly paired (note that you do *not* want to modify <pfamilies> here -- that would be destroying the info you need to decided who's paired with who)
        if tdbg:
            print('  upd: %s' % ' '.join(cpids))
        assert len(cpids) == 2
        for cid, oid in zip(cpids, reversed(cpids)):
            cntn, ic = all_antns[cid], all_antns[cid]['unique_ids'].index(cid)
            if cntn['paired-uids'][ic] == [oid]:  # already ok
                continue
            if tdbg:
                print('      %s %3d --> %3d   %s' % (cid, len(cntn['paired-uids'][ic]), 1, ' '.join(utils.color('red' if p==oid else None, p) for p in cntn['paired-uids'][ic])))  # ok i know it should always just decrease by one, but maybe something else could break and you get duplicates?
            cntn['paired-uids'][ic] = [oid]
    # ----------------------------------------------------------------------------------------
    def ptn_clean(uloc, cline, cluster, remove_uncertain_pids=False):  # use information from the [clonal family] partitions to decide which of several potential paired uids is the correct one
        # ----------------------------------------------------------------------------------------
        def fidstr(fid):
            if fid == '':
                return fid
            else:
                if fid in name_ids:
                    fstr = name_ids[fid]
                else:
                    fstr, name_dict['potential'], name_dict['used'] = utils.choose_new_uid(name_dict['potential'], name_dict['used'])
                    name_ids[fid] = fstr
                return utils.color(utils.cyclecolor(fid, clist=['red', 'yellow', 'blue_bkg', 'reverse_video', 'green_bkg', 'purple', 'green', 'blue']), fstr)
        # ----------------------------------------------------------------------------------------
        def get_pfamily_dict(old_pfam_dict=None):  # see what others in its family are paired with
            pfdict = {}  # counts how many different uids in <cline> had paired ids from each potential paired family (well, how many pids total, but i think they'll always be the same)
            fid_counter = 0  # give a unique id to each family, for dbg printing purposes
            for uid, pids in zip(cline['unique_ids'], cline['paired-uids']):
                for pid in pids:
                    fline = all_antns[pid]  # family for this paired id
                    fkey = ':'.join(fline['unique_ids'])
                    if fkey not in pfdict:
                        pfdict[fkey] = {'locus' : gval(pid, 'loci'), 'count' : 0, 'id' : fid_counter if old_pfam_dict is None else old_pfam_dict[fkey]['id']}
                        fid_counter += 1
                    pfdict[fkey]['count'] += 1  # NOTE that if two cells from this family are in the same droplet, we can see the same pid here more than once (i.e. two counts here can be from the same uid)
            return pfdict
        # ----------------------------------------------------------------------------------------
        def pfkey(p):
            return ':'.join(all_antns[p]['unique_ids'])  # keystr for the family of paired id <p>
        # ----------------------------------------------------------------------------------------
        def print_dbg(old_pfams, new_pfams):
            # ----------------------------------------------------------------------------------------
            def prpfd(pfdict, extra_str=''):
                print('        %7s votes size  id  cdr3' % extra_str)
                for fkey, fdct in sorted(list(pfdict.items()), key=lambda x: x[1]['count'], reverse=True):
                    if utils.samechain(fdct['locus'], uloc):
                        continue
                    print('           %s    %3d   %3d  %2d %s  %3d' % (utils.locstr(fdct['locus']), fdct['count'], len(antn_dicts[fdct['locus']][fkey]['unique_ids']), fdct['id'], fidstr(fdct['id']), antn_dicts[fdct['locus']][fkey]['cdr3_length']))
            # ----------------------------------------------------------------------------------------
            def lcstr(pids, pfcs, pfids):  # returns string summarizing the families of the paired uids for a uid, e.g. 'k 51  l 3  h 1' if the uid has three potential pids, one from k with which 50 other uids in <cline> are paired, etc.
                if len(pids) == 0: return ''
                spids, spfcs, spfids = zip(*sorted(zip(pids, pfcs, pfids), key=operator.itemgetter(1), reverse=True))
                plocs = [getloc(u) for u in spids]
                return '  '.join('%s %s %s'%(utils.locstr(sl), '%2d'%sp, fidstr(fid)) for sl, sp, fid in zip(plocs, spfcs, spfids) if not utils.samechain(sl, uloc))
            # ----------------------------------------------------------------------------------------
            def truestr(uid, pids):
                if uid is None or args.is_data:
                    return ''
                is_correct = len(pids)==1 and utils.is_correctly_paired(uid, pids[0])
                cfam = crct_fam(true_partitions, uid, pids)
                return '%s ' % (utils.color('green', 'y') if is_correct else (utils.color('yellow', '*') if cfam else utils.color('red', 'n')))
            # ----------------------------------------------------------------------------------------
            # prpfd(old_pfams, extra_str='before:')
            # prpfd(new_pfams, extra_str='after:')
            new_pfcounts = [[new_pfams[pfkey(p)]['count'] for p in pids] for pids in cline['paired-uids']]
            new_pfids = [[new_pfams[pfkey(p)]['id'] for p in pids] for pids in cline['paired-uids']]
            new_estrs = [lcstr(pids, pfcs, pfd) for u, pids, pfcs, pfd in zip(cline['unique_ids'], cline['paired-uids'], new_pfcounts, new_pfids)]
            old_pfids = [[old_pfams[pfkey(p)]['id'] for p in pids] for pids in old_pids]
            old_pfcounts = [[old_pfams[pfkey(p)]['count'] for p in pids] for pids in old_pids]
            old_estrs = [lcstr(pids, pfcs, pfd) for pids, pfcs, pfd in zip(old_pids, old_pfcounts, old_pfids)]
            ulen, olen = [max(utils.len_excluding_colors(s) for s in slist) for slist in [cline['unique_ids'], old_estrs]]
            print('       %s%s     after       %s         pid(s)' % (utils.wfmt('uid', ulen, jfmt='-'), '' if args.is_data else 'correct', utils.wfmt('before', olen, jfmt='-')))
            for uid, pids, oldstr, newstr in zip(cline['unique_ids'], cline['paired-uids'], old_estrs, new_estrs):
                print('      %s    %s%s%s      %s       %s' % (utils.wfmt(uid, ulen), truestr(uid, pids), '' if args.is_data else '      ', utils.color('blue', '-' if len(pids)==0 else '', width=6) if len(pids)!=1 else newstr, oldstr + (olen - utils.len_excluding_colors(oldstr))*' ', ' '.join(pids)))
            print('')
        # ----------------------------------------------------------------------------------------
        if debug:
            name_dict, name_ids = {'potential' : None, 'used' : None}, {}  # infrastructure for color code strs
            old_pids = copy.deepcopy(cline['paired-uids'])
        old_pfams = get_pfamily_dict()  # map from each potential paired family to the number of uids in <cluster> that are potentially paired with it (i.e. the number of uids that are voting for it)

        # for each uid, choose the pid that's of opposite chain, and has the most other uids voting for it (as long as some criteria are met)
        finished_ids = []
        for iseq, uid in enumerate(cline['unique_ids']):
            cline['paired-uids'][iseq] = [p for p in cline['paired-uids'][iseq] if p not in finished_ids and not utils.samechain(getloc(p), getloc(uid))]
            pidfcs = [(p, old_pfams[pfkey(p)]['count']) for p in cline['paired-uids'][iseq]]  # (pid, pcount) for all opposite-chain pids, where <pcount> is the number of votes for <pid>'s family (note that the 'paired-uids' get modified as we go through the loop)
            pid_to_keep = None
            if len(pidfcs) > 0:
                sorted_pids, sorted_pfcs = zip(*sorted(pidfcs, key=operator.itemgetter(1), reverse=True))
                # note that even if there's only one ochain choice, there can be other same-chain ones that we still want to drop (hence the <2 below)
                if len(sorted_pfcs) < 2 or sorted_pfcs[0] > sorted_pfcs[1] or pfkey(sorted_pids[0]) == pfkey(sorted_pids[1]):  # in order to drop the later ones, the first one either has to have more counts, or at least the second one has to be from the same family (in the latter case we still don't know which is the right one, but for the purposes of clustering resolution we just need to know what the family is)
                    pid_to_keep = sorted_pids[0]
                    update_pid_info([uid, pid_to_keep])
                    finished_ids += [uid, pid_to_keep]
            if pid_to_keep is None and remove_uncertain_pids:
                cline['paired-uids'][iseq] = []  # if we didn't decide on an opposite-chain pid, remove all pairing info

        update_all_pid_info(finished_ids)  # have to update everybody before recalculating pfamily dicts

        if debug: # and len(cline['unique_ids']) > 1:  # NOTE it's annoying printing the singletons, but it's way worse when they're just missing and you can't figure out where a sequence went
            print_dbg(old_pfams, get_pfamily_dict(old_pfam_dict=old_pfams))
    # ----------------------------------------------------------------------------------------
    def clean_with_partition_info():  # use information from the [clonal family] partitions to decide which of several potential paired uids is the correct one
        # orig_pids, orig_pfams = None, None
        # if debug:
        #     orig_pids = {l : {':'.join(c) : copy.deepcopy(antn_dicts[l][':'.join(c)]['paired-uids']) for c in cpaths[l].best()} for l in cpaths}  # keep track of original pids
        #     orig_pfams = get_pfamily_dict()
        for ltmp in sorted(cpaths):
            if debug:
                print('  %s starting pair cleaning with partition info' % utils.color('green', ltmp))
            for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
                ptn_clean(ltmp, antn_dicts[ltmp][':'.join(cluster)], cluster, remove_uncertain_pids=True)
        # NOTE i would think it would help to do this twice, only removing uncertain ones the second time (or, maybe better [faster], loop over uids in order of confidence that we'll get them correct)
        # for ltmp in sorted(cpaths):
        #     if debug:
        #         print '%s starting pair cleaning (second round)' % utils.color('green', ltmp)
        #     for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
        #         ptn_clean(ltmp, antn_dicts[ltmp][':'.join(cluster)], cluster, remove_uncertain_pids=True)
    # ----------------------------------------------------------------------------------------
    if true_outfos is not None:
        true_partitions = {l : true_outfos['merged']['cpaths'][l].best() for l in antn_lists}
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l], cpath=cpaths[l]) for l in antn_lists}
    all_uids = set(u for p in cpaths.values() for c in p.best() for u in c)  # all uids that occur in a partition (should I think be the same as the ones for which we have valid/non-failed annotations)
    print('cleaning pair info for %d seqs' % len(all_uids))
    sys.stdout.flush()

    # first collect some information for later use
    pid_groups = []  # list of pid groups, i.e. each element is the uids from a single droplet (for 10x)
    pid_ids = {}  # map from each uid to the index of its pid group
    all_antns = {}  # map from each individual uid to its annotation
    all_loci = {u : l for l, ants in antn_lists.items() for antn in ants for u in antn['unique_ids']}  # for speed
    n_missing = 0
    if debug:
        print('  %s consolidating info for %d loci with family/sequence counts: %s' % (utils.color('blue', '+'.join(sorted(cpaths))), len(cpaths), '  '.join('%s: %d/%d'%(l, len(cpaths[l].best()), sum(len(c) for c in cpaths[l].best())) for l in sorted(cpaths))))
    for ltmp in sorted(cpaths):
        for cluster in cpaths[ltmp].best():
            if ':'.join(cluster) not in antn_dicts[ltmp]:
                print('  %s no annotation for %s' % (utils.color('yellow', 'warning'), ':'.join(cluster)))
                continue
            cline = antn_dicts[ltmp][':'.join(cluster)]
            if 'paired-uids' not in cline:
                print('  %s no paired uids in line' % utils.color('yellow', 'warning'))
                continue  # maybe should still add to all_antns?
            for uid, pids in zip(cline['unique_ids'], cline['paired-uids']):
                if pids is None:
                    raise Exception('None type \'paired-uids\' for %s' % uid)
                missing_ids = set(pids) - all_uids
                n_missing += len(missing_ids)
                pset = set([uid] + pids) - missing_ids
                found = False
                for pd in pids:  # even if we assume the pair info from all seqs is consistent (which this more or less does), we need to do this to get the <ipg>
                    if pd in pid_ids:
                        ipg = pid_ids[pd]
                        found = True
                        pid_groups[ipg] |= pset  # don't really need this, but i guess it catches one possible way in which pair info could be inconsistent? Although we should really probably crash if so
                        break
                if not found:
                    pid_groups.append(pset)
                    ipg = len(pid_groups) - 1
                assert ipg is not None
                for pid in pset:
                    pid_ids[pid] = ipg

            cline['loci'] = [ltmp for _ in cline['unique_ids']]  # could maybe add this somewhere else, like in partitiondriver? (eh, maybe not? the locus is always available in each file from the germline info anyway)
            for uid in cline['unique_ids']:
                all_antns[uid] = cline
    if n_missing > 0:
        print('   %d/%d (%.2f) missing uids when cleaning pair info' % (n_missing, len(all_uids), n_missing / float(len(all_uids))))  # NOTE at least for now we're skipping invalid queries when reading output
    for ipg, pg in enumerate(pid_groups):
        #     print '  %3d %s' % (ipg, ' '.join(pg))
        pid_groups[ipg] = sorted(pg)  # need to sort for replicability

    # idg_ok = check_droplet_id_groups(pid_groups, all_uids)  # NOTE not using the return value here, but I may need to in the future UPDATE turning it off, it's too hard to figure out when we want it on vs not
    fnames = None
    if plotdir is not None:
        fnames = [[], [], []]
        initial_flcounts = plot_uids_before(plotdir, pid_groups, all_antns)
        initial_seqs_per_seq = plot_n_pseqs_per_seq('before')

    # then go through each group trying to remove as many crappy/suspicously similar ones as possible (this step has a fairly minor effect compared to the partition-based step below)
    if debug:
        print('  cleaning %d pid groups (by removing crappy/too similar seqs):' % len(pid_groups))
        def tmpincr(pgroup, cdict):
            if lgstr(pgroup) not in cdict:
                cdict[lgstr(pgroup)] = 0
            cdict[lgstr(pgroup)] += 1
    ok_groups, tried_to_fix_groups, id_removed_groups = {}, {}, {}
    for ipg, pgroup in enumerate(pid_groups):
        pgroup = [u for u in pgroup if getloc(u) != '?']  # maybe need to figure out something different to do with missing ones?
        hids = [u for u in pgroup if utils.has_d_gene(getloc(u))]
        lids = [u for u in pgroup if u not in hids]
        if len(hids) < 2 and len(lids) < 2:  # all good, nothing to do
            pid_groups[ipg] = pgroup
            if debug: tmpincr(pgroup, ok_groups)
            continue
        if debug > 1:
            print('    %s' % lgstr(pgroup), end=' ')
        for chain, idlist in zip(utils.chains, [hids, lids]):
            if len(idlist) < 2:  # skip whichever of the chains has only one id
                continue
            if debug > 1:
                print('\n      too many %s chains: %s' % (chain, lgstr(idlist)))
            ids_to_remove = choose_seqs_to_remove(idlist, tdbg=debug>1)
            if debug and len(ids_to_remove) > 0:
                tmpincr(pgroup, id_removed_groups)
            for rid in ids_to_remove:
                pgroup.remove(rid)
                idlist.remove(rid)
                pid_groups.append(set([rid]))  # add the removed id to a new pid group of its own at the end (so it'll show up as unpaired)
                pid_ids[rid] = len(pid_groups) - 1
            if debug > 1:
                print('      %s: removed %d, leaving %d%s' % (utils.color('green', 'fixed') if len(idlist)==1 else utils.color('red', 'still too many'), len(ids_to_remove), len(idlist), ':' if len(idlist)>1 else ''))
                if len(idlist) > 1:
                    for uid in idlist:
                        prutils.print_seq_in_reco_event(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), one_line=True, extra_str='        ', uid_extra_str=utils.locstr(getloc(uid)))

        pid_groups[ipg] = pgroup
        if debug: tmpincr(pgroup, tried_to_fix_groups)

    if debug:
        def prcgrps(cdict, prestr):
            print('    %s%s' % (prestr, '' if len(cdict)>0 else utils.color('blue', ' none')))
            for lstr, count in sorted(list(cdict.items()), key=operator.itemgetter(1), reverse=True):
                print('      %3d  %s' % (count, lstr))
        prcgrps(ok_groups, 'ok (less than 2 of both h and l ids) to start with:')
        prcgrps(id_removed_groups, 'removed ids from:')
        prcgrps(tried_to_fix_groups, 'after trying to fix:')

    # re-set 'paired-uids' key in <cline> to the pid group (which was modified in the previous cleaning steps)
    for ltmp in sorted(cpaths):
        for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
            cline = antn_dicts[ltmp][':'.join(cluster)]  # if this kicks a key error, the annotation probably failed at some weird point (e.g. after sw, maybe in hmm)
            if any(u not in pid_groups[pid_ids[u]] for u in cline['unique_ids']):  # shouldn't be able to happen any more, but that was a really bad/dumb bug
                raise Exception('one of unique ids %s not in its own pid group' % cline['unique_ids'])
            cline['paired-uids'] = [[p for p in pid_groups[pid_ids[u]] if p != u] for u in cline['unique_ids']]

    # then go through again using cluster/family information to try to cut everyone down to one paired id (and if we can't get down to one, we remove all of them) NOTE also re-sets the actual 'paired-uids' keys
    clean_with_partition_info()

    # partially synchronize info from different loci: for any uid with a single pid for which that pid has no pair info, set the latter to [uid]
    n_fixed = {l : 0 for l in cpaths}
    for ltmp in sorted(cpaths):
        for iclust, cluster in enumerate(cpaths[ltmp].best()):
            cline = antn_dicts[ltmp][':'.join(cluster)]
            for iseq, uid in enumerate(cline['unique_ids']):
                pids = cline['paired-uids'][iseq]
                if len(pids) != 1:  # only look at cases where <uid> has a single unique pid
                    continue
                pantn = all_antns[pids[0]]
                ipid = pantn['unique_ids'].index(pids[0])
                if len(pantn['paired-uids'][ipid]) == 0:
                    pantn['paired-uids'][ipid] = [uid]
                    n_fixed[getloc(pids[0])] += 1
    if len(n_fixed) > 0:
        print('     synchronized/fixed %d pairs where one had no pair info after cleaning: %s' % (sum(n for n in n_fixed.values()), '  '.join('%s %d'%(utils.locstr(l), n_fixed[l]) for l in sorted(n_fixed))))

    if not args.is_data and (performance_outdir is not None or plotdir is not None):
        plot_fraction_correctly_paired(cpaths, antn_dicts, true_outfos, performance_outdir=performance_outdir, plotdir=plotdir, fnames=fnames)
    if plotdir is not None:
        make_final_plots(initial_seqs_per_seq, initial_flcounts)

# ----------------------------------------------------------------------------------------
def compare_partition_pair(cfpart, refpart, remove_from_ref=False, add_all_to_ref=False, antn_list=None, unpids=None, seed_unique_id=None, dbg_str=None, cf_label='inferred', ref_label='true', fail_frac=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def incorporate_duplicates(tpart):  # take the map from uid to list of its duplicates (dup_dict), and add the duplicates to any clusters in partition tpart that contain that uid
        for tclust in tpart:
            for uid in tclust:
                if uid in dup_dict:
                    tclust += dup_dict[uid]
    # ----------------------------------------------------------------------------------------
    dup_dict = {}
    if antn_list is not None:
        dup_dict = {u : l['duplicates'][i] for l in antn_list for i, u in enumerate(l['unique_ids']) if len(l['duplicates'][i]) > 0}  # now that I've turned off default duplicate removal, this won't happen much any more
        if len(dup_dict) > 0:
            print('%s duplicate sequences in joint partition evaluation' % utils.color('yellow', 'warning'))

    cfpart = copy.deepcopy(cfpart)  # don't want to modify opject pointed to by <cfpart)
    if len(dup_dict) > 0:
        incorporate_duplicates(cfpart)
    # these remove/add shenanigans are mostly necessary for heavy chain on the h+k/h+l (not merged) partitions just because h seqs will get put with the wrong light chain
    if remove_from_ref:
        refpart = utils.remove_missing_uids_from_partition(refpart, cfpart, ref_label=ref_label, miss_label=cf_label, fail_frac=fail_frac, warn=True, dbgstr=' (probably either seqs that failed inference [if light chain] or mispaired seqs [if heavy chain])')  # returns a new/copied partition, doesn't modify original
    if unpids is not None:
        refpart = utils.add_missing_uids_to_partition(refpart, cfpart, miss_label=ref_label, ref_label=cf_label, allowed_uids=unpids, warn=True, dbgstr=' (restricting to unpaired_uids, i.e. probably these were unpaired seqs that were re-added to the inferred partition, but in the true partition were paired with the other chain)')  # NOTE order of args (ref/miss) is reversed here
    if add_all_to_ref:
        refpart = utils.add_missing_uids_to_partition(refpart, cfpart, miss_label=ref_label, ref_label=cf_label, fail_frac=fail_frac, warn=True, dbgstr=' (giving up and adding whatever is needed, which is sometimes/maybe just from mispaired seqs)')  # NOTE order of args (ref/miss) is reversed here
    return utils.per_seq_correct_cluster_fractions(cfpart, refpart, seed_unique_id=seed_unique_id, dbg_str=dbg_str, inf_label=cf_label, true_label=ref_label, debug=debug)

# ----------------------------------------------------------------------------------------
def evaluate_joint_partitions(ploci, true_outfos, init_partitions, joint_partitions, antn_lists, seed_unique_ids=None, unpaired_seqs=None, fail_frac=None, debug=False):
    if seed_unique_ids is not None and set(ploci.values()) != set(seed_unique_ids):  # skip non-seed light locus pair
        return
    # NOTE that <joint_partitions> can have many fewer seqs than <init_partitions> since in making <joint_partitions> we remove seqs paired with the other light chain (the weighted average ccfs over the h joint partitions corresponding to both light chains would be exactly comparable to the <init_partitions>, but I think this is fine as it is)
    lpair = tuple(ploci[c] for c in 'hl')  # ick
    ccfs = {}
    for tch in utils.chains:
        ltmp = ploci[tch]
        lcfs = {}
        for dstr, cfpart, true_cpaths in [('single', init_partitions[tch], true_outfos['merged']), ('joint', joint_partitions[tch], true_outfos['lpairs'][lpair])]:
            lcfs[dstr] = compare_partition_pair(cfpart, true_cpaths['cpaths'][ploci[tch]].best(), seed_unique_id=None if seed_unique_ids is None else seed_unique_ids[ltmp], remove_from_ref=True,  add_all_to_ref=tch=='h',  # removes from the true ptn any uids that are missing from the inferred ptn
                                                antn_list=antn_lists[ltmp], unpids=unpaired_seqs[ltmp] if tch=='h' else None, dbg_str='%s %s '%(utils.locstr(ltmp), dstr), fail_frac=None if tch=='h' else fail_frac, debug=debug)  # only apply fail_frac to light chain, since we can't guarantee we put the right h chain seqs to k vs l if pair info is super fucked (should apply fail_frac to heavy chain when we calculate on merged partition)
        ccfs[ltmp] = lcfs
        print('  %s ccfs:     purity  completeness' % utils.locstr(ltmp))
        print('      single  %6.3f %6.3f' % (lcfs['single'][0], lcfs['single'][1]))
        print('       joint  %6.3f %6.3f' % (lcfs['joint'][0], lcfs['joint'][1]))
    return ccfs

# ----------------------------------------------------------------------------------------
# cartoon explaining algorithm here https://github.com/psathyrella/partis/commit/ede140d76ff47383e0478c25fae8a9a9fa129afa#commitcomment-40981229
def merge_chains(ploci, cpaths, antn_lists, unpaired_seqs=None, iparts=None, check_partitions=False, true_outfos=None, input_cpaths=None, input_antn_lists=None, seed_unique_ids=None, overmerge=False, naive_hamming_bound_type=None, fail_frac=None, debug=False):  # NOTE the clusters in the resulting partition generally have the uids in a totally different order to in either of the original partitions
    dbgids = None #['1437084736471665213-igh']  # None
    # ----------------------------------------------------------------------------------------
    def akey(klist):
        return ':'.join(klist)
    # ----------------------------------------------------------------------------------------
    def any_in_common(l1, l2):  # true if any uids in any cluster in l1 are found in any clusters in l2 NOTE kind of dumb but atm only calling this with <l1> of len 1
        # l2ids = set(itertools.chain.from_iterable(l2))
        # return any(len(set(c) & l2ids) > 0 for c in l1)
        l1ids, l2ids = [set(itertools.chain.from_iterable(l)) for l in [l1, l2]]
        return len(l1ids & l2ids) > 0
    # ----------------------------------------------------------------------------------------
    def common_clusters(tclust, tlist, return_indices=False):  # return all clusters in tlist that have uids in common with tclust
        tset = set(tclust)
        return [(i if return_indices else c) for i, c in enumerate(tlist) if len(set(c) & tset) > 0]
    # ----------------------------------------------------------------------------------------
    def is_clean_partition(putative_partition):  # make sure the list of clusters is actually disjoint
        return not any(len(set(c1) & set(c2)) > 0 for c1, c2 in itertools.combinations(putative_partition, 2))
    # ----------------------------------------------------------------------------------------
    # Starting with <single_cluster> (from one chain) and <cluster_list> (all clusters in the other chain that overlap with <single_cluster>), decide which of the "splits" (i.e. cluster boundaries) in <cluster_list> should be applied to <single_cluster>.
    # Reapportions all uids from <single_cluster> and <cluster_list> into <return_clusts>, splitting definitely/first by cdr3, and then (if over some threshold) by naive hamming distance.
    def resolve_discordant_clusters(single_cluster, single_annotation, cluster_list, annotation_list, tdbg=False):
        if dbgids is not None:
            tdbg = len(set(dbgids) & set(u for c in [single_cluster] + cluster_list for u in c)) > 0
        # NOTE single_cluster and cluster_list in general have quite different sets of uids, and that's fine. All that matters here is we're trying to find all the clusters that should be split from one another (without doing some all against all horror)
        if len(cluster_list) < 2:  # nothing to do
            return [single_cluster]  # NOTE <single_cluster> doesn't get used after here
        adict = utils.get_annotation_dict(annotation_list)
        cdr3_groups = utils.group_seqs_by_value(cluster_list, lambda c: adict[akey(c)]['cdr3_length'])  # group the together clusters in <cluster_list> that have the same cdr3 (there's already utils.split_clusters_by_cdr3(), but it uses different inputs (e.g. sw_info) so i think it makes sense to not use it here)
        if tdbg:
            print('   %s one cluster size: %3d  %s' % (utils.color('blue', 'syncing'), len(single_cluster), ':'.join(single_cluster)))
            jstrs = ['           %s %3d  %s' % ('vs %2d with sizes:'%len(cluster_list) if i==0 else '                 ', len(c), ':'.join(c)) for i, c in enumerate(cluster_list)]
            print('\n'.join(jstrs))
            print('     split into %d cdr3 group%s' % (len(cdr3_groups), utils.plural(len(cdr3_groups))))
        if not overmerge:  # default/normal
            assert naive_hamming_bound_type is not None  # ugly, but i want to only have the default set in bin/partis
            _, hi_hbound = utils.get_naive_hamming_bounds(naive_hamming_bound_type, overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))
        else:  # don't do any naive hamming splitting (correct only for --n-final-clusters 1)
            hi_hbound = 1.
        return_clusts = []
        for icdr, cdrgroup in enumerate(cdr3_groups):  # within each cdr3 group, split (i.e. use the cluster boundaries from cluster_list rather than single_cluster) if naive hfrac is > hi_hbound (but then there's shenanigans to adjudicate between different possibilities)
            if tdbg: print('      %s' % utils.color('purple', 'icdr %d' % icdr))

            # first figure out who needs to be split from whom
            clusters_to_split = {akey(c) : [] for c in cdrgroup}  # map from each cluster ('s key) to a list of clusters from which it should be split
            for c1, c2 in itertools.combinations(cdrgroup, 2):  # we could take account of the hfrac of both chains at this point, but looking at only the "split" one rather than the "merged" one, as we do here, is i think equivalent to assuming the merged one has zero hfrac, which is probably fine, since we only split if the split chain is very strongly suggesting we split
                hfrac = utils.hamming_fraction(adict[akey(c1)]['naive_seq'], adict[akey(c2)]['naive_seq'], align_if_necessary=True)  # all clusters with the same cdr3 len have been padded in waterer so their naive seqs are the same length
                if hfrac > hi_hbound:
                    clusters_to_split[akey(c1)].append(c2)
                    clusters_to_split[akey(c2)].append(c1)
                    if tdbg: print('         hfrac split %.3f > %.3f  %3d %3d  %s   %s' % (hfrac, hi_hbound, len(c1), len(c2), ':'.join(c1), ':'.join(c2)))

            # then do the splitting, which is accomplished by merging each cluster in <cdrgroup> with every other cluster in <cdrgroup> from which we aren't supposed to split it (i.e. that aren't in its <clusters_to_split>)
            if tdbg:
                print('                  N to     new')
                print('          size    split   cluster?')
            tmpclusts_for_return = []  # final (return) clusters for this cdr3 class
            for cclust in cdrgroup:
                split_clusts = clusters_to_split[akey(cclust)]
                if tdbg: print('         %4d    %3d' % (len(cclust), len(split_clusts)), end=' ')
                found_one = False
                for rclust in tmpclusts_for_return:  # look for an existing return cluster to which we can merge cclust, i.e. that doesn't have any uids from which we want to split
                    if any_in_common([rclust], split_clusts):  # if any uid in rclust is in a cluster from which we want to be split, skip it, i.e. don't merge with that cluster (note that we have to do it by uid because the rclusts are already merged so don't necessarily correspond to any existing cluster)
                        continue
                    # if found_one: print 'it happened!'  # can't happen any more since I switched to 'break' (although see note below)
                    if tdbg: print('     merge with size %3d   %s' % (len(rclust), ':'.join(cclust)))
                    rclust += cclust
                    found_one = True
                    break  # i.e. we just merge with the first one we find and stop looking; if there's more than one, it means we could merge all three together if we wanted (triangle inequality-ish, see diagram linked at top of fcn), but i doubt it'll matter either way, and this is easier
                if not found_one:
                    if tdbg: print('      y                    %s' % ':'.join(cclust))
                    tmpclusts_for_return.append(cclust)  # if we didn't find an existing cluster that we can add it to, add it as a new cluster

            return_clusts += tmpclusts_for_return

        if tdbg:
            print('      resolved clusters:')
            for tclust in return_clusts:
                print('          %s' % ':'.join(tclust))
        return return_clusts
    # ----------------------------------------------------------------------------------------
    def incorporate_rclusts(final_partition, fclust_sets, fclust_indices, resolved_clusters):  # incorporate <resolved_clusters> into <final_partition>
        if debug:
            dbgheader = ['    adding %d resolved cluster%s to %d clusters in final partition' % (len(resolved_clusters), utils.plural(len(resolved_clusters)), len(final_partition)), '      ifclust  fset   rset   common  after: fset  rset',]
            # ----------------------------------------------------------------------------------------
            def appdbg(new_fset, rset, common_uids, xdbg=None):
                def cstr(tclust): return '(empty)' if len(tclust)==0 else ':'.join(utils.color('light_blue_bkg' if u in common_uids else None, u) for u in sorted(tclust))
                dbgstr.append(' %3d %s%3d %s  %3d    %s   %s' % (len(new_fset | common_uids), ('%-s'%utils.color('red', str(len(new_fset)), width=3, padside='right')) if len(common_uids&new_fset)==0 else '   ',
                                                                 len(rset | common_uids), ('%-s'%utils.color('red', str(len(rset)), width=3, padside='right')) if len(common_uids&rset)==0 else '   ',
                                                                 len(common_uids), cstr(new_fset), cstr(rset)))
                if xdbg is not None:
                    dbgstr.append(xdbg)
            # ----------------------------------------------------------------------------------------
            def prdbg(ifclust, dbgheader, dbgstr):
                if len(dbgheader) > 0:
                    print('\n'.join(dbgheader))
                    dbgheader = []
                print('       %4d  %s' % (ifclust, '\n             '.join(dbgstr)))
                return dbgheader
        # ----------------------------------------------------------------------------------------
        rc_sets, rc_ids = [set(c) for c in resolved_clusters], set(itertools.chain.from_iterable(resolved_clusters))  # just for speed (note that rc_ids doesn't get updated, which i think is ok)
        def get_overlap_indices():  # separate fcn just to make profiling easier (should maybe remove fcn at some point)
            return sorted(set(fclust_indices.get(u) for u in rc_ids) - set([None]))  # indices of final clusters that overlap with <resolved_clusters> (None thing is weird, but i want to avoid having to both see if the uid is in there, and then access it) (sort is just so results are identical to previous code that worked differently, can eventually be removed)
        ifovrlps = get_overlap_indices()
        n_clean = len(fclust_sets) - len(ifovrlps)
        # for each cluster that's already in <final_partition> that has uids in common with a cluster in <resolved_clusters>, decide how to apportion the common uids (basically we remove them from the larger of the two clusters)
        for ifclust in ifovrlps:  # iteration/<ifclust> won't get as far as any clusters that we're just adding (to the end of <final_partition>), which is what we want
            fclust, old_fset = final_partition[ifclust], fclust_sets[ifclust]
            def get_rc_indices():  # separate fcn just to make profiling easier (should maybe remove fcn at some point)
                return [i for i, c in enumerate(rc_sets) if len(c & old_fset) > 0]  # indices of any resolved_clusters that overlap with this fclust
            irclusts = get_rc_indices()
            if debug: dbgstr = []
            new_fset = set(fclust)  # we'll remove uids from this, and then replace fclust with its remains
            for irclust in irclusts:  # resolve any discrepancies between these newly-resolved clusters and fclust
                rset = set(resolved_clusters[irclust])
                common_uids = new_fset & rset
                if irclusts.index(irclust) == 0:
                    if len(new_fset) > len(rset):  # remove the common ids from the larger one (effectively splitting according to the splittier one)
                        new_fset -= common_uids
                    else:
                        rset -= common_uids
                    if debug: xdbg = None
                else:  # if this is the second or greater rclust, we need to apply the splits between rclusts (otherwise we can end up merging uids from different rclusts)
                    new_fset -= common_uids
                    rset -= common_uids
                    resolved_clusters.append(sorted(common_uids))  # this adds a cluster at the end, which of course gets ignored in this loop over irclusts, but will get considered in the next fclust
                    rc_sets.append(common_uids)
                    if debug: xdbg = '                  %s  %s' % (utils.color('red', '+%-3d'%len(resolved_clusters[-1])), ':'.join(resolved_clusters[-1]))
                resolved_clusters[irclust] = sorted(rset)  # replace this resolved cluster with a copy of itself that may have had any common uids removed (if it was bigger than fclust)
                rc_sets[irclust] = rset
                if debug: appdbg(new_fset, rset, common_uids, xdbg)
            final_partition[ifclust] = sorted(new_fset)  # replace <fclust> (even if nothing was removed, which shuffles the order of unchanged clusters, but oh well)
            fclust_sets[ifclust] = set(final_partition[ifclust])
            for fid in final_partition[ifclust]:
                fclust_indices[fid] = ifclust
            if debug: dbgheader = prdbg(ifclust, dbgheader, dbgstr)
        if check_partitions:
            assert is_clean_partition(resolved_clusters)
        final_partition += resolved_clusters  # add the (potentially modified) resolved clusters
        fclust_sets += [set(c) for c in resolved_clusters]
        for iftmp in range(len(final_partition) - len(resolved_clusters), len(final_partition)):
            for fid in final_partition[iftmp]:
                fclust_indices[fid] = iftmp
    # ----------------------------------------------------------------------------------------
    def re_add_unpaired(joint_partitions, unpaired_seqs, tdbg=False):
        if huge_dbg or tdbg:
            print('      re-adding %d unpaired' % sum(len(lseqs) for lseqs in unpaired_seqs.values()))
        if huge_dbg:
            print('        finished:', end=' ')
            sys.stdout.flush()
            ihuge = 0
        n_added = {tch : {'singleton' : 0, 'new-cluster' : 0, 'existing-cluster' : 0} for tch in ploci}
        for tch, ltmp in ploci.items():
            jp_indices = {u : i for i, c in enumerate(joint_partitions[tch]) for u in c}  # just for speed
            for upid, nearfo in unpaired_seqs[ltmp].items():  # <upid> is uid of seq with bad/no pair info, <nearfo['nearest']> is uid of nearest seq in <upid>'s original family
                if huge_dbg:
                    if ihuge % (n_huge_dbg // 5) == 0:
                        print(' %d' % ihuge, end=' ')
                        sys.stdout.flush()
                    ihuge += 1
                if nearfo['nearest'] is None and len(nearfo['single-chain-family']) == 1:  # first clause: no paired seqs in its family, second clause: actually a singleton
                    joint_partitions[tch].append([upid])
                    jp_indices[upid] = len(joint_partitions[tch]) - 1
                    n_added[tch]['singleton'] += 1
                    continue
                nearids = set([nearfo['nearest-paired']] if nearfo['nearest-paired'] is not None else nearfo['single-chain-family'])  # if there's any paired seqs in its single chain family, attach it to the nearest one of those; otherwise try to keep all the unpaired seqs from the family together (note that the old method, of always attaching to the 'nearest' id whether it was paired or not, had the effect of splitting in some cases, which we don't want)
                ijclusts = sorted(set(jp_indices.get(u) for u in nearids) - set([None]))  # sort is just so results are identical to previous code that worked differently, can eventually be removed
                if len(ijclusts) < 1:  # it didn't have a 'nearest-paired' (i.e. no paired seqs in its single chain cluster), and we haven't gotten to any of the other unpaired seqs from its single chain cluster (when we do get to them, they'll get added to this cluster)
                    joint_partitions[tch].append([upid])
                    jp_indices[upid] = len(joint_partitions[tch]) - 1
                    n_added[tch]['new-cluster'] += 1
                    continue
                if len(ijclusts) > 1:
                    print('  %s multiple jclusts for %s: %s' % (utils.wrnstr(), upid, ijclusts))
                joint_partitions[tch][ijclusts[0]].append(upid)
                jp_indices[upid] = ijclusts[0]
                n_added[tch]['existing-cluster'] += 1
            if tdbg:
                print('          %s joint clusters with re-added unpaired seqs in %s:' % (utils.locstr(ltmp), utils.color('red', 'red')))
                for ijoint, uipairs in itertools.groupby(sorted(jp_indices.items(), key=operator.itemgetter(1)), key=operator.itemgetter(1)):
                    print('                %s' % ' '.join(utils.color('red' if u in unpaired_seqs[ltmp] else None, u) for u in joint_partitions[tch][ijoint]))
        if huge_dbg:
            print('')
        totstr = '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in joint_partitions[tch])) for tch in sorted(ploci))
        print('    re-added unpaired seqs (%s) to give total seqs in joint partitions: %s' % (', '.join('%s %d'%(utils.locstr(ploci[tch]), sum(nfo.values())) for tch, nfo in n_added.items()), totstr))
        sys.stdout.flush()
        print('                             new        existing')
        print('            singleton      cluster       cluster')
        for tch in 'hl':
            print('       %s      %7d     %7d     %7d' % (utils.locstr(ploci[tch]), n_added[tch]['singleton'], n_added[tch]['new-cluster'], n_added[tch]['existing-cluster']))
    # ----------------------------------------------------------------------------------------
    print('    merging %s partitions' % '+'.join(list(ploci.values())))
    sys.stdout.flush()
    init_partitions = {}
    for tch in utils.chains:
        if iparts is None or ploci[tch] not in iparts:
            init_partitions[tch] = cpaths[ploci[tch]].best()  # <cpaths> (and thus <init_partitions>) are after the badly paired seqs were removed, while <input_cpaths> are the real initial ones (before anything was removed)
        else:
            init_partitions[tch] = cpaths[ploci[tch]].partitions[iparts[ploci[tch]]]
            print('  %s using non-best partition index %d for %s (best is %d)' % (utils.color('red', 'note'), iparts[ploci[tch]], tch, cpaths[ploci[tch]].i_best))

    l_translations = translate_paired_uids(ploci, init_partitions, antn_lists)
    if debug:
        for tstr, tpart in [('heavy', init_partitions['h']), ('light', init_partitions['l'])]:
            ptnprint(tpart, extrastr=utils.color('blue', '%s  '%tstr), print_partition_indices=True, n_to_print=1, sort_by_size=False, print_header=tstr=='heavy')

    common_uids, _, _ = utils.check_intersection_and_complement(init_partitions['h'], init_partitions['l'], only_warn=True, a_label='heavy', b_label='light')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
    if len(common_uids) == 0:
        if all(len(init_partitions[c]) > 0 for c in 'hl'):
            print('  %s no uids in common between heavy (%d uids) and light (%d uids) partitions' % (utils.color('yellow', 'warning'), len(init_partitions['h']), len(init_partitions['l'])))
        return {ploci[ch] : [] for ch in init_partitions}

    antn_dict = {ch : utils.get_annotation_dict(antn_lists[ploci[ch]]) for ch in ploci}

    final_partition = []
    fclust_sets, fclust_indices = [], {} # just for speed
    initp_sets = {ch : [set(c) for c in ptn] for ch, ptn in init_partitions.items()}  # just for speed
    if debug:
        hdbg = ['    N        N       hclusts     lclusts       h/l',
                '  hclusts  lclusts    sizes       sizes      overlaps']
    n_huge_dbg = 10000  # this sucks and is ugly, but it's nice to get visual confirmation that things are still running (although it really only matters on like millions of sequences, so maybe not worth having it?)
    n_total_seqs = sum(len(c) for p in init_partitions.values() for c in p) + (0 if unpaired_seqs is None else len(unpaired_seqs))
    huge_dbg = n_total_seqs > n_huge_dbg  # print some extra dbg on larger samples
    if huge_dbg:
        print('      resolving %d init clusts' % sum(len(p) for p in init_partitions.values()))
        print('        finished:', end=' ')
        sys.stdout.flush()
        ihuge = 0
    # For each single cluster in each partition, get a list of the clusters in the other partition that have common uids
    # Pass this cluster + list to a fcn to resolve discrepancies by splitting on the cluster boundaries in <cluster_list> that we're sure of (i.e. that have different cdr3, or very different naive hamming fraction)
    for h_initclust, l_initclust in [(c, None) for c in init_partitions['h']] + [(None, c) for c in init_partitions['l']]:  # just loops over each single cluster in h and l partitions, but in a way that we know whether the single cluster is from h or l
        single_chain, list_chain = 'h' if l_initclust is None else 'l', 'l' if l_initclust is None else 'h'
        single_cluster = h_initclust if single_chain == 'h' else l_initclust
        single_cset = set(single_cluster)
        cluster_list = [init_partitions[list_chain][i] for i, c in enumerate(initp_sets[list_chain]) if len(single_cset & c) > 0]
        single_annotation = antn_dict[single_chain][akey(single_cluster)]
        annotation_list = [antn_dict[list_chain][akey(c)] for c in cluster_list]

        if debug and len(cluster_list) > 1:
            print('\n'.join(hdbg))
            hclusts, lclusts = ([single_cluster], cluster_list) if single_chain == 'h' else (cluster_list, [single_cluster])
            overlaps = [[len(set(hc) & set(lc)) for lc in lclusts] for hc in hclusts]
            overlapstr = '   '.join([' '.join(str(ov) for ov in ovlist) for ovlist in overlaps])
            def getcstr(clist): return ' '.join(str(len(c)) for c in clist)
            hcstr, lcstr = getcstr(hclusts), getcstr(lclusts)
            cw = 10
            if len(hcstr) < cw and len(lcstr) < cw:  # fits on a single line
                print(('    %2d      %2d         %-'+str(cw)+'s  %-'+str(cw)+'s  %s') % (len(hclusts), len(lclusts), hcstr, lcstr, overlapstr))
            else:  # split the last few columns over multiple lines
                print(('    %2d      %2d         %-s') % (len(hclusts), len(lclusts), hcstr))
                print(('    %2s      %2s         %-'+str(cw)+'s%-s') % ('', '', '', lcstr))
                print(('    %2s      %2s         %-'+str(cw)+'s%-'+str(cw)+'s   %s') % ('', '', '', '', overlapstr))
        if huge_dbg:
            if ihuge % (n_huge_dbg // 5) == 0:
                print(' %d' % ihuge, end=' ')
                sys.stdout.flush()
            ihuge += 1

        resolved_clusters = resolve_discordant_clusters([u for u in single_cluster], single_annotation, [[u for u in c] for c in cluster_list], annotation_list, tdbg=debug)
        if check_partitions:
            assert is_clean_partition(resolved_clusters)
        incorporate_rclusts(final_partition, fclust_sets, fclust_indices, resolved_clusters)
    if huge_dbg:
        print('')
    if debug:
        print('    removing %d/%d empty clusters' % (final_partition.count([]), len(final_partition)))
    final_partition = [c for c in final_partition if len(c) > 0]
    if debug:
        print('    final: %s' % utils.cluster_size_str(final_partition))
    tmpstrs = ['   N clusters without bad/unpaired seqs:'] \
              + ['%s %4d --> %-4d%s'  % (utils.locstr(ploci[tch]), len(init_partitions[tch]), len(final_partition), chstr(len(init_partitions[tch]), len(final_partition))) for tch in utils.chains]
    print('\n        '.join(tmpstrs))
    sys.stdout.flush()

    # ptnprint(final_partition, sort_by_size=False) #extrastr=utils.color('blue', '%s  '%tstr), print_partition_indices=True, n_to_print=1, sort_by_size=False, print_header=tstr=='heavy')

    if check_partitions:
        assert is_clean_partition(final_partition)
        for tch, initpart in init_partitions.items():
            _, _, _ = utils.check_intersection_and_complement(initpart, final_partition, only_warn=True, a_label=tch, b_label='joint')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
            init_ids, final_ids = set([u for c in initpart for u in c]), set([u for c in final_partition for u in c])
            if len(init_ids - final_ids) > 0:  # if any seqs from the initial partition aren't in the final partition
                raise Exception('missing %d uids from joint partition: %s' % (len(init_ids - final_ids), ' '.join(init_ids - final_ids)))
            assert len(set([u for c in initpart for u in c]) - set([u for c in final_partition for u in c])) == 0  # everybody from both initial partitions is in final_partition
        assert len(set([u for c in final_partition for u in c]) - set([u for c in init_partitions['h'] for u in c]) - set([u for c in init_partitions['l'] for u in c])) == 0  # nobody extra got added (i don't see how this could happen, but maybe it's just checking that I didnt' modify the initial partitions)

    joint_partitions = {ch : copy.deepcopy(final_partition) for ch in utils.chains}
    if len(l_translations) > 0:
        untranslate_pids(ploci, init_partitions, antn_lists, l_translations, joint_partitions, antn_dict)  # NOTE code after here (at least randomly_pair_unpaired_seqs()) assumes that corresponding h/l clusters are in same order in each partition

    if unpaired_seqs is not None:  # it might be cleaner to have this elsewhere, but I want it to happen before we evaluate, and it's also nice to have evaluation in here
        re_add_unpaired(joint_partitions, unpaired_seqs, tdbg=debug)

    ccfs = None
    if true_outfos is not None:
        assert iparts is None  # just for now
        ccfs = evaluate_joint_partitions(ploci, true_outfos, {tch : input_cpaths[ploci[tch]].best() for tch in utils.chains}, joint_partitions, antn_lists, seed_unique_ids=seed_unique_ids, unpaired_seqs=unpaired_seqs, fail_frac=fail_frac, debug=debug)

    cluster_pairs = list(zip(*[joint_partitions[c] for c in 'hl']))
    return {ploci[ch] : jp for ch, jp in joint_partitions.items()}, cluster_pairs, ccfs
