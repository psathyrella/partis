import collections
import copy
import itertools
import numpy
import sys
import operator
import string
import os
import math
import json

import hutils
import utils
import prutils
from clusterpath import ptnprint, ClusterPath
from hist import Hist
import glutils
import treeutils

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
    print '    cleaning paired files in %s/' % bdir
    fnames = paired_dir_fnames(bdir, suffix=suffix, ig_or_tr=ig_or_tr)
    if extra_files is not None:
        fnames = extra_files + fnames  # put 'em at the start, since presumably they're actual files, not dirs
    utils.clean_files(fnames, expect_missing=expect_missing)

# ----------------------------------------------------------------------------------------
# add_selection_metrics: list of selection metrics to add (plotdir_fcn is atm only if adding selection metrics)
def read_locus_output_files(tmploci, ofn_fcn, lpair=None, read_selection_metrics=False, add_selection_metrics=None, plotdir_fcn=None, seed_unique_id=None, dont_add_implicit_info=False, dbgstr='', debug=False):
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
                        print '  %s unexpected seq(s) paired with seed id %s: %s' % (utils.color('yellow', 'warning'), sid, pids)
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
                print '%s: no %s %s output file, skipping: %s' % (utils.color('blue', '+'.join(lpair) if lpair is not None else ltmp), ltmp, dbgstr, ofn)
            continue
        lpfos['glfos'][ltmp], lpfos['antn_lists'][ltmp], lpfos['cpaths'][ltmp] = utils.read_output(ofn, dont_add_implicit_info=dont_add_implicit_info, skip_failed_queries=True)
        if debug:
            print '    read %d %s annotations with %d seqs from %s' % (len(lpfos['antn_lists'][ltmp]), ltmp, sum(len(l['unique_ids']) for l in lpfos['antn_lists'][ltmp]), ofn)
        if read_selection_metrics and os.path.exists(treeutils.smetric_fname(ofn)):  # if it doesn't exist, the info should be in the regular output file
            read_smetrics(ofn, lpfos['antn_lists'][ltmp])
        if add_selection_metrics is not None:
            for smetric in add_selection_metrics:
                print '            adding selection metrics to annotations read from %s%s' % (ofn_fcn(ltmp, lpair=lpair), '' if plotdir_fcn is None else ' and plotting to %s' % plotdir_fcn(ltmp, lpair=lpair))
                treeutils.calculate_individual_tree_metrics(smetric, lpfos['antn_lists'][ltmp], base_plotdir=None if plotdir_fcn is None else plotdir_fcn(ltmp, lpair=lpair))
        parse_pairing_info(ltmp, lpfos['antn_lists'][ltmp])
    if all(lpfos['glfos'][l] is None for l in tmploci):  # if there was no info for *any* of the loci, set Nones one level up (it's just easier to have the Nones there)
        lpfos = {k : None for k in lpfos}
    return lpfos

# ----------------------------------------------------------------------------------------
def read_lpair_output_files(lpairs, ofn_fcn, read_selection_metrics=False, add_selection_metrics=None, plotdir_fcn=None, seed_unique_id=None, dont_add_implicit_info=False, dbgstr='', debug=False):
    lp_infos = {}
    for lpair in lpairs:
        lp_infos[tuple(lpair)] = read_locus_output_files(lpair, ofn_fcn, lpair=lpair, read_selection_metrics=read_selection_metrics, add_selection_metrics=add_selection_metrics, plotdir_fcn=plotdir_fcn,
                                                         seed_unique_id=seed_unique_id, dont_add_implicit_info=dont_add_implicit_info, dbgstr=dbgstr, debug=debug)
    return lp_infos

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
def write_concatd_output_files(glfos, antn_lists, ofn_fcn, headers, use_pyyaml=False, work_fnames=None, cpaths=None, true_partitions=None, dont_write_git_info=False):
    for ltmp in sorted(glfos):  # not really a reason to write igh first, but i guess it's nice to be consistent
        ofn = ofn_fcn(ltmp, joint=True)
        if utils.has_d_gene(ltmp):
            cp = ClusterPath(partition=utils.get_partition_from_annotation_list(antn_lists[ltmp])) if cpaths is None else cpaths[ltmp]
            partition_lines = cp.get_partition_lines(true_partition=None if true_partitions is None else true_partitions[ltmp], calc_missing_values='best')
            utils.write_annotations(ofn, glfos[ltmp], antn_lists[ltmp], headers, partition_lines=partition_lines, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info)
        else:
            utils.makelink(os.path.dirname(ofn), ofn_fcn(ltmp, lpair=utils.getlpair(ltmp)), ofn)
        if work_fnames is not None:
            work_fnames.append(ofn)

# ----------------------------------------------------------------------------------------
def get_antn_pairs(lpair, lpfos):  # return list of (hline, lline) pairs
    if None in lpfos.values():
        return []
    assert len(set(len(lpfos['antn_lists'][l]) for l in lpair)) == 1  # make sure the lists for both loci are the same length
    return zip(*[lpfos['antn_lists'][l] for l in lpair])

# ----------------------------------------------------------------------------------------
def get_both_lpair_antn_pairs(lpairs, lp_infos):  # ok this name sucks, but this lumps together both h+k and h+l, whereas get_antn_pairs() just gives you one or the other
    antn_pairs = []
    for lpair in lpairs:
        antn_pairs += get_antn_pairs(lpair, lp_infos[tuple(lpair)])
    return antn_pairs

# ----------------------------------------------------------------------------------------
# NOTE the deep copies are really important, since we later deduplicate the concatd partitions + annotations (they didn't used to be there, which caused a bug/crash) [well even if we didn't deduplicate, it was fucking stupid to not deep copy them)
# NOTE also that you'll probably have duplicates in both the partitions and annotations after this
def concat_heavy_chain(lpairs, lp_infos):  # yeah yeah this name sucks but i want it to be different to the one in the calling scripts
    # ----------------------------------------------------------------------------------------
    def glpf(p, k, l):  # short for "get key (k) from lp_infos for lpair (p) and locus (l) NOTE duplicates code in write_lpair_output_files
        if tuple(p) not in lp_infos or lp_infos[tuple(p)][k] is None:
            return None
        return lp_infos[tuple(p)][k].get(l)
    # ----------------------------------------------------------------------------------------
    glfos, antn_lists, joint_cpaths = {}, {}, {}
    for lpair in lpairs:
        for ltmp in lpair:
            if glpf(lpair, 'glfos', ltmp) is None:  # this lpair's output files were empty
                continue
            if ltmp in glfos:  # heavy chain, second time through: merge chain glfos from those paired with igk and with igl
                glfos[ltmp] = glutils.get_merged_glfo(glfos[ltmp], glpf(lpair, 'glfos', ltmp))
                antn_lists[ltmp] += copy.deepcopy(glpf(lpair, 'antn_lists', ltmp))
                if glpf(lpair, 'cpaths', ltmp) is not None:
                    assert len(joint_cpaths[ltmp].partitions) == 1  # they should always be 1 anyway, and if they weren't, it'd make it more complicated to concatenate them
                    joint_cpaths[ltmp] = ClusterPath(partition=joint_cpaths[ltmp].best() + copy.deepcopy(glpf(lpair, 'cpaths', ltmp).best()))
            else:  # light + heavy chain first time through
                glfos[ltmp] = copy.deepcopy(glpf(lpair, 'glfos', ltmp))
                antn_lists[ltmp] = copy.deepcopy(glpf(lpair, 'antn_lists', ltmp))
                if glpf(lpair, 'cpaths', ltmp) is not None:
                    joint_cpaths[ltmp] = copy.deepcopy(glpf(lpair, 'cpaths', ltmp))
    return glfos, antn_lists, joint_cpaths

# ----------------------------------------------------------------------------------------
# NOTE i'm really not sure that this fcn needs to exist -- don't the joint partitions for h and l end up ordered i.e. with each cluster tied to its partner? Or at least I should be able to keep track of who goes with who so I don't need to reconstruct it here
def find_cluster_pairs(lp_infos, lpair, antn_lists=None, required_keys=None, quiet=False, debug=False):  # the annotation lists should just be in the same order, but after adding back in all the unpaired sequences to each chain they could be a bit wonky
    # can set <antn_lists> if you don't already have lp_infos (and don't need glfos)
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
        return all_ids
    # ----------------------------------------------------------------------------------------
    if antn_lists is not None:
        assert lp_infos is None
        lp_infos = {tuple(lpair) : {'antn_lists' : {l : antn_lists[l] for l in lpair}, 'cpaths' : {ltmp : ClusterPath(partition=[l['unique_ids'] for l in antn_lists[ltmp]]) for ltmp in lpair}}}
    if required_keys is None:
        required_keys = []
    if 'paired-uids' not in required_keys:
        required_keys.append('paired-uids')

    lp_antn_pairs = []
    lpk = tuple(lpair)
    if None in lp_infos[lpk].values():
        return lp_antn_pairs
    h_part, l_part = [sorted(lp_infos[lpk]['cpaths'][l].best(), key=len, reverse=True) for l in lpair]
    h_atn_dict, l_atn_dict = [utils.get_annotation_dict(lp_infos[lpk]['antn_lists'][l], cpath=lp_infos[lpk]['cpaths'][l]) for l in lpair]
    if debug:
        print '  finding cluster pairs for %s partitions with cluster sizes:\n     %s: %s\n     %s: %s' % ('+'.join(lpair), lpair[0], ' '.join(str(len(c)) for c in h_part), lpair[1], ' '.join(str(len(c)) for c in l_part))
        print '          sizes'
        print '          h   l    l index'
    n_skipped = {k : 0 for k in required_keys + ['zero-len-paired-uids']}
    unpaired_l_clusts = [c for c in l_part]
    for h_clust in h_part:
        h_atn = h_atn_dict[':'.join(h_clust)]

        if any(k not in h_atn for k in required_keys):  # skip any annotations that are missing any of these keys (atm, only used to skip ones without 'tree-info', which usually means clusters that were smaller than min selection metric cluster size
            for rk in set(required_keys) - set(h_atn):
                n_skipped[rk] += 1
            continue
        if len(getpids(h_atn)) == 0:
            n_skipped['zero-len-paired-uids'] += 1
            continue

        l_clusts = [c for c in l_part if len(set(getpids(h_atn)) & set(c)) > 0]
        if len(l_clusts) != 1:
            if not quiet:
                print '  %s couldn\'t find a unique light cluster (found %d, looked in %d) for heavy cluster with size %d and %d paired ids (heavy: %s  pids: %s)' % (utils.color('yellow', 'warning'), len(l_clusts), len(l_part), len(h_clust), len(getpids(h_atn)), ':'.join(h_clust), ':'.join(getpids(h_atn)))
            continue
        assert len(l_clusts) == 1
        if ':'.join(l_clusts[0]) not in l_atn_dict:
            print '      %s missing annotation for light chain when finding cluster pairs: %s (paired with %s)' % (utils.color('yellow', 'warning'), ':'.join(l_clusts[0]), ':'.join(h_clust))
            unpaired_l_clusts.remove(l_clusts[0])  # i guess i want to remove it from here? i guess we know who it's paired with, but there's no annotation so we can't do anything with it
            continue
        l_atn = l_atn_dict[':'.join(l_clusts[0])]
        h_atn['loci'] = [lpair[0] for _ in h_atn['unique_ids']]  # this kind of sucks, but it seems like the best option a.t.m. (see note in event.py)
        l_atn['loci'] = [lpair[1] for _ in l_atn['unique_ids']]
        lp_antn_pairs.append((h_atn, l_atn))
        unpaired_l_clusts.remove(l_clusts[0])
        if debug:
            print '        %3d %3d   %3d' % (len(h_clust), len(l_clusts[0]), l_part.index(l_clusts[0]))
    if len(unpaired_l_clusts) > 0:
        print '    %s: %d unpaired light cluster%s after finding h/l cluster pairs' % ('+'.join(lpair), len(unpaired_l_clusts), utils.plural(len(unpaired_l_clusts)))
        # this is just too verbose atm (and hopefully not necessary?)
        # for lc in unpaired_l_clusts:
        #     if ':'.join(lc) not in l_atn_dict:
        #         print '        %s missing annotation for unpaired light chain when finding cluster pairs: %s' % (utils.color('yellow', 'warning'), ':'.join(lc))
        #         continue
        #     lpids = getpids(l_atn_dict[':'.join(lc)])
        #     hpclusts = [c for c in h_part if len(set(lpids) & set(c)) > 0]0
        #     if len(hpclusts) > 0:  # i think this would mean that the pairing info was non-reciprocal, which probably isn't really possible?
        #         print '       %s unpaired light cluster with size %d overlaps with heavy cluster(s): %s' % (utils.color('yellow', 'warning'), len(lc), ' '.join(str(len(c)) for c in hpclusts))
    if any(n > 0 for k, n in n_skipped.items() if k!='zero-len-paired-uids'):
        print '    %s: skipped %d annotations missing required keys: %s' % ('+'.join(lpair), sum(n_skipped.values()), '  '.join('%s: %d'%(k, n) for k, n in sorted(n_skipped.items()) if n>0 and k!='zero-len-paired-uids'))
    if n_skipped['zero-len-paired-uids'] > 0:
            print '    %s: skipped %d annotations with zero length paired uids' % ('+'.join(lpair), n_skipped['zero-len-paired-uids'])
    if debug:
        print '  '
    return lp_antn_pairs

# ----------------------------------------------------------------------------------------
def apportion_cells_to_droplets(outfos, metafos, mean_cells_per_droplet):
    n_droplets = int(0.5 * float(len(outfos)) / mean_cells_per_droplet)  # (randomly) apportion cells among this many droplets (0.5 is because <outfos> includes both heavy and light sequences)
    droplet_ids = [[] for _ in range(n_droplets)]  # list of sequence ids for each droplet
    sfo_dict = {s['name'] : s for s in outfos}  # temp, to keep track of who still needs apportioning (but we do modify its sfos, which are shared with <outfos>)
    while len(sfo_dict) > 0:
        tid = next(iter(sfo_dict))
        idrop = numpy.random.choice(range(len(droplet_ids)))
        droplet_ids[idrop] += [tid] + metafos[tid]['paired-uids']  # add <tid> plus its paired ids to this drop (note that these are the original/correct paired ids, which is what we want)
        for uid in [tid] + metafos[tid]['paired-uids']:
            sfo_dict[uid]['droplet-ids'] = droplet_ids[idrop]
            del sfo_dict[uid]
    for sfo in outfos:
        metafos[sfo['name']]['paired-uids'] = [u for u in sfo['droplet-ids'] if u != sfo['name']]
    print '  apportioned %d seqs among %d droplets (mean/2 %.1f): %s' % (len(outfos), n_droplets, numpy.mean([len(d) for d in droplet_ids]) / 2, ' '.join(str(len(d)) for d in droplet_ids))
# ----------------------------------------------------------------------------------------
def remove_reads_from_droplets(outfos, metafos, fraction_of_reads_to_remove):
    n_to_remove = int(fraction_of_reads_to_remove * len(outfos))
    ifos_to_remove = numpy.random.choice(range(len(outfos)), size=n_to_remove, replace=False)
    for ifo in ifos_to_remove:
        del metafos[outfos[ifo]['name']]
    outfos = [outfos[ifo] for ifo in range(len(outfos)) if ifo not in ifos_to_remove]
    print '  removed %d / %d = %.2f seqs from outfos' % (n_to_remove, len(outfos) + n_to_remove, n_to_remove / float(len(outfos) + n_to_remove))
    return outfos
# ----------------------------------------------------------------------------------------
# write fasta and meta file with all simulation loci together
def write_merged_simu(antn_lists, fastafname, metafname, mean_cells_per_droplet=None, fraction_of_reads_to_remove=None):  # NOTE that this writes a new input meta info file, which is where partis will then get the paird uid info if --input-metfname is set, but does *not* modify the 'paired-uids' key in the original simulation files (since we want those to be correct even if we're adding extra/removing cells from droplets)
    # merge together info from all loci into <outfos> and <metafos>
    outfos, metafos = [], {}
    for ltmp in antn_lists:
        for tline in antn_lists[ltmp]:
            for uid, seq, pids in zip(tline['unique_ids'], tline['input_seqs'], tline['paired-uids']):
                outfos.append({'name' : uid, 'seq' : seq})
                metafos[uid] = {'locus' : ltmp, 'paired-uids' : pids}

    if mean_cells_per_droplet is not None:
        apportion_cells_to_droplets(outfos, metafos, mean_cells_per_droplet)
    if fraction_of_reads_to_remove is not None:
        outfos = remove_reads_from_droplets(outfos, metafos, fraction_of_reads_to_remove)

    # write merged fasta and input meta files
    with open(fastafname, 'w') as outfile:
        for sfo in outfos:
            outfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
    with open(metafname, 'w') as mfile:
        json.dump(metafos, mfile)

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

    for ch in antn_dict:  # atm this only gets used for dbg, but still nice to properly fix it
        antn_dict[ch] = utils.get_annotation_dict(antn_lists[ploci[ch]])

# ----------------------------------------------------------------------------------------
def remove_badly_paired_seqs(ploci, outfos, debug=False):  # remove seqs paired with the other/wrong light chain, as well as those with no pairing info (the latter we keep track of so we can insert them later into the right final cluster)
    # ----------------------------------------------------------------------------------------
    def add_unpaired(cline, iseq, uid):
        sorted_hdists = sorted([(u, utils.hamming_distance(cline['seqs'][i], cline['seqs'][iseq])) for i, u in enumerate(cline['unique_ids']) if i != iseq], key=operator.itemgetter(1))
        nearest_uid = sorted_hdists[0][0] if len(sorted_hdists) > 0 else None
        unpaired_seqs[cline['loci'][iseq]][uid] = nearest_uid  # unless there's no other seqs in the cluster, attach it to the nearest seq by hamming distance
    # ----------------------------------------------------------------------------------------
    cpaths, antn_lists, glfos = [outfos[k] for k in ['cpaths', 'antn_lists', 'glfos']]
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l]) for l in antn_lists}
    all_loci = {u : l for l, ants in antn_lists.items() for antn in ants for u in antn['unique_ids']}  # this includes the heavy ones, which we don't need, but oh well
    all_pids = {u : pids[0] for alist in antn_lists.values() for l in alist for u, pids in zip(l['unique_ids'], l['paired-uids']) if len(pids)==1}  # uid : pid for all uid's that have a single unique pid (which should be all of them, since we just ran pair cleaning -- otherwise we crash below) (I'm pretty sure that the partition implied by the annotations is identical to the one in <cpaths>, and it's nice to loop over annotations for this)
    unpaired_seqs = {l : {} for l in ploci.values()}  # map for each locus from the uid of each seq with no (or non-reciprocal) pairing info to the nearest sequence in its family (after merging partitions we'll insert it into the family that this nearest seq ended up in)
    lp_cpaths, lp_antn_lists = {}, {}
    if debug:
        print '  removing bad/un-paired seqs'
        print '          N       N      no   other  non-'
        print '        before removed  info  light recip'
    for tch in sorted(ploci):
        new_partition, new_antn_list = [], []
        for iclust, cluster in enumerate(cpaths[ploci[tch]].best()):
            cline = antn_dicts[ploci[tch]][':'.join(cluster)]
            iseqs_to_remove = []
            n_no_info, n_other_light, n_non_reciprocal = 0, 0, 0  # just for dbg NOTE n_other_light are the only ones we *really* want to remove, since they're h seqs paired with the other light chain, whereas the other two categories we eventually want to re-add since we're not sure who they're paired with
            for iseq, uid in enumerate(cline['unique_ids']):
                pids = cline['paired-uids'][iseq]
                if len(pids) == 0:  # no pairing info
                    iseqs_to_remove.append(iseq)
                    add_unpaired(cline, iseq, uid)
                    n_no_info += 1
                elif len(pids) > 1:  # shouldn've all been removed by pair info cleaning
                    raise Exception('multiple paired uids for \'%s\': %s' % (uid, pids))
                else:
                    if tch == 'h' and all_loci[utils.get_single_entry(pids)] != ploci['l']:  # if it's the other light chain
                        iseqs_to_remove.append(iseq)
                        n_other_light += 1
                    else:  # also remove any non-reciprocal pairings (I think this will still miss any whose partner was removed) NOTE it would be nice to enforce reciprocal pairings in clean_pair_info(), but atm i think we can't look at both chains at once in that fcn
                        if all_pids[uid] not in all_pids or all_pids[all_pids[uid]] != uid:  # if uid's pid isn't in all_pids, or if it is but it's a different uid
                            iseqs_to_remove.append(iseq)
                            add_unpaired(cline, iseq, uid)
                            n_non_reciprocal += 1
            iseqs_to_keep = [i for i in range(len(cline['unique_ids'])) if i not in iseqs_to_remove]
            if len(iseqs_to_keep) > 0:
                new_partition.append([cluster[i] for i in iseqs_to_keep])
                new_cline = utils.get_non_implicit_copy(cline)
                utils.restrict_to_iseqs(new_cline, iseqs_to_keep, glfos[ploci[tch]])
                new_antn_list.append(new_cline)
            if debug:
                def fstr(v): return '' if v==0 else '%d'%v
                print '    %s   %3d    %3s     %3s   %3s    %3s' % (utils.locstr(ploci[tch]) if iclust==0 else ' ', len(cline['unique_ids']), fstr(len(iseqs_to_remove)), fstr(n_no_info), fstr(n_other_light), fstr(n_non_reciprocal))
        lp_cpaths[ploci[tch]] = ClusterPath(seed_unique_id=cpaths[ploci[tch]].seed_unique_id, partition=new_partition)
        lp_antn_lists[ploci[tch]] = new_antn_list

    if debug:
        print '    totals before: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in cpaths[ploci[tch]].best())) for tch in sorted(ploci))
        print '    totals after: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in lp_cpaths[ploci[tch]].best())) for tch in sorted(ploci))

    return lp_cpaths, lp_antn_lists, unpaired_seqs

# ----------------------------------------------------------------------------------------
def clean_pair_info(cpaths, antn_lists, is_data=False, plotdir=None, performance_outdir=None, paired_data_type='10x', collapse_similar_paired_seqs=False, max_hdist=4, debug=False):
    # ----------------------------------------------------------------------------------------
    def check_droplet_id_groups(pid_groups, all_uids, tdbg=False):
        if not is_data:
            # print '  note: couldn\'t get droplet id from \'%s\', so assuming this isn\'t 10x data' % next(iter(all_uids))  # NOTE i'm not sure that this gives the same one as the previous line
            return False
        if len(set(len(g) for g in pid_groups)) == 1:
            print '    %s all pid groups are length 1 in check_droplet_id_groups(). Maybe you\'re missing pairing info?' % utils.color('yellow', 'warning')
            return False
        # check against the droplet id method (we could just do it this way, but it would only work for 10x, and only until they change their naming convention)
        pgroup_strs = set(':'.join(sorted(pg)) for pg in pid_groups)  # <pid_groups>: list of pid groups, i.e. each element is the uids from a single droplet (for 10x), so <pgroup_strs> converts each group into an indexable key
        n_not_found = 0
        if tdbg:
            print '              found?   drop id           contigs     overlaps (with any non-identical groups)'
        def kfcn(u): return utils.get_droplet_id(u, dtype=paired_data_type)
        for dropid, drop_queries in itertools.groupby(sorted(all_uids, key=kfcn), key=kfcn):  # group all queries into droplets (assuming they conform to 10x convention for droplet ids)
            dqlist = list(drop_queries)
            found = ':'.join(sorted(dqlist)) in pgroup_strs  # was this exact combination of queries in pid_groups?
            if not found:  # if not, see if any pid groups have some of these queries
                overlaps = [g for g in pgroup_strs if dropid in g]  # this should essentially always be length 1, except when we're missing pairing info in simulation (in which case i know we don't even want to be running this simulation, but whatever I'm testing edge cases)
                n_not_found += 1
            if tdbg or not found:
                ostr = ' '.join(sorted(utils.get_contig_id(q) for q in overlaps[0].split(':'))) if len(overlaps)==1 else 'multiple'
                print '  %25s    %s   %-8s   %s' % (utils.color('green', '-') if found else utils.color('red', 'x'), dropid, ' '.join(sorted(utils.get_contig_id(q) for q in dqlist)),
                                                    utils.color('red', ostr if not found else ''))
        if n_not_found > 0:
            print '  %s droplet id group check failed for %d groups, i.e. droplet ids parsed from uids don\'t match pair info: either pairing info is messed up or missing, or this is simulation and you didn\'t set --is-simu (if the latter, ignore this)' % (utils.color('red', 'error'), n_not_found)
        return True
    # ----------------------------------------------------------------------------------------
    def plot_uids_before(plotdir, pid_groups, all_antns):
        # ----------------------------------------------------------------------------------------
        def fnfplot(logstr, fhists, n_max_bins=15):
            import plotting
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
        binlabels, _ = zip(*sorted(flfracs.items(), key=operator.itemgetter(1), reverse=True))

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
        pidlengths = {}
        for ltmp in sorted(cpaths):
            for cluster in cpaths[ltmp].best():
                atn = antn_dicts[ltmp][':'.join(cluster)]
                pidlengths.update({u : len(set(pids) - set([u])) for u, pids in zip(atn['unique_ids'], atn['paired-uids'])})
        ahist = Hist(value_list=pidlengths.values(), init_int_bins=True)
        fn = ahist.fullplot(plotdir, 'paired-seqs-per-seq-%s'%pstr, pargs={'remove_empty_bins' : True}, fargs={'xlabel' : 'N paired seqs per seq', 'ylabel' : 'counts', 'title' : pstr, 'xbounds' : (-0.05, 1.05*ahist.xmax), 'xticks' : [i for i in range(0, int(ahist.xmax+1), 1 if 'after' in pstr else 2)]})
        fnames[1].append(fn)
        return pidlengths
    # ----------------------------------------------------------------------------------------
    def make_fraction_correct_plot():
        # ----------------------------------------------------------------------------------------
        def splid(utmp):
            # ustr, ltmp = utmp.split('-')
            # assert ltmp in utils.loci
            # return ustr
            ltmp = utils.get_single_entry([l for l in utils.loci if '-'+l in utmp])  # ICK (works also for bcr-phylo)
            return utmp.replace('-'+ltmp, '-LOCUS')
        # ----------------------------------------------------------------------------------------
        def gpt(uid, pids):
            if len(pids) == 0:
                return 'unpaired'
            elif len(pids) > 1:
                return 'multiple'
            elif splid(uid) == splid(pids[0]):
                return 'correct'
            else:
                return 'mispaired'
        # ----------------------------------------------------------------------------------------
        fcinfo = collections.OrderedDict([('correct', 0), ('mispaired', 0), ('unpaired', 0), ('multiple', 0)])
        for ltmp in sorted(cpaths):
            for cluster in cpaths[ltmp].best():
                atn = antn_dicts[ltmp][':'.join(cluster)]
                for uid, pids in zip(atn['unique_ids'], atn['paired-uids']):
                    fcinfo[gpt(uid, pids)] += 1
        fchist = hutils.make_hist_from_dict_of_counts(fcinfo, 'string', 'pair cleaning performance', no_sort=True)
        fchist.normalize()
        fcplname = 'true-pair-clean-performance'
        if performance_outdir is not None:
            fchist.write('%s/%s.csv'%(performance_outdir, fcplname))
        if plotdir is not None:
            fn = fchist.fullplot(plotdir, fcplname, pargs={'ignore_overflows' : True}, fargs={'xbounds' : (0.95, 1.05*len(fcinfo)), 'ybounds' : (0., 1.05), 'xticklabelsize' : 15, 'ylabel' : 'fraction of seqs'})
            fnames.append([fn])
    # ----------------------------------------------------------------------------------------
    def make_final_plots(initial_seqs_per_seq, initial_flcounts):
        final_seqs_per_seq = plot_n_pseqs_per_seq('after')
        import plotting
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
        if uid not in all_antns:
            return '?'
        return utils.per_seq_val(all_antns[uid], 'loci', uid)
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
            print '      removed %d with missing annotations' % len(ids_to_remove)

        # among any pairs of sequences that are [almost] identical at all non-ambiguous position, keep only the longest one (note that this is really preprocessing/error correction, so probably shouldn't really be here)
        if collapse_similar_paired_seqs:
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
                print '        %d pair%s equivalent with hdists %s' % (n_equivalent, utils.plural(n_equivalent), ' '.join(dbgstr))

        # if specified, remove unproductive (only on real data, since simulation usually has lots of stop codons)
        if is_data and remove_unproductive:
            dbgstr = []
            unproductive_ids = []
            for uid in chain_ids:
                if not utils.is_functional(all_antns[uid], all_antns[uid]['unique_ids'].index(uid)):
                    unproductive_ids.append(uid)
                    if tdbg:
                        dbgstr.append(utils.is_functional_dbg_str(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), sep='+'))
            ids_to_remove |= set(unproductive_ids)
            if tdbg and len(unproductive_ids) > 0:
                print '        %d unproductive  %s' % (len(unproductive_ids), ',  '.join(dbgstr))

        return ids_to_remove

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
    def clean_with_partition_info(cluster):  # use information from the [clonal family] partitions to decide which of several potential paired uids is the correct one
        # NOTE it would probably be better to structure this so that within a droplet, we first do seqs that we think we'll be able to correctly pair (e.g. the ones without ties). Maybe just keep looping over iseqs until we don't add any new pair info?
        # ----------------------------------------------------------------------------------------
        def get_pfamily_dict(cline, extra_str='', old_pfam_dict=None):  # see what others in its family are paired with
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
            if debug and len(cline['unique_ids']) > 1:
                print '    %7s votes size  id  cdr3' % extra_str
                for fkey, fdict in sorted(pfdict.items(), key=lambda x: x[1]['count'], reverse=True):
                    print '       %s    %3d   %3d  %2d %s  %3d' % (utils.locstr(fdict['locus']), fdict['count'], len(antn_dicts[fdict['locus']][fkey]['unique_ids']), fdict['id'], fidstr(fdict['id']), antn_dicts[fdict['locus']][fkey]['cdr3_length'])

            return pfdict
        # ----------------------------------------------------------------------------------------
        def update_pid_info(cpids, tdbg=False):  # remove the two ids <cpids> from everybody's paired ids, since we've just decided they're properly paired (note that you do *not* want to modify <pfamilies> here -- that would be destroying the info you need to decided who's paired with who)
            if tdbg:
                print '  upd: %s' % ' '.join(cpids)
            for cid in cpids:
                for iseq, uid in enumerate(cline['unique_ids']):  # remove <cid> from all paired uid lists (except the other one in <cpids>)
                    if cid not in cline['paired-uids'][iseq]:
                        continue
                    if uid in cpids:
                        continue
                    if tdbg:
                        print '    %s %3d --> %3d   %s' % (uid, len(cline['paired-uids'][iseq]), len([p for p in cline['paired-uids'][iseq] if p != cid]), ' '.join(utils.color('red' if p==cid else None, p) for p in cline['paired-uids'][iseq]))  # ok i know it should always just decrease by one, but maybe something else could break and you get duplicates?
                    cline['paired-uids'][iseq] = [p for p in cline['paired-uids'][iseq] if p != cid]
                    # if len(cline['paired-uids'][iseq]) == 0:
                    #     raise Exception('removed all paired ids')
        # ----------------------------------------------------------------------------------------
        def lcstr(pids, pfcs, pfids=None):  # returns string summarizing the families of the paired uids for a uid, e.g. 'k 51  l 3  h 1' if the uid has three potential pids, one from k with which 50 other uids in <cline> are paired, etc.
            if len(pids) == 0: return ''
            if pfids is None:
                spids, spfcs = zip(*sorted(zip(pids, pfcs), key=operator.itemgetter(1), reverse=True))
                spfids = ['' for _ in spids]
            else:
                assert len(pids) == 1  # if we passed in <pfids> (an id for each paired family), this should be after cleaning, so there should only be one of them
                spids, spfcs, spfids = pids, pfcs, pfids
            return '  '.join('%s %s %s'%(lg, '%1d'%sp if pfids is None else '%2d'%sp, fidstr(fid)) for lg, sp, fid in zip(lgstr(spids, dont_sort=True).split(' '), spfcs, spfids))
        # ----------------------------------------------------------------------------------------
        def print_dbg():
            new_pfams = get_pfamily_dict(cline, extra_str='after:', old_pfam_dict=old_pfams)
            pfcounts = [[new_pfams[pfkey(p)]['count'] for p in pids] for pids in cline['paired-uids']]
            pfids = [[new_pfams[pfkey(p)]['id'] for p in pids] for pids in cline['paired-uids']]
            uid_extra_strs = ['%s: %s'%(utils.locstr(l), lcstr(pids, pfcs, pfids)) for l, pids, pfcs, pfids in zip(cline['loci'], cline['paired-uids'], pfcounts, pfids)]
            old_pfcounts = [[old_pfams[pfkey(p)]['count'] for p in pids] for pids in old_pids]
            old_estrs = ['%s: %s'%(utils.locstr(l), lcstr(pids, pfcs)) for l, pids, pfcs in zip(cline['loci'], old_pids, old_pfcounts)]
            for istr, (oldstr, newstr) in enumerate(zip(old_estrs, uid_extra_strs)):
                if newstr != oldstr:
                    uid_extra_strs[istr] = '%s%s (%s)' % (newstr, ' '*(12 - utils.len_excluding_colors(newstr)), oldstr)
            utils.print_reco_event(cline, uid_extra_strs=uid_extra_strs, extra_str='      ')
            print ''
        # ----------------------------------------------------------------------------------------
        def pfkey(p): return ':'.join(all_antns[p]['unique_ids'])  # keystr for the family of paired id <p>
        # ----------------------------------------------------------------------------------------
        cline = antn_dicts[ltmp][':'.join(cluster)]
        if any(u not in pid_groups[pid_ids[u]] for u in cline['unique_ids']):  # shouldn't be able to happen any more, but that was a really bad/dumb bug
            raise Exception('one of unique ids %s not in its own pid group' % cline['unique_ids'])
        cline['paired-uids'] = [[p for p in pid_groups[pid_ids[u]] if p != u] for u in cline['unique_ids']]  # re-set 'paired-uids' key in <cline> to the pid group (which was modified in the previous cleaning steps)

        if debug:
            old_pids = copy.deepcopy(cline['paired-uids'])
        old_pfams = get_pfamily_dict(cline, extra_str='before:')  # map from each potential paired family to the number of uids in <cluster> that are potentially paired with it (i.e. the number of uids that are voting for it)

        # for each uid, choose the pid that's of opposite chain, and has the most other uids voting for it (as long as some criteria are met)
        for iseq, uid in enumerate(cline['unique_ids']):
            pid_to_keep = None
            ochain_pidfcs = [(p, old_pfams[pfkey(p)]['count']) for p in cline['paired-uids'][iseq] if not utils.samechain(getloc(p), getloc(uid))]  # (pid, pcount) for all opposite-chain pids, where <pcount> is the number of votes for <pid>'s family (note that the 'paired-uids' get modified as we go through the loop)
            if len(ochain_pidfcs) > 0:
                sorted_pids, sorted_pfcs = zip(*sorted(ochain_pidfcs, key=operator.itemgetter(1), reverse=True))
                # note that even if there's only one ochain choice, there can be other same-chain ones that we still want to drop (hence the <2 below)
                if len(sorted_pfcs) < 2 or sorted_pfcs[0] > sorted_pfcs[1] or pfkey(sorted_pids[0]) == pfkey(sorted_pids[1]):  # in order to drop the later ones, the first one either has to have more counts, or at least the second one has to be from the same family (in the latter case we still don't know which is the right one, but for the purposes of clustering resolution we just need to know what the family is)
                    pid_to_keep = sorted_pids[0]
                    update_pid_info([uid, pid_to_keep])
            cline['paired-uids'][iseq] = [pid_to_keep] if pid_to_keep is not None else []  # if we didn't decide on an opposite-chain pid, remove all pairing info

        if debug: # and len(cline['unique_ids']) > 1:  # NOTE it's annoying printing the singletons, but it's way worse when they're just missing and you can't figure out where a sequence went
            print_dbg()

    # ----------------------------------------------------------------------------------------
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l], cpath=cpaths[l]) for l in antn_lists}
    all_uids = set(u for p in cpaths.values() for c in p.best() for u in c)  # all uids that occur in a partition (should I think be the same as the ones for which we have valid/non-failed annotations)

    # first collect some information for later use
    pid_groups = []  # list of pid groups, i.e. each element is the uids from a single droplet (for 10x)
    pid_ids = {}  # map from each uid to the index of its pid group
    all_antns = {}  # map from each individual uid to its annotation
    n_missing = 0
    if debug:
        print '  %s consolidating info for %d loci with family/sequence counts: %s' % (utils.color('blue', '+'.join(sorted(cpaths))), len(cpaths), '  '.join('%s: %d/%d'%(l, len(cpaths[l].best()), sum(len(c) for c in cpaths[l].best())) for l in sorted(cpaths)))
    for ltmp in sorted(cpaths):
        for cluster in cpaths[ltmp].best():
            cline = antn_dicts[ltmp][':'.join(cluster)]
            if 'paired-uids' not in cline:
                print '  %s no paired uids in line' % utils.color('yellow', 'warning')
                continue  # maybe should still add to all_antns?
            for uid, pids in zip(cline['unique_ids'], cline['paired-uids']):
                if pids is None:
                    raise Exception('None type \'paired-uids\' for %s' % uid)
                missing_ids = set(pids) - all_uids
                n_missing += len(missing_ids)
                pset = set([uid] + pids) - missing_ids
                found = False
                for ipg, pgroup in enumerate(pid_groups):  # look for an existing pid group with some overlap
                    if any(p in pgroup for p in pset):  # could maybe check for consistency if some of them are already in there (i.e. from reciprocal info in another chain)?
                        found = True
                        pgroup |= pset
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
        print '   %d/%d missing uids when cleaning pair info' % (n_missing, len(all_uids))  # NOTE at least for now we're skipping invalid queries when reading output
    # for ipg, pg in enumerate(pid_groups):
    #     print '  %3d %s' % (ipg, ' '.join(pg))

    idg_ok = check_droplet_id_groups(pid_groups, all_uids)  # NOTE not using the return value here, but I may need to in the future
    if plotdir is not None:
        fnames = [[], [], []]
        initial_flcounts = plot_uids_before(plotdir, pid_groups, all_antns)
        initial_seqs_per_seq = plot_n_pseqs_per_seq('before')

    # then go through each group trying to remove as many crappy/suspicously similar ones as possible (this step has a fairly minor effect compared to the partition-based step below)
    if debug:
        print '  cleaning %d pid groups:' % len(pid_groups)
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
            print '    %s' % lgstr(pgroup),
        for chain, idlist in zip(utils.chains, [hids, lids]):
            if len(idlist) < 2:  # skip whichever of the chains has only one id
                continue
            if debug > 1:
                print '\n      too many %s chains: %s' % (chain, lgstr(idlist))
            ids_to_remove = choose_seqs_to_remove(idlist, tdbg=debug>1)
            if debug and len(ids_to_remove) > 0:
                tmpincr(pgroup, id_removed_groups)
            for rid in ids_to_remove:
                pgroup.remove(rid)
                idlist.remove(rid)
                pid_groups.append(set([rid]))  # add the removed id to a new pid group of its own at the end (so it'll show up as unpaired)
                pid_ids[rid] = len(pid_groups) - 1
            if debug > 1:
                print '      %s: removed %d, leaving %d%s' % (utils.color('green', 'fixed') if len(idlist)==1 else utils.color('red', 'still too many'), len(ids_to_remove), len(idlist), ':' if len(idlist)>1 else '')
                if len(idlist) > 1:
                    for uid in idlist:
                        prutils.print_seq_in_reco_event(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), one_line=True, extra_str='        ', uid_extra_str=utils.locstr(getloc(uid)))

        pid_groups[ipg] = pgroup
        if debug: tmpincr(pgroup, tried_to_fix_groups)

    if debug:
        def prcgrps(cdict, prestr):
            print '    %s' % prestr
            for lstr, count in sorted(cdict.items(), key=operator.itemgetter(1), reverse=True):
                print '      %3d  %s' % (count, lstr)
        prcgrps(ok_groups, 'ok to start with:')
        prcgrps(id_removed_groups, 'removed ids from:')
        prcgrps(tried_to_fix_groups, 'after tring to fix:')

    # then go through again using cluster/family information to try to cut everyone down to one paired id (and if we can't get down to one, we remove all of them) NOTE also re-sets the actual 'paired-uids' keys
    for ltmp in sorted(cpaths):
        if debug:
            print '%s' % utils.color('green', ltmp)
            name_dict, name_ids = {'potential' : None, 'used' : None}, {}
        for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
            clean_with_partition_info(cluster)

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
    if n_fixed > 0:
        print '     synchronized/fixed %d pairs where one had no pair info after cleaning: %s' % (sum(n for n in n_fixed.values()), '  '.join('%s %d'%(utils.locstr(l), n_fixed[l]) for l in sorted(n_fixed)))

    if not is_data and (performance_outdir is not None or plotdir is not None):
        make_fraction_correct_plot()
    if plotdir is not None:
        make_final_plots(initial_seqs_per_seq, initial_flcounts)

# ----------------------------------------------------------------------------------------
def compare_partition_pair(cfpart, refpart, remove_from_ref=False, add_to_ref=False, antn_list=None, seed_unique_id=None, dbg_str=None, cf_label='inferred', ref_label='true', debug=False):
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
            print '%s duplicate sequences in joint partition evaluation' % utils.color('yellow', 'warning')

    cfpart = copy.deepcopy(cfpart)  # don't want to modify opject pointed to by <cfpart)
    if len(dup_dict) > 0:
        incorporate_duplicates(cfpart)
    if remove_from_ref:
        refpart = utils.remove_missing_uids_from_partition(refpart, cfpart, debug=debug)  # returns a new/copied partition, doesn't modify original
    if add_to_ref:
        refpart = utils.add_missing_uids_to_partition(refpart, cfpart, miss_label=ref_label, ref_label=cf_label)  # NOTE order is reversed here
    return utils.per_seq_correct_cluster_fractions(cfpart, refpart, seed_unique_id=seed_unique_id, dbg_str=dbg_str, inf_label=cf_label, true_label=ref_label, debug=debug)
    # TODO figure out which cases of 'missing' uids should really be removed, and which should be singletons

# ----------------------------------------------------------------------------------------
def evaluate_joint_partitions(ploci, true_partitions, init_partitions, joint_partitions, antn_lists, seed_unique_ids=None, debug=False):
    if seed_unique_ids is not None and set(ploci.values()) != set(seed_unique_ids):  # skip non-seed light locus pair
        return
    # NOTE that <joint_partitions> can have many fewer seqs than <init_partitions> since in making <joint_partitions> we remove seqs paired with the other light chain (the weighted average ccfs over the h joint partitions corresponding to both light chains would be exactly comparable to the <init_partitions>, but I think this is fine as it is)
    for tch in utils.chains:
        ltmp = ploci[tch]
        ccfs = {}
        for dstr, cfpart in [('single', init_partitions[tch]), ('joint', joint_partitions[tch])]:
            ccfs[dstr] = compare_partition_pair(cfpart, true_partitions[ltmp], seed_unique_id=None if seed_unique_ids is None else seed_unique_ids[ltmp], remove_from_ref=True,  # removes from the true ptn any uids that are missing from the inferred ptn
                                                antn_list=antn_lists[ltmp], dbg_str='%s %s '%(utils.locstr(ltmp), dstr), debug=debug)
        print '  %s ccfs:     purity  completeness' % utils.locstr(ltmp)
        print '      single  %6.3f %6.3f' % (ccfs['single'][0], ccfs['single'][1])
        print '       joint  %6.3f %6.3f' % (ccfs['joint'][0], ccfs['joint'][1])

# ----------------------------------------------------------------------------------------
# cartoon explaining algorithm here https://github.com/psathyrella/partis/commit/ede140d76ff47383e0478c25fae8a9a9fa129afa#commitcomment-40981229
def merge_chains(ploci, cpaths, antn_lists, unpaired_seqs=None, iparts=None, check_partitions=False, true_partitions=None, input_cpaths=None, input_antn_lists=None, seed_unique_ids=None, overmerge=False, naive_hamming_bound_type=None, debug=False):  # NOTE the clusters in the resulting partition generally have the uids in a totally different order to in either of the original partitions
    dbgids = None #['1437084736471665213-igh']  # None
    # ----------------------------------------------------------------------------------------
    def akey(klist):
        return ':'.join(klist)
    # ----------------------------------------------------------------------------------------
    def any_in_common(l1, l2):  # true if any uids in any cluster in l1 are found in any clusters in l2
        for tclust in l1:
            tset = set(tclust)
            if any(len(tset & set(tc)) > 0 for tc in l2):
                return True
        return False
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
            print '   %s one cluster size: %3d  %s' % (utils.color('blue', 'syncing'), len(single_cluster), ':'.join(single_cluster))
            jstrs = ['           %s %3d  %s' % ('vs %2d with sizes:'%len(cluster_list) if i==0 else '                 ', len(c), ':'.join(c)) for i, c in enumerate(cluster_list)]
            print '\n'.join(jstrs)
            print '     split into %d cdr3 group%s' % (len(cdr3_groups), utils.plural(len(cdr3_groups)))
        if not overmerge:  # default/normal
            assert naive_hamming_bound_type is not None  # ugly, but i want to only have the default set in bin/partis
            _, hi_hbound = utils.get_naive_hamming_bounds(naive_hamming_bound_type, overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))
        else:  # don't do any naive hamming splitting (correct only for --n-final-clusters 1)
            hi_hbound = 1.
        return_clusts = []
        for icdr, cdrgroup in enumerate(cdr3_groups):  # within each cdr3 group, split (i.e. use the cluster boundaries from cluster_list rather than single_cluster) if naive hfrac is > hi_hbound (but then there's shenanigans to adjudicate between different possibilities)
            if tdbg: print '      %s' % utils.color('purple', 'icdr %d' % icdr)

            # first figure out who needs to be split from whom
            clusters_to_split = {akey(c) : [] for c in cdrgroup}  # map from each cluster ('s key) to a list of clusters from which it should be split
            for c1, c2 in itertools.combinations(cdrgroup, 2):  # we could take account of the hfrac of both chains at this point, but looking at only the "split" one rather than the "merged" one, as we do here, is i think equivalent to assuming the merged one has zero hfrac, which is probably fine, since we only split if the split chain is very strongly suggesting we split
                hfrac = utils.hamming_fraction(adict[akey(c1)]['naive_seq'], adict[akey(c2)]['naive_seq'], align_if_necessary=True)  # all clusters with the same cdr3 len have been padded in waterer so their naive seqs are the same length
                if hfrac > hi_hbound:
                    clusters_to_split[akey(c1)].append(c2)
                    clusters_to_split[akey(c2)].append(c1)
                    if tdbg: print '         hfrac split %.3f > %.3f  %3d %3d  %s   %s' % (hfrac, hi_hbound, len(c1), len(c2), ':'.join(c1), ':'.join(c2))

            # then do the splitting, which is accomplished by merging each cluster in <cdrgroup> with every other cluster in <cdrgroup> from which we aren't supposed to split it (i.e. that aren't in its <clusters_to_split>)
            if tdbg:
                print '                  N to     new'
                print '          size    split   cluster?'
            tmpclusts_for_return = []  # final (return) clusters for this cdr3 class
            for cclust in cdrgroup:
                split_clusts = clusters_to_split[akey(cclust)]
                if tdbg: print '         %4d    %3d' % (len(cclust), len(split_clusts)),
                found_one = False
                for rclust in tmpclusts_for_return:  # look for an existing return cluster to which we can merge cclust, i.e. that doesn't have any uids from which we want to split
                    if any_in_common([rclust], split_clusts):  # if any uid in rclust is in a cluster from which we want to be split, skip it, i.e. don't merge with that cluster (note that we have to do it by uid because the rclusts are already merged so don't necessarily correspond to any existing cluster)
                        continue
                    # if found_one: print 'it happened!'  # can't happen any more since I switched to 'break' (although see note below)
                    if tdbg: print '     merge with size %3d   %s' % (len(rclust), ':'.join(cclust))
                    rclust += cclust
                    found_one = True
                    break  # i.e. we just merge with the first one we find and stop looking; if there's more than one, it means we could merge all three together if we wanted (triangle inequality-ish, see diagram linked at top of fcn), but i doubt it'll matter either way, and this is easier
                if not found_one:
                    if tdbg: print '      y                    %s' % ':'.join(cclust)
                    tmpclusts_for_return.append(cclust)  # if we didn't find an existing cluster that we can add it to, add it as a new cluster

            return_clusts += tmpclusts_for_return

        if tdbg:
            print '      resolved clusters:'
            for tclust in return_clusts:
                print '          %s' % ':'.join(tclust)
        return return_clusts

    # ----------------------------------------------------------------------------------------
    init_partitions = {}
    for tch in utils.chains:
        if iparts is None or ploci[tch] not in iparts:
            init_partitions[tch] = cpaths[ploci[tch]].best()  # <cpaths> (and thus <init_partitions>) are after the badly paired seqs were removed, while <input_cpaths> are the real initial ones (before anything was removed)
        else:
            init_partitions[tch] = cpaths[ploci[tch]].partitions[iparts[ploci[tch]]]
            print '  %s using non-best partition index %d for %s (best is %d)' % (utils.color('red', 'note'), iparts[ploci[tch]], tch, cpaths[ploci[tch]].i_best)

    l_translations = translate_paired_uids(ploci, init_partitions, antn_lists)
    if debug:
        for tstr, tpart in [('heavy', init_partitions['h']), ('light', init_partitions['l'])]:
            ptnprint(tpart, extrastr=utils.color('blue', '%s  '%tstr), print_partition_indices=True, n_to_print=1, sort_by_size=False, print_header=tstr=='heavy')

    common_uids, _, _ = utils.check_intersection_and_complement(init_partitions['h'], init_partitions['l'], only_warn=True, a_label='heavy', b_label='light', debug=True)  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
    if len(common_uids) == 0:
        if all(len(init_partitions[c]) > 0 for c in 'hl'):
            print '  %s no uids in common between heavy (%d uids) and light (%d uids) partitions' % (utils.color('yellow', 'warning'), len(init_partitions['h']), len(init_partitions['l']))
        return {ploci[ch] : [] for ch in init_partitions}

    antn_dict = {ch : utils.get_annotation_dict(antn_lists[ploci[ch]]) for ch in ploci}

    final_partition = []
    if debug:
        hdbg = ['    N        N       hclusts     lclusts       h/l',
                '  hclusts  lclusts    sizes       sizes      overlaps']
    # For each single cluster in each partition, get a list of the clusters in the other partition that have common uids
    # Pass this cluster + list to a fcn to resolve discrepancies by splitting on the cluster boundaries in <cluster_list> that we're sure of (i.e. that have different cdr3, or very different naive hamming fraction)
    for h_initclust, l_initclust in [(c, None) for c in init_partitions['h']] + [(None, c) for c in init_partitions['l']]:  # just loops over each single cluster in h and l partitions, but in a way that we know whether the single cluster is from h or l
        single_chain, list_chain = 'h' if l_initclust is None else 'l', 'l' if l_initclust is None else 'h'
        single_cluster = h_initclust if single_chain == 'h' else l_initclust
        cluster_list = common_clusters(single_cluster, init_partitions[list_chain])
        single_annotation = antn_dict[single_chain][akey(single_cluster)]
        annotation_list = [antn_dict[list_chain][akey(c)] for c in cluster_list]

        if debug and len(cluster_list) > 1:
            print '\n'.join(hdbg)
            hclusts, lclusts = ([single_cluster], cluster_list) if single_chain == 'h' else (cluster_list, [single_cluster])
            overlaps = [[len(set(hc) & set(lc)) for lc in lclusts] for hc in hclusts]
            overlapstr = '   '.join([' '.join(str(ov) for ov in ovlist) for ovlist in overlaps])
            def getcstr(clist): return ' '.join(str(len(c)) for c in clist)
            hcstr, lcstr = getcstr(hclusts), getcstr(lclusts)
            cw = 10
            if len(hcstr) < cw and len(lcstr) < cw:  # fits on a single line
                print ('    %2d      %2d         %-'+str(cw)+'s  %-'+str(cw)+'s  %s') % (len(hclusts), len(lclusts), hcstr, lcstr, overlapstr)
            else:  # split the last few columns over multiple lines
                print ('    %2d      %2d         %-s') % (len(hclusts), len(lclusts), hcstr)
                print ('    %2s      %2s         %-'+str(cw)+'s%-s') % ('', '', '', lcstr)
                print ('    %2s      %2s         %-'+str(cw)+'s%-'+str(cw)+'s   %s') % ('', '', '', '', overlapstr)

        resolved_clusters = resolve_discordant_clusters(copy.deepcopy(single_cluster), single_annotation, copy.deepcopy(cluster_list), annotation_list, tdbg=debug)
        if check_partitions:
            assert is_clean_partition(resolved_clusters)

        if debug:
            dbgheader = ['    adding %d resolved cluster%s to %d clusters in final partition' % (len(resolved_clusters), utils.plural(len(resolved_clusters)), len(final_partition)), '      ifclust  fset   rset   common  after: fset  rset',]
            # ----------------------------------------------------------------------------------------
            def appdbg(new_fset, rset, common_uids, xdbg=None):
                def cstr(tclust): return '(empty)' if len(tclust)==0 else ':'.join(utils.color('light_blue_bkg' if u in common_uids else None, u) for u in tclust)
                dbgstr.append(' %3d %s%3d %s  %3d    %s   %s' % (len(new_fset | common_uids), ('%-s'%utils.color('red', str(len(new_fset)), width=3, padside='right')) if len(common_uids&new_fset)==0 else '   ',
                                                                 len(rset | common_uids), ('%-s'%utils.color('red', str(len(rset)), width=3, padside='right')) if len(common_uids&rset)==0 else '   ',
                                                                 len(common_uids), cstr(new_fset), cstr(rset)))
                if xdbg is not None:
                    dbgstr.append(xdbg)
            # ----------------------------------------------------------------------------------------
            def prdbg(ifclust, dbgheader, dbgstr):
                if len(dbgheader) > 0:
                    print '\n'.join(dbgheader)
                    dbgheader = []
                print '       %4d  %s' % (ifclust, '\n             '.join(dbgstr))
                return dbgheader
        # ----------------------------------------------------------------------------------------
        n_clean = 0
        # for each cluster that's already in <final_partition> that has uids in common with a cluster in <resolved_clusters>, decide how to apportion the common uids (basically we remove them from the larger of the two clusters)
        for ifclust in range(len(final_partition)):  # iteration/<ifclust> won't get as far as any clusters that we're just adding (to the end of <final_partition>), which is what we want
            fclust = final_partition[ifclust]
            if not any_in_common([fclust], resolved_clusters):  # this is probably faster than combining it with getting the common cluster indices below, but maybe not
                n_clean += 1
                continue
            irclusts = common_clusters(fclust, resolved_clusters, return_indices=True)  # indices of any resolved_clusters that overlap with this fclust
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
                    resolved_clusters.append(list(common_uids))  # this adds a cluster at the end, which of course gets ignored in this loop over irclusts, but will get considered in the next fclust
                    if debug: xdbg = '                  %s  %s' % (utils.color('red', '+%-3d'%len(resolved_clusters[-1])), ':'.join(resolved_clusters[-1]))
                resolved_clusters[irclust] = list(rset)  # replace this resolved cluster with a copy of itself that may have had any common uids removed (if it was bigger than fclust)
                if debug: appdbg(new_fset, rset, common_uids, xdbg)
            final_partition[ifclust] = list(new_fset)  # replace <fclust> (even if nothing was removed, which shuffles the order of unchanged clusters, but oh well)
            if debug: dbgheader = prdbg(ifclust, dbgheader, dbgstr)
        assert is_clean_partition(resolved_clusters)
        final_partition += resolved_clusters  # add the (potentially modified) resolved clusters

    if debug:
        print '    removing %d/%d empty clusters' % (final_partition.count([]), len(final_partition))
    final_partition = [c for c in final_partition if len(c) > 0]
    if debug:
        print '    final: %s' % utils.cluster_size_str(final_partition)
    tmpstrs = ['   N clusters without bad/unpaired seqs:'] \
              + ['%s %4d --> %-4d%s'  % (utils.locstr(ploci[tch]), len(init_partitions[tch]), len(final_partition), chstr(len(init_partitions[tch]), len(final_partition))) for tch in utils.chains]
    print '\n        '.join(tmpstrs)

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
        untranslate_pids(ploci, init_partitions, antn_lists, l_translations, joint_partitions, antn_dict)

    if unpaired_seqs is not None:  # it might be cleaner to have this elsewhere, but I want it to happen before we evaluate, and it's also nice to have evaluation in here
        n_added = {tch : 0 for tch in ploci}
        for tch, ltmp in ploci.items():
            for upid, nearid in unpaired_seqs[ltmp].items():  # <upid> is uid of seq with bad/no pair info, <nearid> is uid of nearest seq in <upid>'s original family
                if nearid is None:  # it was a singleton, so keep it one
                    joint_partitions[tch].append([upid])
                    n_added[tch] += 1
                    continue
                jclusts = [c for c in joint_partitions[tch] if nearid in c]
                if len(jclusts) < 1:
                    joint_partitions[tch].append([upid])  # this should mean that its <nearid> was also missing pairing info, so will also be in <unpaired_seqs>, so add it as a singleton, and then when the <nearid> comes up it should get added to this new cluster
                    n_added[tch] += 1
                    continue
                jclust = utils.get_single_entry(jclusts)  # if <nearid> is in more than one cluster in the partition, it isn't a partition (which I think will only happen if it's an unfinished/uncleaned seed unique id partition)
                jclust.append(upid)
                n_added[tch] += 1
        totstr = '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in joint_partitions[tch])) for tch in sorted(ploci))
        print '    re-added unpaired seqs (%s) to give total seqs in joint partitions: %s' % (', '.join('%s %d'%(utils.locstr(ploci[tch]), n) for tch, n in n_added.items()), totstr)

    if true_partitions is not None:
        assert iparts is None  # just for now
        evaluate_joint_partitions(ploci, true_partitions, {tch : input_cpaths[ploci[tch]].best() for tch in utils.chains}, joint_partitions, antn_lists, seed_unique_ids=seed_unique_ids, debug=debug)

    return {ploci[ch] : jp for ch, jp in joint_partitions.items()}
