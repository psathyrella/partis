import copy
import itertools
import numpy
import sys
import operator
import string
import os
import math

import utils
import prutils
from clusterpath import ptnprint, ClusterPath
from hist import Hist

naive_hamming_bound_type = 'naive-hamming' #'likelihood'

# ----------------------------------------------------------------------------------------
# return standardized file name (including subdirs) in directory structure that we use for paired heavy/light i/o
def paired_fn(bdir, locus, lpair=None, suffix='.fa', ig_or_tr='ig'):  # if set, only file(s) for this <locus>, and/or only files for this <lpair> of loci. If <lpair> is set but <locus> is None, returns subdir name
    if lpair is not None:
        bdir = '%s/%s' % (bdir, '+'.join(lpair))
        if locus is None:
            return bdir
    return '%s/%s%s' % (bdir, locus, suffix)

# ----------------------------------------------------------------------------------------
def paired_dir_fnames(bdir, no_pairing_info=False, only_paired=False, suffix='.fa', ig_or_tr='ig', include_failed=False, include_meta=False):  # return all files + dirs from previous fcn
    fnames = []
    if not only_paired:
        fnames += [paired_fn(bdir, l, suffix=suffix) for l in utils.sub_loci(ig_or_tr)]  # single-locus files
    if include_failed:
        fnames += [paired_fn(bdir, 'failed', suffix=suffix)]  # kind of hackey, but bin/split-loci.py does this, so we need to be able to clean it up
    if include_meta:
        fnames += [paired_fn(bdir, 'meta', suffix='.yaml')]  # this is also kind of hackey
    if not no_pairing_info:
        fnames += [paired_fn(bdir, l, lpair=lp, suffix=suffix) for lp in utils.locus_pairs[ig_or_tr] for l in lp]  # paired single-locus files
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
    # then go through the light chain annotations + partition swapping names
    l_translations = {}
    for lline in antn_lists[ploci['l']]:
        for iseq, l_id in enumerate(lline['unique_ids']):
            if l_id not in h_paired_uids:
                raise Exception('no paired uids for %s' % l_id)  # everybody has to have exactly one paired id at this point
            lline['unique_ids'][iseq] = h_paired_uids[l_id]
            l_translations[h_paired_uids[l_id]] = l_id  # so we can go back to <l_id> afterwards
    if len(h_paired_uids) > 0:
        init_partitions['l'] = [[h_paired_uids.get(u, u) for u in c] for c in init_partitions['l']]
    return l_translations

# ----------------------------------------------------------------------------------------
# reverse action of previous fcn
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
def remove_badly_paired_seqs(ploci, cpaths, antn_lists, glfos, debug=False):  # remove seqs paired with the wrong light chain, as well as those with no pairing info (the latter we keep track of so we can insert them later into the right final cluster)
    # ----------------------------------------------------------------------------------------
    def add_unpaired(cline, iseq, uid):
        sorted_hdists = sorted([(u, utils.hamming_distance(cline['seqs'][i], cline['seqs'][iseq])) for i, u in enumerate(cline['unique_ids']) if i != iseq], key=operator.itemgetter(1))
        nearest_uid = sorted_hdists[0][0] if len(sorted_hdists) > 0 else None
        unpaired_seqs[cline['loci'][iseq]][uid] = nearest_uid  # unless there's no other seqs in the cluster, attach it to the nearest seq by hamming distance
    # ----------------------------------------------------------------------------------------
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l]) for l in antn_lists}
    all_loci = {u : l for ants in antn_lists.values() for antn in ants for u, l in zip(antn['unique_ids'], antn['loci'])}  # this includes the heavy ones, which we don't need, but oh well
    all_pids = {u : pids[0] for alist in antn_lists.values() for l in alist for u, pids in zip(l['unique_ids'], l['paired-uids']) if len(pids)==1}  # I'm pretty sure that the partition implied by the annotations is identical to the one in <cpaths>, and it's nice to loop over annotations for this
    unpaired_seqs = {l : {} for l in ploci.values()}  # map for each locus from the uid of each seq with no (or non-reciprocal) pairing info to the nearest sequence in its family (after merging partitions we'll insert it into the family that this nearest seq ended up in)
    lp_cpaths, lp_antn_lists = {}, {}
    if debug:
        print '  removing bad/un-paired seqs'
        print '          N       N      no   wrong  non-'
        print '        before removed  info  light recip'
    for tch in sorted(ploci):
        new_partition, new_antn_list = [], []
        for iclust, cluster in enumerate(cpaths[ploci[tch]].best()):
            cline = antn_dicts[ploci[tch]][':'.join(cluster)]
            iseqs_to_remove = []
            n_no_info, n_bad_paired, n_non_reciprocal = 0, 0, 0  # just for dbg NOTE n_bad_paired are the only ones we *really* want to remove, since we know they're paired with the long light chain, whereas the other two categories we eventually want to re-add since we're not sure who they're paired with
            for iseq, uid in enumerate(cline['unique_ids']):
                pids = cline['paired-uids'][iseq]
                if len(pids) == 0:  # no pairing info
                    iseqs_to_remove.append(iseq)
                    add_unpaired(cline, iseq, uid)
                    n_no_info += 1
                elif len(pids) > 1:
                    raise Exception('multiple paired uids for \'%s\': %s' % (uid, pids))
                else:
                    if tch == 'h' and all_loci[utils.get_single_entry(pids)] != ploci['l']:  # if it's the other light chain
                        iseqs_to_remove.append(iseq)
                        n_bad_paired += 1
                    else:
                        if all_pids[uid] not in all_pids or all_pids[all_pids[uid]] != uid:  # also remove any non-reciprocal pairings (I think this will still miss any whose partner was removed) NOTE it would be nice to enforce reciprocal pairings in clean_pair_info(), but atm i think we can't look at both chains at once in that fcn
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
                print '    %s   %3d    %3s     %3s   %3s    %3s' % (utils.locstr(ploci[tch]) if iclust==0 else ' ', len(cline['unique_ids']), fstr(len(iseqs_to_remove)), fstr(n_no_info), fstr(n_bad_paired), fstr(n_non_reciprocal))
        lp_cpaths[ploci[tch]] = ClusterPath(seed_unique_id=cpaths[ploci[tch]].seed_unique_id, partition=new_partition)
        lp_antn_lists[ploci[tch]] = new_antn_list

    if debug:
        print '    totals before: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in cpaths[ploci[tch]].best())) for tch in sorted(ploci))
        print '    totals after: %s' % '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in lp_cpaths[ploci[tch]].best())) for tch in sorted(ploci))

    return lp_cpaths, lp_antn_lists, unpaired_seqs

# ----------------------------------------------------------------------------------------
def clean_pair_info(cpaths, antn_lists, max_hdist=4, is_data=False, plotdir=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def check_droplet_id_groups(all_uids, tdbg=False):
        try:
            utils.get_droplet_id(next(iter(all_uids)))
        except:
            print '  note: couldn\'t get droplet id from \'%s\', so assuming this isn\'t 10x data' % next(iter(all_uids))  # NOTE i'm not sure that this gives the same one as the previous line
            return False
        # check against the droplet id method (we could just do it this way, but it would only work for 10x, and only until they change their naming convention)
        pgroup_strs = set(':'.join(sorted(pg)) for pg in pid_groups)
        n_not_found = 0
        for dropid, drop_queries in itertools.groupby(sorted(all_uids, key=utils.get_droplet_id), key=utils.get_droplet_id):
            dqlist = list(drop_queries)
            found = ':'.join(sorted(dqlist)) in pgroup_strs
            if not found:
                overlaps = [g for g in pgroup_strs if dropid in g]
                overlaps = utils.get_single_entry(overlaps)
                n_not_found += 1
            if tdbg or not found:
                print '  %25s %s  %s  %s' % (utils.color('green', '-') if found else utils.color('red', 'x'), dropid, ' '.join(sorted(utils.get_contig_id(q) for q in dqlist)), utils.color('red', ' '.join(sorted(utils.get_contig_id(q) for q in overlaps.split(':'))) if not found else ''))
        if n_not_found > 0:
            print '  %s droplet id group check failed for %d groups' % (utils.color('red', 'error'), n_not_found)
        return True
    # ----------------------------------------------------------------------------------------
    def plot_uids_before(plotdir, pid_groups, all_antns):
        # ----------------------------------------------------------------------------------------
        def fnfplot(logstr, fhists, n_max_bins=15):
            import plotting
            fklabels = {'func' : 'all func.', 'nonfunc' : 'any non.'}
            fig, ax = plotting.mpl_init()
            for fk, fcolor in zip(fhists, plotting.default_colors):
                fhists[fk].mpl_plot(ax, label=fklabels[fk], color=fcolor)
                if logstr == '':
                    fhists[fk].write('%s/%s.csv'%(plotdir, fk + '-per-drop'))
                xticks = fhists[fk].get_bin_centers()
                xticklabels = fhists[fk].bin_labels
                xbounds = None if fhists[fk].n_bins<n_max_bins else (0, n_max_bins)
            plotting.mpl_finish(ax, plotdir, 'func-non-func-per-drop'+logstr, xlabel='N seqs per droplet', ylabel='counts', title='before', log='' if logstr=='' else 'y', leg_loc=(0.6, 0.7), xticks=xticks, xticklabels=xticklabels, xbounds=xbounds)
        # ----------------------------------------------------------------------------------------
        bhist = Hist(value_list=[len(pg) for pg in pid_groups], init_int_bins=True)
        bhist.fullplot(plotdir, 'seqs-per-droplet', fargs={'xlabel' : 'seqs per droplet', 'ylabel' : 'counts', 'title' : 'before'})
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
        ahist.fullplot(plotdir, 'paired-seqs-per-seq-%s'%pstr, pargs={'remove_empty_bins' : True}, fargs={'xlabel' : 'N paired seqs per seq', 'ylabel' : 'counts', 'title' : pstr, 'xbounds' : (-0.05, 1.05*ahist.xmax), 'xticks' : [i for i in range(0, int(ahist.xmax+1))]})
        return pidlengths
    # ----------------------------------------------------------------------------------------
    def make_final_plots(initial_seqs_per_seq, initial_flcounts):
        final_seqs_per_seq = plot_n_pseqs_per_seq('after')
        import plotting
        plotting.plot_smatrix(plotdir, 'pseq-matrix', final_seqs_per_seq, initial_seqs_per_seq, n_max_bins=12, xlabel='after', ylabel='before', lfcn=lambda x: 'miss.' if x==-1 else str(x), title='N paired seqs per seq')
        final_flcounts = {}  # note that this has to be per seq (even though that kind of double counts) since otherwise we wouldn't have a way to determine correspondence between initial and final
        for ltmp in sorted(cpaths):
            for cluster in cpaths[ltmp].best():
                atn = antn_dicts[ltmp][':'.join(cluster)]
                final_flcounts.update({u : lgstr(set([u] + pids), for_plot=True) for u, pids in zip(atn['unique_ids'], atn['paired-uids'])})  # have to make sure <u> is included in <pids> (as well as that there's no duplicates)
        plotting.plot_smatrix(plotdir, 'flcount-matrix', final_flcounts, initial_flcounts, kfcn=len, n_max_bins=15,
                              lfcn=lambda x: 'miss.' if x==-1 else ('none' if x=='' else str(x)), xlabel='after', ylabel='before', title='pair combo (per seq)', tdbg=2 if debug else False)
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
        return ' '.join(lgstrs)
    # ----------------------------------------------------------------------------------------
    def choose_seqs_to_remove(chain_ids, tdbg=False):  # choose one of <chain_ids> to eliminate (based on identical/similar seq collapse and productivity)
        ids_to_remove = set(u for u in chain_ids if getloc(u)=='?')  # remove any with missing annotations
        if tdbg and len(ids_to_remove) > 0:  # i think this actually can't happen a.t.m.
            print '      removed %d with missing annotations' % len(ids_to_remove)

        # among any pairs of sequences that are [almost] identical at all non-ambiguous position, keep only the longest one
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

        # remove unproductive (only on real data, since simulation usually has lots of stop codons)
        dbgstr = []
        unproductive_ids = []
        for uid in chain_ids:
            if is_data and not utils.is_functional(all_antns[uid], all_antns[uid]['unique_ids'].index(uid)):
                unproductive_ids.append(uid)
                if tdbg:
                    dbgstr.append(utils.is_functional_dbg_str(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), sep='+'))
        # unproductive_ids = [u for u in chain_ids if not utils.is_functional(all_antns[u], all_antns[u]['unique_ids'].index(u))]  # this way is only one line, which may or may not be nicer
        if tdbg and len(unproductive_ids) > 0:
            print '        %d unproductive  %s' % (len(unproductive_ids), ',  '.join(dbgstr))
            ids_to_remove |= set(unproductive_ids)

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
            new_pfamilies = get_pfamily_dict(cline, extra_str='after:')
            pfcounts = [[new_pfamilies[pfkey(p)]['count'] for p in pids] for pids in cline['paired-uids']]
            pfids = [[new_pfamilies[pfkey(p)]['id'] for p in pids] for pids in cline['paired-uids']]
            uid_extra_strs = ['%s: %s'%(utils.locstr(l), lcstr(pids, pfcs, pfids)) for l, pids, pfcs, pfids in zip(cline['loci'], cline['paired-uids'], pfcounts, pfids)]
            old_pfcounts = [[pfamilies[pfkey(p)]['count'] for p in pids] for pids in old_pids]
            old_estrs = ['%s: %s'%(utils.locstr(l), lcstr(pids, pfcs)) for l, pids, pfcs in zip(cline['loci'], old_pids, old_pfcounts)]
            for istr, (oldstr, newstr) in enumerate(zip(old_estrs, uid_extra_strs)):
                if newstr != oldstr:
                    uid_extra_strs[istr] = '%s%s (%s)' % (newstr, ' '*(12 - utils.len_excluding_colors(newstr)), oldstr)
            utils.print_reco_event(cline, uid_extra_strs=uid_extra_strs, extra_str='      ')
        # ----------------------------------------------------------------------------------------
        cline = antn_dicts[ltmp][':'.join(cluster)]
        if any(u not in pid_groups[pid_ids[u]] for u in cline['unique_ids']):  # shouldn't be able to happen any more, but that was a really bad/dumb bug
            raise Exception('one of unique ids %s not in its own pid group' % cline['unique_ids'])
        cline['paired-uids'] = [[p for p in pid_groups[pid_ids[u]] if p != u] for u in cline['unique_ids']]

        pfamilies = get_pfamily_dict(cline, extra_str='before:')

        def pfkey(p): return ':'.join(all_antns[p]['unique_ids'])  # keystr for the family of paired id <p>
        old_pids = copy.deepcopy(cline['paired-uids'])  # just for dbg

        # for each uid, choose the pid that's of opposite chain, and has the most other uids voting for it
        for iseq, (uid, pids) in enumerate(zip(cline['unique_ids'], cline['paired-uids'])):
            pid_to_keep = None
            ochain_pidfcs = [(p, pfamilies[pfkey(p)]['count']) for p in pids if not utils.samechain(getloc(p), getloc(uid))]
            if len(ochain_pidfcs) > 0:
                sorted_pids, sorted_pfcs = zip(*sorted(ochain_pidfcs, key=operator.itemgetter(1), reverse=True))
                # note that even if there's only one ochain choice, there can be other same-chain ones that we still want to drop (hence the <2 below)
                if len(sorted_pfcs) < 2 or sorted_pfcs[0] > sorted_pfcs[1] or pfkey(sorted_pids[0]) == pfkey(sorted_pids[1]):  # in order to drop the later ones, the first one either has to have more counts, or at least the second one has to be from the same family (in the latter case we still don't know which is the right one, but for the purposes of clustering resolution we just need to know what the family is)
                    pid_to_keep = sorted_pids[0]
            cline['paired-uids'][iseq] = [pid_to_keep] if pid_to_keep is not None else []  # if we didn't decide on an opposite-chain pid, remove all pairing info

        if debug: # and len(cline['unique_ids']) > 1:  # NOTE it's annoying printing the singletons, but it's way worse when they're just missing and you can't figure out where a sequence went
            print_dbg()

    # ----------------------------------------------------------------------------------------
    def get_pfamily_dict(cline, extra_str=''):  # see what others in its family are paired with
        pfamilies = {}  # counts how many different uids in <cline> had paired ids from each potential paired family (well, how many pids total, but i think they'll always be the same)
        fid_counter = 0  # give a unique id to each family, for dbg printing purposes
        for uid, pids in zip(cline['unique_ids'], cline['paired-uids']):
            for pid in pids:
                fline = all_antns[pid]  # family for this paired id
                fkey = ':'.join(fline['unique_ids'])
                if fkey not in pfamilies:
                    pfamilies[fkey] = {'locus' : gval(pid, 'loci'), 'count' : 0, 'id' : fid_counter}
                    fid_counter += 1
                pfamilies[fkey]['count'] += 1
        if debug and len(cline['unique_ids']) > 1:
            print '    %6s N  size  id  cdr3' % extra_str
            for fkey, fdict in sorted(pfamilies.items(), key=lambda x: x[1]['count'], reverse=True):
                print '       %s %3d  %3d   %2d %s  %3d' % (utils.locstr(fdict['locus']), fdict['count'], len(antn_dicts[fdict['locus']][fkey]['unique_ids']), fdict['id'], fidstr(fdict['id']), antn_dicts[fdict['locus']][fkey]['cdr3_length'])

        return pfamilies

    # ----------------------------------------------------------------------------------------
    antn_dicts = {l : utils.get_annotation_dict(antn_lists[l]) for l in antn_lists}
    all_uids = set(u for p in cpaths.values() for c in p.best() for u in c)  # all uids that occur in a partition (should I think be the same as the ones for which we have valid/non-failed annotations)

    # first make a map from each uid (for all loci) to its annotation
    pid_groups = []  # list of pid groups, i.e. each element is the uids from a single droplet (for 10x)
    pid_ids = {}  # map from each uid to the index of its pid group
    all_antns = {}
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
                missing_ids = set(pids) - all_uids
                n_missing += len(missing_ids)
                pset = set([uid] + pids) - missing_ids
                found = False
                for ipg, pgroup in enumerate(pid_groups):
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
        print '   %d missing uids' % n_missing # NOTE at least for now we're skipping invalid queries when reading output
    # for ipg, pg in enumerate(pid_groups):
    #     print '  %3d %s' % (ipg, ' '.join(pg))

    idg_ok = check_droplet_id_groups(all_uids)  # NOTE not using the return value here, but I may need to in the future
    if plotdir is not None:
        initial_flcounts = plot_uids_before(plotdir, pid_groups, all_antns)
        initial_seqs_per_seq = plot_n_pseqs_per_seq('before')

    # then go through each group trying to remove as many crappy/suspicously similar ones as possible (this step has a fairly minor effect compared to the partition-based step below)
    if debug:
        print '  cleaning %d pid groups:' % len(pid_groups)
        def tmpincr(pgroup, cdict):
            if lgstr(pgroup) not in cdict:
                cdict[lgstr(pgroup)] = 0
            cdict[lgstr(pgroup)] += 1
    ok_groups, tried_to_fix_groups = {}, {}
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
            ids_to_remove = choose_seqs_to_remove(idlist)
            for rid in ids_to_remove:
                pgroup.remove(rid)
                idlist.remove(rid)
                pid_groups.append(set([rid]))  # add the removed id to a new pid group of its own at the end (so it'll show up as unpaired)
                pid_ids[rid] = len(pid_groups) - 1
            if debug > 1:
                print '      %s: removed %d, leaving %d' % (utils.color('green', 'fixed') if len(idlist)==1 else utils.color('red', 'nope'), len(ids_to_remove), len(idlist))
                if len(idlist) > 1:
                    for uid in idlist:
                        prutils.print_seq_in_reco_event(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), one_line=True, extra_str='        ', uid_extra_str=utils.locstr(getloc(uid)))

        pid_groups[ipg] = pgroup
        if debug: tmpincr(pgroup, tried_to_fix_groups)

    if debug:
        print '    ok to start with:'
        for lstr, count in sorted(ok_groups.items(), key=operator.itemgetter(1), reverse=True):
            print '      %3d  %s' % (count, lstr)
        print '    tried to fix (afterwards):'
        for lstr, count in sorted(tried_to_fix_groups.items(), key=operator.itemgetter(1), reverse=True):
            print '      %3d  %s' % (count, lstr)

    # then go through again using cluster/family information to try to cut everyone down to one paired id (and if we can't get down to one, we remove all of them) NOTE also re-sets the actual 'paired-uids' keys
    for ltmp in sorted(cpaths):
        if debug:
            print '%s' % utils.color('green', ltmp)
            name_dict, name_ids = {'potential' : None, 'used' : None}, {}
        for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
            clean_with_partition_info(cluster)

    if plotdir is not None:
        make_final_plots(initial_seqs_per_seq, initial_flcounts)

# ----------------------------------------------------------------------------------------
def compare_partition_pair(cfpart, refpart, remove_from_ref=False, antn_list=None, dbg_str=None, cf_label='inferred', ref_label='true', debug=False):
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
    if remove_from_ref:  # we could also use add_missing_uids_to_partition()
        refpart = utils.remove_missing_uids_from_ref_partition(refpart, cfpart, debug=debug)  # returns a new/copied partition, doesn't modify original
    return utils.per_seq_correct_cluster_fractions(cfpart, refpart, dbg_str=dbg_str, inf_label=cf_label, true_label=ref_label, debug=debug)
    # TODO figure out which cases of 'missing' uids should really be removed, and which should be singletons

# ----------------------------------------------------------------------------------------
def evaluate_joint_partitions(ploci, true_partitions, init_partitions, joint_partitions, antn_lists, debug=False):
    # NOTE that <joint_partitions> can have many fewer seqs than <init_partitions> since in making <joint_partitions> we remove seqs paired with the other light chain (the weighted average ccfs over the h joint partitions corresponding to both light chains would be exactly comparable to the <init_partitions>, but I think this is fine as it is)
    for tch in utils.chains:
        ltmp = ploci[tch]
        ccfs = {}
        for dstr, cfpart in [('single', init_partitions[tch]), ('joint', joint_partitions[tch])]:
            ccfs[dstr] = compare_partition_pair(cfpart, true_partitions[ltmp], remove_from_ref=True,  # removes from the true ptn any uids that are missing from the inferred ptn
                                                antn_list=antn_lists[ltmp], dbg_str='%s %s '%(utils.locstr(ltmp), dstr), debug=debug)
        print '  %s ccfs:     purity  completeness' % utils.locstr(ltmp)
        print '      single  %6.3f %6.3f' % (ccfs['single'][0], ccfs['single'][1])
        print '       joint  %6.3f %6.3f' % (ccfs['joint'][0], ccfs['joint'][1])

# ----------------------------------------------------------------------------------------
# cartoon explaining algorithm here https://github.com/psathyrella/partis/commit/ede140d76ff47383e0478c25fae8a9a9fa129afa#commitcomment-40981229
def merge_chains(ploci, cpaths, antn_lists, unpaired_seqs=None, iparts=None, check_partitions=False, true_partitions=None, input_cpaths=None, input_antn_lists=None, debug=False):  # NOTE the clusters in the resulting partition generally have the uids in a totally different order to in either of the original partitions
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
        lo_hbound, hi_hbound = utils.get_naive_hamming_bounds(naive_hamming_bound_type, overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))  # these are the wider bounds, so < lo is almost certainly clonal, > hi is almost certainly not
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
        print '  %s no uids in common between heavy and light' % utils.color('yellow', 'warning')

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
        print '    final: %s' % ' '.join([str(len(c)) for c in final_partition])
    def chstr(n_before, n_after):
        if n_before == n_after: return ''
        else: return ' ' + utils.color('red', '%+d' % (n_after - n_before))
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
        if debug:
            totstr = '  '.join('%s %d'%(utils.locstr(ploci[tch]), sum(len(c) for c in joint_partitions[tch])) for tch in sorted(ploci))
            print '    re-added unpaired seqs (%s) to give total seqs in joint partitions: %s' % (', '.join('%s %d'%(utils.locstr(ploci[tch]), n) for tch, n in n_added.items()), totstr)

    if true_partitions is not None:
        assert iparts is None  # just for now
        evaluate_joint_partitions(ploci, true_partitions, {tch : input_cpaths[ploci[tch]].best() for tch in utils.chains}, joint_partitions, antn_lists, debug=debug)

    tmpstrs = ['   N clusters with all seqs:'] \
              + ['%s %4d --> %-4d%s'  % (utils.locstr(ploci[tch]), len(input_cpaths[ploci[tch]].best()), len(joint_partitions[tch]), chstr(len(input_cpaths[ploci[tch]].best()), len(joint_partitions[tch]))) for tch in utils.chains]
    print '\n        '.join(tmpstrs)
    if debug:
        for tch in utils.chains:
            input_antn_dict = utils.get_annotation_dict(input_antn_lists[ploci[tch]])
            print '%s' % utils.color('green', ploci[tch])
            assert iparts is None  # just for now
            for tclust in input_cpaths[ploci[tch]].best():  # loop over clusters in the initial partition
                jfamilies = [c for c in joint_partitions[tch] if len(set(tclust) & set(c)) > 0]  # clusters in the joint partition that overlap with this cluster
                uid_extra_strs = None
                if len(jfamilies) == 0:  # couldn't find it (it was probably paired with the other light chain)
                    uid_extra_strs = [utils.color('blue', '?') for _ in tclust]
                else:
                    def getjcstr(u):  # str for length of <u>'s cluster in <jfamilies>
                        jfs = [f for f in jfamilies if u in f]
                        return utils.color('blue', '?') if len(jfs)==0 else ('%d: %s'%(len(utils.get_single_entry(jfs)), ' '.join(jfs[0])))
                    uid_extra_strs = [getjcstr(u) for u in tclust]
                utils.print_reco_event(input_antn_dict[':'.join(tclust)], uid_extra_strs=uid_extra_strs, extra_str='      ')

    return {ploci[ch] : jp for ch, jp in joint_partitions.items()}
