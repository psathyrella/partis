import copy
import itertools
import numpy
import sys
import operator
import string

import utils
import prutils
from clusterpath import ptnprint

naive_hamming_bound_type = 'naive-hamming' #'likelihood'

# ----------------------------------------------------------------------------------------
# rename all uids in the light chain partition, and annotations that are paired with a heavy chain uid, to that heavy chain uid (pairings must, at this stage, be unique)
def translate_paired_uids(ploci, init_partitions, antn_lists):
    h_paired_uids = {}  # map ot each heavy chain uid <u> from its paired light chain uid <pids[0]>
    for hline in antn_lists[ploci['h']]:
        for h_id, pids in zip(hline['unique_ids'], hline['paired-uids']):
            if len(pids) == 0:
                continue
            elif len(pids) > 1:
                raise Exception('multiple paired uids %s for %s sequence %s' % (' '.join(pids), ploci['h'], h_id))
            h_paired_uids[pids[0]] = h_id
    l_translations = {}
    for lline in antn_lists[ploci['l']]:
        for iseq in range(len(lline['unique_ids'])):
            l_id = lline['unique_ids'][iseq]
            if l_id not in h_paired_uids:  # this <l_id> wasn't paired with any heavy chain ids
                continue
            lline['unique_ids'][iseq] = h_paired_uids[l_id]
            l_translations[h_paired_uids[l_id]] = l_id  # so we can go back to <l_id> afterwards
    if len(h_paired_uids) > 0:
        init_partitions['l'] = [[h_paired_uids.get(u, u) for u in c] for c in init_partitions['l']]
    return l_translations

# ----------------------------------------------------------------------------------------
# reverse action of previous fcn
def untranslate_pids(ploci, init_partitions, antn_lists, l_translations, joint_partitions):
    for lline in antn_lists[ploci['l']]:
        lline['unique_ids'] = [l_translations.get(u, u) for u in lline['unique_ids']]
    init_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in init_partitions['l']]
    joint_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in joint_partitions['l']]
    # also need to remove h ids from the l joint partition and vice versa
    all_loci = {u : l for ants in antn_lists.values() for antn in ants for u, l in zip(antn['unique_ids'], antn['loci'])}  # it seems like i should have this info somewhere already, but maybe not?
    for tch in joint_partitions:
        joint_partitions[tch] = [[u for u in c if all_loci[u]==ploci[tch]] for c in joint_partitions[tch]]  # i think they should all be in <all_loci>
    antn_dict = {ch : utils.get_annotation_dict(antn_lists[ploci[ch]]) for ch in ploci}  # atm this only gets used for dbg, but still nice to properly fix it

# ----------------------------------------------------------------------------------------
def clean_pair_info(cpaths, antn_lists, max_hdist=4, n_max_clusters=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def check_droplet_id_groups(all_uids, tdbg=False):
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
    def lgstr(lgroup, sort=True):
        return ' '.join(utils.locstr(l) for l in (sorted if sort else utils.pass_fcn)([getloc(u) for u in lgroup]))
    # ----------------------------------------------------------------------------------------
    def choose_seqs_to_remove(chain_ids, tdbg=False):  # choose one of <chain_ids> to eliminate
        # look for pairs with the same locus that
        ids_to_remove = set(u for u in chain_ids if getloc(u)=='?')
        if tdbg and len(ids_to_remove) > 0:  # i think this actually can't happen a.t.m. TODO maybe remove it
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
            if hdist <= max_hdist:  # TODO would be nice to be able to combine their sequences, but I think propagating the resulting annotation modifications would be hard
                better_id, worse_id = sorted(tpair, key=lambda q: utils.ambig_frac(gval(q, 'seqs')))  # TODO if we're tossing one with hdist > 0, maybe should take the lower-shm one if they're the same length?
                ids_to_remove.add(worse_id)
                n_equivalent += 1
        if tdbg and len(dbgstr) > 0:
            print '        %d pair%s equivalent with hdists %s' % (n_equivalent, utils.plural(n_equivalent), ' '.join(dbgstr))

        # remove unproductive
        dbgstr = []
        unproductive_ids = []
        for uid in chain_ids:
            if not utils.is_functional(all_antns[uid], all_antns[uid]['unique_ids'].index(uid)):
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
        print '  %s consolidating info for %d loci with family/sequence counts: %s' % (utils.color('blue', '+'.join(cpaths)), len(cpaths), '  '.join('%s: %d/%d'%(l, len(cpaths[l].best()), sum(len(c) for c in cpaths[l].best())) for l in sorted(cpaths)))
    for ltmp in sorted(cpaths):
        for cluster in cpaths[ltmp].best():
            cline = antn_dicts[ltmp][':'.join(cluster)]
            if 'paired-uids' not in cline:
                print '  %s no paired-uids in line' % utils.color('yellow', 'warning')
                continue  # maybe should still add to all_antns?
            for uid, pids in zip(cline['unique_ids'], cline['paired-uids']):
                missing_ids = set(pids) - all_uids
                n_missing += len(missing_ids)
                pset = set([uid] + pids) - missing_ids
                found = False
                for ipg, pgroup in enumerate(pid_groups):
                    if any(p in pgroup for p in pset):  # TODO should maybe check for consistency if some of them are already in there (i.e. from reciprocal info in another chain)?
                        found = True
                        pgroup |= pset
                        break
                if not found:
                    pid_groups.append(pset)
                    ipg = len(pid_groups) - 1
                assert ipg is not None
                for pid in pset:
                    pid_ids[pid] = ipg

            cline['loci'] = [ltmp for _ in cline['unique_ids']]  # TODO maybe should add this somewhere else, like in partitiondriver? (eh, maybe not? the locus is always available in each file from the germline info anyway)
            for uid in cline['unique_ids']:
                all_antns[uid] = cline
    if n_missing > 0:
        print '   %d missing uids' % n_missing # NOTE at least for now we're skipping invalid queries when reading output
    # for ipg, pg in enumerate(pid_groups):
    #     print '  %3d %s' % (ipg, ' '.join(pg))

    check_droplet_id_groups(all_uids)
    # TODO handle/keep better track of failures
    # TODO make sure the missing ids are actually really gone

    # then go through each group trying to remove as many crappy/suspicous ones as possible
    print '  cleaning %d pid groups:' % len(pid_groups)
    ok_groups, tried_to_fix_groups = {}, {}
    for ipg, pgroup in enumerate(pid_groups):
        pgroup = [u for u in pgroup if getloc(u) != '?']  # TODO figure out what to do with missing ones
        hids = [u for u in pgroup if utils.has_d_gene(getloc(u))]
        lids = [u for u in pgroup if u not in hids]
        if len(hids) < 2 and len(lids) < 2:
            if lgstr(pgroup) not in ok_groups:
                ok_groups[lgstr(pgroup)] = 0
            ok_groups[lgstr(pgroup)] += 1
            pid_groups[ipg] = pgroup
            continue
        if debug > 1:
            print '    %s' % lgstr(pgroup),
        for chain, idlist in zip(utils.chains, [hids, lids]):
            if len(idlist) < 2:
                continue
            if debug > 1:
                print '\n      too many %s chains: %s' % (chain, lgstr(idlist))
            ids_to_remove = choose_seqs_to_remove(idlist)
            for rid in ids_to_remove:
                pgroup.remove(rid)
                idlist.remove(rid)
            if debug > 1:
                print '      %s: removed %d, leaving %d' % (utils.color('green', 'fixed') if len(idlist)==1 else utils.color('red', 'nope'), len(ids_to_remove), len(idlist))
                if len(idlist) > 1:
                    for uid in idlist:
                        prutils.print_seq_in_reco_event(all_antns[uid], all_antns[uid]['unique_ids'].index(uid), one_line=True, extra_str='        ', uid_extra_str=utils.locstr(getloc(uid)))

        pid_groups[ipg] = pgroup
        if lgstr(pgroup) not in tried_to_fix_groups:
            tried_to_fix_groups[lgstr(pgroup)] = 0
        tried_to_fix_groups[lgstr(pgroup)] += 1

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
            cpaths[ltmp].print_partitions()
            name_dict, name_ids = {'potential' : None, 'used' : None}, {}
        for iclust, cluster in enumerate(sorted(cpaths[ltmp].best(), key=len, reverse=True)):
            if n_max_clusters is not None and iclust >= n_max_clusters:
                break
            cline = antn_dicts[ltmp][':'.join(cluster)]
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
                def lcstr(pids, pfcs, pfids=None):  # returns string summarizing the families of the paired uids for a uid, e.g. 'k 51  l 3  h 1' if the uid has three potential pids, one from k with which 50 other uids in <cline> are paired, etc.
                    if len(pids) == 0: return ''
                    if pfids is None:
                        spids, spfcs = zip(*sorted(zip(pids, pfcs), key=operator.itemgetter(1), reverse=True))
                        spfids = ['' for _ in spids]
                    else:
                        assert len(pids) == 1  # if we passed in <pfids> (an id for each paired family), this should be after cleaning, so there should only be one of them
                        spids, spfcs, spfids = pids, pfcs, pfids
                    return '  '.join('%s %s %s'%(lg, '%1d'%sp if pfids is None else '%2d'%sp, fidstr(fid)) for lg, sp, fid in zip(lgstr(spids, sort=False).split(' '), spfcs, spfids))
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
def evaluate_joint_partitions(ploci, true_partitions, init_partitions, joint_partitions, antn_lists):
    # ----------------------------------------------------------------------------------------
    def incorporate_duplicates(tpart, dup_dict):  # take the map from uid to list of its duplicates (dup_dict), and add the duplicates to any clusters in partition tpart that contain that uid
        for tclust in tpart:
            for uid in tclust:
                if uid in dup_dict:
                    tclust += dup_dict[uid]
    # ----------------------------------------------------------------------------------------
    cmp_partitions = {}  # (potentially) modified versions of the initial heavy/light partitions
    ccfs = {}
    for chain in utils.chains:
        cmp_partitions[chain] = copy.deepcopy(init_partitions[chain])
        true_partitions[chain] = utils.remove_missing_uids_from_true_partition(true_partitions[chain], cmp_partitions[chain], debug=False)  # NOTE it would probably be better to not modify the true partition, since it's getting passed in from outside
        dup_dict = {u : l['duplicates'][i] for l in antn_lists[ploci[chain]] for i, u in enumerate(l['unique_ids']) if len(l['duplicates'][i]) > 0}
        if len(dup_dict) > 0:
            incorporate_duplicates(cmp_partitions[chain], dup_dict)
        ccfs[chain] = {'before' : utils.new_ccfs_that_need_better_names(cmp_partitions[chain], true_partitions[chain])}

        if len(dup_dict) > 0:
            incorporate_duplicates(joint_partitions[chain], dup_dict)  # NOTE this modifies the joint partition
        j_part = utils.get_deduplicated_partitions([joint_partitions[chain]])[0]  # TODO why do i need this?
        j_part = utils.remove_missing_uids_from_true_partition(j_part, true_partitions[chain], debug=False)  # we already removed failed queries from each individual chain's partition, but then if the other chain didn't fail it'll still be in the joint partition
        ccfs[chain]['joint'] = utils.new_ccfs_that_need_better_names(j_part, true_partitions[chain])

    print '             purity  completeness'
    for chain in utils.chains:
        print '   %s before  %6.3f %6.3f' % (chain, ccfs[chain]['before'][0], ccfs[chain]['before'][1])
    for chain in utils.chains:
        print '    joint    %6.3f %6.3f   (%s true)' % (ccfs[chain]['joint'][0], ccfs[chain]['joint'][1], chain)

# ----------------------------------------------------------------------------------------
# cartoon explaining algorithm here https://github.com/psathyrella/partis/commit/ede140d76ff47383e0478c25fae8a9a9fa129afa#commitcomment-40981229
def merge_chains(ploci, cpaths, antn_lists, iparts=None, check_partitions=False, true_partitions=None, debug=False):  # NOTE the clusters in the resulting partition generally have the uids in a totally different order to in either of the original partitions
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
        # NOTE single_cluster and cluster_list in general have quite different sets of uids, and that's fine. All that matters here is we're trying to find all the clusters that should be split from one another (without doing some all against all horror)
        if len(cluster_list) < 2:  # nothing to do
            return [single_cluster]  # NOTE <single_cluster> doesn't get used after here
        adict = utils.get_annotation_dict(annotation_list)
        cdr3_groups = utils.group_seqs_by_value(cluster_list, lambda c: adict[akey(c)]['cdr3_length'])  # group the together clusters in <cluster_list> that have the same cdr3 (there's already utils.split_clusters_by_cdr3(), but it uses different inputs (e.g. sw_info) so i think it makes sense to not use it here)
        if tdbg:
            print '   %s one cluster vs %d clusters' % (utils.color('blue', 'syncing'), len(cluster_list))
            print '     split into %d cdr3 groups' % len(cdr3_groups)
        lo_hbound, hi_hbound = utils.get_naive_hamming_bounds(naive_hamming_bound_type, overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))  # these are the wider bounds, so < lo is almost certainly clonal, > hi is almost certainly not
        return_clusts = []
        for icdr, cdrgroup in enumerate(cdr3_groups):  # within each cdr3 group, split (i.e. use the cluster boundaries from cluster_list rather than single_cluster) if naive hfrac is > hi_hbound (but then there's shenanigans to adjudicate between different possibilities)
            if tdbg: print '      %s hfrac bound %.2f' % (utils.color('purple', 'icdr %d' % icdr), hi_hbound)

            # first figure out who needs to be split from whom
            clusters_to_split = {akey(c) : [] for c in cdrgroup}  # map from each cluster ('s key) to a list of clusters from which it should be split
            for c1, c2 in itertools.combinations(cdrgroup, 2):  # we could take account of the hfrac of both chains at this point, but looking at only the "split" one rather than the "merged" one, as we do here, is i think equivalent to assuming the merged one has zero hfrac, which is probably fine, since we only split if the split chain is very strongly suggesting we split
                hfrac = utils.hamming_fraction(adict[akey(c1)]['naive_seq'], adict[akey(c2)]['naive_seq'], align_if_necessary=True)  # all clusters with the same cdr3 len have been padded in waterer so their naive seqs are the same length
                if hfrac > hi_hbound:
                    clusters_to_split[akey(c1)].append(c2)
                    clusters_to_split[akey(c2)].append(c1)

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
                    if tdbg: print '     merging with size %d' % len(rclust)
                    rclust += cclust
                    found_one = True
                    break  # i.e. we just merge with the first one we find and stop looking; if there's more than one, it means we could merge all three together if we wanted (triangle inequality-ish, see diagram linked at top of fcn), but i doubt it'll matter either way, and this is easier
                if not found_one:
                    if tdbg: print '      y'
                    tmpclusts_for_return.append(cclust)  # if we didn't find an existing cluster that we can add it to, add it as a new cluster

            return_clusts += tmpclusts_for_return

        # if debug:
        #     print '      returning: %s' % ' '.join([str(len(c)) for c in return_clusts])
        #     # ptnprint(return_clusts)
        return return_clusts

    # ----------------------------------------------------------------------------------------
    init_partitions = {}
    for tch in utils.chains:
        if iparts is None or ploci[tch] not in iparts:
            init_partitions[tch] = cpaths[ploci[tch]].best()
        else:
            init_partitions[tch] = cpaths[ploci[tch]].partitions[iparts[ploci[tch]]]
            print '  %s using non-best partition index %d for %s (best is %d)' % (utils.color('red', 'note'), iparts[ploci[tch]], tch, cpaths[ploci[tch]].i_best)

    l_translations = translate_paired_uids(ploci, init_partitions, antn_lists)
    if debug:
        for tstr, tpart in [('heavy', init_partitions['h']), ('light', init_partitions['l'])]:
            ptnprint(tpart, extrastr=utils.color('blue', '%s  '%tstr), print_partition_indices=True, n_to_print=1, sort_by_size=False, print_header=tstr=='heavy')

    common_uids, _, _ = utils.check_intersection_and_complement(init_partitions['h'], init_partitions['l'], only_warn=True, a_label='heavy', b_label='light')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
    if len(common_uids) == 0:
        raise Exception('no uids in common between heavy and light')

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

        resolvedbg = False
        if debug:
            dbgheader = ['    adding %d resolved cluster%s to %d clusters in final partition' % (len(resolved_clusters), utils.plural(len(resolved_clusters)), len(final_partition)), '      ifclust N rclusts',]
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
                if len(new_fset) > len(rset):  # remove the common ids from the larger one (effectively splitting according to the splittier one)
                    new_fset -= common_uids
                    if debug: dbgstr.append('  fclust %d --> %d' % (len(new_fset) + len(common_uids), len(new_fset)))
                else:
                    rset -= common_uids
                    if debug: dbgstr.append('  rclust %d --> %d' % (len(rset) + len(common_uids), len(rset)))
                resolved_clusters[irclust] = list(rset)
            if resolvedbg:
                if len(dbgheader) > 0:
                    print '\n'.join(dbgheader)
                    dbgheader = []
                print '       %4d  %4d  %s' % (ifclust, len(irclusts), ''.join(dbgstr))
            final_partition[ifclust] = list(new_fset)
        # if debug:
        #     print '       %d fclusts clean' % n_clean
        assert is_clean_partition(resolved_clusters)
        final_partition += resolved_clusters

    if debug:
        print '    removing %d/%d empty clusters' % (final_partition.count([]), len(final_partition))
    final_partition = [c for c in final_partition if len(c) > 0]
    if debug:
        print '    final: %s' % ' '.join([str(len(c)) for c in final_partition])
    def chstr(n_before, n_after):
        if n_before == n_after: return ''
        else: return ' ' + utils.color('red', '%+d' % (n_after - n_before))
    print '   N clusters:\n        h %4d --> %-4d%s\n        l %4d --> %-4d%s'  % (len(init_partitions['h']), len(final_partition), chstr(len(init_partitions['h']), len(final_partition)),
                                                                                   len(init_partitions['l']), len(final_partition), chstr(len(init_partitions['l']), len(final_partition)))

    # ptnprint(final_partition, sort_by_size=False) #extrastr=utils.color('blue', '%s  '%tstr), print_partition_indices=True, n_to_print=1, sort_by_size=False, print_header=tstr=='heavy')

    if check_partitions:  # TODO work on this
        assert is_clean_partition(final_partition)
        for tch, initpart in init_partitions.items():
            _, _, _ = utils.check_intersection_and_complement(initpart, final_partition, only_warn=True, a_label=tch, b_label='joint')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
            assert len(set([u for c in initpart for u in c]) - set([u for c in final_partition for u in c])) == 0  # everybody from both initial partitions is in final_partition
        assert len(set([u for c in final_partition for u in c]) - set([u for c in init_partitions['h'] for u in c]) - set([u for c in init_partitions['l'] for u in c])) == 0  # nobody extra got added (i don't see how this could happen, but maybe it's just checking that I didnt' modify the initial partitions)

    joint_partitions = {ch : copy.deepcopy(final_partition) for ch in utils.chains}
    if len(l_translations) > 0:
        untranslate_pids(ploci, init_partitions, antn_lists, l_translations, joint_partitions)
    if true_partitions is not None:
        evaluate_joint_partitions(ploci, true_partitions, init_partitions, joint_partitions, antn_lists)

    if debug:
        for tch in utils.chains:
            ltmp = ploci[tch]
            print '%s' % utils.color('green', ltmp)
            for tclust in cpaths[ltmp].best():  # loop over clusters in the initial partition
                if len(tclust) == 1:
                    continue
                jfamilies = [c for c in joint_partitions[tch] if len(set(tclust) & set(c)) > 0]  # clusters in the joint partition that overlap with this cluster
                uid_extra_strs = None
                if len(jfamilies) == 0:
                    uid_extra_strs = [utils.color('blue', '?') for _ in tclust]
                elif len(jfamilies) == 1:  # everyone in <tclust> is also in the same joint partition family
                    uid_extra_strs = [utils.color('blue', '-') for _ in tclust]
                else:
                    def getjcstr(u):  # str for length of <u>'s cluster in <jfamilies>
                        jfs = [f for f in jfamilies if u in f]
                        return utils.color('blue', '?') if len(jfs)==0 else ('%d: %s'%(len(utils.get_single_entry(jfs)), ' '.join(jfs[0])))
                    uid_extra_strs = [getjcstr(u) for u in tclust]
                utils.print_reco_event(antn_dict[tch][':'.join(tclust)], uid_extra_strs=uid_extra_strs, extra_str='      ')

    return {ploci[ch] : jp for ch, jp in joint_partitions.items()}
