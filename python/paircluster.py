import copy
import itertools
import numpy
import sys

import utils
from clusterpath import ptnprint

# ----------------------------------------------------------------------------------------
# rename all uids in the light chain partition and annotations that are paired with a heavy chain uid to that heavy chain uid (pairings must, at this stage, be unique)
def translate_paired_uids(ploci, init_partitions, antn_lists):
    h_paired_uids = {}  # map ot each heavy chain uid <u> from its paired light chain uid <pids[0]>
    for hline in antn_lists[ploci['h']]:
        for h_id, pids in zip(hline['unique_ids'], hline['paired-uids']):
            if len(pids) > 1:
                raise Exception('multiple paired uids %s for %s sequece %s' % (' '.join(pids), ploci['h'], h_id))
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
def untranslate_pids(ploci, init_partitions, antn_lists, l_translations):
    for lline in antn_lists[ploci['l']]:
        lline['unique_ids'] = [l_translations.get(u, u) for u in lline['unique_ids']]
    init_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in init_partitions['l']]

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
        if len(cluster_list) == 1:  # nothing to do
            return [single_cluster]  # NOTE <single_cluster> doesn't get used after here
        adict = utils.get_annotation_dict(annotation_list)
        cdr3_groups = utils.group_seqs_by_value(cluster_list, lambda c: adict[akey(c)]['cdr3_length'])  # group the together clusters in <cluster_list> that have the same cdr3 (there's already utils.split_clusters_by_cdr3(), but it uses different inputs (e.g. sw_info) so i think it makes sense to not use it here)
        if tdbg:
            print '   %s one cluster vs %d clusters' % (utils.color('blue', 'syncing'), len(cluster_list))
            print '     split into %d cdr3 groups' % len(cdr3_groups)
        lo_hbound, hi_hbound = utils.get_naive_hamming_bounds('likelihood', overall_mute_freq=numpy.mean([f for l in annotation_list for f in l['mut_freqs']]))  # these are the wider bounds, so < lo is almost certainly clonal, > hi is almost certainly not
        return_clusts = []
        for icdr, cdrgroup in enumerate(cdr3_groups):  # within each cdr3 group, split (i.e. use the cluster boundaries from cluster_list rather than single_cluster) if naive hfrac is > hi_hbound (but then there's shenanigans to adjudicate between different possibilities)
            if tdbg: print '      %s hfrac bound %.2f' % (utils.color('purple', 'icdr %d' % icdr), hi_hbound)

            # first figure out who needs to be split from whom
            clusters_to_split = {akey(c) : [] for c in cdrgroup}  # map from each cluster ('s key) to a list of clusters from which it should be split
            for c1, c2 in itertools.combinations(cdrgroup, 2):  # we could take account of the hfrac of both chains at this point, but looking at only the "split" one rather than the "merged" one, as we do here, is i think equivalent to assuming the merged one has zero hfrac, which is probably fine, since we only split if the split chain is very strongly suggesting we split
                hfrac = utils.hamming_fraction(adict[akey(c1)]['naive_seq'], adict[akey(c2)]['naive_seq'])  # all clusters with the same cdr3 len have been padded in waterer so their naive seqs are the same length
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

        if debug:
            print '      returning: %s' % ' '.join([str(len(c)) for c in return_clusts])
            # ptnprint(return_clusts)
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
        print '    N        N       hclusts     lclusts       h/l'
        print '  hclusts  lclusts    sizes       sizes      overlaps'
    # For each single cluster in each partition, get a list of the clusters in the other partition that have common uids
    # Pass this cluster + list to a fcn to resolve discrepancies by splitting on the cluster boundaries in <cluster_list> that we're sure of (i.e. that have different cdr3, or very different naive hamming fraction)
    for h_initclust, l_initclust in [(c, None) for c in init_partitions['h']] + [(None, c) for c in init_partitions['l']]:  # just loops over each single cluster in h and l partitions, but in a way that we know whether the single cluster is from h or l
        single_chain, list_chain = 'h' if l_initclust is None else 'l', 'l' if l_initclust is None else 'h'
        single_cluster = h_initclust if single_chain == 'h' else l_initclust
        cluster_list = common_clusters(single_cluster, init_partitions[list_chain])
        single_annotation = antn_dict[single_chain][akey(single_cluster)]
        annotation_list = [antn_dict[list_chain][akey(c)] for c in cluster_list]

        if debug:
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

        resolved_clusters = resolve_discordant_clusters(copy.deepcopy(single_cluster), single_annotation, copy.deepcopy(cluster_list), annotation_list)
        if check_partitions:
            assert is_clean_partition(resolved_clusters)
        if debug:
            print '    adding %d resolved cluster%s to %d clusters in final partition' % (len(resolved_clusters), utils.plural(len(resolved_clusters)), len(final_partition))
            print '      ifclust N rclusts'
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
            if debug:
                print '       %4d  %4d  %s' % (ifclust, len(irclusts), ''.join(dbgstr))
            final_partition[ifclust] = list(new_fset)
        if debug:
            print '       %d fclusts clean' % n_clean
        assert is_clean_partition(resolved_clusters)
        final_partition += resolved_clusters

    if debug:
        print '    removing %d/%d empty clusters' % (final_partition.count([]), len(final_partition))
    final_partition = [c for c in final_partition if len(c) > 0]
    # if debug:
    #     print '    final: %s' % ' '.join([str(len(c)) for c in final_partition])
    def chstr(n_before, n_after):
        if n_before == n_after: return ''
        else: return ' ' + utils.color('red', '%+d' % (n_after - n_before))
    print '   N clusters:\n        h %4d --> %-4d%s\n        l %4d --> %-4d%s'  % (len(init_partitions['h']), len(final_partition), chstr(len(init_partitions['h']), len(final_partition)),
                                                                                   len(init_partitions['l']), len(final_partition), chstr(len(init_partitions['l']), len(final_partition)))

    if check_partitions:
        assert is_clean_partition(final_partition)
        for tch, initpart in init_partitions.items():
            _, _, _ = utils.check_intersection_and_complement(initpart, final_partition, only_warn=True, a_label=tch, b_label='joint')  # check that h and l partitions have the same uids (they're expected to be somewhat different because of either failed queries or duplicates [note that this is why i just turned off default duplicate removal])
            assert len(set([u for c in initpart for u in c]) - set([u for c in final_partition for u in c])) == 0  # everybody from both initial partitions is in final_partition
        assert len(set([u for c in final_partition for u in c]) - set([u for c in init_partitions['h'] for u in c]) - set([u for c in init_partitions['l'] for u in c])) == 0  # nobody extra got added (i don't see how this could happen, but maybe it's just checking that I didnt' modify the initial partitions)

    joint_partitions = {ch : copy.deepcopy(final_partition) for ch in utils.chains}
    if len(l_translations) > 0:
        untranslate_pids(ploci, init_partitions, antn_lists, l_translations)
        joint_partitions['l'] = [[l_translations.get(u, u) for u in c] for c in joint_partitions['l']]
    if true_partitions is not None:
        evaluate_joint_partitions(ploci, true_partitions, init_partitions, joint_partitions, antn_lists)

    return {ploci[ch] : jp for ch, jp in joint_partitions.items()}
