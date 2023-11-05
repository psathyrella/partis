#!/usr/bin/env python
from __future__ import absolute_import
import json
import csv
import os
import sys
import argparse
import operator
import colored_traceback.always
import collections
import random
import numpy

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import paircluster
import glutils
from clusterpath import ClusterPath
import seqfileopener

# ----------------------------------------------------------------------------------------
dstr = """
Uses vsearch (or the \'locus\' key in --input-metfname) to split the sequences in <fname> according to their loci, writing each locus to its own fasta file <locus>.fa.
If \'paired-uids\' are available in --input-metafname, also splits the heavy sequences according to the light chain locus with which they\'re paired, resulting in subdirectories e.g. igh+igk/ and igh+igl/.
Use --reverse-negative-strands to check both senses for each input sequence.
"""
parser = argparse.ArgumentParser(description=dstr,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('fname', help='fasta input file')
parser.add_argument('--outdir', help='directory to which to write output files (if not set, output is written to directory of <fname>)')
parser.add_argument('--reverse-negative-strands', action='store_true', help='align every sequence both forwards and revcomp\'d, then for each sequence keep the sense with better alignment.')
parser.add_argument('--species', default='human', choices=('human', 'macaque', 'mouse'), help='Which species?')
parser.add_argument('--germline-dir', default=partis_dir + '/data/germlines', help='doesn\'t need to be the germlines corresponding to this sample since it\'s just so it can figure out which is igh vs igk vs igl, so the default is probably fine')
parser.add_argument('--workdir', default=utils.choose_random_subdir('/tmp/%s/partis' % os.getenv('USER', default='partis-work')), help='working directory for vsearch')
parser.add_argument('--vsearch-binary', help='Path to vsearch binary (vsearch binaries for linux and darwin are included in partis/bin/, so leaving this unset should work, but for other systems you need to get your own)')
parser.add_argument('--vsearch-threshold', type=float, default=0.4, help='default identity threshold for vsearch')
parser.add_argument('--debug', type=int, default=1)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--random-seed', type=int, default=1)
parser.add_argument('--guess-pairing-info', action='store_true', help=utils.did_help['guess'])
parser.add_argument('--droplet-id-separators', help=utils.did_help['seps'])
parser.add_argument('--droplet-id-indices', help=utils.did_help['indices'])
parser.add_argument('--fasta-info-index', type=int, help='zero-based index in fasta info/meta string of sequence name/uid (e.g. if name line is \'>stuff more-stuff NAME extra-stuff\' the index should be 2)')
parser.add_argument('--input-metafname', help='yaml file with meta information keyed by sequence id. See same argument in main partis help, and https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#input-meta-info for an example.')
parser.add_argument('--for-testing-n-max-queries', type=int, default=-1, help='only for testing, applied when reading initial fasta file, just in case it\'s huge and you want to run quickly without having to read the whole file')
parser.add_argument('--n-max-queries', type=int, default=-1, help='see partis help (although here it applies to droplets, not individual seqs)')
parser.add_argument('--n-random-queries', type=int, help='see partis help (although here it applies to droplets, not individual seqs)')
parser.add_argument('--ig-or-tr', default='ig', choices=utils.locus_pairs.keys(), help='antibodies or TCRs?')

# ----------------------------------------------------------------------------------------
def use_rev_comp(pline, rline):  # decide whether positive sense <pline> or negative sense <rline> has better alignment
    assert pline['unique_ids'][0] == rline['unique_ids'][0]
    if rline.get('invalid', False):
        return False
    elif pline.get('invalid', False):
        return True
    elif rline['score'] > pline['score']:
        return True
    else:
        return False

# ----------------------------------------------------------------------------------------
def run_vsearch(seqfos):  # run vsearch to see if you can get a match for each locus for every sequence
    print '  running vsearch on %d sequences:' % len(seqfos)
    n_rev_compd, n_total = 0, 0
    for locus in utils.sub_loci(args.ig_or_tr):
        lglfo = glutils.read_glfo(args.germline_dir, locus)
        annotations = utils.run_vsearch_with_duplicate_uids('search', seqfos, args.workdir + '/vsearch', args.vsearch_threshold, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True, expect_failure=True, extra_str='   %s  fwd:'%utils.color('blue', locus) if args.reverse_negative_strands else '    %s: '%locus)
        assert len(annotations) == len(seqfos)
        if args.reverse_negative_strands:  # it might be nicer to user vsearch options to run on both senses at once, but otoh this might be nicer
            revnotations = utils.run_vsearch_with_duplicate_uids('search', revfos, args.workdir + '/vsearch', args.vsearch_threshold, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True, expect_failure=True, extra_str='        rev:')
            assert len(revnotations) == len(seqfos)
        for il, (sfo, line) in enumerate(zip(seqfos, annotations)):
            assert sfo['name'] == line['unique_ids'][0]  # note that they're not full annotations, they just have a couple keys
            if args.reverse_negative_strands and use_rev_comp(line, revnotations[il]):
                sfo['seq'] = revfos[il]['seq']
                line = revnotations[il]
                n_rev_compd += 1
            sfo[locus] = line  # add info for each locus to the input seqfos
            n_total += 1
    if args.reverse_negative_strands:
        print '  used rev comp for %d/%d locus results (for %d seqs)' % (n_rev_compd, n_total, len(seqfos))

# ----------------------------------------------------------------------------------------
def write_locus_file(locus, ofos, lpair=None, extra_str='  '):
    ofn = paircluster.paired_fn(args.outdir, locus=locus, lpair=lpair)
    if utils.output_exists(args, ofn, leave_zero_len=len(ofos)==0, offset=4):  # NOTE not really sure this does anything (or if i want it) now that I'm cleaning/looking for the whole dir at the start of this script
        return
    if not os.path.exists(os.path.dirname(ofn)):
        os.makedirs(os.path.dirname(ofn))
    if len(ofos) == 0:
        # print '%s%s: nothing to write' % (extra_str, locus)
        open(ofn, 'w').close()
        return
    print '%s%s: %d to %s/%s' % (extra_str, locus, len(ofos), os.path.basename(os.path.dirname(ofn)), os.path.basename(ofn))
    with open(ofn, 'w') as lfile:
        for sfo in ofos:
            lfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))

# ----------------------------------------------------------------------------------------
def read_meta_info(seqfos):  # read all input meta info, and add pairing info (if present) to <paired_uids>
    dummy_annotation_list = [{'unique_ids' : [sfo['name']]} for sfo in seqfos]
    seqfileopener.read_input_metafo([args.input_metafname], dummy_annotation_list)  # , required_keys=['paired-uids'])
    for line in dummy_annotation_list:
        uid = utils.get_single_entry(line['unique_ids'])
        if 'loci' in line:
            meta_loci[uid] = line['loci'][0]
        if 'paired-uids' in line:
            paired_uids[uid] = line['paired-uids'][0]
    if len(paired_uids) > 0:
        print '    read pairing info for %d seqs from input meta file' % len(paired_uids)
        if len(paired_uids) < len(seqfos):
            print '      %s only read pairing info for %d/%d seqfos' % (utils.color('yellow', 'warning'), len(paired_uids), len(seqfos))
    if len(meta_loci) > 0:
        print '    read loci for %d sequences from input meta file (so not running vsearch)' % len(meta_loci)
        if len(meta_loci) < len(seqfos):
            print '      %s only read locus info for %d/%d seqfos' % (utils.color('yellow', 'warning'), len(meta_loci), len(seqfos))
    input_metafos = utils.read_json_yaml(args.input_metafname)
    for uid in input_metafos:  # we want to copy over any additional meta info (not paired uids or loci) to the output meta info file (since if we're guessing pair info, the uid names will change, so the original one is no good)
        additional_mfo = {k : v for k, v in input_metafos[uid].items() if k not in ['loci', 'paired-uids']}
        if len(additional_mfo) > 0:
            input_metafos[uid] = additional_mfo
    return input_metafos

# ----------------------------------------------------------------------------------------
def print_pairing_info(outfos, paired_uids):
    loci_by_uid = {sfo['name'] : l for l in outfos for sfo in outfos[l]}  # locus of each sequence, just for counting below
    print_cutoff = 0.01
    print '            count  frac  paired with'
    for locus in utils.sub_loci(args.ig_or_tr):
        plocicounts = {}
        for sfo in outfos[locus]:
            plstr = ' '.join(utils.locstr(l) for l in sorted([loci_by_uid.get(pid, '?') for pid in paired_uids[sfo['name']]]))
            if plstr not in plocicounts:
                plocicounts[plstr] = 0
            plocicounts[plstr] += 1
        total = sum(plocicounts.values())
        n_skipped = 0
        for ipl, (plstr, counts) in enumerate(sorted(plocicounts.items(), key=operator.itemgetter(1), reverse=True)):
            if counts / float(total) < print_cutoff:
                n_skipped += counts
                continue
            print '       %s  %6d  %5.2f   %s' % (utils.locstr(locus) if ipl==0 else ' ', counts, counts / float(total), plstr)
        if n_skipped > 0:
            print '                +%d counts skipped with <%.3f each' % (n_skipped , print_cutoff) # utils.color('yellow', 'note

# ----------------------------------------------------------------------------------------
args = parser.parse_args()
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if os.path.dirname(args.fname) == '':
    args.fname = '%s/%s' % (os.getcwd(), args.fname)
if args.outdir is None:
    args.outdir = utils.getprefix(args.fname)
args.droplet_id_indices = utils.get_arg_list(args.droplet_id_indices, intify=True)
args.input_metafname = utils.fpath(args.input_metafname)

if any(os.path.exists(ofn) for ofn in paircluster.paired_dir_fnames(args.outdir)):
    if args.overwrite:
        paircluster.clean_paired_dir(args.outdir)
    else:
        print '  split-loci.py output exists and --overwrite was not set, so not doing anything: %s' % args.outdir
        sys.exit(0)

seqfos = utils.read_fastx(args.fname, n_max_queries=args.for_testing_n_max_queries)
if args.n_max_queries != -1 or args.n_random_queries is not None:
    seqfos = utils.subset_paired_queries(seqfos, args.droplet_id_separators, args.droplet_id_indices, n_max_queries=args.n_max_queries, n_random_queries=args.n_random_queries)
if args.fasta_info_index is not None:
    for sfo in seqfos:
        sfo['name'] = sfo['infostrs'][args.fasta_info_index]
if args.reverse_negative_strands:
    revfos = [{'name' : s['name'], 'seq' : utils.revcomp(s['seq'])} for s in seqfos]  # NOTE this is not on an equal footing with <seqfos>, since we add all the vsearch info to <seqfos>, then use it do decide on locus, and then to write output

if os.path.exists(args.germline_dir + '/' + args.species):  # ick that is hackey
    args.germline_dir += '/' + args.species

# read input meta file and/or run vsearch
paired_uids, meta_loci, input_metafos = {}, {}, {}
if args.input_metafname is not None:
    input_metafos = read_meta_info(seqfos)
if len(meta_loci) == 0:  # default: no input locus info
    run_vsearch(seqfos)

# then, for each sequence, choose the locus with the best-scoring match (in practice i doubt you ever really get multiple loci with matches)
outfos = collections.OrderedDict(((l, []) for l in utils.sub_loci(args.ig_or_tr)))
failed_seqs = []
if args.debug > 1:
    print '    printing scores for locus determination:'
    n_skipped = 0
for sfo in seqfos:
    if len(meta_loci) == 0:  # default: use vsearch match scores
        lscores = {l : sfo[l]['score'] if 'invalid' not in sfo[l] else 0 for l in utils.sub_loci(args.ig_or_tr)}
        locus, max_score = sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)[0]
        if max_score == 0:
            failed_seqs.append(sfo)
            continue
    else:  # if we were passed input locus info
        locus = meta_loci[sfo['name']]
    outfos[locus].append(sfo)
    if args.debug > 1:
        def lpstr(spair):
            l, s = spair
            return '%s %s' % (utils.locstr(l) if l==locus else l.replace('ig', ''), utils.color('red' if s!=0 else None, '%3d'%s))
        if lscores.values().count(0) == 2:
            n_skipped += 1
        else:
            print '       %s   %s' % ('  '.join(lpstr(s) for s in sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)), sfo['name'])
if args.debug > 1 and n_skipped > 0:
    print '      skipped %d seqs with non-zero scores from only one locus' % n_skipped

print 'totals: %s%s' % (' '.join(('%s %d'%(l, len(sfos))) for l, sfos in outfos.items()), '' if len(failed_seqs) == 0 else ' (%s: %d)'%(utils.color('yellow', 'failed'), len(failed_seqs)))
assert sum(len(ofo) for ofo in outfos.values()) + len(failed_seqs) == len(seqfos)

if args.guess_pairing_info:
    if len(paired_uids) > 0:
        raise Exception('can\'t/shouldn\'t guess pairing info if we already have it from elsewhere')
    for locus in outfos:
        for ofo in outfos[locus]:
            new_name = ofo['name'] + '-' + locus
            if ofo['name'] in input_metafos:
                input_metafos[new_name] = input_metafos[ofo['name']]
                del input_metafos[ofo['name']]
            ofo['name'] = new_name
    guessed_metafos = utils.extract_pairing_info(seqfos, droplet_id_separators=args.droplet_id_separators, droplet_id_indices=args.droplet_id_indices, debug=max(1, args.debug))
    for uid in set(guessed_metafos) & set(input_metafos):
        guessed_metafos[uid].update(input_metafos[uid])
    for uid, mfo in guessed_metafos.items():
        paired_uids[uid] = mfo['paired-uids']

# remove failed uids from paired_uids
failed_uids = set(s['name'] for s in failed_seqs)
n_removed = 0
for fid in failed_uids:
    if fid in paired_uids:
        del paired_uids[fid]
        n_removed += 1
paired_uids = {uid : list(set(paired_uids[uid]) - failed_uids) for uid in paired_uids}
print '  removed %d failed uids from paired_uids' % n_removed

if args.debug and len(paired_uids) > 0:
    print_pairing_info(outfos, paired_uids)

print 'writing to %s/' % args.outdir
if len(failed_seqs) > 0:
    write_locus_file('failed', failed_seqs)

for locus in outfos:  # first write the single files with all seqs for each locus
    write_locus_file(locus, outfos[locus])

omfname = '%s/meta.yaml' % args.outdir
if args.guess_pairing_info:
    with open(omfname, 'w') as outfile:  # NOTE file name duplicates code in bin/partis
        json.dump(guessed_metafos, outfile)
elif args.input_metafname is not None and not os.path.exists(omfname):
    utils.makelink(os.path.dirname(omfname), args.input_metafname, omfname)

if len(paired_uids) == 0:
    print 'no pairing info'
else:
    print 'writing to paired subdirs'
    for lpair in utils.locus_pairs[args.ig_or_tr]:
        h_locus, l_locus = lpair
        all_paired_uids = set(pid for s in outfos[l_locus] for pid in paired_uids[s['name']])  # all uids that are paired with any <l_locus> uid (note that in general not all of these are heavy chain)
        h_outfo = [sfo for sfo in outfos[h_locus] if sfo['name'] in all_paired_uids]  # any h seq that was in a droplet with any of <l_uids>
        # # old way of pairing with identical uids (should really remove this but too chicken atm):
        #     l_uids = set(sfo['name'] for sfo in outfos[l_locus])
        #     h_outfo = [sfo for sfo in outfos[h_locus] if sfo['name'] in l_uids]  # if no meta file specified, just pair up the ones with identical uids

        print '    %d/%d %s seqs pair with %s%s' % (len(h_outfo), len(outfos[h_locus]), h_locus, l_locus, utils.color('yellow', ' (warning)') if len(h_outfo)==0 else '')
        if len(h_outfo) == 0:
            outfos[l_locus] = []  # makes more sense to not write them if there's no heavy paired with them, since it'd just duplicate the unpaired ones in the parent dir
        write_locus_file(h_locus, h_outfo, lpair=lpair, extra_str='      ')
        write_locus_file(l_locus, outfos[l_locus], lpair=lpair, extra_str='      ')
