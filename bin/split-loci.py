#!/usr/bin/env python
import csv
import os
import sys
import argparse
import operator
import colored_traceback.always
import collections

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

dstr = """
use vsearch to split the sequences in <fname> according to their ig loci, writing each locus to its own fasta file <fname>-<locus>.fa.
If --split-heavy-seqs is set, also splits the igh sequences according to the light chain locus with which they\'re paired, thus resulting in two directories igh+igk/ and igh+igl/.
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
parser.add_argument('--debug', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--fasta-info-index', type=int, help='zero-based index in fasta info/meta string of sequence name/uid (e.g. if name line is \'>stuff more-stuff NAME extra-stuff\' the index should be 2)')
parser.add_argument('--split-heavy-seqs', action='store_true', help='In addition to writing a fasta file for each locus with every sequence from that locus, also split the igh sequences into separate subdirs for those sequences paired with igk and those paired with igl (so each of these two subdirs is ready for partis --paired-loci)')
parser.add_argument('--ig-or-tr', default='ig', choices=utils.locus_pairs.keys(), help='antibodies or TCRs?')
args = parser.parse_args()

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
if os.path.dirname(args.fname) == '':
    args.fname = '%s/%s' % (os.getcwd(), args.fname)
seqfos = utils.read_fastx(args.fname)
if args.fasta_info_index is not None:
    for sfo in seqfos:
        sfo['name'] = sfo['infostrs'][args.fasta_info_index]
if args.reverse_negative_strands:
    revfos = [{'name' : s['name'], 'seq' : utils.revcomp(s['seq'])} for s in seqfos]  # NOTE this is not on an equal footing with <seqfos>, since we add all the vsearch info to <seqfos>, then use it do decide on locus, and then to write output

if os.path.exists(args.germline_dir + '/' + args.species):  # ick that is hackey
    args.germline_dir += '/' + args.species

tmploci = [l for l in utils.loci if args.ig_or_tr in l]

# run vsearch to see if you can get a match for each locus for every sequence
n_rev_compd, n_total = 0, 0
for locus in tmploci:
    lglfo = glutils.read_glfo(args.germline_dir, locus)
    annotations = utils.run_vsearch_with_duplicate_uids('search', seqfos, args.workdir + '/vsearch', 0.3, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True, expect_failure=True, extra_str='   %s  fwd:'%utils.color('blue', locus) if args.reverse_negative_strands else '   %s: '%locus)
    assert len(annotations) == len(seqfos)
    if args.reverse_negative_strands:  # it might be nicer to user vsearch options to run on both senses at once, but otoh this might be nicer
        revnotations = utils.run_vsearch_with_duplicate_uids('search', revfos, args.workdir + '/vsearch', 0.3, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True, expect_failure=True, extra_str='        rev:')
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

# then, for each sequence, choose the locus with the best-scoring match (in practice i doubt you ever really get multiple loci with matches)
outfos = collections.OrderedDict(((l, []) for l in tmploci))
failed_ids = set()
for sfo in seqfos:
    lscores = {l : sfo[l]['score'] if 'invalid' not in sfo[l] else 0 for l in tmploci}
    locus, max_score = sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)[0]
    if max_score == 0:
        failed_ids.add(sfo['name'])
    outfos[locus].append(sfo)
    if args.debug:
        def lpstr(spair):
            l, s = spair
            return '%s %d' % (utils.color('blue' if l==locus else None, l), s)
        print '   %s: %s' % (sfo['name'], '  '.join(lpstr(s) for s in sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)))

print 'totals: %s%s' % (' '.join(('%s %d'%(l, len(sfos))) for l, sfos in outfos.items()), '' if len(failed_ids) == 0 else ' (%s: %d)'%(utils.color('yellow', 'failed'), len(failed_ids)))

# ----------------------------------------------------------------------------------------
def write_locus_file(locus, ofn, ofos, extra_str='  '):
    if len(ofos) == 0:
        print '%s%s: nothing to write' % (extra_str, locus)
        return
    if utils.output_exists(args, ofn, offset=4):
        return
    print '%s%s: %d to %s' % (extra_str, locus, len(ofos), ofn)
    if not os.path.exists(os.path.dirname(ofn)):
        os.makedirs(os.path.dirname(ofn))
    with open(ofn, 'w') as lfile:
        for sfo in ofos:
            lfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))

# ----------------------------------------------------------------------------------------
print 'writing:'
if args.outdir is None:
    outbase = utils.getprefix(args.fname)  # if --outdir wasn't set, we include basename of input file, since they might have multiple files in the same dir
    tmp_sep = '-'
else:
    outbase = args.outdir  # but if they did set --outdir, it's nicer to have plain locus file/path names
    tmp_sep = '/'
for locus in outfos:  # first write the single files with all seqs for each locus
    write_locus_file(locus, '%s%s%s.fa' % (outbase, tmp_sep, locus), outfos[locus])
if args.split_heavy_seqs:  # then, if necessary, write the ones that're split by pairing
    print '  writing to paired subdirs'
    for h_locus, l_locus in utils.locus_pairs[args.ig_or_tr]:
        l_uids = set(sfo['name'] for sfo in outfos[l_locus])
        h_outfo = [sfo for sfo in outfos[h_locus] if sfo['name'] in l_uids]  # heavy chain seqs corresponding to this light chain
        ldir = '%s%s%s+%s' % (outbase, tmp_sep, h_locus, l_locus)
        print '      %d/%d %s seqs pair with %s' % (len(h_outfo), len(outfos[h_locus]), h_locus, l_locus)
        if len(h_outfo) == 0:
            print '      %s no heavy chains paired with %s' % (utils.color('yellow', 'warning'), l_locus)
        write_locus_file(h_locus, '%s/%s.fa' % (ldir, h_locus), h_outfo, extra_str='    ')
        write_locus_file(l_locus, '%s/%s.fa' % (ldir, l_locus), outfos[l_locus], extra_str='    ')
