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
use vsearch to split the sequences in <fname> according to their ig loci, either writing each locus to its own fasta file <fname>-<locus>.fa (by default), 
or additionally splitting the igh sequences according to the light chain locus with which they\'re paired, thus resulting in two directories igh+igk/ and igh+igl/ (if --split-heavy-seqs is set).
"""
parser = argparse.ArgumentParser(description=dstr,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('fname', help='fasta input file')
parser.add_argument('--outdir', help='output directory (if not set, output is written to directory of <fname>)')
parser.add_argument('--species', default='human', choices=('human', 'macaque', 'mouse'), help='Which species?')
parser.add_argument('--germline-dir', default=partis_dir + '/data/germlines', help='doesn\'t need to be the germlines corresponding to this sample since it\'s just so it can figure out which is igh vs igk vs igl, so the default is probably fine')
parser.add_argument('--workdir', default=utils.choose_random_subdir('/tmp/%s/partis' % os.getenv('USER', default='partis-work')), help='working directory for vsearch')
parser.add_argument('--vsearch-binary', help='Path to vsearch binary (vsearch binaries for linux and darwin are included in partis/bin/, so leaving this unset should work, but for other systems you need to get your own)')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--fasta-info-index', type=int, help='zero-based index in fasta info of sequence name string (e.g. if name line is \'>stuff more-stuff NAME extra-stuff\' the index should be 2)')
parser.add_argument('--split-heavy-seqs', action='store_true', help='instead of the default of writing each locus to its own .fa file, split igh into separate dirs for those sequences paired with igk and igl (i.e. so each of the two dirs is ready for partis --paired-loci)')
parser.add_argument('--ig-or-tr', default='ig', choices=utils.locus_pairs.keys(), help='antibodies or TCRs?')
args = parser.parse_args()

seqfos = utils.read_fastx(args.fname)
if args.fasta_info_index is not None:
    for sfo in seqfos:
        sfo['name'] = sfo['infostrs'][args.fasta_info_index]

if os.path.exists(args.germline_dir + '/' + args.species):  # ick that is hackey
    args.germline_dir += '/' + args.species

tmploci = [l for l in utils.loci if args.ig_or_tr in l]

# run vsearch to see if you can get a match for each locus for every sequence
for locus in tmploci:
    lglfo = glutils.read_glfo(args.germline_dir, locus)
    annotations = utils.run_vsearch_with_duplicate_uids('search', seqfos, args.workdir + '/vsearch', 0.3, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True, expect_failure=True, extra_str='   %s: '%locus)
    assert len(annotations) == len(seqfos)
    for sfo, line in zip(seqfos, annotations):
        assert sfo['name'] == line['unique_ids'][0]  # note that they're not full annotations, they just have a couple keys
        sfo[locus] = line

# then, for each sequence, choose the locus with the best-scoring match (in practice i doubt you ever really get multiple loci with matches)
outfos = collections.OrderedDict(((l, []) for l in tmploci))
for sfo in seqfos:
    lscores = {l : sfo[l]['score'] if 'invalid' not in sfo[l] else 0 for l in tmploci}
    locus, max_score = sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)[0]
    if max_score == 0:
        raise Exception('max score for locus %s is zero' % locus)
    outfos[locus].append(sfo)
    if args.debug:
        def lpstr(spair):
            l, s = spair
            return '%s %d' % (utils.color('blue' if l==locus else None, l), s)
        print '   %s: %s' % (sfo['name'], '  '.join(lpstr(s) for s in sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)))

print 'totals: %s' % ' '.join(('%s %d'%(l, len(sfos))) for l, sfos in outfos.items())

# ----------------------------------------------------------------------------------------
def write_locus_file(locus, ofn, ofos):
    if len(ofos) == 0:
        return
    if utils.output_exists(args, ofn, offset=4):
        return
    print '    %s: %d to %s' % (locus, len(ofos), ofn)
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
    for h_locus, l_locus in utils.locus_pairs[args.ig_or_tr]:
        l_uids = set(sfo['name'] for sfo in outfos[l_locus])
        h_outfo = [sfo for sfo in outfos[h_locus] if sfo['name'] in l_uids]  # heavy chain seqs corresponding to this light chain
        ldir = '%s%s%s+%s' % (outbase, tmp_sep, h_locus, l_locus)
        print '  %d/%d %s seqs pair with %s' % (len(h_outfo), len(outfos[h_locus]), h_locus, l_locus)
        write_locus_file(h_locus, '%s/%s.fa' % (ldir, h_locus), h_outfo)
        write_locus_file(l_locus, '%s/%s.fa' % (ldir, l_locus), outfos[l_locus])
