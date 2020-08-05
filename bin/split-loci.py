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

parser = argparse.ArgumentParser(description='split the sequences in the input fasta <fname> according to their ig loci, writing to fasta output files <fname>-<locus>.fa', formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('fname', help='fasta input file')
parser.add_argument('--species', default='human', choices=('human', 'macaque', 'mouse'), help='Which species?')
parser.add_argument('--germline-dir', default=partis_dir + '/data/germlines', help='doesn\'t need to be the germlines corresponding to this sample since it\'s just so it can figure out which is igh vs igk vs igl, so the default is probably fine')
parser.add_argument('--workdir', default=utils.choose_random_subdir('/tmp/%s/partis' % os.getenv('USER', default='partis-work')))
parser.add_argument('--vsearch-binary', help='Path to vsearch binary (vsearch binaries for linux and darwin are pre-installed in bin/, so leaving this unset should work, but for other systems you need to get your own)')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--fasta-info-index', type=int, help='index in fasta info of sequence name string')
args = parser.parse_args()

seqfos = utils.read_fastx(args.fname)
if args.fasta_info_index is not None:
    for sfo in seqfos:
        sfo['name'] = sfo['infostrs'][args.fasta_info_index]

if os.path.exists(args.germline_dir + '/' + args.species):  # ick that is hackey
    args.germline_dir += '/' + args.species

# run vsearch to see if you can get a match for each locus for every sequence
tmploci = [l for l in utils.loci if 'ig' in l]
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

print 'writing:'
for locus in outfos:
    if len(outfos[locus]) == 0:
        continue
    lfname = '%s-%s.fa' % (utils.getprefix(args.fname), locus)
    print '    %s: %d to %s' % (locus, len(outfos[locus]), lfname)
    with open(lfname, 'w') as lfile:
        for sfo in outfos[locus]:
            lfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
