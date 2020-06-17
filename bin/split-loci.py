#!/usr/bin/env python
import csv
import os
import sys
import argparse
import operator
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser(description='split the sequences in the input fasta <fname> according to their ig loci, writing to fasta output files <fname>-<locus>.fa')
parser.add_argument('fname', help='fasta input file')
parser.add_argument('--species', default='human', choices=('human', 'macaque', 'mouse'), help='Which species?')
parser.add_argument('--germline-dir', default=partis_dir + '/data/germlines')
parser.add_argument('--workdir', default=utils.choose_random_subdir('/tmp/%s/partis' % os.getenv('USER', default='partis-work')))
parser.add_argument('--vsearch-binary', help='Path to vsearch binary (vsearch binaries for linux and darwin are pre-installed in bin/, but for other systems you need to get your own)')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--fasta-info-index', type=int, help='index in fasta info of sequence name string')
args = parser.parse_args()

seqfos = utils.read_fastx(args.fname)
if args.fasta_info_index is not None:
    for sfo in seqfos:
        sfo['name'] = sfo['infostrs'][args.fasta_info_index]

if os.path.exists(args.germline_dir + '/' + args.species):  # ick that is hackey
    args.germline_dir += '/' + args.species

seq_dict = {sfo['name'] : sfo['seq'] for sfo in seqfos}
vs_infos = {}
for locus in [l for l in utils.loci if 'ig' in l]:
    lglfo = glutils.read_glfo(args.germline_dir, locus)
    vs_infos[locus] = utils.run_vsearch('search', seq_dict, args.workdir + '/vsearch', threshold=0.3, glfo=lglfo, print_time=True, vsearch_binary=args.vsearch_binary, get_annotations=True)

outfos = {l : [] for l in vs_infos}
for sfo in seqfos:
    lscores = {l : vs_infos[l]['annotations'][sfo['name']]['score'] if sfo['name'] in vs_infos[l]['annotations'] else 0 for l in vs_infos}
    locus, max_score = sorted(lscores.items(), key=operator.itemgetter(1), reverse=True)[0]
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
