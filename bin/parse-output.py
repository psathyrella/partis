#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

helpstr = """
Extract sequences from a partis output file and write them to a fasta, csv, or tsv file, optionally with a limited amount of extra information for each sequence.
To write airr-standard tsv files, use the partis --airr-output option.
For details of partis output files, see the manual.
To view the partitions and annotations in a partis output file, use the partis \'view-output\' action.
Example usage:
    ./bin/parse-output.py test/reference-results/partition-new-simu.yaml out.fa
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('infile', help='partis output file from which to read input')
parser.add_argument('outfile', help='file to which to write output extracted from <infile> (fasta or csv/tsv)')
parser.add_argument('--extra-columns', help='colon-separated list of additional partis output columns (beyond sequences), to write to the output file. If writing to a fasta file, the column values are appended after the sequence name, separated by --fasta-info-separator. If writing to csv/tsv, they\'re written as proper, labeled columns.')
parser.add_argument('--partition-index', type=int, help='if set, use the partition at this index in the cluster path, rather than the default of using the best partition')
parser.add_argument('--seed-unique-id', help='if set, take sequences only from the cluster containing this seed sequence, rather than the default of taking all sequences from all clusters')
parser.add_argument('--cluster-index', type=int, help='if set, take sequences only from the cluster at this index in the partition, rather than the default of taking all sequences from all clusters. This index is with respect to the cluster order found in the file (which, in contrast to plots made by --plotdir, is *not* sorted by size)')
parser.add_argument('--indel-reversed-seqs', action='store_true', help='if set, take sequences that have had any shm indels "reversed" (i.e. insertions are reversed, and deletions are replaced with the germline bases) rather than the default of using sequences from the original input file. Indel-reversed sequences can be convenient because they are by definition the same length as and aligned to the naive sequence.')
parser.add_argument('--glfo-dir', help='Directory with germline info. Only necessary for old-style csv output files. Equivalent to a parameter dir with \'/hmm/germline-sets\' appended.')
parser.add_argument('--locus', default='igh', help='only used for old-style csv output files')
parser.add_argument('--plotdir', help='if set, plot annotation parameters from infile to --plotdir and exit (you still have to set outfile, sorry, it\'s nice having it be a positional arg, but it doesn\'t get used). To add e.g. per-gene-per-position plots comment/uncomment args in the call below.')
parser.add_argument('--fasta-info-separator', default=' ', help='character to use ')

if 'extract-fasta.py' in sys.argv[0]:  # if they're trying to run this old script, which is now just a link to this one, print a warning and rejigger the arguments so it still works
    print '  note: running deprecated script %s, which currently is just a link pointing to %s' % (os.path.basename(sys.argv[0]), os.path.basename(os.path.realpath( __file__)))
    print '  note: transferring deprecated arguments --input-file and --fasta-output-file to the first two positional arguments (this will continue to work, you only need to change things if you want this warning to go away)'
    utils.insert_in_arglist(sys.argv, [utils.get_val_from_arglist(sys.argv, '--input-file'), utils.get_val_from_arglist(sys.argv, '--fasta-output-file')], sys.argv[0])
    utils.remove_from_arglist(sys.argv, '--input-file', has_arg=True)
    utils.remove_from_arglist(sys.argv, '--fasta-output-file', has_arg=True)

args = parser.parse_args()
args.extra_columns = utils.get_arg_list(args.extra_columns)
assert utils.getsuffix(args.outfile) in ['.csv', '.tsv', '.fa', '.fasta']

default_glfo_dir = partis_dir + '/data/germlines/human'
if utils.getsuffix(args.infile) == '.csv' and args.glfo_dir is None:
    print '  note: reading deprecated csv format, so need to get germline info from a separate directory; --glfo-dir was not set, so using default %s. If it doesn\'t crash, it\'s probably ok.' % default_glfo_dir
    args.glfo_dir = default_glfo_dir
glfo, annotation_list, cpath = utils.read_output(args.infile, glfo_dir=args.glfo_dir, locus=args.locus)

if args.plotdir is not None:
    from parametercounter import ParameterCounter
    setattr(args, 'region_end_exclusions', {r : [0 for e in ['5p', '3p']] for r in utils.regions})  # hackity hackity hackity
    pcounter = ParameterCounter(glfo, args)
    for line in annotation_list:
        pcounter.increment(line)
    pcounter.plot(args.plotdir) #, make_per_base_plots=True) #, only_overall=True, make_per_base_plots=True
    sys.exit(0)

if cpath is None or cpath.i_best is None:
    clusters_to_use = [l['unique_ids'] for l in annotation_list]
    print '  no cluster path in input file, so just using all %d sequences (in %d clusters) in annotations' % (sum(len(c) for c in clusters_to_use), len(clusters_to_use))
else:
    ipartition = cpath.i_best if args.partition_index is None else args.partition_index
    print '  found %d clusters in %s' % (len(cpath.partitions[ipartition]), 'best partition' if args.partition_index is None else 'partition at index %d (of %d)' % (ipartition, len(cpath.partitions)))
    if args.cluster_index is None:
        clusters_to_use = cpath.partitions[ipartition]
        print '    taking all %d clusters' % len(clusters_to_use)
    else:
        clusters_to_use = [cpath.partitions[ipartition][args.cluster_index]]
        print '    taking cluster at index %d with size %d' % (args.cluster_index, len(clusters_to_use[0]))
    if args.seed_unique_id is not None:
        clusters_to_use = [c for c in clusters_to_use if args.seed_unique_id in c]  # NOTE can result in more than one cluster with the seed sequence (e.g. if this file contains intermediate annotations from seed partitioning))
        print '    removing clusters not containing sequence \'%s\' (leaving %d)' % (args.seed_unique_id, len(clusters_to_use))

seqfos = []
annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access
for cluster in clusters_to_use:
    if ':'.join(cluster) not in annotations:
        print '  %s cluster with size %d not in annotations, so skipping it' % (utils.color('red', 'warning'), len(cluster))
        continue
    cluster_annotation = annotations[':'.join(cluster)]
    newfos = [{'name' : u, 'seq' : s} for u, s in zip(cluster_annotation['unique_ids'], cluster_annotation['seqs' if args.indel_reversed_seqs else 'input_seqs'])]
    if args.extra_columns is not None:
        for ecol in args.extra_columns:
            if ecol not in cluster_annotation:
                raise Exception('column \'%s\' not found in annotations' % ecol)
            for iseq in range(len(newfos)):
                ival = cluster_annotation[ecol]
                if ecol in utils.linekeys['per_seq']:
                    ival = ival[iseq]
                newfos[iseq][ecol] = ival
    seqfos += newfos

if not os.path.exists(os.path.dirname(os.path.abspath(args.outfile))):
    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)))
print '  writing %d sequences to %s' % (len(seqfos), args.outfile)
with open(args.outfile, 'w') as ofile:
    if utils.getsuffix(args.outfile) in ['.csv', '.tsv']:
        writer = csv.DictWriter(ofile, seqfos[0].keys(), delimiter=',' if utils.getsuffix(args.outfile)=='.csv' else '\t')
        writer.writeheader()
        for sfo in seqfos:
            writer.writerow(sfo)
    elif utils.getsuffix(args.outfile) in ['.fa', '.fasta']:
        for sfo in seqfos:
            estr = ''
            if args.extra_columns is not None:
                estr = args.fasta_info_separator
                estr += args.fasta_info_separator.join(str(sfo[c]) for c in args.extra_columns)
            ofile.write('>%s%s\n%s\n' % (sfo['name'], estr, sfo['seq']))
    else:
        assert False
