import bz2
import gzip
import copy
import os
import sys
import csv
from collections import OrderedDict
import random
import re
from Bio import SeqIO
import string
import itertools

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
def translate_columns(line, translations):  # NOTE similar to code in utils.get_arg_list()
    for key in line.keys():
        if key in translations and key != translations[key]:
            line[translations[key]] = line[key]
            del line[key]

# ----------------------------------------------------------------------------------------
def get_more_names(potential_names):
    potential_names += [''.join(ab) for ab in itertools.combinations(potential_names, 2)]

# ----------------------------------------------------------------------------------------
def abbreviate(used_names, potential_names, unique_id):
    if len(used_names) >= len(potential_names):
        get_more_names(potential_names)
    ilet = 0
    new_id = potential_names[ilet]
    while new_id in used_names:  # NOTE this is kind of wasteful, since they're both ordered I could just keep track of which one to use next
        new_id = potential_names[ilet]
        ilet += 1
    used_names.add(new_id)
    return new_id

# ----------------------------------------------------------------------------------------
def get_seqfile_info(args, glfo=None):
    """ return list of sequence info from files of several types """

    # WARNING defaults for <name_column> and <seq_column> also set in partis (since we call this from places other than partis, but we also want people to be able set them from the partis command line)
    internal_name_column = 'unique_id'  # key we use in the internal dictionaries
    internal_seq_column = 'seq'
    name_column = args.name_column
    if name_column is None:  # header we expect in the file
        name_column = internal_name_column
    seq_column = args.seq_column
    if seq_column is None:
        seq_column = internal_seq_column

    if not args.is_data and glfo is None:
        print '  WARNING glfo is None, so not adding implicit info'

    suffix = os.path.splitext(args.infname)[1]
    if len(re.findall('\.[ct]sv', suffix)) > 0:
        if suffix == '.csv':
            delimiter = ','
        elif suffix == '.tsv':
            delimiter = '\t'
        else:
            assert False
        seqfile = opener('r')(args.infname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    else:
        if suffix == '.fasta' or suffix == '.fa':
            ftype = 'fasta'
        elif suffix == '.fastq' or suffix == '.fq':
             ftype = 'fastq'
        else:
            raise Exception('couldn\'t handle file extension for %s' % args.infname)
        reader = []
        n_fasta_queries = 0
        for seq_record in SeqIO.parse(args.infname, ftype):

            # if command line specified query or reco ids, skip other ones (can't have/don't allow simulation info in a fast[aq])
            if args.queries is not None and seq_record.name not in args.queries:
                continue

            reader.append({})
            reader[-1][name_column] = seq_record.name
            reader[-1][seq_column] = str(seq_record.seq).upper()
            n_fasta_queries += 1
            if args.n_max_queries > 0 and n_fasta_queries >= args.n_max_queries:
                break

    input_info = OrderedDict()
    reco_info = None
    if not args.is_data:
        reco_info = OrderedDict()
    n_queries = 0
    found_seed = False
    used_names = set()  # for abbreviating
    if args.abbreviate:
        potential_names = list(string.ascii_lowercase)
    iname = None  # line number -- used as sequence id if there isn't a <name_column>
    for line in reader:
        if seq_column not in line:
            raise Exception('mandatory header \'%s\' not present in %s (you can set column names with --name-column and --seq-column)' % (seq_column, args.infname))
        if name_column not in line and iname is None:
            iname = 0
        if name_column != internal_name_column or seq_column != internal_seq_column:
            if iname is not None:
                line[internal_name_column] = '%09d' % iname
                iname += 1
            translate_columns(line, {name_column : internal_name_column, seq_column: internal_seq_column})
        utils.process_input_line(line)
        unique_id = line[internal_name_column]
        if any(fc in unique_id for fc in utils.forbidden_characters):
            raise Exception('found a forbidden character (one of %s) in sequence id \'%s\' -- sorry, you\'ll have to replace it with something else' % (' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), unique_id))

        if args.abbreviate:
            unique_id = abbreviate(used_names, potential_names, unique_id)

        # if command line specified query or reco ids, skip other ones
        if args.queries is not None and unique_id not in args.queries:
            continue
        if args.reco_ids is not None and line['reco_id'] not in args.reco_ids:
            continue

        if unique_id in input_info:
            raise Exception('found id %s twice in file %s' % (unique_id, args.infname))

        if args.seed_unique_id is not None and unique_id == args.seed_unique_id:
            found_seed = True

        input_info[unique_id] = {'unique_ids' : [unique_id, ], 'seqs' : [line[internal_seq_column], ]}

        if n_queries == 0 and args.is_data and 'v_gene' in line:
            print '  note: found simulation info in %s -- are you sure you didn\'t mean to set --is-simu?' % args.infname

        if not args.is_data:
            if 'v_gene' not in line:
                raise Exception('simulation info not found in %s' % args.infname)
            reco_info[unique_id] = copy.deepcopy(line)
            reco_info[unique_id]['unique_ids'] = [unique_id, ]
            reco_info[unique_id]['seqs'] = [line[internal_seq_column], ]
            reco_info[unique_id]['indelfos'] = [line['indelfo'], ]
            for key in ['unique_id', 'seq', 'indelfo']:
                del reco_info[unique_id][key]
            if glfo is not None:
                utils.add_implicit_info(glfo, reco_info[unique_id], existing_implicit_keys=('cdr3_length', ))

        n_queries += 1
        if args.n_max_queries > 0 and n_queries >= args.n_max_queries:
            break

    if len(input_info) == 0:
        raise Exception('didn\'t end up pulling any input info out of %s while looking for queries: %s reco_ids: %s\n' % (args.infname, str(args.queries), str(args.reco_ids)))
    if args.seed_unique_id is not None and not found_seed:
        raise Exception('couldn\'t find seed %s in %s' % (args.seed_unique_id, args.infname))

    return (input_info, reco_info)
