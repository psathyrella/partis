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
def add_seed_seq(args, input_info, reco_info, is_data):
    input_info[args.seed_unique_id] = {'unique_ids' : [args.seed_unique_id, ], 'seqs' : [args.seed_seq, ]}
    if not is_data:
        reco_info[args.seed_unique_id] = 'unknown!'  # hopefully more obvious than a key error

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
def get_seqfile_info(infname, is_data, n_max_queries=-1, args=None, glfo=None, simglfo=None):
    """ return list of sequence info from files of several types """

    # WARNING defaults for <name_column> and <seq_column> also set in partis (since we call this from places other than partis, but we also want people to be able set them from the partis command line)
    internal_name_column = 'unique_id'  # key we use in the internal dictionaries
    internal_seq_column = 'seq'

    name_column = internal_name_column
    seq_column = internal_seq_column
    if args is not None:
        if args.name_column is not None:
            name_column = args.name_column
        if args.seq_column is not None:
            seq_column = args.seq_column

    if not is_data and glfo is None:
        print '  WARNING glfo is None, so not adding implicit info'

    suffix = os.path.splitext(infname)[1]
    if len(re.findall('\.[ct]sv', suffix)) > 0:
        if suffix == '.csv':
            delimiter = ','
        elif suffix == '.tsv':
            delimiter = '\t'
        else:
            assert False
        seqfile = opener('r')(infname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    else:
        if suffix == '.fasta' or suffix == '.fa':
            ftype = 'fasta'
        elif suffix == '.fastq' or suffix == '.fq':
             ftype = 'fastq'
        else:
            raise Exception('couldn\'t handle file extension for %s' % infname)
        reader = []
        n_fasta_queries = 0
        for seq_record in SeqIO.parse(infname, ftype):

            # if command line specified query or reco ids, skip other ones (can't have/don't allow simulation info in a fast[aq])
            if args is not None and args.queries is not None and seq_record.name not in args.queries:
                continue

            reader.append({})
            reader[-1][name_column] = seq_record.name
            reader[-1][seq_column] = str(seq_record.seq).upper()
            n_fasta_queries += 1
            if n_max_queries > 0 and n_fasta_queries >= n_max_queries:
                break

    input_info = OrderedDict()
    reco_info = None
    if not is_data:
        reco_info = OrderedDict()
    already_printed_forbidden_character_warning = False
    n_queries_added = 0
    found_seed = False
    used_names = set()  # for abbreviating
    if args is not None and args.abbreviate:
        potential_names = list(string.ascii_lowercase)
    iname = None  # line number -- used as sequence id if there isn't a <name_column>
    iline = -1
    for line in reader:
        iline += 1
        if args is not None and args.istartstop is not None:
            if iline < args.istartstop[0]:
                continue
            if iline >= args.istartstop[1]:
                break
        if seq_column not in line:
            raise Exception('mandatory header \'%s\' not present in %s (you can set column names with --name-column and --seq-column)' % (seq_column, infname))
        if name_column not in line and iname is None:
            iname = 0
        if iname is not None:
            line[internal_name_column] = '%09d' % iname
            iname += 1
        if name_column != internal_name_column or seq_column != internal_seq_column:
            translate_columns(line, {name_column : internal_name_column, seq_column: internal_seq_column})
        utils.process_input_line(line)
        unique_id = line[internal_name_column]
        if any(fc in unique_id for fc in utils.forbidden_characters):
            if not already_printed_forbidden_character_warning:
                print '  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (utils.color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), unique_id, unique_id.translate(utils.forbidden_character_translations))
                already_printed_forbidden_character_warning = True
            unique_id = unique_id.translate(utils.forbidden_character_translations)
        if args is not None and args.abbreviate:
            unique_id = abbreviate(used_names, potential_names, unique_id)

        # if command line specified query or reco ids, skip other ones
        if args is not None:
            if args.queries is not None and unique_id not in args.queries:
                continue
            if args.reco_ids is not None and line['reco_id'] not in args.reco_ids:
                continue
            if args.seed_unique_id is not None and unique_id == args.seed_unique_id:
                found_seed = True

        if unique_id in input_info:
            raise Exception('found id %s twice in file %s' % (unique_id, infname))

        if len(line[internal_seq_column].translate(None, ''.join(utils.alphabet))) > 0:
            raise Exception('unexpected character (not among %s) in input sequence with id %s:\n  %s' % (utils.nukes + utils.ambiguous_bases, unique_id, line[internal_seq_column]))

        input_info[unique_id] = {'unique_ids' : [unique_id, ], 'seqs' : [line[internal_seq_column], ]}

        if n_queries_added == 0 and is_data and 'v_gene' in line:
            print '  note: found simulation info in %s -- are you sure you didn\'t mean to set --is-simu?' % infname

        if not is_data:
            if 'v_gene' not in line:
                raise Exception('simulation info not found in %s' % infname)
            reco_info[unique_id] = copy.deepcopy(line)
            reco_info[unique_id]['unique_ids'] = [unique_id, ]
            reco_info[unique_id]['seqs'] = [line[internal_seq_column], ]
            reco_info[unique_id]['indelfos'] = [line['indelfo'], ]
            for key in ['unique_id', 'seq', 'indelfo']:
                del reco_info[unique_id][key]
            if simglfo is not None:
                utils.add_implicit_info(simglfo, reco_info[unique_id])

        n_queries_added += 1
        if n_max_queries > 0 and n_queries_added >= n_max_queries:
            break

    if args is not None:
        if args.istartstop is not None:
            n_lines_in_file = iline + 1
            if n_lines_in_file < args.istartstop[1]:
                raise Exception('--istartstop upper bound %d larger than number of lines in file %d' % (args.istartstop[1], n_lines_in_file))
        if len(input_info) == 0:
            raise Exception('didn\'t find the specified --queries (%s) or --reco-ids (%s) in %s' % (str(args.queries), str(args.reco_ids), infname))
        if args.queries is not None:
            missing_queries = set(args.queries) - set(input_info)
            extra_queries = set(input_info) - set(args.queries)  # this is just checking for a bug in the code just above here...
            if len(missing_queries) > 0:
                raise Exception('didn\'t find some of the specified --queries: %s' % ' '.join(missing_queries))
            if len(extra_queries) > 0:
                raise Exception('extracted uids %s that weren\'t specified with --queries' % ' '.join(extra_queries))
        if args.seed_unique_id is not None:
            if found_seed:
                if args.seed_seq is not None:  # and input_info[args.seed_unique_id]['seqs'][0] != args.seed_seq:
                    # raise Exception('incompatible --seed-unique-id and --seed-seq (i.e. the sequence in %s corresponding to %s wasn\'t %s)' % (infname, args.seed_unique_id, args.seed_seq))
                    raise Exception('--seed-seq was specified, but --seed-unique-id was also present in input file')
            else:
                if args.seed_seq is None:
                    raise Exception('couldn\'t find seed unique id %s in %s' % (args.seed_unique_id, infname))
                add_seed_seq(args, input_info, reco_info, is_data)
        elif args.seed_seq is not None:
            args.seed_unique_id = 'seed-seq'
            add_seed_seq(args, input_info, reco_info, is_data)
        elif args.random_seed_seq:  # already checked (in bin/partis) that other seed args aren't set
            args.seed_unique_id = random.choice(input_info.keys())
            print '    chose random seed unique id %s' % args.seed_unique_id

    return input_info, reco_info
