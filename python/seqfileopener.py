import bz2
import gzip
import copy
import os
import sys
import csv
from collections import OrderedDict
import random
import re
import string
import itertools

import utils

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
        seqfile = open(infname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    else:
        reader = []
        n_fasta_queries = 0
        already_printed_forbidden_character_warning = False
        for seqinfo in utils.read_fastx(infname):
            uid = seqinfo['name']

            # if command line specified query or reco ids, skip other ones (can't have/don't allow simulation info in a fast[aq])
            if args is not None and args.queries is not None and uid not in args.queries:
                continue

            reader.append({})

            if any(fc in uid for fc in utils.forbidden_characters):
                if not already_printed_forbidden_character_warning:
                    print '  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (utils.color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), uid, uid.translate(utils.forbidden_character_translations))
                    already_printed_forbidden_character_warning = True
                uid = uid.translate(utils.forbidden_character_translations)

            reader[-1]['unique_ids'] = uid
            reader[-1]['input_seqs'] = str(seqinfo['seq']).upper()
            n_fasta_queries += 1
            if n_max_queries > 0 and n_fasta_queries >= n_max_queries:
                break

    input_info = OrderedDict()
    reco_info = None
    if not is_data:
        reco_info = OrderedDict()
    # already_printed_forbidden_character_warning = False
    n_queries_added = 0
    found_seed = False
    used_names = set()  # for abbreviating
    if args is not None and args.abbreviate:
        potential_names = list(string.ascii_lowercase)
    iname = None  # line number -- used as sequence id if there isn't a name column in the file
    iline = -1
    for line in reader:
        iline += 1
        if args is not None:
            if args.istartstop is not None:
                if iline < args.istartstop[0]:
                    continue
                if iline >= args.istartstop[1]:
                    break
            if args.name_column is not None:
                line['unique_ids'] = line[args.name_column]
                del line[args.name_column]
            if args.seq_column is not None:
                line['input_seqs'] = line[args.seq_column]
                if args.seq_column != 'seqs':  # stupid god damn weird backwards compatibility edge case bullshit
                    del line[args.seq_column]
        if iname is None and 'unique_ids' not in line and 'unique_id' not in line:
            print '  %s: couldn\'t find a name (unique id) column, so using line number as the sequence label (you can set the name column with --name-column)' % (utils.color('yellow', 'warning'))
            iname = 0
        if iname is not None:
            line['unique_ids'] = '%09d' % iname
            iname += 1
        if 'input_seqs' not in line and 'seq' not in line:
            raise Exception('couldn\'t find a sequence column in %s (you can set this with --seq-column)' % infname)
        utils.process_input_line(line)
        if len(line['unique_ids']) > 1:
            raise Exception('can\'t yet handle multi-seq csv input files')
        uid = line['unique_ids'][0]
        inseq = line['input_seqs'][0]

        # NOTE I just moved this to the .fa loop, since otherwise we have no way of knowing how to interpret special characters... nevertheless if someone passesin a csv with special characters as part of a uid this will break
        # if any(fc in uid for fc in utils.forbidden_characters):
        #     if not already_printed_forbidden_character_warning:
        #         print '  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (utils.color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), uid, uid.translate(utils.forbidden_character_translations))
        #         already_printed_forbidden_character_warning = True
        #     uid = uid.translate(utils.forbidden_character_translations)
        if args is not None:
            if args.abbreviate:  # note that this changes <uid>, but doesn't modify <line>
                uid = abbreviate(used_names, potential_names, uid)
            if args.queries is not None and uid not in args.queries:
                continue
            if args.reco_ids is not None and line['reco_id'] not in args.reco_ids:
                continue
            if args.seed_unique_id is not None and uid == args.seed_unique_id:
                found_seed = True

        if uid in input_info:
            raise Exception('found uid \'%s\' twice in input file %s' % (uid, infname))

        if len(inseq.translate(None, ''.join(utils.alphabet))) > 0:
            raise Exception('unexpected character (not among %s) in input sequence with id %s:\n  %s' % (utils.nukes + utils.ambiguous_bases, uid, inseq))

        input_info[uid] = {'unique_ids' : [uid, ], 'seqs' : [inseq, ]}

        if n_queries_added == 0 and is_data and 'v_gene' in line:
            print '  note: found simulation info in %s -- are you sure you didn\'t mean to set --is-simu?' % infname

        if not is_data:
            if 'v_gene' not in line:
                raise Exception('simulation info not found in %s' % infname)
            reco_info[uid] = copy.deepcopy(line)
            if simglfo is not None:
                utils.add_implicit_info(simglfo, reco_info[uid])

        n_queries_added += 1
        if n_max_queries > 0 and n_queries_added >= n_max_queries:
            break

    if args is not None:
        if args.istartstop is not None:
            n_lines_in_file = iline + 1
            if n_lines_in_file < args.istartstop[1]:
                raise Exception('--istartstop upper bound %d larger than number of lines in file %d' % (args.istartstop[1], n_lines_in_file))
        if len(input_info) == 0:
            if args.queries is not None:
                raise Exception('didn\'t find the specified --queries (%s) in %s' % (str(args.queries), infname))
            if args.reco_ids is not None:
                raise Exception('didn\'t find the specified --reco-ids (%s) in %s' % (str(args.reco_ids), infname))
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

    if len(input_info) == 0:
        raise Exception('didn\'t read any sequences from %s' % infname)

    return input_info, reco_info
