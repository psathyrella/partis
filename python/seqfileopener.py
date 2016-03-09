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
def abbreviate(used_names, potential_names, unique_id):
    ilet = 0
    new_id = potential_names[ilet]
    while new_id in used_names:  # NOTE this is kind of wasteful, since they're both ordered I could just keep track of which one to use next
        new_id = potential_names[ilet]
        ilet += 1
    used_names.add(new_id)
    return new_id

# ----------------------------------------------------------------------------------------
def get_seqfile_info(fname, is_data, glfo=None, n_max_queries=-1, queries=None, reco_ids=None, name_column=None, seq_column=None, seed_unique_id=None, abbreviate_names=False):
    """ return list of sequence info from files of several types """

    # WARNING defaults for <name_column> and <seq_column> also set in partis.py (since we call this from places other than partis.py, but we also want people to be able set them from the partis.py command line)
    internal_name_column = 'unique_id'  # key we use in the internal dictionaries
    internal_seq_column = 'seq'
    if name_column is None:  # header we expect in the file
        name_column = internal_name_column
    if seq_column is None:
        seq_column = internal_seq_column

    if not is_data and glfo is None:
        print '  WARNING glfo is None, so not adding implicit info'

    suffix = os.path.splitext(fname)[1]
    if len(re.findall('\.[ct]sv', suffix)) > 0:
        if suffix == '.csv':
            delimiter = ','
        elif suffix == '.tsv':
            delimiter = '\t'
        else:
            assert False
        seqfile = opener('r')(fname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    else:
        if suffix == '.fasta' or suffix == '.fa':
            ftype = 'fasta'
        elif suffix == '.fastq' or suffix == '.fq':
             ftype = 'fastq'
        else:
            raise Exception('couldn\'t handle file extension for %s' % fname)
        reader = []
        n_fasta_queries = 0
        for seq_record in SeqIO.parse(fname, ftype):

            # if command line specified query or reco ids, skip other ones (can't have/don't allow simulation info in a fast[aq])
            if queries is not None and seq_record.name not in queries:
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
    n_queries = 0
    found_seed = False
    used_names = set()  # for abbreviating
    if abbreviate_names:
        potential_names = list(string.ascii_lowercase)
        while n_max_queries > len(potential_names):
            potential_names += [''.join(ab) for ab in itertools.combinations(potential_names, 2)]
    for line in reader:
        if name_column not in line or seq_column not in line:
            raise Exception('mandatory headers \'%s\' and \'%s\' not both present in %s (set with --name-column and --seq-column)' % (name_column, seq_column, fname))
        if name_column != internal_name_column or seq_column != internal_seq_column:
            translate_columns(line, {name_column : internal_name_column, seq_column: internal_seq_column})
        utils.process_input_line(line)
        unique_id = line[internal_name_column]
        if abbreviate_names:
            unique_id = abbreviate(used_names, potential_names, unique_id)
        if any(fc in unique_id for fc in utils.forbidden_characters):
            raise Exception('found a forbidden character (one of %s) in sequence id \'%s\' -- sorry, you\'ll have to replace it with something else' % (' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), unique_id))

        # if command line specified query or reco ids, skip other ones
        if queries is not None and unique_id not in queries:
            continue
        if reco_ids is not None and line['reco_id'] not in reco_ids:
            continue

        if unique_id in input_info:
            raise Exception('found id %s twice in file %s' % (unique_id, fname))

        if seed_unique_id is not None and unique_id == seed_unique_id:
            found_seed = True

        input_info[unique_id] = {'unique_id' : unique_id, 'seq' : line[internal_seq_column]}

        if n_queries == 0 and is_data and 'v_gene' in line:
            print 'WARNING found simulation info in %s -- are you sure you didn\'t mean to set --is-simu?' % fname

        if not is_data:
            if 'v_gene' not in line:
                raise Exception('simulation info not found in %s' % fname)
            reco_info[unique_id] = copy.deepcopy(line)
            reco_info[unique_id]['unique_id'] = unique_id  # in case we're abbreviating
            if glfo is not None:
                utils.add_implicit_info(glfo, reco_info[unique_id], multi_seq=False, existing_implicit_keys=('cdr3_length', ))  # single seqs, since each seq is on its own line in the file

        n_queries += 1
        if n_max_queries > 0 and n_queries >= n_max_queries:
            break

    if len(input_info) == 0:
        raise Exception('didn\'t end up pulling any input info out of %s while looking for queries: %s reco_ids: %s\n' % (fname, str(queries), str(reco_ids)))
    if seed_unique_id is not None and not found_seed:
        raise Exception('couldn\'t find seed %s in %s' % (seed_unique_id, fname))
    
    return (input_info, reco_info)
