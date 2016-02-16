import bz2
import gzip
import os
import sys
import csv
from collections import OrderedDict
import random
from Bio import SeqIO

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
def get_seqfile_info(fname, is_data, glfo=None, n_max_queries=-1, queries=None, reco_ids=None):
    """ return list of sequence info from files of several types """

    suffix = os.path.splitext(fname)[1]
    if suffix == '.csv':
        delimiter = ','
        name_column = 'unique_id'
        seq_column = 'seq'
        seqfile = opener('r')(fname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    elif suffix == '.tsv':
        delimiter = '\t'
        name_column = 'name'
        seq_column = 'nucleotide'
        seqfile = opener('r')(fname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    else:
        if suffix == '.fasta' or suffix == '.fa':
            ftype = 'fasta'
        elif suffix == '.fastq' or suffix == '.fq':
             ftype = 'fastq'
        else:
            raise Exception('couldn\'t handle file extension for %s' % fname)
        name_column = 'unique_id'
        seq_column = 'seq'
        reader = []
        n_fasta_queries = 0
        for seq_record in SeqIO.parse(fname, ftype):

            # if command line specified query or reco ids, skip other ones
            if queries is not None and seq_record.name not in queries:
                continue
            # if reco_ids is not None and line['reco_id'] not in reco_ids:  # probably no reco ids in a fasta file
            #     continue

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
    for line in reader:
        if '.csv' in fname and name_column not in line:  # hackey hackey hackey
            name_column = 'name'
            seq_column = 'nucleotide'
        utils.process_input_line(line)
        unique_id = line[name_column]
        if ':' in unique_id:
            raise Exception('found a \':\' in sequence id \'%s\' -- you\'ll have to replace it with something else, as we use \':\'s internally to concatenate sequence ids' % unique_id)

        # if command line specified query or reco ids, skip other ones
        if queries is not None and unique_id not in queries:
            continue
        if reco_ids is not None and line['reco_id'] not in reco_ids:
            continue

        input_info[unique_id] = {'unique_id' : unique_id, 'seq' : line[seq_column]}
        if not is_data:
            if 'v_gene' not in line:
                raise Exception('simulation info not found in %s -- if this is data add option --is-data' % fname)
            reco_info[unique_id] = dict(line)
            if 'indels' in line and line['indels']['reversed_seq'] != '':  # TODO unhackify this
                reco_info[unique_id]['seq'] = line['indels']['reversed_seq']
            if 'indels' not in line:  # TODO unhackify this
                reco_info[unique_id]['indels'] = None
            if glfo is not None:
                utils.remove_implicit_info(reco_info[unique_id], multi_seq=False)
                utils.add_implicit_info(glfo, reco_info[unique_id], multi_seq=False)  # each seq is on its own line in the file
        n_queries += 1
        if n_max_queries > 0 and n_queries >= n_max_queries:
            break

    if len(input_info) == 0:
        raise Exception('didn\'t end up pulling any input info out of %s while looking for queries: %s reco_ids: %s\n' % (fname, str(queries), str(reco_ids)))
    
    return (input_info, reco_info)
