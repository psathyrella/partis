import bz2
import gzip
import os.path
import sys
import csv
from collections import OrderedDict
import random
from Bio import SeqIO

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
def get_seqfile_info(fname, is_data, germline_seqs=None, cyst_positions=None, tryp_positions=None, n_max_queries=-1, queries=None, reco_ids=None, randomize_order=False, replace_N_with=None):
    """ return list of sequence info from files of several types """
    if replace_N_with is not None:
        print 'WARNING replacing N with %s' % replace_N_with
    if is_data:
        assert not randomize_order  # only really makes sense to randomize simulation
    else:
        assert germline_seqs is not None
        assert cyst_positions is not None
        assert tryp_positions is not None

    if '.csv' in fname:
        delimiter = ','
        name_column = 'unique_id'
        seq_column = 'seq'
        seqfile = opener('r')(fname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    elif '.tsv' in fname:
        delimiter = '\t'
        name_column = 'name'
        seq_column = 'nucleotide'
        seqfile = opener('r')(fname)
        reader = csv.DictReader(seqfile, delimiter=delimiter)
    elif '.fasta' in fname or '.fa' in fname or '.fastq' in fname or '.fq' in fname:
        name_column = 'unique_id'
        seq_column = 'seq'
        reader = []
        n_fasta_queries = 0
        ftype = 'fasta' if ('.fasta' in fname or '.fa' in fname) else 'fastq'
        for seq_record in SeqIO.parse(fname, ftype):
            reader.append({})
            reader[-1][name_column] = seq_record.name
            reader[-1][seq_column] = str(seq_record.seq).upper()
            n_fasta_queries += 1
            if n_max_queries > 0 and n_fasta_queries >= n_max_queries:
                break
    else:
        raise Exception('unrecognized file format %s' % fname)

    input_info = OrderedDict()
    reco_info = None
    if not is_data:
        reco_info = OrderedDict()
    n_queries = 0
    namelist = []  # only used for randomization
    for line in reader:
        utils.process_input_line(line, int_columns=('v_5p_del', 'd_5p_del', 'cdr3_length', 'j_5p_del', 'j_3p_del', 'd_3p_del', 'v_3p_del'))
        # if command line specified query or reco ids, skip other ones
        if queries is not None and line[name_column] not in queries:
            continue
        if reco_ids is not None and line['reco_id'] not in reco_ids:
            continue

        if replace_N_with is not None:
            line[seq_column] = line[seq_column].replace('N', replace_N_with)

        input_info[line[name_column]] = {'unique_id':line[name_column], 'seq':line[seq_column]}
        if not is_data:
            reco_info[line['unique_id']] = line
            utils.add_match_info(germline_seqs, line, cyst_positions, tryp_positions)
        if randomize_order:
            namelist.append(line[name_column])
        n_queries += 1
        if n_max_queries > 0 and n_queries >= n_max_queries:
            break

    if len(input_info) == 0:
        raise Exception('didn\'t end up pulling any input info out of %s while looking for queries: %s reco_ids: %s\n' % (fname, str(queries), str(reco_ids)))
    
    # for k in reco_info.keys():
    #     print reco_info[k]['reco_id']
    if randomize_order:
        rand_input_info, rand_reco_info = OrderedDict(), OrderedDict()
        while len(namelist) > 0:
            irand = random.randint(0, len(namelist) - 1)  # NOTE interval is inclusive
            uid = namelist[irand]
            rand_input_info[uid] = input_info[uid]
            rand_reco_info[uid] = reco_info[uid]
            namelist.remove(uid)

        input_info, reco_info = rand_input_info, rand_reco_info
    # print '---'
    # for k in reco_info.keys():
    #     print reco_info[k]['reco_id']

    return (input_info, reco_info)
