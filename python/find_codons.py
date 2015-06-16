#!/usr/bin/env python
import sys
import json
from Bio import SeqIO
import utils
from opener import opener

""" find conserved cysteine and tryptophan in new alleles/genes """


dirname = 'data/imgt'
align_fname = dirname + '/ighv-aligned.fasta'

def get_info(fname):
    info = {}
    for seq_record in SeqIO.parse(fname, 'fasta'):
        info[seq_record.name] = str(seq_record.seq)
    return info

align_info = get_info(align_fname)

with opener('r')(dirname + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
    cyst_positions = json.load(json_file)

align_cpos_aa = 104 - 1  # cysteine position in amino acids (-1 to convert to zero-indexing)
align_cpos = align_cpos_aa * 3

bad_genes = []
for name, align_seq in align_info.items():  # NOTE there's a whole bunch of sequences in seqinfo that aren't in aliinfo. If imgt doesn't have alignment info for 'em, though, I'm ignoring 'em
    print '%-20s' % name,
    # check for unexpected characters
    for pos in align_seq:
        if pos not in utils.nukes + ['.', ]:  #'ACGT.':
            print 'ERROR unexpected character %s in %s from %s' % (pos, name, align_fname)
            sys.exit()

    # see if it's too short (WTF?!?!)
    if align_cpos >= len(align_seq):
        print 'too short!'
        bad_genes.append(name)
        continue
    try:
        utils.check_conserved_cysteine(align_seq, align_cpos, debug=True)
    except:
        bad_genes.append(name)
        continue

    # remove dots
    n_dots = align_seq.count('.')
    real_cpos = align_cpos - n_dots
    utils.check_conserved_cysteine(align_seq.replace('.', ''), real_cpos, debug=True)
    if name in cyst_positions and real_cpos == cyst_positions[name]['cysteine-position']:
        print 'ok'
    else:
        if name in cyst_positions and real_cpos != cyst_positions[name]['cysteine-position']:
            print 'not the same, new: %d old: %s' % (real_cpos, cyst_positions[name]['cysteine-position']),
            print '  switching to the new one'
        else:
            print 'new: %d' % (real_cpos)
        cyst_positions[name] = {}
        cyst_positions[name]['cysteine-position'] = real_cpos

print 'bad genes:'
print '\\|'.join(bad_genes).replace('*', '\\*')
print ''

outfname = 'out.json'
outfile = open(outfname, 'w')
json.dump(cyst_positions, outfile, indent=0)
outfile.close()
print 'wrote output to', outfname
