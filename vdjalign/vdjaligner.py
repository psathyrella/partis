#!/usr/bin/env python
import argparse
import contextlib
import logging
import shutil
import subprocess

from pkg_resources import resource_stream

log = logging.getLogger('vdjalign')

from . import bwa, util
from imgt import igh

@contextlib.contextmanager
def bwa_index_of_package_imgt_fasta(file_name):
    with util.tempdir(prefix='bwa-') as td, \
         resource_stream('vdjalign.imgt.data', file_name) as in_fasta:
        with open(td(file_name), 'w') as fp:
            shutil.copyfileobj(in_fasta, fp)
        cmd = ['bwa', 'index', file_name]
        log.info('Running "%s" in %s', ' '.join(cmd), td())
        subprocess.check_call(cmd, cwd=td())
        yield bwa.load_index(td(file_name))

def identify_v(bwa_index, cysteine_map, s):
    v_alignments = bwa_index.align(s)

    for idx, v_alignment in enumerate(v_alignments):
        cigar = v_alignment.cigar
        qstart = 0
        if cigar[0][1] == 'S':
            qstart = cigar[0][0]
        if cigar[-1][1] == 'S':
            cigar.pop()

        qend = sum(c for c, op in cigar if op != 'D')

        rend = v_alignment.pos + sum(c for c, op in cigar if op not in 'IS')

        cysteine_position = cysteine_map[v_alignment.reference]
        if not cysteine_position:
            if idx == 0:
                log.warn("No Cysteine position for %s", v_alignment.reference)
            continue

        if rend < cysteine_position + 2:
            #if not idx == 0:
                #log.warn('Alignment does not cover Cysteine: (%d < %d) %s|%s', rend, cysteine_position + 2,
                        #s[:qend], s[qend:])
            continue

        cdr3_start = cysteine_position - v_alignment.pos + qstart

        yield {'v_start': qstart, 'v_end': qend,
               'cdr3_start': cdr3_start,
               'v_alignment': v_alignment,
               'frame': (cdr3_start % 3) + 1}

def identify_j(bwa_index, s, seq_start=0, frame=1):
    alignments = bwa_index.align(s[seq_start:])

    for alignment in alignments:
        cigar = alignment.cigar
        qstart = seq_start
        if cigar[0][1] == 'S':
            qstart += cigar[0][0]

        codon_start = qstart
        rem = (codon_start - frame + 1) % 3
        if rem:
            codon_start += 3 - rem

        for i in xrange(codon_start, len(s), 3):
            if s[i:i+3] ==  'TGG':
                cdr3_end = i + 3
                break
        else:
            continue

        yield {'j_start': qstart, 'cdr3_end': cdr3_end, 'j_alignment': alignment}

def process_reads(sequence_iter, v_index, j_index):
    cysteine_map = igh.cysteine_map()
    for name, sequence in sequence_iter:
        v_alignment = next(identify_v(v_index, cysteine_map, sequence), None)
        r = {'name': name, 'sequence': sequence}

        if v_alignment is None:
            r['error'] = 'no V-gene alignments'
            yield r
            continue

        j_alignment = next(identify_j(j_index, sequence, seq_start=v_alignment['v_end'], frame=v_alignment['frame']), None)
        if j_alignment is None:
            r['error'] = 'no J-gene alignments'
            yield r
            continue

        r.update(v_alignment)
        for k, v in j_alignment.iteritems():
            assert k not in r
            r[k] = v
        yield r

        #log.info("%s\n%s\n%s",
                #name,
                #sequence,
                #'{0}{1}{2}'.format('V' * alignment['v_end'],
                    #'N' * (alignment['j_start'] - alignment['v_end']),
                    #'J' * (len(sequence) - alignment['j_start'])))
        #log.info('%08s %s', name, codons.translate(cdr3))
        #log.info("CDR3: L=%d, %s", len(cdr3), cdr3)
        #log.info(alignment)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fastx_file', type=argparse.FileType('r'))
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    indexed = bwa_index_of_package_imgt_fasta
    with indexed('ighv.fasta') as v_index, indexed('ighj.fasta') as j_index, \
         a.fastx_file as ifp:
        j_index.min_seed_len = 5
        sequences = ((n, s) for n, s, _ in util.readfq(ifp))

        result = next(process_reads(sequences, v_index, j_index))
        print list(result['v_alignment'].identify_variations())

if __name__ == '__main__':
    main()
