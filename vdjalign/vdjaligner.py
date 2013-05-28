#!/usr/bin/env python
import argparse
import collections
import contextlib
import logging
import shutil
import subprocess

from pkg_resources import resource_stream

log = logging.getLogger('vdjalign')

from . import bwa, util

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

def main():
    p = argparse.ArgumentParser()
    p.add_argument('fastx_file', type=argparse.FileType('r'))
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    indexed = bwa_index_of_package_imgt_fasta
    with indexed('ighv.fasta') as v_index, indexed('ighj.fasta') as j_index:
        with a.fastx_file as fp:
            sequences = util.readfq(fp)

            c = collections.Counter()
            for name, sequence, _ in sequences:
                v_alignments = v_index.align(sequence)
                if not v_alignments:
                    logging.warn('No V alignments for %s', name)
                    continue

                best_v = v_alignments[0]
                best_v_cigar = best_v['cigar']
                v_qstart = 0
                if best_v_cigar[0][1] == 'S':
                    v_qstart = best_v_cigar[0][0]
                if  best_v_cigar[-1][1] == 'S':
                    best_v_cigar.pop()

                v_qend = sum(c for c, op in best_v['cigar'] if op != 'D')

                v_rend = best_v['pos'] + sum(c for c, op in best_v['cigar'] if op not in 'IS')

                ## TODO: should be able to determine cyteine location from sequence
                #try:
                    #v_qend = max(sequence.rindex('TGT', 0, v_qend), sequence.rindex('TGC', 0, v_qend)) + 2
                #except ValueError:
                    #logging.warn('No Cysteine')
                    #continue
                #assert v_qend > 0

                j_alignments = j_index.align(sequence[v_qend:])
                if not j_alignments:
                    logging.warn('No J alignments for %s', name)
                    continue

                best_j = j_alignments[0]
                best_j_cigar = best_j['cigar']
                j_qstart = v_qend
                if best_j_cigar[0][1] == 'S':
                    j_qstart += best_j_cigar[0][0]

                #try:
                    #j_qstart = sequence.index('TGG', j_qstart)
                #except ValueError:
                    #logging.warn('No Tryptophan')
                    #continue

                c[best_j['reference'].split('*', 1)[0], j_qstart - v_qend] += 1

                print sequence
                print '{0}{1}{2}{3}'.format('.' * v_qstart,
                        v_index.fetch_reference(best_v['rid'], best_v['pos'], v_rend).lower(),
                        '.' * (j_qstart - v_qend),
                        j_index.fetch_reference(best_j['rid'], best_j['pos'], best_j['pos'] + sum(c for c, op in best_j_cigar if op not in 'IS')).lower())

        logging.info(c)

if __name__ == '__main__':
    main()
