#!/usr/bin/env python

import sys
import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_alphabet

class Recombinator(object):
    """ Simulates the process of VDJ recombination """ 
    
    seqs = {}
    def __init__(self):
        """ Initialize from files """ 
        self.read_file('v', 'ighv.fasta')
        self.read_file('d', 'ighd.fasta')
        self.read_file('j', 'ighd.fasta')

    def combine(self):
        """ Run the combination """
        self.v_seq = self.choose_seq('v')
        self.d_seq = self.choose_seq('d')
        self.j_seq = self.choose_seq('j')
        self.erode('right', self.v_seq, 'v')
        self.erode('left', self.d_seq, 'd')
        self.erode('right', self.d_seq, 'd')
        self.erode('left', self.j_seq, 'j')
        vd_insertion = self.get_insertion()
        vd_seq = self.join(self.v_seq, vd_insertion, self.d_seq)
        dj_insertion = self.get_insertion()
        self.join(vd_seq, dj_insertion, self.j_seq)

    def read_file(self, region, fname):
        """ Read the various germline variants from file and store as
        Seq objects """
        self.seqs[region] = []
        for seq_record in SeqIO.parse(fname, "fasta"):
            self.seqs[region].append(seq_record.seq)

    def choose_seq(self, region):
        """ Choose which of the germline variants to use """
        i_seq = random.randint(0, len(self.seqs[region])-1)
        print ('choosing %s'
               'sequence number %d (of %d)') % (region,
                                                          i_seq,
                                                          len(self.seqs[region]))
        return self.seqs[region][i_seq]

    def erode(self, location, seq, region):
        """ Erode (delete) some number of letters from the <location> side of
        <seq>, where location is 'left' or 'right' """
        n_to_remove = random.randint(0,3) # number of bases to remove
        fragment_before = ''
        fragment_after = ''
        if location == 'left':
            fragment_before = seq[:n_to_remove + 3]
            seq = seq[n_to_remove:len(seq)]
            fragment_after = seq[:n_to_remove + 3]
        elif location == 'right':
            fragment_before = seq[len(seq) - n_to_remove - 3 :]
            seq = seq[0:len(seq)-n_to_remove]
            fragment_after = seq[len(seq) - n_to_remove - 3 :]
        else:
            print 'no way bub'
            sys.exit()
        tail_after = seq[len(seq) - n_to_remove - 3 : ]
        print 'eroding %3d from %5s of %s region: %15s --> %-15s' % (n_to_remove, location, region, fragment_before, fragment_after)

    def int_to_nucleotide(self, number):
        """ Convert between (0,1,2,3) and (A,C,G,T) """
        if number == 0:
            return 'A'
        elif number == 1:
            return 'C'
        elif number == 2:
            return 'G'
        elif number == 3:
            return 'T'
        else:
            print 'ERROR bad number'
            sys.exit()

    def get_insertion(self):
        """ Get the non-templated sequence to insert between
        templated regions """
        length = random.randint(0,5)
        insert_seq_str = ''
        for ic in range(0,length):
            insert_seq_str += self.int_to_nucleotide(random.randint(0,3))
        print 'inserting %d: %s' % (length, insert_seq_str)
        return Seq(insert_seq_str, generic_alphabet)

    def join(self, left_seq, insert_sequence, right_seq):
        print 'joining: %s + %s + %s' % (left_seq, insert_sequence, right_seq)
        new_seq = left_seq + insert_sequence + right_seq
        print new_seq
        return new_seq

if __name__ == '__main__':
    reco = Recombinator()
    reco.combine()
