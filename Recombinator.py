#!/usr/bin/env python
import sys
import random
from Bio.Seq import Seq
from Bio.Alphabet import generic_alphabet

class Recombinator(object):
    seqs = {}
    def __init__(self):
        self.read_file('v', 'igv.txt')
        self.read_file('d', 'igd.txt')
        self.v_seq = self.seqs['v'][random.randint(0,len(self.seqs['v'])-1)]
        self.d_seq = self.seqs['d'][random.randint(0,len(self.seqs['d'])-1)]
        self.erode('right', self.v_seq)
        self.erode('left', self.d_seq)
        self.insert_seq = self.get_insertion()
        self.join(self.v_seq, self.insert_seq, self.d_seq)

    def read_file(self, region, fname):
        infile = open(fname)
        lines = infile.readlines()
        self.seqs[region] = []
        for line in lines:
            self.seqs[region].append(Seq(line.rstrip(), generic_alphabet))

    def erode(self, location, seq):
        print 'eroding'
        n_to_remove = random.randint(0,3) # number of bases to remove
        print 
        print '  ',n_to_remove,'-->',seq
        if location == 'left':
            seq = seq[n_to_remove:len(seq)]
        elif location == 'right':
            seq = seq[0:len(seq)-n_to_remove]
        else:
            print 'no way bub'
            sys.exit()
        print '  ',seq

    def int_to_nucleotide(self, number):
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
        length = random.randint(0,5)
        print 'inserting ',length
        insert_seq_str = ''
        for ic in range(0,length):
            insert_seq_str += self.int_to_nucleotide(random.randint(0,3))
        return Seq(insert_seq_str, generic_alphabet)

    def join(self, left_seq, insert_sequence, right_seq):
        print left_seq,insert_sequence,right_seq
        new_seq = left_seq + insert_sequence + right_seq
        print new_seq

reco = Recombinator()
