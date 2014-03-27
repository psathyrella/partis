#!/usr/bin/env python
""" Simulates the process of VDJ recombination """ 
import sys
import random
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_alphabet

def int_to_nucleotide(number):
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
        print 'ERROR nucleotide number not in [0,4]'
        sys.exit()

class Recombinator(object):
    """ Simulates the process of VDJ recombination """ 
    all_seqs = {}  # all the Vs, all the Ds...
    seqs = {}  # one of each
    regions = ['v', 'd', 'j']
    def __init__(self):
        """ Initialize from files """ 
        for region in self.regions:
            self.read_file(region, 'igh'+region+'.fasta')

    def combine(self):
        """ Run the combination """
        for region in self.regions:
            self.seqs[region] = self.choose_seq(region)
        self.erode('right', self.seqs['v'], 'v')
        self.erode('left', self.seqs['d'], 'd')
        self.erode('right', self.seqs['d'], 'd')
        self.erode('left', self.seqs['j'], 'j')
        vd_seq = self.join(self.seqs['v'], self.get_insertion(), self.seqs['d'])
        self.join(vd_seq, self.get_insertion(), self.seqs['j'])

    def read_file(self, region, fname):
        """ Read the various germline variants from file and store as
        Seq objects
        """
        self.all_seqs[region] = []
        for seq_record in SeqIO.parse(fname, "fasta"):
            self.all_seqs[region].append(seq_record.seq)

    def choose_seq(self, region):
        """ Choose which of the germline variants to use """
        i_seq = random.randint(0, len(self.all_seqs[region])-1)
        print 'choosing %s' % region,
        print 'sequence number %d (of %d)' % (i_seq, len(self.all_seqs[region]))
        return self.all_seqs[region][i_seq]

    def get_n_to_erode(self):
        """ Number of bases to remove before joining to the neighboring
        sequence
        """
        # NOTE: casting to an int reduces the mean somewhat
        return int(numpy.random.exponential(4))

    def erode(self, location, seq, region):
        """ Erode (delete) some number of letters from the <location> side of
        <seq>, where location is 'left' or 'right'
        """
        n_to_erode = self.get_n_to_erode()
        fragment_before = ''
        fragment_after = ''
        if location == 'left':
            fragment_before = seq[:n_to_erode + 3]
            seq = seq[n_to_erode:len(seq)]
            fragment_after = seq[:n_to_erode + 3]
        elif location == 'right':
            fragment_before = seq[len(seq) - n_to_erode - 3 :]
            seq = seq[0:len(seq)-n_to_erode]
            fragment_after = seq[len(seq) - n_to_erode - 3 :]
        else:
            print 'ERROR location must be \"left\" or \"right\"'
            sys.exit()
        print 'eroding %3d from %5s' % (n_to_erode, location)
        print ' of %s region: %15s' % (region, fragment_before)
        print ' --> %-15s' % fragment_after

    def get_insertion(self):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        length = int(numpy.random.exponential(6))
        insert_seq_str = ''
        for _ in range(0, length):
            insert_seq_str += int_to_nucleotide(random.randint(0, 3))
        print 'inserting %d: %s' % (length, insert_seq_str)
        return Seq(insert_seq_str, generic_alphabet)

    def join(self, left_seq, insert_sequence, right_seq):
        """ Join three sequences in the indicated order """
        new_seq = left_seq + insert_sequence + right_seq
        print 'joining: %s + %s + %s' % (left_seq, insert_sequence, right_seq)
        print 'gives: %s' % new_seq
        return new_seq

if __name__ == '__main__':
    reco = Recombinator()
    reco.combine()
