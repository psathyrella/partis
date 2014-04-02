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
        print 'ERROR nucleotide number not in [0,3]'
        sys.exit()

class Recombinator(object):
    """ Simulates the process of VDJ recombination """

    # parameters that control recombination, erosion, and whatnot
    mean_nukes_eroded = 4  # Mean number of nucleotides to remove at each side of the boundary.
                           # NOTE actual mean is a bit less than this because I round to an integer
    mute_rate = 0.06  # average number of point mutations per base
    mean_insertion_length = 6  # mean length of the non-templated insertion

    all_seqs = {}  # all the Vs, all the Ds...
    seqs = {}   # one of each
    regions = ['v', 'd', 'j']
    version_freq_table = {}

    def __init__(self):
        """ Initialize from files """ 
        for region in self.regions:
            self.read_vdj_versions(region, 'igh'+region+'.fasta')
        self.read_vdj_version_freqs('version-counter/filtered-vdj-probs.txt')

    def combine(self):
        """ Run the combination. Returns the new sequence. """
        vdj_combo_label = self.choose_vdj_combo()  # a tuple with the names of the chosen versions
        for region in self.regions:
            self.seqs[region] = self.get_seq(region, vdj_combo_label)
        print '  eroding'
        self.erode('right', self.seqs['v'], 'v')
        self.erode('left', self.seqs['d'], 'd')
        self.erode('right', self.seqs['d'], 'd')
        self.erode('left', self.seqs['j'], 'j')
        vd_seq = self.join(self.seqs['v'], self.get_insertion(), self.seqs['d'])
        recombined_seq = self.join(vd_seq, self.get_insertion(), self.seqs['j'])
        final_seq = self.mutate(recombined_seq)
        print '  before mute:',recombined_seq
        print '  after mute: ',final_seq
        return final_seq

    def read_vdj_versions(self, region, fname):
        """ Read the various germline variants from file and store as
        Seq objects
        """
        self.all_seqs[region] = {}
        for seq_record in SeqIO.parse(fname, "fasta"):
            self.all_seqs[region][seq_record.name] = seq_record.seq

    def read_vdj_version_freqs(self, fname):
        """ Read the frequencies at which various VDJ combinations appeared
        in data. This file was created with the code in version-counter/
        """
        with open(fname) as infile:
            for line in infile:
                if line[0]=='#':
                    continue
                line_fragments = line.split()
                v_gene = str(line_fragments[0])
                d_gene = str(line_fragments[1])
                j_gene = str(line_fragments[2])
                cdr3_length  = int(line_fragments[3])
                prob = float(line_fragments[4])
                # TODO add some input sanitization here
                self.version_freq_table[(v_gene, d_gene, j_gene, cdr3_length)] = prob

    def choose_vdj_combo(self):
        """ Choose which combination germline variants to use """
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        print '  choosing combination',
        for vdj_choice in self.version_freq_table:
#            print '  += %.2e for %60s' % (self.version_freq_table[vdj_choice], vdj_choice),
            sum_prob += self.version_freq_table[vdj_choice]
#            print ' %.3e < %.3e ? %d' % (iprob, sum_prob, iprob < sum_prob)
            if iprob < sum_prob:
                print ' ---> ',vdj_choice
                return vdj_choice
        print 'ERROR I shouldn\'t be here'
        sys.exit()  # shouldn't get here!

    def get_seq(self, region, vdj_combo_label):
        gene_name = ''
        if region == 'v':
            gene_name = vdj_combo_label[0]
        elif region == 'd':
            gene_name = vdj_combo_label[1]
        elif region == 'j':
            gene_name = vdj_combo_label[2]
        else:
            print 'ERROR that don\'t make no sense'
            sys.exit()
        print '    returning %s: %s' % (region, gene_name)
        return self.all_seqs[region][gene_name]

    def get_n_to_erode(self):
        """ Number of bases to remove before joining to the neighboring
        sequence
        """
        # NOTE casting to an int reduces the mean somewhat
        return int(numpy.random.exponential(self.mean_nukes_eroded))

    def erode(self, location, seq, region):
        """ Erode (delete) some number of letters from the <location> side of
        <seq>, where location is 'left' or 'right'
        """
        n_to_erode = self.get_n_to_erode()
        fragment_before = ''
        fragment_after = ''
        if location == 'left':
            fragment_before = seq[:n_to_erode + 3] + '...'
            seq = seq[n_to_erode:len(seq)]
            fragment_after = seq[:n_to_erode + 3] + '...'
        elif location == 'right':
            fragment_before = '...' + seq[len(seq) - n_to_erode - 3 :]
            seq = seq[0:len(seq)-n_to_erode]
            fragment_after = '...' + seq[len(seq) - n_to_erode - 3 :]
        else:
            print 'ERROR location must be \"left\" or \"right\"'
            sys.exit()
        print '    %3d from %5s' % (n_to_erode, location),
        print ' of %s region: %15s' % (region, fragment_before),
        print ' --> %-15s' % fragment_after

    def get_insertion(self):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        length = int(numpy.random.exponential(self.mean_insertion_length))
        insert_seq_str = ''
        for _ in range(0, length):
            insert_seq_str += int_to_nucleotide(random.randint(0, 3))
        print '  inserting %d: %s' % (length, insert_seq_str)
        return Seq(insert_seq_str, generic_alphabet)

    def join(self, left_seq, insert_sequence, right_seq):
        """ Join three sequences in the indicated order.
        Return the new sequence.
        """
        new_seq = left_seq + insert_sequence + right_seq
        print '  joining: %s + %s + %s' % (left_seq, insert_sequence, right_seq)
        print '  gives: %s' % new_seq
        return new_seq

    def mutate(self, seq):
        """ Apply point mutations to <seq>, then return it """
        n_mutes = int(numpy.random.poisson(self.mute_rate*len(seq)))
        original_seq = seq
        print '  adding %3d point mutations' % n_mutes
        for _ in range(n_mutes):
            position = random.randint(0, len(seq)-1)
            old_nuke = seq[position]
            new_nuke = old_nuke
            while new_nuke == old_nuke:
                new_nuke = int_to_nucleotide(random.randint(0,3))
            print '    muting %s to %s at position %3d' % (old_nuke, new_nuke, position)
            seq = seq[:position] + new_nuke + seq[position+1:]

        return seq
