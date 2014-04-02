#!/usr/bin/env python
""" Simulates the process of VDJ recombination """ 
import sys
import csv
import json
import random
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_alphabet

#----------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """

    # parameters that control recombination, erosion, and whatnot
    mean_nukes_eroded = 4  # Mean number of nucleotides to remove at each side of the boundary.
                           # NOTE actual mean is a bit less than this because I round to an integer
    mute_rate = 0.06  # average number of point mutations per base
    mean_insertion_length = 6  # mean length of the non-templated insertion

    all_seqs = {}  # all the Vs, all the Ds...
    regions = ['v', 'd', 'j']
    version_freq_table = {}  # list of the probabilities with which each VDJ combo appears in data

    def __init__(self):
        """ Initialize from files """
        for region in self.regions:
            self.read_vdj_versions(region, 'data/igh'+region+'.fasta')
        self.read_vdj_version_freqs('version-counter/filtered-vdj-probs.txt')
        with open('data/v-meta.json','r') as jfile:
            self.json_data = json.load(jfile)

    def combine(self, outfile=''):
        """ Run the combination. Returns the new sequence. """
        vdj_combo_label = self.choose_vdj_combo()  # a tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length)
        reco_info = {}  # collect information about the recombination process for output to csv file
        reco_info['vdj_combo'] = vdj_combo_label
        chosen_seqs = {}
        for region in self.regions:
            chosen_seqs[region] = self.get_seq(region, vdj_combo_label)
        cysteine_position = self.json_data[vdj_combo_label[0]]['cysteine-position']
        print '  cysteine position: %d' % cysteine_position
        print '  eroding'
        chosen_seqs['v'] = self.erode(reco_info, 'right', chosen_seqs['v'], 'v', cysteine_position)
        chosen_seqs['d'] = self.erode(reco_info, 'left', chosen_seqs['d'], 'd')
        chosen_seqs['d'] = self.erode(reco_info, 'right', chosen_seqs['d'], 'd')
        chosen_seqs['j'] = self.erode(reco_info, 'left', chosen_seqs['j'], 'j')
        reco_info['vd_insertion'] = self.get_insertion()
        vd_seq = self.join(chosen_seqs['v'], reco_info['vd_insertion'], chosen_seqs['d'])
        reco_info['dj_insertion'] = self.get_insertion()
        recombined_seq = self.join(vd_seq, reco_info['dj_insertion'], chosen_seqs['j'])
        reco_info['mutations'] = ''
        final_seq = self.mutate(reco_info, recombined_seq, cysteine_position)
        reco_info['seq'] = final_seq
        print '  before mute:',recombined_seq
        print '  after mute: ',final_seq
        cysteine_word = str(final_seq[cysteine_position:cysteine_position+3])
        if cysteine_word != 'TGT' and cysteine_word != 'TGC':
            print 'ERROR cysteine in V is messed up, you probably eroded it off: %s' % cysteine_word
            sys.exit()
        if outfile != '':
            self.write_csv(outfile, reco_info)
        return final_seq

    def read_vdj_versions(self, region, fname):
        """ Read the various germline variants from file and store as
        Seq objects
        """
        self.all_seqs[region] = {}
        for seq_record in SeqIO.parse(fname, "fasta"):
            # This line works equally well without the string conversion
            # ... er, but it seems simpler to not use the fancy class if'n
            # I don't need it? It's the same speed both ways.
            self.all_seqs[region][seq_record.name] = str(seq_record.seq)

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
            sum_prob += self.version_freq_table[vdj_choice]
            if iprob < sum_prob:
                print ' ---> ',vdj_choice
                return vdj_choice
        print 'ERROR I shouldn\'t be here'
        sys.exit()

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

    def erode(self, reco_info, location, seq, region, cysteine_position=-1):
        """ Erode (delete) some number of letters from the <location> side of
        <seq>, where location is 'left' or 'right'
        """
        n_to_erode = self.get_n_to_erode()
        if cysteine_position > 0:
            if location != 'right' or region != 'v':
                print 'ERROR that doesn\'t make sense'
                sys.exit()
            while len(seq) - n_to_erode - 3 <= cysteine_position:
                print '    about to erode off cysteine (%d), try again' % n_to_erode
                n_to_erode = self.get_n_to_erode()

        reco_info[region + '_' + location] = n_to_erode
        fragment_before = ''
        fragment_after = ''
        if location == 'left':
            fragment_before = seq[:n_to_erode + 3] + '...'
            new_seq = seq[n_to_erode:len(seq)]
            fragment_after = new_seq[:n_to_erode + 3] + '...'
        elif location == 'right':
            fragment_before = '...' + seq[len(seq) - n_to_erode - 3 :]
            new_seq = seq[0:len(seq)-n_to_erode]
            fragment_after = '...' + new_seq[len(new_seq) - n_to_erode - 3 :]
        else:
            print 'ERROR location must be \"left\" or \"right\"'
            sys.exit()
        print '    %3d from %5s' % (n_to_erode, location),
        print ' of %s region: %15s' % (region, fragment_before),
        print ' --> %-15s' % fragment_after
        if len(fragment_after) == 0:
            print '    NOTE eroded away entire sequence'
        return new_seq

    def get_insertion(self):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        length = int(numpy.random.exponential(self.mean_insertion_length))
        insert_seq_str = ''
        for _ in range(0, length):
            insert_seq_str += int_to_nucleotide(random.randint(0, 3))
        print '  inserting %d: %s' % (length, insert_seq_str)
        return insert_seq_str

    def join(self, left_seq, insert_sequence, right_seq):
        """ Join three sequences in the indicated order.
        Return the new sequence.
        """
        new_seq = left_seq + insert_sequence + right_seq
        print '  joining: %s + %s + %s' % (left_seq, insert_sequence, right_seq)
        print '  gives: %s' % new_seq
        return new_seq

    def mutate(self, reco_info, seq, cysteine_position):
        """ Apply point mutations to <seq>, then return it """
        n_mutes = int(numpy.random.poisson(self.mute_rate*len(seq)))
        original_seq = seq
        print '  adding %d point mutations' % n_mutes
        for _ in range(n_mutes):
            position = random.randint(0, len(seq)-1)
            while (position == cysteine_position or
                   position == (cysteine_position + 1) or
                   position == (cysteine_position + 2)):
                print '  position is in cysteine codon (%d), trying again.' % position
                position = random.randint(0, len(seq)-1)
            old_nuke = seq[position]
            new_nuke = old_nuke
            while new_nuke == old_nuke:
                new_nuke = int_to_nucleotide(random.randint(0,3))
            print '    muting %s to %s at position %3d' % (old_nuke, new_nuke, position)
            reco_info['mutations'] += '%d%s%s-' % (position, old_nuke, new_nuke)
            seq = seq[:position] + new_nuke + seq[position+1:]

        reco_info['mutations'] = reco_info['mutations'].rstrip('-')
        return seq

    def write_csv(self, outfile, reco_info):
        with open(outfile, 'ab') as csvfile:
            csv_writer = csv.writer(csvfile)
            # NOTE this is the *original* cdr3_length, i.e. from Connor's file, *not* the length in the sequence I'm spitting out
            # TODO fix that shit
            columns = ('vdj_combo', 'v_right', 'd_left', 'd_right', 'j_left', 'vd_insertion', 'dj_insertion', 'mutations', 'seq')
            csv_row = []
            for column in columns:
                print column,
                if column == 'vdj_combo':
                    csv_row.append(reco_info[column][0])
                    csv_row.append(reco_info[column][1])
                    csv_row.append(reco_info[column][2])
                else:
                    csv_row.append(reco_info[column])
                    
            csv_writer.writerow(csv_row)
#            print csv_row
                    
