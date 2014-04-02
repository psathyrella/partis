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
def region_to_int(region):
    """ Convert region: ('v','d','j') --> 0,1,2 """
    if region == 'v':
        return 0
    elif region == 'd':
        return 1
    elif region == 'j':
        return 2
    else:
        print 'ERROR bad region'
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
        with open('data/v-meta.json') as json_file:  # get location of <begin> cysteine in v regions
            self.cyst_positions = json.load(json_file)
        with open('data/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in j regions (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

    def combine(self, outfile=''):
        """ Run the combination. Returns the new sequence. """
        vdj_combo_label = self.choose_vdj_combo()  # a tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length)
        reco_info = {}  # collect information about the recombination process for output to csv file
        reco_info['vdj_combo'] = vdj_combo_label
        chosen_seqs = {}
        cyst_position = self.cyst_positions[vdj_combo_label[0]]['cysteine-position']
        tryp_position = int(self.tryp_positions[vdj_combo_label[2]])
        for region in self.regions:
            chosen_seqs[region] = self.get_seq(region, vdj_combo_label)
            print '    found %s: %s' % (region, vdj_combo_label[region_to_int(region)]),
            if region == 'v':
                print ' (cysteine: %d)' % cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % tryp_position
            else:
                print ''
        print '  eroding'
        chosen_seqs['v'] = self.erode(reco_info, 'right', chosen_seqs['v'], 'v', cyst_position)
        chosen_seqs['d'] = self.erode(reco_info, 'left', chosen_seqs['d'], 'd')
        chosen_seqs['d'] = self.erode(reco_info, 'right', chosen_seqs['d'], 'd')
        chosen_seqs['j'] = self.erode(reco_info, 'left', chosen_seqs['j'], 'j', tryp_position)
        reco_info['vd_insertion'] = self.get_insertion()
        reco_info['dj_insertion'] = self.get_insertion()
        print '  joining'
        print '         v: %s' % chosen_seqs['v']
        print '    insert: %s' % reco_info['vd_insertion']
        print '         d: %s' % chosen_seqs['d']
        print '    insert: %s' % reco_info['dj_insertion']
        print '         j: %s' % chosen_seqs['j']
        recombined_seq = chosen_seqs['v'] + reco_info['vd_insertion'] + chosen_seqs['d'] + reco_info['dj_insertion'] + chosen_seqs['j']
        reco_info['mutations'] = ''
        # figure out the position the tryptophan's at in the final sequence
        length_to_left_of_j = len(chosen_seqs['v'] + reco_info['vd_insertion'] + chosen_seqs['d'] + reco_info['dj_insertion'])
        final_tryp_position = tryp_position - reco_info['j_left'] + length_to_left_of_j
        print '  final tryptophan position: %d' % final_tryp_position
        final_seq = self.mutate(reco_info, recombined_seq, (cyst_position, final_tryp_position))
        reco_info['seq'] = final_seq
        self.check_conserved_codons(final_seq, cyst_position, final_tryp_position)
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
        for vdj_choice in self.version_freq_table:
            sum_prob += self.version_freq_table[vdj_choice]
            if iprob < sum_prob:
                print '  chose combination: ',vdj_choice
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
        return self.all_seqs[region][gene_name]

    def get_n_to_erode(self):
        """ Number of bases to remove before joining to the neighboring
        sequence
        """
        # NOTE casting to an int reduces the mean somewhat
        return int(numpy.random.exponential(self.mean_nukes_eroded))

    def erode(self, reco_info, location, seq, region, protected_position=-1):
        """ Erode some number of letters from seq.

        Nucleotides are removed from the <location> ('left' or 'right') side of
        <seq>. The codon beginning at index <protected_position> is protected
        from removal.
        """
        n_to_erode = self.get_n_to_erode()
        if protected_position > 0:
            if location == 'right' and region == 'v':
                while len(seq) - n_to_erode <= protected_position + 2:
                    print '    about to erode cysteine (%d), try again' % n_to_erode
                    n_to_erode = self.get_n_to_erode()
            elif location == 'left' and region == 'j':
                while n_to_erode - 1 >= protected_position:
                    print '    about to erode tryptophan (%d), try again' % n_to_erode
                    n_to_erode = self.get_n_to_erode()
            else:
                print 'ERROR unanticipated protection'
                sys.exit()

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
        return insert_seq_str

    def is_position_protected(self, protected_positions, prospective_position):
        for position in protected_positions:
            if (prospective_position == position or
                prospective_position == (position + 1) or
                prospective_position == (position + 2)):
                return True
        return False

    def mutate(self, reco_info, seq, protected_positions):
        """ Apply point mutations to <seq>, then return it """
        n_mutes = int(numpy.random.poisson(self.mute_rate*len(seq)))
        original_seq = seq
        print '  making %d point mutations: ' % n_mutes,
        mute_locations = []
        for _ in range(n_mutes):
            position = random.randint(0, len(seq)-1)
            while self.is_position_protected(protected_positions, position):
                position = random.randint(0, len(seq)-1)
            mute_locations.append(position)
            old_nuke = seq[position]
            new_nuke = old_nuke
            while new_nuke == old_nuke:
                new_nuke = int_to_nucleotide(random.randint(0,3))
            print '%d%s%s' % (position, old_nuke, new_nuke),
            reco_info['mutations'] += '%d%s%s-' % (position, old_nuke, new_nuke)
            seq = seq[:position] + new_nuke + seq[position+1:]
        print ''
        reco_info['mutations'] = reco_info['mutations'].rstrip('-')

        # print out mutation locations
        mute_locations = set(mute_locations)
        mute_print_str = ''
        for ich in range(len(seq)):
            if ich in mute_locations:
                if self.is_position_protected(protected_positions, ich):
                    print 'ERROR mutated a protected position'
                mute_print_str += '|'
            elif self.is_position_protected(protected_positions, ich):
                mute_print_str += 'o'
            else:
                mute_print_str += ' '
        print '  before mute:',original_seq
        print '  mutations:  ',mute_print_str,'    (o: protected codons)'
        print '  after mute: ',seq

        return seq

    def write_csv(self, outfile, reco_info):
        with open(outfile, 'ab') as csvfile:
            csv_writer = csv.writer(csvfile)
            # NOTE this is the *original* cdr3_length, i.e. from Connor's file, *not* the length in the sequence I'm spitting out
            # TODO fix that shit
            columns = ('vdj_combo', 'v_right', 'd_left', 'd_right', 'j_left', 'vd_insertion', 'dj_insertion', 'mutations', 'seq')
            csv_row = []
            for column in columns:
                if column == 'vdj_combo':
                    csv_row.append(reco_info[column][0])
                    csv_row.append(reco_info[column][1])
                    csv_row.append(reco_info[column][2])
                else:
                    csv_row.append(reco_info[column])
                    
            csv_writer.writerow(csv_row)
                    
    def check_conserved_codons(self, seq, cyst_position, tryp_position):
        """ Double check that we conserved the cysteine and the tryptophan. """
        cyst_word = str(seq[cyst_position:cyst_position+3])
        if cyst_word != 'TGT' and cyst_word != 'TGC':
            print 'ERROR cysteine in V is messed up: %s' % cyst_word
            sys.exit()
        tryp_word = str(seq[tryp_position:tryp_position+3])
        if tryp_word != 'TGG':
            print 'ERROR tryptophan in J is messed up: %s' % tryp_word
            sys.exit()

