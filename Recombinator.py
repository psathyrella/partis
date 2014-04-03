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
    """ Convert ('v','d','j') --> (0,1,2) """
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
def check_conserved_cysteine(seq, cyst_position):
    """ Ensure there's a cysteine at <cyst_position> in <seq>. """
    cyst_word = str(seq[cyst_position:cyst_position+3])
    if cyst_word != 'TGT' and cyst_word != 'TGC':
        print 'ERROR cysteine in V is messed up: %s' % cyst_word
        sys.exit()
def check_conserved_tryptophan(seq, tryp_position):
    """ Ensure there's a tryptophan at <tryp_position> in <seq>. """
    tryp_word = str(seq[tryp_position:tryp_position+3])
    if tryp_word != 'TGG':
        print 'ERROR tryptophan in J is messed up: %s' % tryp_word
        sys.exit()
def check_conserved_codons(seq, cyst_position, tryp_position):
    """ Double check that we conserved the cysteine and the tryptophan. """
    check_conserved_cysteine(seq, cyst_position)
    check_conserved_tryptophan(seq, tryp_position)

#----------------------------------------------------------------------------------------
def is_position_protected(protected_positions, prospective_position):
    """ Would a mutation at <prospective_position> screw up a protected codon? """
    for position in protected_positions:
        if (prospective_position == position or
            prospective_position == (position + 1) or
            prospective_position == (position + 2)):
            return True
    return False

#----------------------------------------------------------------------------------------
def would_erode_conserved_codon(lengths, seqs, cyst_position, tryp_position):
    """ Would any of the erosion <lengths> delete a conserved codon? """
    if len(lengths) == 0:  # i.e. if we haven't yet filled the lengths
        return True
    # check conserved cysteine
    if len(seqs['v']) - lengths['v_right'] <= cyst_position + 2:
        print '      about to erode cysteine (%d), try again' % lengths['v_right']
        return True  # i.e. it *would* screw it up
    # check conserved tryptophan
    if lengths['j_left'] - 1 >= tryp_position:
        print '      about to erode tryptophan (%d), try again' % lengths['j_left']
        return True

    return False  # *whew*, it won't erode either of 'em

#----------------------------------------------------------------------------------------
def is_erosion_longer_than_seq(lengths, seqs):
    """ Are any of the proposed erosion <lengths> longer than the seq to be eroded? """
    if lengths['v_right'] > len(seqs['v']):  # NOTE not actually possible since we already know we didn't erode the cysteine
        print '      v_right erosion too long (%d)' % lengths['v_right']
        return True
    if lengths['d_left'] + lengths['d_right'] > len(seqs['d']):
        print '      d erosions too long (%d)' % (lengths['d_left'] + lengths['d_right'])
        return True
    if lengths['j_left'] > len(seqs['j']):  # NOTE also not possible for the same reason
        print '      j_left erosion too long (%d)' % lengths['j_left']
        return True
    return False

#----------------------------------------------------------------------------------------
def find_tryp_in_joined_seq(tryp_position_in_j, v_seq, vd_insertion, d_seq, dj_insertion, j_seq, j_erosion):
    """ Find the <end> tryptophan in a joined sequence.

    Given local tryptophan position in the j region, figure
    out what position it's at in the final sequence.
    NOTE tryp_position_in_j is the position *before* the j was eroded,
    but this fcn assumes that the j *has* been eroded.
    """
    check_conserved_tryptophan(j_seq, tryp_position_in_j - j_erosion)  # make sure tryp is where it's supposed to be
    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion)
    return tryp_position_in_j - j_erosion + length_to_left_of_j

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """

    # parameters that control recombination, erosion, and whatnot
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
        with open('data/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with open('data/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
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

        # choose germline seqs
        for region in self.regions:
            chosen_seqs[region] = self.all_seqs[region][vdj_combo_label[region_to_int(region)]]
            print '    found %s: %s' % (region, vdj_combo_label[region_to_int(region)]),
            if region == 'v':
                print ' (cysteine: %d)' % cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % tryp_position
            else:
                print ''

        # erode, insert, and combine
        self.erode_and_insert(reco_info, vdj_combo_label, chosen_seqs, cyst_position, tryp_position)
        print '  joining'
        print '         v: %s' % chosen_seqs['v']
        print '    insert: %s' % reco_info['vd_insertion']
        print '         d: %s' % chosen_seqs['d']
        print '    insert: %s' % reco_info['dj_insertion']
        print '         j: %s' % chosen_seqs['j']
        recombined_seq = chosen_seqs['v'] + reco_info['vd_insertion'] + chosen_seqs['d'] + reco_info['dj_insertion'] + chosen_seqs['j']
        final_tryp_position = find_tryp_in_joined_seq(tryp_position,  # NOTE remember this is the tryp position, in j alone, *before* erosion
                                                      chosen_seqs['v'],
                                                      reco_info['vd_insertion'],
                                                      chosen_seqs['d'],
                                                      reco_info['dj_insertion'],
                                                      chosen_seqs['j'],
                                                      reco_info['j_left'])
        print '  final tryptophan position: %d' % final_tryp_position

        # add point mutations
        reco_info['mutations'] = ''  # string to keep track of all mutations for output to csv
        final_seq = self.mutate(reco_info, recombined_seq, (cyst_position, final_tryp_position))
        check_conserved_codons(final_seq, cyst_position, final_tryp_position)
        reco_info['seq'] = final_seq

        # make sure cdr3 length matches the desired length in vdj_combo_label
        final_cdr3_length = final_tryp_position - cyst_position + 3
        assert final_cdr3_length == vdj_combo_label[3]

        if outfile != '':
            self.write_csv(outfile, reco_info)

        return final_seq

    def read_vdj_versions(self, region, fname):
        """ Read the various germline variants from file. """
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
        print 'ERROR shouldn\'t have fallen through'
        sys.exit()

    def get_insert_delete_lengths(self, cyst_position, tryp_position, desired_cdr3_length, current_cdr3_length):
        """ Choose random insertion and deletion lengths consisent with desired cdr3 length. """
        lengths = {}  # return values
        net_length_change = desired_cdr3_length - current_cdr3_length

        # first choose a total insertion length (vd + dj)
        total_insertion_length = int(numpy.random.exponential(2*self.mean_insertion_length))
        while total_insertion_length <= net_length_change:  # have to insert at least this many
            total_insertion_length = int(numpy.random.exponential(2*self.mean_insertion_length))
        total_deletion_length = total_insertion_length - net_length_change
        print '    trying totals: insert %d delete %d' % (total_insertion_length, total_deletion_length)
        assert total_insertion_length >= 0 and total_deletion_length >= 0

        # divide total_insertion_length into vd_insertion and dj_insertion
        partition_point = numpy.random.uniform(0, total_insertion_length)
        lengths['vd_insertion'] = int(round(partition_point))
        lengths['dj_insertion'] = total_insertion_length - lengths['vd_insertion']
        print '      insertion lengths: %d %d' % (lengths['vd_insertion'], lengths['dj_insertion'])
        assert lengths['vd_insertion'] + lengths['dj_insertion'] == total_insertion_length  # check for rounding problems

        # divide total_deletion_length into the four erosions
        partition_point_2 = int(round(numpy.random.uniform(0,total_deletion_length)))  # 'middle' point of the interval
        partition_point_1 = int(round(numpy.random.uniform(0, partition_point_2)))  # leftmost
        partition_point_3 = int(round(numpy.random.uniform(partition_point_2, total_deletion_length)))  # rightmost
        lengths['v_right'] = partition_point_1
        lengths['d_left'] = partition_point_2 - partition_point_1
        lengths['d_right'] = partition_point_3 - partition_point_2
        lengths['j_left'] = total_deletion_length - partition_point_3
        print '      erosion lengths: %d %d %d %d' % (lengths['v_right'], lengths['d_left'], lengths['d_right'], lengths['j_left'])
        assert lengths['v_right'] + lengths['d_left'] + lengths['d_right'] + lengths['j_left'] == total_deletion_length

        return lengths

    def erode(self, region, location, seqs, lengths, protected_position=-1):
        """ Erode some number of letters from seq.

        Nucleotides are removed from the <location> ('left' or 'right') side of
        <seq>. The codon beginning at index <protected_position> is optionally
        protected from removal.
        """
        seq = seqs[region]
        n_to_erode = lengths[region + '_' + location]
        if protected_position > 0:  # this check is redundant at this point
            if location == 'right' and region == 'v':
                if len(seq) - n_to_erode <= protected_position + 2:
                    assert False
            elif location == 'left' and region == 'j':
                if n_to_erode - 1 >= protected_position:
                    assert False
            else:
                print 'ERROR unanticipated protection'
                sys.exit()

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

    def get_insertion(self, length):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        insert_seq_str = ''
        for _ in range(0, length):
            insert_seq_str += int_to_nucleotide(random.randint(0, 3))
        return insert_seq_str
        
    def erode_and_insert(self, reco_info, vdj_combo_label, seqs, cyst_position, tryp_position):
        """ Erode and insert based on the cdr3 length in vdj_combo_label. """
        desired_cdr3_length = vdj_combo_label[3]
        # find 'current' cdr3 length (i.e. with no insertions or erosions)
        tryp_position_in_joined_seq = find_tryp_in_joined_seq(tryp_position, seqs['v'],  '', seqs['d'], '', seqs['j'], 0)
        check_conserved_codons(seqs['v'] + seqs['d'] + seqs['j'], cyst_position, tryp_position_in_joined_seq)
        current_cdr3_length = tryp_position_in_joined_seq - cyst_position + 3
        print '  choosing lengths'
        print '    cdr3: %d --> %d (%+d)' % (current_cdr3_length, desired_cdr3_length, desired_cdr3_length - current_cdr3_length)

        # Choose the deletion and insertion lengths, rechoosing if the original
        #   choice would have eroded a conserved codon
        lengths = {}
        while (would_erode_conserved_codon(lengths, seqs, cyst_position, tryp_position) or
               is_erosion_longer_than_seq(lengths, seqs)):
            lengths = self.get_insert_delete_lengths(cyst_position, tryp_position, desired_cdr3_length, current_cdr3_length)

        # copy lengths over to reco_info
        for key,val in lengths.iteritems():
            reco_info[key] = lengths[key]

        print '  eroding'
        seqs['v'] = self.erode('v', 'right', seqs, lengths, cyst_position)
        seqs['d'] = self.erode('d', 'left', seqs, lengths)
        seqs['d'] = self.erode('d', 'right', seqs, lengths)
        seqs['j'] = self.erode('j', 'left', seqs, lengths, tryp_position)
        for boundary in ('vd', 'dj'):
            reco_info[boundary + '_insertion'] = self.get_insertion(lengths[boundary + '_insertion'])

    def mutate(self, reco_info, seq, protected_positions):
        """ Apply point mutations to <seq>, then return it. """
        n_mutes = int(numpy.random.poisson(self.mute_rate*len(seq)))
        original_seq = seq
        print '  making %d point mutations: ' % n_mutes,
        mute_locations = []
        for _ in range(n_mutes):
            position = random.randint(0, len(seq)-1)
            while is_position_protected(protected_positions, position):
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
                if is_position_protected(protected_positions, ich):
                    print 'ERROR mutated a protected position'
                mute_print_str += '|'
            elif is_position_protected(protected_positions, ich):
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
            columns = ('vdj_combo', 'vd_insertion', 'dj_insertion', 'v_right', 'd_left', 'd_right', 'j_left', 'mutations', 'seq')
            csv_row = []
            for column in columns:
                if column == 'vdj_combo':
                    csv_row.append(reco_info[column][0])
                    csv_row.append(reco_info[column][1])
                    csv_row.append(reco_info[column][2])
                    csv_row.append(reco_info[column][3])
                else:
                    csv_row.append(reco_info[column])
                    
            csv_writer.writerow(csv_row)
