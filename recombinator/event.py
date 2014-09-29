""" Container to hold the information for a single recombination event. """
import csv
import random
import numpy
import os

from utils.opener import opener
from utils import utils

#----------------------------------------------------------------------------------------
class RecombinationEvent(object):
    """ Container to hold the information for a single recombination event. """
    def __init__(self, germlines):
        self.germlines = germlines
        self.vdj_combo_label = ()  # A tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length, <erosion lengths>)
                                   # NOTE I leave the lengths in here as strings
        self.genes = {}
        self.original_seqs = {}
        self.eroded_seqs = {}
        self.cyst_position = -1
        self.tryp_position = -1  # NOTE this is the position *within* the j gene *only*
        self.final_tryp_position = -1  # while *this* is the tryp position in the final recombined sequence
        self.erosions = {}  # erosion lengths for the event
        self.cdr3_length = 0  # NOTE this is the *desired* cdr3_length, i.e. after erosion and insertion
        self.insertion_lengths = {}
        self.insertions = {}
        self.recombined_seq = ''  # combined sequence *before* mutations
        self.final_seqs = []
        self.original_cyst_word = ''
        self.original_tryp_word = ''

    # ----------------------------------------------------------------------------------------
    def set_vdj_combo(self, vdj_combo_label, cyst_positions, tryp_positions, all_seqs, debug=False):
        """ Set the label which labels the gene/length choice (a tuple of strings) as well as it's constituent parts """
        self.vdj_combo_label = vdj_combo_label
        for region in utils.regions:
            self.genes[region] = vdj_combo_label[utils.index_keys[region + '_gene']]
            self.original_seqs[region] = all_seqs[region][self.genes[region]]
            self.original_seqs[region] = self.original_seqs[region].replace('N', utils.int_to_nucleotide(random.randint(0, 3)))  # replace any Ns with a random nuke (a.t.m. use the same nuke for all Ns in a given seq)
        self.cyst_position = cyst_positions[self.genes['v']]['cysteine-position']
        self.tryp_position = int(tryp_positions[self.genes['j']])  # tryp position in j
        self.cdr3_length = int(vdj_combo_label[utils.index_keys['cdr3_length']])
        for boundary in utils.boundaries:
            self.insertion_lengths[boundary] = int(vdj_combo_label[utils.index_keys[boundary + '_insertion']])
        for erosion_location in utils.erosions:
            self.erosions[erosion_location] = int(vdj_combo_label[utils.index_keys[erosion_location + '_del']])

        # set the original conserved codon words, so we can revert them if they get mutated
        self.original_cyst_word = str(self.original_seqs['v'][self.cyst_position : self.cyst_position + 3 ])
        self.original_tryp_word = str(self.original_seqs['j'][self.tryp_position : self.tryp_position + 3 ])

        if debug:
            self.print_gene_choice()

    # ----------------------------------------------------------------------------------------
    def set_final_tryp_position(self, debug=False):
        """ Set tryp position in the final, combined sequence. """
        self.final_tryp_position = utils.find_tryp_in_joined_seq(self.tryp_position,
                                                                self.eroded_seqs['v'],
                                                                self.insertions['vd'],
                                                                self.eroded_seqs['d'],
                                                                self.insertions['dj'],
                                                                self.eroded_seqs['j'],
                                                                self.erosions['j_5p'])
        if debug:
            print '  final tryptophan position: %d' % self.final_tryp_position
        # make sure cdr3 length matches the desired length in vdj_combo_label
        final_cdr3_length = self.final_tryp_position - self.cyst_position + 3
        if debug:
            print '  final_tryp_position - cyst_position + 3 = %d - %d + 3 = %d (should be %d)' % (self.final_tryp_position, self.cyst_position, final_cdr3_length, self.cdr3_length)
        utils.check_conserved_codons(self.eroded_seqs['v'] + self.insertions['vd'] + self.eroded_seqs['d'] + self.insertions['dj'] + self.eroded_seqs['j'], self.cyst_position, self.final_tryp_position)
        assert final_cdr3_length == int(self.cdr3_length)

    # ----------------------------------------------------------------------------------------
    def write_event(self, outfile, total_length_from_right=0):
        """ Write out all info to csv file. """
        columns = ('unique_id', 'reco_id', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'vd_insertion', 'dj_insertion', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'seq')
        mode = ''
        if os.path.isfile(outfile):
            mode = 'ab'
        else:
            mode = 'wb'
        with opener(mode)(outfile) as csvfile:
            writer = csv.DictWriter(csvfile, columns)
            if mode == 'wb':  # write the header if file wasn't there before
                writer.writeheader()
            # fill the row with values
            row = {}
            # first the stuff that's common to the whole recombination event
            row['cdr3_length'] = self.cdr3_length
            for region in utils.regions:
                row[region + '_gene'] = self.genes[region]
            for boundary in utils.boundaries:
                row[boundary + '_insertion'] = self.insertions[boundary]
            for erosion_location in utils.erosions:
                row[erosion_location + '_del'] = self.erosions[erosion_location]
            # hash the information that uniquely identifies each recombination event
            reco_id = ''
            for column in row:
                assert 'unique_id' not in row
                assert 'seq' not in row
                reco_id += str(row[column])
            row['reco_id'] = hash(reco_id)
            # then the stuff that's particular to each mutant/clone
            for imute in range(len(self.final_seqs)):
                row['seq'] = self.final_seqs[imute]
                if total_length_from_right > 0:
                    row['seq'] = row['seq'][len(row['seq'])-total_length_from_right : ]
                unique_id = ''  # Hash to uniquely identify the sequence.
                for column in row:
                    unique_id += str(row[column])
                unique_id += str(numpy.random.uniform())
                row['unique_id'] = hash(unique_id)
                writer.writerow(row)

    # ----------------------------------------------------------------------------------------
    def print_event(self, total_length_from_right=0):
        line = {}  # collect some information into a form that print_reco_event understands
        line['cdr3_length'] = self.cdr3_length
        for region in utils.regions:
            line[region + '_gene'] = self.genes[region]
        for boundary in utils.boundaries:
            line[boundary + '_insertion'] = self.insertions[boundary]
        for erosion_location in utils.erosions:
            line[erosion_location + '_del'] = self.erosions[erosion_location]
        for imute in range(len(self.final_seqs)):
            line['seq'] = self.final_seqs[imute]
            if total_length_from_right > 0:
                line['v_5p_del'] = len(line['seq']) - total_length_from_right
                line['seq'] = line['seq'][len(line['seq'])-total_length_from_right : ]
            utils.print_reco_event(self.germlines, line, self.cyst_position, self.final_tryp_position, one_line=(imute!=0))

    # ----------------------------------------------------------------------------------------
    def print_gene_choice(self):
        print '    chose:  gene             length'
        for region in utils.regions:
            print '        %s  %-18s %-3d' % (region, self.genes[region], len(self.original_seqs[region])),
            if region == 'v':
                print ' (cysteine: %d)' % self.cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % self.tryp_position
            else:
                print ''

    # ----------------------------------------------------------------------------------------
    def revert_conserved_codons(self, seq):
        """ revert conserved cysteine and tryptophan to their original bases, eg if they were messed up by s.h.m. """
        # TODO how badly does this screw up the tree you can infer from the seqs?
        cpos = self.cyst_position
        if seq[cpos : cpos + 3] != self.original_cyst_word:
            seq = seq[:cpos] + self.original_cyst_word + seq[cpos+3:]
        tpos = self.final_tryp_position
        if seq[tpos : tpos + 3] != self.original_tryp_word:
            seq = seq[:tpos] + self.original_tryp_word + seq[tpos+3:]

        return seq
