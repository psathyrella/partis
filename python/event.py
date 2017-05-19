""" Container to hold the information for a single recombination event. """
import csv
import sys
import random
import numpy
import os
import copy

import utils

#----------------------------------------------------------------------------------------
class RecombinationEvent(object):
    """ Container to hold the information for a single recombination event. """
    def __init__(self, glfo):
        self.glfo = glfo
        self.vdj_combo_label = ()  # A tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length, <erosion lengths>)
                                   # NOTE I leave the lengths in here as strings
        self.genes = {}
        self.original_seqs = {}
        self.eroded_seqs = {}
        self.local_codon_positions = {r : -1 for r in utils.conserved_codons[glfo['locus']]}  # local: without the v left erosion
        self.final_codon_positions = {r : -1 for r in utils.conserved_codons[glfo['locus']]}  # final: with it
        self.erosions = {}  # erosion lengths for the event
        self.effective_erosions = {}  # v left and j right erosions
        self.cdr3_length = 0  # NOTE this is the *desired* cdr3_length, i.e. after erosion and insertion
        self.insertion_lengths = {}
        self.insertions = {}
        self.recombined_seq = ''  # combined sequence *before* mutations
        self.final_seqs, self.indelfos = [], []
        self.unmutated_codons = None

        self.line = None  # dict with info in format of utils.py/output files

    # ----------------------------------------------------------------------------------------
    def set_vdj_combo(self, vdj_combo_label, glfo, debug=False, mimic_data_read_length=False):
        """ Set the label which labels the gene/length choice (a tuple of strings) as well as it's constituent parts """
        self.vdj_combo_label = vdj_combo_label
        for region in utils.regions:
            self.genes[region] = vdj_combo_label[utils.index_keys[region + '_gene']]
            self.original_seqs[region] = glfo['seqs'][region][self.genes[region]]
            self.original_seqs[region] = self.original_seqs[region].replace('N', utils.int_to_nucleotide(random.randint(0, 3)))  # replace any Ns with a random nuke (a.t.m. use the same nuke for all Ns in a given seq)
        for region, codon in utils.conserved_codons[glfo['locus']].items():
            self.local_codon_positions[region] = glfo[codon + '-positions'][self.genes[region]]  # position in uneroded germline gene
        for boundary in utils.boundaries:
            self.insertion_lengths[boundary] = int(vdj_combo_label[utils.index_keys[boundary + '_insertion']])
        for erosion in utils.real_erosions:
            self.erosions[erosion] = int(vdj_combo_label[utils.index_keys[erosion + '_del']])
        for erosion in utils.effective_erosions:
            if mimic_data_read_length:  # use v left and j right erosions from data?
                self.effective_erosions[erosion] = int(vdj_combo_label[utils.index_keys[erosion + '_del']])
            else:  # otherwise ignore data, and keep the entire v and j genes
                self.effective_erosions[erosion] = 0

        if debug:
            self.print_gene_choice()

    # ----------------------------------------------------------------------------------------
    def set_final_codon_positions(self):
        """ Set tryp position in the final, combined sequence. """
        self.final_codon_positions['v'] = self.local_codon_positions['v'] - self.effective_erosions['v_5p']
        length_to_left_of_j = len(self.eroded_seqs['v'] + self.insertions['vd'] + self.eroded_seqs['d'] + self.insertions['dj'])
        self.final_codon_positions['j'] = self.local_codon_positions['j'] - self.erosions['j_5p'] + length_to_left_of_j
        self.cdr3_length = self.final_codon_positions['j'] - self.final_codon_positions['v'] + 3

    # ----------------------------------------------------------------------------------------
    def write_event(self, outfile, irandom=None):
        """ 
        Write out all info to csv file.
        NOTE/RANT so, in calculating each sequence's unique id, we need to hash more than the information about the rearrangement
            event and mutation, because if we create identical events and sequences in independent recombinator threads, we *need* them
            to have different unique ids (otherwise all hell will break loose when you try to analyze them). The easy way to avoid this is
            to add a random number to the information before you hash it... but then you have no way to reproduce that random number when 
            you want to run again with a set random seed to get identical output. The FIX for this at the moment is to pass in <irandom>, i.e.
            the calling proc tells write_event() that we're writing the <irandom>th event that that calling event is working on. Which effectively
            means we (drastically) reduce the period of our random number generator for hashing in exchange for reproducibility. Should be ok...
        """
        columns = ('unique_ids', 'reco_id') + utils.index_columns + ('cdr3_length', 'input_seqs', 'indel_reversed_seqs', 'indelfos')
        mode = ''
        if os.path.isfile(outfile):
            mode = 'ab'
        else:
            mode = 'wb'
        with open(outfile, mode) as csvfile:
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
            for erosion in utils.real_erosions:
                row[erosion + '_del'] = self.erosions[erosion]
            for erosion in utils.effective_erosions:
                row[erosion + '_del'] = self.effective_erosions[erosion]
            # hash the information that uniquely identifies each recombination event
            str_for_reco_id = ''
            for column in row:
                assert 'unique_ids' not in row
                assert 'seqs' not in row
                str_for_reco_id += str(row[column])
            row['reco_id'] = hash(str_for_reco_id)  # note that this gives the same reco id for the same rearrangement parameters, even if they come from a separate rearrangement event
            assert 'fv_insertion' not in row  # well, in principle it's ok if they're there, but in that case I'll need to at least think about updating some things
            assert 'jf_insertion' not in row
            row['fv_insertion'] = ''
            row['jf_insertion'] = ''
            # then the stuff that's particular to each mutant/clone
            for imute in range(len(self.final_seqs)):
                row['seqs'] = [self.indelfos[imute]['reversed_seq'], ]  # add this as 'seqs' (instead of 'indel_reversed_seqs'), since we want get_line_for_output() to use empty strings if there's no indels
                row['input_seqs'] = [self.final_seqs[imute], ]
                str_for_unique_id = ''  # Hash to uniquely identify the sequence.
                for column in row:
                    str_for_unique_id += str(row[column])
                if irandom is None:  # NOTE see note above
                    str_for_unique_id += str(numpy.random.uniform())
                else:
                    str_for_unique_id += str(irandom)
                row['unique_ids'] = [hash(str_for_unique_id), ]
                row['indelfos'] = [self.indelfos[imute], ]
                writer.writerow(utils.get_line_for_output(row))

    # ----------------------------------------------------------------------------------------
    def getline(self):  # don't access <self.line> directly
        if self.line is not None:
            return self.line

        line = {}  # collect some information into a form that the print fcn understands
        for region in utils.regions:
            line[region + '_gene'] = self.genes[region]
        for boundary in utils.boundaries:
            line[boundary + '_insertion'] = self.insertions[boundary]
        for erosion in utils.real_erosions:
            line[erosion + '_del'] = self.erosions[erosion]
        for erosion in utils.effective_erosions:
            line[erosion + '_del'] = self.effective_erosions[erosion]
        assert 'fv_insertion' not in line  # well, in principle it's ok if they're there, but in that case I'll need to at least think about updating some things
        assert 'jf_insertion' not in line
        line['fv_insertion'] = ''
        line['jf_insertion'] = ''
        line['input_seqs'] = self.final_seqs
        line['indel_reversed_seqs'] = []
        for iseq in range(len(self.indelfos)):
            if self.indelfos[iseq]['reversed_seq'] != '':
                line['indel_reversed_seqs'].append(self.indelfos[iseq]['reversed_seq'])
            else:
                line['indel_reversed_seqs'].append(line['input_seqs'][iseq])
        line['seqs'] = line['indel_reversed_seqs']
        line['indelfos'] = self.indelfos
        line['unique_ids'] = [str(i) for i in range(len(self.final_seqs))]
        line['cdr3_length'] = self.cdr3_length
        line['codon_positions'] = copy.deepcopy(self.final_codon_positions)
        utils.add_implicit_info(self.glfo, line)

        self.line = line
        return self.line

    # ----------------------------------------------------------------------------------------
    def print_event(self):
        utils.print_reco_event(self.getline(), extra_str='    ')  # don't access <self.line> directly

    # ----------------------------------------------------------------------------------------
    def print_gene_choice(self):
        print '    chose:  gene             length'
        for region in utils.regions:
            print '        %s  %-18s %-3d' % (region, self.genes[region], len(self.original_seqs[region])),
            if region in self.local_codon_positions:
                print ' (%s: %d)' % (utils.conserved_codons[self.glfo['locus']][region], self.local_codon_positions[region])
            else:
                print ''

    # ----------------------------------------------------------------------------------------
    def revert_conserved_codons(self, seq):
        """ revert conserved cysteine and tryptophan to their original bases, eg if they were messed up by s.h.m. """
        for region, pos in self.final_codon_positions.items():
            if seq[pos : pos + 3] != self.unmutated_codons[region]:
                assert len(self.unmutated_codons[region]) == 3
                seq = seq[:pos] + self.unmutated_codons[region] + seq[pos + 3 :]
            assert utils.codon_unmutated(utils.conserved_codons[self.glfo['locus']][region], seq, pos)
        return seq
