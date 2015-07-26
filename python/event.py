""" Container to hold the information for a single recombination event. """
import csv
import random
import numpy
import os

from opener import opener
import utils

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
        self.local_cyst_position, self.final_cyst_position = -1, -1  # the 'local' one is without the v left erosion
        self.local_tryp_position, self.final_tryp_position = -1, -1  # NOTE the first is the position *within* the j gene *only*, while the second is the tryp position in the final recombined sequence
        self.erosions = {}  # erosion lengths for the event
        self.effective_erosions = {}  # v left and j right erosions
        self.cdr3_length = 0  # NOTE this is the *desired* cdr3_length, i.e. after erosion and insertion
        self.insertion_lengths = {}
        self.insertions = {}
        self.recombined_seq = ''  # combined sequence *before* mutations
        self.final_seqs, self.indelfo = [], []
        self.original_cyst_word = ''
        self.original_tryp_word = ''

    # ----------------------------------------------------------------------------------------
    def set_vdj_combo(self, vdj_combo_label, cyst_positions, tryp_positions, all_seqs, debug=False, dont_mimic_data_read_length=False):
        """ Set the label which labels the gene/length choice (a tuple of strings) as well as it's constituent parts """
        self.vdj_combo_label = vdj_combo_label
        for region in utils.regions:
            self.genes[region] = vdj_combo_label[utils.index_keys[region + '_gene']]
            self.original_seqs[region] = all_seqs[region][self.genes[region]]
            self.original_seqs[region] = self.original_seqs[region].replace('N', utils.int_to_nucleotide(random.randint(0, 3)))  # replace any Ns with a random nuke (a.t.m. use the same nuke for all Ns in a given seq)
        self.local_cyst_position = cyst_positions[self.genes['v']]['cysteine-position']  # cyst position in uneroded v
        self.local_tryp_position = int(tryp_positions[self.genes['j']])  # tryp position within j only
        self.cdr3_length = int(vdj_combo_label[utils.index_keys['cdr3_length']])
        for boundary in utils.boundaries:
            self.insertion_lengths[boundary] = int(vdj_combo_label[utils.index_keys[boundary + '_insertion']])
        for erosion in utils.real_erosions:
            self.erosions[erosion] = int(vdj_combo_label[utils.index_keys[erosion + '_del']])
        for erosion in utils.effective_erosions:
            if not dont_mimic_data_read_length:  # use v left and j right erosions from data?
                self.effective_erosions[erosion] = int(vdj_combo_label[utils.index_keys[erosion + '_del']])
            else:  # otherwise ignore data, and keep the entire v and j genes
                self.effective_erosions[erosion] = 0

        # set the original conserved codon words, so we can revert them if they get mutated
        self.original_cyst_word = str(self.original_seqs['v'][self.local_cyst_position : self.local_cyst_position + 3 ])
        self.original_tryp_word = str(self.original_seqs['j'][self.local_tryp_position : self.local_tryp_position + 3 ])

        if debug:
            self.print_gene_choice()

    # ----------------------------------------------------------------------------------------
    def set_final_cyst_tryp_positions(self, debug=False, total_length_from_right=-1):
        """ Set tryp position in the final, combined sequence. """
        self.final_cyst_position = self.local_cyst_position - self.effective_erosions['v_5p']
        self.final_tryp_position = utils.find_tryp_in_joined_seq(self.local_tryp_position,
                                                                self.eroded_seqs['v'],
                                                                self.insertions['vd'],
                                                                self.eroded_seqs['d'],
                                                                self.insertions['dj'],
                                                                self.eroded_seqs['j'],
                                                                self.erosions['j_5p'])
        if debug:
            print '  final tryptophan position: %d' % self.final_tryp_position
        # make sure cdr3 length matches the desired length in vdj_combo_label
        final_cdr3_length = self.final_tryp_position - self.final_cyst_position + 3
        if debug:
            print '  final_tryp_position - final_cyst_position + 3 = %d - %d + 3 = %d (should be %d)' % (self.final_tryp_position, self.final_cyst_position, final_cdr3_length, self.cdr3_length)
        utils.check_both_conserved_codons(self.eroded_seqs['v'] + self.insertions['vd'] + self.eroded_seqs['d'] + self.insertions['dj'] + self.eroded_seqs['j'], self.final_cyst_position, self.final_tryp_position)
            
        assert final_cdr3_length == int(self.cdr3_length)

        assert total_length_from_right == -1  # deprecated (I think) now that I'm adding the mimic_data_read_length option
        # if total_length_from_right >= 0:
        #     self.effective_erosions['v_5p'] = len(self.recombined_seq) - total_length_from_right

    # ----------------------------------------------------------------------------------------
    def write_event(self, outfile, total_length_from_right=0, irandom=None):
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
        columns = ('unique_id', 'reco_id') + utils.index_columns + ('seq', 'indels')
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
            for erosion in utils.real_erosions:
                row[erosion + '_del'] = self.erosions[erosion]
            for erosion in utils.effective_erosions:
                row[erosion + '_del'] = self.effective_erosions[erosion]
            # hash the information that uniquely identifies each recombination event
            reco_id = ''
            for column in row:
                assert 'unique_id' not in row
                assert 'seq' not in row
                reco_id += str(row[column])
            row['reco_id'] = hash(reco_id)
            assert 'fv_insertion' not in row  # well, in principle it's ok if they're there, but in that case I'll need to at least think about updating some things
            assert 'jf_insertion' not in row
            row['fv_insertion'] = ''
            row['jf_insertion'] = ''
            # then the stuff that's particular to each mutant/clone
            for imute in range(len(self.final_seqs)):
                row['seq'] = self.final_seqs[imute]
                if total_length_from_right > 0:
                    row['seq'] = row['seq'][len(row['seq'])-total_length_from_right : ]
                unique_id = ''  # Hash to uniquely identify the sequence.
                for column in row:
                    unique_id += str(row[column])
                if irandom is None:  # NOTE see note above
                    unique_id += str(numpy.random.uniform())
                else:
                    unique_id += str(irandom)
                row['unique_id'] = hash(unique_id)
                row['indels'] = self.indelfo[imute]
                writer.writerow(row)

    # ----------------------------------------------------------------------------------------
    def print_event(self, total_length_from_right=0):
        line = {}  # collect some information into a form that print_reco_event understands
        line['cdr3_length'] = self.cdr3_length
        for region in utils.regions:
            line[region + '_gene'] = self.genes[region]
        for boundary in utils.boundaries:
            line[boundary + '_insertion'] = self.insertions[boundary]
        for erosion in utils.real_erosions:
            line[erosion + '_del'] = self.erosions[erosion]
        for erosion in utils.effective_erosions:
            line[erosion + '_del'] = self.effective_erosions[erosion]
        line['cyst_position'] = self.final_cyst_position
        line['tryp_position'] = self.final_tryp_position
        assert 'fv_insertion' not in line  # well, in principle it's ok if they're there, but in that case I'll need to at least think about updating some things
        assert 'jf_insertion' not in line
        line['fv_insertion'] = ''
        line['jf_insertion'] = ''
        line['seqs'] = self.final_seqs
        line['unique_ids'] = [i for i in range(len(self.final_seqs))]
        utils.print_reco_event(self.germlines, line, indelfos=self.indelfo)
        # for imute in range(len(self.final_seqs)):
        #     line['seq'] = self.final_seqs[imute]
        #     if total_length_from_right > 0:
        #         line['seq'] = line['seq'][len(line['seq'])-total_length_from_right : ]
        #     line['indels'] = self.indelfo[imute]
        #     utils.print_reco_event(self.germlines, line, one_line=(imute!=0 and len(self.indelfo[imute]['indels']) == 0))

    # ----------------------------------------------------------------------------------------
    def print_gene_choice(self):
        print '    chose:  gene             length'
        for region in utils.regions:
            print '        %s  %-18s %-3d' % (region, self.genes[region], len(self.original_seqs[region])),
            if region == 'v':
                print ' (cysteine: %d)' % self.local_cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % self.local_tryp_position
            else:
                print ''

    # ----------------------------------------------------------------------------------------
    def revert_conserved_codons(self, seq):
        """ revert conserved cysteine and tryptophan to their original bases, eg if they were messed up by s.h.m. """
        cpos = self.final_cyst_position
        if seq[cpos : cpos + 3] != self.original_cyst_word:
            seq = seq[:cpos] + self.original_cyst_word + seq[cpos+3:]
        tpos = self.final_tryp_position
        if seq[tpos : tpos + 3] != self.original_tryp_word:
            seq = seq[:tpos] + self.original_tryp_word + seq[tpos+3:]

        return seq
