""" Container to hold the information for a single recombination event. """
import csv
import sys
import random
import numpy
import os
import copy

import utils
import indelutils
import treeutils

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

        # per-rearrangement even codon positions (i.e. per-sequence, i.e. *pre*-shm indels)
        self.pre_erosion_codon_positions = {r : -1 for r in utils.conserved_codons[glfo['locus']]}  # local: without the v left erosion
        self.post_erosion_codon_positions = {r : -1 for r in utils.conserved_codons[glfo['locus']]}  # final: with it

        self.erosions = {}  # erosion lengths for the event
        self.effective_erosions = {}  # v left and j right erosions
        self.cdr3_length = 0  # NOTE this is the *desired* cdr3_length, i.e. after erosion and insertion
        self.insertion_lengths = {}
        self.insertions = {}
        self.recombined_seq = ''  # combined sequence *before* mutations
        self.final_seqs, self.indelfos, self.final_codon_positions = [], [], []
        self.unmutated_codons = None
        self.leaf_names = None  # keeps track of leaf names (of form t<n>) so that when we get to deciding on final uids, we still know which leaf in <self.tree> corresponds to which sequence (order in <self.leaf_names> is same as in <self.final_seqs>)
        self.tree = None

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
            self.pre_erosion_codon_positions[region] = glfo[codon + '-positions'][self.genes[region]]  # position in uneroded germline gene
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
    def set_naive_seq(self, use_dummy_insertions=False):
        if use_dummy_insertions:
            vd_str, dj_str = utils.ambiguous_bases[0] * self.insertion_lengths['vd'], utils.ambiguous_bases[0] * self.insertion_lengths['dj']
        else:
            vd_str, dj_str = self.insertions['vd'], self.insertions['dj']
        self.recombined_seq = self.eroded_seqs['v'] + vd_str + self.eroded_seqs['d'] + dj_str + self.eroded_seqs['j']

    # ----------------------------------------------------------------------------------------
    def is_there_a_stop_codon(self):
        return utils.is_there_a_stop_codon(self.recombined_seq, '', '', self.effective_erosions['v_5p'])  # fv and jf insertions are always empty (not even defined in reco_event, at least at this point), I don't remember why, but I think it's on purpose

    # ----------------------------------------------------------------------------------------
    def set_post_erosion_codon_positions(self):
        """ Set tryp position in the final, combined sequence. """
        self.post_erosion_codon_positions['v'] = self.pre_erosion_codon_positions['v'] - self.effective_erosions['v_5p']
        length_to_left_of_j = len(self.eroded_seqs['v'] + self.insertions['vd'] + self.eroded_seqs['d'] + self.insertions['dj'])
        self.post_erosion_codon_positions['j'] = self.pre_erosion_codon_positions['j'] - self.erosions['j_5p'] + length_to_left_of_j
        self.cdr3_length = self.post_erosion_codon_positions['j'] - self.post_erosion_codon_positions['v'] + 3

    # ----------------------------------------------------------------------------------------
    def set_reco_id(self, line, irandom=None):
        # NOTE maybe this rant is deprecated?
        """ 
        NOTE/RANT so, in calculating each sequence's unique id, we need to hash more than the information about the rearrangement
            event and mutation, because if we create identical events and sequences in independent recombinator threads, we *need* them
            to have different unique ids (otherwise all hell will break loose when you try to analyze them). The easy way to avoid this is
            to add a random number to the information before you hash it... but then you have no way to reproduce that random number when 
            you want to run again with a set random seed to get identical output. The FIX for this at the moment is to pass in <irandom>, i.e.
            the calling proc tells write_event() that we're writing the <irandom>th event that that calling event is working on. Which effectively
            means we (drastically) reduce the period of our random number generator for hashing in exchange for reproducibility. Should be ok...
        """
        reco_id_columns = [r + '_gene' for r in utils.regions] + [b + '_insertion' for b in utils.boundaries] + [e + '_del' for e in utils.all_erosions]
        def randstr():  # gah, duplicated in set_unique_ids()
            return str(numpy.random.uniform() if irandom is None else irandom)
        reco_id_str = ''.join([str(line[c]) for c in reco_id_columns])
        line['reco_id'] = hash(reco_id_str)  # note that this gives the same reco id for the same rearrangement parameters, even if they come from a separate rearrangement event
        return reco_id_str  # this is pretty hackey, but I want to split up the reco and unique id setting so I can call only the former from bin/bcr-phylo-run.py (UPDATE: but then I had to start changing the uids as well, so, oh, well)

    # ----------------------------------------------------------------------------------------
    def set_unique_ids(self, line, reco_id_str, irandom=None):
        unique_id_columns = ['seqs', 'input_seqs']
        def randstr():  # gah, duplicated in set_reco_id()
            return str(numpy.random.uniform() if irandom is None else irandom)
        uidstrs = [''.join([str(line[c][iseq]) for c in unique_id_columns]) for iseq in range(len(line['input_seqs']))]
        uidstrs = [reco_id_str + uidstrs[iseq] + randstr() + str(iseq) for iseq in range(len(uidstrs))]  # NOTE i'm not sure I really like having the str(iseq), but it mimics the way things used to be by accident/bug (i.e. identical sequences in the same simulated rearrangement event get different uids), so I'm leaving it in for the moment to ease transition after a rewrite
        line['unique_ids'] = [str(hash(ustr)) for ustr in uidstrs]

    # ----------------------------------------------------------------------------------------
    def set_ids(self, line, irandom=None):
        reco_id_str = self.set_reco_id(line, irandom=irandom)
        self.set_unique_ids(line, reco_id_str, irandom=irandom)

    # ----------------------------------------------------------------------------------------
    def set_tree(self, treestr):
        self.tree = treeutils.get_dendro_tree(treestr=treestr)

    # ----------------------------------------------------------------------------------------
    def setline(self, irandom=None):  # don't access <self.line> directly
        if self.line is not None:
            return self.line

        line = {}
        for region in utils.regions:
            line[region + '_gene'] = self.genes[region]
        for boundary in utils.boundaries:
            line[boundary + '_insertion'] = self.insertions[boundary]
        for boundary in utils.effective_boundaries:
            line[boundary + '_insertion'] = ''  # NOTE 'fv' and 'jf' insertions are hereby hardcoded to zero (I'm just writing this here to make it easily searchable -- I don't remember why it's set up that way)
        for erosion in utils.real_erosions:
            line[erosion + '_del'] = self.erosions[erosion]
        for erosion in utils.effective_erosions:
            line[erosion + '_del'] = self.effective_erosions[erosion]
        line['input_seqs'] = self.final_seqs
        line['indelfos'] = self.indelfos
        line['seqs'] = [self.indelfos[iseq]['reversed_seq'] if indelutils.has_indels(self.indelfos[iseq]) else line['input_seqs'][iseq] for iseq in range(len(line['input_seqs']))]
        self.set_ids(line, irandom=irandom)
        treeutils.translate_labels(self.tree, zip(self.leaf_names, line['unique_ids']))  # ordering in <self.leaf_names> is set in recombinator.add_mutants()
        line['tree'] = self.tree.as_string(schema='newick')

        utils.add_implicit_info(self.glfo, line)

        self.line = line

    # ----------------------------------------------------------------------------------------
    def print_gene_choice(self):
        print '    chose:  gene             length'
        for region in utils.regions:
            print '        %s  %-18s %-3d' % (region, self.genes[region], len(self.original_seqs[region])),
            if region in self.pre_erosion_codon_positions:
                print ' (%s: %d)' % (utils.conserved_codons[self.glfo['locus']][region], self.pre_erosion_codon_positions[region])
            else:
                print ''

    # ----------------------------------------------------------------------------------------
    def revert_conserved_codons(self, seq, debug=False):
        """ revert conserved cysteine and tryptophan to their original bases, eg if they were messed up by s.h.m. """
        for region, pos in self.post_erosion_codon_positions.items():  #  NOTE this happens *before* shm indels, i.e. we use self.post_erosion_codon_positions rather than self.final_codon_positions
            if seq[pos : pos + 3] != self.unmutated_codons[region]:
                assert len(self.unmutated_codons[region]) == 3
                if debug:
                    print '    reverting %s --> %s' % (seq[pos : pos + 3], self.unmutated_codons[region])
                seq = seq[:pos] + self.unmutated_codons[region] + seq[pos + 3 :]
            assert utils.codon_unmutated(utils.conserved_codons[self.glfo['locus']][region], seq, pos)
        return seq
