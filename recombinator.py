""" Simulates the process of VDJ recombination """ 
import sys
import csv
import json
import random
import numpy
import math
import os
from Bio import SeqIO

from opener import opener
import recoutils as util

#----------------------------------------------------------------------------------------
class RecombinationEvent(object):
    """ Container to hold the information for a single recombination event. """

    def __init__(self):
        self.vdj_combo_label = ()  # A tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length, <erosion lengths>)
                                   # NOTE I leave the lengths in here as strings
        self.gene_names = {}
        self.seqs = {}
        self.cyst_position = -1
        self.tryp_position = -1  # NOTE this is the position *within* the j gene *only*
        self.final_tryp_position = -1  # while *this* is the tryp position in the final recombined sequence
        self.erosions = {}  # erosion lengths for the event
        self.cdr3_length = 0  # NOTE this is the *desired* cdr3_length, i.e. after erosion and insertion
        self.current_cdr3_length = 0  # while this is the lenght beforehand
        self.net_length_change = 0  # and this is the difference between the two
        self.total_deletion_length = 0
        self.insertion_lengths = {}
        self.insertions = {}
        self.recombined_seq = ''  # combined sequence *before* mutations
        self.mutations = []
        self.final_seqs = []

    def add_empty_mutant(self):
        """ Prepare to add a new mutant. """
        self.mutations.append('')

    def add_mutation(self, position, old_nuke, new_nuke):
        """ Add a single mutation to the list of mutations. """
        if len(self.mutations[-1]) != 0:  # add a separator between adjacent mutations
            self.mutations[-1] += '-'
        self.mutations[-1] += '%d%s%s' % (position, old_nuke, new_nuke)

    def set_vdj_combo(self, vdj_combo_label, cyst_position, tryp_position, all_seqs):
        """ Set the label which labels the gene/length choice (a tuple of strings)
        as well as it's constituent parts. """
        self.vdj_combo_label = vdj_combo_label
        self.cyst_position = cyst_position
        self.tryp_position = tryp_position
        for region in util.regions:
            self.gene_names[region] = vdj_combo_label[util.index_keys[region + '_gene']]
            self.seqs[region] = all_seqs[region][vdj_combo_label[util.index_keys[region + '_gene']]]
        # set the erosion lengths
        for erosion_location in util.erosions:
            self.erosions[erosion_location] = int(vdj_combo_label[util.index_keys[erosion_location + '_del']])
        # set the desired cdr3_length
        self.cdr3_length = int(vdj_combo_label[util.index_keys['cdr3_length']])
        # find 'current' cdr3 length (i.e. with no insertions or erosions)
        tryp_position_in_joined_seq = util.find_tryp_in_joined_seq(self.tryp_position,
                                                                   self.seqs['v'],
                                                                   '',
                                                                   self.seqs['d'],
                                                                   '',
                                                                   self.seqs['j'],
                                                                   0)
        util.check_conserved_codons(self.seqs['v'] + self.seqs['d'] + self.seqs['j'], self.cyst_position, tryp_position_in_joined_seq)
        self.current_cdr3_length = tryp_position_in_joined_seq - self.cyst_position + 3
        self.net_length_change = self.cdr3_length - self.current_cdr3_length
        
        self.total_deletion_length = 0
        for erosion_location in util.erosions:
            self.total_deletion_length += self.erosions[erosion_location]

    def set_final_tryp_position(self):
        """ Set tryp position in the final, combined sequence. """
        self.final_tryp_position = util.find_tryp_in_joined_seq(self.tryp_position,
                                                                self.seqs['v'],
                                                                self.insertions['vd'],
                                                                self.seqs['d'],
                                                                self.insertions['dj'],
                                                                self.seqs['j'],
                                                                self.erosions['j_5p'])
        print '  final tryptophan position: %d' % self.final_tryp_position
        # make sure cdr3 length matches the desired length in vdj_combo_label
        final_cdr3_length = self.final_tryp_position - self.cyst_position + 3
        assert final_cdr3_length == int(self.cdr3_length)


    def write_event(self, outfile):
        """ Write out all info to csv file. """
        columns = ('unique_id', 'reco_id', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'vd_insertion', 'dj_insertion', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'mutations', 'seq')
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
            for region in util.regions:
                row[region + '_gene'] = self.gene_names[region]
            for boundary in util.boundaries:
                row[boundary + '_insertion'] = self.insertions[boundary]
            for erosion_location in util.erosions:
                row[erosion_location + '_del'] = self.erosions[erosion_location]
            # hash the information that uniquely identifies each recombination event
            reco_id = ''
            for column in row:
                assert 'unique_id' not in row
                assert 'mutations' not in row
                assert 'seq' not in row
                reco_id += str(row[column])
            row['reco_id'] = hash(reco_id)
            # then the stuff that's particular to each mutant/clone
            for imute in range(len(self.final_seqs)):
                row['seq'] = self.final_seqs[imute]
                row['mutations'] = self.mutations[imute]
                unique_id = ''  # Hash to uniquely identify the sequence.
                for column in row:
                    unique_id += str(row[column])
                unique_id += str(numpy.random.uniform())
                row['unique_id'] = hash(unique_id)
                
                # write the row
                writer.writerow(row)

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """
    def __init__(self):
        # parameters that control recombination, erosion, and whatnot
        self.mute_rate = 0.06  # average number of point mutations per base
        self.mean_insertion_length = 6  # mean length of the non-templated insertion
        self.mean_n_clones = 5  # mean number of sequences to toss from each rearrangement event
    
        self.all_seqs = {}  # all the Vs, all the Ds...
        self.index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
        self.version_freq_table = {}  # list of the probabilities with which each VDJ combo appears in data

        print '  init'
        print '    reading vdj versions'
        for region in util.regions:
            self.read_vdj_versions(region, 'data/igh'+region+'.fasta')
        print '    reading version freqs'
        self.read_vdj_version_freqs('data/human-beings/C/N/01-C-N_filtered.vdjcdr3.probs.csv.bz2')
        print '    reading cyst and tryp positions'
        with opener('r')('data/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with opener('r')('data/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

    def combine(self, outfile=''):
        """ Run the combination. """
        reco_event = RecombinationEvent()
        print '      choosing genes %45s %10s %10s %10s %10s' % (' ', 'cdr3', 'deletions', 'net', 'ok?')
        while self.are_erosion_lengths_inconsistent(reco_event):
            self.choose_vdj_combo(reco_event)  # set a vdj/erosion choice in reco_event
        self.get_insertion_lengths(reco_event)

        print '    chose:  gene             length'
        for region in util.regions:
            print '        %s  %-18s %-3d' % (region, reco_event.gene_names[region], len(reco_event.seqs[region])),
            if region == 'v':
                print ' (cysteine: %d)' % reco_event.cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % reco_event.tryp_position
            else:
                print ''

        assert not util.is_erosion_longer_than_seq(reco_event)
        assert not util.would_erode_conserved_codon(reco_event)
        # erode, insert, and combine
        self.erode_and_insert(reco_event)
        print '  joining'
        print '         v: %s' % reco_event.seqs['v']
        print '    insert: %s' % reco_event.insertions['vd']
        print '         d: %s' % reco_event.seqs['d']
        print '    insert: %s' % reco_event.insertions['dj']
        print '         j: %s' % reco_event.seqs['j']
        reco_event.recombined_seq = reco_event.seqs['v'] + reco_event.insertions['vd'] + reco_event.seqs['d'] + reco_event.insertions['dj'] + reco_event.seqs['j']
        reco_event.set_final_tryp_position()

        # toss a bunch of clones: add point mutations
        for iclone in range(int(round(numpy.random.poisson(self.mean_n_clones)))):
            print '  clone %d' % iclone
            self.add_mutant(reco_event)
            util.check_conserved_codons(reco_event.final_seqs[-1], reco_event.cyst_position, reco_event.final_tryp_position)

#        # write some stuff that can be used by hmmer for training profiles
#        # NOTE at the moment I have this *appending* to the files
#        self.write_final_vdj(reco_event)

        # write final output to csv
        if outfile != '':
            print '  writing'
            reco_event.write_event(outfile)

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
        in data. This file was created with versioncounter.py
        """
#        with bz2.BZ2File(fname) as infile:
        with opener('r')(fname) as infile:
            in_data = csv.DictReader(infile)
            total = 0.0  # check that the probs sum to 1.0
            for line in in_data:
                total += float(line['prob'])
                index = tuple(line[column] for column in util.index_columns)
                assert index not in self.version_freq_table
                self.version_freq_table[index] = float(line['prob'])
            assert math.fabs(total - 1.0) < 1e-8

    def choose_vdj_combo(self, reco_event):
        """ Choose which combination germline variants to use """
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        for vdj_choice in self.version_freq_table:
            sum_prob += self.version_freq_table[vdj_choice]
            if iprob < sum_prob:
                reco_event.set_vdj_combo(vdj_choice,
                                         self.cyst_positions[vdj_choice[util.index_keys['v_gene']]]['cysteine-position'],
                                         int(self.tryp_positions[vdj_choice[util.index_keys['j_gene']]]),
                                         self.all_seqs)
                return

        assert False  # shouldn't fall through to here

    def get_insertion_lengths(self, reco_event):
        """ Partition the necessary insertion length between the vd and dj boundaries. """
        # first get total insertion length
        total_insertion_length = reco_event.total_deletion_length + reco_event.net_length_change
        assert total_insertion_length >= 0

        # then divide total_insertion_length into vd_insertion and dj_insertion
        partition_point = numpy.random.uniform(0, total_insertion_length)
        reco_event.insertion_lengths['vd'] = int(round(partition_point))
        reco_event.insertion_lengths['dj'] = total_insertion_length - reco_event.insertion_lengths['vd']
        print '      insertion lengths: %d %d' % (reco_event.insertion_lengths['vd'], reco_event.insertion_lengths['dj'])
        assert reco_event.insertion_lengths['vd'] + reco_event.insertion_lengths['dj'] == total_insertion_length  # check for rounding problems

    def erode(self, region, location, reco_event, protected_position=-1):
        """ Erode some number of letters from seq.

        Nucleotides are removed from the <location> ('5p' or '3p') side of
        <seq>. The codon beginning at index <protected_position> is optionally
        protected from removal.
        """
        seq = reco_event.seqs[region]
        n_to_erode = reco_event.erosions[region + '_' + location]
        if protected_position > 0:  # this check is redundant at this point
            if location == '3p' and region == 'v':
                if len(seq) - n_to_erode <= protected_position + 2:
                    assert False
            elif location == '5p' and region == 'j':
                if n_to_erode - 1 >= protected_position:
                    assert False
            else:
                print 'ERROR unanticipated protection'
                sys.exit()

        fragment_before = ''
        fragment_after = ''
        if location == '5p':
            fragment_before = seq[:n_to_erode + 3] + '...'
            new_seq = seq[n_to_erode:len(seq)]
            fragment_after = new_seq[:n_to_erode + 3] + '...'
        elif location == '3p':
            fragment_before = '...' + seq[len(seq) - n_to_erode - 3 :]
            new_seq = seq[0:len(seq)-n_to_erode]
            fragment_after = '...' + new_seq[len(new_seq) - n_to_erode - 3 :]
        else:
            print 'ERROR location must be \"5p\" or \"3p\"'
            sys.exit()
        print '    %3d from %s' % (n_to_erode, location),
        print 'of %s: %15s' % (region, fragment_before),
        print ' --> %-15s' % fragment_after
        if len(fragment_after) == 0:
            print '    NOTE eroded away entire sequence'

        reco_event.seqs[region] = new_seq

    def get_insertion(self, boundary, reco_event):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        insert_seq_str = ''
        for _ in range(0, reco_event.insertion_lengths[boundary]):
            insert_seq_str += util.int_to_nucleotide(random.randint(0, 3))
        reco_event.insertions[boundary] = insert_seq_str
        
    def erode_and_insert(self, reco_event):
        """ Erode and insert based on the lengths in reco_event. """
        print '  eroding'
        self.erode('v', '3p', reco_event, reco_event.cyst_position)
        self.erode('d', '5p', reco_event)
        self.erode('d', '3p', reco_event)
        self.erode('j', '5p', reco_event, reco_event.tryp_position)

        # then insert
        for boundary in util.boundaries:
            self.get_insertion(boundary, reco_event)

    def add_mutant(self, reco_event):
        """ Apply point mutations to a copy of reco_event.recombined_seq
        and append the result to reco_event.final_seqs.
        """
        seq = reco_event.recombined_seq[:]  # make sure it's a copy, not a ref, because it gets mutated
        protected_positions = (reco_event.cyst_position, reco_event.final_tryp_position)
        n_mutes = int(round(numpy.random.poisson(self.mute_rate*len(seq))))
        print '    making %d point mutations: ' % n_mutes,
        reco_event.add_empty_mutant()
        mute_locations = []  # list to help printing below
        for _ in range(n_mutes):
            position = random.randint(0, len(seq)-1)
            while util.is_position_protected(protected_positions, position):
                position = random.randint(0, len(seq)-1)
            mute_locations.append(position)
            old_nuke = seq[position]
            new_nuke = old_nuke
            while new_nuke == old_nuke:
                new_nuke = util.int_to_nucleotide(random.randint(0,3))
            print '%d%s%s' % (position, old_nuke, new_nuke),
            reco_event.add_mutation(position, old_nuke, new_nuke)
            seq = seq[:position] + new_nuke + seq[position+1:]

        print ''
        reco_event.final_seqs.append(seq)

        # print out mutation locations
        mute_locations = set(mute_locations)
        mute_print_str = ''
        for ich in range(len(reco_event.final_seqs[-1])):
            if ich in mute_locations:
                if util.is_position_protected(protected_positions, ich):
                    print 'ERROR mutated a protected position'
                    sys.exit()
                mute_print_str += '|'
            elif util.is_position_protected(protected_positions, ich):
                mute_print_str += 'o'
            else:
                mute_print_str += ' '
        print '    before mute:', reco_event.recombined_seq
        print '    mutations:  ', mute_print_str,'    (o: protected codons)'
        print '    after mute: ', reco_event.final_seqs[-1]

    def write_final_vdj(self, reco_event):
        """ Write the eroded and mutated v, d, and j regions to file. """
        # first do info for the whole reco event
        original_seqs = {}
        for region in util.regions:
            original_seqs[region] = self.all_seqs[region][reco_event.gene_names[region]]
        # work out the final starting positions and lengths
        v_start = 0
        v_length = len(original_seqs['v']) - reco_event.erosions['v_3p']
        d_start = v_length + len(reco_event.insertions['vd'])
        d_length = len(original_seqs['d']) - reco_event.erosions['d_5p'] - reco_event.erosions['d_3p']
        j_start = v_length + len(reco_event.insertions['vd']) + d_length + len(reco_event.insertions['dj'])
        j_length = len(original_seqs['j']) - reco_event.erosions['j_5p']
        # then do stuff that's particular to each mutant
        for final_seq in reco_event.final_seqs:
            assert len(final_seq) == v_length + len(reco_event.insertions['vd']) + d_length + len(reco_event.insertions['dj']) + j_length
            # get the final seqs (i.e. what v, d, and j look like in the final sequence)
            final_seqs = {}
            final_seqs['v'] = final_seq[v_start:v_start+v_length]
            final_seqs['d'] = final_seq[d_start:d_start+d_length]
            final_seqs['j'] = final_seq[j_start:j_start+j_length]
            # pad with dots so it looks like (ok, is) an m.s.a. file
            final_seqs['v'] = final_seqs['v'] + reco_event.erosions['v_3p'] * '.'
            final_seqs['d'] = reco_event.erosions['d_5p'] * '.' + final_seqs['d'] + reco_event.erosions['d_3p'] * '.'
            final_seqs['j'] = reco_event.erosions['j_5p'] * '.' + final_seqs['j']
            for region in util.regions:
                sanitized_name = reco_event.gene_names[region]  # replace special characters in gene names
                sanitized_name = sanitized_name.replace('*','-')
                sanitized_name = sanitized_name.replace('/','-')
                out_dir = 'data/msa/' + region
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                out_fname = out_dir + '/' + sanitized_name + '.sto'
                with opener('ab')(out_fname) as outfile:
                    outfile.write('%15d   %s\n' % (hash(numpy.random.uniform()), final_seqs[region]))

    def are_erosion_lengths_inconsistent(self, reco_event):
        """ Are the erosion lengths inconsistent with the cdr3 length? """
        if reco_event.vdj_combo_label == ():  # haven't filled it yet
            return True
        # now are the erosion lengths we chose consistent with the cdr3_length we chose?
        total_deletion_length = 0
        for erosion in util.erosions:
            total_deletion_length += int(reco_event.vdj_combo_label[util.index_keys[erosion + '_del']])

        # print some crap
        gene_choices = reco_event.vdj_combo_label[util.index_keys['v_gene']] + ' ' + reco_event.vdj_combo_label[util.index_keys['d_gene']] + ' ' + reco_event.vdj_combo_label[util.index_keys['j_gene']]
        print '               try: %45s %10s %10d %10d' % (gene_choices, reco_event.vdj_combo_label[util.index_keys['cdr3_length']], total_deletion_length, reco_event.net_length_change),
        if -total_deletion_length > reco_event.net_length_change:
            print '%10s' % 'no'
        else:
            print '%10s' % 'yes'

        # i.e. we're *in*consistent if net change is negative and also less than total deletions
        return -total_deletion_length > reco_event.net_length_change
