""" Simulates the process of VDJ recombination """ 
import sys
import csv
import json
import random
import numpy
import bz2
import math
from Bio import SeqIO

from opener import opener
import recoutils as util

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """

    # parameters that control recombination, erosion, and whatnot
    mute_rate = 0.06  # average number of point mutations per base
    mean_insertion_length = 6  # mean length of the non-templated insertion
    mean_n_clones = 5  # mean number of sequences to toss from each rearrangement event

    all_seqs = {}  # all the Vs, all the Ds...
    regions = ['v', 'd', 'j']
    erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
    index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
    version_freq_table = {}  # list of the probabilities with which each VDJ combo appears in data

    def __init__(self):
        """ Initialize from files """
        print '  init'
        print '    reading vdj versions'
        # tuple with the variables that are used to index selection frequencies
        self.index_columns = ('v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del')
        for i in range(len(self.index_columns)):  # dict so we can access them by name instead of by index number
            self.index_keys[self.index_columns[i]] = i
        for region in self.regions:
            self.read_vdj_versions(region, 'data/igh'+region+'.fasta')
        print '    reading version freqs'
        self.read_vdj_version_freqs('data/human-beings/C/N/01-C-N_filtered.vdjcdr3.probs.csv.bz2')
        print '    reading cyst and tryp positions'
        with open('data/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with open('data/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

    def combine(self, outfile=''):
        """ Run the combination. """
        vdj_combo_label = ()  # a tuple with the names of the chosen versions (v_gene, d_gene, j_gene, cdr3_length, <erosion lengths>)
        chosen_seqs = {}
        cyst_position = -1
        tryp_position = -1
        net_length_change = 0
        print '    choosing genes: %45s %10s %10s %10s %10s' % (' ', 'cdr3_length', 'deletions', 'net', 'ok?')
        while self.are_erosion_lengths_inconsistent(vdj_combo_label, net_length_change):
            # first choose a combination
            vdj_combo_label = self.choose_vdj_combo()
            for region in self.regions:
                chosen_seqs[region] = self.all_seqs[region][vdj_combo_label[self.index_keys[region + '_gene']]]

            # then update some parameters we need in order to figure out whether it's consistent
            cyst_position = self.cyst_positions[vdj_combo_label[self.index_keys['v_gene']]]['cysteine-position']
            tryp_position = int(self.tryp_positions[vdj_combo_label[self.index_keys['j_gene']]])
            desired_cdr3_length = int(vdj_combo_label[self.index_keys['cdr3_length']])
            # find 'current' cdr3 length (i.e. with no insertions or erosions)
            tryp_position_in_joined_seq = util.find_tryp_in_joined_seq(tryp_position, chosen_seqs['v'],  '', chosen_seqs['d'], '', chosen_seqs['j'], 0)
            util.check_conserved_codons(chosen_seqs['v'] + chosen_seqs['d'] + chosen_seqs['j'], cyst_position, tryp_position_in_joined_seq)
            current_cdr3_length = tryp_position_in_joined_seq - cyst_position + 3
            net_length_change = desired_cdr3_length - current_cdr3_length
    
        reco_info = {}  # collect information about the recombination process for output to csv file
        reco_info['vdj_combo'] = vdj_combo_label
        total_deletion_length = 0
        for erosion in self.erosions:
            deletion = int(vdj_combo_label[self.index_keys[erosion + '_del']])
            total_deletion_length += deletion
            reco_info[erosion + '_del'] = deletion
        self.get_insertion_lengths(reco_info, total_deletion_length, net_length_change)

        print '    chose:  gene             length'
        for region in self.regions:
            print '        %s  %-18s %-3d' % (region, vdj_combo_label[self.index_keys[region + '_gene']], len(chosen_seqs[region])),
            if region == 'v':
                print ' (cysteine: %d)' % cyst_position
            elif region == 'j':
                print ' (tryptophan: %d)' % tryp_position
            else:
                print ''

        assert not util.is_erosion_longer_than_seq(reco_info, chosen_seqs)
        assert not util.would_erode_conserved_codon(reco_info, chosen_seqs, cyst_position, tryp_position)
        # erode, insert, and combine
        self.erode_and_insert(reco_info, chosen_seqs, cyst_position, tryp_position)
        print '  joining'
        print '         v: %s' % chosen_seqs['v']
        print '    insert: %s' % reco_info['vd_insertion']
        print '         d: %s' % chosen_seqs['d']
        print '    insert: %s' % reco_info['dj_insertion']
        print '         j: %s' % chosen_seqs['j']
        recombined_seq = chosen_seqs['v'] + reco_info['vd_insertion'] + chosen_seqs['d'] + reco_info['dj_insertion'] + chosen_seqs['j']
        final_tryp_position = util.find_tryp_in_joined_seq(tryp_position,  # NOTE remember this is the tryp position, in j alone, *before* erosion
                                                      chosen_seqs['v'],
                                                      reco_info['vd_insertion'],
                                                      chosen_seqs['d'],
                                                      reco_info['dj_insertion'],
                                                      chosen_seqs['j'],
                                                      reco_info['j_5p_del'])
        print '  final tryptophan position: %d' % final_tryp_position
        # make sure cdr3 length matches the desired length in vdj_combo_label
        final_cdr3_length = final_tryp_position - cyst_position + 3
        assert final_cdr3_length == int(vdj_combo_label[self.index_keys['cdr3_length']])

        # toss a bunch of clones: add point mutations
        for ic in range(int(round(numpy.random.poisson(self.mean_n_clones)))):
            print '  clone %d' % ic
            reco_info['seq'] = self.mutate(reco_info, recombined_seq, (cyst_position, final_tryp_position))
            util.check_conserved_codons(reco_info['seq'], cyst_position, final_tryp_position)

            # write some stuff that can be used by hmmer for training profiles
#            self.write_final_vdj(vdj_combo_label, reco_info)

            if outfile != '':
                self.write_csv(outfile, reco_info)

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
        with bz2.BZ2File(fname) as infile:
            in_data = csv.DictReader(infile)
            total = 0.0  # check that the probs sum to 1.0
            for line in in_data:
                total += float(line['prob'])
                index = tuple(line[column] for column in self.index_columns)
                assert index not in self.version_freq_table
                self.version_freq_table[index] = float(line['prob'])
            assert math.fabs(total - 1.0) < 1e-8

    def choose_vdj_combo(self):
        """ Choose which combination germline variants to use """
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        for vdj_choice in self.version_freq_table:
            sum_prob += self.version_freq_table[vdj_choice]
            if iprob < sum_prob:
                return vdj_choice
        assert False  # shouldn't fall through to here

    def get_insertion_lengths(self, reco_info, total_deletion_length, net_length_change):
        """ Partition the necessary insertion length between the vd and dj boundaries. """
        # first get total insertion length
        total_insertion_length = total_deletion_length + net_length_change
        assert total_insertion_length >= 0

        # then divide total_insertion_length into vd_insertion and dj_insertion
        partition_point = numpy.random.uniform(0, total_insertion_length)
        reco_info['vd_insertion_length'] = int(round(partition_point))
        reco_info['dj_insertion_length'] = total_insertion_length - reco_info['vd_insertion_length']
        print '      insertion lengths: %d %d' % (reco_info['vd_insertion_length'], reco_info['dj_insertion_length'])
        assert reco_info['vd_insertion_length'] + reco_info['dj_insertion_length'] == total_insertion_length  # check for rounding problems

    def erode(self, region, location, seqs, lengths, protected_position=-1):
        """ Erode some number of letters from seq.

        Nucleotides are removed from the <location> ('5p' or '3p') side of
        <seq>. The codon beginning at index <protected_position> is optionally
        protected from removal.
        """
        seq = seqs[region]
        n_to_erode = lengths[region + '_' + location + '_del']
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
        return new_seq

    def get_insertion(self, length):
        """ Get the non-templated sequence to insert between
        templated regions
        """
        insert_seq_str = ''
        for _ in range(0, length):
            insert_seq_str += util.int_to_nucleotide(random.randint(0, 3))
        return insert_seq_str
        
    def erode_and_insert(self, reco_info, seqs, cyst_position, tryp_position):
        """ Erode and insert based on the lengths in reco_info. """
        print '  eroding'
        seqs['v'] = self.erode('v', '3p', seqs, reco_info, cyst_position)
        seqs['d'] = self.erode('d', '5p', seqs, reco_info)
        seqs['d'] = self.erode('d', '3p', seqs, reco_info)
        seqs['j'] = self.erode('j', '5p', seqs, reco_info, tryp_position)
        for boundary in ('vd', 'dj'):
            reco_info[boundary + '_insertion'] = self.get_insertion(reco_info[boundary + '_insertion_length'])

    def mutate(self, reco_info, seq, protected_positions):
        """ Apply point mutations to <seq>, then return it. """
        n_mutes = int(round(numpy.random.poisson(self.mute_rate*len(seq))))
        original_seq = seq
        print '    making %d point mutations: ' % n_mutes,
        mute_locations = []
        reco_info['mutations'] = ''
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
            reco_info['mutations'] += '%d%s%s-' % (position, old_nuke, new_nuke)
            seq = seq[:position] + new_nuke + seq[position+1:]
        print ''
        reco_info['mutations'] = reco_info['mutations'].rstrip('-')

        # print out mutation locations
        mute_locations = set(mute_locations)
        mute_print_str = ''
        for ich in range(len(seq)):
            if ich in mute_locations:
                if util.is_position_protected(protected_positions, ich):
                    print 'ERROR mutated a protected position'
                    sys.exit()
                mute_print_str += '|'
            elif util.is_position_protected(protected_positions, ich):
                mute_print_str += 'o'
            else:
                mute_print_str += ' '
        print '    before mute:', original_seq
        print '    mutations:  ', mute_print_str,'    (o: protected codons)'
        print '    after mute: ', seq

        return seq

    def write_final_vdj(self, vdj_combo_label, reco_info):
        """ Write the eroded and mutated v, d, and j regions to file. """
        original_seqs = {}
        for region in self.regions:
            original_seqs[region] = self.all_seqs[region][vdj_combo_label[self.index_keys[region]]]
        # work out the final starting positions and lengths
        v_start = 0
        v_length = len(original_seqs['v']) - reco_info['v_3p_del']
        d_start = v_length + len(reco_info['vd_insertion'])
        d_length = len(original_seqs['d']) - reco_info['d_5p_del'] - reco_info['d_3p_del']
        j_start = v_length + len(reco_info['vd_insertion']) + d_length + len(reco_info['dj_insertion'])
        j_length = len(original_seqs['j']) - reco_info['j_5p_del']
        assert len(reco_info['seq']) == v_length + len(reco_info['vd_insertion']) + d_length + len(reco_info['dj_insertion']) + j_length
        # get the final seqs (i.e. what v, d, and j look like in the final sequence)
        final_seqs = {}
        final_seqs['v'] = reco_info['seq'][v_start:v_start+v_length]
        final_seqs['d'] = reco_info['seq'][d_start:d_start+d_length]
        final_seqs['j'] = reco_info['seq'][j_start:j_start+j_length]
        # pad with dots so it looks like (ok, is) an m.s.a. file
        final_seqs['v'] = final_seqs['v'] + reco_info['v_3p_del'] * '.'
        final_seqs['d'] = reco_info['d_5p_del'] * '.' + final_seqs['d'] + reco_info['d_3p_del'] * '.'
        final_seqs['j'] = reco_info['j_5p_del'] * '.' + final_seqs['j']
        for region in self.regions:
            sanitized_name = vdj_combo_label[self.index_keys[region]]  # replace special characters in gene names
            sanitized_name = sanitized_name.replace('*','-')
            sanitized_name = sanitized_name.replace('/','-')
            out_fname = 'data/msa/' + sanitized_name + '.sto'
            with open(out_fname, 'ab') as outfile:
                outfile.write('%15d   %s\n' % (hash(numpy.random.uniform()), final_seqs[region]))

    def write_csv(self, outfile, reco_info):
        """ Write out all info to csv file. """
        with open(outfile, 'ab') as csvfile:
            csv_writer = csv.writer(csvfile)
            columns = ('vdj_combo', 'vd_insertion', 'dj_insertion', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'mutations', 'seq')
            string_to_hash = ''  # hash the information that uniquely identifies each recombination event
            string_to_hash_unique = ''  # Hash to *uniquely* identify the sequence.
                                        #   It's just strings of all the column values mashed together
                                        #   with a random number tossed on top.
            csv_row = []
            for column in columns:
                if column == 'vdj_combo':
                    for i in range(4):
                        csv_row.append(reco_info[column][i])
                        string_to_hash += str(reco_info[column][i])
                        string_to_hash_unique += str(reco_info[column][i])
                else:
                    csv_row.append(reco_info[column])
                    string_to_hash_unique += str(reco_info[column])
                    if column != 'mutations' and column != 'seq':
                        string_to_hash += str(reco_info[column])

            csv_row.insert(0, hash(string_to_hash))
            string_to_hash_unique += str(numpy.random.uniform())
            csv_row.insert(0, hash(string_to_hash_unique))
            csv_writer.writerow(csv_row)

    def are_erosion_lengths_inconsistent(self, vdj_combo_label, net_length_change):
        """ Are the erosion lengths inconsistent with the cdr3 length? """
        if vdj_combo_label == ():  # haven't filled it yet
            return True
        # now are the erosion lengths we chose consistent with the cdr3_length we chose?
        total_deletion_length = 0
        for erosion in self.erosions:
            total_deletion_length += int(vdj_combo_label[self.index_keys[erosion + '_del']])

        # print some crap
        gene_choices = vdj_combo_label[self.index_keys['v_gene']] + ' ' + vdj_combo_label[self.index_keys['d_gene']] + ' ' + vdj_combo_label[self.index_keys['j_gene']]
        print '            trying: %45s %10s %10d %10d' % (gene_choices, vdj_combo_label[self.index_keys['cdr3_length']], total_deletion_length, net_length_change),
        if -total_deletion_length > net_length_change:
            print '%10s' % 'no'
        else:
            print '%10s' % 'yes'

        # i.e. we're *in*consistent if net change is negative and also less than total deletions
        return -total_deletion_length > net_length_change
