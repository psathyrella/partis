""" Simulates the process of VDJ recombination """
import sys
import csv
import time
import json
import random
from cStringIO import StringIO
import treegenerator
import numpy
import os
import re
from subprocess import check_output

from Bio import SeqIO, Phylo
import dendropy

from opener import opener
import paramutils
import utils
from event import RecombinationEvent

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """
    def __init__(self, args, seed, sublabel=None):
        self.args = args

        if sublabel == None:
            self.workdir = self.args.workdir + '/recombinator'
            self.outfname = self.args.outfname
        else:  # need a separate workdir for each subprocess
            self.workdir = self.args.workdir + '/recombinator-' + sublabel
            self.outfname = self.workdir + '/' + os.path.basename(self.args.outfname)
        utils.prep_dir(self.workdir)

        if not self.args.simulate_partially_from_scratch:
            parameter_dir = self.args.parameter_dir
        else:  # we start from scratch, except for the mute freq stuff
            parameter_dir = self.args.scratch_mute_freq_dir

        if parameter_dir is None or not os.path.exists(parameter_dir):
            raise Exception('parameter dir ' + parameter_dir + ' d.n.e')

        self.index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
        self.mute_models = {}
        # self.treeinfo = []  # list of newick-formatted tree strings with region-specific branch info tacked at the end
        for region in utils.regions:
            self.mute_models[region] = {}
            for model in ['gtr', 'gamma']:
                self.mute_models[region][model] = {}

        self.glfo = utils.read_germline_set(self.args.datadir)

        self.allowed_genes = self.get_allowed_genes(parameter_dir)  # set of genes a) for which we read per-position mutation information and b) from which we choose when running partially from scratch
        self.version_freq_table = self.read_vdj_version_freqs(parameter_dir)  # list of the probabilities with which each VDJ combo (plus other rearrangement parameters) appears in data
        self.insertion_content_probs = self.read_insertion_content(parameter_dir)
        self.all_mute_freqs = {}
        self.parameter_dir = parameter_dir  # damnit, I guess I do need to save this in self

        # read shm info NOTE I'm not inferring the gtr parameters a.t.m., so I'm just (very wrongly) using the same ones for all individuals
        with opener('r')(self.args.gtrfname) as gtrfile:  # read gtr parameters
            reader = csv.DictReader(gtrfile)
            for line in reader:
                parameters = line['parameter'].split('.')
                region = parameters[0][3].lower()
                assert region == 'v' or region == 'd' or region == 'j'
                model = parameters[1].lower()
                parameter_name = parameters[2]
                assert model in self.mute_models[region]
                self.mute_models[region][model][parameter_name] = line['value']
        treegen = treegenerator.TreeGenerator(args, parameter_dir, seed=seed)
        self.treefname = self.workdir + '/trees.tre'
        treegen.generate_trees(seed, self.treefname)
        with opener('r')(self.treefname) as treefile:  # read in the trees (and other info) that we just generated
            self.treeinfo = treefile.readlines()
        if not self.args.no_clean:
            os.remove(self.treefname)

        if os.path.exists(self.outfname):
            os.remove(self.outfname)
        elif not os.path.exists(os.path.dirname(os.path.abspath(self.outfname))):
            os.makedirs(os.path.dirname(os.path.abspath(self.outfname)))

    # ----------------------------------------------------------------------------------------
    def read_insertion_content(self, parameter_dir):
        if self.args.simulate_partially_from_scratch:
            return {b : {n : 1./len(utils.nukes) for n in utils.nukes} for b in utils.boundaries}

        insertion_content_probs = {}
        for bound in utils.boundaries:
            insertion_content_probs[bound] = {}
            with opener('r')(parameter_dir + '/' + bound + '_insertion_content.csv') as icfile:
                reader = csv.DictReader(icfile)
                total = 0
                for line in reader:
                    insertion_content_probs[bound][line[bound + '_insertion_content']] = int(line['count'])
                    total += int(line['count'])
                for nuke in utils.nukes:
                    if nuke not in insertion_content_probs[bound]:
                        print '    %s not in insertion content probs, adding with zero' % nuke
                        insertion_content_probs[bound][nuke] = 0
                    insertion_content_probs[bound][nuke] /= float(total)

            assert utils.is_normed(insertion_content_probs[bound])

        return insertion_content_probs


    # ----------------------------------------------------------------------------------------
    def get_mute_freqs(self, gene_or_insert_name):
        if gene_or_insert_name not in self.all_mute_freqs:
            self.read_mute_freq_stuff(gene_or_insert_name)
        return self.all_mute_freqs[gene_or_insert_name]

    # ----------------------------------------------------------------------------------------
    def read_mute_freq_stuff(self, gene_or_insert_name):
        if gene_or_insert_name[:2] in utils.boundaries:
            replacement_genes = utils.find_replacement_genes(self.parameter_dir, min_counts=-1, all_from_region='v')
            self.all_mute_freqs[gene_or_insert_name], _ = paramutils.read_mute_info(self.parameter_dir, this_gene=gene_or_insert_name, approved_genes=replacement_genes)
        else:
            gene_counts = utils.read_overall_gene_probs(self.parameter_dir, only_gene=gene_or_insert_name, normalize=False, expect_zero_counts=True)
            replacement_genes = None
            if gene_counts < self.args.min_observations_to_write:  # if we didn't see it enough, average over all the genes that find_replacement_genes() gives us NOTE if <gene_or_insert_name> isn't in the dict, it's because it's <args.datadir> but not in the parameter dir
                replacement_genes = utils.find_replacement_genes(self.parameter_dir, min_counts=self.args.min_observations_to_write, gene_name=gene_or_insert_name, single_gene=False)
            self.all_mute_freqs[gene_or_insert_name], _ = paramutils.read_mute_info(self.parameter_dir, this_gene=gene_or_insert_name, approved_genes=replacement_genes)

    # ----------------------------------------------------------------------------------------
    def combine(self, irandom):
        """
        Create a recombination event and write it to disk
        <irandom> is used as the seed for the myriad random number calls.
        If combine() is called with the same <irandom>, it will find the same event, i.e. it should be a random number, not just a seed
        """
        if self.args.debug:
            print 'combine (seed %d)' % irandom
        numpy.random.seed(irandom)
        random.seed(irandom)

        reco_event = RecombinationEvent(self.glfo)
        self.choose_vdj_combo(reco_event)
        self.erode_and_insert(reco_event)
        if self.args.debug:
            print '  joining eroded seqs'
            print '         v: %s' % reco_event.eroded_seqs['v']
            print '    insert: %s' % reco_event.insertions['vd']
            print '         d: %s' % reco_event.eroded_seqs['d']
            print '    insert: %s' % reco_event.insertions['dj']
            print '         j: %s' % reco_event.eroded_seqs['j']
        reco_event.recombined_seq = reco_event.eroded_seqs['v'] + reco_event.insertions['vd'] + reco_event.eroded_seqs['d'] + reco_event.insertions['dj'] + reco_event.eroded_seqs['j']
        codons_ok = reco_event.set_final_cyst_tryp_positions(debug=self.args.debug)
        if not codons_ok:
            print '    unproductive event -- rerunning'  # probably a weirdly long v_3p or j_5p deletion
            return False

        self.add_mutants(reco_event, irandom)  # toss a bunch of clones: add point mutations

        if self.args.debug:
            reco_event.print_event()

        # write output to csv
        reco_event.write_event(self.outfname, irandom=irandom)

        return True

    # ----------------------------------------------------------------------------------------
    def get_parameter_dir_genes(self, parameter_dir):
        parameter_dir_genes = set()
        for region in utils.regions:
            col = region + '_gene'
            column_and_deps = [col, ] + utils.column_dependencies[col]
            with open(parameter_dir + '/' + utils.get_parameter_fname(column_and_deps=column_and_deps)) as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    parameter_dir_genes.add(line[region + '_gene'])
        return parameter_dir_genes

    # ----------------------------------------------------------------------------------------
    def get_allowed_genes(self, parameter_dir):
        # first get all the genes that are available
        if self.args.simulate_partially_from_scratch:  # start with all the ones in args.datadir
            tmplist = [self.glfo['seqs'][r].keys() for r in utils.regions]
            allowed_set = set([g for glist in tmplist for g in glist])
        else:  # start with all the ones in the parameter directory
            allowed_set = self.get_parameter_dir_genes(parameter_dir)

        # then, if specified, require that they're also in args.only_genes
        if self.args.only_genes is not None:
            allowed_set = allowed_set & set(self.args.only_genes)

        allowed_genes = {r : [] for r in utils.regions}
        for gene in allowed_set:
            allowed_genes[utils.get_region(gene)].append(gene)

        return allowed_genes

    # ----------------------------------------------------------------------------------------
    def freqtable_index(self, line):
        return tuple(line[column] for column in utils.index_columns)

    # ----------------------------------------------------------------------------------------
    def read_vdj_version_freqs(self, parameter_dir):
        """ Read the frequencies at which various VDJ combinations appeared in data """
        if self.args.simulate_partially_from_scratch:
            return None

        version_freq_table = {}
        with opener('r')(parameter_dir + '/' + utils.get_parameter_fname('all')) as infile:
            in_data = csv.DictReader(infile)
            total = 0.0
            for line in in_data:  # NOTE do *not* assume the file is sorted
                skip = False
                for region in utils.regions:
                    if line[region + '_gene'] not in self.allowed_genes[region]:
                        skip = True
                        break
                if skip:
                    continue
                total += float(line['count'])
                index = self.freqtable_index(line)
                assert index not in version_freq_table
                version_freq_table[index] = float(line['count'])

        if len(version_freq_table) == 0:
            raise Exception('didn\'t find any gene combinations in %s' % fname)

        # then normalize
        test_total = 0.0
        for index in version_freq_table:
            version_freq_table[index] /= total
            test_total += version_freq_table[index]
        assert utils.is_normed(test_total, this_eps=1e-8)
        assert len(version_freq_table) < 1e8  # if it gets *too* large, choose_vdj_combo() below isn't going to work because of numerical underflow. Note there's nothing special about 1e8, it's just that I'm pretty sure we're fine *up* to that point, and once we get beyond it we should think about doing things differently
        return version_freq_table

    # ----------------------------------------------------------------------------------------
    def choose_vdj_combo(self, reco_event):
        """ Choose the set of rearrangement parameters """

        vdj_choice = None
        if self.args.simulate_partially_from_scratch:  # generate an event without using the parameter directory
            tmpline = {}
            for region in utils.regions:
                tmpline[region + '_gene'] = numpy.random.choice(self.allowed_genes[region])
            for effrode in utils.effective_erosions:
                tmpline[effrode + '_del'] = 0
            for effbound in utils.effective_boundaries:
                tmpline[effbound + '_insertion'] = 0
            for erosion in utils.real_erosions:
                tmpline[erosion + '_del'] = numpy.random.geometric(1. / utils.scratch_mean_erosion_lengths[erosion])
            for bound in utils.boundaries:
                tmpline[bound + '_insertion'] = numpy.random.geometric(1. / utils.scratch_mean_insertion_lengths[bound])
            vdj_choice = self.freqtable_index(tmpline)
        else:  # use real parameters from a directory
            iprob = numpy.random.uniform(0, 1)
            sum_prob = 0.0
            for tmpchoice in self.version_freq_table:  # assign each vdj choice a segment of the interval [0,1], and choose the one which contains <iprob>
                sum_prob += self.version_freq_table[tmpchoice]
                if iprob < sum_prob:
                    vdj_choice = tmpchoice
                    break

            assert vdj_choice is not None  # shouldn't fall through to here

        reco_event.set_vdj_combo(vdj_choice, self.glfo, debug=self.args.debug, mimic_data_read_length=self.args.mimic_data_read_length)

    # ----------------------------------------------------------------------------------------
    def erode(self, erosion, reco_event):
        """ apply <erosion> to the germline seqs in <reco_event> """
        seq = reco_event.eroded_seqs[erosion[0]]  # <erosion> looks like [vdj]_[35]p
        n_to_erode = reco_event.erosions[erosion] if erosion in utils.real_erosions else reco_event.effective_erosions[erosion]
        fragment_before = ''  # fragments to print
        fragment_after = ''
        if '5p' in erosion:
            fragment_before = seq[:n_to_erode + 3] + '...'
            new_seq = seq[n_to_erode:len(seq)]
            fragment_after = new_seq[:n_to_erode + 3] + '...'
        else:
            assert '3p' in erosion
            fragment_before = '...' + seq[len(seq) - n_to_erode - 3 :]
            new_seq = seq[0:len(seq)-n_to_erode]
            fragment_after = '...' + new_seq[len(new_seq) - n_to_erode - 3 :]

        if self.args.debug:
            print '    %3d from %s' % (n_to_erode, erosion[2:]),
            print 'of %s: %15s' % (erosion[0], fragment_before),
            print ' --> %-15s' % fragment_after
        if len(fragment_after) == 0:
            print '    NOTE eroded away entire sequence'

        reco_event.eroded_seqs[erosion[0]] = new_seq

    # ----------------------------------------------------------------------------------------
    def insert(self, boundary, reco_event):
        insert_seq_str = ''
        probs = self.insertion_content_probs[boundary]
        for _ in range(0, reco_event.insertion_lengths[boundary]):
            iprob = numpy.random.uniform(0, 1)
            sum_prob = 0.0
            new_nuke = ''  # this is just to make sure I don't fall through the loop over nukes
            for nuke in utils.nukes:  # assign each nucleotide a segment of the interval [0,1], and choose the one which contains <iprob>
                sum_prob += probs[nuke]
                if iprob < sum_prob:
                    new_nuke = nuke
                    break
            assert new_nuke != ''
            insert_seq_str += new_nuke

        reco_event.insertions[boundary] = insert_seq_str

    # ----------------------------------------------------------------------------------------
    def erode_and_insert(self, reco_event):
        """ Erode the germline seqs, and add insertions, based on the info in <reco_event> """
        if self.args.debug:
            print '  eroding'
        for region in utils.regions:
            reco_event.eroded_seqs[region] = reco_event.original_seqs[region]
        for erosion in utils.real_erosions:
            self.erode(erosion, reco_event)
        if self.args.mimic_data_read_length:
            for erosion in utils.effective_erosions:
                self.erode(erosion, reco_event)
        for boundary in utils.boundaries:
            self.insert(boundary, reco_event)

    # ----------------------------------------------------------------------------------------
    def write_mute_freqs(self, region, gene_or_insert_name, seq, reco_event, reco_seq_fname, is_insertion=False):
        """ Read position-by-position mute freqs from disk for <gene_or_insert_name>, renormalize, then write to a file for bppseqgen. """
        mute_freqs = self.get_mute_freqs(gene_or_insert_name)

        rates = []  # list with a relative mutation rate for each position in <seq>
        total = 0.0
        # assert len(mute_freqs) == len(seq)  # only equal length if no erosions NO oh right but mute_freqs only covers areas we could align to...
        left_erosion_length = dict(reco_event.erosions.items() + reco_event.effective_erosions.items())[region + '_5p']
        for inuke in range(len(seq)):  # append a freq for each nuke
            position = inuke + left_erosion_length
            freq = 0.0
            if position in mute_freqs:
                freq = mute_freqs[position]
            else:
                freq = mute_freqs['overall_mean']
            rates.append(freq)
            total += freq

        # normalize to the number of sites (i.e. so an average site is given value 1.0)
        assert total != 0.0  # I am not hip enough to divide by zero
        for inuke in range(len(seq)):
            rates[inuke] *= float(len(seq)) / total
        total = 0.0

        # and... double check it, just for shits and giggles
        for inuke in range(len(seq)):
            total += rates[inuke]
        assert utils.is_normed(total / float(len(seq)))
        assert len(rates) == len(seq)  # you just can't be too careful. what if gremlins ate a few while python wasn't looking?

        # write the input file for bppseqgen, one base per line
        with opener('w')(reco_seq_fname) as reco_seq_file:
            reco_seq_file.write('state\trate\n')
            for inuke in range(len(seq)):
                reco_seq_file.write('%s\t%.15f\n' % (seq[inuke], rates[inuke]))

        # NOTE I need to find a tool to give me the total branch length of the chosen tree, so I can compare to the number of mutations I see

    # ----------------------------------------------------------------------------------------
    def run_bppseqgen(self, seq, chosen_tree, gene_or_insert_name, reco_event, seed, is_insertion=False):
        """ Run bppseqgen on sequence

        Note that this is in general a piece of the full sequence (say, the V region), since
        we have different mutation models for different regions. Returns a list of mutated
        sequences.
        """
        region = ''
        if is_insertion:
            region = 'v'  # NOTE should really do something other than just use the v model for insertion mutations
        else:
            region = utils.get_region(gene_or_insert_name)

        treg = re.compile('t[0-9][0-9]*')  # find number of leaf nodes
        n_leaf_nodes = len(treg.findall(chosen_tree))
        if len(seq) == 0:  # zero length insertion (or d)
            return ['' for _ in range(n_leaf_nodes)]  # return an empty string for each leaf node

        # write the tree to a tmp file
        if is_insertion:
            label = gene_or_insert_name[:2]
        else:
            label = utils.get_region(gene_or_insert_name)
        treefname = self.workdir + '/' + label + '-tree.tre'
        reco_seq_fname = self.workdir + '/' + label + '-start-seq.txt'
        leaf_seq_fname = self.workdir + '/' + label + '-leaf-seqs.fa'
        if n_leaf_nodes == 1:  # add an extra leaf to one-leaf trees so bppseqgen doesn't barf (we'll ignore it later)
            lreg = re.compile('t1:[0-9]\.[0-9][0-9]*')
            leafstr = lreg.findall(chosen_tree)
            assert len(leafstr) == 1
            leafstr = leafstr[0]
            chosen_tree = chosen_tree.replace(leafstr, '(' + leafstr + ',' + leafstr + '):0.0')
        with opener('w')(treefname) as treefile:
            treefile.write(chosen_tree)
        self.write_mute_freqs(region, gene_or_insert_name, seq, reco_event, reco_seq_fname, is_insertion=is_insertion)

        # build up the command line
        # docs: http://biopp.univ-montp2.fr/apidoc/bpp-phyl/html/classbpp_1_1GTR.html that page is too darn hard to google
        bpp_binary = self.args.partis_dir + '/packages/bpp/bin/bppseqgen'
        if not os.path.exists(bpp_binary):
            raise Exception('bpp not found in %s' % os.path.dirname(bpp_binary))
        command = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:' + self.args.partis_dir + '/packages/bpp/lib\n'
        command += bpp_binary
        command += ' input.tree.file=' + treefname
        command += ' output.sequence.file=' + leaf_seq_fname
        command += ' number_of_sites=' + str(len(seq))
        command += ' input.tree.format=Newick'
        command += ' output.sequence.format=Fasta\(\)'
        command += ' alphabet=DNA'
        command += ' --seed=' + str(seed)
        command += ' model=GTR\('
        for par in self.mute_models[region]['gtr']:
            val = self.mute_models[region]['gtr'][par]
            command += par + '=' + val + ','
        command = command.rstrip(',')
        command += '\)'
        # NOTE should I use the "equilibrium frequencies" option?
        command += ' rate_distribution=\'Gamma(n=4,alpha=' + self.mute_models[region]['gamma']['alpha']+ ')\''
        command += ' input.infos.states=state'
        command += ' input.infos=' + reco_seq_fname
        command += ' input.infos.rates=rate'
        check_output(command, shell=True)

        mutated_seqs = []
        for seq_record in SeqIO.parse(leaf_seq_fname, "fasta"):  # get the leaf node sequences from the file that bppseqgen wrote
            mutated_seqs.append(str(seq_record.seq))
            if n_leaf_nodes == 1:  # skip the extra leaf we added above
                break
        assert n_leaf_nodes == len(mutated_seqs)

        # self.check_tree_simulation(leaf_seq_fname, chosen_tree)

        if not self.args.no_clean:
            os.remove(reco_seq_fname)  # clean up temp files
            os.remove(treefname)
            os.remove(leaf_seq_fname)

        return mutated_seqs

    # ----------------------------------------------------------------------------------------
    def get_rescaled_trees(self, treestr, branch_length_ratios, debug=False):
        """
        Trees are generated with the mean branch length observed in data over the whole sequence, because we want to use topologically
        the same tree for the whole sequence. But we observe different branch lengths for each region, so we need to rescale the tree for
        v, d, and j
        """
        rescaled_trees = {}
        if debug:
            print '      rescaling tree:'
        for region in utils.regions:
            # rescale the tree
            rescaled_trees[region] = treegenerator.rescale_tree(treestr, branch_length_ratios[region])
            if debug:
                print '         %s by %f (new depth %f): %s -> %s' % (region, branch_length_ratios[region], treegenerator.get_leaf_node_depths(rescaled_trees[region])['t1'], treestr, rescaled_trees[region])

            # and then check it NOTE can remove this eventually
            initial_depths = {}
            for node, depth in treegenerator.get_leaf_node_depths(treestr).items():
                initial_depths[node] = depth
            for node, depth in treegenerator.get_leaf_node_depths(rescaled_trees[region]).items():
                depth_ratio = depth / initial_depths[node]
                assert utils.is_normed(depth_ratio / branch_length_ratios[region], this_eps=1e-6)
        return rescaled_trees

    # ----------------------------------------------------------------------------------------
    def add_shm_insertion(self, reco_event, seq, pos, length):
        """ insert a random sequence with <length> beginning at <pos> """
        inserted_sequence = ''
        for ipos in range(length):
            inuke = random.randint(0, len(utils.nukes) - 1)  # inclusive
            inserted_sequence += utils.nukes[inuke]
        return_seq = seq[ : pos] + inserted_sequence + seq[pos : ]
        reco_event.indelfos[-1]['indels'].append({'type' : 'insertion', 'pos' : pos, 'len' : length, 'seqstr' : inserted_sequence})
        if self.args.debug:
            print '          inserting %s at %d' % (inserted_sequence, pos)
        return return_seq

    # ----------------------------------------------------------------------------------------
    def add_single_indel(self, seq, reco_event):
        if self.args.indel_location == None:  # uniform over entire sequence
            pos = random.randint(0, len(seq) - 1)  # this will actually exclude either before the first index or after the last index. No, I don't care.
        elif self.args.indel_location == 'v':  # within the meat of the v
            pos = random.randint(10, reco_event.final_cyst_position)
        elif self.args.indel_location == 'cdr3':  # inside cdr3
            pos = random.randint(reco_event.final_cyst_position, reco_event.final_tryp_position)
        else:
            assert False

        length = numpy.random.geometric(1. / self.args.mean_indel_length)

        if numpy.random.uniform(0, 1) < 0.5:  # fifty-fifty chance of insertion and deletion
            new_seq = self.add_shm_insertion(reco_event, seq, pos, length)
        else:
            deleted_seq = seq[ : pos] + seq[pos + length : ]  # delete <length> bases beginning with <pos>
            reco_event.indelfos[-1]['indels'].append({'type' : 'deletion', 'pos' : pos, 'len' : length, 'seqstr' : seq[pos : pos + length]})
            if self.args.debug:
                print '          deleting %d bases at %d' % (length, pos)
            new_seq = deleted_seq

        return new_seq

    # ----------------------------------------------------------------------------------------
    def add_shm_indels(self, reco_event):
        if self.args.debug and self.args.indel_frequency > 0.:
            print '      indels'
        for iseq in range(len(reco_event.final_seqs)):
            reco_event.indelfos.append(utils.get_empty_indel())
            if self.args.indel_frequency == 0.:  # no indels at all
                continue
            if numpy.random.uniform(0, 1) > self.args.indel_frequency:  # no indels for this sequence
                if self.args.debug:
                    print '        0'
                continue
            seq = reco_event.final_seqs[iseq]
            reco_event.indelfos[-1]['reversed_seq'] = seq  # set the original sequence (i.e. with all the indels reversed)
            n_indels = 1  #numpy.random.geometric(1. / self.args.mean_n_indels)
            if self.args.debug:
                print '        %d' % n_indels
            for _ in range(n_indels):
                seq = self.add_single_indel(seq, reco_event)
            reco_event.final_seqs[iseq] = seq

    # ----------------------------------------------------------------------------------------
    def add_mutants(self, reco_event, irandom):
        chosen_treeinfo = self.treeinfo[random.randint(0, len(self.treeinfo)-1)]
        chosen_tree = chosen_treeinfo.split(';')[0] + ';'
        branch_length_ratios = {}  # NOTE a.t.m (and probably permanently) the mean branch lengths for each region are the *same* for all the trees in the file, I just don't have a better place to put them while I'm passing from TreeGenerator to here than at the end of each line in the file
        for tmpstr in chosen_treeinfo.split(';')[1].split(','):  # looks like e.g.: (t2:0.003751736951,t1:0.003751736951):0.001248262937;v:0.98,d:1.8,j:0.87, where the newick trees has branch lengths corresponding to the whole sequence  (i.e. the weighted mean of v, d, and j)
            region = tmpstr.split(':')[0]
            assert region in utils.regions
            ratio = float(tmpstr.split(':')[1])
            if self.args.mutation_multiplier is not None:  # multiply the branch lengths by some factor
                # if self.args.debug:
                # print '    adding branch length factor %f ' % self.args.mutation_multiplier
                ratio *= self.args.mutation_multiplier
            branch_length_ratios[region] = ratio

        if self.args.debug:  # NOTE should be the same for t[0-9]... but I guess I should check at some point
            print '  using tree with total depth %f' % treegenerator.get_leaf_node_depths(chosen_tree)['t1']  # kind of hackey to just look at t1, but they're all the same anyway and it's just for printing purposes...
            if len(re.findall('t', chosen_tree)) > 1:  # if more than one leaf
                Phylo.draw_ascii(Phylo.read(StringIO(chosen_tree), 'newick'))
            else:
                print '    one leaf'
            print '    with branch length ratios ', ', '.join(['%s %f' % (region, branch_length_ratios[region]) for region in utils.regions])

        scaled_trees = self.get_rescaled_trees(chosen_tree, branch_length_ratios)
        # NOTE would be nice to parallelize this
        mutes = {}
        for region in utils.regions:
            mutes[region] = self.run_bppseqgen(reco_event.eroded_seqs[region], scaled_trees[region], reco_event.genes[region], reco_event, seed=irandom, is_insertion=False)
        mutes['vd'] = self.run_bppseqgen(reco_event.insertions['vd'], scaled_trees['v'], 'vd_insert', reco_event, seed=irandom, is_insertion=True)  # NOTE would be nice to use a better mutation model for the insertions
        mutes['dj'] = self.run_bppseqgen(reco_event.insertions['dj'], scaled_trees['j'], 'dj_insert', reco_event, seed=irandom, is_insertion=True)

        assert len(reco_event.final_seqs) == 0
        for iseq in range(len(mutes['v'])):
            seq = mutes['v'][iseq] + mutes['vd'][iseq] + mutes['d'][iseq] + mutes['dj'][iseq] + mutes['j'][iseq]  # build final sequence
            seq = reco_event.revert_conserved_codons(seq)  # if mutation screwed up the conserved codons, just switch 'em back to what they were to start with
            reco_event.final_seqs.append(seq)  # set final sequnce in reco_event

        assert not utils.are_conserved_codons_screwed_up(reco_event)

        self.add_shm_indels(reco_event)

    # ----------------------------------------------------------------------------------------
    def check_tree_simulation(self, leaf_seq_fname, chosen_tree_str, reco_event=None):
        """ See how well we can reconstruct the true tree """
        clean_up = False
        if leaf_seq_fname == '':  # we need to make the leaf seq file based on info in reco_event
            clean_up = True
            leaf_seq_fname = self.workdir + '/leaf-seqs.fa'
            with opener('w')(leaf_seq_fname) as leafseqfile:
                for iseq in range(len(reco_event.final_seqs)):
                    leafseqfile.write('>t' + str(iseq+1) + '\n')  # NOTE the *order* of the seqs doesn't correspond to the tN number. does it matter?
                    leafseqfile.write(reco_event.final_seqs[iseq] + '\n')

        with opener('w')(os.devnull) as fnull:
            inferred_tree_str = check_output('FastTree -gtr -nt ' + leaf_seq_fname, shell=True, stderr=fnull)
        if clean_up and not self.args.no_clean:
            os.remove(leaf_seq_fname)
        chosen_tree = dendropy.Tree.get_from_string(chosen_tree_str, 'newick')
        inferred_tree = dendropy.Tree.get_from_string(inferred_tree_str, 'newick')
        if self.args.debug:
            print '        tree diff -- symmetric %d   euke %f   rf %f' % (chosen_tree.symmetric_difference(inferred_tree), chosen_tree.euclidean_distance(inferred_tree), chosen_tree.robinson_foulds_distance(inferred_tree))
