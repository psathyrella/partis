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
import glutils
from event import RecombinationEvent

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """
    def __init__(self, args, glfo, seed, workdir, outfname):  # NOTE <gldir> is not in general the same as <args.initial_germline_dir>
        self.args = args
        self.glfo = glfo

        # NOTE in general *not* the same as <self.args.workdir> and <self.args.outfname>
        self.workdir = workdir
        self.outfname = outfname
        utils.prep_dir(self.workdir)

        if self.args.rearrange_from_scratch:
            parameter_dir = self.args.scratch_mute_freq_dir  # if you make up mute freqs from scratch, unless you're really careful you tend to get nonsense results for a lot of things (e.g. allele finding)
        else:
            parameter_dir = self.args.parameter_dir + '/' + self.args.parameter_type

        if parameter_dir is None or not os.path.exists(parameter_dir):
            raise Exception('parameter dir ' + parameter_dir + ' d.n.e')

        self.index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
        self.mute_models = {}
        # self.treeinfo = []  # list of newick-formatted tree strings with region-specific branch info tacked at the end
        for region in utils.regions:
            self.mute_models[region] = {}
            for model in ['gtr', 'gamma']:
                self.mute_models[region][model] = {}

        self.allele_prevalence_freqs = glutils.read_allele_prevalence_freqs(args.allele_prevalence_fname) if args.allele_prevalence_fname is not None else {}
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
        os.remove(self.treefname)

        if os.path.exists(self.outfname):
            os.remove(self.outfname)
        elif not os.path.exists(os.path.dirname(os.path.abspath(self.outfname))):
            os.makedirs(os.path.dirname(os.path.abspath(self.outfname)))

    # ----------------------------------------------------------------------------------------
    def read_insertion_content(self, parameter_dir):
        if self.args.rearrange_from_scratch:
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
        if self.args.mutate_from_scratch:  # XXX GODDAMMIT i remember putting this 'xxx' here for a reason and I have no fucking clue what it was
            self.all_mute_freqs[gene_or_insert_name] = {'overall_mean' : self.args.flat_mute_freq}
        elif gene_or_insert_name[:2] in utils.boundaries:
            replacement_genes = utils.find_replacement_genes(self.parameter_dir, min_counts=-1, all_from_region='v')
            self.all_mute_freqs[gene_or_insert_name], _ = paramutils.read_mute_info(self.parameter_dir, this_gene=gene_or_insert_name, chain=self.args.chain, approved_genes=replacement_genes)
        else:
            gene_counts = utils.read_overall_gene_probs(self.parameter_dir, only_gene=gene_or_insert_name, normalize=False, expect_zero_counts=True)
            replacement_genes = None
            if gene_counts < self.args.min_observations_to_write:  # if we didn't see it enough, average over all the genes that find_replacement_genes() gives us NOTE if <gene_or_insert_name> isn't in the dict, it's because it's <args.datadir> but not in the parameter dir UPDATE not using datadir like this any more, so previous statement may not be true
                replacement_genes = utils.find_replacement_genes(self.parameter_dir, min_counts=self.args.min_observations_to_write, gene_name=gene_or_insert_name)
            self.all_mute_freqs[gene_or_insert_name], _ = paramutils.read_mute_info(self.parameter_dir, this_gene=gene_or_insert_name, chain=self.args.chain, approved_genes=replacement_genes)

    # ----------------------------------------------------------------------------------------
    def combine(self, initial_irandom):
        """ keep running self.try_to_combine() until you get a good event """
        failed = True
        itry = 0
        while failed:
            if itry > 0 and self.args.debug:
                print '    unproductive event -- rerunning (try %d)  ' % itry  # probably a weirdly long v_3p or j_5p deletion
            failed = not self.try_to_combine(initial_irandom + itry)
            itry += 1

    # ----------------------------------------------------------------------------------------
    def try_to_combine(self, irandom):
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
        reco_event.set_final_codon_positions()

        # set the original conserved codon words, so we can revert them if they get mutated NOTE we do it here, *after* setting the full recombined sequence, so the germline Vs that don't extend through the cysteine don't screw us over
        reco_event.unmutated_codons = {}
        for region, codon in utils.conserved_codons[self.args.chain].items():
            fpos = reco_event.final_codon_positions[region]
            original_codon = reco_event.recombined_seq[fpos : fpos + 3]
            reco_event.unmutated_codons[region] = reco_event.recombined_seq[fpos : fpos + 3]
            # print fpos, original_codon, utils.codon_ok(codon, reco_event.recombined_seq, fpos)

        codons_ok = utils.both_codons_ok(self.glfo['chain'], reco_event.recombined_seq, reco_event.final_codon_positions, extra_str='      ', debug=self.args.debug)
        if not codons_ok:
            if self.args.rearrange_from_scratch and self.args.generate_germline_set:  # if you let it try more than once, it screws up the desired allele prevalence ratios
                raise Exception('arg')
            return False
        in_frame = reco_event.cdr3_length % 3 == 0
        if self.args.rearrange_from_scratch and not in_frame:
            raise Exception('arg 2')  # if you let it try more than once, it screws up the desired allele prevalence ratios
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
        if self.args.rearrange_from_scratch:  # start with all of 'em
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
        if self.args.rearrange_from_scratch:
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
    def get_scratchline(self):
        tmpline = {}
        for region in utils.regions:
            probs = None
            if region in self.allele_prevalence_freqs and len(self.allele_prevalence_freqs[region]) > 0:  # should really change it so it has to be the one or the other
                probs = [self.allele_prevalence_freqs[region][g] for g in self.allowed_genes[region]]
            tmpline[region + '_gene'] = numpy.random.choice(self.allowed_genes[region], p=probs)
        for effrode in utils.effective_erosions:
            tmpline[effrode + '_del'] = 0
        for effbound in utils.effective_boundaries:
            tmpline[effbound + '_insertion'] = ''

        # ----------------------------------------------------------------------------------------
        def try_scratch_erode_insert(tmpline):
            utils.remove_all_implicit_info(tmpline)
            for erosion in utils.real_erosions:  # includes various contortions to avoid eroding the entire gene
                region = erosion[0]
                gene_length = len(self.glfo['seqs'][region][tmpline[region + '_gene']])
                if self.args.chain != 'h' and region == 'd':  # light chains dummy d treatment
                    assert gene_length == 1 and tmpline['d_gene'] == glutils.dummy_d_genes[self.args.chain]
                    tmpline[erosion + '_del'] = 1 if '5p' in erosion else 0  # always erode the whole dummy d from the left
                else:
                    max_erosion = max(0, gene_length/2 - 2)  # now that, son, is a heuristic
                    if region in utils.conserved_codons[self.args.chain]:
                        codon_pos = self.glfo[utils.conserved_codons[self.args.chain][region] + '-positions'][tmpline[region + '_gene']]
                        if '3p' in erosion:
                            n_bases_to_codon = gene_length - codon_pos - 3
                        elif '5p' in erosion:
                            n_bases_to_codon = codon_pos
                        max_erosion = min(max_erosion, n_bases_to_codon)
                    tmpline[erosion + '_del'] = min(max_erosion, numpy.random.geometric(1. / utils.scratch_mean_erosion_lengths[erosion]) - 1)
            for bound in utils.boundaries:
                mean_length = utils.scratch_mean_insertion_lengths[self.args.chain][bound]
                length = 0 if mean_length == 0 else numpy.random.geometric(1. / mean_length) - 1
                probs = [self.insertion_content_probs[bound][n] for n in utils.nukes]
                tmpline[bound + '_insertion'] = ''.join(numpy.random.choice(utils.nukes, size=length, p=probs))

            # have to add the 'seqs' by hand so utils.add_implicit_info doesn't barf (this duplicates code later on in recombinator)
            gl_seqs = {r : self.glfo['seqs'][r][tmpline[r + '_gene']] for r in utils.regions}
            for erosion in utils.real_erosions:
                region = erosion[0]
                e_length = tmpline[erosion + '_del']
                if '5p' in erosion:
                    gl_seqs[region] = gl_seqs[region][e_length:]
                elif '3p' in erosion:
                    gl_seqs[region] = gl_seqs[region][:len(gl_seqs[region]) - e_length]
            tmpline['seqs'] = [gl_seqs['v'] + tmpline['vd_insertion'] + gl_seqs['d'] + tmpline['dj_insertion'] + gl_seqs['j'], ]
            utils.add_implicit_info(self.glfo, tmpline)
            assert len(tmpline['in_frames']) == 1

        while 'in_frames' not in tmpline or not tmpline['in_frames'][0]:
            try_scratch_erode_insert(tmpline)

        # convert insertions back to lengths
        for bound in utils.boundaries + utils.effective_boundaries:
            tmpline[bound + '_insertion'] = len(tmpline[bound + '_insertion'])

        return tmpline

    # ----------------------------------------------------------------------------------------
    def choose_vdj_combo(self, reco_event):
        """ Choose the set of rearrangement parameters """

        vdj_choice = None
        if self.args.rearrange_from_scratch:  # generate an event without using the parameter directory
            vdj_choice = self.freqtable_index(self.get_scratchline())
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
    def write_mute_freqs(self, gene, seq, reco_event, reco_seq_fname):
        """ Read position-by-position mute freqs from disk for <gene>, renormalize, then write to a file for bppseqgen. """
        mute_freqs = self.get_mute_freqs(gene)

        rates = []  # list with a relative mutation rate for each position in <seq>
        total = 0.0
        # assert len(mute_freqs) == len(seq)  # only equal length if no erosions NO oh right but mute_freqs only covers areas we could align to...
        left_erosion_length = dict(reco_event.erosions.items() + reco_event.effective_erosions.items())[utils.get_region(gene) + '_5p']
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
            # NOTE really not sure why this doesn't really [seems to require an "extra" column] work with csv.DictWriter, but it doesn't -- bppseqgen barfs (I think maybe it expects a different newline character? don't feel like working it out)
            headstr = 'state'
            if not self.args.mutate_from_scratch:
                headstr += '\trate'
            reco_seq_file.write(headstr + '\n')
            for inuke in range(len(seq)):
                linestr = seq[inuke]
                if not self.args.mutate_from_scratch:
                    linestr += '\t%f' % rates[inuke]
                reco_seq_file.write(linestr + '\n')

    # ----------------------------------------------------------------------------------------
    def prepare_bppseqgen(self, seq, chosen_tree, n_leaf_nodes, gene, reco_event, seed):
        """ write input files and get command line options necessary to run bppseqgen on <seq> (which is a part of the full query sequence) """
        if len(seq) == 0:
            return None

        # write the tree to a tmp file
        workdir = self.workdir + '/' + utils.get_region(gene)
        os.makedirs(workdir)
        treefname = workdir + '/tree.tre'
        reco_seq_fname = workdir + '/start-seq.txt'
        leaf_seq_fname = workdir + '/leaf-seqs.fa'
        if n_leaf_nodes == 1:  # add an extra leaf to one-leaf trees so bppseqgen doesn't barf (when we read the output, we ignore the second leaf)
            lreg = re.compile('t1:[0-9]\.[0-9][0-9]*')
            leafstr = lreg.findall(chosen_tree)
            assert len(leafstr) == 1
            leafstr = leafstr[0]
            chosen_tree = chosen_tree.replace(leafstr, '(' + leafstr + ',' + leafstr + '):0.0')
        with opener('w')(treefname) as treefile:
            treefile.write(chosen_tree)
        self.write_mute_freqs(gene, seq, reco_event, reco_seq_fname)

        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] += ':' + self.args.partis_dir + '/packages/bpp/lib'

        # build up the command line
        # docs: http://biopp.univ-montp2.fr/apidoc/bpp-phyl/html/classbpp_1_1GTR.html that page is too darn hard to google
        bpp_binary = self.args.partis_dir + '/packages/bpp/bin/bppseqgen'
        if not os.path.exists(bpp_binary):
            raise Exception('bpp not found in %s' % os.path.dirname(bpp_binary))

        command = bpp_binary  # NOTE should I use the "equilibrium frequencies" option?
        command += ' alphabet=DNA'
        command += ' --seed=' + str(seed)
        command += ' input.infos=' + reco_seq_fname  # input file (specifies initial "state" for each position, and possibly also the mutation rate at that position)
        command += ' input.infos.states=state'  # column name in input file BEWARE bio++ undocumented defaults (i.e. look in the source code)
        command += ' input.tree.file=' + treefname
        command += ' input.tree.format=Newick'
        command += ' output.sequence.file=' + leaf_seq_fname
        command += ' output.sequence.format=Fasta'
        if self.args.mutate_from_scratch:
            command += ' model=JC69'
            command += ' input.infos.rates=none'  # BEWARE bio++ undocumented defaults (i.e. look in the source code)
            if self.args.flat_mute_freq is not None:
                command += ' rate_distribution=Constant'
            else:
                command += ' rate_distribution=Gamma(n=4,alpha=' + self.mute_models[utils.get_region(gene)]['gamma']['alpha']+ ')'
        else:
            command += ' input.infos.rates=rate'  # column name in input file
            pvpairs = [p + '=' + v for p, v in self.mute_models[utils.get_region(gene)]['gtr'].items()]
            command += ' model=GTR(' + ','.join(pvpairs) + ')'

        return {'cmd_str' : command, 'outfname' : leaf_seq_fname, 'workdir' : workdir, 'other-files' : [reco_seq_fname, treefname], 'env' : env}

    # ----------------------------------------------------------------------------------------
    def read_bppseqgen_output(self, cmdfo, n_leaf_nodes):
        mutated_seqs = []
        for seq_record in SeqIO.parse(cmdfo['outfname'], "fasta"):  # get the leaf node sequences from the file that bppseqgen wrote
            mutated_seqs.append(str(seq_record.seq))
            if n_leaf_nodes == 1:  # skip the extra leaf we added earlier
                break
        assert n_leaf_nodes == len(mutated_seqs)
        # self.check_tree_simulation(leaf_seq_fname, chosen_tree)
        os.remove(cmdfo['outfname'])
        for otherfname in cmdfo['other-files']:
            os.remove(otherfname)
        os.rmdir(cmdfo['workdir'])
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
            pos = random.randint(10, reco_event.final_codon_positions['v'])
        elif self.args.indel_location == 'cdr3':  # inside cdr3
            pos = random.randint(reco_event.final_codon_positions['v'], reco_event.final_codon_positions['j'])
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
        treg = re.compile('t[0-9][0-9]*')
        n_leaf_nodes = len(treg.findall(chosen_tree))
        cmdfos = []
        for region in utils.regions:
            simstr = reco_event.eroded_seqs[region]
            if region == 'd':
                simstr = reco_event.insertions['vd'] + simstr + reco_event.insertions['dj']
            cmdfos.append(self.prepare_bppseqgen(simstr, scaled_trees[region], n_leaf_nodes, reco_event.genes[region], reco_event, seed=irandom))

        utils.run_cmds([cfo for cfo in cmdfos if cfo is not None], sleep=False)  # shenanigan is to handle zero-length regional seqs

        mseqs = {}
        for ireg in range(len(utils.regions)):
            if cmdfos[ireg] is None:
                mseqs[utils.regions[ireg]] = ['' for _ in range(n_leaf_nodes)]  # return an empty string for each leaf node
            else:
                mseqs[utils.regions[ireg]] = self.read_bppseqgen_output(cmdfos[ireg], n_leaf_nodes)

        assert len(reco_event.final_seqs) == 0
        for iseq in range(n_leaf_nodes):
            seq = mseqs['v'][iseq] + mseqs['d'][iseq] + mseqs['j'][iseq]
            seq = reco_event.revert_conserved_codons(seq)  # if mutation screwed up the conserved codons, just switch 'em back to what they were to start with
            reco_event.final_seqs.append(seq)  # set final sequnce in reco_event

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
        os.remove(leaf_seq_fname)
        chosen_tree = dendropy.Tree.get_from_string(chosen_tree_str, 'newick')
        inferred_tree = dendropy.Tree.get_from_string(inferred_tree_str, 'newick')
        if self.args.debug:
            print '        tree diff -- symmetric %d   euke %f   rf %f' % (chosen_tree.symmetric_difference(inferred_tree), chosen_tree.euclidean_distance(inferred_tree), chosen_tree.robinson_foulds_distance(inferred_tree))
