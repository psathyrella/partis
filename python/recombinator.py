import copy
import operator
import tempfile
import sys
import csv
import time
import json
import random
import treegenerator
import numpy
import os
import re
import subprocess

import paramutils
import utils
import glutils
import indelutils
from event import RecombinationEvent

dummy_name_so_bppseqgen_doesnt_break = 'xxx'  # bppseqgen ignores branch length before mrca, so we add a spurious leaf with this name and the same total depth as the rest of the tree, then remove it after getting bppseqgen's output

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """
    def __init__(self, args, glfo, seed, workdir, outfname):  # NOTE <gldir> is not in general the same as <args.initial_germline_dir> # rm workdir
        self.args = args
        self.glfo = glfo

        # NOTE in general *not* the same as <self.args.workdir> and <self.args.outfname>
        self.workdir = tempfile.mkdtemp()
        self.outfname = outfname
        utils.prep_dir(self.workdir)

        # set <self.parameter_dir> (note that this is in general *not* the same as self.args.parameter_dir)
        if self.args.rearrange_from_scratch:  # currently not allowed to mutate from scratch without also rearranging from scratch (enforced in bin/partis)
            if self.args.mutate_from_scratch:
                self.parameter_dir = None
            else:
                self.parameter_dir = self.args.scratch_mute_freq_dir  # if you make up mute freqs from scratch, unless you're really careful you tend to get nonsense results for a lot of things (e.g. allele finding). So it's easier to copy over a reasonable set of mut freq parameters from somewhere.
        else:
            self.parameter_dir = self.args.parameter_dir + '/' + self.args.parameter_type

        self.index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
        self.mute_models = {}
        # self.treeinfo = []  # list of newick-formatted tree strings with region-specific branch info tacked at the end
        for region in utils.regions:
            self.mute_models[region] = {}
            for model in ['gtr', 'gamma']:
                self.mute_models[region][model] = {}

        self.allele_prevalence_freqs = glutils.read_allele_prevalence_freqs(args.allele_prevalence_fname) if args.allele_prevalence_fname is not None else {}
        self.version_freq_table = self.read_vdj_version_freqs()  # list of the probabilities with which each VDJ combo (plus other rearrangement parameters) appears in data (none if rearranging from scratch)
        self.insertion_content_probs = self.read_insertion_content()  # dummy/uniform if rearranging from scratch
        self.all_mute_freqs = {}

        # read shm info NOTE I'm not inferring the gtr parameters a.t.m., so I'm just (very wrongly) using the same ones for all individuals
        with open(self.args.gtrfname, 'r') as gtrfile:  # read gtr parameters
            reader = csv.DictReader(gtrfile)
            for line in reader:
                parameters = line['parameter'].split('.')
                region = parameters[0][3].lower()
                assert region == 'v' or region == 'd' or region == 'j'
                model = parameters[1].lower()
                parameter_name = parameters[2]
                assert model in self.mute_models[region]
                self.mute_models[region][model][parameter_name] = line['value']
        treegen = treegenerator.TreeGenerator(args, self.parameter_dir, seed=seed)
        self.treefname = self.workdir + '/trees.tre'
        treegen.generate_trees(seed, self.treefname)  # NOTE not really a newick file, since I hack on the per-region branch length info at the end of each line
        with open(self.treefname, 'r') as treefile:  # read in the trees (and other info) that we just generated
            self.treeinfo = treefile.readlines()
        os.remove(self.treefname)

        if os.path.exists(self.outfname):
            os.remove(self.outfname)
        elif not os.path.exists(os.path.dirname(os.path.abspath(self.outfname))):
            os.makedirs(os.path.dirname(os.path.abspath(self.outfname)))

        self.validation_values = {'heights' : {t : {'in' : [], 'out' : []} for t in ['all'] + utils.regions}}

    # ----------------------------------------------------------------------------------------
    def __del__(self):
        os.rmdir(self.workdir)

    # ----------------------------------------------------------------------------------------
    def read_insertion_content(self):
        if self.args.rearrange_from_scratch:
            return {b : {n : 1./len(utils.nukes) for n in utils.nukes} for b in utils.boundaries}

        insertion_content_probs = {}
        for bound in utils.boundaries:
            insertion_content_probs[bound] = {}
            with open(self.parameter_dir + '/' + bound + '_insertion_content.csv', 'r') as icfile:
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
    def get_mute_freqs(self, gene):
        if gene not in self.all_mute_freqs:
            self.read_mute_freq_stuff(gene)
        return self.all_mute_freqs[gene]

    # ----------------------------------------------------------------------------------------
    def read_mute_freq_stuff(self, gene):
        assert gene[:2] not in utils.boundaries  # make sure <gene> isn't actually an insertion (we used to pass insertions in here separately, but now they're smooshed onto either end of d)
        if self.args.mutate_from_scratch:
            self.all_mute_freqs[gene] = {'overall_mean' : self.args.default_scratch_mute_freq if self.args.flat_mute_freq is None else self.args.flat_mute_freq}
        else:
            gene_counts = utils.read_overall_gene_probs(self.parameter_dir, only_gene=gene, normalize=False, expect_zero_counts=True)
            approved_genes = [gene]
            if gene_counts < self.args.min_observations_to_write:  # if we didn't see it enough, average over all the genes that find_replacement_genes() gives us NOTE if <gene> isn't in the dict, it's because it's <args.datadir> but not in the parameter dir UPDATE not using datadir like this any more, so previous statement may not be true
                approved_genes += utils.find_replacement_genes(self.parameter_dir, min_counts=self.args.min_observations_to_write, gene_name=gene)
            self.all_mute_freqs[gene] = paramutils.read_mute_freqs_with_weights(self.parameter_dir, approved_genes)

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
        reco_event.set_post_erosion_codon_positions()

        # set the original conserved codon words, so we can revert them if they get mutated NOTE we do it here, *after* setting the full recombined sequence, so the germline Vs that don't extend through the cysteine don't screw us over
        reco_event.unmutated_codons = {}
        for region, codon in utils.conserved_codons[self.args.locus].items():
            fpos = reco_event.post_erosion_codon_positions[region]
            original_codon = reco_event.recombined_seq[fpos : fpos + 3]
            reco_event.unmutated_codons[region] = reco_event.recombined_seq[fpos : fpos + 3]
            # print fpos, original_codon, utils.codon_unmutated(codon, reco_event.recombined_seq, fpos)

        codons_ok = utils.both_codons_unmutated(self.glfo['locus'], reco_event.recombined_seq, reco_event.post_erosion_codon_positions, extra_str='      ', debug=self.args.debug)
        if not codons_ok:
            if self.args.rearrange_from_scratch and self.args.generate_germline_set:  # if you let it try more than once, it screws up the desired allele prevalence ratios
                raise Exception('arg')
            return False
        in_frame = utils.in_frame(reco_event.recombined_seq, reco_event.post_erosion_codon_positions, '', reco_event.effective_erosions['v_5p'])  # NOTE empty string is the fv insertion, which is hard coded to zero in event.py. I no longer recall the details of that decision, but I have a large amount of confidence that it's more sensible than it looks
        if self.args.rearrange_from_scratch and not in_frame:
            raise Exception('arg 2')  # if you let it try more than once, it screws up the desired allele prevalence ratios
            return False

        self.add_mutants(reco_event, irandom)

        reco_event.write_event(self.outfname, reco_event.line)

        return True

    # ----------------------------------------------------------------------------------------
    def freqtable_index(self, line):
        return tuple(line[column] for column in utils.index_columns)

    # ----------------------------------------------------------------------------------------
    def read_vdj_version_freqs(self):
        """ Read the frequencies at which various VDJ combinations appeared in data """
        if self.args.rearrange_from_scratch:
            return None

        version_freq_table = {}
        with open(self.parameter_dir + '/' + utils.get_parameter_fname('all', 'r')) as infile:
            in_data = csv.DictReader(infile)
            total = 0.0
            for line in in_data:  # NOTE do *not* assume the file is sorted
                skip = False
                for region in utils.regions:
                    if line[region + '_gene'] not in self.glfo['seqs'][region]:
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
    def try_scratch_erode_insert(self, tmpline, debug=False):
        utils.remove_all_implicit_info(tmpline)
        for erosion in utils.real_erosions:  # includes various contortions to avoid eroding the entire gene
            region = erosion[0]
            gene_length = len(self.glfo['seqs'][region][tmpline[region + '_gene']])
            if region == 'd' and not utils.has_d_gene(self.args.locus):  # dummy d genes: always erode the whole thing from the left
                assert gene_length == 1 and tmpline['d_gene'] == glutils.dummy_d_genes[self.args.locus]
                tmpline[erosion + '_del'] = 1 if '5p' in erosion else 0
            else:
                max_erosion = max(0, gene_length/2 - 2)  # heuristic
                if region in utils.conserved_codons[self.args.locus]:  # make sure not to erode a conserved codon
                    codon_pos = utils.cdn_pos(self.glfo, region, tmpline[region + '_gene'])
                    if '3p' in erosion:
                        n_bases_to_codon = gene_length - codon_pos - 3
                    elif '5p' in erosion:
                        n_bases_to_codon = codon_pos
                    max_erosion = min(max_erosion, n_bases_to_codon)
                tmpline[erosion + '_del'] = min(max_erosion, numpy.random.geometric(1. / utils.scratch_mean_erosion_lengths[erosion]) - 1)
        for bound in utils.boundaries:
            mean_length = utils.scratch_mean_insertion_lengths[self.args.locus][bound]
            length = 0 if mean_length == 0 else numpy.random.geometric(1. / mean_length) - 1
            probs = [self.insertion_content_probs[bound][n] for n in utils.nukes]
            tmpline[bound + '_insertion'] = ''.join(numpy.random.choice(utils.nukes, size=length, p=probs))

        if debug:
            print '    erosions:  %s' % ('   '.join([('%s %d' % (e, tmpline[e + '_del'])) for e in utils.real_erosions]))
            print '    insertions:  %s' % ('   '.join([('%s %s' % (b, tmpline[b + '_insertion'])) for b in utils.boundaries]))

        # have to add some things by hand so utils.add_implicit_info() doesn't barf (this duplicates code later on in recombinator)
        gl_seqs = {r : self.glfo['seqs'][r][tmpline[r + '_gene']] for r in utils.regions}
        for erosion in utils.real_erosions:
            region = erosion[0]
            e_length = tmpline[erosion + '_del']
            if '5p' in erosion:
                gl_seqs[region] = gl_seqs[region][e_length:]
            elif '3p' in erosion:
                gl_seqs[region] = gl_seqs[region][:len(gl_seqs[region]) - e_length]
        tmpline['seqs'] = [gl_seqs['v'] + tmpline['vd_insertion'] + gl_seqs['d'] + tmpline['dj_insertion'] + gl_seqs['j'], ]
        tmpline['unique_ids'] = [None]  # this is kind of hackey, but some things in the implicit info adder use it to get the number of sequences
        tmpline['input_seqs'] = copy.deepcopy(tmpline['seqs'])
        tmpline['indelfos'] = [indelutils.get_empty_indel(), ]
        utils.add_implicit_info(self.glfo, tmpline)
        assert len(tmpline['in_frames']) == 1

    # ----------------------------------------------------------------------------------------
    def get_scratchline(self):
        tmpline = {}

        # first choose the things that we'll only need to try choosing once (genes and effective (non-physical) deletions/insertions)
        for region in utils.regions:
            if len(self.glfo['seqs'][region]) == 0:
                raise Exception('no genes to choose from for %s' % region)
            probs = None  # it would make more sense to only do this prob calculation once, rather than for each event
            if region in self.allele_prevalence_freqs and len(self.allele_prevalence_freqs[region]) > 0:  # should really change it so it has to be the one or the other
                probs = [self.allele_prevalence_freqs[region][g] for g in self.glfo['seqs'][region].keys()]
            tmpline[region + '_gene'] = numpy.random.choice(self.glfo['seqs'][region].keys(), p=probs)  # order is arbitrary, but guaranteed to be the same as the previous line (https://docs.python.org/2/library/stdtypes.html#dict.items)
        for effrode in utils.effective_erosions:
            tmpline[effrode + '_del'] = 0
        for effbound in utils.effective_boundaries:
            tmpline[effbound + '_insertion'] = ''

        # then choose the things that we may need to try a few times (physical deletions/insertions)
        n_tries = 0
        while 'in_frames' not in tmpline or not tmpline['in_frames'][0]:
            self.try_scratch_erode_insert(tmpline)
            n_tries +=1
            if n_tries > 500:
                print tmpline['v_gene']
                print '%s tried to get an in-frame rearrangement %d times already' % (utils.color('yellow', 'warning'), n_tries)

        # convert insertions back to lengths (hoo boy this shouldn't need to be done)
        for bound in utils.all_boundaries:
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
    def write_mute_freqs(self, gene, seq, reco_event, reco_seq_fname):  # TODO unsurprisingly, this function profiles out to be kind of a dumb way to do it, in terms of run time
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
        with open(reco_seq_fname, 'w') as reco_seq_file:
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
        # add dummy leaf that we'll subsequently ignore (such are the vagaries of bppseqgen)
        chosen_tree = '(%s,%s:%.15f):0.0;' % (chosen_tree.rstrip(';'), dummy_name_so_bppseqgen_doesnt_break, treegenerator.get_mean_height(chosen_tree))
        with open(treefname, 'w') as treefile:
            treefile.write(chosen_tree)
        self.write_mute_freqs(gene, seq, reco_event, reco_seq_fname)

        env = os.environ.copy()
        env['LD_LIBRARY_PATH'] = env.get('LD_LIBRARY_PATH', '') + ':' + self.args.partis_dir + '/packages/bpp/lib'

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
        mutated_seqs = {}
        for seqfo in utils.read_fastx(cmdfo['outfname']):  # get the leaf node sequences from the file that bppseqgen wrote
            if seqfo['name'] == dummy_name_so_bppseqgen_doesnt_break:  # in the unlikely (impossible unless we change tree generators and don't tell them to use the same leaf names) event that we get a non-dummy leaf with this name, it'll fail at the assertion just below
                continue
            mutated_seqs[seqfo['name'].strip('\'')] = seqfo['seq']
        try:
            mutated_seqs = [mutated_seqs['t' + str(iseq + 1)] for iseq in range(len(mutated_seqs))]
        except KeyError as ke:
            raise Exception('leaf name %s not as expected in bppseqgen output %s' % (ke, cmdfo['outfname']))
        assert n_leaf_nodes == len(mutated_seqs)
        os.remove(cmdfo['outfname'])
        for otherfname in cmdfo['other-files']:
            os.remove(otherfname)
        os.rmdir(cmdfo['workdir'])
        return mutated_seqs

    # ----------------------------------------------------------------------------------------
    def add_shm_indels(self, reco_event):
        # NOTE that it will eventually make sense to add shared indel mutation according to the chosen tree -- i.e., probably, with some probability apply an indel instead of a point mutation
        if self.args.debug and self.args.indel_frequency > 0.:
            print '      indels'
        reco_event.indelfos = [indelutils.get_empty_indel() for _ in range(len(reco_event.final_seqs))]
        for iseq in range(len(reco_event.final_seqs)):
            if self.args.indel_frequency == 0.:  # no indels at all
                continue
            if numpy.random.uniform(0, 1) > self.args.indel_frequency:  # no indels for this sequence
                if self.args.debug:
                    print '        0'
                continue
            n_indels = numpy.random.choice(self.args.n_indels_per_indeld_seq)
            input_seq, indelfo = indelutils.add_indels(n_indels, reco_event.final_seqs[iseq], reco_event.recombined_seq,  # NOTE modifies <indelfo> and <codon_positions>
                                                       self.args.mean_indel_length, reco_event.final_codon_positions[iseq], indel_location=self.args.indel_location, dbg_pad=8, debug=self.args.debug)
            reco_event.final_seqs[iseq] = input_seq
            indelfo['v_gene'] = reco_event.genes['v']
            indelfo['j_gene'] = reco_event.genes['j']
            reco_event.indelfos[iseq] = indelfo
            # TODO need to update other things, at least input_seqs?

    # ----------------------------------------------------------------------------------------
    def add_mutants(self, reco_event, irandom):
        if self.args.mutation_multiplier is not None and self.args.mutation_multiplier == 0.:  # some of the stuff below fails if mut mult is actually 0.
            reco_event.final_seqs.append(reco_event.recombined_seq)  # set final sequnce in reco_event
            reco_event.indelfos = [indelutils.get_empty_indel() for _ in range(len(reco_event.final_seqs))]
            return

        # When generating trees, each tree's number of leaves and total depth are chosen from the specified distributions (a.t.m., by default n-leaves is from a geometric/zipf, and depth is from data)
        # This chosen depth corresponds to the sequence-wide mutation frequency.
        # In order to account for varying mutation rates in v, d, and j we simulate these regions separately, by appropriately rescaling the tree for each region.
        # i.e.: here we get the sequence-wide mute freq from the tree, and rescale it by the repertoire-wide ratios from data (which are stored in the tree file).
        # looks like e.g.: (t2:0.003751736951,t1:0.003751736951):0.001248262937;v:0.98,d:1.8,j:0.87, where the newick trees has branch lengths corresponding to the whole sequence  (i.e. the weighted mean of v, d, and j)
        # NOTE a.t.m (and probably permanently) the mean branch lengths for each region are the same for all the trees in the file, I just don't have a better place to put them while I'm passing from TreeGenerator to here than at the end of each line in the file
        treefostr = self.treeinfo[random.randint(0, len(self.treeinfo)-1)]  # per-region mutation info is tacked on after the tree... sigh. kind of hackey but works ok.
        assert treefostr.count(';') == 1
        isplit = treefostr.find(';') + 1
        chosen_tree = treefostr[:isplit]  # includes semi-colon
        mutefo = [rstr for rstr in treefostr[isplit:].split(',')]
        mean_total_height = treegenerator.get_mean_height(chosen_tree)
        regional_heights = {}  # per-region height, including <self.args.mutation_multiplier>
        for tmpstr in mutefo:
            region, ratio = tmpstr.split(':')
            assert region in utils.regions
            ratio = float(ratio)
            if self.args.mutation_multiplier is not None:  # multiply the branch lengths by some factor
                ratio *= self.args.mutation_multiplier
            regional_heights[region] = mean_total_height * ratio

        scaled_trees = {r : treegenerator.rescale_tree(chosen_tree, regional_heights[r]) for r in utils.regions}

        if self.args.debug:
            print '  chose tree with total height %f' % treegenerator.get_mean_height(chosen_tree)
            print '    regional trees rescaled to heights:  %s' % ('   '.join(['%s %.3f  (expected %.3f)' % (region, treegenerator.get_mean_height(scaled_trees[region]), regional_heights[region]) for region in utils.regions]))
            print treegenerator.get_ascii_tree(chosen_tree, extra_str='    ')

        n_leaves = treegenerator.get_n_leaves(chosen_tree)
        cmdfos = []
        for region in utils.regions:
            simstr = reco_event.eroded_seqs[region]
            if region == 'd':
                simstr = reco_event.insertions['vd'] + simstr + reco_event.insertions['dj']
            cmdfos.append(self.prepare_bppseqgen(simstr, scaled_trees[region], n_leaves, reco_event.genes[region], reco_event, seed=irandom))

        utils.run_cmds([cfo for cfo in cmdfos if cfo is not None], sleep=False)  # shenanigan is to handle zero-length regional seqs

        mseqs = {}
        for ireg in range(len(utils.regions)):  # NOTE kind of sketchy just using index in <utils.regions> (although it just depends on the loop immediately above a.t.m.)
            if cmdfos[ireg] is None:
                mseqs[utils.regions[ireg]] = ['' for _ in range(n_leaves)]  # return an empty string for each leaf node
            else:
                mseqs[utils.regions[ireg]] = self.read_bppseqgen_output(cmdfos[ireg], n_leaves)

        assert len(reco_event.final_seqs) == 0

        for iseq in range(n_leaves):
            seq = mseqs['v'][iseq] + mseqs['d'][iseq] + mseqs['j'][iseq]
            seq = reco_event.revert_conserved_codons(seq, debug=self.args.debug)  # if mutation screwed up the conserved codons, just switch 'em back to what they were to start with
            reco_event.final_seqs.append(seq)  # set final sequnce in reco_event
            reco_event.final_codon_positions.append(copy.deepcopy(reco_event.post_erosion_codon_positions))  # separate codon positions for each sequence, because of shm indels

        self.add_shm_indels(reco_event)
        reco_event.setline(irandom)  # set the line here because we use it when checking tree simulation, and want to make sure the uids are always set at the same point in the workflow

        self.check_tree_simulation(mean_total_height, regional_heights, scaled_trees, mseqs, reco_event)

        if self.args.debug:
            utils.print_reco_event(reco_event.line, extra_str='    ')

    # ----------------------------------------------------------------------------------------
    def infer_tree_from_leaves(self, region, in_tree, leafseqs):
        with tempfile.NamedTemporaryFile() as tmpfile:
            for iseq in range(len(leafseqs)):
                tmpfile.write('>t%s\n%s\n' % (iseq+1, leafseqs[iseq]))  # NOTE the order of the leaves/names is checked when reading bppseqgen output
            tmpfile.flush()  # BEWARE if you forget this you are fucked
            with open(os.devnull, 'w') as fnull:
                out_tree = subprocess.check_output('./bin/FastTree -gtr -nt ' + tmpfile.name, shell=True, stderr=fnull)

        in_height = treegenerator.get_mean_height(in_tree)
        out_height = treegenerator.get_mean_height(out_tree)
        base_width = 100
        print '  %s trees: input/output' % region
        print treegenerator.get_ascii_tree(in_tree, extra_str='        ', width=base_width)
        print treegenerator.get_ascii_tree(out_tree, extra_str='        ', width=int(base_width*out_height/in_height))

        if 'dendropy' not in sys.modules:
            import dendropy
        in_dtree = sys.modules['dendropy'].Tree.get_from_string(in_tree, 'newick')
        out_dtree = sys.modules['dendropy'].Tree.get_from_string(out_tree, 'newick')
        if self.args.debug:
            print '                   heights: %.3f   %.3f' % (in_height, out_height)
            print '      symmetric difference: %d' % in_dtree.symmetric_difference(out_dtree)
            print '        euclidean distance: %f' % in_dtree.euclidean_distance(out_dtree)
            print '              r-f distance: %f' % in_dtree.robinson_foulds_distance(out_dtree)

    # ----------------------------------------------------------------------------------------
    def check_tree_simulation(self, mean_total_height, regional_heights, scaled_trees, mseqs, reco_event, debug=False):
        assert reco_event.line is not None  # make sure we already set it
        mean_observed = {n : 0.0 for n in ['all'] + utils.regions}
        for iseq in range(len(reco_event.final_seqs)):
            if debug:
                print '   %4d  %.3f  %.3f' % (iseq, reco_event.line['mut_freqs'][iseq], mean_total_height)
            mean_observed['all'] += reco_event.line['mut_freqs'][iseq]
            for region in utils.regions:  # NOTE for simulating, we mash the insertions in with the D, but this isn't accounted for here
                rrate = utils.get_mutation_rate(reco_event.line, iseq=iseq, restrict_to_region=region)
                if debug:
                    print '       %.3f  %.3f' % (rrate, regional_heights[region])
                mean_observed[region] += rrate
        if debug:
            print '             in          out'
        for rname in ['all'] + utils.regions:
            mean_observed[rname] /= float(len(reco_event.final_seqs))
            if rname == 'all':
                input_height = mean_total_height
                if self.args.mutation_multiplier is not None:  # multiply the branch lengths by some factor
                    input_height *= self.args.mutation_multiplier
            else:
                input_height = regional_heights[rname]
            self.validation_values['heights'][rname]['in'].append(input_height)
            self.validation_values['heights'][rname]['out'].append(mean_observed[rname])
            if debug:
                print '  %4s    %7.3f     %7.3f' % (rname, input_height, mean_observed[rname])

        # need to root them
        # for region in utils.regions:
        #     self.infer_tree_from_leaves(region, scaled_trees[region], mseqs[region])  # NOTE haven't yet reverted codons

    # ----------------------------------------------------------------------------------------
    def print_validation_values(self):
        # NOTE the v, d, and all are systematically low, while j is high. Don't feel like continuing to figure out all the contributors a.t.m. though
        print '  tree heights:'
        print '        in      out          diff'
        for vtype in ['all'] + utils.regions:
            vvals = self.validation_values['heights'][vtype]
            deltas = [(vvals['out'][i] - vvals['in'][i]) for i in range(len(vvals['in']))]
            print '      %.3f   %.3f    %+.3f +/- %.3f      %s' % (numpy.mean(vvals['in']), numpy.mean(vvals['out']),
                                                                   numpy.mean(deltas), numpy.std(deltas) / len(deltas), vtype)  # NOTE each delta is already the mean of <n_leaves> independent measurements
