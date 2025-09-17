from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import copy
import operator
import tempfile
import sys
import csv
import time
import json
import random
import numpy
import os
import re
import math
import subprocess
import scipy.stats
from collections import defaultdict

from . import paramutils
from . import utils
from . import glutils
from . import treeutils
from . import indelutils
from . import treegenerator
from .event import RecombinationEvent
from io import open

dummy_name_so_bppseqgen_doesnt_break = 'xxx'  # bppseqgen ignores branch length before mrca, so we add a spurious leaf with this name and the same total depth as the rest of the tree, then remove it after getting bppseqgen's output

# ----------------------------------------------------------------------------------------
class SimuGiveUpError(Exception):
    pass

#----------------------------------------------------------------------------------------
class Recombinator(object):
    """ Simulates the process of VDJ recombination """
    def __init__(self, args, glfo, seed, workdir, heavy_chain_events=None):  # NOTE <gldir> is not in general the same as <args.initial_germline_dir> # rm workdir
        self.args = args
        self.glfo = glfo
        if len(glfo['seqs']['v']) > 100:  # this is kind of a shitty criterion, but I don't know what would be better (we basically just want to warn people if they're simulating from data/germlines/human)
            print('  note: simulating with a very large number (%d) of V genes (the use of realistic diploid sets can be controlled either by using inferred germline sets that you\'ve got lying around (--reco-parameter-dir), or with --generate-germline-set)' % len(glfo['seqs']['v']))

        self.workdir = tempfile.mkdtemp()
        utils.prep_dir(self.workdir)

        assert self.args.parameter_dir is None  # we don't use the overall dir anywhere here, only reco_ and shm_
        self.reco_parameter_dir = utils.parameter_type_subdir(self.args, self.args.reco_parameter_dir) if self.args.reco_parameter_dir is not None else None  # only used if not rearranging from scratch
        self.shm_parameter_dir = utils.parameter_type_subdir(self.args, self.args.shm_parameter_dir) if self.args.shm_parameter_dir is not None else None  # only used if not mutating from scratch
        print('    reco params: %s     shm params: %s' % (utils.non_none([self.reco_parameter_dir, 'scratch']), utils.non_none([self.shm_parameter_dir, 'scratch'])))

        self.index_keys = {}  # this is kind of hackey, but I suspect indexing my huge table of freqs with a tuple is better than a dict
        self.mute_models = {}
        # self.treeinfo = []  # list of newick-formatted tree strings plus region-specific branch length info (initialized below)
        for region in utils.regions:
            self.mute_models[region] = {}
            for model in ['gtr', 'gamma']:
                self.mute_models[region][model] = {}

        self.allele_prevalence_freqs = None if args.allele_prevalence_fname is None else glutils.read_allele_prevalence_freqs(args.allele_prevalence_fname)
        self.version_freq_table = self.read_vdj_version_freqs()  # list of the probabilities with which each VDJ combo (plus other rearrangement parameters) appears in data (none if rearranging from scratch)
        self.insertion_content_probs = self.read_insertion_content()  # dummy/uniform if rearranging from scratch
        self.all_mute_freqs = {}  # NOTE see description of the difference in hmmwriter.py
        self.all_mute_counts = {}

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
        treegen = treegenerator.TreeGenerator(args, self.shm_parameter_dir)
        self.treefname = self.workdir + '/trees.yaml'  # yaml file contains a list of tree strings (like a multi-line newick file), as well as some mutation info
        treegen.generate_trees(seed, self.treefname, self.workdir)
        with open(self.treefname, 'r') as treefile:  # read in the trees (and other info) that we just generated
            self.treeinfo = json.load(treefile)
        os.remove(self.treefname)

        # because when writing bppseqgen input we carefully adjust the per-base rates based on the germline base (so that the prob to go from germline to germline is zero, i.e. we're not at equilibrium) bppseqgen makes too many mutations (since it assumes we're at equilibrium/assumes the rate matrix and initial state are uncorrelated)
        # we correct for this by multiplying by this heuristic factor, although andy has some corrections that we could calculate, see data/recombinator/andy-bpp-per-base-correction.txt
        self.per_base_mutation_multiplier = 0.6

        self.validation_values = {'heights' : {t : {'in' : [], 'out' : []} for t in ['all'] + utils.regions}}
        self.validation_values['bpp-times'] = []

        self.heavy_chain_events = heavy_chain_events

        self.n_max_tries = 10000  # applies in two different (recursive) places, which is a little weird, but whatever

    # ----------------------------------------------------------------------------------------
    def __del__(self):
        if len(os.listdir(self.workdir)) == 0:
            os.rmdir(self.workdir)
        else:
            print('  couldn\'t exit cleanly, workdir %s not empty' % self.workdir)

    # ----------------------------------------------------------------------------------------
    def read_insertion_content(self):
        # ----------------------------------------------------------------------------------------
        def default_content():
            return {n : 1./len(utils.nukes) for n in utils.nukes}

        # ----------------------------------------------------------------------------------------
        if self.args.rearrange_from_scratch:
            return {b : default_content() for b in utils.boundaries}

        insertion_content_probs = {}
        for bound in utils.boundaries:
            insertion_content_probs[bound] = {}
            with open(self.reco_parameter_dir + '/' + bound + '_insertion_content.csv', 'r') as icfile:
                reader = csv.DictReader(icfile)
                total = 0
                for line in reader:
                    insertion_content_probs[bound][line[bound + '_insertion_content']] = int(line['count'])
                    total += int(line['count'])
                if total > 0:
                    for nuke in utils.nukes:
                        if nuke not in insertion_content_probs[bound]:
                            print('    %s not in insertion content probs, adding with zero' % nuke)
                            insertion_content_probs[bound][nuke] = 0
                        insertion_content_probs[bound][nuke] /= float(total)
                else:  # i think this will only happen for light chain (i.e. when one of the bounds has all-zero counts)
                    insertion_content_probs[bound] = default_content()

            assert utils.is_normed(insertion_content_probs[bound])

        return insertion_content_probs


    # ----------------------------------------------------------------------------------------
    def get_mute_freqs(self, gene):
        if gene not in self.all_mute_freqs:
            self.read_mute_freq_stuff(gene)
        return self.all_mute_freqs[gene]

    # ----------------------------------------------------------------------------------------
    def get_mute_counts(self, gene):
        if gene not in self.all_mute_counts:
            self.read_mute_freq_stuff(gene)
        return self.all_mute_counts[gene]

    # ----------------------------------------------------------------------------------------
    def read_mute_freq_stuff(self, gene):
        assert gene[:2] not in utils.boundaries  # make sure <gene> isn't actually an insertion (we used to pass insertions in here separately, but now they're smooshed onto either end of d)
        if self.args.mutate_from_scratch:
            self.all_mute_freqs[gene] = {'overall_mean' : self.args.scratch_mute_freq}
            if not self.args.no_per_base_mutation:
                raise Exception('can\'t yet mutate from scratch with per-base mutation')
            # self.all_mute_counts[gene] = {'overall_mean' : self.args.scratch_mute_freq} # TODO see TODOs further down, but at the moment we don't use these if --mutate-from-scratch is set
        else:
            extra_genes = []
            gene_counts = utils.read_overall_gene_probs(self.shm_parameter_dir, only_gene=gene, normalize=False, expect_zero_counts=True)
            if gene_counts < self.args.min_observations_per_gene:  # if we didn't see it enough, average over all the genes that find_replacement_genes() gives us NOTE if <gene> isn't in the dict, it's because it's in <args.datadir> but not in the parameter dir UPDATE not using datadir like this any more, so previous statement may not be true
                extra_genes = utils.find_replacement_genes(self.shm_parameter_dir, min_counts=self.args.min_observations_per_gene, gene_name=gene)

            self.all_mute_freqs[gene] = paramutils.read_mute_freqs_with_weights(self.shm_parameter_dir, list(set([gene] + extra_genes)))  # NOTE these fcns do quite different things as far as smoothing (see comments elsewhere)
            per_base_extra_genes = None
            if len(extra_genes) > 0 and gene_counts == 0:  # if we have even one count for the gene we want, we want to use it, since the per-base counts will be wrong for even other alleles at any positions at which they differe
                per_base_extra_genes = list(set(extra_genes) - set([gene]))
            self.all_mute_counts[gene] = paramutils.read_mute_counts(self.shm_parameter_dir, gene, utils.get_locus(gene), extra_genes=per_base_extra_genes)

    # ----------------------------------------------------------------------------------------
    def combine(self, initial_irandom, i_choose_tree=None, i_heavy_event=None):
        """ keep running self.try_to_combine() until you get a good event """
        line = None
        itry = 0
        while line is None:
            if itry > 0 and self.args.debug:
                print('    unproductive event -- rerunning (try %d)  ' % itry)  # probably a weirdly long v_3p or j_5p deletion
            line = self.try_to_combine(initial_irandom + itry, i_choose_tree=i_choose_tree, i_heavy_event=i_heavy_event)
            itry += 1
            if itry > self.n_max_tries:
                raise Exception('too many tries %d in recombinator' % itry)
        return line

    # ----------------------------------------------------------------------------------------
    def try_to_combine(self, irandom, i_choose_tree=None, i_heavy_event=None):
        """
        Create a recombination event and write it to disk
        <irandom> is used as the seed for the myriad random number calls.
        If combine() is called with the same <irandom>, it will find the same event, i.e. it should be a random number, not just a seed
        """
        if self.args.debug:
            print('combine (seed %d)' % irandom)
        numpy.random.seed(irandom)
        random.seed(irandom)

        reco_event = RecombinationEvent(self.glfo)
        try:
            self.choose_vdj_combo(reco_event, i_heavy_event=i_heavy_event)  # it's kind of ugly to raise an exception here, rather than having the fcn[s] return None as below, but the control flow gets complicated that way. Maybe would be better to switch all of them to exceptions
        except SimuGiveUpError as err:
            print(err)
            return None

        self.erode_and_insert(reco_event)

        # set the original conserved codon words, so we can revert them if they get mutated NOTE we do it here, *after* setting the full recombined sequence, so the germline Vs that don't extend through the cysteine don't screw us over (update: we should no longer ever encounter Vs that're screwed up like this)
        reco_event.unmutated_codons = {}
        for region, codon in utils.conserved_codons[self.args.locus].items():
            fpos = reco_event.post_erosion_codon_positions[region]
            original_codon = reco_event.recombined_seq[fpos : fpos + 3]
            reco_event.unmutated_codons[region] = reco_event.recombined_seq[fpos : fpos + 3]
            # print fpos, original_codon, utils.codon_unmutated(codon, reco_event.recombined_seq, fpos)

        codons_ok = utils.both_codons_unmutated(self.glfo['locus'], reco_event.recombined_seq, reco_event.post_erosion_codon_positions, extra_str='      ', debug=self.args.debug)
        if not codons_ok:
            if self.args.rearrange_from_scratch and self.args.generate_germline_set:
                raise Exception('mutated invariant codons, but since --rearrange-from-scratch and --generate-germline-set are set, we can\'t retry, since it would screw up the prevalence ratios')  # if you let it try more than once, it screws up the desired allele prevalence ratios
            return None
        in_frame = utils.in_frame(reco_event.recombined_seq, reco_event.post_erosion_codon_positions, '', reco_event.effective_erosions['v_5p'])  # NOTE empty string is the fv insertion, which is hard coded to zero in event.py. I no longer recall the details of that decision, but I have a large amount of confidence that it's more sensible than it looks
        if not self.args.allow_nonfunctional_scratch_seqs and self.args.rearrange_from_scratch and not in_frame:
            raise Exception('out of frame rearrangement, but since --rearrange-from-scratch is set we can\'t retry (it would screw up the prevalence ratios)')  # if you let it try more than once, it screws up the desired allele prevalence ratios
            # return None

        self.add_mutants(reco_event, irandom, i_choose_tree=i_choose_tree)

        line = reco_event.line
        # NOTE don't use reco_event after here, since we don't modify it when we remove non-functional sequences (as noted elsewhere, it would be nice to eventually update to just using <line>s instead of <reco_event> now that that's possible)
        if self.args.remove_nonfunctional_seqs:
            functional_iseqs = [iseq for iseq in range(len(line['unique_ids'])) if utils.is_functional(line, iseq)]
            if len(functional_iseqs) == 0:  # none are functional -- try again
                return None
            if len(functional_iseqs) < len(line['unique_ids']):  # it's generally very rare for them to all be functional
                if self.args.debug:
                    print('      removing %d nonfunctional seqs (of %d)' % (len(line['unique_ids']) - len(functional_iseqs), len(line['unique_ids'])))
                utils.restrict_to_iseqs(line, functional_iseqs, self.glfo)

        return line

    # ----------------------------------------------------------------------------------------
    def freqtable_index(self, line):
        return tuple(line[column] for column in utils.index_columns)

    # ----------------------------------------------------------------------------------------
    def read_vdj_version_freqs(self):
        """ Read the frequencies at which various rearrangement events (VDJ combinations + insertion/deletion lengths) appeared in data """
        if self.args.rearrange_from_scratch:
            return None

        version_freq_table = {}
        n_skipped_gene, n_skipped_cdr3, n_used = 0., 0., 0.
        cdr3_counts = {}
        with open(self.reco_parameter_dir + '/' + utils.get_parameter_fname('all', 'r')) as infile:
            in_data = csv.DictReader(infile)
            total = 0.0
            for line in in_data:  # NOTE do *not* assume the file is sorted
                if any(line[region + '_gene'] not in self.glfo['seqs'][region] for region in utils.regions):
                    n_skipped_gene += float(line['count'])
                    continue
                if self.args.allowed_cdr3_lengths is not None and int(line['cdr3_length']) not in self.args.allowed_cdr3_lengths:
                    n_skipped_cdr3 += float(line['count'])  #  NOTE this isn't really right if we're also skipping genes above, but whatever (well it's not really wrong, but it'd be different if we did it before the gene skipping rather than after)
                    continue
                total += float(line['count'])
                index = self.freqtable_index(line)
                assert index not in version_freq_table
                version_freq_table[index] = float(line['count'])
                n_used += float(line['count'])
                if self.args.allowed_cdr3_lengths is not None:  # maybe it's worth printing what the final cdr3 length breakdown is? (so it's obvious if e.g. you asked for 27, 30, 33 but only got 33)
                    if int(line['cdr3_length']) not in cdr3_counts:
                        cdr3_counts[int(line['cdr3_length'])] = 0
                    cdr3_counts[int(line['cdr3_length'])] += float(line['count'])
            if n_skipped_gene > 0:
                print('    skipped %.0f / %.0f (%.3f) vdj freq counts with genes not in glfo (used %.0f)' % (n_skipped_gene, n_skipped_gene + n_used, n_skipped_gene / (n_skipped_gene + n_used), n_used))
            if n_skipped_cdr3 > 0:
                print('    skipped %.0f / %.0f (%.3f) vdj freq counts with cdr3 lengths not among [%s] (used %.0f)' % (n_skipped_cdr3, n_skipped_cdr3 + n_used, n_skipped_cdr3 / (n_skipped_cdr3 + n_used), ' '.join(str(c) for c in self.args.allowed_cdr3_lengths), n_used))
                print('       final cdr3 lengths: %s' % '   '.join(('%d: %.0f'%(l, cdr3_counts[l])) for l in sorted(cdr3_counts)))

        if len(version_freq_table) == 0:
            raise Exception('didn\'t find any gene combinations in %s' % self.reco_parameter_dir + '/' + utils.get_parameter_fname('all', 'r'))

        # then normalize
        test_total = 0.0
        for index in version_freq_table:
            version_freq_table[index] /= total
            test_total += version_freq_table[index]
        assert utils.is_normed(test_total, this_eps=1e-8)
        assert len(version_freq_table) < 1e8  # if it gets *too* large, choose_vdj_combo() below isn't going to work because of numerical underflow. Note there's nothing special about 1e8, it's just that I'm pretty sure we're fine *up* to that point, and once we get beyond it we should think about doing things differently
        return version_freq_table

    # ----------------------------------------------------------------------------------------
    # keep track of allowed options, and if specified also restrict options according to the specified correlation value (corr of 1 means 1 option, of 0 means take all options)
    def handle_options_for_correlation(self, pthis, this_options, corr_vals, allowed_vals, parent_line, probs=None, mean_max=None):  # return <this_options> (allowed values for parameter <pthis>) that has been restricted according to any correlations specified in self.correlation_values that are for keys already in <parent_line> (<allowed_vals> contains the options for other ones)
        # NOTE this needs to be called even for parameters we don't want to restrict, in order to add <this_options> to allowed_vals['current']

        if mean_max is not None:  # need to init list of possible values + probs
            assert probs is None and this_options is None
            mean_len, max_len = mean_max
            this_options, probs = zip(*[(l - 1, scipy.stats.geom.pmf(l, 1. / mean_len)) for l in range(1, max_len + 1)])  # NOTE has to match numpy.random.geometric below NOTE also probs aren't normalized here, but that shouldn't matter NOTE also l-1 is to match similar subtractions in non-correlation code

        allowed_vals['current'][pthis] = this_options  # note that this means all the modifications to <this_options> below are also modifying the list in allowed_vals['current']
        if corr_vals is None:
            return this_options, None if probs is None else utils.normalize(probs)

        for (param_pair), corr_val in corr_vals.items():
            if pthis != param_pair[1]:
                continue
            pother = param_pair[0]
            n_before = len(this_options)
            iother = allowed_vals['parent'][pother].index(parent_line[pother])  # pother should (now) always already be in parent_line, since we're enforcing the ordering of the allowed correlation pairs
            istart = iother % len(this_options)
            n_to_take = max(1, int((1. - corr_val) * float(len(this_options))))
            i_taken = [i % len(this_options) for i in range(istart, istart + n_to_take)]  # start from <istart> and add <n_to_take>, wrapping around to index 0 if necessary
            if self.args.debug > 1:
                def cstr(i, cfcn): return utils.color('red' if cfcn(i) else None, str(i), width=3)
                pair_lstrs = ('', '') if self.args.paired_correlation_values is None else [' (%s)  '%utils.color('blue', l) for l in (utils.get_locus(parent_line['v_gene']), self.args.locus)]
                print('     %s%s %s%s correlation: %.2f, parent val %s iother %d, restricting to max(1, (1 - %.2f) * %d) = %d values starting from index %d %% %d = %d' % (param_pair[0], pair_lstrs[0], param_pair[1], pair_lstrs[1],
                                                                                                                                                            corr_val, parent_line[pother], iother, corr_val, len(this_options), n_to_take, iother, len(this_options), istart))
                print('      %10s: %s' % (param_pair[0], ' '.join(cstr(i, lambda x: x==iother) for i in range(len(allowed_vals['parent'][pother])))))
                print('      %10s: %s' % (param_pair[1], ' '.join(cstr(i, lambda x: x in i_taken) for i in range(len(this_options)))))
            this_options = [this_options[i] for i in i_taken]
            if probs is not None:
                assert len(probs) == n_before
                probs = [probs[i] for i in i_taken]
            if self.args.debug:
                print('        reducing options for %s (given previous choice of %s) from %d to %d' % (pthis, pother, n_before, len(this_options)))

        return this_options, None if probs is None else utils.normalize(probs)

    # ----------------------------------------------------------------------------------------
    def try_scratch_erode_insert(self, tmpline, corr_vals=None, allowed_vals=None, parent_line=None, debug=False):  # non-None corr_vals determines if we're applying correlations
        # ----------------------------------------------------------------------------------------
        def get_heavy_del(region, erosion, gene_length):
            if self.args.no_insertions_or_deletions:
                return 0
            max_erosion = max(0, gene_length//2 - 2)  # heuristic
            if region in utils.conserved_codons[self.args.locus]:  # make sure not to erode a conserved codon
                codon_pos = utils.cdn_pos(self.glfo, region, tmpline[region + '_gene'])
                if '3p' in erosion:
                    n_bases_to_codon = gene_length - codon_pos - 3
                elif '5p' in erosion:
                    n_bases_to_codon = codon_pos
                max_erosion = min(max_erosion, n_bases_to_codon)
            mean_len = utils.scratch_mean_erosion_lengths[self.args.locus][erosion]
            if corr_vals is None and allowed_vals is None:  # the case where they're different is heavy chain for paired correlation, when corr_vals is None so we don't apply any correlations, but allowed_vals is *not* None so we can keep track of allowed values
                e_len = min(max_erosion, len_fcn(1. / mean_len) - 1)
            else:
                lens, probs = self.handle_options_for_correlation(erosion+'_del', None, corr_vals, allowed_vals, parent_line, mean_max=(mean_len, max_erosion))
                e_len = numpy.random.choice(lens, p=probs)
            return e_len
        # ----------------------------------------------------------------------------------------
        len_fcn = numpy.random.geometric  # has to match scipy.stats.geom above NOTE also the -1 in several places
        utils.remove_all_implicit_info(tmpline)
        for erosion in utils.real_erosions:  # includes various contortions to avoid eroding the entire gene
            region = erosion[0]
            gene_length = len(self.glfo['seqs'][region][tmpline[region + '_gene']])
            if region == 'd' and not utils.has_d_gene(self.args.locus):  # dummy d genes: always erode the whole thing from the left
                assert gene_length == 1 and tmpline['d_gene'] == glutils.dummy_d_genes[self.args.locus]
                tmpline[erosion + '_del'] = 1 if '5p' in erosion else 0
            else:
                tmpline[erosion + '_del'] = get_heavy_del(region, erosion, gene_length)
        for bound in utils.boundaries:
            mean_len = utils.scratch_mean_insertion_lengths[self.args.locus][bound]
            if mean_len == 0 or self.args.no_insertions_or_deletions:
                i_len = 0  # mean_len if 0 means it *needs* to be zero, e.g. vd insertion for light chain
            else:
                if corr_vals is None and allowed_vals is None:
                    i_len = len_fcn(1. / mean_len) - 1
                else:
                    lens, probs = self.handle_options_for_correlation(bound+'_insertion', None, corr_vals, allowed_vals, parent_line, mean_max=(mean_len, 2*int(mean_len)))  # it's kind of weird to just limit it to twice the mean length here, but since handle_options_for_correlation() doesn't account for probs when choosing which to keep, if we make it bigger we choose those super large values too frequently
                    i_len = numpy.random.choice(lens, p=probs)  # maybe this should also always be 0 if mean_len is 0?
            cnt_probs = [self.insertion_content_probs[bound][n] for n in utils.nukes]
            tmpline[bound + '_insertion'] = ''.join(numpy.random.choice(utils.nukes, size=i_len, p=cnt_probs, replace=True))

        if debug:
            print('    erosions:  %s' % ('   '.join([('%s %d' % (e, tmpline[e + '_del'])) for e in utils.real_erosions])))
            print('    insertions:  %s' % ('   '.join([('%s %s' % (b, tmpline[b + '_insertion'])) for b in utils.boundaries])))

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
        tmpline['input_seqs'] = copy.deepcopy(tmpline['seqs'])  # NOTE has to be updated _immediately_ so seqs and input_seqs don't get out of sync
        tmpline['indelfos'] = [indelutils.get_empty_indel(), ]
        utils.add_implicit_info(self.glfo, tmpline)
        assert len(tmpline['in_frames']) == 1
        if self.args.paired_correlation_values is not None:
            if utils.has_d_gene(self.args.locus):  # if we're doing heavy/light correlations we need to write some extra info to the output file: either allowed vals (if this is heavy chain, so we can use them to determine light chain rearrangement parameters) or heavy chain uids (if this is light chain, so we can later make sure we pair the correct h/l events together)
                tmpline['heavy-chain-correlation-info'] = allowed_vals['current']
            else:
                tmpline['heavy-chain-correlation-info'] = None if parent_line is None else {'heavy-chain-uids' : parent_line['unique_ids']}

    # ----------------------------------------------------------------------------------------
    def get_scratchline(self, i_heavy_event=None):
        # ----------------------------------------------------------------------------------------
        def keep_trying(tmpline):
            if not self.args.allow_nonfunctional_scratch_seqs and not tmpline['in_frames'][0]:
                failcounters['out-of-frame'] += 1
                return True
            if not self.args.allow_nonfunctional_scratch_seqs and tmpline['stops'][0]:
                failcounters['stop'] += 1
                return True
            if self.args.allowed_cdr3_lengths is not None and tmpline['cdr3_length'] not in self.args.allowed_cdr3_lengths:
                failcounters['allowed-cdr3'] += 1
                return True
            return False
        # ----------------------------------------------------------------------------------------
        def fcstrs():
            return ['%s %d'%(f, c) for f, c in sorted(failcounters.items(), key=operator.itemgetter(1), reverse=True)]

        # ----------------------------------------------------------------------------------------
        tmpline = {}

        assert self.args.correlation_values is None or self.args.paired_correlation_values is None  # this is already checked elsewhere, but having it here makes it clearer how things work
        corr_vals, allowed_vals, parent_line = None, None, None  # allowed_vals is to keep track of all the options we had for each parameter, for when a later parameter wants to be correlated with it (note: could probably include allowed_genes and whatnot in allowed_vals now)
        if self.args.correlation_values is not None:
            corr_vals = self.args.correlation_values
            allowed_vals = {'current' : {}}
            allowed_vals['parent'] = allowed_vals['current']
            parent_line = tmpline
        elif self.args.paired_correlation_values is not None:
            allowed_vals = {'current' : {}, 'parent' : None}
            if not utils.has_d_gene(self.args.locus):  # we only *apply* paired correlations in light chain (whereas for heavy chain we need to keep track of original options but without applying any correlations)
                corr_vals = self.args.paired_correlation_values
                assert i_heavy_event is not None
                allowed_vals['parent'] = self.heavy_chain_events[i_heavy_event]['heavy-chain-correlation-info']
                parent_line = self.heavy_chain_events[i_heavy_event]
                if self.args.debug:
                    print('    taking parent values for correlation from heavy chain event with: %s  %s  %s  cdr3: %d  uids: %s' % (utils.color_gene(parent_line['v_gene']), utils.color_gene(parent_line['d_gene']), utils.color_gene(parent_line['j_gene']), parent_line['cdr3_length'], ' '.join(parent_line['unique_ids'])))

        # first choose the things that we'll only need to try choosing once (genes and effective (non-physical) deletions/insertions)
        for region in utils.regions:
            if len(self.glfo['seqs'][region]) == 0:
                raise Exception('no genes to choose from for %s (if light chain d, you need to explicitly add the dummy d to --only-genes)' % region)
            probs = None  # it would make more sense to only do this prob calculation once, rather than for each event
            allowed_genes = sorted(self.glfo['seqs'][region].keys())  # sorted is really just to make the correlation stuff feel more repeatable (.keys() should have repeatable order within a version: https://docs.python.org/2/library/stdtypes.html#dict.items)
            if self.allele_prevalence_freqs is not None:
                if region in self.allele_prevalence_freqs and len(self.allele_prevalence_freqs[region]) > 0:  # should really change it so it has to be the one or the other
                    probs = [self.allele_prevalence_freqs[region][g] for g in allowed_genes]
            if corr_vals is not None or allowed_vals is not None:
                allowed_genes, probs = self.handle_options_for_correlation(region+'_gene', allowed_genes, corr_vals, allowed_vals, parent_line, probs=probs)
            tmpline[region + '_gene'] = str(numpy.random.choice(allowed_genes, p=probs))
        for effrode in utils.effective_erosions:
            tmpline[effrode + '_del'] = 0
        for effbound in utils.effective_boundaries:
            tmpline[effbound + '_insertion'] = ''

        # then choose the things that we may need to try a few times (physical deletions/insertions)
        itry, failcounters = 0, defaultdict(int)
        while itry == 0 or keep_trying(tmpline):  # keep trying until it's both in frame and has no stop codons
            if self.args.debug and itry > 0:
                print('    %s: retrying scratch rearrangement (fail counters: %s)' % (utils.color('blue', 'itry %d'%itry), '  '.join(fcstrs())))
            self.try_scratch_erode_insert(tmpline, corr_vals=corr_vals, allowed_vals=allowed_vals, parent_line=parent_line)  # NOTE the content of these insertions doesn't get used. They're converted to lengths just below (we make up new ones in self.erode_and_insert())
            itry += 1
            if itry % 500 == 0:
                print('      %s finding an in-frame and stop-less %srearrangement is taking an oddly large number of tries (%d so far with failcounters: %s)' % ('note:', '' if self.args.allowed_cdr3_lengths is None else '(and with --allowed-cdr3-length) ', itry, '  '.join(fcstrs())))
            if itry > self.n_max_tries:
                raise SimuGiveUpError('    %s too many tries (%d > %d) when trying to get scratch line so giving up (well, probably retrying from further up, i.e. with new/iterated seed)' % (utils.wrnstr(), itry, self.n_max_tries))

        # convert insertions back to lengths (hoo boy this shouldn't need to be done)
        for bound in utils.all_boundaries:
            tmpline[bound + '_insertion'] = len(tmpline[bound + '_insertion'])

        return tmpline

    # ----------------------------------------------------------------------------------------
    def choose_vdj_combo(self, reco_event, i_heavy_event=None):  # NOTE similarity to hist.sample()
        """ Choose the set of rearrangement parameters """

        vdj_choice = None
        h_corr_line = None
        if self.args.rearrange_from_scratch:  # generate an event without using the parameter directory
            tmpline = self.get_scratchline(i_heavy_event=i_heavy_event)
            vdj_choice = self.freqtable_index(tmpline)
            if self.args.paired_correlation_values:
                h_corr_line = tmpline
        else:  # use real parameters from a directory
            iprob = numpy.random.uniform(0, 1)
            sum_prob = 0.0
            for tmpchoice in self.version_freq_table:  # assign each vdj choice a segment of the interval [0,1], and choose the one which contains <iprob>
                sum_prob += self.version_freq_table[tmpchoice]
                if iprob < sum_prob:
                    vdj_choice = tmpchoice
                    break

            assert vdj_choice is not None  # shouldn't fall through to here

        reco_event.set_vdj_combo(vdj_choice, self.glfo, debug=self.args.debug, mimic_data_read_length=self.args.mimic_data_read_length, h_corr_line=h_corr_line)

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
            print('    %3d from %s' % (n_to_erode, erosion[2:]), end=' ')
            print('of %s: %15s' % (erosion[0], fragment_before), end=' ')
            print(' --> %-15s' % fragment_after)
        if len(fragment_after) == 0:
            print('    NOTE eroded away entire sequence')

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
            print('  eroding')
        for region in utils.regions:
            reco_event.eroded_seqs[region] = reco_event.original_seqs[region]
        for erosion in utils.real_erosions:
            self.erode(erosion, reco_event)
        if self.args.mimic_data_read_length:
            for erosion in utils.effective_erosions:
                self.erode(erosion, reco_event)

        itry = 0
        reco_event.set_naive_seq(use_dummy_insertions=True)  # see if there's a stop due to stuff other than the new insertions, in which case we can't do anything about it here
        pre_existing_stop = reco_event.is_there_a_stop_codon()
        while itry == 0 or (not pre_existing_stop and reco_event.is_there_a_stop_codon()):  # note that if there's already a stop codon in the non-insert bits, this lets us add additional stop codons in the insertions (but that makes sense, since we don't have a way to tell where the stop codons are [and we don't care])
            for boundary in utils.boundaries:
                self.insert(boundary, reco_event)
            reco_event.set_naive_seq()
            itry += 1
            if itry % 50 == 0:
                print('%s adding insertions is taking an oddly large number of tries (%d so far)' % (utils.color('yellow', 'warning'), itry))

        if self.args.debug:
            print('  joining eroded seqs')
            print('         v: %s' % reco_event.eroded_seqs['v'])
            print('    insert: %s' % reco_event.insertions['vd'])
            print('         d: %s' % reco_event.eroded_seqs['d'])
            print('    insert: %s' % reco_event.insertions['dj'])
            print('         j: %s' % reco_event.eroded_seqs['j'])
        reco_event.set_post_erosion_codon_positions()

    # ----------------------------------------------------------------------------------------
    def write_mute_freqs(self, reco_event, reco_seq_fname, per_base_freqs=None, debug=False):  # unsurprisingly, this function profiles out to be kind of a dumb way to do it, in terms of run time
        # write per-position mute freqs, but also collects per-base (per-ACGT) freqs (<per_base_freqs>) and returns them if they're needed
        # ----------------------------------------------------------------------------------------
        def get_pbfreqs(naive_base, pbcounts=None, inuke=None, rgene=None):
            def def_count(base, count=None):  # default, i.e. if we have no other information (this is used twice, first to set all bases if we have no info, and second [if <count> is set] to set pseudocount values if we don't have enough counts for some/all bases)
                if base == naive_base:
                    return 0
                else:  # but for the other three we just want to set 1 if there's no info or zero counts
                    return 1 if count is None else max(count, 1)
            def min_count(count, tot_counts):  # bppseqgen barfs if any count is too small (and it gets normalized after this so ends up smaller): ParameterException: ConstraintException: Parameter::setValue(0)]1e-06; 0.999999[(HKY85.theta1)
                min_fraction = 0.01  # this is an important parameter -- it determines if it's possible to mutate back to the original naive base. If it's 0, then it's not possible, which I think is what we want, since while that isn't really right, if it's greater than 0 then we'll get the much more common occurrence that we *really* don't want of the original/naive base "mutating" to itself with an initial mutation (whereas if it's 0, that's only wrong like 1/4 of the times the position mutations twice, which we don't care about at all)
                return max(count, min_fraction * tot_counts)  # if count/tot_counts is less than min_fraction, return min_fraction
            if debug and pbcounts is not None:  # originally just to check that we have the right position in the gene and counts (but it happens too much just from random stuff to be worth printing unless debug is one)
                if sum(pbcounts.values()) > 10 and any(c > pbcounts[naive_base] for n, c in pbcounts.items() if n != naive_base):  # ok now that i've actually run with this check, it picks up quite a few cases where I'm presuming we have the wrong germline gene, in which case it's probably better that it's "wrong"? jeez i dunno, doesn't matter
                    print('    %s non-germline base has more counts than germline base (%s) at ipos %s in %s: %s' % (utils.color('red', 'warning'), naive_base, inuke, utils.color_gene(rgene), pbcounts))  # formatting inuke as string on the off chance we get here when the calling fcn doesn't pass it
            if pbcounts is None:
                pbcounts = {n : def_count(n) for n in utils.nukes}
            pbcounts = {n : def_count(n, count=c) for n, c in pbcounts.items()}  # add pseudocounts (NOTE this is quite a bit less involved than in hmmwriter.py process_mutation_info() and get_emission_prob())
            pbcounts = {n : min_count(c, sum(pbcounts.values())) for n, c in pbcounts.items()}  # make sure none of them are too small
            tmptot = sum(pbcounts.values())
            pbcounts = {n : pbcounts[n] / float(tmptot) for n, c in pbcounts.items()}
            return pbcounts
        # ----------------------------------------------------------------------------------------
        def get_region_freqs(region):
            # ----------------------------------------------------------------------------------------
            rfreqs[region] = []  # list with a relative mutation rate for each position in <rseq>
            rtotals[region] = 0.
            rseq = reco_event.eroded_seqs[region]
            rgene = reco_event.genes[region]
            if not self.args.no_per_base_mutation:
                per_base_freqs[region] = []
            if len(rseq) == 0:  # i think this is how it handles light chain d? not checking right now, just copying how it did it below
                return
            mute_freqs = self.get_mute_freqs(rgene)
            if not self.args.no_per_base_mutation:
                mute_counts = self.get_mute_counts(rgene)
            all_erosions = dict(list(reco_event.erosions.items()) + list(reco_event.effective_erosions.items()))  # arg, this is hackey, but I don't want to change Event right now
            cposlist = None
            if not self.args.mutate_conserved_codons and region in utils.conserved_codons[self.args.locus]:
                codon = utils.conserved_codons[self.args.locus][region]
                cpos = self.glfo[codon + '-positions'][rgene]
                cposlist = list(range(cpos, cpos + 3))

            # NOTE <inuke> is position/index in the *eroded* sequence that we're dealing with, while <position> is in the uneroded germline gene
            for inuke in range(len(rseq)):  # append a freq for each nuke
                position = inuke + all_erosions[region + '_5p']
                freq = 0.0
                if position in mute_freqs:
                    freq = mute_freqs[position]
                else:
                    freq = mute_freqs['overall_mean']
                if cposlist is not None and position in cposlist:  # set freq for conserved codons to zero
                    freq = 0.0
                rfreqs[region].append(freq)
                rtotals[region] += freq
                if not self.args.no_per_base_mutation:
                    per_base_freqs[region].append(get_pbfreqs(rseq[inuke], pbcounts=mute_counts.get(position), inuke=inuke, rgene=rgene))

            # normalize to the number of sites in this region, i.e. so an average site is given value 1.0 (I'm not sure that this really needs to be done, but it might not be exactly normalized before this)
            assert rtotals[region] != 0.0
            for inuke in range(len(rseq)):
                rfreqs[region][inuke] *= float(len(rseq)) / rtotals[region]
            rtotals[region] = float(len(rseq))  # dammit i keep forgetting this: if you're going to use the total later, you obviously have to reset it after normalizing

            # then rescale by regional branch lengths
            rlength_ratio = self.treeinfo['branch-length-ratios'][region]
            for inuke in range(len(rseq)):
                rfreqs[region][inuke] *= rlength_ratio  # so it's no longer normalized after this -- we renormalize over all regions afterward
            rtotals[region] *= rlength_ratio

            if debug:
                print('    %s: read freqs for %d positions: %d to %d' % (region, len(rfreqs[region]), all_erosions[region + '_5p'], len(rseq) - 1 + all_erosions[region + '_5p']))
                print('          normalized to 1, then multiplied by %.3f (total %.3f)' % (rlength_ratio, rtotals[region]))

        # ----------------------------------------------------------------------------------------
        rfreqs, rtotals = {}, {}
        # first add each region
        for region in utils.regions:
            get_region_freqs(region)
        # then add insertions
        mean_freq = numpy.mean(rfreqs['d' if utils.has_d_gene(self.args.locus) else 'v'])  # use mute freq from d for heavy chain, v for light chain (i wish i had the cdr3 mute freq here, but i don't)
        for bound in utils.boundaries:
            rfreqs[bound] = [mean_freq for _ in range(len(reco_event.insertions[bound]))]
            rtotals[bound] = mean_freq * len(reco_event.insertions[bound])
            if debug:
                print('    %s: added %d positions with freq %.3f from %s' % (bound, len(rfreqs[bound]), mean_freq, 'd' if utils.has_d_gene(self.args.locus) else 'v'))
            if not self.args.no_per_base_mutation:
                per_base_freqs[bound] = [get_pbfreqs(nb) for nb in reco_event.insertions[bound]]

        final_freqs = [f for r in ['v', 'vd', 'd', 'dj', 'j'] for f in rfreqs[r]]
        final_seq = reco_event.recombined_seq
        assert len(final_freqs) == len(final_seq)
        if not self.args.no_per_base_mutation:
            per_base_freqs['final'] = [f for r in ['v', 'vd', 'd', 'dj', 'j'] for f in per_base_freqs[r]]
            assert len(per_base_freqs['final']) == len(final_seq)

        # normalize to the number of sites (i.e. so an average site is given value 1.0)
        final_total = sum(rtotals.values())
        for inuke in range(len(final_seq)):
            final_freqs[inuke] *= float(len(final_seq)) / final_total

        # you might think you can remove this, but i can tell you from experience that you really, really shouldn't
        final_total = 0.0
        for inuke in range(len(final_seq)):
            final_total += final_freqs[inuke]
        assert utils.is_normed(final_total / float(len(final_seq)))

        # write the input file for bppseqgen, one base per line
        with open(reco_seq_fname, 'w') as reco_seq_file:
            # NOTE really not sure why this doesn't really [seems to require an "extra" column] work with csv DictWriter, but it doesn't -- bppseqgen barfs (I think maybe it expects a different newline character? don't feel like working it out)
            headstr = 'state'
            if not self.args.mutate_from_scratch:
                headstr += '\trate'
            reco_seq_file.write(headstr + '\n')
            for inuke in range(len(final_seq)):
                linestr = final_seq[inuke]
                if not self.args.mutate_from_scratch:
                    linestr += '\t%f' % final_freqs[inuke]
                reco_seq_file.write(linestr + '\n')
        if debug:
            print('    wrote %d positions%s to %s' % (len(final_seq), '' if self.args.mutate_from_scratch else ' with per-position rates', reco_seq_fname))

    # ----------------------------------------------------------------------------------------
    def prepare_bppseqgen(self, cmdfos, reco_event, seed):
        """ write input files and get command line options necessary to run bppseqgen on <seq> (which is a part of the full query sequence) """
        tmptree = copy.deepcopy(reco_event.tree)  # we really don't want to modify the tree in the event
        if not self.args.no_per_base_mutation:
            if self.args.debug:
                print('      rescaling tree by %.2f to account for non-equilibrium per-base mutation' % self.per_base_mutation_multiplier)
            tmptree.scale_edges(self.per_base_mutation_multiplier)
        for node in tmptree.preorder_internal_node_iter():  # bppseqgen barfs if any node labels aren't of form t<N>, so we have to de-label all the internal nodes, which have been labelled by the code in treeutils
            node.taxon = None
        chosen_tree = tmptree.as_string(schema='newick').strip()
        chosen_tree = '(%s,%s:%.15f):0.0;' % (chosen_tree.rstrip(';'), dummy_name_so_bppseqgen_doesnt_break, treeutils.get_mean_leaf_height(treestr=chosen_tree))  # add dummy leaf that we'll subsequently ignore (such are the vagaries of bppseqgen; see https://github.com/BioPP/bppsuite/issues/3)

        workdir = self.workdir + '/bpp'
        treefname = workdir + '/tree.tre'
        reco_seq_fname = workdir + '/start-seq.txt'
        leaf_seq_fname = '%s/leaf-seqs.fa' % self.workdir
        workfnames = [reco_seq_fname, treefname]

        os.makedirs(workdir)
        with open(treefname, 'w') as treefile:
            treefile.write(chosen_tree)
        per_base_freqs = {} if not self.args.no_per_base_mutation else None
        self.write_mute_freqs(reco_event, reco_seq_fname, per_base_freqs=per_base_freqs)

        bpp_path = '%s/packages/%s' % (self.args.partis_dir, 'bpp-newlik/_build' if not self.args.no_per_base_mutation else 'bpp')
        bpp_binary = '%s/bin/bppseqgen' % bpp_path
        if not os.path.exists(bpp_binary):
            raise Exception('bppseqgen binary not found: %s' % bpp_binary)
        env = os.environ.copy()
        env['LD_LIBRARY_PATH'] = '%s/lib%s' % (bpp_path, (':' + env.get('LD_LIBRARY_PATH')) if env.get('LD_LIBRARY_PATH') is not None else '')

        full_seq = reco_event.recombined_seq
        if not self.args.no_per_base_mutation:
            if self.args.debug:
                print('      using bpp per-base mutation with path %s' % bpp_path)
            paramfname = workdir + '/cfg.bpp'  # docs: http://biopp.univ-montp2.fr/apidoc/bpp-phyl/html/index.html that page is too darn hard to google
            workfnames.append(paramfname)
            plines = ['alphabet = DNA']
            plines += ['number_of_sites = %d' % len(full_seq)]
            plines += ['input.tree1 = user(file=%s)' % treefname]
            plines += ['input.tree1.format = Newick']  # this is the default, but adding this keeps it from printing a warning
            plines += ['input.infos = %s' % reco_seq_fname]  # init state (i.e. naive base) and rate for each position
            plines += ['input.infos.states = state']  # column name for reco_seq_fname

            if self.args.mutate_from_scratch:  # this isn't per-base mutation, it's non-per-base but using the newlik branch, but i can't get it to work (it's crashing because i'm not quite specifying parameters right, but there's no damn docs, or i can't find them, and i'm tired of guessing), so I'm just going back to the old version. It sucks to carry two bpp versions, but oh well
                # NOTE if you implement this, you'll have to check all the places where not self.args.no_per_base_mutation is used to see if it should be self.args.newlik or something
                raise Exception('can\'t yet mutate from scratch with per-base mutation')
                plines += ['input.infos.rates = none']  # column name for reco_seq_fname  # BEWARE bio++ undocumented defaults (i.e. look in the source code)
                plines += ['model1 = JC69']
                if self.args.flat_mute_freq:
                    plines += ['rate_distribution1 = Constant']
                else:
                    plines += ['rate_distribution1 = Gamma(n=4,alpha=' + self.mute_models['v']['gamma']['alpha'] + ')']  # eh, maybe just use v, it was probably inferred more accurately, i don't think this really varies over regions, and I don't really care about this bit much any more since the per-base stuff will soon be default. Plus, this inference was just taken from that old connor paper so i have no idea how good it is
                plines += ['process1 = Homogeneous(model=model1, tree=1, rate=1)']  # not really sure what the rate does here
                plines += ['simul1 = Single(process=%d, output.sequence.file=%s)' % (len(full_seq) + 1, leaf_seq_fname)]
            else:  # default: per-position read from parameter file
                # NOTE/TODO this successfully gets us per-base mutation rates, but the overall mutation isn't right -- the more asymmetric the rates to the four bases, the higher the tree depth bppseqgen gives back. Not sure why yet
                # updating bio++ to (mostly) newlik branch, maybe roughly the same as v2.4.0?
                # https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/biopp-devel-forum/cX-Q9aks7CA/mW1uQfCy-nUJ
                # general philosophy here: download + compile a specific version of bio++, test it *thoroughly* to make sure it's doing *exactly* the mutation we need, and then only allow use of that specific version (this probably screws people on other operating systems, or at least they'll have to recompile)
                #  - this is needed cause bpp does lots of weird undocumented stuff (e.g. undocumented default values, if you pass it a nonsense/unknown parameter it just ignores it silently instead of telling you)
                plines += ['input.infos.rates = rate']  # column name for reco_seq_fname
                plines += ['rate_distribution1 = Constant()']  # every site has a totally different model, so no point in having a non-constant rates-across-sites
                plines += ['']

                # make a model for each site, with per-base rates (equilibrium/init freqs) as the observed fraction of times we saw each base at that position
                for inuke in range(len(full_seq)):
                    plines += ['model%d = HKY85(kappa=1., initFreqs=values(%s))' % (inuke + 1, ', '.join(['%f' % per_base_freqs['final'][inuke][n] for n in sorted(utils.nukes)]))]  # NOTE bio++ manual says the alphabet is always in alphabetical order, so we assume here it's ACGT (if it isn't, this is all wrong)
                plines += ['']

                # make a homogeneous (same rate for entire tree) process for each site, all with the same tree and rate
                for inuke in range(len(full_seq)):
                    plines += ['process%d = Homogeneous(model=%d, tree=1, rate=1)' % tuple(inuke + 1 for _ in range(2))]  # NOTE I"m not really sure that the rate does anything, since I"m passing input.infos (homogeneous: same rate for whole tree)
                plines += ['']

                # make a final combination process incorporating all the previous ones
                plines += ['process%d = Partition( \\' % (len(full_seq) + 1)]
                for inuke in range(len(full_seq)):
                    plines += ['                     process%d=%d, process%d.sites=%d, \\' % tuple(inuke + 1 for _ in range(4))]
                plines += [')', '']

                # simulation with the final process
                plines += ['simul1 = Single(process=%d, output.sequence.file=%s)' % (len(full_seq) + 1, leaf_seq_fname)]

            with open(paramfname, 'w') as pfile:
                pfile.write('\n'.join(plines))
            command = '%s param=%s --seed=%d' % (bpp_binary, paramfname, seed)  # not sure how to set the seed in the param file, but this works, so oh well
        else:
            if self.args.debug:
                print('      using bpp uniform mutation with path %s' % bpp_path)
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
                if self.args.flat_mute_freq:
                    command += ' rate_distribution=Constant'
                else:
                    command += ' rate_distribution=Gamma(n=4,alpha=' + self.mute_models['v']['gamma']['alpha'] + ')'  # eh, maybe just use v, it was probably inferred more accurately, i don't think this really varies over regions, and I don't really care about this bit much any more since the per-base stuff will soon be default. Plus, this inference was just taken from that old connor paper so i have no idea how good it is
            else:
                command += ' input.infos.rates=rate'  # column name in input file
                pvpairs = [p + '=' + v for p, v in self.mute_models['v']['gtr'].items()]  # see note about using only v a couple lines above
                command += ' model=GTR(' + ','.join(pvpairs) + ')'

        cmdfos.append({'cmd_str' : command, 'outfname' : leaf_seq_fname, 'workdir' : workdir, 'workfnames' : workfnames, 'env' : env})  # used to run each region separately so it made more sense as a list

    # ----------------------------------------------------------------------------------------
    def read_bppseqgen_output(self, cmdfo, reco_event):
        mutated_seqs = {}
        for seqfo in utils.read_fastx(cmdfo['outfname']):  # get the leaf node sequences from the file that bppseqgen wrote
            if seqfo['name'] == dummy_name_so_bppseqgen_doesnt_break:  # in the unlikely (impossible unless we change tree generators and don't tell them to use the same leaf names) event that we get a non-dummy leaf with this name, it'll fail at the assertion just below
                continue
            mutated_seqs[seqfo['name'].strip('\'')] = seqfo['seq']
        try:  # make sure names are all of form t<n>, and keep track of which sequences goes with which name (have to keep around the t<n> labels so we can translate the tree labels, in event.py)
            names_seqs = [('t' + str(iseq + 1), mutated_seqs['t' + str(iseq + 1)]) for iseq in range(len(mutated_seqs))]
        except KeyError as ke:
            raise Exception('leaf name %s not as expected in bppseqgen output %s' % (ke, cmdfo['outfname']))
        assert treeutils.get_n_leaves(reco_event.tree) == len(names_seqs)
        os.remove(cmdfo['outfname'])

        assert len(reco_event.final_seqs) == 0
        reco_event.leaf_names = []  # i'm pretty sure there's a good reason this starts as None but the seqs start as a length zero list, but I don't remember it
        for name, seq in names_seqs:
            if not self.args.mutate_conserved_codons:
                seq = reco_event.revert_conserved_codons(seq, debug=self.args.debug)  # if mutation screwed up the conserved codons, just switch 'em back to what they were to start with UPDATE should really remove this now that we're setting the rates for these positions to zero
            if self.args.mutate_stop_codons:
                seq = reco_event.mutate_away_stop_codons(seq, debug=self.args.debug)
            reco_event.final_seqs.append(seq)
            reco_event.leaf_names.append(name)
            reco_event.final_codon_positions.append(copy.deepcopy(reco_event.post_erosion_codon_positions))  # separate codon positions for each sequence, because of shm indels

    # ----------------------------------------------------------------------------------------
    def add_shm_indels(self, reco_event):
        # NOTE that it will eventually make sense to add shared indel mutation according to the chosen tree -- i.e., probably, with some probability apply an indel instead of a point mutation
        if self.args.debug and self.args.indel_frequency > 0.:
            print('      indels')
        reco_event.indelfos = [indelutils.get_empty_indel() for _ in range(len(reco_event.final_seqs))]
        for iseq in range(len(reco_event.final_seqs)):
            if self.args.indel_frequency == 0.:  # no indels at all
                continue
            if numpy.random.uniform(0, 1) > self.args.indel_frequency:  # no indels for this sequence
                if self.args.debug:
                    print('        0')
                continue
            n_indels = numpy.random.choice(self.args.n_indels_per_indeld_seq)
            input_seq, indelfo = indelutils.add_indels(n_indels, reco_event.final_seqs[iseq], reco_event.recombined_seq,  # NOTE modifies <indelfo> and <codon_positions>
                                                       self.args.mean_indel_length, reco_event.final_codon_positions[iseq], indel_location=self.args.indel_location, dbg_pad=8, debug=self.args.debug)
            reco_event.final_seqs[iseq] = input_seq
            indelfo['genes'] = {r : reco_event.genes[r] for r in utils.regions}
            reco_event.indelfos[iseq] = indelfo

    # ----------------------------------------------------------------------------------------
    def add_mutants(self, reco_event, irandom, i_choose_tree=None):
        if self.args.mutation_multiplier is not None and self.args.mutation_multiplier == 0.:  # some of the stuff below fails if mut mult is actually 0.
            reco_event.final_seqs.append(reco_event.recombined_seq)  # set final sequnce in reco_event
            reco_event.indelfos = [indelutils.get_empty_indel() for _ in range(len(reco_event.final_seqs))]
            # reco_event.setline(irandom)  # doesn't work
            if self.args.rearrange_from_scratch:
                print('  %s not setting reco_event.line, which may cause a crash. If you\'re setting --mutation-multiplier to exactly zero you can fix this by setting it to arbitrary very small value e.g. 0.00000001 (sorry)' % utils.wrnstr())
            return

        # When generating trees, each tree's number of leaves and total depth are chosen from the specified distributions (a.t.m., by default n-leaves is from a geometric/zipf, and depth is from data)
        # This chosen depth corresponds to the sequence-wide mutation frequency (the newick trees have branch lengths corresponding to the whole sequence  (i.e. the weighted mean of v, d, and j))
        # In order to account for varying mutation rates in v, d, and j we also get the repertoire-wide ratio of mutation freqs for each region from treegenerator
        # We used to make a separate tree for each region, and rescale that tree by the appropriate ratio and simulate with three separate bppseqgen processes, but now we use one tree for the whole sequence, and do the rescaling when we write the per-position mutation rates for each region
        if i_choose_tree is not None:
            if i_choose_tree >= len(self.treeinfo['trees']):
                if self.args.debug:  # it should have warned about this (well, more events than trees) already when first dealing with trees, before running anything
                    print('      i_choose_tree %d larger than n trees %d, %%\'ing down to %d' % (i_choose_tree, len(self.treeinfo['trees']), i_choose_tree % len(self.treeinfo['trees'])))
                i_choose_tree = i_choose_tree % len(self.treeinfo['trees'])
            itree = i_choose_tree
        else:
            itree = random.randint(0, len(self.treeinfo['trees'])-1)
        chosen_treestr = self.treeinfo['trees'][itree]
        if self.args.debug:
            print('    choosing itree %d (of %d)' % (itree, len(self.treeinfo['trees'])))
        reco_event.set_tree(chosen_treestr)  # leaf names are still just like t<n>
        if self.args.mutation_multiplier is not None:
            reco_event.tree.scale_edges(self.args.mutation_multiplier)

        if self.args.debug:
            mheight = treeutils.get_mean_leaf_height(tree=reco_event.tree)
            print('  chose tree with total height %f%s' % (mheight, (' (includes factor %.2f from --mutation-multiplier)' % self.args.mutation_multiplier) if self.args.mutation_multiplier is not None else ''))
            print('    regional heights:  %s' % ('   '.join(['%s %.3f' % (r, mheight * self.treeinfo['branch-length-ratios'][r]) for r in utils.regions])))

        cmdfos = []
        self.prepare_bppseqgen(cmdfos, reco_event, seed=irandom)
        assert len(cmdfos) == 1  # used to be one cmd for each region

        start = time.time()
        utils.run_cmds(cmdfos, sleep=False, clean_on_success=True)
        self.validation_values['bpp-times'].append(time.time()-start)

        self.read_bppseqgen_output(cmdfos[0], reco_event)

        self.add_shm_indels(reco_event)
        reco_event.setline(irandom)  # set the line here because we use it when checking tree simulation, and want to make sure the uids are always set at the same point in the workflow
        if self.args.check_tree_depths or self.args.debug:
            self.check_tree_simulation(reco_event)

        if self.args.debug:
            print('  bppseqgen ran on the following tree (mean depth %.3f, imbalance %.4f) in %.2fs:' % (treeutils.get_mean_leaf_height(tree=reco_event.tree), treeutils.get_imbalance(reco_event.tree), self.validation_values['bpp-times'][-1]))
            print(treeutils.get_ascii_tree(dendro_tree=reco_event.tree, extra_str='      '))
            utils.print_reco_event(reco_event.line, extra_str='    ')

    # ----------------------------------------------------------------------------------------
    def check_tree_simulation(self, reco_event, debug=False):  # also adds validation values for this event, so you can later print them for all of the events
        # NOTE turning on this debug just tells you the values for this event, but if you want to average over lots of events you turn on the print_validation_values() call in bin/partis
        ltmp = reco_event.line
        mheight = treeutils.get_mean_leaf_height(tree=reco_event.tree)
        if debug:
            print('          in       out    (input tree heights vs output fraction of positions mutated)')
        for rname in ['all'] + utils.regions:
            input_height = mheight * (1 if rname=='all' else self.treeinfo['branch-length-ratios'][rname])
            mean_obs = numpy.mean([utils.get_mutation_rate(ltmp, iseq=i, restrict_to_region='' if rname=='all' else rname) for i in range(len(ltmp['unique_ids']))])  # could use the 'mut_freqs' key to avoid recalculatation, but this is a bit cleaner
            self.validation_values['heights'][rname]['in'].append(input_height)
            self.validation_values['heights'][rname]['out'].append(mean_obs)
            if debug:
                print('  %4s %7.3f  %7.3f' % (rname, input_height, mean_obs))

        if self.args.debug:  # <debug> above is for if we're debugging this function, whereas we want to print this stuff when --debug is set
            if any(indelutils.has_indels(reco_event.line['indelfos'][i]) for i in range(len(reco_event.line['unique_ids']))):
                print('    skipping tree difference metrics since there\'s shm indels (would have to align before passing to fasttree)')
            else:
                treeutils.compare_input_tree_to_leaf_seqs('all', reco_event.line['tree'], reco_event.final_seqs, ltmp['naive_seq'])  # NOTE we want to pass in the treestr rather than the dendro tree because a) the stupid dendropy functions modify the dendro tree you give them b) they also need the two trees to have the same taxon namespace so we'd have to make a new one anyway

    # ----------------------------------------------------------------------------------------
    def print_validation_values(self):  # the point of having this separate from the previous function is that it's run after you've simulated lots of different events (unlike everything else in this class, the validation values aggregate over different events)
        print('  input tree height vs output mut frac (means over leaves/seqs) averaged over %d events:' % len(self.validation_values['heights']['all']['in']))
        print('             in      out       diff    std err   std dev')
        for vtype in ['all'] + utils.regions:
            vvals = self.validation_values['heights'][vtype]
            deltas = [(vvals['out'][i] - vvals['in'][i]) for i in range(len(vvals['in']))]
            std_val = numpy.std(deltas, ddof=1) if len(deltas) > 1 else float('nan')
            print('       %3s  %.3f   %.3f    %+.3f +/- %.3f     %.3f' % (vtype, numpy.mean(vvals['in']), numpy.mean(vvals['out']), numpy.mean(deltas), std_val / math.sqrt(len(deltas)), std_val))  # NOTE each delta is already the mean of <n_leaves> (non-independent) measurements
