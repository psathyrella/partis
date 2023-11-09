from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import csv
import time
import sys
import itertools

from . import utils
from . import glutils
from .hist import Hist
from . import hutils
from . import plotconfig
from .mutefreqer import MuteFreqer
from .corrcounter import CorrCounter
from . import indelutils
from io import open

# ----------------------------------------------------------------------------------------
class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self, glfo, args, count_correlations=False):
        self.glfo = glfo
        self.args = args
        self.mfreqer = MuteFreqer(self.glfo, exclusions=args.region_end_exclusions, args=args)
        self.corrcounter = CorrCounter() if count_correlations else None  # NOTE if you want paired h/l correlations you have to use corr counter outside of the parameter counter
        self.reco_total = 0  # total number of recombination events
        self.mute_total = 0  # total number of sequences
        self.counts = {}
        self.counts['all'] = {}
        for column in utils.column_dependencies:
            self.counts[column] = {}
        self.string_columns = set([r + '_gene' for r in utils.regions])
        for bound in utils.boundaries:
            self.counts[bound + '_insertion_content'] = {n : 0 for n in utils.nukes}  # base content of each insertion
            self.string_columns.add(bound + '_insertion_content')
        self.counts['aa_cdr3_length'] = {}
        self.counts['non_vj_length'] = {}
        self.counts['seq_content'] = {n : 0 for n in utils.nukes}  # now I'm adding the aa content, I wish this had nucleotide in the name, but I don't want to change it since it corresponds to a million existing file paths
        self.counts['cluster_size'] = {}
        self.counts['net_indel_length'] = {}
        self.init_aa_stuff()
        self.counts['seq_aa_content'] = {a : 0 for a in self.all_aa}
        self.string_columns.add('seq_content')
        self.string_columns.add('seq_aa_content')

        self.no_write_columns = ['aa_cdr3_length', 'non_vj_length', 'seq_aa_content', 'net_indel_length']  # don't write these to the parameter dir, since a) cdr3 length is better viewed as an output of more fundamental parameters (gene choice, insertion + deletion lengths) and b) I"m adding them waaay long after the others, and I don't want to add a new file to the established parameter directory structure. (I'm adding these because I want them plotted)

        self.columns_to_subset_by_gene = [e + '_del' for e in utils.all_erosions] + [b + '_insertion' for b in utils.boundaries]
        self.mean_columns = ['aa_cdr3_length', 'non_vj_length']
        self.log_columns = ['cluster_size', 'net_indel_length']

    # ----------------------------------------------------------------------------------------
    def init_aa_stuff(self):
        codons = itertools.product(utils.nukes + ['N'], repeat=3)  # I cannot for the life of me find anything in Bio that will give me the list of amino acids, wtf, but I'm tired of googling, this will be fine
        self.all_aa = set([utils.ltranslate(''.join(c)) for c in codons])

    # ----------------------------------------------------------------------------------------
    def get_index(self, info, deps):
        index = []
        for ic in deps:
            if ic[2:] == '_insertion':  # insertion length
                index.append(len(info[ic]))
            else:
                assert 'insertion' not in ic
                assert 'content' not in ic
                index.append(info[ic])
        return tuple(index)

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        self.increment_per_family_params(info)
        for iseq in range(len(info['seqs'])):
            self.increment_per_sequence_params(info, iseq)

    # ----------------------------------------------------------------------------------------
    def sub_increment(self, column, index):
        if index not in self.counts[column]:
            self.counts[column][index] = 0
        self.counts[column][index] += 1

    # ----------------------------------------------------------------------------------------
    def increment_per_sequence_params(self, info, iseq):
        """ increment parameters that differ for each sequence within the clonal family """
        self.mute_total += 1
        self.mfreqer.increment(info, iseq)
        for nuke in utils.nukes:
            self.counts['seq_content'][nuke] += info['seqs'][iseq].count(nuke)

        # aa seq content stuff
        nseq = info['seqs'][iseq]
        if info['v_5p_del'] > 0:
            nseq = info['v_5p_del'] * utils.ambig_base + nseq
        if len(info['fv_insertion']) > 0:
            nseq = nseq[len(info['fv_insertion']) :]
        if len(nseq) % 3 != 0:
            nseq += utils.ambig_base * (3 - (len(nseq) % 3))  # I think I could replace this with the new utils.pad_nuc_seq()
        aaseq = utils.ltranslate(nseq)
        for aa in self.all_aa:
            self.counts['seq_aa_content'][aa] += aaseq.count(aa)

        self.sub_increment('net_indel_length', (indelutils.net_length(info['indelfos'][iseq]), ))

    # ----------------------------------------------------------------------------------------
    def increment_per_family_params(self, info):
        """ increment parameters that are the same for the entire clonal family """
        self.reco_total += 1

        all_index = self.get_index(info, tuple(list(utils.index_columns) + ['cdr3_length', ]))  # NOTE this cdr3_length is for getting a unique index for the rearrangement event parameters, and is thus unrelated to the key aa_cdr3_length for plotting
        if all_index not in self.counts['all']:
            self.counts['all'][all_index] = 0
        self.counts['all'][all_index] += 1

        for deps in utils.column_dependency_tuples:
            column = deps[0]
            index = self.get_index(info, deps)
            self.sub_increment(column, index)

        # have to be done separately, since they're not index columns (and we don't want them to be, since they're better viewed as derivative -- see note in self.write())
        self.sub_increment('aa_cdr3_length', (info['cdr3_length'] // 3, ))  # oh, jeez, this has to be a tuple to match the index columns, that's ugly
        self.sub_increment('non_vj_length', (utils.get_non_vj_len(info), ))
        self.sub_increment('cluster_size', (len(info['unique_ids']), ))

        for bound in utils.boundaries:
            for nuke in info[bound + '_insertion']:
                if nuke == utils.ambig_base:
                    continue
                self.counts[bound + '_insertion_content'][nuke] += 1

        if self.corrcounter is not None:
            self.corrcounter.increment(info)

    # ----------------------------------------------------------------------------------------
    def clean_plots(self, plotdir):
        self.mfreqer.clean_plots(plotdir + '/mute-freqs')
        if self.corrcounter is not None:
            self.corrcounter.clean_plots(plotdir + '/correlations')
        utils.prep_dir(plotdir + '/overall', wildlings=('*.csv', '*.svg'))
        for column in self.counts:
            if column in self.columns_to_subset_by_gene:
                thisplotdir = plotdir + '/' + column
                utils.prep_dir(thisplotdir, wildlings=['*.csv', '*.svg'])

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False, only_overall=False, make_per_base_plots=False):  # NOTE most of the time in here is taken up by mutefrequer.finalize() (if it write() wasn't called first, that is)
        from . import plotting
        print('  plotting parameters in %s' % plotdir, end=' ')
        sys.stdout.flush()
        start = time.time()

        self.clean_plots(plotdir)

        self.mfreqer.plot(plotdir + '/mute-freqs', only_csv=only_csv, only_overall=only_overall, make_per_base_plots=make_per_base_plots)
        if self.corrcounter is not None:
            self.corrcounter.plot(plotdir + '/correlations', only_csv=only_csv)

        overall_plotdir = plotdir + '/overall'

        for column in self.counts:
            if column == 'all':
                continue
            values, gene_values = {}, {}
            for index, count in self.counts[column].items():
                column_val = index[0]

                if column_val not in values:
                    values[column_val] = 0.0
                values[column_val] += count

                if column in self.columns_to_subset_by_gene:
                    gene = index[1]  # NOTE this is hackey, but it works find now and will fail obviously if I ever change the correlations to be incompatible. so screw it
                    utils.split_gene(gene)  # checks validity of gene
                    if gene not in gene_values:
                        gene_values[gene] = {}
                    if column_val not in gene_values[gene]:
                        gene_values[gene][column_val] = 0.0
                    gene_values[gene][column_val] += count

            var_type = 'string' if column in self.string_columns else 'int'

            hist = hutils.make_hist_from_dict_of_counts(values, var_type, column)
            plotting.draw_no_root(hist, plotname=column, plotdir=overall_plotdir, xtitle=plotconfig.xtitles.get(column, column), plottitle=plotconfig.plot_titles.get(column, column),
                                  errors=True, write_csv=True, only_csv=only_csv, stats='mean' if column in self.mean_columns else None, log='y' if column in self.log_columns else '', normalize=True)

            if column in self.columns_to_subset_by_gene and not only_overall:
                thisplotdir = plotdir + '/' + column
                for gene in gene_values:
                    plotname = utils.sanitize_name(gene) + '-' + column
                    hist = hutils.make_hist_from_dict_of_counts(gene_values[gene], var_type, plotname)
                    plotting.draw_no_root(hist, plotname=plotname, plotdir=thisplotdir, xtitle=plotconfig.plot_titles.get(column, column), plottitle=gene, errors=True, write_csv=True, only_csv=only_csv)
                if not only_csv:
                    plotting.make_html(thisplotdir)

        if not only_csv:
            plotting.make_html(overall_plotdir)

        print('(%.1f sec)' % (time.time()-start))

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir):  # NOTE most of the time in here is taken up by mutefrequer.finalize() (if it plot() wasn't called first, that is)
        print('    writing parameters to %s' % base_outdir, end=' ')
        sys.stdout.flush()
        start = time.time()

        if os.path.exists(base_outdir + '/' + glutils.glfo_dir):
            for tmploc in [l for l in utils.loci if os.path.exists(base_outdir + '/' + glutils.glfo_dir + '/' + l)]:
                glutils.remove_glfo_files(base_outdir + '/' + glutils.glfo_dir, tmploc, print_warning=False)
        subdirs = ['hmms', 'mute-freqs', 'correlations', glutils.glfo_dir]  # need to clean correlations even if we're not writing it, since they might already be there from a previous run
        utils.prep_dir(base_outdir, subdirs=subdirs, wildlings=('*.csv', '*.yaml', '*.fasta'))  # it's kind of hackey to specify the /hmms dir here, but as soon as we write the parameters below, the previous yamels are out of date, so it's pretty much necessary
        if self.corrcounter is None and os.path.isdir(base_outdir + '/correlations'):  # but if it is there, we want to rm it
            os.rmdir(base_outdir + '/correlations')

        self.mfreqer.write(base_outdir, base_outdir + '/mute-freqs')
        genes_with_counts = [g[0] for r in utils.regions for g in self.counts[r + '_gene'].keys()]
        glutils.write_glfo(base_outdir + '/' + glutils.glfo_dir, self.glfo, only_genes=genes_with_counts, debug=False)

        for column in self.counts:
            index = None
            outfname = None
            if column in self.no_write_columns:
                continue
            elif column == 'all':
                index = tuple(list(utils.index_columns) + ['cdr3_length', ])
                outfname = base_outdir + '/' + utils.get_parameter_fname(column='all')
            elif '_content' in column or column == 'cluster_size':
                index = [column,]
                outfname = base_outdir + '/' + column + '.csv'
            else:
                index = [column,] + utils.column_dependencies[column]
                outfname = base_outdir + '/' + utils.get_parameter_fname(column_and_deps=index)
            if os.path.isfile(outfname):
                os.remove(outfname)
            elif not os.path.exists(base_outdir):
                os.makedirs(base_outdir)
            with open(outfname, utils.csv_wmode()) as outfile:
                out_fieldnames = list(index)
                out_fieldnames.append('count')
                out_data = csv.DictWriter(outfile, out_fieldnames)
                out_data.writeheader()
                # NOTE this will in general not be sorted
                for key, count in self.counts[column].items():
                    line = {}
                    for ic in range(len(key)):
                        line[index[ic]] = key[ic]
                    line['count'] = count
                    out_data.writerow(line)

        print('(%.1f sec)' % (time.time()-start))
