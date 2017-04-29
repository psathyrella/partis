import os
import csv
import time
import sys

import utils
import glutils
from hist import Hist
import plotconfig
from mutefreqer import MuteFreqer

# ----------------------------------------------------------------------------------------
class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args
        self.mfreqer = MuteFreqer(self.glfo)
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
        self.counts['seq_content'] = {n : 0 for n in utils.nukes}
        self.string_columns.add('seq_content')

        self.columns_to_subset_by_gene = [e + '_del' for e in utils.all_erosions] + [b + '_insertion' for b in utils.boundaries]

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
    def increment_per_sequence_params(self, info, iseq):
        """ increment parameters that differ for each sequence within the clonal family """
        self.mute_total += 1
        self.mfreqer.increment(info, iseq)
        for nuke in info['seqs'][iseq]:
            if nuke in utils.ambiguous_bases:
                continue
            self.counts['seq_content'][nuke] += 1

    # ----------------------------------------------------------------------------------------
    def increment_per_family_params(self, info):
        """ increment parameters that are the same for the entire clonal family """
        self.reco_total += 1

        all_index = self.get_index(info, tuple(list(utils.index_columns) + ['cdr3_length', ]))
        if all_index not in self.counts['all']:
            self.counts['all'][all_index] = 0
        self.counts['all'][all_index] += 1

        for deps in utils.column_dependency_tuples:
            column = deps[0]
            index = self.get_index(info, deps)
            if index not in self.counts[column]:
                self.counts[column][index] = 0
            self.counts[column][index] += 1

        for bound in utils.boundaries:
            for nuke in info[bound + '_insertion']:
                if nuke in utils.ambiguous_bases:
                    continue
                self.counts[bound + '_insertion_content'][nuke] += 1

    # ----------------------------------------------------------------------------------------
    def clean_plots(self, plotdir):
        self.mfreqer.clean_plots(plotdir + '/mute-freqs')
        utils.prep_dir(plotdir + '/overall')  #, multilings=('*.csv', '*.svg'))
        for column in self.counts:
            if column in self.columns_to_subset_by_gene:
                thisplotdir = plotdir + '/' + column
                utils.prep_dir(thisplotdir, wildlings=['*.csv', '*.svg'])

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False, only_overall=False):
        import plotting
        print '  plotting parameters',
        sys.stdout.flush()
        start = time.time()

        self.clean_plots(plotdir)

        self.mfreqer.plot(plotdir + '/mute-freqs', only_csv=only_csv, only_overall=only_overall)

        overall_plotdir = plotdir + '/overall'

        for column in self.counts:
            if column == 'all':
                continue
            values, gene_values = {}, {}
            for index, count in self.counts[column].iteritems():
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

            hist = plotting.make_hist_from_dict_of_counts(values, var_type, column, sort=True)
            plotting.draw_no_root(hist, plotname=column, plotdir=overall_plotdir, xtitle=plotconfig.xtitles.get(column, column), plottitle=plotconfig.plot_titles.get(column, column), errors=True, write_csv=True, only_csv=only_csv)

            if column in self.columns_to_subset_by_gene and not only_overall:
                thisplotdir = plotdir + '/' + column
                for gene in gene_values:
                    plotname = utils.sanitize_name(gene) + '-' + column
                    hist = plotting.make_hist_from_dict_of_counts(gene_values[gene], var_type, plotname, sort=True)
                    plotting.draw_no_root(hist, plotname=plotname, plotdir=thisplotdir, xtitle=plotconfig.plot_titles.get(column, column), plottitle=gene, errors=True, write_csv=True, only_csv=only_csv)
                if not only_csv:
                    plotting.make_html(thisplotdir)

        if not only_csv:
            plotting.make_html(overall_plotdir)

        print '(%.1f sec)' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir):
        print '    writing parameters',
        sys.stdout.flush()
        start = time.time()

        if os.path.exists(base_outdir + '/' + glutils.glfo_dir):
            glutils.remove_glfo_files(base_outdir + '/' + glutils.glfo_dir, self.glfo['locus'])  # NOTE I think this will fail if I ever start having multiple loci in one dir
        utils.prep_dir(base_outdir, subdirs=('hmms', 'mute-freqs'), wildlings=('*.csv', '*.yaml', '*.fasta'))  # it's kind of hackey to specify the /hmms dir here, but as soon as we write the parameters below, the previous yamels are out of date, so it's pretty much necessary

        self.mfreqer.write(base_outdir + '/mute-freqs', mean_freq_outfname=base_outdir + '/REGION-mean-mute-freqs.csv')  # REGION is replace by each region in the three output files)
        genes_with_counts = [g[0] for r in utils.regions for g in self.counts[r + '_gene'].keys()]
        glutils.write_glfo(base_outdir + '/' + glutils.glfo_dir, self.glfo, only_genes=genes_with_counts, debug=False)

        for column in self.counts:
            index = None
            outfname = None
            if column == 'all':
                index = tuple(list(utils.index_columns) + ['cdr3_length', ])
                outfname = base_outdir + '/' + utils.get_parameter_fname(column='all')
            elif '_content' in column:
                index = [column,]
                outfname = base_outdir + '/' + column + '.csv'
            else:
                index = [column,] + utils.column_dependencies[column]
                outfname = base_outdir + '/' + utils.get_parameter_fname(column_and_deps=index)
            if os.path.isfile(outfname):
                os.remove(outfname)
            elif not os.path.exists(base_outdir):
                os.makedirs(base_outdir)
            with open(outfname, 'w') as outfile:
                out_fieldnames = list(index)
                out_fieldnames.append('count')
                out_data = csv.DictWriter(outfile, out_fieldnames)
                out_data.writeheader()
                # NOTE this will in general not be sorted
                for key, count in self.counts[column].iteritems():
                    line = {}
                    for ic in range(len(key)):
                        line[index[ic]] = key[ic]
                    line['count'] = count
                    out_data.writerow(line)

        print '(%.1f sec)' % (time.time()-start)
