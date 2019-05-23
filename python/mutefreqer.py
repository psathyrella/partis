import sys
import multiprocessing
import csv
import os

from hist import Hist
import utils
import glutils
# import paramutils

# ----------------------------------------------------------------------------------------
class MuteFreqer(object):
    def __init__(self, glfo, exclusions, calculate_uncertainty=True):
        self.glfo = glfo
        self.exclusions = exclusions
        self.calculate_uncertainty = calculate_uncertainty

        self.counts, self.freqs = {}, {}  # per-gene, per-position counts/rates
        self.n_bins, self.xmin, self.xmax = 40, 0., 0.4
        self.mean_rates = {n : Hist(self.n_bins, self.xmin, self.xmax, xtitle='mut freq', ytitle='freq', title='full seq' if n == 'all' else n.upper())
                           for n in ['all', 'cdr3'] + utils.regions}  # kind of annoying that it's 'all' here, but '' in plotconfig.rstrings (and therefore in performanceplotter) but it's too much trouble to change to '' here
        nmaxes = {'v' : 55, 'd' : 14, 'j' : 15, 'cdr3' : 20, 'all' : 90}
        # note that we want the mean rates to go up really high, so there aren't many overflows, but the n-muted we can't get away with having a huge range because there has to be a bin for each n-muted (to avoid aliasing)
        self.mean_n_muted = {n : Hist(nmaxes[n] + 1, -0.5, nmaxes[n] + 0.5, xtitle='n mutated', ytitle='counts', title='full seq' if n == 'all' else n.upper())
                             for n in ['all', 'cdr3'] + utils.regions}
        self.per_gene_mean_rates = {}

        self.finalized = False
        self.n_cached, self.n_not_cached = 0, 0

        self.subplotdirs = ['overall', ] + ['per-gene/' + r for r in utils.regions] + ['per-gene-per-position/' + r for r in utils.regions]  # + ['per-gene-per-position-per-base/' + r for r in utils.regions]

    # ----------------------------------------------------------------------------------------
    def increment(self, info, iseq):
        singlefo = utils.synthesize_single_seq_line(info, iseq)
        dummy_iseq = 0

        freq, n_muted = utils.get_mutation_rate_and_n_muted(singlefo, dummy_iseq)
        self.mean_rates['all'].fill(freq)  # mean freq over whole sequence (excluding insertions)
        self.mean_n_muted['all'].fill(n_muted)

        freq, n_muted = utils.get_mutation_rate_and_n_muted(singlefo, dummy_iseq, restrict_to_region='cdr3')
        self.mean_rates['cdr3'].fill(freq)
        self.mean_n_muted['cdr3'].fill(n_muted)

        for region in utils.regions:
            # first do mean freqs
            regional_freq, regional_n_muted  = utils.get_mutation_rate_and_n_muted(singlefo, dummy_iseq, restrict_to_region=region)
            self.mean_rates[region].fill(regional_freq)  # per-region mean freq
            self.mean_n_muted[region].fill(regional_n_muted)

            # then do per-gene and per-gene-per-position freqs
            gene = singlefo[region + '_gene']
            if gene not in self.counts:
                self.counts[gene] = {}
                self.per_gene_mean_rates[gene] = Hist(self.n_bins, self.xmin, self.xmax, xtitle='mut freq', ytitle='freq', title=gene)
            self.per_gene_mean_rates[gene].fill(regional_freq)

            gcts = self.counts[gene]  # shorthand name

            assert len(singlefo[region + '_qr_seqs']) == 1  # don't really need this anymore since we're using utils.synthesize_single_seq_line(), but it's nice to have some explicit check immediatley before using [0] below
            germline_seq = singlefo[region + '_gl_seq']
            query_seq = singlefo[region + '_qr_seqs'][0]
            assert len(germline_seq) == len(query_seq)

            istart = self.exclusions[region][0]
            istop = len(germline_seq) - self.exclusions[region][1]
            for ipos in range(istart, istop):  # NOTE this is similar to the stuff in allelefinder, except in allelefinder we need every single sequence to be the same length (so they go in the correct [comparable?] bin), whereas here we do not
                igl = ipos + int(singlefo[region + '_5p_del'])  # account for left-side deletions in the indexing

                if germline_seq[ipos] in utils.ambiguous_bases or query_seq[ipos] in utils.ambiguous_bases:  # skip if either germline or query sequence is ambiguous at this position
                    continue

                if igl not in gcts:  # if we have not yet observed this position in a query sequence, initialize it
                    gcts[igl] = {n : 0 for n in utils.nukes + ['total', ]}
                    gcts[igl]['gl_nuke'] = germline_seq[ipos]

                gcts[igl]['total'] += 1
                gcts[igl][query_seq[ipos]] += 1  # note that if <query_seq[ipos]> isn't among <utils.nukes>, this will toss a key error

    # ----------------------------------------------------------------------------------------
    def get_uncertainty(self, obs, total):
        import fraction_uncertainty
        if self.calculate_uncertainty:  # it's kinda slow
            errs = fraction_uncertainty.err(obs, total)
            # if errs[2]:
            #     self.n_cached += 1
            # else:
            #     self.n_not_cached += 1
        else:
            errs = 0., 1.

        return errs[0], errs[1]

    # ----------------------------------------------------------------------------------------
    def finalize(self):
        """ convert from counts to mut freqs """
        assert not self.finalized

        for gene in self.counts:
            freqs = {position : {} for position in self.counts[gene]}
            for position in sorted(self.counts[gene].keys()):
                n_conserved, n_mutated = 0, 0
                for nuke in utils.nukes:
                    ncount, total = self.counts[gene][position][nuke], self.counts[gene][position]['total']
                    nuke_freq = float(ncount) / total
                    freqs[position][nuke] = nuke_freq
                    freqs[position][nuke + '_lo_err'], freqs[position][nuke + '_hi_err'] = self.get_uncertainty(ncount, total)
                    if nuke == self.counts[gene][position]['gl_nuke']:
                        n_conserved += ncount
                    else:
                        n_mutated += ncount  # sum over A,C,G,T
                freqs[position]['freq'] = float(n_mutated) / total
                freqs[position]['freq_lo_err'], freqs[position]['freq_hi_err'] = self.get_uncertainty(n_mutated, total)

            self.freqs[gene] = freqs

        for hist in self.mean_rates.values():
            hist.normalize()
        for hist in self.per_gene_mean_rates.values():
            hist.normalize()

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def clean_plots(self, plotdir):
        for substr in self.subplotdirs:
            utils.prep_dir(plotdir + '/' + substr, wildlings=('*.csv', '*.svg'))

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False, only_overall=False, make_per_base_plots=False):
        import plotting
        if not self.finalized:
            self.finalize()

        overall_plotdir = plotdir + '/overall'

        for gene in self.freqs:
            if only_overall:
                continue
            freqs = self.freqs[gene]
            if len(freqs) == 0:
                if gene not in glutils.dummy_d_genes.values():
                    print '    %s no mutefreqer obs for %s' % (utils.color('red', 'warning'), utils.color_gene(gene))
                continue
            sorted_positions = sorted(freqs.keys())
            genehist = Hist(sorted_positions[-1] - sorted_positions[0] + 1, sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5, xtitle='position', ytitle='mut freq', title=gene)
            for position in sorted_positions:
                hi_diff = abs(freqs[position]['freq'] - freqs[position]['freq_hi_err'])
                lo_diff = abs(freqs[position]['freq'] - freqs[position]['freq_lo_err'])
                err = 0.5*(hi_diff + lo_diff)
                genehist.set_ibin(genehist.find_bin(position), freqs[position]['freq'], error=err)
            xline = None
            figsize = [7, 4]
            if utils.get_region(gene) in utils.conserved_codons[self.glfo['locus']]:
                xline = utils.cdn_pos(self.glfo, utils.get_region(gene), gene)
            if utils.get_region(gene) == 'v':
                figsize[0] *= 3.5
            elif utils.get_region(gene) == 'j':
                figsize[0] *= 2
            plotting.draw_no_root(self.per_gene_mean_rates[gene], plotdir=plotdir + '/per-gene/' + utils.get_region(gene), plotname=utils.sanitize_name(gene), errors=True, write_csv=True, only_csv=only_csv, shift_overflows=True)
            # per-position plots:
            plotting.draw_no_root(genehist, plotdir=plotdir + '/per-gene-per-position/' + utils.get_region(gene), plotname=utils.sanitize_name(gene), errors=True, write_csv=True, xline=xline, figsize=figsize, only_csv=only_csv, shift_overflows=True)

            if make_per_base_plots:
                # per-position, per-base plots: (super slow, so commented by default)
                import paramutils
                plotting_info = []
                for pos in sorted_positions:
                    plotting_info.append({
                        'name' : str(pos),  # not sure if this is right, but I think so
                        'nuke_freqs' : {n : freqs[pos][n] for n in utils.nukes},
                    })
                paramutils.make_mutefreq_plot(plotdir + '/' + utils.get_region(gene) + '-per-base', utils.sanitize_name(gene), plotting_info)

        # make mean mute freq hists
        for rstr in ['all', 'cdr3'] + utils.regions:
            if rstr == 'all':
                bounds = (0.0, 0.4)
            else:
                bounds = (0.0, 0.6 if rstr == 'd' else 0.4)
            plotting.draw_no_root(self.mean_rates[rstr], plotname=rstr+'_mean-freq', plotdir=overall_plotdir, stats='mean', bounds=bounds, write_csv=True, only_csv=only_csv, shift_overflows=True)
            plotting.draw_no_root(self.mean_n_muted[rstr], plotname=rstr+'_mean-n-muted', plotdir=overall_plotdir, stats='mean', write_csv=True, only_csv=only_csv, shift_overflows=True)

        if not only_csv:  # write html file and fix permissiions
            for substr in self.subplotdirs:
                plotting.make_html(plotdir + '/' + substr)

    # ----------------------------------------------------------------------------------------
    def write_single_gene(self, gene, outfname):
        gcounts, freqs = self.counts[gene], self.freqs[gene]
        with open(outfname, 'w') as outfile:
            nuke_header = [n + xtra for n in utils.nukes for xtra in ('', '_obs', '_lo_err', '_hi_err')]
            writer = csv.DictWriter(outfile, ('position', 'mute_freq', 'lo_err', 'hi_err') + tuple(nuke_header))
            writer.writeheader()
            for position in sorted(gcounts.keys()):
                row = {'position':position,
                       'mute_freq':freqs[position]['freq'],
                       'lo_err':freqs[position]['freq_lo_err'],
                       'hi_err':freqs[position]['freq_hi_err']}
                for nuke in utils.nukes:
                    row[nuke] = freqs[position][nuke]
                    row[nuke + '_obs'] = gcounts[position][nuke]
                    row[nuke + '_lo_err'] = freqs[position][nuke + '_lo_err']
                    row[nuke + '_hi_err'] = freqs[position][nuke + '_hi_err']
                writer.writerow(row)

    # ----------------------------------------------------------------------------------------
    def write(self, outdir, mean_freq_outfname):
        if not self.finalized:
            self.finalize()

        for gene in self.counts:  # don't bother parallelizing this, all the time is spent finalizing
            self.write_single_gene(gene, outdir + '/' + utils.sanitize_name(gene) + '.csv')

        assert 'REGION' in mean_freq_outfname
        self.mean_rates['all'].write(mean_freq_outfname.replace('REGION', 'all'))  # hackey hackey hackey replacement... *sigh*
        for region in utils.regions:
            self.mean_rates[region].write(mean_freq_outfname.replace('REGION', region))
