""" Process observed counts of recombination parameters into a form useful for sampling. """

import bz2
import csv
import subprocess
import math
import os
import operator

class VersionCounter:
    """ Process observed counts of recombination parameters into a form useful for sampling.

    Input: csv with counts for each vdj gene choice, cdr3_length, and erosion lengths.
    Output: 1) probabilities for each (vdj gene choice + cdr3_length combo), eg: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42) : 0.003
            2) for each gene and end (5p or 3p), probability distribution of erosion lengths, eg: ('IGHV3-23*04', 3p) : {0:0.5, 1:0.3, 2:0.2}
    """
    min_counts = 3  # if a combo has fewer counts than this, ignore it completely

    version_freqs = {}  # Dict with entries like: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42) : 360, i.e. we saw this combination of genes with cdr3_length 42 a total of 360 times
    erosion_counts = {}  # Dict of erosion side, gene name, and number eroded to frequency, e.g.: erosion_counts['5p']['IGHV3-64*04']['3'] = 35
                         #   means that for that gene, we saw 5p erosions of length 3 35 times
    erosion_probs = {}  # Same as erosion_counts, but normalized
    regions = ['v', 'd', 'j']
    erosions = ['5p', '3p']
    gene_names = {}

    def __init__(self, infname, base_outdir):
        self.infname = infname
        self.base_outdir = base_outdir
        if not os.path.exists(self.base_outdir):
            os.makedirs(self.base_outdir)
        for region in self.regions:
            self.gene_names[region] = set()
            self.erosion_counts[region] = {}
        for erosion in self.erosions:
            self.erosion_counts[erosion] = {}

    def count(self):
        """ Run the whole counting procedure. """
        total = self.parse_input()
        self.write_gene_choice_probs(total)
        self.normalize_erosion_counts()
        self.check_normalization(self.erosion_probs)
        self.write_erosion_probs()
 
    def parse_input(self):
        """ Convert input csv file to dicts. Also sum totals for gene freqs. """
        print '  parsing input file %s' % self.infname
        total = 0.0
        with bz2.BZ2File(self.infname) as infile:
            in_data = csv.DictReader(infile)
    #        progress_count = 0
    #        n_lines = 5000000 #subprocess.check_output(['bzcat ', infname, ' | wc'])
            for line in in_data:
    #            if progress_count % 100000 == 0:
    #                print '  %3.1e / %3.0e = %5.2f  %s' % (progress_count, n_lines, float(progress_count) / n_lines, line['count'])
    #            progress_count += 1
                if int(line['count']) < self.min_counts:
                    print '    breaking because count is less than %d.' % self.min_counts
                    print '    NOTE that this assumes input file is sorted!'
                    break
                total += int(line['count'])
                index = (line['v_gene'], line['d_gene'], line['j_gene'], line['cdr3_length'])
                if index not in self.version_freqs:
                    self.version_freqs[index] = 0
                self.version_freqs[index] += int(line['count'])
                for region in self.regions:
                    if line[region + '_gene'] not in self.gene_names[region]:
                        self.gene_names[region].add(line[region + '_gene'])
                    for erosion in self.erosions:
                        erosion_name = region + '_' + erosion + '_del'
                        gene_name = line[region + '_gene']
                        if gene_name not in self.erosion_counts[erosion]:
                            self.erosion_counts[erosion][gene_name] = {}
                        n_eroded = int(line[erosion_name])
                        if n_eroded not in self.erosion_counts[erosion][gene_name]:  # e.g.: if we haven't already seen IGHV3-64*04 with an erosion of length 4 on the 3p end
                            self.erosion_counts[erosion][gene_name][n_eroded] = 0
                        self.erosion_counts[erosion][gene_name][n_eroded] += int(line['count'])
    
        return total
    
    def write_gene_choice_probs(self, total):
        """ Write vdj combo freqs to file. """
        outfname = self.base_outdir + '/' + os.path.basename(self.infname.replace('.csv.bz2','.probs.csv.bz2'))
        print '  writing gene choice probabilities'
        with bz2.BZ2File(outfname, 'w') as outfile:
            out_fieldnames = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'prob']
            out_data = csv.DictWriter(outfile, out_fieldnames)
            out_data.writeheader()
            check_tot = 0.0
            for index,count in self.version_freqs.iteritems():
                prob = float(count) / total
                line = {}
                line['v_gene'] = index[0]
                line['d_gene'] = index[1]
                line['j_gene'] = index[2]
                line['cdr3_length'] = index[3]
                line['prob'] = prob
                check_tot += prob
                out_data.writerow(line)
    
    def normalize_erosion_counts(self):
        """ Normalize the counts in erosion_counts. """
        print '  normalizing erosion counts'
        for erosion in self.erosions:
            self.erosion_probs[erosion] = {}
            for region in self.regions:
                for gene_name in self.gene_names[region]:
                    self.erosion_probs[erosion][gene_name] = {}
                    assert gene_name in self.erosion_counts[erosion]
                    total = 0.0  # total for this gene and this side (e.g. 3p of IGHV3-64*04)
                    for n_eroded,count in self.erosion_counts[erosion][gene_name].iteritems():
                        total += count
                    for n_eroded,count in self.erosion_counts[erosion][gene_name].iteritems():
                        self.erosion_probs[erosion][gene_name][n_eroded] = float(count) / total
    
    def check_normalization(self, erosion_probs):
        """ Silly little double check. """
        for erosion in self.erosions:
            for region in self.regions:
                for gene_name in self.gene_names[region]:
                    total = 0.0
                    for n_eroded,prob in self.erosion_probs[erosion][gene_name].iteritems():
                        total += prob
                    assert math.fabs(total - 1.0) < 1e-5
    
    def write_erosion_probs(self):
        """ Write erosion_probs to disk. """
        print '  writing erosion probs'
        for region in self.regions:
            for erosion in self.erosions:
                outdir = self.base_outdir + '/' + region + '/' + erosion
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                for gene_name in self.gene_names[region]:
                    # replace * and / in the gene names so they can be made into file names
                    sanitized_gene_name = gene_name.replace('*',"_star_").replace('/','_slash_')
                    outfname = outdir + '/' + sanitized_gene_name + '.csv.bz2'
                    with bz2.BZ2File(outfname, 'w') as outfile:  # this compression is kind of pointless a.t.m.
                        out_fieldnames = [region + '_gene', 'n_eroded', 'prob']
                        writer = csv.writer(outfile, out_fieldnames)
                        # sort by n_eroded so file is more human readable
                        sorted_values = sorted(self.erosion_probs[erosion][gene_name].iteritems(),
                                               key=operator.itemgetter(0))
                        for n_eroded,prob in sorted_values:
                             writer.writerow([gene_name, str(n_eroded), str(prob)])
