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
    Output: probabilities for each (vdj gene choice + cdr3_length combo + erosion length), eg: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42, 0, 2, 5, 2) : 3.4e-5
    """
    min_counts = 10  # if a combo has fewer counts than this, ignore it completely

    version_freqs = {}  # Dict with entries like: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42) : 360, i.e. we saw this combination of genes with cdr3_length 42 a total of 360 times
    regions = ['v', 'd', 'j']
    erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']

    def __init__(self, infname, base_outdir):
        self.infname = infname
        self.base_outdir = base_outdir

    def count(self):
        """ Run the whole counting procedure. """
        total = self.parse_input()
        self.write_gene_choice_probs(total)
 
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
                index = (line['v_gene'], line['d_gene'], line['j_gene'], line['cdr3_length'], line['v_3p_del'], line['d_5p_del'], line['d_3p_del'], line['j_5p_del'])
                if index not in self.version_freqs:  # duplicates are from non-sensical 'v_5p_del' and 'j_3p_del' erosions
                    self.version_freqs[index] = 0
                self.version_freqs[index] += int(line['count'])
    
        return total
    
    def write_gene_choice_probs(self, total):
        """ Write vdj combo freqs to file. """
        outfname = self.base_outdir + '/' + os.path.basename(self.infname.replace('.csv.bz2','.probs.csv.bz2'))
        print '  writing gene choice probabilities to %s' % outfname
        with bz2.BZ2File(outfname, 'w') as outfile:
            out_fieldnames = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'prob']
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
                line['v_3p_del'] = index[4]
                line['d_5p_del'] = index[5]
                line['d_3p_del'] = index[6]
                line['j_5p_del'] = index[7]
                line['prob'] = prob
                check_tot += prob
                out_data.writerow(line)
