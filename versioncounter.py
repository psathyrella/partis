#!/usr/bin/env python

""" Read observed counts of the vdj gene chosen and deletion lengths.

Store them in a similar file, but with the counts normalized.
"""

import bz2
import csv
import subprocess
import math
import os

#erosions = ['v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del']
regions = ['v', 'd', 'j']
erosions = ['5p', '3p']
gene_names = {}

min_counts = 3  # if a combo has fewer counts than this, ignore it completely

def parse_input(infname, version_freqs, erosion_counts):
    """ get sum of the \'count\' column in <infname>. """
    print '  parsing input file %s' % infname
    total = 0.0
    with bz2.BZ2File(infname) as infile:
        in_data = csv.DictReader(infile)
#        progress_count = 0
#        n_lines = 5000000 #subprocess.check_output(['bzcat ', infname, ' | wc'])
        for line in in_data:
#            if progress_count % 100000 == 0:
#                print '  %3.1e / %3.0e = %5.2f  %s' % (progress_count, n_lines, float(progress_count) / n_lines, line['count'])
#            progress_count += 1
            if int(line['count']) < min_counts:
                print '    breaking because count is less than %d. \n    NOTE that this assumes input file is sorted!' % min_counts
                break
            total += int(line['count'])
            index = (line['v_gene'], line['d_gene'], line['j_gene'], line['cdr3_length'])
            if index not in version_freqs:
                version_freqs[index] = 0
            version_freqs[index] += int(line['count'])
            for region in regions:
                if line[region + '_gene'] not in gene_names[region]:
                    gene_names[region].add(line[region + '_gene'])
                for erosion in erosions:
                    erosion_name = region + '_' + erosion + '_del'
                    gene_name = line[region + '_gene']
                    if gene_name not in erosion_counts[erosion]:
                        erosion_counts[erosion][gene_name] = {}
                    if line[erosion_name] not in erosion_counts[erosion][gene_name]:  # e.g.: if we haven't already seen IGHV3-64*04 with an erosion of length 4 on the 3p end
                        erosion_counts[erosion][gene_name][line[erosion_name]] = 0
                    erosion_counts[erosion][gene_name][line[erosion_name]] += int(line['count'])

    return total

def write_gene_choice_probs(infname, outfname, total, version_freqs):
    print '  writing gene choice probabilities'
    with bz2.BZ2File(outfname, 'w') as outfile:
        out_fieldnames = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'prob']
        out_data = csv.DictWriter(outfile, out_fieldnames)
        out_data.writeheader()
        check_tot = 0.0
        for index,count in version_freqs.iteritems():
            prob = float(count) / total
            line = {}
            line['v_gene'] = index[0]
            line['d_gene'] = index[1]
            line['j_gene'] = index[2]
            line['cdr3_length'] = index[3]
            line['prob'] = prob
            check_tot += prob
            out_data.writerow(line)
        print 'check: %f' % check_tot

def normalize_erosion_counts(erosion_counts):
    print '  normalizing erosion counts'
    erosion_probs = {}
    for erosion in erosions:
        erosion_probs[erosion] = {}
        for region in regions:
            for gene_name in gene_names[region]:
                erosion_probs[erosion][gene_name] = {}
                assert gene_name in erosion_counts[erosion]
                total = 0.0  # total for this gene and this side (e.g. 3p of IGHV3-64*04)
                for n_eroded,count in erosion_counts[erosion][gene_name].iteritems():
                    total += count
                for n_eroded,count in erosion_counts[erosion][gene_name].iteritems():
                    erosion_probs[erosion][gene_name][n_eroded] = float(count) / total
    return erosion_probs

def check_normalization(erosion_probs):                    
    # SILLY little double check
    for erosion in erosions:
        for region in regions:
            for gene_name in gene_names[region]:
#                print region,erosion,gene_name
                total = 0.0
                for n_eroded,prob in erosion_probs[erosion][gene_name].iteritems():
                    total += prob
#                    print '  ',n_eroded,prob,total
                assert math.fabs(total - 1.0) < 1e-5

def write_erosion_probs(basedir, erosion_probs):
    print '  writing erosion probs'
    for region in regions:
        for erosion in erosions:
            outdir = basedir + '/' + region + '/' + erosion
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            for gene_name in gene_names[region]:
                sanitized_gene_name = gene_name.replace('*',"_star_").replace('/','_slash_')
                outfname = outdir + '/' + sanitized_gene_name + '.csv.bz2'
                with bz2.BZ2File(outfname, 'w') as outfile:  # this compression is super pointless
                    out_fieldnames = [region + '_gene', 'n_eroded', 'prob']
                    writer = csv.writer(outfile, out_fieldnames)
                    for n_eroded,prob in erosion_probs[erosion][gene_name].iteritems():
                        writer.writerow([gene_name, str(n_eroded), str(prob)])

version_freqs = {}
erosion_counts = {}  # Dict of erosion side, gene name, and number eroded to frequency, e.g.: erosion_counts['5p']['IGHV3-64*04']['3'] = 35
                     #   means that for that gene, we saw 5p erosions of length 3 35 times
for region in regions:
    gene_names[region] = set()
    erosion_counts[region] = {}
for erosion in erosions:
    erosion_counts[erosion] = {}

infname = 'data/human-beings/01-C-N_filtered.vdjcdr3.csv.bz2'
total = parse_input(infname, version_freqs, erosion_counts)
write_gene_choice_probs(infname, infname.replace('.csv.bz2','.probs.csv.bz2'), total, version_freqs)
erosion_probs = normalize_erosion_counts(erosion_counts)
check_normalization(erosion_probs)
write_erosion_probs('data/human-beings', erosion_counts)
