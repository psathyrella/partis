#!/usr/bin/env python
import csv
from collections import OrderedDict
import re
import os
from subprocess import check_output
import glob
import sys

import utils
from opener import opener
from performanceplotter import PerformancePlotter

label='new-imgt'
indir = 'caches/recombinator/performance/' + label
simfname = indir + '/simu.csv'
datadir = 'data/imgt'
plotdir = os.getenv('www') + '/partis/ihhhmmm-performance/' + label

# first word    value name     position   required?
line_order = \
[
    ('State',     '',            -1,   True,    None),
    ('IGHV',      'v_gene',       0,   True,    None),
    ('input:',    'v_qr_seq',     1,   True,    None),
    ('germline:', 'v_gl_seq',     1,   True,    None),
    ('match:',    '',            -1,   True,    None),
    ('V-D',       '',            -1,   False,   None),
    ('emitted',   'vd_insertion', 2,   False,   ''  ),
    ('IGHD',      'd_gene',       0,   True,    None),
    ('input:',    'd_qr_seq',     1,   True,    None),
    ('germline:', 'd_gl_seq',     1,   True,    None),
    ('match:',    '',            -1,   True,    None),
    ('D-J',       '',            -1,   False,   None),
    ('emitted',   'dj_insertion', 2,   False,   ''  ),
    ('IGHJ',      'j_gene',       0,   True,    None),
    ('input:',    'j_qr_seq',     1,   True,    None),
    ('germline:', 'j_gl_seq',     1,   True,    None),
    ('match:',    '',            -1,   True,    None)
]

# ----------------------------------------------------------------------------------------
def clean_value(column, value):
    if column == 'v_gene':
        if '_HM855577' in value:
            value = value[ : value.find('_HM855577')]
        if '_HM855939' in value:
            value = value[ : value.find('_HM855939')]  # what the hell, dude?
            value = value.replace('NL*', 'NL1*')  # and again I ask, what the hell?
        if '(' in value:
            value = value[ : value.find('(')]
        return value
    # elif column == 'j_gene':
    #     if 'P' in value:
    #         return value + '_P'
    #     else:
    #         return value + '_F'
    elif '_seq' in column or '_insertion' in column:
      return value.upper()
    else:
      return value

# ----------------------------------------------------------------------------------------
class FileKeeper(object):
    def __init__(self, lines):
        self.lines = lines
        self.iline = 0
        self.line = self.lines[self.iline].strip().split()
        self.eof = False
    def increment(self):  # I feel like this would be easier with python's existing file class, but I haven't figured how to do it that way
        self.iline += 1
        if self.iline >= len(self.lines):
            self.eof = True
            return
            # print 'ERROR iline %d too large (length %d)' % (self.iline, len(self.lines))
            # assert False
        self.line = self.lines[self.iline].strip().split()

# ----------------------------------------------------------------------------------------
class IhhhmmmParser(object):
    def __init__(self):
        self.debug = 1
        n_max_queries = 10
        queries = []

        self.germline_seqs = utils.read_germlines(datadir)  #, add_fp=True)
        self.perfplotter = PerformancePlotter(self.germline_seqs, plotdir, 'ihhhmmm')

        self.details = OrderedDict()

        # get sequence info that was passed to ihhhmmm
        self.siminfo = OrderedDict()
        with opener('r')(simfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if len(queries) > 0 and line['unique_id'] not in queries:
                    continue
                self.siminfo[line['unique_id']] = line
                iline += 1
                if n_max_queries > 0 and iline >= n_max_queries:
                    break

        # first get info for the queries on which it succeeded
        fostream_names = glob.glob(indir + '/*.txt.fostream')
        fostream_names.sort()
        for infname in fostream_names:
            print infname
            with opener('r')(infname) as infile:
                self.parse_file(infile)

        # then try to get whatever you can for the failures
        txtfnames = [foname.replace('.fostream', '') for foname in fostream_names]
        failtails = {}
        n_partially_failed = 0
        for fname in txtfnames:
            for line in open(fname).readlines():
                if len(line.strip()) == 0:  # skip blank lines
                    continue
                if 'NA' not in line:  # skip lines that were ok
                    continue
                line = line.replace('"', '')
                line = line.split(';')
                info = {}
                unique_id = line[0]
                info['unique_id'] = unique_id
                for stuff in line:
                    for region in utils.regions:  # add the first instance of IGH[VDJ] (if it's there at all)
                        if 'IGH'+region.upper() in stuff and region+'_gene' not in info:
                            genes = re.findall('IGH' + region.upper() + '[^ ][^ ]*', stuff)
                            if len(genes) == 0:
                                print 'ERROR no %s genes in %s' % (region, stuff)
                            gene = genes[0]
                            if gene not in self.germline_seqs[region]:
                                print 'ERROR bad gene %s for %s' % (gene, unique_id)
                                sys.exit()
                            info[region + '_gene'] = gene
                self.perfplotter.add_partial_fail(self.siminfo[unique_id], info)
                if self.debug:
                    print '%-20s  partial fail %s %s %s' % (unique_id,
                                                         utils.color_gene(info['v_gene']) if 'v_gene' in info else '',
                                                         utils.color_gene(info['d_gene']) if 'd_gene' in info else '',
                                                         utils.color_gene(info['j_gene']) if 'j_gene' in info else ''),
                    print '  (true %s %s %s)' % tuple([self.siminfo[unique_id][region + '_gene'] for region in utils.regions])
                failtails[unique_id] = info
                n_partially_failed += 1

        # now check that we got results for all the queries we wanted
        n_failed = 0
        for unique_id in self.siminfo:
            if unique_id not in self.details and unique_id not in failtails:
                print '%-20s  no info' % unique_id
                self.perfplotter.add_fail()
                n_failed += 1

        print ''
        print 'partially failed: %d / %d = %f' % (n_partially_failed, len(self.siminfo), float(n_partially_failed) / len(self.siminfo))
        print 'failed: %d / %d = %f' % (n_failed, len(self.siminfo), float(n_failed) / len(self.siminfo))
        print ''

        self.perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def parse_detail(self, fk):
        assert fk.iline < len(fk.lines)

        while fk.line[1] != 'Details':
            fk.increment()
            if fk.eof:
                return

        fk.increment()
        info = {}
        for begin_line, column, index, required, default in line_order:
            if fk.line[0].find(begin_line) != 0:
                if required:
                    print 'oop', begin_line, fk.line
                    sys.exit()
                else:
                    info[column] = default
                    continue
            if column != '':
                info[column] = clean_value(column, fk.line[index])
                if column.find('_gene') == 1:
                    region = column[0]
                    info[region + '_5p_del'] = int(fk.line[fk.line.index('start:') + 1]) - 1  # NOTE their indices are 1-based
                    gl_length = int(fk.line[fk.line.index('gene:') + 1]) - 1
                    match_end = int(fk.line[fk.line.index('end:') + 1]) - 1
                    assert gl_length >= match_end
                    info[region + '_3p_del'] = gl_length - match_end

            fk.increment()

        info['fv_insertion'] = ''
        info['jf_insertion'] = ''
        info['seq'] = info['v_qr_seq'] + info['vd_insertion'] + info['d_qr_seq'] + info['dj_insertion'] + info['j_qr_seq']

        for unique_id in self.siminfo:
            if self.siminfo[unique_id]['seq'] == info['seq']:
                info['unique_id'] = unique_id
                break

        if 'unique_id' not in info:  # arg. probably a j right or v left erosion made it not match
            grep_match = check_output(['grep', '-Hrn', info['seq'], simfname])
            info['unique_id'] = grep_match.split(':')[2].split(',')[0]
            if info['unique_id'] not in self.siminfo:
                print 'ERROR %s not in siminfo' % info['unique_id']

        if self.debug:
            print info['unique_id']
            utils.print_reco_event(self.germline_seqs, self.siminfo[info['unique_id']], label='true:', extra_str='    ')
            utils.print_reco_event(self.germline_seqs, info, label='inferred:', extra_str='    ')

        for region in utils.regions:
            if info[region + '_gene'] not in self.germline_seqs[region]:
                print 'ERROR %s not in germlines' % info[region + '_gene']
                sys.exit()
            if info[region + '_gl_seq'] not in self.germline_seqs[region][info[region + '_gene']]:
                print 'ERROR gl match not found in gl for %s' % info[region + '_gene']
                print '  ', info[region + '_gl_seq']
                print '  ', self.germline_seqs[region][info[region + '_gene']]                

        self.perfplotter.evaluate(self.siminfo[info['unique_id']], info)
        self.details[info['unique_id']] = info
        
    # ----------------------------------------------------------------------------------------
    def parse_file(self, infile):
        fk = FileKeeper(infile.readlines())
        while fk.iline < len(fk.lines):
            self.parse_detail(fk)
        
iparser = IhhhmmmParser()
