#!/usr/bin/env python
import csv
import os
import glob
import sys

import utils
from opener import opener
from performanceplotter import PerformancePlotter

indir = 'data/performance/ihmmune/'
seqfname = 'caches/recombinator/longer-reads/simu.csv'
datadir = 'data/imgt'
plotdir = os.getenv('www') + '/partis/ihhhmmm-performance'
# header=('unique_id', 'v_gene', 'd_gene', 'j_gene', 'v_seq', 'd_seq', 'j_seq')

# first word    value name       pos of value
line_order = \
[
    ('State',     '',             4),
    ('IGHV',      'v_gene',       0),
    ('input:',    'v_qr_seq',     1),
    ('germline:', 'v_gl_seq',     1),
    ('match:',    '',            -1),
    ('V-D',       '',            -1),
    ('emitted',   'vd_insertion', 2),
    ('IGHD',      'd_gene',       0),
    ('input:',    'd_qr_seq',     1),
    ('germline:', 'd_gl_seq',     1),
    ('match:',    '',            -1),
    ('D-J',       '',            -1),
    ('emitted',   'dj_insertion', 2),
    ('IGHJ',      'j_gene',       0),
    ('input:',    'j_qr_seq',     1),
    ('germline:', 'j_gl_seq',     1),
    ('match:',    '',            -1)
]

# ----------------------------------------------------------------------------------------
def clean_value(column, value):
    if column == 'v_gene':
      return value[ : value.find('(')]
    elif column == 'j_gene':
        if 'P' in value:
            return value + '_P'
        else:
            return value + '_F'
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
    def increment(self):  # I feel like this would be easier with python's existing file class, but I haven't figured how to do it that way
        self.iline += 1
        self.line = self.lines[self.iline].strip().split()

# ----------------------------------------------------------------------------------------
class IhhhmmmParser(object):
    def __init__(self):
        self.debug = 0
        n_max_queries = 1
        queries = []

        self.germline_seqs = utils.read_germlines(datadir, add_fp=True)
        perfplotter = PerformancePlotter(self.germline_seqs, plotdir, 'ihhhmmm')

        # get sequence info that was passed to ihhhmmm
        self.seqinfo = {}
        with opener('r')(seqfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if len(queries) > 0 and line['unique_id'] not in queries:
                    continue
                self.seqinfo[line['unique_id']] = line
                iline += 1
                if n_max_queries > 0 and iline >= n_max_queries:
                    break

        infnames = glob.glob(indir + '/*.txt.fostream')
        for unique_id in self.seqinfo:
            utils.print_reco_event(self.germline_seqs, self.seqinfo[unique_id])
            for infname in infnames:
                print infname
                with opener('r')(infname) as infile:
                    details = self.parse_file(infile)
                    # reader = csv.DictReader(infile, delimiter=';')
                    # for line in reader:
                    #     line = line[None]
                    #     for ihead in range(len(header)):
                    #         print header[ihead], line[ihead]
                    #     sys.exit()

    # ----------------------------------------------------------------------------------------
    def parse_detail(self, fk):
        assert fk.iline < len(fk.lines)

        while fk.line[1] != 'Details':
            print '  skipping', fk.line
            fk.increment()

        fk.increment()
        info = {}
        for begin_line, column, index in line_order:
            if fk.line[0].find(begin_line) != 0:
                print 'oop', begin_line, fk.line
                sys.exit()
            if column != '':
                info[column] = clean_value(column, fk.line[index])
                print column, info[column]
            fk.increment()

        for region in utils.regions:
            if info[region + '_gene'] not in self.germline_seqs[region]:
                print 'ERROR %s not in germlines' % info[region + '_gene']
                sys.exit()
            if info[region + '_gl_seq'] not in self.germline_seqs[region][info[region + '_gene']]:
                print 'ERROR gl match not found in gl for %s' % info[region + '_gene']
                print '  ', info[region + '_gl_seq']
                print '  ', self.germline_seqs[region][info[region + '_gene']]                

        return info
        
    # ----------------------------------------------------------------------------------------
    def parse_file(self, infile):
        details = []
        fk = FileKeeper(infile.readlines())
        while fk.iline < len(fk.lines):
            details.append(self.parse_detail(fk))
            sys.exit()
        
iparser = IhhhmmmParser()
