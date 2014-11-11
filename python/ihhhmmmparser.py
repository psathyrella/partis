#!/usr/bin/env python
import csv
from collections import OrderedDict
import os
from subprocess import check_call
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
        n_max_queries = 100
        queries = []

        self.germline_seqs = utils.read_germlines(datadir, add_fp=True)
        perfplotter = PerformancePlotter(self.germline_seqs, plotdir, 'ihhhmmm')

        # get sequence info that was passed to ihhhmmm
        self.seqinfo = OrderedDict()
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
        infnames.sort()
        for infname in infnames:
            print infname
            with opener('r')(infname) as infile:
                details = self.parse_file(infile)

    # ----------------------------------------------------------------------------------------
    def parse_detail(self, fk):
        assert fk.iline < len(fk.lines)

        while fk.line[1] != 'Details':
            print '  skipping', fk.line
            fk.increment()

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
                if column == 'j_gene'):
                    gl_length = int(fk[fk.line.index('gene:') + 1])
                    match_end = int(fk[fk.line.index('index:') + 1])
                    if match_end != gl_length:
                        print '    j match not to end, expanding %d --> %d' % (match_end, gl_end)
                    asdfjkl;asdfjkl;sdgasghil;
                    
                print column, info[column]
            fk.increment()

        info['seq'] = info['v_qr_seq'] + info['vd_insertion'] + info['d_qr_seq'] + info['dj_insertion'] + info['j_qr_seq']

        for unique_id in self.seqinfo:
            if self.seqinfo[unique_id]['seq'] == info['seq']:
                info['unique_id'] = unique_id
                break

        if 'unique_id' not in info:
            print 'seq not found'
            check_call(['grep', '-Hrn', info['seq'], seqfname])
            print '  ', info['seq']
            sys.exit()

        # for unique_id in self.seqinfo:
        #     utils.print_reco_event(self.germline_seqs, self.seqinfo[unique_id])

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
        
iparser = IhhhmmmParser()
