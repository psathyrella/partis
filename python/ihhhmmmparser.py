#!/usr/bin/env python
import csv
from collections import OrderedDict
import os
from subprocess import check_output
import glob
import sys

import utils
from opener import opener
from performanceplotter import PerformancePlotter

indir = 'data/performance/ihmmune/'
seqfname = 'caches/recombinator/longer-reads/simu.csv'
datadir = 'data/imgt'
plotdir = os.getenv('www') + '/partis/ihhhmmm-performance'

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
        if '_HM855939' in value:
            value = value[ : value.find('_HM855939')]  # what the hell, dude?
            value = value.replace('NL*', 'NL1*')  # and again I ask, what the hell?
        if '(' in value:
            value = value[ : value.find('(')]
        return value
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

        self.germline_seqs = utils.read_germlines(datadir)  #, add_fp=True)
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
            # if self.debug:
            #     print '  skipping', fk.line
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

        for unique_id in self.seqinfo:
            if self.seqinfo[unique_id]['seq'] == info['seq']:
                info['unique_id'] = unique_id
                break

        if 'unique_id' not in info:  # arg. probabl a j right or v left erosion made it not match
            grep_match = check_output(['grep', '-Hrn', info['seq'], seqfname])
            info['unique_id'] = grep_match.split(':')[2].split(',')[0]
            if info['unique_id'] not in self.seqinfo:
                print 'ERROR %s not in seqinfo' % info['unique_id']

        utils.print_reco_event(self.germline_seqs, self.seqinfo[info['unique_id']], label='true:')
        utils.print_reco_event(self.germline_seqs, info, label='inferred:')

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
