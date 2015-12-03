#!/usr/bin/env python
import csv
import argparse
from collections import OrderedDict
import re
import os
from subprocess import check_call
import glob
import sys

import utils
from opener import opener
from performanceplotter import PerformancePlotter

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
    def __init__(self, args):
        self.args = args

        self.germline_seqs = utils.read_germlines(self.args.datadir, remove_N_nukes=True)
        self.perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir, 'ihhhmmm')

        self.details = OrderedDict()
        self.failtails = {}
        self.n_partially_failed = 0

        # get sequence info that was passed to ihhhmmm
        self.siminfo = OrderedDict()
        self.sim_need = []  # list of queries that we still need to find
        with opener('r')(self.args.simfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if self.args.queries != None and line['unique_id'] not in self.args.queries:
                    continue
                self.siminfo[line['unique_id']] = line
                self.sim_need.append(line['unique_id'])
                iline += 1
                if args.n_queries > 0 and iline >= args.n_queries:
                    break

        fostream_names = glob.glob(self.args.indir + '/*.fostream')
        if len(fostream_names) == 0:
            raise Exception('no fostreams found in %s' % args.indir)
        fostream_names.sort()  # maybe already sorted?
        for infname in fostream_names:
            if len(self.sim_need) == 0:
                break

            # try to get whatever you can for the failures
            unique_ids = self.find_partial_failures(infname)  # returns list of unique ids in this file

            with opener('r')(infname) as infile:
                self.parse_file(infile, unique_ids)

        # now check that we got results for all the queries we wanted
        n_failed = 0
        for unique_id in self.siminfo:
            if unique_id not in self.details and unique_id not in self.failtails:
                print '%-20s  no info' % unique_id
                self.perfplotter.add_fail()
                n_failed += 1

        print ''
        print 'partially failed: %d / %d = %.2f' % (self.n_partially_failed, len(self.siminfo), float(self.n_partially_failed) / len(self.siminfo))
        print 'failed:           %d / %d = %.2f' % (n_failed, len(self.siminfo), float(n_failed) / len(self.siminfo))
        print ''

        self.perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def parse_file(self, infile, unique_ids):
        fk = FileKeeper(infile.readlines())
        i_id = 0
        while not fk.eof and len(self.sim_need) > 0:
            self.parse_detail(fk, unique_ids[i_id])
            i_id += 1
        
    # ----------------------------------------------------------------------------------------
    def parse_detail(self, fk, unique_id):
        assert fk.iline < len(fk.lines)

        while fk.line[1] != 'Details':
            fk.increment()
            if fk.eof:
                return

        fk.increment()
        info = {}
        info['unique_id'] = unique_id
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
                # if '[' in info[column]:
                #     print 'added', column, clean_value(column, fk.line[index])
                if column.find('_gene') == 1:
                    region = column[0]
                    info[region + '_5p_del'] = int(fk.line[fk.line.index('start:') + 1]) - 1  # NOTE their indices are 1-based
                    gl_length = int(fk.line[fk.line.index('gene:') + 1]) - 1
                    match_end = int(fk.line[fk.line.index('end:') + 1]) - 1
                    assert gl_length >= match_end
                    info[region + '_3p_del'] = gl_length - match_end

            fk.increment()

        if unique_id not in self.sim_need:
            while not fk.eof and fk.line[1] != 'Details':  # skip stuff until start of next Detail block
                fk.increment()
            return

        info['fv_insertion'] = ''
        info['jf_insertion'] = ''
        info['seq'] = info['v_qr_seq'] + info['vd_insertion'] + info['d_qr_seq'] + info['dj_insertion'] + info['j_qr_seq']

        if '-' in info['seq']:
            print 'ERROR found a dash in %s, returning failure' % unique_id
            while not fk.eof and fk.line[1] != 'Details':  # skip stuff until start of next Detail block
                fk.increment()
            return

        if info['seq'] not in self.siminfo[unique_id]['seq']:  # arg. I can't do != because it tacks on v left and j right deletions
            print 'ERROR didn\'t find the right sequence for %s' % unique_id
            print '  ', info['seq']
            print '  ', self.siminfo[unique_id]['seq']
            sys.exit()

        if self.args.debug:
            print unique_id
            for region in utils.regions:
                infer_gene = info[region + '_gene']
                true_gene = self.siminfo[unique_id][region + '_gene']
                if utils.are_alleles(infer_gene, true_gene):
                    regionstr = utils.color('bold', utils.color('blue', region))
                    truestr = ''  #'(originally %s)' % match_name
                else:
                    regionstr = utils.color('bold', utils.color('red', region))
                    truestr = '(true: %s)' % utils.color_gene(true_gene).replace(region, '')
                print '  %s %s %s' % (regionstr, utils.color_gene(infer_gene).replace(region, ''), truestr)

            utils.print_reco_event(self.germline_seqs, self.siminfo[unique_id], label='true:', extra_str='    ')
            utils.print_reco_event(self.germline_seqs, info, label='inferred:', extra_str='    ')

        for region in utils.regions:
            if info[region + '_gene'] not in self.germline_seqs[region]:
                print 'ERROR %s not in germlines' % info[region + '_gene']
                assert False

            gl_seq = info[region + '_gl_seq']
            if '[' in gl_seq:  # ambiguous
                for nuke in utils.nukes:
                    gl_seq = gl_seq.replace('[', nuke)
                    if gl_seq in self.germline_seqs[region][info[region + '_gene']]:
                        print '  replaced [ with %s' % nuke
                        break
                info[region + '_gl_seq'] = gl_seq

            if info[region + '_gl_seq'] not in self.germline_seqs[region][info[region + '_gene']]:
                print 'ERROR gl match not found for %s in %s' % (info[region + '_gene'], unique_id)
                print '  ', info[region + '_gl_seq']
                print '  ', self.germline_seqs[region][info[region + '_gene']]                
                self.perfplotter.add_partial_fail(self.siminfo[unique_id], info)
                while not fk.eof and fk.line[1] != 'Details':  # skip stuff until start of next Detail block
                    fk.increment()
                return

        self.perfplotter.evaluate(self.siminfo[unique_id], info)
        self.details[unique_id] = info
        self.sim_need.remove(unique_id)

        while not fk.eof and fk.line[1] != 'Details':  # skip stuff until start of next Detail block
            fk.increment()
        
    # ----------------------------------------------------------------------------------------
    def find_partial_failures(self, fostream_name):
        unique_ids = []
        for line in open(fostream_name.replace('.fostream', '')).readlines():
            if len(self.sim_need) == 0:
                return
            if len(line.strip()) == 0:  # skip blank lines
                continue

            line = line.replace('"', '')
            line = line.split(';')

            unique_id = line[0]
            
            if 'NA' not in line:  # skip lines that were ok
                unique_ids.append(unique_id)
                continue
            if unique_id not in self.sim_need:
                continue
            if unique_id not in self.siminfo:
                continue  # not looking for this <unique_id> a.t.m.

            info = {}
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
            if self.args.debug:
                print '%-20s  partial fail %s %s %s' % (unique_id,
                                                     utils.color_gene(info['v_gene']) if 'v_gene' in info else '',
                                                     utils.color_gene(info['d_gene']) if 'd_gene' in info else '',
                                                     utils.color_gene(info['j_gene']) if 'j_gene' in info else ''),
                print '  (true %s %s %s)' % tuple([self.siminfo[unique_id][region + '_gene'] for region in utils.regions])
            self.failtails[unique_id] = info
            self.n_partially_failed += 1
            self.sim_need.remove(unique_id)

        return unique_ids

# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-queries', type=int, default=-1)
    parser.add_argument('--queries')
    parser.add_argument('--plotdir', required=True)
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
    parser.add_argument('--datadir', default='data/imgt')
    parser.add_argument('--indir', required=True)  # data/performance/ihhhmmm
    parser.add_argument('--simfname')
    args = parser.parse_args()
    args.queries = utils.get_arg_list(args.queries)
    
    # check_call(['tar', 'xzf', args.indir + '.tgz', '-C', 'data/performance/'])  # untar the ihmmune-align output
    ihhhmmmparser = IhhhmmmParser(args)
