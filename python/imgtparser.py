#!/usr/bin/env python
import csv
import sys
import os
import re
from bs4 import BeautifulSoup

from opener import opener
import utils
import joinparser

# from performanceplotter import PerformancePlotter

class IMGTParser(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, seqfname, infname, datadir):
        self.debug = 0
        n_max_queries = 10
        queries = []

        self.germline_seqs = utils.read_germlines(datadir, remove_N_nukes=False)
        # perfplotter = PerformancePlotter(self.germline_seqs, os.getenv('www') + '/partis/imgt_performance', 'imgt')

        # get info that was passed to joinsolver
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

        self.info = {}
        n_failed = 0
        with opener('r')(infname) as infile:
            soup = BeautifulSoup(infile)
            for unique_id in self.seqinfo:
                print unique_id
                imgtinfo = []
                for pre in soup.find_all('pre'):  # NOTE this loops over everything an awful lot of times. shouldn't really matter for now
                    if unique_id in pre.text:
                        imgtinfo.append(pre.text)
                self.info[unique_id] = self.parse_query(unique_id, imgtinfo)
                if len(self.info[unique_id]) == 0:
                    print '    giving up'
                    n_failed += 1

    # ----------------------------------------------------------------------------------------
    def parse_query(self, unique_id, query_info):
        if len(query_info) == 0:  # one for the query sequence, then one for v, d, and j
            print 'no info for',unique_id
            return {}
        elif len(query_info) < 4:
            regions_ok = ''
            for info in query_info:
                for region in utils.regions:
                    if 'IGH' + region.upper() in info:
                        regions_ok += region
            for region in utils.regions:
                if region not in regions_ok:
                    print '    ERROR no %s matches' % region
                    return {}
            assert False  # shouldn't get here
        elif len(query_info) != 4:
            print 'info for', unique_id, 'all messed up'
            for info in query_info:
                print info
            sys.exit()

        query_seq = query_info[0].replace('>', '').replace(unique_id, '')  # strip off the unique id
        query_seq = ''.join(query_seq.split()).upper()  # strip off white space and uppercase it
        assert query_seq == self.seqinfo[unique_id]['seq']

        qr_seqs, match_names, gl_seqs = {}, {}, {}
        for ireg in range(len(utils.regions)):
            region = utils.regions[ireg]
            print '  ', region
            info = query_info[ireg + 1].splitlines()
            if unique_id not in info[0]:  # remove the line marking cdr3 and framework regions
                info.pop(0)
            qr_seq = info[0].split()[1]  # this line should be '<unique_id> .............<query_seq>'
            match_names[region] = info[1].split()[2]
            gl_seq = info[1].split()[4]
            assert qr_seq.replace('.', '').upper() in self.seqinfo[unique_id]['seq']

            # work out the erosions
            print '    qr', qr_seq
            print '      ', gl_seq
            qr_ldots = qr_seq.rfind('.') + 1  # first strip off any dots on the left of query seq
            qr_seq = qr_seq[qr_ldots : ]
            gl_seq = gl_seq[qr_ldots : ]
            print '    qr', qr_seq
            print '      ', gl_seq
            gl_ldots = gl_seq.rfind('.') + 1  # then remove dots on the left of the germline seq
            qr_seq = qr_seq[gl_ldots : ]
            gl_seq = gl_seq[gl_ldots : ]
            del_5p = qr_ldots + gl_ldots
            print '    qr', qr_seq
            print '      ', gl_seq
            qr_seq = qr_seq[ : len(gl_seq)]  # then strip the right-hand portion of the query sequence that isn't aligned to the germline
            print '    qr', qr_seq
            print '      ', gl_seq
            del_3p = len(gl_seq) - len(qr_seq)  # then do the same for the germline overhanging on the right of the query
            gl_seq = gl_seq[ : len(qr_seq)]
            assert len(gl_seq) == len(qr_seq)
            print '    qr', qr_seq
            print '      ', gl_seq, del_5p, del_3p

            gl_seqs[region] = gl_seq
            qr_seqs[region] = qr_seq

        line = {}
        line['seq'] = query_seq
        return line
# joinparser.figure_out_which_damn_gene(self.germline_seqs, 
#             if match_names[region] not in self.germline_seqs[region]:
#                 print 'ERROR %s not found in germline file' % match_names[region]
#                 sys.exit()

iparser = IMGTParser('caches/recombinator/simu.csv', '/home/dralph/Dropbox/imgtvquest.html', datadir='./data')
