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
        self.debug = 1
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

        n_failed = 0
        with opener('r')(infname) as infile:
            soup = BeautifulSoup(infile)
            for unique_id in self.seqinfo:
                print unique_id
                imgtinfo = []
                for pre in soup.find_all('pre'):  # NOTE this loops over everything an awful lot of times. shouldn't really matter for now
                    if unique_id in pre.text:
                        imgtinfo.append(pre.text)
                line = self.parse_query(unique_id, imgtinfo)
                if len(line) == 0:
                    print '    giving up'
                    n_failed += 1
                joinparser.add_insertions(line)
                joinparser.resolve_overlapping_matches(line, debug=True)
                # perfplotter.evaluate(self.seqinfo[unique_id], line, unique_id)
                if self.debug:
                    utils.print_reco_event(self.germline_seqs, line)

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

        full_qr_seq = query_info[0].replace('>', '').replace(unique_id, '')  # strip off the unique id
        full_qr_seq = ''.join(full_qr_seq.split()).upper()  # strip off white space and uppercase it
        assert full_qr_seq == self.seqinfo[unique_id]['seq']

        line = {}
        line['unique_id'] = unique_id
        line['seq'] = full_qr_seq
        # qrbounds = {}
        for ireg in range(len(utils.regions)):
            region = utils.regions[ireg]
            info = query_info[ireg + 1].splitlines()
            if unique_id not in info[0]:  # remove the line marking cdr3 and framework regions
                info.pop(0)
            qr_seq = info[0].split()[1].upper()  # this line should be '<unique_id> .............<query_seq>'
            match_name = str(info[1].split()[2])
            gl_seq = info[1].split()[4].upper()
            assert qr_seq.replace('.', '') in self.seqinfo[unique_id]['seq']

            if self.debug:
                print '  ', region, match_name
                print '    qr', qr_seq
                print '      ', gl_seq

            # replace the dots (gaps) in the gl match
            new_qr_seq, new_gl_seq = [], []
            for inuke in range(min(len(qr_seq), len(gl_seq))):
                if gl_seq[inuke] == '.':
                    pass
                else:
                    new_qr_seq.append(qr_seq[inuke])  # this should only be out of range if the v match extends through the whole query sequence, i.e. friggin never
                    new_gl_seq.append(gl_seq[inuke])
            for inuke in range(len(gl_seq), len(qr_seq)):
                new_qr_seq.append(qr_seq[inuke])
            for inuke in range(len(qr_seq), len(gl_seq)):
                new_gl_seq.append(gl_seq[inuke])
            qr_seq = ''.join(new_qr_seq)
            gl_seq = ''.join(new_gl_seq)

            # work out the erosions
            qr_ldots = qr_seq.rfind('.') + 1  # first strip off any dots on the left of query seq
            qr_seq = qr_seq[qr_ldots : ]
            gl_seq = gl_seq[qr_ldots : ]
            gl_ldots = gl_seq.rfind('.') + 1  # then remove dots on the left of the germline seq
            qr_seq = qr_seq[gl_ldots : ]
            gl_seq = gl_seq[gl_ldots : ]
            del_5p = qr_ldots + gl_ldots
            qr_seq = qr_seq[ : len(gl_seq)]  # then strip the right-hand portion of the query sequence that isn't aligned to the germline
            del_3p = len(gl_seq) - len(qr_seq)  # then do the same for the germline overhanging on the right of the query
            gl_seq = gl_seq[ : len(qr_seq)]
            assert len(gl_seq) == len(qr_seq)
            new_gl_seq = []
            for inuke in range(len(gl_seq)):  # replace dashes (matched bases)
                assert gl_seq[inuke] != '.'  # hoping there's no gaps in here
                if gl_seq[inuke] == '-':
                    new_gl_seq.append(qr_seq[inuke])
                else:
                    new_gl_seq.append(gl_seq[inuke])
            gl_seq = ''.join(new_gl_seq)
            if self.debug:
                print '    qr', qr_seq
                print '      ', gl_seq, del_5p, del_3p

            assert len(re.findall(qr_seq, full_qr_seq)) == 1
            qr_start = full_qr_seq.find(qr_seq)
            assert qr_start >= 0
            # qrbounds[region] = (qr_start, qr_start + len(qr_seq))
            match_name = joinparser.figure_out_which_damn_gene(self.germline_seqs, match_name, gl_seq, debug=self.debug)

            adaptive_gl_seq = self.germline_seqs[region][match_name]
            if region == 'j':  # remove the extra righthand bases in the imgt version
                assert adaptive_gl_seq[del_5p : ].find(gl_seq) == 0  # left hand side of the two should be the same now
                assert len(adaptive_gl_seq[del_5p : ]) == len(gl_seq)  # should be ok for now
                del_3p = 0
            if gl_seq != adaptive_gl_seq[del_5p : len(adaptive_gl_seq) - del_3p]:
                print 'ERROR doesn\'t match adaptive gl version'
                print 'imgt              ', gl_seq
                print 'adaptive          ', adaptive_gl_seq[del_5p : len(adaptive_gl_seq) - del_3p]
                print 'adaptive untrimmed', adaptive_gl_seq
                sys.exit()
            line[region + '_gene'] = match_name
            line[region + '_qr_seq'] = qr_seq
            line[region + '_gl_seq'] = gl_seq
            line[region + '_5p_del'] = del_5p
            line[region + '_3p_del'] = del_3p
            
        return line
# joinparser.figure_out_which_damn_gene(self.germline_seqs, 
#             if match_names[region] not in self.germline_seqs[region]:
#                 print 'ERROR %s not found in germline file' % match_names[region]
#                 sys.exit()

iparser = IMGTParser('caches/recombinator/simu.csv', '/home/dralph/Dropbox/imgtvquest.html', datadir='./data')
